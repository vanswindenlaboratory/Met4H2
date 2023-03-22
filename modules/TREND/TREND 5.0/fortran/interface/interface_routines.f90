

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


    ! module for file interface_routines_module.f90
    module interface_routines_module
    !global use inclusion

    use module_all_types
    use controlling
    use calc_functions
    use module_parameters
    use module_all_types
    use omp_lib
    use main_Data
    use module_eqn_props
    use module_matrices
    use setup_module
    use interface_support_module
    use rhomix_pt_module
    use ancillary_equations_mix_module
    use constraints
    use unit_convertion
    contains

    !double precision function TREND_EOS_STDCALL (calctype, input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "TREND_EOS_STDCALL" :: TREND_EOS_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: TREND_EOS_STDCALL
    !
    !
    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !
    !!define input variables
    !character (12), intent(in) :: input
    !character(12), intent(in) :: calctype
    !character(20) :: unitdefinition
    !double precision, intent(in) :: prop1, prop2
    !character(30), dimension(30), intent(in) :: fluids
    !double precision, dimension (30), intent(in) :: moles
    !character (255), intent(in) :: path
    !integer, dimension (30), intent(in) :: EOS_indicator
    !integer, intent(in) :: MIX_indicator
    !integer :: errorflag
    !
    !double precision :: TREND_EOS
    !
    !
    !TREND_EOS_STDCALL = TREND_EOS (calctype, input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, errorflag, handle)
    !
    !end function
    !
    !
    !
    !
    !!##################################################################
    !!TREND_EOS function calculates a thermophysical property depending on the input
    !!calctype       character(12)           specification of the property that has to be calculated possible calculations are: see if-block below
    !!input          character(12)           specifies the given state point where the property has to be calculated, possible inputs are:
    !!                                       TP (prop1 is temperature, prop2 is pressure)
    !!                                       TD (prop1 is temperature, prop2 is density)
    !!                                       PS(prop1 is pressure, prop2 is entropy)
    !!                                       PH (prop1 is pressure, prop2 is enthalpy)
    !!                                       TLIQ (prop1 is temperature on bubble line, prop2 is ignored)
    !!                                       TVAP(prop1 is temperature on dew line, prop2 is ignored)
    !!                                       PLIQ (prop1 is pressure on bubble line, prop2 is ignored)
    !!                                       PVAP (prop1 is pressure on dew line, prop2 is ignored)
    !!prop1, prop2   double                  values for state point depending on input
    !!fluids         character(30),dim(30)   fluid/mixture calculations are conducted with
    !!moles          double,dim(30)          composition  calculations are conducted with
    !!EOS_indicator  integer,dim(30)         equation of state indicator calculations are conducted with
    !!                                       1: Helmholtz
    !!                                       2: SRK
    !!                                       3: PR
    !!                                       4: LKP
    !!                                       5: Generalized EOS (51: based on Alexandrov 52: based on Span/Wagner 53: based on Sun/Ely)
    !!                                       6: PCSAFT
    !!                                       7:
    !!                                       8: RKM
    !!                                       9: COSTALD
    !!MIX_indicator  integer                 indicator for mixing rule
    !!                                       1:      Helmholtz default,
    !!                                       110:    All components are mixed according to the linear mixing rule
    !!                                       111:    check if binary mix files are available, if not use linear mixing rule
    !!                                       120:    All components are mixed according to the Lorentz-Berthelot mixing rule.
    !!                                       121:    check if binary mix files are available, if not use Lorentz-Berthelot mixing rule.
    !!                                       2:      mixing rule according to SRK model
    !!                                       3:      mixing rule according to PR model
    !!                                       4:      mixing rule according to LKP model
    !!                                       5:      mixing rule according to Gen EOS model
    !!                                       6:      mixing rule according to PC SAFT model
    !!                                       7:
    !!                                       8:      mixing rule according to RKM model
    !!                                       9:      mixing rule according to COSTALD model
    !!path           chracter(255)           main directory of TREND
    !!unitdefintion  character(20)           choose molar/specific input and output
    !!handle         integer(c_intptr_t)     memory adress of gl type
    !!##################################################################
    !double precision function TREND_EOS (calctype, input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, errorflag, handle)!, unit)
    !!DEC$ ATTRIBUTES DLLEXPORT :: TREND_EOS
    !
    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !!Define input variables
    !character (12) :: input, unit
    !character(12) :: calctype
    !double precision :: prop1, prop2, PROP_EOS, wm_mix
    !character (30), dimension(30) :: fluids
    !double precision, dimension (30) :: moles
    !character (255) :: path
    !integer, dimension (30) :: EOS_indicator
    !integer :: MIX_indicator
    !character (1) :: propID
    !integer :: calctype_internal
    !integer :: ididit, errorflag
    !character (255) :: errorcodes
    !character (20) :: unitdefinition
    !logical :: fluid_present
    !integer :: i,nr_fld
    !integer :: IS_ERROR
    !
    !
    !!-------------------------------------------------------------------------------------------------------------------------------------------------
    !
    !!Reset internal calctype
    !calctype_internal = 0
    !errorflag = 0
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !i = sizeof(gl)
    !!default settings for all calculations:
    !gl%calc_ref = .false.  ! calculate reference state
    !gl%transport = 0       ! needed to check the range of validity for transport properties; 1 for transport properties, 0 for every other property
    !gl%VLE_needed = .false.
    !
    !if(trim(unitdefinition) == "specific") then
    !    gl%unitin=2
    !elseif(trim(unitdefinition) == "molar") then
    !    gl%unitin=1
    !else
    !    errorflag=-123456
    !    return
    !end if
    !
    !!Check for errors:
    !if(errorflag.ne.0) then
    !    TREND_EOS = errorflag
    !    return
    !end if
    !
    !!calctype to lower case
    !call uppertolower_char(calctype,len(calctype))
    !
    !!input to lower case
    !call uppertolower_char(input,len(input))
    !
    !!
    !
    !!-------------------------------------------------------------------------------------------------------------------------------------------------
    !
    !
    !!call uppertolower(gl, unit)
    !
    !!if(unitin == "mass")                         !define input for molspec routine: 1=molar input; 2=specific input
    !!    gl%unitin=2
    !!elseif(unitin == "mole")
    !!    gl%unitin=1
    !!end if
    !!
    !!!gl%molexp = 0
    !
    !
    !!***********************************************************************************************************************
    !if(trim(calctype) .eq. "t") then		![K]
    !    !Check if input includes T
    !    propID = "t"
    !    Call check_triv_input(gl,input, propID, ididit)
    !    if (ididit == 1) then
    !        if (prop1 .lt. 1.d-12) then     !temperature input less than 0 --> not defined
    !            TREND_EOS = -9911.d0
    !        else
    !            TREND_EOS = prop1
    !        end if
    !        return
    !    else if (ididit == 2) then
    !        if (prop2 .lt. 1.d-12) then     !temperature input less than 0 --> not defined
    !            TREND_EOS = -9911.d0
    !        else
    !            TREND_EOS = prop2
    !        end if
    !        return
    !    end if
    !    calctype_internal = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "d") then		![mol/m^3]
    !    !Check if input includes T
    !    propID = "d"
    !    Call check_triv_input(gl,input, propID, ididit)
    !    if (ididit == 1) then
    !        if (prop1 .lt. 1.d-12) then     !density input less than 0 --> not defined
    !            TREND_EOS = -9922.d0
    !        else
    !            TREND_EOS = prop1
    !        end if
    !        return
    !    else if (ididit == 2) then
    !        if (prop2 .lt. 1.d-12) then     !density input less than 0 --> not defined
    !            TREND_EOS = -9922.d0
    !        else
    !            TREND_EOS = prop2
    !        end if
    !        return
    !    end if
    !    calctype_internal = 2
    !    gl%converttype = 2
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "p") then		![MPa]
    !    !Check if input includes T
    !    propID = "p"
    !    Call check_triv_input(gl,input, propID, ididit)
    !    if (ididit == 1) then
    !        if (prop1 .lt. 1.d-12) then     !pressure input less than 0 --> not defined
    !            TREND_EOS = -9933.d0
    !        else
    !            TREND_EOS = prop1
    !        end if
    !        return
    !    else if (ididit == 2) then
    !        if (prop2 .lt. 1.d-12) then     !pressure input less than 0 --> not defined
    !            TREND_EOS = -9933.d0
    !        else
    !            TREND_EOS = prop2
    !        end if
    !        return
    !    end if
    !    calctype_internal = 3
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "u") then		![J/mol]
    !    gl%calc_ref = .true.  ! calculate reference state
    !    calctype_internal = 4
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "h") then		![J/mol]
    !    !Check if input includes T
    !    propID = "H"
    !    Call check_triv_input(gl,input, propID, ididit)
    !    if (ididit == 1) then
    !        TREND_EOS = prop1
    !        return
    !    else if (ididit == 2) then
    !        TREND_EOS = prop2
    !        return
    !    end if
    !    gl%calc_ref = .true.  ! calculate reference state
    !    calctype_internal = 5
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "s") then		![J/(mol K)]
    !    !Check if input includes T
    !    propID = "S"
    !    Call check_triv_input(gl,input, propID, ididit)
    !    if (ididit == 1) then
    !        TREND_EOS = prop1
    !        return
    !    else if (ididit == 2) then
    !        TREND_EOS = prop2
    !        return
    !    end if
    !    gl%calc_ref = .true.  ! calculate reference state
    !    calctype_internal = 6
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "g") then		![J/mol]
    !    gl%calc_ref = .true.  ! calculate reference state
    !    calctype_internal = 7
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a") then		![J/mol]
    !    gl%calc_ref = .true.  ! calculate reference state
    !    calctype_internal = 8
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "cp") then		![J/(mol K)]
    !    calctype_internal = 9
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "cv") then		![J/(mol K)]
    !    calctype_internal = 10
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "ws") then		![m/s]
    !    calctype_internal = 11
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "bvir") then		![m^3/mol]
    !    calctype_internal = 12
    !    gl%converttype = 3
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "cvir") then		![m^6/mol^2]
    !    calctype_internal = 13
    !    gl%converttype = 4
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "cp0") then      ![J/(mol K)]
    !    calctype_internal = 14
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "q") then		![mol/mol]
    !    calctype_internal = 15
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "z") then		![-]
    !    calctype_internal = 16
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "pnum") then		![MPa]
    !    calctype_internal = 17
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "wsnum") then	![m/s]
    !    calctype_internal = 18
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a00") then		![-]
    !    calctype_internal = 19
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a01") then		![-]
    !    calctype_internal = 20
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a02") then		![-]
    !    calctype_internal = 21
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a03") then		![-]
    !    calctype_internal = 22
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a10") then		![-]
    !    calctype_internal = 23
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a11") then		![-]
    !    calctype_internal = 24
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a12") then		![-]
    !    calctype_internal = 25
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a20") then		![-]
    !    calctype_internal = 26
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a21") then		![-]
    !    calctype_internal = 27
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "a30") then		![-]
    !    calctype_internal = 28
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "ur") then		![J/mol]
    !    calctype_internal = 29
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "hr") then		![J/mol]
    !    calctype_internal = 30
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "sr") then		![J/(mol K)]
    !    calctype_internal = 31
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "cvr") then		![J/(mol K)]
    !    calctype_internal = 32
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "cpotr") then		![J/mol]
    !    calctype_internal = 33
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "gruen") then		![-]
    !    calctype_internal = 34
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "pip") then		![-]
    !    calctype_internal = 35
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "de") then		![-]
    !    gl%transport = 3    !if transport=3, range of validity of the dielectricity EOS will be evaluated (if no "&" occurs in the input code)
    !    calctype_internal = 36
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "st") then		![mN/m]
    !    calctype_internal = 37
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "eta") then		![nu Pa /s]
    !    gl%transport = 1       !if transport=1, range of validity of the viscosity EOS will be evaluated (if no "&" occurs in the input code)
    !    calctype_internal = 38
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "tcx") then		![W/(m K)]
    !    gl%transport = 2       !if transport=2, range of validity of the thermal conductivity EOS will be evaluated (if no "&" occurs in the input code)
    !    calctype_internal = 39
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "vexp") then     ![1 / K]
    !    calctype_internal = 40
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "dvir") then		![m^9/mol^4]
    !    calctype_internal = 41
    !    gl%converttype = 5
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "riem") then         ![-]
    !    calctype_internal = 42
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "compt") then		![1/MPa]
    !    calctype_internal = 43
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "comps") then		![1/MPa]
    !    calctype_internal = 44
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "throt") then		![m^3/mol]
    !    calctype_internal = 45
    !    gl%converttype = 3
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "exps") then         ![-]
    !    calctype_internal = 46
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "expt") then         ![-]
    !    calctype_internal = 47
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "dudv") then         ![J / (mol m^3)]
    !    calctype_internal = 48
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "jtco") then			![]
    !    calctype_internal = 49
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "dpdt") then			![MPa/K]
    !    calctype_internal = 50
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "ve") then			![m^3/mol]
    !    calctype_internal = 51
    !    gl%converttype = 3
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "he") then			![J/mol]
    !    gl%calc_ref = .true.  ! calculate reference state
    !    calctype_internal = 52
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "ge") then			![J/mol]
    !    gl%calc_ref = .true.  ! calculate reference state
    !    calctype_internal = 53
    !    gl%converttype = 1
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "b12") then			![m^3/mol]
    !    calctype_internal = 54
    !    gl%converttype = 3
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "gammagd") then		![-]
    !    calctype_internal = 55
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "pr") then			![MPa]
    !    calctype_internal = 56
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "dbdt") then			![m^3/(mol T)]
    !    calctype_internal = 57
    !    gl%converttype = 3
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "dcdt") then			![m^6/(mol^2 T)]
    !    gl%converttype = 4
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "dpdd") then			![(Mpa m^3) / mol]
    !    calctype_internal = 59
    !    gl%converttype = 3
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "d2pdd2") then		![Mpa/(mol/m^3)^2]
    !    calctype_internal = 60
    !    gl%converttype = 4
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "d2pddt") then		![(Mpa m^3) / (mol K)]
    !    calctype_internal = 61
    !    gl%converttype = 3
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "d2pdt2") then		![Mpa / K^2]
    !    calctype_internal = 62
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "dddt") then         !
    !    calctype_internal = 63
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "dspin") then		![mol/m^3]
    !    calctype_internal = 64
    !    gl%converttype = 2
    !    !***********************************************************************************************************************
    !elseif(trim(calctype) .eq. "phase") then
    !    calctype_internal = 65
    !    !***********************************************************************************************************************
    !endif
    !
    !
    !if (calctype_internal .eq. 0) then
    !    TREND_EOS = -9935.d0
    !    return
    !end if
    !
    !
    !
    !!!  character(12) :: specific
    !!!
    !!!  specific = "specific"
    !!!
    !!!  if (specific == "specific")
    !!!      x_specific = moles
    !!!      call masstomolar(x_specific, fluids, x_moles, mw_array)
    !!!  end if
    !
    !!   if (gl%unitin == 2) then                                !muss noch verschoben werden.... die moles sind ja noch nicht gesetzt. irgendwo nach read_from_fluidfile
    !!
    !!   call wm_mix_calc(gl, wm_mix)
    !!
    !!   gl%molefractions = gl%molefractions * wm_mix
    !!
    !!   end if
    !!
    !
    !TREND_EOS = PROP_EOS(gl, calctype_internal, input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path)
    !
    !if(unitdefinition=="specific") then
    !    gl%unitin=2
    !end if
    !
    !
    !if ((IS_ERROR(int(TREND_EOS)) == 1).and.((TREND_EOS < -1000.d0) .and. (dabs(TREND_EOS - int(TREND_EOS)) < 1.d-14))) then
    !    TREND_EOS = TREND_EOS
    !else
    !    !molesv=0.d0
    !    call moles_incheck(gl,moles, errorflag)
    !    if(errorflag ==0) then
    !        gl%molfractions = moles
    !
    !        if(gl%unitin == 1)   then
    !
    !            TREND_EOS =  TREND_EOS
    !
    !        elseif((gl%unitin == 2) .and. (gl%converttype == 1)) then
    !
    !            call wm_mix_calc(gl, wm_mix)
    !
    !            TREND_EOS =  TREND_EOS / wm_mix !*1000
    !
    !        elseif((gl%unitin == 2) .and. (gl%converttype ==2))then
    !
    !            call wm_mix_calc(gl, wm_mix)
    !
    !            TREND_EOS =  TREND_EOS * wm_mix !* 1000
    !
    !        elseif((gl%unitin == 2) .and. (gl%converttype ==3))then
    !
    !            call wm_mix_calc(gl, wm_mix)
    !
    !            TREND_EOS =  TREND_EOS / wm_mix !* 1000
    !
    !        elseif((gl%unitin == 2) .and. (gl%converttype ==4))then
    !
    !            call wm_mix_calc(gl, wm_mix)
    !
    !            TREND_EOS =  TREND_EOS * (wm_mix)**(-2.d0) !* 1e6
    !
    !        elseif((gl%unitin == 2) .and. (gl%converttype ==5))then
    !
    !            call wm_mix_calc(gl, wm_mix)
    !
    !            TREND_EOS =  TREND_EOS * (wm_mix)**(-4.d0) !* 1e12
    !
    !        else
    !            TREND_EOS = TREND_EOS              !error in mass / mole specific handling
    !
    !        end if
    !
    !    else
    !        TREND_EOS = TREND_EOS
    !    end if
    !
    !
    !end if
    !end function TREND_EOS
    !
    !subroutine TREND_SPEC_EOS_STDCALL (fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, limits_text, limits_values, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "TREND_SPEC_EOS_STDCALL" :: TREND_SPEC_EOS_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: TREND_SPEC_EOS_STDCALL

    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !
    !
    !!define input variables
    !character(30), dimension(30), intent(in) :: fluids, limits_text
    !double precision, dimension(30), intent(in) :: moles
    !integer, dimension (30), intent(in) :: EOS_indicator
    !integer, intent(in) :: MIX_indicator
    !character (255), intent(in) :: path
    !character(20), intent(in) :: unitdefinition
    !double precision, dimension(30,31) :: limits_values
    !integer :: errorflag
    !
    !
    !call TREND_SPEC_EOS (fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, limits_text, limits_values, errorflag, handle)
    !
    !end subroutine
    !
    !!##################################################################
    !!TREND_SPEC_EOS function collects fluid/mixture specific information
    !!fluids         character(30),dim(30)   fluid/mixture calculations are conducted with
    !!moles          double,dim(30)          composition  calculations are conducted with
    !!EOS_indicator  integer,dim(30)         equation of state indicator calculations are conducted with
    !!                                       1: Helmholtz
    !!                                       2: SRK
    !!                                       3: PR
    !!                                       4: LKP
    !!                                       5: Generalized EOS (51: based on Alexandrov 52: based on Span/Wagner 53: based on Sun/Ely)
    !!                                       6: PCSAFT
    !!                                       7:
    !!                                       8: RKM
    !!                                       9: COSTALD
    !!MIX_indicator  integer                 indicator for mixing rule
    !!                                       1:      Helmholtz default,
    !!                                       110:    All components are mixed according to the linear mixing rule
    !!                                       111:    check if binary mix files are available, if not use linear mixing rule
    !!                                       120:    All components are mixed according to the Lorentz-Berthelot mixing rule.
    !!                                       121:    check if binary mix files are available, if not use Lorentz-Berthelot mixing rule.
    !!                                       2:      mixing rule according to SRK model
    !!                                       3:      mixing rule according to PR model
    !!                                       4:      mixing rule according to LKP model
    !!                                       5:      mixing rule according to Gen EOS model
    !!                                       6:      mixing rule according to PC SAFT model
    !!                                       7:
    !!                                       8:      mixing rule according to RKM model
    !!                                       9:      mixing rule according to COSTALD model
    !!path           chracter(255)           main directory of TREND
    !!unitdefintion  character(20)           choose molar/specific input and output
    !!limits_text    character(30),dim(30)   desciption of the fluid property or limit:
    !!                                       1:  molar mass
    !!                                       2:  temperature at triple point
    !!                                       3:  pressure at triple point
    !!                                       4:  temperature at crtical point
    !!                                       5:  pressure at critical point
    !!                                       6:  density at critical point
    !!                                       7:  acentric factor
    !!                                       8:  minimum specified temperature of eos
    !!                                       9:  maximum specified temperature of eos
    !!                                       10: maximum specified pressure of eos
    !!                                       11: maximum specified density of eos
    !!                                       12: cas-nr -> value: -1000
    !!limits_values  double,dim(30)          value of the corresponding property or limi
    !!##################################################################
    !subroutine TREND_SPEC_EOS (fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, limits_text, limits_values, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT :: TREND_SPEC_EOS
    !

    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !
    !integer:: errorflag,e
    !character(255):: path
    !character(20) :: unitdefinition
    !integer :: nf, ne, nm, i
    !character (30), dimension (30) :: fluids
    !double precision, dimension (30) :: moles
    !character (30), dimension (30) :: limits_text
    !integer, dimension (30) :: EOS_indicator
    !integer :: MIX_indicator
    !double precision, dimension(30,31) :: limits_values
    !
    !character(12):: input
    !double precision:: prop1,prop2
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !input = ''
    !input = 'TP'
    !prop1 = 300.d0
    !prop2 = 1.d0
    !errorflag = 0
    !if(trim(unitdefinition) == "specific") then
    !    gl%unitin=2
    !elseif(trim(unitdefinition) == "molar") then
    !    gl%unitin=1
    !else
    !    limits_values=-123456
    !    return
    !end if
    !!                 1     2       3       4       5      6              7           8       9    10    11     12
    !limits_text = (/'MW','Ttrip','ptrip','Tcrit','pcrit','Dcrit','AF','Tmin','Tmax','pmax','Dmax','','','','','','','','','','','','','','','','','','',''/)
    !limits_values = 0.d0
    !
    !!dummy setup call to guarantte memory allocation
    !call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator, errorflag)
    !
    !
    !if (errorflag == 0) then
    !
    !    if (gl%ncomp == 1) then
    !        limits_values(1,1) = gl%wm(1)
    !    else
    !        call wm_mix_calc(gl,limits_values(1,gl%ncomp+1))
    !    end if
    !    do i = 1, gl%ncomp
    !        limits_values(1,i) = gl%wm(i)
    !        limits_values(2,i) = gl%ttp(i)
    !        limits_values(3,i) = gl%ptp(i)
    !        limits_values(4,i) = gl%tc(i)
    !        limits_values(5,i) = gl%pc(i)
    !
    !        if (gl%unitin == 2) then
    !            limits_values(6,i) = gl%rhoc(i) * gl%wm(i)
    !        else
    !            limits_values(6,i) = gl%rhoc(i)
    !        endif
    !        if (gl%substcasnr(i) == trim('999-999-999')) limits_values(6,i) = gl%rhoc(i) / gl%Factor
    !
    !        limits_values(7,i) = gl%accen(i)
    !        limits_values(8,i) = gl%tminfluid(i)
    !        limits_values(9,i) = gl%tmaxfluid(i)
    !        limits_values(10,i) = gl%pmaxfluid(i)
    !        if (gl%unitin == 2) then
    !            limits_values(11,i) = gl%rhomaxfluid(i) * gl%wm(i)
    !        else
    !            limits_values(11,i) = gl%rhomaxfluid(i)
    !        endif
    !        if (gl%substcasnr(1) == trim('999-999-999')) limits_values(11,i) = gl%rhomaxfluid(i) / gl%Factor
    !    end do
    !
    !else
    !    limits_values = errorflag
    !end if
    !
    !
    !end subroutine TREND_SPEC_EOS



    !!****************************************************************************************
    !! Fluidname to CASNR
    !!****************************************************************************************
    !subroutine CAS_FROM_FLUID_STDCALL (fluids, EOS_indicator, MIX_indicator, path, CASNR, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "CAS_FROM_FLUID_STDCALL" :: CAS_FROM_FLUID_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: CAS_FROM_FLUID_STDCALL
    !

    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !
    !
    !!define input variables
    !character (30), dimension(30), intent(inout) :: fluids
    !integer, dimension(30), intent(inout) :: EOS_indicator
    !integer, intent(inout) :: MIX_indicator
    !character (255), intent(inout) :: path
    !!define output variables
    !character (12), intent (out) :: CASNR
    !integer :: errorflag
    !
    !
    !
    !call CAS_FROM_FLUID (fluids, EOS_indicator, MIX_indicator, path, CASNR, errorflag, handle)
    !
    !end subroutine
    !
    !subroutine CAS_FROM_FLUID (fluids, EOS_indicator, MIX_indicator, path, CASNR, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT :: CAS_FROM_FLUID
    !


    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !
    !!define input variables
    !character (30), dimension(30), intent(inout) :: fluids
    !integer, dimension(30), intent(inout) :: EOS_indicator
    !integer, intent(inout) :: MIX_indicator
    !character (255), intent(inout) :: path
    !!define output variables
    !character (12), intent (out) :: CASNR
    !
    !character (12) :: input
    !double precision :: prop1, prop2
    !double precision, dimension(30) :: moles
    !
    !integer :: errorflag
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !gl%calc_ref = .false.  ! no need to calculate reference state
    !gl%transport = 0       ! needed to check the range of validity for transport properties; 0 for every other property
    !
    !input = 'tp'
    !prop1 = 300.d0
    !prop2 = 1.d0
    !moles = 0.d0
    !moles(1:count(fluids /= '')) = 1.d0 / count(fluids /= '')
    !
    !call setup (gl,input, prop1, prop2, fluids, moles, path,EOS_indicator, mix_indicator, errorflag)
    !if (errorflag /= 0) then
    !    write(CASNR,*)errorflag
    !else
    !    CASNR = gl%substcasnr(1)
    !end if
    !
    !
    !end subroutine CAS_FROM_FLUID



    !!------------------------------------------------
    !!Function for composition convertion
    !!------------------------------------------------
    !subroutine COMP_STDCALL(fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, composition, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "COMP_STDCALL" :: COMP_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: COMP_STDCALL

    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !
    !
    !!definie input:
    !character(30), dimension(30), intent(inout) :: fluids
    !double precision, dimension(30), intent(inout) :: moles
    !integer, dimension(30), intent(inout) :: EOS_indicator!, unitdefinition
    !integer, intent(inout) :: MIX_indicator
    !character (255), intent(inout) :: path
    !character(20) :: unitdefinition
    !double precision, dimension (30), intent (out) :: composition
    !integer :: errorflag
    !
    !call COMP (fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, composition, errorflag, handle)
    !
    !end subroutine
    !
    !subroutine COMP(fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, composition, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT::COMP
    !
    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !
    !
    !!definie input:
    !character(30), dimension(30), intent(inout) :: fluids
    !double precision, dimension(30), intent(inout) :: moles
    !integer, dimension(30), intent(inout) :: EOS_indicator!, unitdefinition
    !integer, intent(inout) :: MIX_indicator
    !character (255), intent(inout) :: path
    !character(20) :: unitdefinition
    !character(12) :: input
    !double precision :: prop1, prop2, wm_mix
    !integer :: errorflag, converttype            !dummy for setup
    !double precision, dimension (30), intent(out) :: composition
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !prop1 = 300.d0
    !prop2 = 1.d0
    !input = trim('tp')
    !
    !!gl%unitin = 1
    !
    !call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator,errorflag)
    !
    !unitdefinition = trim(unitdefinition)
    !
    !call uppertolower_char(unitdefinition,len(unitdefinition))
    !
    !if(unitdefinition == "molartospecific") then
    !    converttype = 2
    !elseif(unitdefinition == "specifictomolar") then
    !    converttype = 1
    !else
    !    errorflag = -123456
    !    return
    !end if
    !
    !call wm_mix_calc(gl, wm_mix)
    !
    !call convert_fractions(gl, converttype, wm_mix, moles)
    !
    !
    !composition = moles
    !
    !end subroutine COMP



    !!*******************************************************
    !! Triple point of pure substances
    !!*******************************************************
    !! Andreas, July 2014
    !!
    !subroutine TRIPLE_POINT_STDCALL (Temp, press, fluids, path, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "TRIPLE_POINT_STDCALL" :: TRIPLE_POINT_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: TRIPLE_POINT_STDCALL
    !

    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !
    !
    !!define input variables
    !character (30), dimension(30), intent(in) :: fluids
    !character (255), intent(in) :: path
    !!define output variables
    !double precision, intent (out) :: Temp, press
    !integer :: errorflag
    !
    !
    !
    !call TRIPLE_POINT (Temp, press, fluids, path, errorflag, handle)
    !
    !end subroutine

    !subroutine TRIPLE_POINT(Temp, press, fluids, path, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT :: TRIPLE_POINT
    !




    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !
    !character (12) :: input
    !double precision, intent (out) :: Temp, press
    !character (30), dimension(30) :: fluids
    !double precision, dimension(30) :: moles
    !character (255) :: path
    !integer, dimension(30) :: EOS_indicator
    !integer:: MIX_indicator
    !double precision, dimension(30) ::  x_vap, x_liq, x_sol, x_hyd, x_known
    !double precision, dimension (3) :: phasefrac, rho
    !
    !double precision :: rhovap_est, rholiq_est
    !
    !integer:: errorflag, iphase, iter, iFlash
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !gl%transport = 0       ! needed to check the range of validity for transport properties; 0 for every other property
    !
    !!Call setup to read fluid files and check inputs
    !Temp = 300.D0
    !press = 1.D0
    !input = "tp+"
    !moles = 0.d0
    !moles(1) = 1.d0
    !EOS_indicator = 0
    !MIX_indicator = 0
    !EOS_indicator(1) = 1
    !MIX_indicator = 1
    !
    !call setup (gl,input, Temp, press, fluids, moles, path, EOS_indicator, MIX_indicator, errorflag)
    !
    !! irregular input, or problems with fluid file or reference state
    !if (errorflag /= 0) then
    !    Temp = errorflag
    !    press = errorflag
    !    return
    !end if
    !
    !if (gl%ncomp > 1) then
    !    errorflag = -9902
    !    Temp = errorflag
    !    press = errorflag
    !    return
    !end if
    !
    !!Composition of all phases is known (pure substance)
    !x_known(1) = 1.D0
    !x_vap(1) = 1.D0
    !x_liq(1) = 1.D0
    !x_hyd = 0.D0
    !x_sol(1) = 1.D0
    !
    !rhovap_est = 0.D0
    !rholiq_est = 0.D0
    !
    !iFlash = 8
    !iphase = 2
    !gl%solidtype(2) = 0
    !gl%solid_pos = 1
    !
    !if (gl%Fluidlist_hydrate(1) == "co2") then
    !    !CO2
    !    Temp = 216.5D0
    !    press = 0.517D0
    !    gl%solidtype(1) = 2
    !elseif (gl%Fluidlist_hydrate(1) == "water") then
    !    !H2O
    !    Temp = 273.D0
    !    press = 611.D-6
    !    gl%solidtype(1) = 1
    !else
    !    !error, fluid not yet implemented
    !    errorflag = -12902
    !    Temp = errorflag
    !    press = errorflag
    !    return
    !end if
    !
    !call ptflash_solid_NC_3P(gl,press, Temp, x_known, rho, x_vap, x_liq, x_sol, x_hyd, rhovap_est, &
    !    & rholiq_est, Phasefrac, iFlash, iphase, iter, errorflag)
    !
    !if (errorflag /= 0) then
    !    Temp = errorflag
    !    press = errorflag
    !    return
    !end if
    !
    !end subroutine TRIPLE_POINT





    !!__________________________________________________________________________________________________________________________________________________________________________________________________________
    !
    !!***************************************************************************************
    !!Subroutine for the estimation of uncertainties of the equation of state used.
    !!At the moment This routine works only for PURE SUBSTANCES and with input T and p !!
    !!***************************************************************************************
    !subroutine ALL_UNCTY_STDCALL (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, UNC_D, UNC_U, UNC_H, UNC_S, UNC_CP, UNC_CV, UNC_WS, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "ALL_UNCTY_STDCALL" :: ALL_UNCTY_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: ALL_UNCTY_STDCALL

    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !
    !
    !!define input variables
    !character (12), intent(in) :: input
    !double precision, intent(in) :: prop1, prop2
    !character (30), dimension(30), intent(in) :: fluids
    !double precision, dimension(30), intent(in) :: moles
    !integer, dimension(30), intent(in) :: EOS_indicator
    !character (255), intent(in) :: path
    !integer:: MIX_indicator
    !
    !!define output variables
    !double precision, intent(out):: UNC_D, UNC_U, UNC_H, UNC_S, UNC_CP, UNC_CV, UNC_WS, p
    !integer, intent(out) :: errorflag
    !
    !
    !
    !call ALL_UNCTY (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, UNC_D, UNC_U, UNC_H, UNC_S, UNC_CP, UNC_CV, UNC_WS, errorflag, handle)
    !
    !end subroutine

    !SUBROUTINE ALL_UNCTY (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, UNC_D, UNC_U, UNC_H, UNC_S, UNC_CP, UNC_CV, UNC_WS, errorflag, handle)
    !
    !
    !!   NOT AVAILABLE IN TREND 1.1
    !!  !DEC$ ATTRIBUTES DLLEXPORT :: ALL_UNCTY
    !




    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !
    !double precision, intent(out):: UNC_D, UNC_U, UNC_H, UNC_S, UNC_CP, UNC_CV, UNC_WS
    !
    !double precision :: prop1, prop2
    !character (30), dimension(30) :: fluids
    !character (255) :: path
    !integer, dimension(30) :: EOS_indicator
    !integer:: MIX_indicator
    !character (12) :: input
    !
    !double precision, dimension(30) :: moles
    !
    !integer, intent(out) :: errorflag
    !integer :: ILIMITS
    !
    !double precision:: x_Phase(30,5)    ! vector conaining the compositions of all phases
    !integer:: nrofphases                ! at this point: only 1 or 2 implemented
    !!Indicate which phases are present
    !integer:: phasetype(5)              !phasetype contains the phase indicator number
    !!define variables for calculating
    !double precision :: t, d, p
    !double precision :: dvap
    !double precision :: dliq
    !double precision :: vapfrac
    !
    !double precision :: h, s
    !double precision:: rho(5)           ! vector containing the densities of all phases
    !double precision :: x_liq(30)
    !double precision :: x_vap(30)
    !
    !
    !!Uncertainty variables
    !integer:: uncty_call_prop, iter, nrsubst, IFlash, phase
    !LOGICAL :: twophase
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !gl%transport = 0       ! needed to check the range of validity for transport properties; 0 for every other property
    !
    !
    !ILIMITS = 0
    !errorflag = 0
    !UNC_D = 0.D0
    !UNC_U = 0.D0
    !UNC_H = 0.D0
    !UNC_S = 0.D0
    !UNC_CP = 0.D0
    !UNC_CV = 0.D0
    !UNC_WS = 0.D0
    !
    !twophase = .false.
    !vapfrac = 0.D0 !Not needed since input is only T and p
    !
    !!trim input
    !call uppertolower_char(input,len(input))
    !
    !!if (EOS_ind == 2) then
    !!    UNC_D = -2908
    !!    UNC_U = -2908
    !!    UNC_H = -2908
    !!    UNC_S = -2908
    !!    UNC_CP = -2908
    !!    UNC_CV = -2908
    !!    UNC_WS = -2908
    !!    return
    !!end if
    !if (.not.allocated(gl%uncty)) allocate(gl%uncty)
    !!!call setup - >  initializing, reading fluid files
    !call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator, errorflag)
    !
    !! irregular input, or problems with fluid file or reference state
    !if (errorflag /= 0) then
    !    UNC_D = errorflag
    !    UNC_U = errorflag
    !    UNC_H = errorflag
    !    UNC_S = errorflag
    !    UNC_CP = errorflag
    !    UNC_CV = errorflag
    !    UNC_WS = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !
    !!Catch wrong inputs
    !if (gl%ncomp > 1) then
    !    errorflag = -7000      !No mixtures allowed
    !    UNC_D = errorflag
    !    UNC_U = errorflag
    !    UNC_H = errorflag
    !    UNC_S = errorflag
    !    UNC_CP = errorflag
    !    UNC_CV = errorflag
    !    UNC_WS = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !if (gl%Eq_type(1) > 1) then
    !    errorflag = -7002
    !    UNC_D = errorflag
    !    UNC_U = errorflag
    !    UNC_H = errorflag
    !    UNC_S = errorflag
    !    UNC_CP = errorflag
    !    UNC_CV = errorflag
    !    UNC_WS = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!calculate CALC
    !!______________________________________________________________________________
    !if (input == 'td') then
    !    t = prop1
    !    d = prop2
    !
    !    !check whether pure fluid or mixture
    !    if (gl%ncomp == 1) then !pure fluid
    !        nrsubst = 1
    !        call PhaseDet_pure(gl,p, T, d, dvap, dliq, phase, vapfrac, nrofphases, errorflag)
    !    else  ! mixture
    !        nrsubst = 0
    !        call PhaseDet_td(gl,p, t, moles, rho, x_Phase, phasetype, vapfrac, d, nrofphases, errorflag)
    !    end if
    !    !______________________________________________________________________________
    !    !calculate TP
    !else if (input == 'tp') then
    !    t = prop1
    !    p = prop2
    !    !check whether pure fluid or mixture
    !    if (gl%ncomp == 1) then !pure fluid
    !        nrsubst = 1
    !        call PhaseDet_pure_tp(gl,p, t, d, dvap, dliq, phase, vapfrac, nrofphases, errorflag)
    !    else  !mixture
    !        nrsubst = 0
    !        call PhaseDet(gl,p, t, moles, rho, x_Phase, phasetype, vapfrac, nrofphases, errorflag)
    !    end if
    !
    !    !______________________________________________________________________________
    !    !calculate PH
    !else if (input == 'ph') then
    !    p = prop1
    !    h = prop2
    !    !check whether pure fluid or mixture
    !    if (gl%ncomp == 1) then !pure fluid
    !        nrsubst = 1
    !        call PhaseDet_ph_pure(gl,p, t, d, dvap, dliq, phase, vapfrac, h, nrofphases, errorflag)
    !    else    !mixture
    !        nrsubst = 0
    !        call PhaseDet_ph(gl,p, t, moles, rho, x_phase, phasetype, vapfrac, h, nrofphases, errorflag)
    !    end if
    !    !______________________________________________________________________________
    !    !calculate PS
    !else if (input == 'ps') then
    !    p = prop1
    !    s = prop2
    !    !check whether pure fluid or mixture
    !    if (gl%ncomp == 1) then !pure fluid
    !        nrsubst = 1
    !        call PhaseDet_ps_pure(gl,p, t, d, dvap, dliq, phase, vapfrac, s, nrofphases, errorflag)
    !    else    !mixture
    !        nrsubst = 0
    !        call PhaseDet_ps(gl,p, t, moles, rho, x_phase, phasetype, vapfrac, s, nrofphases, errorflag)
    !    end if
    !    !______________________________________________________________________________
    !else if (input == 'tliq') then
    !    t = prop1
    !    ! check whether pure fluid or mixture
    !    if (gl%ncomp == 1) then !pure fluid
    !        nrsubst = 1
    !        call Flash_Pure_PhaseBoundary(gl,p, t, dvap, dliq, 1, errorflag, iter, nrsubst)
    !    else    !mixture
    !        nrsubst = 0
    !        iFlash = 3
    !        p = 0.d0
    !        call Flash_PhaseBoundary(gl,p, t, moles, x_vap, x_liq, 0.d0, 0.d0, vapfrac, iFlash,&
    !            &  0, errorflag, iter)
    !    end if
    !    vapfrac = 0.d0
    !    nrofphases = 2
    !    !______________________________________________________________________________
    !else if (input == 'tvap') then
    !    t = prop1
    !    ! check whether pure fluid or mixture
    !    if (gl%ncomp == 1) then !pure fluid
    !        nrsubst = 1
    !        call Flash_Pure_PhaseBoundary(gl,p, t, dvap, dliq, 1, errorflag, iter, nrsubst)
    !    else    !mixture
    !        nrsubst = 0
    !        iFlash = 4
    !        p = 0.d0
    !        call Flash_PhaseBoundary(gl,p, t, moles, x_vap, x_liq, 0.d0, 0.d0, vapfrac, iFlash,&
    !            &  0, errorflag, iter)
    !    end if
    !    vapfrac = 1.d0
    !    nrofphases = 2
    !    !______________________________________________________________________________
    !else if (input == 'pliq') then
    !    p = prop1
    !    ! check whether pure fluid or mixture
    !    if (gl%ncomp == 1) then !pure fluid
    !        nrsubst = 1
    !        call Flash_Pure_PhaseBoundary(gl,p, t, dvap, dliq, 2, errorflag, iter, nrsubst)
    !    else    !mixture
    !        nrsubst = 0
    !        iFlash = 1
    !        t = 0.d0
    !        call Flash_PhaseBoundary(gl,p, t, moles, x_vap, x_liq, 0.d0, 0.d0, vapfrac, iFlash,&
    !            &  0, errorflag, iter)
    !    end if
    !    vapfrac = 0.d0
    !    nrofphases = 2
    !    !______________________________________________________________________________
    !else if (input == 'pvap') then
    !    p = prop1
    !    ! check whether pure fluid or mixture
    !    if (gl%ncomp == 1) then !pure fluid
    !        nrsubst = 1
    !        call Flash_Pure_PhaseBoundary(gl,p, t, dvap, dliq, 2, errorflag, iter, nrsubst)
    !
    !    else    !mixture
    !        nrsubst = 0
    !        iFlash = 2
    !        t = 0.d0
    !        call Flash_PhaseBoundary(gl,p, t, moles, x_vap, x_liq, 0.d0, 0.d0, vapfrac, iFlash,&
    !            &  0, errorflag, iter)
    !    end if
    !    vapfrac = 1.d0
    !    nrofphases = 2
    !end if
    !
    !!errorflag =! 0
    !if (errorflag /= 0) then
    !    UNC_D = errorflag
    !    UNC_U = errorflag
    !    UNC_H = errorflag
    !    UNC_S = errorflag
    !    UNC_CP = errorflag
    !    UNC_CV = errorflag
    !    UNC_WS = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !
    !!---end of Phasedet---
    !
    !
    !!errorflag = 0
    !if (nrofphases == 2) then    !2phase                                                                       !2phase
    !    if (gl%ncomp == 1) then                                                                                !2phase, pure fluid
    !        errorflag = -7003
    !        UNC_D = errorflag
    !        UNC_U = errorflag
    !        UNC_H = errorflag
    !        UNC_S = errorflag
    !        UNC_CP = errorflag
    !        UNC_CV = errorflag
    !        UNC_WS = errorflag
    !        input=gl%inptorig
    !        return
    !    else
    !        !2phase, mixture
    !    end if
    !
    !else        !1phase
    !    !Code kommt unten
    !
    !end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !
    !!Calculate the uncertainties of the following properties:
    !!H  : 1
    !!U  : 2
    !!CV : 3
    !!S  : 4
    !!CP : 5
    !!WS : 6
    !!t = prop1
    !!p = prop2
    !!d = rhomix_calc(gl,t, p, 0.D0, 0, 1)
    !!-----------------------------------------------------
    !!TEST OF UNCERTAINTY ROUTINES!!
    !!Get density and enthalpy uncertainty (The uncertainty in density is always given by the authors and updated no matter which uncertainty is calculated)
    !Uncty_call_prop = 1
    !call callUncty(gl,uncty_call_prop, t, d, twophase, vapfrac, errorflag, path)
    !if (errorflag == 0) then
    !    UNC_H = gl%uncty%uncty_h
    !    UNC_D = gl%uncty%uncty_d
    !else
    !    UNC_D = errorflag
    !    UNC_U = errorflag
    !    UNC_H = errorflag
    !    UNC_S = errorflag
    !    UNC_CP = errorflag
    !    UNC_CV = errorflag
    !    UNC_WS = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !!Get internal energy uncertainty
    !Uncty_call_prop = 2
    !call callUncty(gl,uncty_call_prop, t, d, twophase, vapfrac, errorflag, path)
    !if (errorflag == 0) then
    !    UNC_U = gl%uncty%uncty_u
    !else
    !    UNC_D = errorflag
    !    UNC_U = errorflag
    !    UNC_H = errorflag
    !    UNC_S = errorflag
    !    UNC_CP = errorflag
    !    UNC_CV = errorflag
    !    UNC_WS = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !!Get isochoric heat capacity uncertainty
    !Uncty_call_prop = 3
    !call callUncty(gl,uncty_call_prop, t, d, twophase, vapfrac, errorflag, path)
    !if (errorflag == 0) then
    !    UNC_CV = gl%uncty%uncty_cv
    !else
    !    UNC_D = errorflag
    !    UNC_U = errorflag
    !    UNC_H = errorflag
    !    UNC_S = errorflag
    !    UNC_CP = errorflag
    !    UNC_CV = errorflag
    !    UNC_WS = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !!Get entropy uncertainty
    !Uncty_call_prop = 4
    !call callUncty(gl,uncty_call_prop, t, d, twophase, vapfrac, errorflag, path)
    !if (errorflag == 0) then
    !    UNC_S = gl%uncty%uncty_s
    !else
    !    UNC_D = errorflag
    !    UNC_U = errorflag
    !    UNC_H = errorflag
    !    UNC_S = errorflag
    !    UNC_CP = errorflag
    !    UNC_CV = errorflag
    !    UNC_WS = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !!Get isobaric heat capacity uncertainty
    !Uncty_call_prop = 5
    !call callUncty(gl,uncty_call_prop, t, d, twophase, vapfrac, errorflag, path)
    !if (errorflag == 0) then
    !    UNC_CP = gl%uncty%uncty_cp
    !else
    !    UNC_D = errorflag
    !    UNC_U = errorflag
    !    UNC_H = errorflag
    !    UNC_S = errorflag
    !    UNC_CP = errorflag
    !    UNC_CV = errorflag
    !    UNC_WS = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !!Get sound speed uncertainty
    !Uncty_call_prop = 6
    !call callUncty(gl,uncty_call_prop, t, d, twophase, vapfrac, errorflag, path)
    !if (errorflag == 0) then
    !    UNC_WS = gl%uncty%uncty_w
    !else
    !    UNC_D = errorflag
    !    UNC_U = errorflag
    !    UNC_H = errorflag
    !    UNC_S = errorflag
    !    UNC_CP = errorflag
    !    UNC_CV = errorflag
    !    UNC_WS = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !!-----------------------------------------------------
    !
    !
    !!check limits
    !call check_limits (gl,input, t, p, d, dliq, dvap, nrofphases, x_phase, phasetype, ILIMITS)
    !
    !if (ILIMITS /= 0) then
    !    UNC_D = ILIMITS
    !    UNC_U = ILIMITS
    !    UNC_H = ILIMITS
    !    UNC_S = ILIMITS
    !    UNC_CP = ILIMITS
    !    UNC_CV = ILIMITS
    !    UNC_WS = ILIMITS
    !    errorflag = ILIMITS
    !    input=gl%inptorig
    !    return
    !end if
    !
    !! für Krohne Umrechnung in spezifische Einheiten
    !!call wm_mix_calc(wm_mix)
    !!UNC_D = UNC_D*wm_mix
    !!UNC_U = UNC_U/wm_mix
    !!UNC_H = UNC_H/wm_mix
    !!UNC_S = UNC_S/wm_mix
    !!UNC_CP = UNC_CP/wm_mix
    !!UNC_CV = UNC_CV/wm_mix
    !
    !input=gl%inptorig
    !
    !
    !end subroutine ALL_UNCTY




    !!*******************************************************
    !! subroutine to get literature reference for fluid or mixing rule
    !!*******************************************************
    !subroutine LITERATURE_REF_STDCALL (fluids, EOS_indicator, MIX_indicator, prop, path, reference_out, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "LITERATURE_REF_STDCALL" :: LITERATURE_REF_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: LITERATURE_REF_STDCALL

    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !
    !!define input variables
    !character (30), dimension(30) :: fluids
    !integer, dimension(30) :: EOS_indicator
    !integer, intent(inout):: MIX_indicator
    !character (10), intent(in) :: prop
    !character (255), intent(in) :: path
    !!define output variables
    !character (1000), intent(out) :: reference_out
    !integer, intent(out) :: errorflag
    !
    !
    !
    !call LITERATURE_REF (fluids, eos_indicator, mix_indicator, prop, path, reference_out, errorflag, handle)
    !
    !end subroutine

    !subroutine LITERATURE_REF (fluids, EOS_indicator, MIX_indicator, prop, path, reference_out, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT :: LITERATURE_REF
    !


    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !character (30), dimension(30) :: fluids
    !character (10) :: prop
    !character (255) :: path
    !character (1000) :: reference_out
    !
    !integer, dimension(30) :: EOS_indicator
    !integer, intent(inout):: MIX_indicator
    !character (12) :: input
    !double precision :: prop1, prop2
    !double precision, dimension (30) :: moles
    !
    !integer :: i,j, errorflag, nf
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !if (.not. allocated(gl%litref)) allocate(gl%litref)
    !
    !errorflag = 0
    !input = 'tp'
    !prop1 = 300d0
    !prop2 = 0.1d0
    !nf = 0
    !do i = 1, 2
    !    if (EOS_indicator(i) /= 0) then
    !        nf = nf + 1
    !    endif
    !enddo
    !moles = 0.d0
    !if (nf == 1) then
    !    moles(1) = 1.d0
    !    !EOS_indicator = '1;1'
    !else
    !    moles(1:2) = 0.5d0
    !    !EOS_indicator = '1;1;1'
    !endif
    !
    !call setup (gl, input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator, errorflag)
    !if (errorflag /= 0) then
    !    reference_out = 'Setup error'
    !    return
    !endif
    !call uppertolower_char(prop, len(prop))
    !
    !select case (trim(prop))
    !case ('p','d','t','u','h','s','g','a','cp','cv','ws','bvir','cvir','z')
    !    if (gl%ncomp == 1) then
    !        reference_out = gl%litref%lit_ref_res(1)
    !    elseif (gl%ncomp == 2) then
    !        reference_out = gl%litref%lit_ref_mix
    !    endif
    !case ('de')
    !    reference_out = gl%litref%lit_ref_de(1)
    !case ('tcx')
    !    reference_out = gl%litref%lit_ref_tcx(1)
    !case ('stn')
    !    reference_out = gl%litref%lit_ref_stn(1)
    !case ('eta')
    !    reference_out = gl%litref%lit_ref_eta(1)
    !case ('mlt')
    !    reference_out = gl%litref%lit_ref_mlt(1)
    !case ('sbl')
    !    reference_out = gl%litref%lit_ref_sbl(1)
    !    case default
    !    reference_out = 'Invalid property input'
    !end select
    !
    !end subroutine LITERATURE_REF

    !!*******************************************************************************
    !! "Expert" user functions for the calculation of two and three phase equilibria
    !! Andreas June 2013
    !!*******************************************************************************
    !subroutine FLASH_EXPERT_STDCALL (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, xliq_in, xvap_in, rholiq_est, rhovap_est, iPhase_try, xliq, xvap, propliq, propvap, vapfrac, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "FLASH_EXPERT_STDCALL" :: FLASH_EXPERT_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: FLASH_EXPERT_STDCALL
    !
    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !
    !!define input variables
    !character (12) :: input
    !character(20) :: unitdefinition
    !double precision :: prop1, prop2
    !character (30), dimension(30) :: fluids
    !double precision, dimension(30) :: moles
    !integer, dimension(30) :: EOS_indicator
    !integer :: MIX_indicator
    !character (255) :: path
    !double precision, dimension(30):: xliq_in, xvap_in
    !double precision:: rhovap_est, rholiq_est
    !integer :: iphase_try
    !!define output variables
    !double precision, dimension(30), intent(out) :: xliq
    !double precision, dimension(30), intent(out) :: xvap
    !double precision, dimension(30), intent(out) :: propliq
    !double precision, dimension(30), intent(out) :: propvap
    !double precision, intent(out) :: vapfrac
    !
    !integer :: errorflag
    !
    !
    !
    !call FLASH_EXPERT (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, xliq_in, xvap_in, rholiq_est, rhovap_est, iPhase_try, xliq, xvap, propliq, propvap, vapfrac, errorflag, handle)
    !
    !end subroutine

    !SUBROUTINE FLASH_EXPERT (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, &
    !    & xliq_in, xvap_in, rholiq_est, rhovap_est, iPhase_try, xliq, xvap, propliq, propvap, vapfrac, errorflag, handle)
    !! EXPERT routine for the calculation of a two phase flash
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! ALL INITIAL ESTIMATES NECESSARY HAVE TO BE GIVEN BY THE USER!!!
    !! THE RESULT IS NOT EVALUATED FOR STABILITY, THUS METASTABLE STATES MIGHT BE CALCULATES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! 2-phase flash calculation with t-p, p-h, p-s and t-d as input parameters possible
    !! This routine returns the compositions and all thermodynamic properties of the coexisting phases
    !!---------------------------------------------------
    !! Input parameters:
    !!   prop1, prop2                - according to the chosen flash type: T,p / p,s / p,h
    !!   xliq, xvap                  - initial estimates for the phase composition
    !!   rhovap_est, rholiq_est      - initial estimates for the liquid and vapor densities (optional)
    !!   EOS_indicator               - specifies which equation of state shall be used
    !!   path                        - path to the fluid files
    !!   iPhase_try                  - 0 - >  let the density solver choose the likely density for both phases,  vap: iphase = 0, liq : iphase = 0
    !!                               - 1 - >  LLE assumed, vap: iphase = 1, liq : iphase = 0
    !!                               - 2 - >  VLE assumed, vap: iphase = 2, liq : iphase = 0
    !!                               - 3 - >  LLE assumed, vap: iphase = 0, liq : iphase = 1
    !!                               - 4 - >  LLE assumed, vap: iphase = 1, liq : iphase = 1
    !!                               - 5 - >  VLE assumed, vap: iphase = 2, liq : iphase = 1
    !!                               - 6 - >  VLE assumed, vap: iphase = 1, liq : iphase = 2
    !! Output parameters:
    !!   xliq, xvap          composition vectors, length 30
    !!   propliq, propvap    property vectors of the two phases, properties are positioned in the following order:
    !!                       P, T, D, U, H, S, G, A, CP, CV, WS - lengh 30
    !!   vapfrac             vapor fraction: returns -2 for superheated vapor and -1 for subcooled liquid
    !!DEC$ ATTRIBUTES DLLEXPORT :: FLASH_EXPERT
    !
    !!No specific input, as only for expert users
    !
    !



    !
    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !
    !character (12) :: input
    !character(20) :: unitdefinition
    !double precision :: prop1, prop2
    !character (30), dimension(30) :: fluids
    !double precision, dimension(30) :: moles
    !character (255) :: path
    !integer, dimension(30) :: EOS_indicator
    !integer :: MIX_indicator
    !double precision, dimension(30):: xliq_in, xvap_in
    !
    !double precision, dimension(30), intent(out) :: xliq
    !double precision, dimension(30), intent(out) :: xvap
    !double precision, dimension(30), intent(out) :: propliq
    !double precision, dimension(30), intent(out) :: propvap
    !double precision, intent(out) :: vapfrac
    !
    !double precision:: t, p, h, s, dliq, dvap
    !double precision:: rho(5)
    !double precision:: u_calc, h_calc, s_calc, g_calc, A_CALC, cp_calc, cv_calc, ws_calc
    !integer :: nrsubst, errorflag, phasetype(5), iFlash, iter
    !integer :: i, iphase_try
    !double precision:: rhovap_est, rholiq_est
    !double precision :: wm_mix
    !double precision, dimension(5) :: wm_phase
    !
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !!input to lower case
    !call uppertolower_char (input, len(input))
    !
    !!call setup
    !errorflag = 0
    !xvap = 0.D0
    !xliq = 0.D0
    !gl%transport = 0       ! needed to check the range of validity for transport properties; 0 for every other property
    !gl%unitin = 1
    !
    !
    !call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator, errorflag)
    !
    !! irregular input, or problems with fluid file or reference state
    !if (errorflag /= 0) then
    !    xliq = errorflag
    !    xvap = errorflag
    !    propliq = errorflag
    !    propvap = errorflag
    !    vapfrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !
    !!Check whether all necessary input variables have (proper) values
    !!---------------------------
    !!check input variables (T,p,h,s)
    !!---------------------------
    !if (input == 'tp') then
    !    if (prop1 <= 0.D0) then
    !        errorflag = - 9911
    !    end if
    !    if (prop2 <= 0.D0) then
    !        errorflag = - 9933
    !    end if
    !end if
    !if (input == 'ph') then
    !    if (prop1 <= 0.D0) then
    !        errorflag = - 9933
    !    end if
    !end if
    !if (input == 'ps') then
    !    if (prop1 <= 0.D0) then
    !        errorflag = - 9933
    !    end if
    !end if
    !if (input == 'td') then
    !    !if (prop1 <= 0.D0) then
    !    !    errorflag = - 9911
    !    !    input=inptorig
    !    !    return
    !    !end if
    !    !if (prop1 <= 0.D0) then
    !    !    errorflag = - 9922
    !    !    input=inptorig
    !    !    return
    !    !end if
    !    errorflag = -9955
    !    xliq = errorflag
    !    xvap = errorflag
    !    propliq = errorflag
    !    propvap = errorflag
    !    vapfrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !if (input == 'tvap') then
    !    if (prop1 <= 0.D0) then
    !        errorflag = - 9911
    !    end if
    !    if (prop2 <= 0.D0) then
    !        errorflag = - 9933
    !    end if
    !end if
    !if (input == 'tliq') then
    !    if (prop1 <= 0.D0) then
    !        errorflag = - 9911
    !    end if
    !    if (prop2 <= 0.D0) then
    !        errorflag = - 9933
    !    end if
    !end if
    !if (input == 'pvap') then
    !    if (prop1 <= 0.D0) then
    !        errorflag = - 9933
    !    end if
    !    if (prop2 <= 0.D0) then
    !        errorflag = - 9911
    !    end if
    !end if
    !if (input == 'pliq') then
    !    if (prop1 <= 0.D0) then
    !        errorflag = - 9933
    !    end if
    !    if (prop2 <= 0.D0) then
    !        errorflag = - 9911
    !    end if
    !end if
    !
    !!read moles from the character variable and store it in a double type variable
    !call moles_incheck(gl,xliq_in, errorflag)
    !xliq = xliq_in
    !if (errorflag /= 0) then
    !    xliq = errorflag
    !    xvap = errorflag
    !    propliq = errorflag
    !    propvap = errorflag
    !    vapfrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !call moles_incheck(gl,xvap_in, errorflag)
    !xvap = xvap_in
    !if (errorflag /= 0) then
    !    xliq = errorflag
    !    xvap = errorflag
    !    propliq = errorflag
    !    propvap = errorflag
    !    vapfrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !
    !!---------------------------
    !!check initial estimates for the compositions
    !do i = 1, gl%ncomp
    !    if ((xliq(i) <= 0.D0) .or. (xvap(i) <= 0.D0)) then
    !        errorflag = -9951
    !    end if
    !end do
    !!---------------------------
    !!check iphase_try
    !if (iphase_try < 0 .or. iphase_try > 6) then
    !    !Set the iphase_try such that VLE is assumed (iphase_try = 5)
    !    iphase_try = 5
    !end if
    !!---------------------------
    !
    !if (errorflag /= 0) then
    !    xliq = errorflag
    !    xvap = errorflag
    !    propliq = errorflag
    !    propvap = errorflag
    !    vapfrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !
    !! ---------------------------------------------------------------------
    !!  TD-FLASH
    !! ---------------------------------------------------------------------
    !!if (input == 'td') then
    !!    t = prop1
    !!    d = prop2
    !!    if (ncomp == 1) then
    !!        nrsubst = 1     ! pure fluid
    !!        call PhaseDet_pure(p, T, d, dvap, dliq, phase, vapfrac, nrofphases, errorflag)
    !!        xliq = moles
    !!        xvap = moles
    !!    else
    !!        nrsubst = 0     ! mixture
    !!        errorflag = -9902
    !!        return
    !!    end if
    !!end if
    !
    !
    !! ---------------------------------------------------------------------
    !!  TP-FLASH
    !! ---------------------------------------------------------------------
    !if (input == 'tp') then
    !    t = prop1
    !    p = prop2
    !    if (gl%ncomp == 1) then
    !        errorflag = -9901 !pure fluid
    !    else
    !        nrsubst = 0     ! mixture
    !        !call PhaseDet(p, t, moles, rho, x_Phase, phasetype, vapfrac, nrofphases, errorflag)
    !        call Flash_pT_calc(gl,p, t, moles, xvap, xliq, rhovap_est, rholiq_est, vapfrac, iPhase_try, errorflag, iter)
    !        dvap = gl%rho_vap
    !        dliq = gl%rho_liq
    !    end if
    !end if
    !
    !! ---------------------------------------------------------------------
    !!  PH-FLASH
    !! ---------------------------------------------------------------------
    !if (input == 'ph') then
    !    p = prop1
    !    h = prop2
    !    if (gl%ncomp == 1) then
    !        nrsubst = 1  !pure fluid
    !        errorflag = -9901
    !    else
    !        nrsubst = 0     ! mixture
    !        !call PhaseDet_ph(p, t, moles, rho, x_Phase, phasetype, vapfrac, h, nrofphases, errorflag)
    !        call Flash_ph(gl,p, t, moles, xvap, xliq, rhovap_est, rholiq_est, vapfrac, h, iPhase_try, errorflag, iter)
    !        dvap = rho(phasetype(1))
    !        dliq = rho(phasetype(2))
    !    end if
    !end if
    !
    !! ---------------------------------------------------------------------
    !!  PS-FLASH
    !! ---------------------------------------------------------------------
    !if (input == 'ps') then
    !    p = prop1
    !    s = prop2
    !    if (gl%ncomp == 1) then
    !        nrsubst = 1     ! pure fluid
    !        errorflag = -9901
    !    else
    !        nrsubst = 0     ! mixture
    !        !call PhaseDet_ps(p, t, moles, rho, x_Phase, phasetype, vapfrac, s, nrofphases, errorflag)
    !        call Flash_ps(gl,p, t, moles, xvap, xliq, rhovap_est, rholiq_est, vapfrac, s, iPhase_try, errorflag, iter)
    !        dvap = rho(phasetype(1))
    !        dliq = rho(phasetype(2))
    !    end if
    !end if
    !
    !! ---------------------------------------------------------------------
    !!  OTHER INPUT CODES
    !! ---------------------------------------------------------------------
    !if ((input == 'tliq') .or. (input == 'tvap') .or. (input == 'pliq') .or. (input == 'pvap')) then
    !    if (gl%ncomp == 1) then
    !        nrsubst = 1     ! pure fluid
    !        xliq = moles
    !        xvap = moles
    !        if (input(1:1) == 'P') then   ! vapor pressure given
    !            p = prop1
    !            t = prop2
    !            dvap = rhovap_est
    !            dliq = rholiq_est
    !            call Flash_Pure_PhaseBoundary(gl,p, t, dvap, dliq, 2, errorflag, iter, nrsubst)
    !            !call VLEpurePres(p,t,dvap,dliq,errorflag, 1)
    !        else                        ! Sat. Temperature given
    !            t = prop1
    !            p = prop2
    !            dvap = rhovap_est
    !            dliq = rholiq_est
    !            call Flash_Pure_PhaseBoundary(gl,p, t, dvap, dliq, 1, errorflag, iter, nrsubst)
    !            !call VLEpure(p,t,dvap,dliq,errorflag, 1)
    !        end if
    !        vapfrac = 1.d0
    !        if (input(2:4) == 'LIQ') vapfrac = 0.d0
    !    else
    !        nrsubst = 0     ! mixture
    !        !   - BUBBLE POINT: P AND x' VEXTOR GIVEN   --  iFlash = 1
    !        !   - DEW POINT:    P AND x" VECTOR GIVEN   --  iFlash = 2
    !        !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
    !        !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
    !        iFlash = 0
    !        if (input == 'pliq') then
    !            iFlash = 1
    !            p = prop1
    !            t = prop2
    !            xliq = moles
    !        end if
    !        if (input == 'pvap') then
    !            iFlash = 2
    !            p = prop1
    !            t = prop2
    !            xvap = moles
    !        end if
    !        if (input == 'tliq') then
    !            iFlash = 3
    !            t = prop1
    !            p = prop2
    !            xliq = moles
    !        end if
    !        if (input == 'tvap') then
    !            iFlash = 4
    !            t = prop1
    !            p = prop2
    !            xvap = moles
    !        end if
    !        dvap = rhovap_est
    !        dliq = rholiq_est
    !        !call Flash_PhaseBoundary(p, t, moles, xvap, xliq, dvap, dliq, vapfrac, iFlash, 0, errorflag, iter)
    !        call Flash_PhaseBoundary_calc(gl,p, t, moles, xvap, xliq, dvap, dliq, vapfrac, iFlash,&
    !            & iPhase_try, 1, errorflag, iter)
    !        if (errorflag == 0) then
    !            dliq = gl%rho_liq
    !            dvap = gl%rho_vap
    !        end if
    !    end if
    !end if
    !
    !! ---------------------------------------------------------------------
    !! calculate all thermodynamic properties
    !! ---------------------------------------------------------------------
    !if (errorflag == 0) then
    !    !if ((0.d0 <= vapfrac) .and. (vapfrac <= 1.d0) .and. (xvap(1) /= 0.d0)) then   ! two-phase
    !    ! calculate the properties of the liquid phase
    !    gl%molfractions = xliq
    !    call reduced_parameters_calc(gl,t) !Dummy temperature 300 K for the SRK
    !    propliq(1) = t
    !    propliq(2) = p
    !    propliq(3) = dliq
    !    propliq(4) = u_calc(gl,t, dliq, nrsubst)
    !    propliq(5) = h_calc(gl,t, dliq, nrsubst)
    !    propliq(6) = s_calc(gl,t, dliq, nrsubst)
    !    propliq(7) = g_calc(gl,t, dliq, nrsubst)
    !    propliq(8) = A_CALC(gl,t, dliq, nrsubst)
    !    propliq(9) = cp_calc(gl,t, dliq, nrsubst)
    !    propliq(10) = cv_calc(gl,t, dliq, nrsubst)
    !    propliq(11) = ws_calc(gl,t, dliq, nrsubst)
    !
    !    !call wm_mix_calc(gl, wm_mix)
    !    !wm_phase(1) = wm_mix
    !
    !    ! calculate the properties of the vapor phase
    !    gl%molfractions = xvap
    !    call reduced_parameters_calc(gl,t) !Dummy temperature 300 K for the SRK
    !    propvap(1) = t
    !    propvap(2) = p
    !    propvap(3) = dvap
    !    propvap(4) = u_calc(gl,t, dvap, nrsubst)
    !    propvap(5) = h_calc(gl,t, dvap, nrsubst)
    !    propvap(6) = s_calc(gl,t, dvap, nrsubst)
    !    propvap(7) = g_calc(gl,t, dvap, nrsubst)
    !    propvap(8) = A_CALC(gl,t, dvap, nrsubst)
    !    propvap(9) = cp_calc(gl,t, dvap, nrsubst)
    !    propvap(10) = cv_calc(gl,t, dvap, nrsubst)
    !    propvap(11) = ws_calc(gl,t, dvap, nrsubst)
    !
    !    !call wm_mix_calc(gl, wm_mix)
    !    !wm_phase(2) = wm_mix
    !
    !else
    !    xliq = errorflag
    !    xvap = errorflag
    !    propliq = errorflag
    !    propvap = errorflag
    !    vapfrac = errorflag
    !end if
    !
    !
    !!If the SRK is used, check whether an ideal part exists
    !!Andreas November 2013
    !Do i = 1, gl%ncomp
    !    if(gl%Eq_type(1) == 2) then
    !        !---------------------------------------------------------
    !        !If constants B and C for the SRK cp0 model are 0, no model is implemented
    !        if ((dabs(gl%B_cv0(i)) < 1.D-8) .and. (dabs(gl%C_cv0(i)) < 1.D-8)) then
    !            !No ideal part for SRK exists
    !            errorflag = -2908
    !            propvap(4) = errorflag
    !            propvap(5) = errorflag
    !            propvap(6) = errorflag
    !            propvap(7) = errorflag
    !            propvap(8) = errorflag
    !            propvap(9) = errorflag
    !            propvap(10) = errorflag
    !            propvap(11) = errorflag
    !            propliq(4) = errorflag
    !            propliq(5) = errorflag
    !            propliq(6) = errorflag
    !            propliq(7) = errorflag
    !            propliq(8) = errorflag
    !            propliq(9) = errorflag
    !            propliq(10) = errorflag
    !            propliq(11) = errorflag
    !        end if
    !    end if
    !end do
    !
    !input=gl%inptorig
    !
    !
    !End subroutine FLASH_EXPERT
    !
    !
    !subroutine FLASH_NC_3P_EXPERT_STDCALL (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, xliq1_in, xliq2_in, xvap_in, &
    !    & rholiq1_est, rholiq2_est, rhovap_est, xliq1, xliq2, xvap, propliq1, propliq2, propvap, phasefrac, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "FLASH_NC_3P_EXPERT_STDCALL" :: FLASH_NC_3P_EXPERT_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: FLASH_NC_3P_EXPERT_STDCALL
    !

    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !!define input variables
    !character (12) :: input
    !character(20) :: unitdefinition
    !double precision :: prop1, prop2
    !character(30), dimension(30) :: fluids
    !double precision, dimension(30) :: moles
    !integer, dimension(30) :: EOS_indicator
    !integer :: MIX_indicator
    !character (255) :: path
    !double precision, dimension(30):: xliq1_in, xliq2_in, xvap_in
    !double precision:: rhovap_est, rholiq1_est, rholiq2_est
    !integer :: iphase_try
    !!define output variables
    !double precision, dimension(30), intent(out) :: xliq1, xliq2
    !double precision, dimension(30), intent(out) :: xvap
    !double precision, dimension(30), intent(out) :: propliq1, propliq2
    !double precision, dimension(30), intent(out) :: propvap
    !double precision, dimension(3), intent(out) :: Phasefrac
    !
    !integer :: errorflag
    !
    !
    !
    !call FLASH_NC_3P_EXPERT (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, xliq1_in, xliq2_in, xvap_in, &
    !    & rholiq1_est, rholiq2_est, rhovap_est, xliq1, xliq2, xvap, propliq1, propliq2, propvap, phasefrac, errorflag, handle)
    !
    !end subroutine

    !SUBROUTINE FLASH_NC_3P_EXPERT (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, xliq1_in, xliq2_in, xvap_in, &
    !    & rholiq1_est, rholiq2_est, rhovap_est, xliq1, xliq2, xvap, propliq1, propliq2, propvap, Phasefrac, errorflag, handle)
    !
    !



    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !
    !double precision :: p, t
    !double precision, intent(inout) :: prop1, prop2
    !character(30), dimension(30) :: fluids
    !double precision, dimension(30) :: moles
    !character (255), intent(inout) :: path
    !integer, dimension(30) :: EOS_indicator
    !integer :: MIX_indicator
    !
    !double precision, dimension(30):: xliq1_in, xliq2_in, xvap_in
    !
    !double precision, dimension(30), intent(out) :: xvap, xliq1, xliq2
    !double precision, dimension(3), intent(out) :: Phasefrac
    !
    !double precision, dimension(5) :: rho, wm_phase
    !integer :: errorflag, iter, iFlash, nrsubst
    !
    !double precision, intent(in):: rhovap_est, rholiq1_est, rholiq2_est
    !double precision :: dvap, dliq1, dliq2
    !
    !double precision, dimension(30), intent(out) :: propliq1
    !double precision, dimension(30), intent(out) :: propliq2
    !double precision, dimension(30), intent(out) :: propvap
    !
    !character(12), intent(inout) :: input
    !double precision:: u_calc, h_calc, s_calc, g_calc, A_CALC, cp_calc, cv_calc, ws_calc, wm_mix
    !character(20) :: unitdefinition
    !
    !integer:: i
    !
    !!DEC$ ATTRIBUTES DLLEXPORT :: FLASH_NC_3P_EXPERT
    !
    !!input to lower case
    !call uppertolower_char(input, len(input))
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !gl%transport = 0       ! needed to check the range of validity for transport properties; 0 for every other property
    !
    !!call setup
    !errorflag = 0
    !xvap = 0.D0
    !xliq1 = 0.D0
    !xliq2 = 0.D0
    !gl%unitin = 1
    !
    !
    !call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator, errorflag)
    !
    !! irregular input, or problems with fluid file or reference state
    !if (errorflag /= 0) then
    !    xliq1 = errorflag
    !    xliq2 = errorflag
    !    xvap = errorflag
    !    propliq1 = errorflag
    !    propliq2 = errorflag
    !    propvap = errorflag
    !    Phasefrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !
    !
    !if (input == 'tp') then
    !    iflash = 7
    !    t = prop1
    !    p = prop2
    !    if ((Phasefrac(1) <= 0.D0) .and. (Phasefrac(2) <= 0.D0) .and. (Phasefrac(3) <= 0.D0)) then
    !        Phasefrac = 1.D0 / 3.D0
    !    end if
    !elseif (input == 'tl2') then
    !    iflash = 6
    !    t = prop1
    !    p = prop2
    !    if ((Phasefrac(1) <= 0.D0) .and. (Phasefrac(2) <= 0.D0) .and. (Phasefrac(3) <= 0.D0)) then
    !        Phasefrac = 0.5D0
    !    end if
    !elseif (input == 'tvap') then
    !    iflash = 5
    !    t = prop1
    !    p = prop2
    !    if ((Phasefrac(1) <= 0.D0) .and. (Phasefrac(2) <= 0.D0) .and. (Phasefrac(3) <= 0.D0)) then
    !        Phasefrac = 0.5D0
    !    end if
    !elseif (input == 'tliq') then
    !    iflash = 4
    !    t = prop1
    !    p = prop2
    !    if ((Phasefrac(1) <= 0.D0) .and. (Phasefrac(2) <= 0.D0) .and. (Phasefrac(3) <= 0.D0)) then
    !        Phasefrac = 0.5D0
    !    end if
    !elseif (input == 'pl2') then
    !    iflash = 3
    !    p = prop1
    !    t = prop2
    !    if ((Phasefrac(1) <= 0.D0) .and. (Phasefrac(2) <= 0.D0) .and. (Phasefrac(3) <= 0.D0)) then
    !        Phasefrac = 0.5D0
    !    end if
    !elseif (input == 'pvap') then
    !    iflash = 2
    !    p = prop1
    !    t = prop2
    !    if ((Phasefrac(1) <= 0.D0) .and. (Phasefrac(2) <= 0.D0) .and. (Phasefrac(3) <= 0.D0)) then
    !        Phasefrac = 0.5D0
    !    end if
    !elseif (input == 'pliq') then
    !    iflash = 1
    !    p = prop1
    !    t = prop2
    !    if ((Phasefrac(1) <= 0.D0) .and. (Phasefrac(2) <= 0.D0) .and. (Phasefrac(3) <= 0.D0)) then
    !        Phasefrac = 0.5D0
    !    end if
    !else
    !    errorflag = -9955
    !    xliq1 = errorflag
    !    xliq2 = errorflag
    !    xvap = errorflag
    !    propliq1 = errorflag
    !    propliq2 = errorflag
    !    propvap = errorflag
    !    Phasefrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !
    !
    !if (t <= 0.D0) then
    !    errorflag = - 9911
    !end if
    !if (p <= 0.D0) then
    !    errorflag = - 9933
    !end if
    !
    !!read moles from the character variable and store it in a double type variable
    !call moles_incheck(gl,xliq1_in, errorflag)
    !xliq1 = xliq1_in
    !if (errorflag /= 0) then
    !    xliq1 = errorflag
    !    xliq2 = errorflag
    !    xvap = errorflag
    !    propliq1 = errorflag
    !    propliq2 = errorflag
    !    propvap = errorflag
    !    Phasefrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !call moles_incheck(gl,xliq2_in, errorflag)
    !xliq2 = xliq2_in
    !if (errorflag /= 0) then
    !    xliq1 = errorflag
    !    xliq2 = errorflag
    !    xvap = errorflag
    !    propliq1 = errorflag
    !    propliq2 = errorflag
    !    propvap = errorflag
    !    Phasefrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !call moles_incheck(gl,xvap_in, errorflag)
    !xvap = xvap_in
    !if (errorflag /= 0) then
    !    xliq1 = errorflag
    !    xliq2 = errorflag
    !    xvap = errorflag
    !    propliq1 = errorflag
    !    propliq2 = errorflag
    !    propvap = errorflag
    !    Phasefrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !!---------------------------
    !!check initial estimates for the compositions
    !do i = 1, gl%ncomp
    !    if ((xliq2(i) <= 0.D0) .or. (xliq1(i) <= 0.D0) .or. (xvap(i) <= 0.D0)) then
    !        errorflag = -9951
    !    end if
    !end do
    !!---------------------------
    !
    !if (gl%ncomp < 2) then
    !    errorflag = -9901
    !end if
    !
    !if (errorflag /= 0) then
    !    xliq1 = errorflag
    !    xliq2 = errorflag
    !    xvap = errorflag
    !    propliq1 = errorflag
    !    propliq2 = errorflag
    !    propvap = errorflag
    !    Phasefrac = errorflag
    !    input=gl%inptorig
    !    return
    !end if
    !
    !
    !call ptflash_NC_3P(gl,p, t, moles, rho, xvap, xliq1, xliq2, rhovap_est, &
    !    & rholiq1_est, rholiq2_est, Phasefrac, iFlash, iter, errorflag)
    !
    !dvap = rho(1)
    !dliq1 = rho(2)
    !dliq2 = rho(3)
    !
    !! ---------------------------------------------------------------------
    !! calculate all thermodynamic properties
    !! ---------------------------------------------------------------------
    !if (errorflag == 0) then
    !    !if ((0.d0 <= vapfrac) .and. (vapfrac <= 1.d0) .and. (xvap(1) /= 0.d0)) then   ! two-phase
    !    ! calculate the properties of the liquid phase 1
    !    nrsubst = 0
    !    gl%molfractions = xvap
    !    call reduced_parameters_calc(gl,T)
    !    propvap(1) = t
    !    propvap(2) = p
    !    propvap(3) = dvap
    !    propvap(4) = u_calc(gl,t, dvap, nrsubst)
    !    propvap(5) = h_calc(gl,t, dvap, nrsubst)
    !    propvap(6) = s_calc(gl,t, dvap, nrsubst)
    !    propvap(7) = g_calc(gl,t, dvap, nrsubst)
    !    propvap(8) = A_CALC(gl,t, dvap, nrsubst)
    !    propvap(9) = cp_calc(gl,t, dvap, nrsubst)
    !    propvap(10) = cv_calc(gl,t, dvap, nrsubst)
    !    propvap(11) = ws_calc(gl,t, dvap, nrsubst)
    !
    !    !mol specific handle
    !    ! if (gl%unitin == 2) then
    !    ! call wm_mix_calc(gl, wm_mix)
    !    !
    !    !     wm_phase(1) = wm_mix
    !    ! end if
    !    !
    !    ! calculate the properties of the liquid phase 2
    !    gl%molfractions = xliq1
    !
    !
    !    call reduced_parameters_calc(gl,t)
    !    propliq1(1) = t
    !    propliq1(2) = p
    !    propliq1(3) = dliq1
    !    propliq1(4) = u_calc(gl,t, dliq1, nrsubst)
    !    propliq1(5) = h_calc(gl,t, dliq1, nrsubst)
    !    propliq1(6) = s_calc(gl,t, dliq1, nrsubst)
    !    propliq1(7) = g_calc(gl,t, dliq1, nrsubst)
    !    propliq1(8) = A_CALC(gl,t, dliq1, nrsubst)
    !    propliq1(9) = cp_calc(gl,t, dliq1, nrsubst)
    !    propliq1(10) = cv_calc(gl,t, dliq1, nrsubst)
    !    propliq1(11) = ws_calc(gl,t, dliq1, nrsubst)
    !
    !    !mol specific handle
    !    !if (gl%unitin == 2) then
    !    !call wm_mix_calc(gl, wm_mix)
    !    !
    !    !    wm_phase(2) = wm_mix
    !    !end if
    !
    !
    !    ! calculate the properties of the vapor phase
    !    gl%molfractions = xliq2
    !    call reduced_parameters_calc(gl,t)
    !    propliq2(1) = t
    !    propliq2(2) = p
    !    propliq2(3) = dliq2
    !    propliq2(4) = u_calc(gl,t, dliq2, nrsubst)
    !    propliq2(5) = h_calc(gl,t, dliq2, nrsubst)
    !    propliq2(6) = s_calc(gl,t, dliq2, nrsubst)
    !    propliq2(7) = g_calc(gl,t, dliq2, nrsubst)
    !    propliq2(8) = A_CALC(gl,t, dliq2, nrsubst)
    !    propliq2(9) = cp_calc(gl,t, dliq2, nrsubst)
    !    propliq2(10) = cv_calc(gl,t, dliq2, nrsubst)
    !    propliq2(11) = ws_calc(gl,t, dliq2, nrsubst)
    !
    !    !mol specific handle
    !    ! if (gl%unitin == 2) then
    !    ! call wm_mix_calc(gl, wm_mix)
    !    !     wm_phase(3) = wm_mix
    !    ! end if
    !else
    !    xliq1 = errorflag
    !    xliq2 = errorflag
    !    xvap = errorflag
    !    propliq1 = errorflag
    !    propliq2 = errorflag
    !    propvap = errorflag
    !    Phasefrac = errorflag
    !end if
    !
    !
    !!If the SRK is used, check whether an ideal part exists
    !!Andreas November 2013
    !Do i = 1, gl%ncomp
    !    if(gl%Eq_type(1) == 2) then
    !        !---------------------------------------------------------
    !        !If constants B and C for the SRK cp0 model are 0, no model is implemented
    !        if ((dabs(gl%B_cv0(i)) < 1.D-8) .and. (dabs(gl%C_cv0(i)) < 1.D-8)) then
    !            !No ideal part for SRK exists
    !            errorflag = -2908
    !            propvap(4) = errorflag
    !            propvap(5) = errorflag
    !            propvap(6) = errorflag
    !            propvap(7) = errorflag
    !            propvap(8) = errorflag
    !            propvap(9) = errorflag
    !            propvap(10) = errorflag
    !            propvap(11) = errorflag
    !            propliq1(4) = errorflag
    !            propliq1(5) = errorflag
    !            propliq1(6) = errorflag
    !            propliq1(7) = errorflag
    !            propliq1(8) = errorflag
    !            propliq1(9) = errorflag
    !            propliq1(10) = errorflag
    !            propliq1(11) = errorflag
    !            propliq2(4) = errorflag
    !            propliq2(5) = errorflag
    !            propliq2(6) = errorflag
    !            propliq2(7) = errorflag
    !            propliq2(8) = errorflag
    !            propliq2(9) = errorflag
    !            propliq2(10) = errorflag
    !            propliq2(11) = errorflag
    !        end if
    !    end if
    !end do
    !
    !
    !
    !
    !input=gl%inptorig
    !
    !End subroutine FLASH_NC_3P_EXPERT



    !*******************************************************************************
    ! Subroutine for flash calculations up to 3 phases
    ! Vapor, liquid, solid, and hydrate may form
    ! Andreas, March 2014
    ! Modified for mixed hydrates: Property id in interface, changed output properties -> Langmuir Constants and Cage Occupancies moved to separate routines
    ! Sebastian, November 2017
    !*******************************************************************************
    !subroutine FLASH3_STDCALL (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, phasetype, phase_text, x_phase, prop_phase, &
    !    & prop_overall, lnfug_phase, Chempot_phase, phasefrac, prop_name_unit, errorflag, handle) !pmin, pmax,
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "FLASH3_STDCALL" :: FLASH3_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: FLASH3_STDCALL


    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !!define input variables
    !character (12) :: input
    !character(20) :: unitdefinition
    !double precision :: prop1, prop2
    !character (30), dimension(30) :: fluids
    !double precision, dimension(30) :: moles
    !integer, dimension(30) :: EOS_indicator
    !integer :: MIX_indicator
    !character (255) :: path
    !!define output variables
    !double precision, dimension(30, 5), intent(out) :: x_phase
    !double precision, dimension(30, 5), intent(out) :: prop_phase
    !double precision, dimension(30, 5), intent(out) :: lnfug_phase
    !double precision, dimension(30, 5), intent(out) :: ChemPot_phase
    !double precision, dimension(30), intent(out) :: prop_overall
    !double precision, dimension(5), intent(out) :: phasefrac
    !integer, dimension(5), intent(out) :: phasetype
    !character(4), dimension(5), intent(out) :: phase_text
    !character(30), dimension(flash_chardim,2), intent(out) :: prop_name_unit
    !!    double precision, optional, intent(inout) :: pmin, pmax
    !
    !integer :: errorflag
    !
    !call FLASH3 (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, phasetype, phase_text, x_phase, prop_phase, &
    !    & prop_overall, lnfug_phase, Chempot_phase, phasefrac, prop_name_unit, errorflag, handle) !pmin, pmax,
    !
    !!if ((dabs(pmin)<1.d-14).and.(dabs(pmax)<1.d-14)) then
    !!    pmin = 0.d0
    !!    pmax = 0.d0
    !!endif
    !end subroutine
    !
    !
    !!##################################################################
    !!FLASH3 function performs same calculations as TREND_EOS but gives properties of each present phase
    !!input          character(12)           specifies the given state point where the property has to be calculated, possible inputs are:
    !!                                       TP (prop1 is temperature, prop2 is pressure)
    !!                                       TD (prop1 is temperature, prop2 is density)
    !!                                       PS(prop1 is pressure, prop2 is entropy)
    !!                                       PH (prop1 is pressure, prop2 is enthalpy)
    !!                                       TLIQ (prop1 is temperature on bubble line, prop2 is
    !!                                       TVAP(prop1 is temperature on dew line, prop2 is ignored)
    !!                                       PLIQ (prop1 is pressure on bubble line, prop2 is ignored)
    !!                                       PVAP (prop1 is pressure on dew line, prop2 is ignored)
    !!prop1, prop2   double                  values for state point depending on input
    !!fluids         character(30),dim(30)   fluid/mixture calculations are conducted with
    !!moles          double,dim(30)          composition  calculations are conducted with
    !!EOS_indicator  integer,dim(30)         equation of state indicator calculations are conducted with
    !!                                       1: Helmholtz
    !!                                       2: SRK
    !!                                       3: PR
    !!                                       4: LKP
    !!                                       5: Generalized EOS (51: based on Alexandrov 52: based on Span/Wagner 53: based on Sun/Ely)
    !!                                       6: PCSAFT
    !!                                       7:
    !!                                       8: RKM
    !!                                       9: COSTALD
    !!MIX_indicator  integer                 indicator for mixing rule
    !!                                       1:      Helmholtz default,
    !!                                       110:    All components are mixed according to the linear mixing rule
    !!                                       111:    check if binary mix files are available, if not use linear mixing rule
    !!                                       120:    All components are mixed according to the Lorentz-Berthelot mixing rule.
    !!                                       121:    check if binary mix files are available, if not use Lorentz-Berthelot mixing rule.
    !!                                       2:      mixing rule according to SRK model
    !!                                       3:      mixing rule according to PR model
    !!                                       4:      mixing rule according to LKP model
    !!                                       5:      mixing rule according to Gen EOS model
    !!                                       6:      mixing rule according to PC SAFT model
    !!                                       7:
    !!                                       8:      mixing rule according to RKM model
    !!                                       9:      mixing rule according to COSTALD model
    !!path           chracter(255)           main directory of TREND
    !!unitdefintion  character(20)           choose molar/specific input and output
    !! 3-phase flash calculation with t-p, p-h, and p-s as input parameters possible
    !! This routine returns the compositions and all thermodynamic properties of the coexisting phases
    !!---------------------------------------------------
    !! Output parameters:
    !!   phasetype           Vector  phasetype contains the phase indicator number
    !!   phase_text          Vector  phase_text contains the abbreviations to all possible phases (VAP,LIQ1,LIQ2,SOL,HYD)
    !!E.g.: 2 phases are present: liquid and liquid equilibrium
    !!-- >  nrofphases = 2
    !!-- >  phasetype(1) = 2 (light liquid)
    !!-- >  phasetype(2) = 3 (heavy liquid)
    !!E.g.: 1 phase present: liquid
    !!-- >  nrofphases = 1
    !!-- >  phasetype(1) = 3 (heavy liquid) NOTE: In case of one phase liquid, heavy liquid is used!!
    !!phaseindicators: 1: Vapor, 2: light liquid, 3: heavy liquid, 4: solid, 5: hydrate
    !!   x_phase             composition matrix (5,30)
    !!       1   x_vap       - Vapor phase composition
    !!       2   x_liq1      - (lighter) liquid phase composition
    !!       3   x_liq2      - (heavier) liquid phase composition
    !!       4   x_sol       - solid phase composition (so far only CO2 and H2O are implemented)
    !!       5   x_hyd       - hydrate phase compsition
    !!   prop_phase          matrix (30,5) with all properties of all phases in equilibrium
    !!                       - T, P, D, U, H, S, G, A, CP, CV, WS,      Hyd_Nr, Hydrate_Structure, Overall_Small_Cage_Occup, Overall_Large_Cage_occup, phase_fraction
    !!                       - 1, 2, 3, 4, 5, 6, 7, 8,  9, 10, 11,          12                 13,                       14,                       15,             16,
    !!   prop_overall        - overall properties of phases in equilibrium
    !!   lnfug_phase         matrix (30,5) natural logarithm of the fugacities of all components in the phases
    !!   ChemPot_phase       matrix (30,5) chemical potential of all components in the phases
    !!   phasefrac           vector (5) that contains the phasefractions of all phases
    !!   prop_name_unit      matrix (flash_chardim,2) conatins the property name and the corresponding unit
    !!##################################################################
    !subroutine FLASH3 (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, phasetype, phase_text, x_phase, prop_phase,  prop_overall, lnfug_phase, Chempot_phase, phasefrac, prop_name_unit, errorflag, handle) !pmin, pmax,
    !!DEC$ ATTRIBUTES DLLEXPORT :: FLASH3
    !
    !! Adapted for molar <-> specific handling, Benedikt May 2018
    !




    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !character (12) :: input
    !character(20) :: unitdefinition
    !double precision :: prop1, prop2
    !character (30), dimension(30) :: fluids
    !double precision, dimension(30) :: moles
    !character (255) :: path
    !integer, dimension(30) :: EOS_indicator
    !integer :: MIX_indicator
    !!    double precision, optional, intent(inout) :: pmin, pmax
    !
    !double precision, dimension(30, 5), intent(out) :: x_phase
    !double precision, dimension(30, 5), intent(out) :: prop_phase
    !double precision, dimension(30, 5), intent(out) :: lnfug_phase
    !double precision, dimension(30, 5), intent(out) :: ChemPot_phase
    !double precision, dimension(30), intent(out) :: prop_overall
    !double precision, dimension(5), intent(out) :: phasefrac
    !integer, dimension(5), intent(out) :: phasetype
    !character(4), dimension(5), intent(out) :: phase_text
    !character(30), dimension(flash_chardim,2), intent(out) :: prop_name_unit
    !
    !double precision, dimension(30) :: molesorig, prop_vec
    !double precision, dimension(30) :: xvap, xliq1, xliq2, xsol, xhyd, xhelp, x_spec, x_overall
    !double precision:: d, t, p, dliq, dvap, d_phase
    !double precision, dimension(5) :: rho, wm_phase
    !double precision:: u_calc, h_calc, s_calc, g_calc, A_CALC, cp_calc, cv_calc, ws_calc, wm_mix, wm_overall
    !double precision:: u_DryIce, h_DryIce, s_DryIce, g_DryIce, f_DryIce, cp_DryIce
    !double precision:: u_WaterIce, h_WaterIce, s_WaterIce, g_WaterIce, f_WaterIce, cp_WaterIce
    !integer:: nrsubst, errorflag, nrofphases
    !integer :: i, j, k, converttype
    !
    !!Variables needed for the calculation of hydrate enthalpy
    !double precision, dimension(30) :: chem_pot, fug_g, x_fluid
    !double precision :: dchem_pot_dT_fluid, dchem_pot_dx_fluid, dfug_dT_fluid, dfug_dx_fluid, d_fluid
    !double precision :: h_hyd, s_hyd, u_hyd, f_hyd, g_hyd, ChemPot_hyd
    !double precision :: fug_CO2, fug_DryIce
    !integer :: pos_dryIce
    !
    !!New variables needed to calculate the chemical potentials and fugacities of the fluid phases
    !integer :: oir
    !double precision, dimension (30):: chem_pot_fluid
    !double precision, dimension (30):: lnfi_fluid
    !
    !!double precision, dimension(5) :: wm_phase
    !
    !!Pure Hydrates Model
    !!!Variables needed for the hydrate composition and Langmuir constants
    !!double precision, dimension(2) :: occup, CiJ, xH
    !
    !! Mixed hydrates 2015
    !double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double !<- CiJ and occup are used in the commented part ONLY
    !double precision, dimension(30):: xH, fug_gas!, occup_ls, occup_ld, occup_sms, occup_smd !<- xH and fug_gas are used in the commented part ONLY
    !!call my_allocate(gl)
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !!if ((present(pmin)).and.(present(pmax))) then
    !!    if ((dabs(pmin)>1.d-14).and.(dabs(pmax)>1.d-14)) then
    !!        gl%savebounds_p = .true.
    !!        gl%pmin_old = pmin
    !!        gl%pmax_old = pmax
    !!    elseif((dabs(pmin)<1.d-14).and.(dabs(pmax)<1.d-14)) then
    !!        pmin = 0.d0
    !!        pmax = 0.d0
    !!    end if
    !!end if
    !
    !!Initialize variables
    !x_phase = 0.D0
    !prop_phase = 0.D0
    !phasefrac = 0.D0
    !phasetype = 0
    !xvap = 0.D0
    !xliq1 = 0.D0
    !xliq2 = 0.D0
    !xsol = 0.D0
    !xhyd = 0.D0
    !chem_pot = 0.d0
    !fug_g = 0.D0
    !dchem_pot_dT_fluid = 0.d0
    !dchem_pot_dx_fluid = 0.d0
    !dfug_dT_fluid = 0.d0
    !dfug_dx_fluid = 0.d0
    !u_hyd = 0.D0
    !h_hyd = 0.D0
    !s_hyd = 0.d0
    !g_hyd = 0.d0
    !f_hyd = 0.d0
    !lnfug_phase = 0.D0
    !ChemPot_phase = 0.D0
    !prop_vec = 0.d0
    !errorflag = 0
    !
    !!Needed for the chemical potential --> get overall chemical potential
    !oir = 0
    !
    !!trim input
    !call uppertolower_char(input,len(input))
    !
    !call uppertolower_char(unitdefinition,len(unitdefinition))
    !
    !if(trim(unitdefinition) == "molar") then
    !    gl%unitin = 1
    !elseif(trim(unitdefinition) == "specific") then
    !    gl%unitin =2
    !else
    !    errorflag = -123456
    !    x_phase = errorflag
    !    prop_phase = errorflag
    !    phasefrac = errorflag
    !    phasetype = 0
    !    input=gl%inptorig
    !    return
    !end if
    !
    !!Andreas Jan 2015
    !gl%calc_ref = .true.  ! calculate reference state
    !gl%transport = 0       ! needed to check the range of validity for transport properties; 0 for every other property
    !
    !!call setup
    !errorflag = 0
    !call setup (gl,input, prop1, prop2, fluids, moles, path,EOS_indicator, mix_indicator, errorflag)
    !
    !!if (gl%check_solid .eqv. .true.) then
    !!   prop_id = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,0 /)
    !phase_text = (/'VAP','LIQ1','LIQ2','SOL','HYD'/)
    !prop_name_unit(:,1) = (/'Temperature','Pressure','Density','Int. Energy','Enthalpy','Entropy','Gibbs energy','Helmholtz energy','isob. Heat capacity','isoch. Heat capacity','speed of sound','Hydration number','Hydrate Structure', 'Overall Small Cage Occup','Overall Large Cage Occup','','phase fraction','X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15','X16','X17','X18','X19','X20'/)
    !if (gl%unitin == 1) then
    !    prop_name_unit(:,2) = (/'K','MPa','mol/m³','J/mol','J/mol','J/(mol K)','J/mol','J/mol','J/(mol K)','J/(mol K)','m/s','-','-','-','-','','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol','mol/mol'/)
    !elseif (gl%unitin == 2) then
    !    prop_name_unit(:,2) = (/'K','MPa','kg/m³','J/kg','J/kg','J/(kg K)','J/kg','J/kg','J/(kg K)','J/(kg K)','m/s','-','-','-','-','-','','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg','kg/kg'/)
    !endif
    !!else
    !!    prop_id = (/1,2,3,4,5,6,7,8,9,10,11,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,0,0,0,0,0 /)
    !!endif
    !
    !molesorig = gl%molfractions        ! be aware of gl%molfractions with molar / specific input
    !
    !
    !
    !! irregular input, or problems with fluid file or reference state
    !if (errorflag /= 0) then
    !    x_phase = errorflag
    !    prop_phase = errorflag
    !    phasefrac = errorflag
    !    phasetype = 0
    !    input=gl%inptorig
    !    return
    !end if
    !
    !
    !!Handle different input parameters and do all calculations
    !call inpt_handle (gl,input, prop1, prop2, p, t, d, dvap, dliq, moles, rho, x_Phase, phasetype, &
    !    & phasefrac, nrofphases, nrsubst, errorflag)
    !
    !! ---------------------------------------------------------------------
    !! calculate all thermodynamic properties
    !! ---------------------------------------------------------------------
    !if (errorflag == 0) then
    !
    !    prop_overall = 0.D0
    !    do i = 1, nrofphases
    !        prop_phase(1, phasetype(i)) = t
    !        prop_phase(2, phasetype(i)) = p
    !        prop_overall(1) = t
    !        prop_overall(2) = p
    !        if (phasetype(i) < 4) then  !Fluid phase
    !            ! calculate the properties of the first phase
    !            gl%molfractions = x_phase(:,phasetype(i))
    !            !! if((gl%unitin ==1) .and. (gl%unitstat == 2))
    !            !! call set_molarspecific(convert, input, moles, prop1, prop2, x_spec, wm_phase)
    !            !! gl%molefractions = x_spec
    !            !! end if
    !            call wm_mix_calc(gl, wm_mix)
    !            wm_phase(phasetype(i)) = wm_mix
    !            if (gl%ncomp > 1) then
    !                call reduced_parameters_calc(gl,t) !Dummy temperature 300 K for the SRK
    !            end if
    !            d_phase = rho(phasetype(i))
    !            prop_phase(3, phasetype(i)) = d_phase
    !            Do k = 1, gl%ncomp
    !                if(gl%Eq_type(k) > 1) then
    !                    !---------------------------------------------------------
    !                    !If constants B and C for the SRK cp0 model are 0, no model is implemented
    !                    if (((dabs(gl%B_cv0(k)) < 1.D-8) .and. (dabs(gl%C_cv0(k)) < 1.D-8))) then
    !                        prop_phase(4, phasetype(i)) =-2908
    !                        prop_phase(5, phasetype(i)) = -2908
    !                        prop_phase(6, phasetype(i)) = -2908
    !                        prop_phase(7, phasetype(i)) = -2908
    !                        prop_phase(8, phasetype(i)) = -2908
    !                        prop_phase(9, phasetype(i)) = -2908
    !                        prop_phase(10, phasetype(i)) = -2908
    !                        prop_phase(11, phasetype(i)) = -2908
    !                    else
    !                        prop_phase(4, phasetype(i)) = u_calc(gl,t, d_phase, nrsubst)
    !                        prop_phase(5, phasetype(i)) = h_calc(gl,t, d_phase, nrsubst)
    !                        prop_phase(6, phasetype(i)) = s_calc(gl,t, d_phase, nrsubst)
    !                        prop_phase(7, phasetype(i)) = g_calc(gl,t, d_phase, nrsubst)
    !                        prop_phase(8, phasetype(i)) = A_CALC(gl,t, d_phase, nrsubst)
    !                        prop_phase(9, phasetype(i)) = cp_calc(gl,t, d_phase, nrsubst)
    !                        prop_phase(10, phasetype(i)) = cv_calc(gl,t, d_phase, nrsubst)
    !                        prop_phase(11, phasetype(i)) = ws_calc(gl,t, d_phase, nrsubst)
    !                    end if
    !
    !                else
    !                    prop_phase(4, phasetype(i)) = u_calc(gl,t, d_phase, nrsubst)
    !                    prop_phase(5, phasetype(i)) = h_calc(gl,t, d_phase, nrsubst)
    !                    prop_phase(6, phasetype(i)) = s_calc(gl,t, d_phase, nrsubst)
    !                    prop_phase(7, phasetype(i)) = g_calc(gl,t, d_phase, nrsubst)
    !                    prop_phase(8, phasetype(i)) = A_CALC(gl,t, d_phase, nrsubst)
    !                    prop_phase(9, phasetype(i)) = cp_calc(gl,t, d_phase, nrsubst)
    !                    prop_phase(10, phasetype(i)) = cv_calc(gl,t, d_phase, nrsubst)
    !                    prop_phase(11, phasetype(i)) = ws_calc(gl,t, d_phase, nrsubst)
    !                end if
    !            end do
    !
    !
    !            !Calculate the chemical potential for all components in the fluid phase
    !            call Chempot_CALC(gl,t, d_phase, chem_pot_fluid, oir)
    !            ChemPot_phase(:,phasetype(i)) = chem_pot_fluid
    !
    !            !Calculate the fugacity for all components in the fluid phase
    !            call lnf_mix(gl,t, d_phase, p, lnfi_fluid)
    !            lnfug_phase(:,phasetype(i)) = lnfi_fluid
    !
    !        elseif (phasetype(i) == 4) then     !Solid (pure)
    !            gl%molfractions = x_phase(:,phasetype(i))
    !            if (gl%solidtype_akt_phase == 1) then         !Water
    !                prop_phase(3, phasetype(i)) = rho(phasetype(i))
    !                prop_phase(4, phasetype(i)) = u_WaterIce(gl,t, p)
    !                prop_phase(5, phasetype(i)) = h_WaterIce(gl,t, p)
    !                prop_phase(6, phasetype(i)) = s_WaterIce(gl,t, p)
    !                prop_phase(7, phasetype(i)) = g_WaterIce(gl,t, p)
    !                prop_phase(8, phasetype(i)) = f_WaterIce(gl,t, p)
    !                prop_phase(9, phasetype(i)) = cp_WaterIce(gl,t, p)
    !                prop_phase(10, phasetype(i)) = -12900
    !                prop_phase(11, phasetype(i)) = -12900
    !                wm_phase(phasetype(i)) = gl%wm(1)
    !            elseif (gl%solidtype_akt_phase == 2) then     !CO2
    !                prop_phase(3, phasetype(i)) = rho(phasetype(i))
    !                prop_phase(4, phasetype(i)) = u_DryIce(gl,t, p)
    !                prop_phase(5, phasetype(i)) = h_DryIce(gl,t, p)
    !                prop_phase(6, phasetype(i)) = s_DryIce(gl,t, p)
    !                prop_phase(7, phasetype(i)) = g_DryIce(gl,t, p)
    !                prop_phase(8, phasetype(i)) = f_DryIce(gl,t, p)
    !                prop_phase(9, phasetype(i)) = cp_DryIce(gl,t, p)
    !                prop_phase(10, phasetype(i)) = -12900
    !                prop_phase(11, phasetype(i)) = -12900
    !                if(gl%components(1) == 'water') then
    !                    wm_phase(phasetype(i)) = gl%wm(2)
    !                else
    !                    wm_phase(phasetype(i)) = gl%wm(1)
    !                end if
    !            else                                !No solid equation available
    !                prop_phase(3, phasetype(i)) = -9904
    !                prop_phase(4, phasetype(i)) = -9904
    !                prop_phase(5, phasetype(i)) = -9904
    !                prop_phase(6, phasetype(i)) = -9904
    !                prop_phase(7, phasetype(i)) = -9904
    !                prop_phase(8, phasetype(i)) = -9904
    !                prop_phase(9, phasetype(i)) = -9904
    !                prop_phase(10, phasetype(i)) = -9904
    !                prop_phase(11, phasetype(i)) = -9904
    !            end if
    !        elseif (phasetype(i) == 5) then     !Hydrate
    !            gl%molfractions = x_phase(:,phasetype(i))
    !            prop_phase(3, phasetype(i)) = rho(phasetype(i))
    !            prop_phase(9, phasetype(i)) = -12900        !cp
    !            prop_phase(10, phasetype(i)) = -12900       !cv
    !            prop_phase(11, phasetype(i)) = -12900       !w
    !
    !            !Get the T and x partial derivatives of the fugacity and the chemical potential of the fluid phase to
    !            !calculate the entropy and enthalpy of hydrates
    !            if (phasetype(1) < 4) then !first phase must be a fluid phase in order to calculate the needed derivatives
    !                x_fluid = x_phase(:,phasetype(1))
    !                d_fluid = rho(phasetype(1))
    !                xhyd = x_phase(:,phasetype(i))
    !                call hdrt_ancillary_hs(gl,t, d_fluid, p *1.D6, x_fluid, chem_pot, fug_g, errorflag)
    !
    !                call hdrt_Gibbs_energy(gl,xhyd, chem_pot, g_hyd, errorflag)
    !                if (errorflag == 0) then
    !                    call hdrt_entropy(gl,t, p*1.D6, xhyd, fug_g, s_hyd, errorflag)
    !                end if
    !                if (errorflag == 0) then
    !                    call hdrt_enthalpy(gl,t, p*1.D6, xhyd, chem_pot, fug_g,  h_hyd, errorflag)
    !                end if
    !                if (errorflag == 0) then
    !                    u_hyd = h_hyd - p * 1.D6 / rho(phasetype(i))
    !                    f_hyd = g_hyd - p * 1.D6 / rho(phasetype(i))
    !
    !                    prop_phase(4, phasetype(i)) = u_hyd
    !                    prop_phase(5, phasetype(i)) = h_hyd
    !                    prop_phase(6, phasetype(i)) = s_hyd
    !                    prop_phase(7, phasetype(i)) = g_hyd
    !                    prop_phase(8, phasetype(i)) = f_hyd
    !                else
    !
    !                    prop_phase(4, phasetype(i)) = -12900
    !                    prop_phase(5, phasetype(i)) = -12900
    !                    prop_phase(6, phasetype(i)) = -12900
    !                    prop_phase(7, phasetype(i)) = -12900
    !                    prop_phase(8, phasetype(i)) = -12900
    !
    !                end if
    !
    !                !Get the Langmuir constants, hydration number and cage occupancies for the gas hydrate phase
    !                !Andreas, Jan 2015
    !                call hdrt_mole_fract(gl,t,p*1.d6,fug_g(1),occup,CiJ,xH,occup_single, occup_double)
    !
    !                call wm_mix_calc(gl, wm_mix)
    !                wm_phase(phasetype(i)) = wm_mix
    !
    !                !pure hydrates
    !                !prop_phase(12, phasetype(i)) = xH(1)/xH(2)      !Hyd Nr
    !                !prop_phase(13, phasetype(i)) = occup(1)         !Small cage occup
    !                !prop_phase(14, phasetype(i)) = occup(2)         !Large cage occup
    !                !prop_phase(15, phasetype(i)) = CiJ(1) * 1.D6    !Small cage Langmuir constant 1/MPa
    !                !prop_phase(16, phasetype(i)) = CiJ(2) * 1.D6    !Large cage Langmuir constant 1/MPa
    !                !mixed hydrates
    !                prop_phase(12, phasetype(i)) = xH(1)/sum(xH(2:))      !Hyd Nr
    !                prop_phase(13, phasetype(i)) = gl%hdrt_structure_stable  !hydrate structure
    !                prop_phase(14, phasetype(i)) = sum(occup(1,:))         !Overall Small cage occup
    !                prop_phase(15, phasetype(i)) = sum(occup(2,:))         !Overall Large cage occup
    !                !prop_phase(16, phasetype(i)) = sum(CiJ(1,:)) * 1.D6    !Small cage Langmuir constant 1/MPa
    !                !prop_phase(17, phasetype(i)) = sum(CiJ(2,:)) * 1.D6    !Large cage Langmuir constant 1/MPa
    !
    !            else
    !                if ((gl%solidpos_akt_phase == 2) .and. (gl%ncomp == 2)) then !Try dry ice in equilibrium with hydrate
    !                    pos_DryIce = gl%solidtype_akt_phase
    !                    fug_g(1) = fug_DryIce(gl,t, p, pos_dryIce) *1.D6
    !                    if (fug_g(1) < 1.D-12) then
    !                        errorflag = -7779
    !                    else
    !                        !fug_CO2 = fug_g(1)
    !                        chem_pot(2) = g_DryIce(gl,t, p)
    !                        !call hdrt_chem_potent_w(t, p*1.d6, fug_CO2, ChemPot_hyd)
    !                        call hdrt_chem_potent_w(gl,t, p*1.d6, fug_g, ChemPot_hyd)
    !                        chem_pot(1) = ChemPot_hyd
    !
    !                        xhyd = x_phase(:,phasetype(i))
    !
    !                        call hdrt_enthalpy(gl,t, p*1.D6, xhyd, chem_pot, fug_g, h_hyd, errorflag)
    !                    end if
    !                else
    !                    errorflag = -12900
    !                end if
    !            end if
    !
    !        end if
    !
    !        prop_overall(3) = prop_overall(3) + phasefrac(phasetype(i)) / prop_phase(3,phasetype(i)) !Density
    !        prop_overall(4:8) = prop_overall(4:8) + phasefrac(phasetype(i)) * prop_phase(4:8,phasetype(i))
    !
    !    end do
    !    if (dabs(prop_overall(3)) > 1.D-14) then
    !        prop_overall(3) = 1.D0 / prop_overall(3) !Density
    !    end if
    !else
    !    x_phase = errorflag
    !    phasefrac = errorflag
    !    prop_phase = errorflag
    !    phasetype = errorflag
    !    input = gl%inptorig
    !    gl%molfractions = molesorig
    !    return
    !end if
    !
    !do i = 1, nrofphases
    !    !if solids are checked, map the molfractions back to original input
    !    if (gl%check_solid .and. (gl%ncomp > 1)) then
    !        xhelp = x_phase(:,phasetype(i))
    !        do j = 1, gl%ncomp
    !            x_phase(gl%mapping(j),phasetype(i)) = xhelp(j)
    !        end do
    !    end if
    !end do
    !
    !!If the SRK is used, check whether an ideal part exists
    !!Andreas November 2013
    !Do i = 1, gl%ncomp
    !    if(gl%Eq_type(1) == 2) then
    !        !---------------------------------------------------------
    !        !If constants B and C for the SRK cp0 model are 0, no model is implemented
    !        if ((dabs(gl%B_cv0(i)) < 1.D-8) .and. (dabs(gl%C_cv0(i)) < 1.D-8)) then
    !            !No ideal part for SRK exists
    !            errorflag = -2908
    !            prop_phase(4:11,:) = errorflag
    !        end if
    !    end if
    !end do
    !
    !
    !
    !if((gl%unitin==2) .and. (errorflag == 0)) then
    !
    !    converttype = 2         ! converting moles to weight fractions
    !
    !    do i=1,nrofphases
    !
    !        x_spec = x_phase(:,phasetype(i))
    !
    !        call convert_fractions(gl, converttype, wm_phase(phasetype(i)), x_spec)
    !
    !        x_phase(:,phasetype(i)) = x_spec
    !
    !    end do
    !
    !
    !    gl%molfractions = molesorig
    !    call wm_mix_calc(gl, wm_overall)
    !    call convert_fractions(gl, converttype, wm_overall, molesorig)
    !
    !
    !    do i = 1,nrofphases
    !
    !        prop_vec = prop_phase(:, phasetype(i))
    !        if (phasetype(i) == 5) then
    !            write(*,*) wm_phase(i)
    !            continue
    !        end if
    !        wm_mix = wm_phase(phasetype(i))
    !
    !        call convert_flash_spec(gl, prop_vec,  wm_mix)
    !
    !        prop_phase(:, phasetype(i)) = prop_vec
    !
    !    end do
    !
    !    prop_vec = prop_overall
    !
    !    wm_mix = wm_overall
    !    call convert_flash_spec(gl, prop_vec,  wm_mix)
    !
    !    prop_overall = prop_vec
    !
    !end if
    !
    !
    !input=gl%inptorig
    !gl%molfractions = molesorig         !be aware of molar / specific composition  and phase fractions ...
    !
    !!if ((present(pmin)).and.(present(pmax)).and.(dabs(pmin)>1.d-14).and.(dabs(pmax)>1.d-14)) then
    !!    pmin = gl%pmin_old
    !!    pmax = gl%pmax_old
    !!end if
    !end subroutine FLASH3



    !!*******************************************************************************
    !! Subroutine for the calculation of the Henry constants in binary mixtures
    !!*******************************************************************************
    !subroutine HENRY_EOS_STDCALL (Temp, fluids, moles, EOS_indicator, MIX_indicator, path, H12, H21, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "HENRY_EOS_STDCALL" :: HENRY_EOS_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: HENRY_EOS_STDCALL
    !
    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !
    !!define input variables
    !double precision, intent(in) :: Temp
    !character (30), dimension(30) :: fluids
    !double precision, dimension (30) :: moles
    !integer, dimension (30) :: EOS_indicator
    !integer :: MIX_indicator
    !
    !character (255), intent(in) :: path
    !!define output variables
    !double precision, intent (out) :: H12, H21 !why no vector?
    !
    !integer :: errorflag
    !
    !call HENRY_EOS (Temp, fluids, moles, EOS_indicator, MIX_indicator, path, H12, H21, errorflag, handle)
    !
    !end subroutine
    !
    !subroutine HENRY_EOS (Temp, fluids, moles, EOS_indicator, MIX_indicator, path, H12, H21, errorflag, handle)
    !!-------------------------------------------------------------------------------
    !!Subroutine for the calculation of the Henry constants of fluid 1 in fluid 2 H12
    !!and fluid 2 in fluid 1 (H21) at given temperature:
    !!H12 = lim(x1 -> 0) phi1_L * ps,2
    !!H21 = lim(x2 -> 0) phi2_L * ps,1
    !!The calculation is restricted to binary mixtures
    !!Andreas Jäger, March 2017
    !!-------------------------------------------------------------------------------
    !!DEC$ ATTRIBUTES DLLEXPORT :: HENRY_EOS
    !

    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !
    !!define input/outpit variables
    !double precision :: Temp
    !character (30), dimension(30) :: fluids
    !double precision, dimension (30) :: moles
    !character (255) :: path
    !integer, dimension (30) :: EOS_indicator
    !integer :: MIX_indicator
    !double precision:: H12, H21
    !
    !!Variables for setup
    !double precision :: prop1, prop2        !Dummy variables
    !integer :: errorflag
    !character (12) :: input
    !
    !logical:: H12_calc, H21_calc
    !
    !!Variables for the calculation of the pure component vapor pressure
    !double precision :: press
    !double precision :: rhovap_est, rholiq_est   ![mol/m³]
    !integer :: iFlash, nrsubst
    !integer :: iter
    !
    !!Variables for the calculation of the fugacity coefficient of the infinitly diluted component in the solvent
    !double precision:: dens_liq
    !double precision, dimension(30):: fugcoef_mix
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !H12_calc = .false.
    !H12_calc = .false.
    !H12 = 0.D0
    !H21 = 0.D0
    !
    !!call setup
    !errorflag = 0
    !input = 'tp'    !Dummy
    !prop1 = Temp
    !prop2 = 1.0     !Dummy
    !call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator, errorflag)
    !
    !!Catch errors in setup
    !if (errorflag .ne. 0) then
    !    H12 = errorflag
    !    H21 = errorflag
    !    return
    !end if
    !
    !!Calculation only possible for binary mixtures, throw error in other cases
    !if (gl%ncomp .ne. 2) then
    !    errorflag = -5249
    !    H12 = errorflag
    !    H21 = errorflag
    !    return
    !end if
    !
    !!Check which of the Henry constants can be calculated. In order to be able to calculate a Henry constant, the temperature of the solvent needs to be subcritical but higher than the triple temperature.
    !If ((Temp < gl%tc(2)) .and. (Temp > gl%ttp(2))) then
    !    H12_calc = .true.
    !else
    !    H12 = -4409
    !end if
    !
    !If ((Temp < gl%tc(1)) .and. (Temp > gl%ttp(1))) then
    !    H21_calc = .true.
    !else
    !    H21 = -4409
    !end if
    !
    !!If no Henry constant can be calculated, quit with error
    !If ((H12_calc .eqv. .false.) .and. (H21_calc .eqv. .false.)) then
    !    errorflag = -4409
    !    return
    !end if
    !
    !if (H12_calc) then
    !
    !    !Set composition to pure component 2
    !    gl%molfractions(1) = 0.D0
    !    gl%molfractions(2) = 1.D0
    !    call reduced_parameters_calc(gl,Temp)
    !
    !    !Calculate vapor pressure of component 2
    !    rhovap_est = 0.D0
    !    rholiq_est = 0.D0
    !    iFlash = 1
    !    errorflag = 0
    !    nrsubst = 2
    !    press = 0.d0
    !    call Flash_Pure_PhaseBoundary (gl,press, Temp, rhovap_est, rholiq_est, iFlash, errorflag, iter, nrsubst)
    !    if (errorflag .ne. 0) then
    !        H12 = errorflag
    !        gl%molfractions = moles
    !        call reduced_parameters_calc(gl,Temp)
    !        return
    !    end if
    !
    !    !Calculate the fugacity coefficient of the infinitely diluted component in the solvent
    !    dens_liq = rholiq_est   !Density is the density of the saturated liquid of component 2
    !    call FUGCO_CALC_MIX(gl,Temp, dens_liq, fugcoef_mix)
    !    H12 = fugcoef_mix(1) * press
    !
    !end if
    !
    !if (H21_calc) then
    !
    !    !Set composition to pure component 1
    !    gl%molfractions(1) = 1.D0
    !    gl%molfractions(2) = 0.D0
    !    call reduced_parameters_calc(gl,Temp)
    !
    !    !Calculate vapor pressure of component 1
    !    rhovap_est = 0.D0
    !    rholiq_est = 0.D0
    !    iFlash = 1
    !    errorflag = 0
    !    nrsubst = 1
    !    press = 0.d0
    !    call Flash_Pure_PhaseBoundary (gl,press, Temp, rhovap_est, rholiq_est, iFlash, errorflag, iter, nrsubst)
    !    if (errorflag .ne. 0) then
    !        H21 = errorflag
    !        gl%molfractions = moles
    !        call reduced_parameters_calc(gl,Temp)
    !        return
    !    end if
    !
    !    !Calculate the fugacity coefficient of the infinitely diluted component in the solvent
    !    dens_liq = rholiq_est   !Density is the density of the saturated liquid of component 1
    !    call FUGCO_CALC_MIX(gl,Temp, dens_liq, fugcoef_mix)
    !    H21 = fugcoef_mix(2) * press
    !
    !end if
    !
    !gl%molfractions = moles
    !call reduced_parameters_calc(gl,Temp)
    !
    !end subroutine
    !!*******************************************************************************




    !!*******************************************************************************
    !!subroutine to calculate the thermodynamic factor according to Taylor and
    !!Note that the last molefraction of the mixture is replaced by 1.d0-sum of the
    !!remaining molefractions
    !subroutine GAMMATF_EOS_STDCALL (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, GAMMATF, errorflag, handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "GAMMATF_EOS_STDCALL" :: GAMMATF_EOS_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: GAMMATF_EOS_STDCALL
    !
    !use, intrinsic :: iso_c_binding

    !
    !implicit none
    !
    !!Handle for gl
    !type(c_ptr), intent(inout):: handle
    !
    !
    !!define input variables
    !character (12), intent(in) :: input
    !double precision, intent(in) :: prop1, prop2
    !character (30), dimension(30), intent(inout) :: fluids
    !double precision, dimension(30), intent(inout) :: moles
    !integer, dimension(30), intent(inout) :: EOS_indicator
    !integer :: MIX_indicator
    !character (255), intent(in) :: path
    !!define output variables
    !double precision, dimension(30,30), intent(out):: GAMMATF
    !
    !integer :: errorflag
    !
    !
    !call GAMMATF_EOS (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, GAMMATF, errorflag, handle)
    !
    !end subroutine




    !subroutine GAMMATF_EOS(input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, GAMMATF, errorflag, handle)
    !!subroutine for the calculation of the thermodynamic factor
    !! according to Taylor and Kooijman (1991)
    !!Output: GAMMATF as double precision, dimension(30)
    !!DEC$ ATTRIBUTES DLLEXPORT :: GAMMATF_EOS
    !


    !use, intrinsic :: iso_c_binding


    !
    !implicit none
    !
    !!-------------------------------------------------------------------------
    !!Handle for gl
    !type(type_gl) :: gl
    !type(c_ptr), intent(inout):: handle
    !!-------------------------------------------------------------------------
    !
    !!define input variables
    !character (12) :: input
    !double precision:: prop1, prop2
    !character (30), dimension(30), intent(inout) :: fluids
    !double precision, dimension(30), intent(inout) :: moles
    !integer, dimension(30), intent(inout) :: EOS_indicator
    !integer :: MIX_indicator
    !character (255) :: path
    !
    !double precision, dimension(30,30), intent(out):: GAMMATF
    !
    !!define variables for calculating
    !double precision :: t, d, p, h, s
    !double precision :: dvap
    !double precision :: dliq
    !double precision :: propvap
    !double precision :: propliq
    !double precision :: x_liq(30)
    !double precision :: x_vap(30)
    !
    !
    !integer :: errorflag
    !integer :: nrsubst, i
    !
    !double precision:: rho(5)           ! vector containing the densities of all phases
    !double precision:: x_Phase(30,5)    ! vector conaining the compositions of all phases
    !double precision:: vapfrac             ! molar vapor fraction
    !integer:: nrofphases                ! at this point: only 1 or 2 implemented
    !!Indicate which phases are present
    !integer:: phasetype(5)              !phasetype contains the phase indicator number
    !!E.g.: 2 phases are present: liquid and liquid equilibrium
    !!-- >  nrofphases = 2
    !!-- >  phasetype(1) = 2 (light liquid)
    !!-- >  phasetype(2) = 3 (heavy liquid)
    !
    !double precision, dimension(5) :: phasefrac
    !double precision, dimension(5) :: prop_phase
    !
    !double precision:: d_phase
    !
    !!input to lower case
    !call uppertolower_char(input, len(input))
    !
    !!control handle and fluid input fot memory usage
    !call control_fluids(gl,fluids,EOS_indicator,handle)
    !
    !gl%calc_ref = .false.  ! no need calculate reference state
    !gl%transport = 0       ! needed to check the range of validity for transport properties; 0 for every other property
    !
    !!call setup
    !errorflag = 0
    !
    !call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator, errorflag)
    !
    !!Catch errors in setup
    !if (errorflag .ne. 0) then
    !    GAMMATF = errorflag
    !    return
    !end if
    !
    !call reduced_parameters_calc(gl,300.d0)
    !
    !if (gl%ncomp .gt. 1) then
    !    call inpt_handle (gl,input, prop1, prop2, p, t, d, dvap, dliq, moles, rho, x_Phase, phasetype, &
    !        & phasefrac, nrofphases, nrsubst, errorflag)
    !else
    !    GAMMATF = -5248.d0   !property only available for mixtures
    !    return
    !end if
    !
    !if (errorflag .ne. 0) then
    !    GAMMATF = errorflag
    !    return
    !end if
    !
    !if (nrofphases > 1) then
    !    GAMMATF = -6666.d0
    !    return
    !end if
    !
    !if (phasetype(1) < 4) then  !fluid phase
    !    d_phase = rho(phasetype(1))
    !else
    !    GAMMATF = -9904.d0
    !    return
    !end if
    !
    !Call GAMMATF_CALC(gl,T, d_phase, GAMMATF)
    !
    !end subroutine


    !________________________________________________________________________________________________________________________________________________________________________________________________________________






    !***************************************************************************************
    ! Routine for Calculation of all thermophysical properties without any stability tests
    ! Stefan Herrig, November 2013
    !***************************************************************************************
    subroutine ALLPROP_HOM_STDCALL (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, phase_ind,TEOS,DEOS,PEOS,UEOS,HEOS,SEOS,GEOS,AEOS,CPEOS,CVEOS,WSEOS,BEOS,CEOS,CP0EOS,errorflag,handle)
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "ALLPROP_HOM_STDCALL" :: ALLPROP_HOM_STDCALL
    !DEC$ ATTRIBUTES stdcall, reference :: ALLPROP_HOM_STDCALL

    use, intrinsic :: iso_c_binding


    implicit none

    !Handle for gl
    type(c_ptr), intent(inout):: handle


    !define input variables
    character (12), intent(in) :: input
    character(20) :: unitdefinition
    double precision, intent(in) :: prop1, prop2
    character (30), dimension(30), intent(in) :: fluids
    double precision, dimension(30), intent(in) :: moles
    integer, dimension(30), intent(in) :: EOS_indicator
    integer:: MIX_indicator
    character (255), intent(in) :: path
    integer :: errorflag

    !define output variables
    double precision, intent(out) :: TEOS,DEOS,PEOS,UEOS,HEOS,SEOS,CVEOS,CPEOS,CP0EOS,WSEOS,GEOS,AEOS,BEOS,CEOS
    integer, intent(out) :: phase_ind



    call ALLPROP_HOM (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, phase_ind,TEOS,DEOS,PEOS,UEOS,HEOS,SEOS,GEOS,AEOS,CPEOS,CVEOS,WSEOS,BEOS,CEOS,CP0EOS,errorflag,handle)

    end subroutine
    SUBROUTINE ALLPROP_HOM (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition, phase_ind,&
        & TEOS,DEOS,PEOS,UEOS,HEOS,SEOS,GEOS,AEOS,CPEOS,CVEOS,WSEOS,BEOS,CEOS,CP0EOS,errorflag,gl_handle)
    !DEC$ ATTRIBUTES DLLEXPORT :: ALLPROP_HOM




    use, intrinsic :: iso_c_binding

    implicit none

    !-------------------------------------------------------------------------
    !Handle for gl
    type(type_gl) ,pointer:: gl
    type(c_ptr), intent(inout):: gl_handle
    !-------------------------------------------------------------------------

    !define input variables
    ! double precision :: P_CALC,U_CALC,H_CALC,S_CALC,CV_CALC,CP_CALC,CP0_CALC,WS_CALC,G_CALC,A_CALC,B_CALC,C_CALC
    double precision, intent(out) :: TEOS,DEOS,PEOS,UEOS,HEOS,SEOS,CVEOS,CPEOS,CP0EOS,WSEOS,GEOS,AEOS,BEOS,CEOS

    character (12) :: input!, unitdefinition
    character(20) :: unitdefinition
    double precision :: prop1, prop2, wm_mix
    character(30), dimension(30) :: fluids
    double precision, dimension(30) :: moles
    character (255) :: path
    integer, dimension(30) :: EOS_indicator
    integer:: MIX_indicator
    integer, intent(out) :: phase_ind

    !define variables for calculating
    double precision :: t, d, p


    integer :: errorflag
    integer :: unitin
    integer :: nrsubst
    integer :: ILIMITS
    integer :: i

    if(trim(unitdefinition) == "specific") then
        unitin=2
    elseif(trim(unitdefinition) == "molar") then
        unitin=1
    elseif(trim(unitdefinition) == "reduced") then
        unitin=3
    else
        errorflag=-123456
        return
    end if

    !control handle and fluid input fot memory usage
    call control_fluids(gl,input,fluids,moles,EOS_indicator,unitin,gl_handle,errorflag)
    ! irregular input, or problems with fluid file or reference state
    if (errorflag /= 0) then
        PEOS = errorflag
        DEOS = errorflag
        TEOS = errorflag
        UEOS = errorflag
        HEOS = errorflag
        SEOS = errorflag
        CVEOS = errorflag
        CPEOS = errorflag
        CP0EOS = errorflag
        WSEOS = errorflag
        GEOS = errorflag
        AEOS = errorflag
        BEOS = errorflag
        CEOS = errorflag
        input=gl%inptorig
        return
    end if


    gl%transport = 0       ! needed to check the range of validity for transport properties; 0 for every other property

    ILIMITS = 0

    errorflag = 0

    TEOS = 0.D0
    DEOS = 0.D0
    PEOS = 0.d0
    UEOS = 0.d0
    HEOS = 0.d0
    SEOS = 0.d0
    CVEOS = 0.d0
    CPEOS = 0.d0
    CP0EOS = 0.d0
    WSEOS = 0.d0
    GEOS = 0.d0
    AEOS = 0.d0
    BEOS = 0.d0
    CEOS = 0.d0

    !trim input
    call uppertolower_char(input, len(input))
    unitdefinition = trim(unitdefinition)

    if ((phase_ind < 0) .or. (phase_ind > 2)) then
        errorflag = -9917
        PEOS = errorflag
        DEOS = errorflag
        TEOS = errorflag
        UEOS = errorflag
        HEOS = errorflag
        SEOS = errorflag
        CVEOS = errorflag
        CPEOS = errorflag
        CP0EOS = errorflag
        WSEOS = errorflag
        GEOS = errorflag
        AEOS = errorflag
        BEOS = errorflag
        CEOS = errorflag
        input=gl%inptorig
        return
    end if


    !call setup - >  initializing, reading fluid files
    call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator, errorflag)

    if(gl%seawater .or. gl%el_present) then
        errorflag = -12800
    end if

    !Andreas March 2014
    !Routine not yet implemented for solids
    if (gl%check_solid) then
        errorflag = -12900
    end if

    ! irregular input, or problems with fluid file or reference state
    if (errorflag /= 0) then
        PEOS = errorflag
        DEOS = errorflag
        TEOS = errorflag
        UEOS = errorflag
        HEOS = errorflag
        SEOS = errorflag
        CVEOS = errorflag
        CPEOS = errorflag
        CP0EOS = errorflag
        WSEOS = errorflag
        GEOS = errorflag
        AEOS = errorflag
        BEOS = errorflag
        CEOS = errorflag
        input=gl%inptorig
        return
    end if


    !regular input

    call wm_mix_calc(gl, wm_mix)

    !calculate CALC
    !______________________________________________________________________________
    if (input == 'td') then
        t = prop1
        !if (gl%unitin == 2) then
        !d = prop2 / wm_mix
        !elseif (gl%unitin ==1) then
        d = prop2
        !end if

        !check whether pure fluid or mixture
        if (gl%ncomp == 1) then !pure fluid
            nrsubst = 1
        else  !mixture
            nrsubst = 0
        end if

        !convert specific input to molar
        if(gl%unitin==2) then
            call convert_inputprop(gl, 1, d)
        end if
        !______________________________________________________________________________
        !calculate TP
    else if (input == 'tp') then
        t = prop1
        p = prop2
        !check whether pure fluid or mixture
        if (gl%ncomp == 1) then !pure fluid
            nrsubst = 1
        else  !mixture
            nrsubst = 0
        end if
        d = rhomix_calc(gl,t, p, 0.d0, phase_ind, nrsubst)

        !______________________________________________________________________________
        !calculate PH
    else if (input == 'ph') then
        errorflag = -9955
        !______________________________________________________________________________
        !calculate PS
    else if (input == 'ps') then
        errorflag = -9955
        !______________________________________________________________________________
    else if (input == 'tliq') then
        errorflag = -9955
        !______________________________________________________________________________
    else if (input == 'tvap') then
        errorflag = -9955
        !______________________________________________________________________________
    else if (input == 'pliq') then
        errorflag = -9955
        !______________________________________________________________________________
    else if (input == 'pvap') then
        errorflag = -9955
    end if

    !errorflag =! 0
    if (errorflag /= 0) then
        TEOS = errorflag
        DEOS = errorflag
        PEOS = errorflag
        UEOS = errorflag
        HEOS = errorflag
        SEOS = errorflag
        GEOS = errorflag
        AEOS = errorflag
        CPEOS = errorflag
        CVEOS = errorflag
        WSEOS = errorflag
        BEOS = errorflag
        CEOS = errorflag
        CP0EOS = errorflag
        input=gl%inptorig
        return
    end if

    !errorflag = 0

    UEOS = U_CALC(gl,t, d, nrsubst)
    HEOS = H_CALC(gl,t, d, nrsubst)
    SEOS = S_CALC(gl,t, d, nrsubst)
    GEOS = G_CALC(gl,t, d, nrsubst)
    AEOS = A_CALC(gl,t, d, nrsubst)
    CPEOS = CP_CALC(gl,t, d, nrsubst)
    CVEOS = CV_CALC(gl,t, d, nrsubst)
    WSEOS = WS_CALC(gl,t, d, nrsubst)

    DEOS = d
    PEOS = P_CALC(gl,T, d, nrsubst)


    TEOS = t
    BEOS = B_CALC(gl,T,nrsubst)
    CEOS = C_CALC(gl,T,nrsubst)
    CP0EOS = CP0_CALC(gl,T,nrsubst)


    !If the SRK is used, check whether an ideal part exists
    !Andreas November 2013
    Do i = 1, gl%ncomp
        if(gl%Eq_type(1) == 2) then
            !---------------------------------------------------------
            !If constants B and C for the SRK cp0 model are 0, no model is implemented
            if ((dabs(gl%B_cv0(i)) < 1.D-8) .and. (dabs(gl%C_cv0(i)) < 1.D-8)) then
                !No ideal part for SRK exists
                errorflag = -2908
                UEOS = errorflag
                HEOS = errorflag
                SEOS = errorflag
                GEOS = errorflag
                AEOS = errorflag
                CPEOS = errorflag
                CVEOS = errorflag
                WSEOS = errorflag
                CP0EOS = errorflag
            end if
        end if
    end do

    if (gl%unitin ==2) then
        DEOS = DEOS * wm_mix
        UEOS = UEOS / wm_mix
        HEOS = HEOS / wm_mix
        SEOS = SEOS / wm_mix
        GEOS = GEOS / wm_mix
        AEOS = AEOS / wm_mix
        CPEOS = CPEOS / wm_mix
        CVEOS = CVEOS / wm_mix
        BEOS = BEOS / wm_mix
        CEOS =  CEOS / (wm_mix**2.d0)
        CP0EOS = CP0EOS / wm_mix
    end if



    input=gl%inptorig

    end subroutine  ALLPROP_HOM


    !*******************************************************
    ! subroutine to calculate all thermodynamic properties
    !*******************************************************
    subroutine ALLPROP_STDCALL (input, prop1, prop2, fluids, moles, EOS_indicator,  MIX_indicator, path, unitdefinition,TEOS,DEOS,PEOS,UEOS,HEOS,SEOS,GEOS,AEOS,CPEOS,CVEOS,WSEOS,BEOS,CEOS,CP0EOS,QEOS,errorflag,handle)
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "ALLPROP_STDCALL" :: ALLPROP_STDCALL
    !DEC$ ATTRIBUTES stdcall, reference :: ALLPROP_STDCALL

    use, intrinsic :: iso_c_binding

    implicit none

    !Handle for gl
    type(c_ptr), intent(inout):: handle


    !define input variables
    character (12), intent(in) :: input
    character(20) :: unitdefinition
    double precision, intent(in) :: prop1, prop2
    character (30), dimension(30), intent(in) :: fluids
    double precision, dimension(30), intent(in) :: moles
    integer, dimension(30), intent(in) :: EOS_indicator
    integer:: MIX_indicator
    character (255), intent(in) :: path

    !define output variables
    double precision, intent(out) :: TEOS,DEOS,PEOS,UEOS,HEOS,SEOS,GEOS,AEOS,CPEOS,CVEOS,WSEOS,BEOS,CEOS,CP0EOS,QEOS

    integer :: errorflag


    call ALLPROP (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator,  path, unitdefinition,TEOS,DEOS,PEOS,UEOS,HEOS,SEOS,GEOS,AEOS,CPEOS,CVEOS,WSEOS,BEOS,CEOS,CP0EOS,QEOS,errorflag,handle)
    unitdefinition  =unitdefinition
    end subroutine


    SUBROUTINE ALLPROP (input, prop1, prop2, fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition,&
        & TEOS,DEOS,PEOS,UEOS,HEOS,SEOS,GEOS,AEOS,CPEOS,CVEOS,WSEOS,BEOS,CEOS,CP0EOS,QEOS,errorflag,gl_handle)
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "ALLPROP" :: ALLPROP



    use, intrinsic :: iso_c_binding

    implicit none

    !-------------------------------------------------------------------------
    !Handle for gl
    type(type_gl), pointer :: gl
    type(c_ptr), intent(inout):: gl_handle
    !-------------------------------------------------------------------------



    !define input variables
    !double precision :: P_CALC,U_CALC,H_CALC,S_CALC,CV_CALC,CP_CALC,CP0_CALC,WS_CALC,G_CALC,A_CALC,B_CALC,C_CALC
    DOUBLE PRECISION, intent(out) :: TEOS,DEOS,PEOS,UEOS,HEOS,SEOS,CVEOS,CPEOS,CP0EOS,WSEOS,GEOS,AEOS,BEOS,CEOS
    DOUBLE PRECISION, intent(out) :: QEOS
    double precision::u_phase,h_phase,s_phase,A_phase,g_phase,cp_phase,cv_phase,ws_phase,q_phase

    character (12) :: input
    character(20) :: unitdefinition
    double precision :: prop1, prop2
    character (30), dimension(30) :: fluids
    double precision, dimension(30) :: moles
    character (255) :: path
    integer, dimension(30):: EOS_indicator
    integer:: MIX_indicator

    double precision:: rho(5)           ! vector containing the densities of all phases
    double precision:: x_Phase(30,5)    ! vector conaining the compositions of all phases
    integer:: nrofphases                ! at this point: only 1 or 2 implemented
    !Indicate which phases are present
    integer:: phasetype(5)              !phasetype contains the phase indicator number
    integer :: phase
    integer:: iter, iFlash
    double precision :: x_liq(30)
    double precision :: x_vap(30)

    double precision, dimension(30) :: molesallpropin

    !define variables for calculating
    double precision :: t, d, p, h, s
    double precision :: dvap
    double precision :: dliq, wm_mix, wm_mix_all
    double precision :: vapfrac, d_phase, dens_phase, dens
    double precision, dimension(5) :: phasefrac, wm_phase


    integer :: nrsubst
    integer :: unitin
    integer :: errorflag
    integer :: ILIMITS
    integer :: i

    if(trim(unitdefinition) == "specific") then
        unitin=2
    elseif(trim(unitdefinition) == "molar") then
        unitin=1
    elseif(trim(unitdefinition) == "reduced") then
        unitin=3
    else
        errorflag=-123456
        TEOS = errorflag
        DEOS = errorflag
        PEOS = errorflag
        UEOS = errorflag
        HEOS = errorflag
        SEOS = errorflag
        GEOS = errorflag
        AEOS = errorflag
        CPEOS = errorflag
        CVEOS = errorflag
        WSEOS = errorflag
        BEOS = errorflag
        CEOS = errorflag
        CP0EOS = errorflag
        QEOS = errorflag
        prop1 = prop1
        prop2 = prop2
        unitdefinition = unitdefinition
        input = input
        return
    end if

    !control handle and fluid input fot memory usage
    call control_fluids(gl,input,fluids,moles,EOS_indicator,unitin,gl_handle,errorflag)

    if(errorflag /= 0) return


    dens = 0
    TEOS = 0
    DEOS = 0
    PEOS = 0
    UEOS = 0
    HEOS = 0
    SEOS = 0
    GEOS = 0
    AEOS = 0
    CPEOS = 0
    CVEOS = 0
    WSEOS = 0
    BEOS = 0
    CEOS = 0
    CP0EOS = 0
    QEOS = 0

    ILIMITS = 0

    errorflag = 0

    gl%calc_ref = .true.
    gl%transport = 0       ! needed to check the range of validity for transport properties; 0 for every other property

    call uppertolower_char(input, len(input))

    call uppertolower_char(unitdefinition, len(unitdefinition)) !molar <-> specific handling Benedikt May 2018


    !call setup - >  initializing, reading fluid files
    call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator,errorflag)

    if((gl%seawater) .or. (gl%el_present)) then
        errorflag = -12800
    end if

    molesallpropin = moles

    call wm_mix_calc(gl,wm_mix)

    !Andreas March 2014
    !Routine not yet implemented for solids
    if (gl%check_solid) then
        errorflag = -12900
    end if

    ! irregular input, or problems with fluid file or reference state
    if (errorflag /= 0) then
        TEOS = errorflag
        DEOS = errorflag
        PEOS = errorflag
        UEOS = errorflag
        HEOS = errorflag
        SEOS = errorflag
        GEOS = errorflag
        AEOS = errorflag
        CPEOS = errorflag
        CVEOS = errorflag
        WSEOS = errorflag
        BEOS = errorflag
        CEOS = errorflag
        CP0EOS = errorflag
        QEOS = errorflag
        input=gl%inptorig
        return
    end if


    !If the SRK, PR or LKP is used, check whether a model for the ideal gas exists
    Do i = 1, gl%ncomp
        if(gl%Eq_type(i) > 1) then
            !---------------------------------------------------------
            !If constants B and C for the SRK cp0 model are 0, no model is implemented
            if (((dabs(gl%B_cv0(i)) < 1.D-8) .and. (dabs(gl%C_cv0(i)) < 1.D-8))) then
                errorflag = -2908
                TEOS = errorflag
                DEOS = errorflag
                PEOS = errorflag
                UEOS = errorflag
                HEOS = errorflag
                SEOS = errorflag
                GEOS = errorflag
                AEOS = errorflag
                CPEOS = errorflag
                CVEOS = errorflag
                WSEOS = errorflag
                BEOS = errorflag
                CEOS = errorflag
                CP0EOS = errorflag
                QEOS = errorflag
                input=gl%inptorig
                return
            end if
        end if
    end do

    !Handle different input parameters and do all calculations
    call inpt_handle (gl,input, prop1, prop2, p, t, d, dvap, dliq, moles, rho, x_Phase, phasetype, &
        & phasefrac, nrofphases, nrsubst, errorflag)


    !errorflag =! 0
    if (errorflag /= 0) then
        TEOS = errorflag
        DEOS = errorflag
        PEOS = errorflag
        UEOS = errorflag
        HEOS = errorflag
        SEOS = errorflag
        GEOS = errorflag
        AEOS = errorflag
        CPEOS = errorflag
        CVEOS = errorflag
        WSEOS = errorflag
        BEOS = errorflag
        CEOS = errorflag
        CP0EOS = errorflag
        QEOS = errorflag
        input=gl%inptorig
        return
    end if

    TEOS = t
    PEOS = p
    call reduced_parameters_calc(gl,t)
    BEOS = B_calc(gl,t, nrsubst)
    CEOS = C_calc(gl,t, nrsubst)
    CP0EOS = CP0_calc(gl,t, nrsubst)

    !Start of property calculation
    ! Fluid phase

    do i = 1, nrofphases
        if (phasetype(i) < 4) then

            if (gl%ncomp > 1) then !mixture
                !If mixture, recalculate tredmix and rhoredmix
                gl%molfractions = x_phase(:,phasetype(i))


                call reduced_parameters_calc(gl,t) !Dummy temperature 300 K for the SRK
            end if


            d_phase = rho(phasetype(i)) ! properties will be calculated with this phase specific density
            dens_phase = 1.d0/d_phase   ! needed for the summation of the overall density


            u_phase = U_calc(gl,t, d_phase, nrsubst)
            h_phase = H_calc(gl,t, d_phase, nrsubst)
            s_phase = S_calc(gl,t, d_phase, nrsubst)
            g_phase = G_calc(gl,t, d_phase, nrsubst)
            a_phase = A_calc(gl,t, d_phase, nrsubst)
            cp_phase = CP_calc(gl,t, d_phase, nrsubst)
            cv_phase = CV_calc(gl,t, d_phase, nrsubst)
            ws_phase = WS_calc(gl,t, d_phase, nrsubst)
            q_phase = phasefrac(phasetype(i))


            !Density might be needed afterwards for further calculations
            dens = dens + phasefrac(phasetype(i)) * dens_phase      ! CHECK it this leads to proper results (input / output td)

            UEOS = UEOS + phasefrac(phasetype(i)) * u_phase
            hEOS = hEOS + phasefrac(phasetype(i)) * h_phase
            sEOS = sEOS + phasefrac(phasetype(i)) * s_phase
            gEOS = gEOS + phasefrac(phasetype(i)) * g_phase
            aEOS = aEOS + phasefrac(phasetype(i)) * a_phase

            if ((nrofphases > 1) .and. (.not. gl%phaseboundary)) then
                CPEOS = -6666.d0
                CVEOS = -6666.d0
                WSEOS = -6666.d0
            else
                CPEOS = CPEOS + phasefrac(phasetype(i)) * cp_phase
                CVEOS = CVEOS + phasefrac(phasetype(i)) * cv_phase
                WSEOS = WSEOS + phasefrac(phasetype(i)) * ws_phase
            end if

            QEOS = QEOS + phasefrac(phasetype(i)) * q_phase

        else
            errorflag = -12900
            TEOS = errorflag
            DEOS = errorflag
            PEOS = errorflag
            UEOS = errorflag
            HEOS = errorflag
            SEOS = errorflag
            GEOS = errorflag
            AEOS = errorflag
            CPEOS = errorflag
            CVEOS = errorflag
            WSEOS = errorflag
            BEOS = errorflag
            CEOS = errorflag
            CP0EOS = errorflag
            QEOS = errorflag
            input=gl%inptorig
            return
        end if
    end do


    dens = 1.d0/dens
    DEOS = dens

    ! be aware of density

    !check limits
    call check_limits (gl,input, t, p, d, dliq, dvap, nrofphases, x_phase, phasetype, ILIMITS)


    if (ILIMITS /= 0) then
        TEOS = ILIMITS
        DEOS = ILIMITS
        PEOS = ILIMITS
        UEOS = ILIMITS
        HEOS = ILIMITS
        SEOS = ILIMITS
        GEOS = ILIMITS
        AEOS = ILIMITS
        CPEOS = ILIMITS
        CVEOS = ILIMITS
        WSEOS = ILIMITS
        BEOS = ILIMITS
        CEOS = ILIMITS
        CP0EOS = ILIMITS
        QEOS = ILIMITS
        input=gl%inptorig
        return
    end if

    if (gl%unitin == 2) then
        gl% molfractions = moles
        call wm_mix_calc(gl, wm_mix_all)

        DEOS = DEOS * wm_mix_all
        UEOS = UEOS / wm_mix_all
        HEOS = HEOS / wm_mix_all
        SEOS = SEOS / wm_mix_all
        GEOS = GEOS / wm_mix_all
        AEOS = AEOS / wm_mix_all
        CPEOS = CPEOS / wm_mix_all
        CVEOS = CVEOS / wm_mix_all
        BEOS = BEOS / wm_mix_all
        CEOS =  CEOS / (wm_mix_all**2.d0)
        CP0EOS = CP0EOS / wm_mix_all

    end if



    input=gl%inptorig
    unitdefinition = unitdefinition

    end subroutine ALLPROP


    !!DEC$ IF DEFINED(CON_FIT)
    !!####################################################################################################################
    !!####################################################################################################################
    !!Set the parameter in gl from extern arrays
    !subroutine  set_parameter_from_extern(gl,x,n_i,t_i,d_i,p_i,gauss,add_vals,nterms,nrsubst,fit_mode,fit_red)
    !! Variables
    !use module_all_types
    !implicit none
    !type(type_gl):: gl
    !double precision, dimension(:) :: x
    !double precision, dimension(2):: add_vals
    !integer, dimension(:):: nterms,d_i,p_i
    !double precision,dimension(:):: n_i,t_i
    !double precision,dimension(:,:)::gauss
    !integer:: i,nrsubst,nterm,fit_mode,fit_red,n_red_param
    !double precision, allocatable, dimension(:):: d_tmp
    !integer*4 :: Int_4
    !if(fit_red .eq. 1) then
    !    n_red_param = 2
    !else
    !    n_red_param = 0
    !end if
    !
    !! Body of set_parameter_from_extern
    !gl%EQ_type = 0
    !gl%components(1) = 'ljf'
    !gl%EQ_type(1) =  1 !has to be changed later!!!
    !gl%mix_Type = 1
    !gl%Factor = 1.D0
    !gl%factorpress = 1.D6
    !gl%factortrans=1.D0
    !gl%factorrbwr=1.D0
    !gl%factor_VS5eta=1.D0
    !gl%factor_VS5slj=1.D0
    !gl%ncomp = 1
    !gl%REQ(1)  = 1
    !gl%ttp(1) = 0.66d0
    !gl%tred(1) = add_vals(1)
    !gl%rhored(1) = add_vals(2)
    !gl%tc(1) = gl%tred(1)
    !gl%accen = 0d0
    !gl%rhoredmix = gl%rhored(1)
    !gl%tredmix =  gl%tred(1)
    !if(.not. allocated(gl%eos_coeff)) allocate(gl%eos_coeff)
    !gl%eos_coeff%I_pol(1) = nterms(1)
    !gl%eos_coeff%I_exp(1) = nterms(2)
    !gl%eos_coeff%I_gbs(1) = nterms(3)
    !gl%eos_coeff%nreg(1) =  sum(nterms(1:2))
    !nterm  = sum(nterms(1:3))
    !
    !if(.not.allocated(gl%eos_coeff%ni)) then
    !    allocate(gl%eos_coeff%ni(nterm,1))
    !    allocate(gl%eos_coeff%ti(nterm,1))
    !    allocate(gl%eos_coeff%di(nterm,1))
    !    allocate(gl%eos_coeff%p_i(nterm,1))
    !    allocate(gl%eos_coeff%gama(nterm,1))
    !    allocate(gl%eos_coeff%pli(gl%eos_coeff%I_gbs(1),1))
    !    allocate(gl%eos_coeff%tli(gl%eos_coeff%I_gbs(1),1))
    !    allocate(gl%eos_coeff%eta(gl%eos_coeff%I_gbs(1),1))
    !    allocate(gl%eos_coeff%beta(gl%eos_coeff%I_gbs(1),1))
    !    allocate(gl%eos_coeff%gam(gl%eos_coeff%I_gbs(1),1))
    !    allocate(gl%eos_coeff%eps(gl%eos_coeff%I_gbs(1),1))
    !end if
    !
    !
    !allocate(d_tmp(nterm))
    !
    !
    !
    !! density exponent fitting (doesnt make sense, maybe it makes sense )
    !if(fit_mode .eq. 1) then
    !    gl%eos_coeff%ni(1:nterm,1) = x(1:nterm)
    !    gl%eos_coeff%ti(1:nterm,1) = x(1+nterm:2*nterm)
    !    gl%eos_coeff%di(1:nterm,1) = x(2*nterm+4*gl%eos_coeff%I_gbs(1)+1:3*nterm+4*gl%eos_coeff%I_gbs(1))!d_i(1:nterm) !at the end of the vector
    !    d_tmp(1:nterm) =  gl%eos_coeff%di(1:nterm,1)
    !    if(gl%eos_coeff%I_gbs(1) .gt. 0) then
    !        gl%eos_coeff%pli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
    !        gl%eos_coeff%tli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
    !        gl%eos_coeff%eta(1:gl%eos_coeff%I_gbs(1),1) =  x(3*nterm+1:3*nterm+gl%eos_coeff%I_gbs(1))
    !        gl%eos_coeff%beta(1:gl%eos_coeff%I_gbs(1),1) = x(3*nterm+1+gl%eos_coeff%I_gbs(1):3*nterm+2*gl%eos_coeff%I_gbs(1))
    !        gl%eos_coeff%gam(1:gl%eos_coeff%I_gbs(1),1) =  x(3*nterm+1+2*gl%eos_coeff%I_gbs(1):3*nterm+3*gl%eos_coeff%I_gbs(1))
    !        gl%eos_coeff%eps(1:gl%eos_coeff%I_gbs(1),1) =  x(3*nterm+1+3*gl%eos_coeff%I_gbs(1):3*nterm+4*gl%eos_coeff%I_gbs(1))
    !    end if
    !elseif(fit_mode == 2) then ! all are fitted
    !    gl%eos_coeff%ni(1:nterm,1) = x(1:nterm)
    !    gl%eos_coeff%ti(1:nterm,1) = x(1+nterm:2*nterm)
    !    if(gl%eos_coeff%I_gbs(1) .gt. 0) then
    !        gl%eos_coeff%pli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
    !        gl%eos_coeff%tli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
    !        gl%eos_coeff%eta(1:gl%eos_coeff%I_gbs(1),1) =  x(2*nterm+1:2*nterm+gl%eos_coeff%I_gbs(1))
    !        gl%eos_coeff%beta(1:gl%eos_coeff%I_gbs(1),1) = x(2*nterm+1+gl%eos_coeff%I_gbs(1):2*nterm+2*gl%eos_coeff%I_gbs(1))
    !        gl%eos_coeff%gam(1:gl%eos_coeff%I_gbs(1),1) =  x(2*nterm+1+2*gl%eos_coeff%I_gbs(1):2*nterm+3*gl%eos_coeff%I_gbs(1))
    !        gl%eos_coeff%eps(1:gl%eos_coeff%I_gbs(1),1) =  x(2*nterm+1+3*gl%eos_coeff%I_gbs(1):2*nterm+4*gl%eos_coeff%I_gbs(1))
    !    end if
    !    gl%eos_coeff%di(1:nterm,1) = d_i(1:nterm)
    !elseif(fit_mode .eq. 3) then !only n are fitted
    !    gl%eos_coeff%ni(1:nterm,1) = x(1:nterm)
    !    gl%eos_coeff%ti(1:nterm,1) = t_i
    !    gl%eos_coeff%pli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
    !    gl%eos_coeff%tli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
    !    gl%eos_coeff%eta(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,1)
    !    gl%eos_coeff%beta(1:gl%eos_coeff%I_gbs(1),1) = gauss(:,2)
    !    gl%eos_coeff%gam(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,3)
    !    gl%eos_coeff%eps(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,4)
    !    gl%eos_coeff%di(1:nterm,1) = d_i(1:nterm)
    !end if
    !
    !gl%eos_coeff%p_i(:,1) =  0
    !gl%eos_coeff%gama = 0d0
    !
    !do i = gl%eos_coeff%I_pol(1)+1,gl%eos_coeff%nreg(1)
    !    gl%eos_coeff%p_i(i,1) =  p_i(i)
    !    gl%eos_coeff%gama(i,1) = 1.d0
    !end do
    !
    !!do i=1,gl%eos_coeff%I_exp(1)
    !!    gl%eos_coeff%gama(gl%eos_coeff%I_pol(1)+i,1) = x(size(x,1)-gl%eos_coeff%I_exp(1)-n_red_param+i)
    !!    !write(*,*) gl%eos_coeff%gama(gl%eos_coeff%I_pol(1)+i,1)
    !!end do
    !
    !end subroutine set_parameter_from_extern
    !!####################################################################################################################
    !!####################################################################################################################

    !DEC$ IF DEFINED(IDEAL_CURVES)

    !subroutine ideal_curves(x,d_i,p_i,nterms,add_vals,sizes)
    integer function ideal_cur(path_in,press,temps,resolution)
    use ideal_curves
    !use iso_c_bindings
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "ideal_cur" :: ideal_cur
    implicit none
    !integer, dimension(3):: sizes
    !integer:: mode ! 1: ssq calc, 2: single prop
    !double precision, dimension(sizes(1)) :: x
    !double precision, dimension(2):: add_vals
    !integer, dimension(sizes(2)) :: nterms
    !integer, dimension(sum(nterms)) :: d_i,p_i
    character(255):: path_in
    integer:: resolution,i
    double precision, dimension(resolution,4):: press,temps
    type(type_gl),pointer:: gl
    type(ideal):: id
    character(30), dimension(30):: fluids_in
    double precision, dimension(30):: moles_in
    integer, dimension(30):: EOS_INDICATOR_in
    character(20):: unit_in
    character(12):: input_in
    type(c_ptr):: handle
    integer:: nrsubst
    ! Variables

    !set all necessary variables for the trend calculations
    !call set_parameter_from_extern(gl,x,d_i,p_i,add_vals,nterms,nrsubst)
    handle = c_null_ptr
    fluids_in =''
    fluids_in(1) = 'MIE12'
    unit_in = 'molar'
    input_in  = 'TD'
    nrsubst = 1
    moles_in = 0d0
    moles_in(1) = 1d0
    EOS_INDICATOR_in = 0
    EOS_INDICATOR_in(1) = 1
    id  = ideal(gl,input_in,fluids_in,moles_in,EOS_indicator_in,path_in,unit_in,handle,nrsubst)

    plt_i%resolution = 3
    call id%get_classic(gl,nrsubst)

    do i=1,4
        press(:,i)  = id%crv(i)%p_red(:)
        temps(:,i)  = id%crv(i)%t_red(:)
    end do

    ideal_cur = 0
    ! Body of ideal_curves
    !init_ideal(gl,input_in,fluids_in,moles_in,EOS_indicator_in,path_in,unit_in,gl_handle_in,nrsubst)
    end function ideal_cur
    !####################################################################################################################
    !####################################################################################################################
    !DEC$ END IF

    !
    !!   prop    property
    !!    1       T
    !!    2       D
    !!    3       P
    !!    4       U
    !!    5       H
    !!    6       S
    !!    7       G
    !!    8       A
    !!    9       CP
    !!   10       CV
    !!   11       WS
    !!   12       B
    !!   13       C
    !!   14       CP0
    !!   15       Q
    !!   16       Z
    !!   17       PNUMER
    !!   18       WSNUMER
    !!   19       A00
    !!   20       A01
    !!   21       A02
    !!   22       A03
    !!   23       A10
    !!   24       A11
    !!   25       A12
    !!   26       A20
    !!   27       A21
    !!   28       A30
    !!   29       UR
    !!   30       HR
    !!   31       SR
    !!   32       CVR
    !!   33       CPOTR
    !!   34       GRUEN
    !!   35       PIP
    !!   36       DE
    !!   37       ST
    !!   38       VISDYN
    !!   39       TCX
    !!   40       VOLEXP
    !!   41       DVIR
    !!   42       RIEM
    !!   43       COMPT
    !!   44       COMPS
    !!   45       THROT
    !!   46       EXPANS
    !!   47       EXPANT
    !!   48       DUDV
    !!   49       JTCO
    !!   50       DPDT
    !!   51       VE
    !!   52       HE
    !!   53       GE
    !!   54       B12
    !!   55       GAMMAGD
    !!   56       Residual pressure
    !!   57       dBdT (derivative of the second virial coeffcient)
    !!   58       dCdT (derivative of the third virial coeffcient)
    !!   59       dPdD (1st derivative of the pressure with respect to density at constant temperature)
    !!   60       d2PdD2 (2nd derivative of the pressure with respect to density at constant temperature)
    !!   61       d2PdDdT (1st mixed derivatives of the pressure with respect to density and temperature)
    !!   62       d2PdT2 (2nd derivative of the pressure with respect to temperature at constant density)
    !!   63       dDdT (1st derivative of the density with respect to temperature at constant pressure)
    !!   64       DSPIN (Calculation of Spinodals)
    !!   65       PHASEEOS (Phase Determination)
    !!   66       AF Acentric factor
    !!   67       Density maximum of water and heavy water
    !!   68       Boyle-Temperature
    !!   69       Joule-Thomson-Temperature
    !!   70       Joule-Thomson-Inversion-Temperature
    !!   71       First PIP derivative with respect to density (at constant temperature) (numerical derivative)
    !!   72       Second PIP derivative with respect to density^2 (at constant temperature) (numerical derivative)
    !!   73       Second derivative of second virial coeff. B with respect to T^2
    !!   74       Second derivative of third virial coeff. C with respect to T^2
    !!   75       Second derivative of fourth virial coeff. D with respect to T^2
    !!   76       Third derivative of second virial coeff. B with respect to T^3
    !!   77       Third derivative of third virial coeff. C with respect to T^3
    !!   78       Third derivative of fourth virial coeff. D with respect to T^3
    !!   79       Fourth derivative of second virial coeff. B with respect to T^4
    !!   80       Fourth derivative of third virial coeff. C with respect to T^4
    !!   81       Fourth derivative of fourth virial coeff. D with respect to T^4
    !!#####################################################################################################################################################
    !! INTERFACES ROUTINES FOR THE PARALLEL FITTER (Sven P.)
    !!#####################################################################################################################################################
    !double precision function trend_ssq(x,n_i,t_i,d_i,p_i,gauss,nterms,add_vals,datas,datas_id,sizes,mode,input,fit_mode,fit_red)
    !use module_all_types
    !implicit none
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "trend_ssq" :: trend_ssq
    !! Variables
    !
    !integer, dimension(3):: sizes
    !integer:: mode,fit_mode ! 1: ssq calc, 2: single prop
    !double precision, dimension(sizes(1)) :: x
    !double precision, dimension(2):: add_vals
    !integer, dimension(sizes(2)) :: nterms
    !double  precision, dimension(nterms(3),4)::gauss
    !double precision, dimension(sizes(3),5):: datas
    !integer, dimension(sizes(3)):: datas_id
    !character(2):: input
    !double precision, dimension(sizes(3)):: val_calc
    !integer, dimension(sum(nterms)) :: d_i,p_i
    !double precision, dimension(sum(nterms))::n_i,t_i !these are necessary if the parameters are not fitted
    !! Body of trend_ssq
    !integer:: nterm,i,errval, iter, nrsubst,IPHASE,fit_red
    !double  precision:: T,D,rho_vap,rho_liq,press,RHO_EST_GIVEN
    !logical:: num_virial_mode ! if = true: numerical approximation of virial coeffs because the d_i are no integers (fit process)
    !type(c_ptr):: handle
    !type(type_gl):: gl
    !double precision:: test1,test2
    !
    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !!init phase
    !handle = c_null_ptr
    !call  CREATE_FLUID(handle)
    !CALL C_F_POINTER(handle,gl)
    !nrsubst = 1
    !call set_parameter_from_extern(gl,x,n_i,t_i,d_i,p_i,gauss,add_vals,nterms,nrsubst,fit_mode,fit_red)
    !gl%pc(1) = p_calc(gl,gl%tred(1),gl%rhored(1),1)
    !trend_ssq = 0d0
    !nterm = sum(nterms,1)
    !!all numbers are real integers
    !if(all(int(gl%eos_coeff%di(1:nterm,1)) .eq. gl%eos_coeff%di(1:nterm,1))) then
    !    num_virial_mode = .false.
    !else
    !    num_virial_mode = .true.
    !end if
    !
    !
    !
    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !!sum of squares calculation
    !if(mode.eq.1) then
    !    !!$omp  parallel do
    !    do i=1,sizes(3)
    !        if( datas_id(i) == 19) then
    !            val_calc(i) = (datas(i,3)-A00_calc(gl,datas(i,1),datas(i,2),nrsubst))/datas(i,4)*datas(i,5)
    !        elseif(datas_id(i) == 20) then
    !            val_calc(i) = (datas(i,3)-A01_calc(gl,datas(i,1),datas(i,2),nrsubst))/datas(i,4)*datas(i,5)
    !        elseif(datas_id(i) == 21) then
    !            val_calc(i) = (datas(i,3)-A02_calc(gl,datas(i,1),datas(i,2),nrsubst))/datas(i,4)*datas(i,5)
    !        elseif(datas_id(i) == 22) then
    !            val_calc(i) = (datas(i,3)-A03_calc(gl,datas(i,1),datas(i,2),nrsubst))/datas(i,4)*datas(i,5)
    !        elseif( datas_id(i) == 23) then
    !            val_calc(i) = (datas(i,3)-A10_calc(gl,datas(i,1),datas(i,2),nrsubst))/datas(i,4)*datas(i,5)
    !        elseif( datas_id(i) == 24) then
    !            val_calc(i) = (datas(i,3)-A11_calc(gl,datas(i,1),datas(i,2),nrsubst))/datas(i,4)*datas(i,5)
    !        elseif( datas_id(i) == 25) then
    !            val_calc(i) = (datas(i,3)-A11_calc(gl,datas(i,1),datas(i,2),nrsubst))/datas(i,4)*datas(i,5)
    !        elseif( datas_id(i) == 26) then
    !            val_calc(i) = (datas(i,3)-A20_calc(gl,datas(i,1),datas(i,2),nrsubst))/datas(i,4)*datas(i,5)
    !        elseif( datas_id(i) == 27) then
    !            val_calc(i) = (datas(i,3)-A21_calc(gl,datas(i,1),datas(i,2),nrsubst))/datas(i,4)*datas(i,5)
    !        elseif( datas_id(i) == 28) then
    !            val_calc(i) = (datas(i,3)-A30_calc(gl,datas(i,1),datas(i,2),nrsubst))/datas(i,4)*datas(i,5)
    !        elseif( datas_id(i) == 12) then
    !            if(.not.num_virial_mode) then
    !                val_calc(i) = (datas(i,3)-B_CALC(gl,datas(i,1),nrsubst))/datas(i,4)*datas(i,5)
    !            else
    !                val_calc(i) = (datas(i,3)-B_CALC_num(gl,datas(i,1),nrsubst))/datas(i,4)*datas(i,5)
    !            end if
    !        elseif( datas_id(i) == 13) then
    !            if(.not.num_virial_mode) then
    !                val_calc(i) = (datas(i,3)-C_CALC(gl,datas(i,1),nrsubst))/datas(i,4)*datas(i,5)
    !            else
    !                val_calc(i) = (datas(i,3)-C_CALC_num(gl,datas(i,1),nrsubst))/datas(i,4)*datas(i,5)
    !            end if
    !
    !        elseif( datas_id(i) == 41) then
    !            if(.not.num_virial_mode) then
    !                val_calc(i) = (datas(i,3)-D_CALC(gl,datas(i,1),nrsubst))/datas(i,4)*datas(i,5)
    !            else
    !                val_calc(i) = (datas(i,3)-D_CALC_num(gl,datas(i,1),nrsubst))/datas(i,4)*datas(i,5)
    !            end if
    !        end if
    !
    !
    !
    !    end do
    !    !!$omp end parallel do
    !    trend_ssq = sum(val_calc**2)
    !
    !    !single calc of a property (like TREND_CALC)
    !elseif(mode.eq.2)then
    !    trend_ssq = 0d0
    !    if( datas_id(1) == 19) then
    !        trend_ssq = A00_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif(datas_id(1) == 20) then
    !        trend_ssq = A01_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif(datas_id(1) == 21) then
    !        trend_ssq = A02_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif(datas_id(1) == 22) then
    !        trend_ssq = A03_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif( datas_id(1) == 23) then
    !        trend_ssq = A10_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif( datas_id(1) == 24) then
    !        trend_ssq = A11_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif( datas_id(1) == 25) then
    !        trend_ssq = A11_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif( datas_id(1) == 26) then
    !        trend_ssq = A20_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif( datas_id(1) == 27) then
    !        trend_ssq = A21_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif( datas_id(1) == 28) then
    !        trend_ssq = A30_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif(datas_id(1) == 59) then
    !        trend_ssq = DPDD_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif(datas_id(1) == 60) then
    !        trend_ssq = D2PDD2_CALC(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif(datas_id(1) == 35) then
    !        trend_ssq = PIP_CALC(gl,datas(1,1),datas(1,2),nrsubst)
    !    elseif(datas_id(1) == 12) then
    !        trend_ssq = B_CALC(gl,datas(1,1),nrsubst)
    !    elseif(datas_id(1) == 13) then
    !        trend_ssq = C_CALC(gl,datas(1,1),nrsubst)
    !    elseif(datas_id(1) == 41) then
    !        trend_ssq = D_CALC(gl,datas(1,1),nrsubst)
    !    elseif(datas_id(1) == 58) then
    !        trend_ssq = DCDT_CALC(gl,datas(1,1),nrsubst)
    !    elseif(datas_id(1) == 74) then
    !        trend_ssq = D2CDT2_CALC(gl,datas(1,1),nrsubst)
    !    elseif(datas_id(1) == 75) then
    !        !trend_ssq = D_GRUEN_DT_CALC(gl,datas(1,1),datas(1,2),nrsubst)
    !    end if
    !
    !
    !
    !end if
    !
    !
    !call DESTROY_FLUID(handle)
    !
    !end function trend_ssq
    !
    !!#####################################################################################################################################################
    !
    !!intern trend_Sqq
    !!#####################################################################################################################################################
    !double precision function trendvals(x,d_i,p_i,nterms,add_vals,datas,datas_id,sizes,fit_mode,fit_red,val_calc)
    !use module_all_types
    !implicit none
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "TRENDVAL" :: trendvals
    !! Variables
    !
    !integer, dimension(3):: sizes
    !integer:: mode,fit_mode ! 1: ssq calc, 2: single prop
    !double precision, dimension(sizes(1)) :: x
    !double precision, dimension(2):: add_vals
    !integer, dimension(sizes(2)) :: nterms
    !double  precision, dimension(nterms(3),4)::gauss
    !double precision, dimension(sizes(3),5):: datas
    !integer, dimension(sizes(3)):: datas_id
    !character(2):: input
    !double precision, dimension(sizes(3)):: val_calc
    !integer, dimension(sum(nterms)) :: d_i,p_i
    !double precision, dimension(sum(nterms))::n_i,t_i
    !! Body of trend_ssq
    !integer:: nterm,i,errval, iter, nrsubst,IPHASE,fit_red
    !double  precision:: T,D,rho_vap,rho_liq,press,RHO_EST_GIVEN
    !logical:: num_virial_mode ! if = true: numerical approximation of virial coeffs because the d_i are no integers (fit process)
    !type(c_ptr):: handle
    !type(type_gl):: gl
    !double precision:: test1,test2
    !
    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !!init phase
    !handle = c_null_ptr
    !call  CREATE_FLUID(handle)
    !CALL C_F_POINTER(handle,gl)
    !nrsubst = 1
    !call set_parameter_from_extern(gl,x,n_i,t_i,d_i,p_i,gauss,add_vals,nterms,nrsubst,fit_mode,fit_red)
    !gl%pc(1) = p_calc(gl,gl%tred(1),gl%rhored(1),1)
    !trendvals = 0d0
    !nterm = sum(nterms,1)
    !!all numbers are real integers
    !if(all(int(gl%eos_coeff%di(1:nterm,1)) .eq. gl%eos_coeff%di(1:nterm,1))) then
    !    num_virial_mode = .false.
    !else
    !    num_virial_mode = .true.
    !end if
    !
    !do i=1,sizes(3)
    !    if( datas_id(i) == 19) then
    !        val_calc(i) = A00_calc(gl,datas(i,1),datas(i,2),nrsubst)
    !    elseif( datas_id(i) == 23) then
    !        val_calc(i) = A10_calc(gl,datas(i,1),datas(i,2),nrsubst)
    !    elseif(datas_id(i) == 20) then
    !        val_calc(i) = A01_calc(gl,datas(i,1),datas(i,2),nrsubst)
    !    elseif( datas_id(i) == 26) then
    !        val_calc(i) = A20_calc(gl,datas(i,1),datas(i,2),nrsubst)
    !    elseif( datas_id(i) == 12) then
    !        val_calc(i) = B_CALC(gl,datas(i,1),nrsubst)
    !    elseif( datas_id(i) == 13) then
    !        val_calc(i) = C_CALC(gl,datas(i,1),nrsubst)
    !    elseif( datas_id(i) == 41) then
    !        val_calc(i) = D_CALC(gl,datas(i,1),nrsubst)
    !    elseif(datas_id(1) == 59) then
    !        val_calc(i) = DPDD_calc(gl,datas(i,1),datas(i,2),nrsubst)
    !    elseif(datas_id(1) == 60) then
    !        val_calc(i) = D2PDD2_CALC(gl,datas(i,1),datas(i,2),nrsubst)
    !    end if
    !end do
    !!!$omp end parallel do
    !trendvals = 1234d0
    !
    !call DESTROY_FLUID(handle)
    !
    !end function trendvals
    !
    !! end intern trend_Sqq
    !
    !
    !! calculation of real properties .e.g. argon
    !double precision function trendreals(sizes,x,n,t,d,p,gauss,datas,datas_id,nterms,add_vals)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "trendreals" :: trendreals
    !integer, dimension(4) :: sizes
    !integer, dimension(3) :: nterms
    !double precision, dimension(2) :: add_vals
    !double precision, dimension(sizes(1)) :: x
    !double precision, dimension(sizes(2)) :: n,t,d,p
    !double precision, dimension(sizes(3),4) :: gauss
    !double precision, dimension(sizes(4),4) :: datas
    !integer, dimension(sizes(4)) :: datas_id
    !double precision, allocatable,dimension(:) :: val_calc
    !integer :: i
    !double precision :: val
    !type(c_ptr):: handle
    !type(type_gl) :: gl
    !double precision :: press, Temp
    !double precision :: rhovap_est, rholiq_est
    !integer :: iFlash, nrsubst,IPHASE
    !integer :: errval, iter
    !
    !handle = c_null_ptr
    !call  CREATE_FLUID(handle)
    !CALL C_F_POINTER(handle,gl)
    !trendreals = 0d0
    !! Body of set_parameter_from_extern
    !gl%EQ_type = 0
    !gl%components(1) = 'argon'
    !gl%EQ_type(1) =  1 
    !gl%mix_Type = 1                         
    !gl%Factor = 1.D3            
    !gl%factorpress=1.D0         
    !gl%factortrans=1.D6         
    !gl%factorrbwr=1.D2          
    !gl%factor_VS5eta=1.D5       
    !gl%factor_VS5slj=1.D1       
    !gl%ncomp = 1
    !gl%REQ(1)  = 8.314472d0
    !gl%ttp(1) = 0.66d0
    !gl%tred(1) = add_vals(1)
    !gl%rhored(1) = add_vals(2)
    !gl%tc(1) = gl%tred(1)
    !gl%rhoc(1) = gl%rhored(1)
    !gl%accen =0d0
    !gl%rhoredmix = gl%rhored(1)
    !gl%tredmix =  gl%tred(1)
    !nrsubst = 1
    !gl%molfractions = 0d0
    !gl%molfractions(1) = 1d0
    !gl%rhomaxfluid = 0d0
    !gl%rhomaxfluid(1) = 1d6
    !allocate(val_calc(sizes(4)))
    !
    !if(.not. allocated(gl%eos_coeff)) allocate(gl%eos_coeff)
    !gl%eos_coeff%I_pol(1) = nterms(1)
    !gl%eos_coeff%I_exp(1) = nterms(2)
    !gl%eos_coeff%I_gbs(1) = nterms(3)
    !gl%eos_coeff%nreg(1) =  sum(nterms(1:2))
    !nterm  = sum(nterms(1:3))
    !
    !if(.not.allocated(gl%eos_coeff%ni)) then
    !    allocate(gl%eos_coeff%ni(nterm,1))
    !    allocate(gl%eos_coeff%ti(nterm,1))
    !    allocate(gl%eos_coeff%di(nterm,1))
    !    allocate(gl%eos_coeff%p_i(nterm,1))
    !    allocate(gl%eos_coeff%gama(nterm,1))
    !    allocate(gl%eos_coeff%pli(gl%eos_coeff%I_gbs(1),1))
    !    allocate(gl%eos_coeff%tli(gl%eos_coeff%I_gbs(1),1))
    !    allocate(gl%eos_coeff%eta(gl%eos_coeff%I_gbs(1),1))
    !    allocate(gl%eos_coeff%beta(gl%eos_coeff%I_gbs(1),1))
    !    allocate(gl%eos_coeff%gam(gl%eos_coeff%I_gbs(1),1))
    !    allocate(gl%eos_coeff%eps(gl%eos_coeff%I_gbs(1),1))
    !end if
    !
    !! set the parameters
    !gl%eos_coeff%ni(:,1)   = n
    !gl%eos_coeff%ti(:,1)   = t
    !gl%eos_coeff%di(:,1)   = d
    !gl%eos_coeff%p_i(:,1)  = p
    !gl%eos_coeff%gama(:,1) = 1d0
    !gl%eos_coeff%pli(:,1)  = 2d0
    !gl%eos_coeff%tli(:,1)  = 2d0
    !gl%eos_coeff%eta(:,1)  = gauss(:,1)
    !gl%eos_coeff%beta(:,1) = gauss(:,2)
    !gl%eos_coeff%gam(:,1)  =  gauss(:,3)
    !gl%eos_coeff%eps(:,1)  =  gauss(:,4)
    !
    !gl%pc(1) = p_calc(gl,gl%tred(1),gl%rhored(1),1)
    !!calculation of properties
    !! datas(:;1) -> p
    !! datas(:;2) -> t
    !! datas(:,3) -> d
    !! datas(:;1) -> prop
    !
    !do i=1,sizes(4)
    !    val = 0d0
    !    Temp = 0d0
    !    press = 0d0
    !    rhovap_est = 0d0
    !    rholiq_est = 0d0
    !    if(datas_id(i) .eq. 111) then ! DL calculation
    !        iflash = 1
    !        iter = 0
    !        call Flash_Pure_PhaseBoundary(gl ,press, datas(i,2), rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst)
    !        val_calc(i) = rholiq_est/1d3
    !        !if(isnan(val_calc(i))) val_calc(i) = 0d0
    !        trendreals = datas(i,4)*(datas(i,3) - val_calc(i))**2 + trendreals
    !    elseif(datas_id(i) .eq. 112) then ! DV calculation
    !        iflash = 1
    !        iter = 0
    !        call Flash_Pure_PhaseBoundary (gl,press, datas(i,2), rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst)
    !        val_calc(i) = rhovap_est/1d3
    !        !if(isnan(val_calc(i))) val_calc(i) = 0d0
    !        trendreals = datas(i,4)*(datas(i,3) - val_calc(i))**2 + trendreals
    !    elseif(datas_id(i)  .eq. 113) then !PVT calculation
    !        IPHASE = 0
    !        val = rhomix_calc(gl,datas(i,2),datas(i,1), rhovap_est, IPHASE, nrsubst)
    !        val_calc(i) = val/1d3
    !        !if(isnan(val_calc(i))) val_calc(i) = 0d0
    !        trendreals = datas(i,4)*(datas(i,3) - val_calc(i))**2 + trendreals
    !    elseif(datas_id(i) .eq. 114) then ! PV calculation
    !        iflash = 1
    !        iter = 0
    !        call Flash_Pure_PhaseBoundary (gl, press, datas(i,2), rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst)
    !        val_calc(i) = press
    !        !if(isnan(val_calc(i))) val_calc(i) = 0d0
    !        trendreals = datas(i,4)*(datas(i,1) - val_calc(i))**2 + trendreals
    !    end if
    !    if(isnan(trendreals)) then
    !        write(*,*)  trendreals
    !    end if
    !end do
    !call DESTROY_FLUID(handle)
    !end function
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !subroutine trend(x,d_i,p_i,nterms,add_vals,temps,values,sizes,mode,handle)
    !use module_all_types
    !use flash_pure_module
    !implicit none
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "trend" :: trend
    !!input / output variabels
    !integer, dimension(4):: sizes
    !double precision, dimension(sizes(1),2) :: values
    !double precision, dimension(sizes(1)) :: x
    !double precision, dimension(sizes(2)) :: temps
    !integer, dimension(sizes(3)) :: nterms
    !double precision, dimension(sizes(4)):: add_vals
    !integer, dimension(sum(nterms)) :: d_i,p_i
    !integer:: mode ! 1 ssq calc, 2 sat state calc
    !!end input / output variables
    !integer:: nterm,i,errval, iter, nrsubst,fit_mode
    !type(c_ptr):: handle
    !type(type_gl):: gl
    !
    !!
    !!if(.not.c_associated(handle)) then
    !!    call  CREATE_FLUID(handle)
    !!    CALL C_F_POINTER(handle,gl)
    !!    allocate(gl%eos_coeff)
    !!    !  call init_var(gl)
    !!else
    !!    CALL C_F_POINTER(handle,gl)
    !!end if
    !!
    !!call set_parameter_from_extern(gl,x,d_i,p_i,add_vals,nterms,nrsubst,fit_mode)
    !!
    !!gl%pc(1) = p_calc(gl,gl%tred(1),gl%rhored(1),1)
    !!
    !!
    !!!calc all desired densities
    !!do i = 1,size(temps)
    !!    iter= 0
    !!    errval = 0
    !!    call Flash_Pure_PhaseBoundary (gl,0d0, temps(i), values(i,1), values(i,2), 1, errval, iter, 1)
    !!end do
    !
    !end subroutine
    !-------------------------------------------------------------------------------

    !subroutine TREND_FOR_FITTER(eqn,datap,fit_proc,controller,prop_result,cst,eqn_id)
    !
    !!!DEC$ ATTRIBUTES DLLEXPORT :: TREND_FOR_FITTER
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "TREND_FOR_FITTER" :: TREND_FOR_FITTER
    !use crossover
    !
    !implicit none
    !
    !!variable declaration start
    !!initalize the general type_gl
    !type(type_gl), pointer :: gl
    !type(fitter_data_main), intent(inout) :: datap
    !type(levenberg), intent(inout) :: fit_proc
    !type(props), intent(inout):: eqn
    !integer, intent(in):: controller
    !double precision, intent(inout):: prop_result
    !type(constraint),optional:: cst
    !integer,optional:: eqn_id
    !!type(equation_props), intent(inout) :: equation
    !
    !integer:: i,j,k,pot_count,data_loop,data_length
    !integer:: loop_a00,loop_a10,loop_a01,loop_a11,loop_a20,loop_DPDD,loop_D2PDD2
    !integer:: eq_id,nterm
    !type(c_ptr):: func_ptr
    !
    !!matrix for calculating numerical derivatives
    !double precision, allocatable, dimension(:,:,:):: num_derivs_mem
    !integer:: num_deriv_loop
    !double precision:: num_delta
    !
    !!variables for constraints
    !integer:: loop_pots
    !integer:: loop_points
    !integer:: loop_constraints
    !integer:: loop_temp
    !double precision:: calc_val,calc_val_prev
    !double precision :: kernel_part
    !
    !!variables for constraint derivatives
    !double precision :: val_prev
    !double precision,parameter,dimension(2) :: num_del = (/ 1D-8, -2D-8 /)
    !double precision , dimension(2):: val_num_deriv
    !integer:: L,m,n
    !
    !
    !!variable declaration end
    !
    !!DEC ATTRIBUTE FOR PREPROCESSOR
    !!Code is compilated
    !
    !
    !
    !
    !
    !if(.not.c_associated(eqn%handle)) then
    !    call  CREATE_FLUID(eqn%handle)
    !    CALL C_F_POINTER(eqn%handle,gl)
    !    allocate(gl%eos_coeff)
    !    !  call init_var(gl)
    !else
    !    CALL C_F_POINTER(eqn%handle,gl)
    !end if
    !
    !!crossover model is fitted
    !call init_crs(gl,.false.)
    !gl%crs%param = eqn%crs_param
    !
    !
    !!end crossover
    !
    !
    !gl%EQ_type =  1 !has to be changed later!!!
    !gl%mix_Type = 1
    !gl%Factor = 1.D0
    !gl%factorpress = 1.D6
    !gl%factortrans=1.D0
    !gl%factorrbwr=1.D0
    !gl%factor_VS5eta=1.D0
    !gl%factor_VS5slj=1.D0
    !gl%ncomp = 1
    !
    !
    !gl%tred(1) = datap%fluid_param(1,1)
    !gl%rhored(1) = datap%fluid_param(2,1)
    !
    !if(.not. allocated(gl%eos_coeff)) allocate(gl%eos_coeff)
    !gl%eos_coeff%I_pol(1) = eqn%I_P
    !gl%eos_coeff%I_exp(1) = eqn%I_E
    !gl%eos_coeff%I_gbs(1) = eqn%I_G
    !gl%eos_coeff%nreg(1) =  eqn%n_reg
    !!gl%eos_coeff%I_LOG(1) =  eqn%I_L
    !!gl%eos_coeff%I_LOG_D(1) =  eqn%I_L_D
    !nterm  = eqn%n_reg +    eqn%I_G !+ eqn%I_L + eqn%I_L_D
    !
    !if(.not.allocated(gl%eos_coeff%ni)) then
    !    allocate(gl%eos_coeff%ni(nterm,1))
    !    allocate(gl%eos_coeff%ti(nterm,1))
    !    allocate(gl%eos_coeff%di(nterm,1))
    !    allocate(gl%eos_coeff%p_i(nterm,1))
    !    !allocate(gl%eos_coeff%r_i(nterm,1))
    !    !allocate(gl%eos_coeff%q_i(nterm,1))
    !    allocate(gl%eos_coeff%gama(nterm,1))
    !    allocate(gl%eos_coeff%pli(eqn%I_G,1))
    !    allocate(gl%eos_coeff%tli(eqn%I_G,1))
    !    allocate(gl%eos_coeff%eta(eqn%I_G,1))
    !    allocate(gl%eos_coeff%beta(eqn%I_G,1))
    !    allocate(gl%eos_coeff%gam(eqn%I_G,1))
    !    allocate(gl%eos_coeff%eps(eqn%I_G,1))
    !
    !end if
    !
    !do i = 1,gl%eos_coeff%I_pol(1)
    !    gl%eos_coeff%ni(i,1) = eqn%param(i,1)
    !    gl%eos_coeff%ti(i,1) = eqn%param(i,2)
    !    gl%eos_coeff%di(i,1) = eqn%param(i,3)
    !end do
    !
    !gl%eos_coeff%p_i(:,1) =  0
    !gl%eos_coeff%gama(:,1) = 0.d0
    !
    !do i = gl%eos_coeff%I_pol(1)+1,gl%eos_coeff%nreg(1)
    !    gl%eos_coeff%ni(i,1) =   eqn%param(i,1)
    !    gl%eos_coeff%ti(i,1) =   eqn%param(i,2)
    !    gl%eos_coeff%di(i,1) =   eqn%param(i,3)
    !    gl%eos_coeff%p_i(i,1) =  eqn%param(i-gl%eos_coeff%I_pol(1),4)
    !    gl%eos_coeff%gama(i,1) = eqn%param(i-gl%eos_coeff%I_pol(1),5)
    !end do
    !
    !
    !do i = 1,gl%eos_coeff%I_gbs(1)
    !
    !    gl%eos_coeff%ni(i+gl%eos_coeff%nreg(1),1)  =  eqn%param(i+gl%eos_coeff%nreg(1),1)
    !    gl%eos_coeff%ti(i+gl%eos_coeff%nreg(1),1)  =  eqn%param(i+gl%eos_coeff%nreg(1),2)
    !    gl%eos_coeff%di(i+gl%eos_coeff%nreg(1),1)  =  eqn%param(i+gl%eos_coeff%nreg(1),3)
    !    gl%eos_coeff%pli(i,1) = eqn%param(i,6)
    !    gl%eos_coeff%tli(i,1) = eqn%param(i,7)
    !    gl%eos_coeff%eta(i,1) = eqn%param(i,8)
    !    gl%eos_coeff%beta(i,1) =eqn%param(i,9)
    !    gl%eos_coeff%gam(i,1) = eqn%param(i,10)
    !    gl%eos_coeff%eps(i,1) = eqn%param(i,11)
    !
    !end do
    !
    !!gl%eos_coeff%r_i = 0.d0
    !!gl%eos_coeff%q_i = 0.d0
    !
    !
    !!do i=1,gl%eos_coeff%I_LOG(1)
    !!    gl%eos_coeff%ni(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1),1)  =  eqn%param(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1),1)
    !!    gl%eos_coeff%ti(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1),1)  =  eqn%param(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1),2)
    !!    gl%eos_coeff%di(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1),1)  =  eqn%param(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1),3)
    !!    gl%eos_coeff%r_i(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1),1) =  eqn%param(i+gl%eos_coeff%I_exp(1),4)
    !!end do
    !
    !!do i=1,gl%eos_coeff%I_LOG_D(1)
    !!    gl%eos_coeff%ni(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1)+gl%eos_coeff%I_LOG(1),1)  =  eqn%param(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1)+gl%eos_coeff%I_LOG(1),1)
    !!    gl%eos_coeff%ti(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1)+gl%eos_coeff%I_LOG(1),1)  =  eqn%param(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1)+gl%eos_coeff%I_LOG(1),2)
    !!    gl%eos_coeff%di(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1)+gl%eos_coeff%I_LOG(1),1)  =  eqn%param(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1)+gl%eos_coeff%I_LOG(1),3)
    !!    gl%eos_coeff%q_i(i+gl%eos_coeff%nreg(1)+gl%eos_coeff%I_GBS(1)+gl%eos_coeff%I_LOG(1),1) =  eqn%param(i+gl%eos_coeff%I_exp(1)+gl%eos_coeff%I_LOG(1),5)
    !!end do
    !
    !
    !
    !
    !!kernel term variables
    !!! if(eqn%kernel_active) then
    !! !     gl%cro%a20 = eqn%crs_param(1)
    !!      gl%cro%a21 = eqn%crs_param(1)
    !!      gl%cro%gi = eqn%crs_param(2)
    !!      gl%cro%d1 = eqn%crs_param(3)
    !!      gl%cro%v1 = eqn%crs_param(4)
    !!      gl%cro%mo = eqn%crs_param(5)
    !! end if
    !!
    !! !Sum of Square calculation
    !if(controller==1) then
    !    !!loop over potentials
    !    !!ORDER IS IMPORTANT A00,A10,A01,A11,A20 for each potential!
    !
    !
    !    do j=1,datap%n_pot
    !
    !        if(eqn%fit_indicator(eqn%act_fit_comb).eq.2) then !for each datapoint the value of the coefficient functions must be updated
    !            gl%eos_coeff%ni(1:eqn%n_term,1) = eqn%general_coeffs_result(1:eqn%n_term,1,j)
    !        elseif(eqn%fit_indicator(eqn%act_fit_comb).eq.3)then
    !            gl%eos_coeff%ni(1:eqn%n_term,1) = eqn%general_coeffs_result(1:eqn%n_term,1,j)
    !            gl%eos_coeff%ti(1:eqn%n_term,1) = eqn%general_texp_result(1:eqn%n_term,1,j)
    !        end if
    !
    !        gl%tred(1) = datap%fluid_param(1,j)
    !        gl%rhored(1) = datap%fluid_param(2,j)
    !        gl%REQ(1) = datap%r_const(j)
    !        eqn%p_red(j) = P_CALC(gl,gl%tred(1),gl%rhored(1),1)
    !        !fit_proc%r_vec = 0.d0
    !
    !        kernel_part = 0.d0
    !
    !        !-------------- A00 - normal
    !        if(datap%pos_ptypes(1,j) .ne. 0) then
    !            if(datap%pos_ptypes(2,j).ne.0) then
    !                data_length = datap%pos_ptypes(2,j)
    !            else
    !                data_length = datap%n_ptypes(1,j) + datap%pos_ptypes(1,j)
    !            end if
    !            do loop_a00=datap%pos_ptypes(1,j),data_length-1
    !                ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_a00),datap%dd(loop_a00),1,0)
    !
    !                fit_proc%r_vec(loop_a00) = datap%pprop(loop_a00) - A00_CALC(gl,datap%tt(loop_a00),datap%dd(loop_a00),1)       !+ kernel_part
    !            end do
    !        end if
    !
    !        !-------------- A10 - normal
    !        if(datap%pos_ptypes(2,j) .ne. 0) then
    !            if(datap%pos_ptypes(3,j).ne.0) then
    !                data_length = datap%pos_ptypes(3,j)
    !            else
    !                data_length = datap%n_ptypes(2,j) + datap%pos_ptypes(2,j)
    !            end if
    !            do loop_a10=datap%pos_ptypes(2,j),data_length-1
    !                ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_a10),datap%dd(loop_a10),3,0)
    !                fit_proc%r_vec(loop_a10) = datap%pprop(loop_a10) - A10_CALC(gl,datap%tt(loop_a10),datap%dd(loop_a10),1)      !+ kernel_part
    !            end do
    !        end if
    !
    !        !-------------- A01 - normal
    !        if(datap%pos_ptypes(3,j) .ne. 0) then
    !            if(datap%pos_ptypes(4,j).ne.0) then
    !                data_length = datap%pos_ptypes(4,j)
    !            else
    !                data_length = datap%n_ptypes(3,j) + datap%pos_ptypes(3,j)
    !            end if
    !            do loop_a01=datap%pos_ptypes(3,j),data_length-1
    !                !if(eqn%fit_indicator(eqn%act_fit_comb).eq.2) then !for each datapoint the value of the coefficient functions must be updated
    !                !    gl%eos_coeff%ni(1:eqn%n_term,1) = eqn%general_coeffs_result(1:eqn%n_term,1,j)
    !                !end if
    !                ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_a01),datap%dd(loop_a01),2,0)
    !                fit_proc%r_vec(loop_a01) = datap%pprop(loop_a01) - A01_CALC(gl,datap%tt(loop_a01),datap%dd(loop_a01),1)      !+ kernel_part
    !            end do
    !        end if
    !        !-------------- A11 - normal
    !        if(datap%pos_ptypes(4,j) .ne. 0) then
    !            if(datap%pos_ptypes(5,j).ne.0) then
    !                data_length = datap%pos_ptypes(5,j)
    !            else
    !                data_length = datap%n_ptypes(4,j) + datap%pos_ptypes(4,j)
    !            end if
    !            do loop_a11=datap%pos_ptypes(4,j),data_length-1
    !                ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_a11),datap%dd(loop_a11),4,0)
    !                fit_proc%r_vec(loop_a11) = datap%pprop(loop_a11) -A11_CALC(gl,datap%tt(loop_a11),datap%dd(loop_a11),1)      !+ kernel_part
    !            end do
    !        end if
    !        !-------------- A20 - normal
    !        if(datap%pos_ptypes(5,j) .ne. 0) then
    !            if(datap%pos_ptypes(6,j).ne.0) then
    !                data_length = datap%pos_ptypes(6,j)
    !            else
    !                data_length = datap%n_ptypes(5,j) + datap%pos_ptypes(5,j)
    !            end if
    !            do loop_a20=datap%pos_ptypes(5,j),data_length-1
    !                ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_a20),datap%dd(loop_a20),5,0)
    !                fit_proc%r_vec(loop_a20) = datap%pprop(loop_a20) -  A20_CALC(gl,datap%tt(loop_a20),datap%dd(loop_a20),1)      !+ kernel_part
    !            end do
    !        end if
    !
    !
    !
    !        if(fit_proc%fit_crit.eq.1) then ! fitting of critical datapoints
    !
    !            !-------------- A00 - crit
    !            if(datap%pos_ptypes_crit(1,j) .ne. 0) then
    !                if(datap%pos_ptypes_crit(2,j).ne.0) then
    !                    data_length = datap%pos_ptypes_crit(2,j)
    !                else
    !                    data_length = datap%n_ptypes_crit(1,j) + datap%pos_ptypes_crit(1,j)
    !                end if
    !                do loop_a00=datap%pos_ptypes_crit(1,j),data_length-1
    !                    ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_a00),datap%dd(loop_a00),1,0)
    !                    fit_proc%r_vec(loop_a00) = datap%pprop(loop_a00) - A00_CALC(gl,datap%tt(loop_a00),datap%dd(loop_a00),1)      !+ kernel_part
    !                end do
    !            end if
    !
    !            !-------------- A10 - crit
    !            if(datap%pos_ptypes_crit(2,j) .ne. 0) then
    !                if(datap%pos_ptypes_crit(3,j).ne.0) then
    !                    data_length = datap%pos_ptypes_crit(3,j)
    !                else
    !                    data_length = datap%n_ptypes_crit(2,j) + datap%pos_ptypes_crit(2,j)
    !                end if
    !                do loop_a10=datap%pos_ptypes_crit(2,j),data_length-1
    !                    ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_a10),datap%dd(loop_a10),3,0)
    !                    fit_proc%r_vec(loop_a10) = datap%pprop(loop_a10) - A10_CALC(gl,datap%tt(loop_a10),datap%dd(loop_a10),1)      !+ kernel_part
    !                end do
    !            end if
    !
    !            !-------------- A01 - crit
    !            if(datap%pos_ptypes_crit(3,j) .ne. 0) then
    !                if(datap%pos_ptypes_crit(4,j).ne.0) then
    !                    data_length = datap%pos_ptypes_crit(4,j)
    !                else
    !                    data_length = datap%n_ptypes_crit(3,j) + datap%pos_ptypes_crit(3,j)
    !                end if
    !                do loop_a01=datap%pos_ptypes_crit(3,j),data_length-1
    !                    ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_a01),datap%dd(loop_a01),2,0)
    !                    fit_proc%r_vec(loop_a01) = datap%pprop(loop_a01) - A01_CALC(gl,datap%tt(loop_a01),datap%dd(loop_a01),1)      !+ kernel_part
    !                end do
    !            end if
    !            !-------------- A11 - crit
    !            if(datap%pos_ptypes_crit(4,j) .ne. 0) then
    !                if(datap%pos_ptypes_crit(5,j).ne.0) then
    !                    data_length = datap%pos_ptypes_crit(5,j)
    !                else
    !                    data_length = datap%n_ptypes_crit(4,j) + datap%pos_ptypes_crit(4,j)
    !                end if
    !                do loop_a11=datap%pos_ptypes_crit(4,j),data_length-1
    !                    ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_a11),datap%dd(loop_a11),4,0)
    !                    fit_proc%r_vec(loop_a11) = datap%pprop(loop_a11) -A11_CALC(gl,datap%tt(loop_a11),datap%dd(1),1)      !+ kernel_part
    !                end do
    !            end if
    !            !-------------- A20 - crit
    !            if(datap%pos_ptypes_crit(5,j) .ne. 0) then
    !                if(datap%pos_ptypes_crit(6,j).ne.0) then
    !                    data_length = datap%pos_ptypes_crit(6,j)
    !                else
    !                    data_length = datap%n_ptypes_crit(5,j) + datap%pos_ptypes_crit(5,j)
    !                end if
    !                do loop_a20=datap%pos_ptypes_crit(5,j),data_length-1
    !                    ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_a20),datap%dd(loop_a20),5,0)
    !                    fit_proc%r_vec(loop_a20) = datap%pprop(loop_a20) -  A20_CALC(gl,datap%tt(loop_a20),datap%dd(loop_a20),1)      !+ kernel_part
    !                end do
    !            end if
    !
    !
    !            !-------------- DPD - crit
    !            if(datap%pos_ptypes_crit(6,j) .ne. 0) then
    !                if(datap%pos_ptypes_crit(7,j).ne.0) then
    !                    data_length = datap%pos_ptypes_crit(7,j)
    !                else
    !                    data_length = datap%n_ptypes_crit(6,j) + datap%pos_ptypes_crit(6,j)
    !                end if
    !                do loop_DPDD=datap%pos_ptypes_crit(6,j),data_length-1
    !                    ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_DPDD),datap%dd(loop_DPDD),9,0)
    !                    fit_proc%r_vec(loop_DPDD) = datap%pprop(loop_DPDD) -  DPDD_CALC(gl,datap%tt(loop_DPDD),datap%dd(loop_DPDD),1)      !+ kernel_part
    !                end do
    !            end if
    !
    !            !-------------- D2PD2 - crit
    !            if(datap%pos_ptypes_crit(7,j) .ne. 0) then
    !                if(datap%pos_ptypes_crit(8,j).ne.0) then
    !                    data_length = datap%pos_ptypes_crit(8,j)
    !                else
    !                    data_length = datap%n_ptypes_crit(7,j) + datap%pos_ptypes_crit(7,j)
    !                end if
    !                do loop_D2PDD2=datap%pos_ptypes_crit(7,j),data_length-1
    !                    ! if(eqn%kernel_active)  kernel_part = cross_over_model(gl,datap%tt(loop_D2PDD2),datap%dd(loop_D2PDD2),10,0)
    !                    fit_proc%r_vec(loop_D2PDD2) = datap%pprop(loop_D2PDD2) -  D2PDD2_CALC(gl,datap%tt(loop_D2PDD2),datap%dd(loop_D2PDD2),1)      !+ kernel_part
    !                end do
    !            end if
    !
    !
    !        end if
    !
    !
    !
    !
    !    end do
    !
    !    !Check two phase region for datapoints A20 and A11
    !elseif(controller==2) then
    !
    !
    !elseif(controller==3) then !calculate the derivative for all datapoints with respect to the reducing parameters tc and rhoc
    !
    !    allocate(num_derivs_mem(datap%ndata,4,datap%n_pot))
    !
    !    num_delta = 1.d-8
    !    eqn%red_num_derivs = 0.d0
    !
    !    do j=1,datap%n_pot
    !
    !        gl%tred(1) = datap%fluid_param(1,j)
    !        gl%rhored(1) = datap%fluid_param(2,j)
    !        gl%REQ(1) = datap%r_const(j)
    !        eqn%p_red(j) = P_CALC(gl,gl%tred(1),gl%rhored(1),1)
    !
    !        do num_deriv_loop=1,4 !loop over TC+,TC-,RHOC+,RHOC-
    !
    !            !set the correct delta for numerical differentiation
    !            if(num_deriv_loop.eq.1) then
    !                gl%tred(1) = gl%tred(1) + num_delta
    !                eqn%t_red(j) = eqn%t_red(j) + num_delta
    !            elseif(num_deriv_loop.eq.2) then
    !                gl%tred(1) = gl%tred(1)- num_delta
    !                eqn%t_red(j) = eqn%t_red(j) - num_delta
    !            elseif(num_deriv_loop.eq.3) then
    !                gl%rhored(1) = gl%rhored(1) + num_delta
    !                eqn%r_red(j) = eqn%r_red(j) + num_delta
    !            elseif(num_deriv_loop.eq.4) then
    !                gl%rhored(1) = gl%rhored(1) - num_delta
    !                eqn%r_red(j) = eqn%r_red(j) - num_delta
    !            end if
    !
    !            !IF change, change also in general_variations.f90 !!!!!!!!!!!!!!!
    !            eqn%general_coeff_factors(1,:) = eqn%t_red(:)
    !            eqn%general_coeff_factors(2,:) = eqn%t_red(:)**2
    !            eqn%general_texp_factors(1,:)  = eqn%t_red(:)
    !            eqn%general_texp_factors(2,:)  = eqn%t_red(:)**2
    !            eqn%general_texp_factors(3,:)  = eqn%t_red(:)*  eqn%r_red(:)
    !
    !            !after the reducing parameters where changed, also the corrsponding correlation function has to be changed
    !            do k=1,eqn%n_term
    !                !coeff update derivatives
    !                eqn%general_coeffs_result(k,1,j) = eqn%general_coeffs(k,1) &
    !                    &+eqn%general_coeffs(k,2)*eqn%general_coeff_factors(1,j)&
    !                    &+eqn%general_coeffs(k,3)*eqn%general_coeff_factors(2,j)
    !
    !                !texp update for numerical derivatives
    !                eqn%general_texp_result(k,1,j) = eqn%general_texp(k,1) &
    !                    &+eqn%general_texp(k,2)*eqn%general_texp_factors(1,j) &
    !                    &+eqn%general_texp(k,3)*eqn%general_texp_factors(2,j) &
    !                    &+eqn%general_texp(k,4)*eqn%general_texp_factors(3,j)
    !            end do
    !
    !            !update trend intern variables
    !            gl%eos_coeff%ni(1:eqn%n_term,1) = eqn%general_coeffs_result(1:eqn%n_term,1,j)
    !            gl%eos_coeff%ti(1:eqn%n_term,1) = eqn%general_texp_result(1:eqn%n_term,1,j)
    !
    !            !-------------- A00 - normal
    !            if(datap%pos_ptypes(1,j) .ne. 0) then
    !                if(datap%pos_ptypes(2,j).ne.0) then
    !                    data_length = datap%pos_ptypes(2,j)
    !                else
    !                    data_length = datap%n_ptypes(1,j) + datap%pos_ptypes(1,j)
    !                end if
    !                do loop_a00=datap%pos_ptypes(1,j),data_length-1
    !                    num_derivs_mem(loop_a00,num_deriv_loop,j) = A00_CALC(gl,datap%tt(loop_a00),datap%dd(loop_a00),1)
    !                end do
    !            end if
    !
    !            !-------------- A10 - normal
    !            if(datap%pos_ptypes(2,j) .ne. 0) then
    !                if(datap%pos_ptypes(3,j).ne.0) then
    !                    data_length = datap%pos_ptypes(3,j)
    !                else
    !                    data_length = datap%n_ptypes(2,j) + datap%pos_ptypes(2,j)
    !                end if
    !                do loop_a10=datap%pos_ptypes(2,j),data_length-1
    !                    num_derivs_mem(loop_a10,num_deriv_loop,j) = A10_CALC(gl,datap%tt(loop_a10),datap%dd(loop_a10),1)
    !                end do
    !            end if
    !
    !            !-------------- A01 - normal
    !            if(datap%pos_ptypes(3,j) .ne. 0) then
    !                if(datap%pos_ptypes(4,j).ne.0) then
    !                    data_length = datap%pos_ptypes(4,j)
    !                else
    !                    data_length = datap%n_ptypes(3,j) + datap%pos_ptypes(3,j)
    !                end if
    !                do loop_a01=datap%pos_ptypes(3,j),data_length-1
    !                    num_derivs_mem(loop_a01,num_deriv_loop,j) = A01_CALC(gl,datap%tt(loop_a01),datap%dd(loop_a01),1)
    !                end do
    !            end if
    !            !-------------- A11 - normal
    !            if(datap%pos_ptypes(4,j) .ne. 0) then
    !                if(datap%pos_ptypes(5,j).ne.0) then
    !                    data_length = datap%pos_ptypes(5,j)
    !                else
    !                    data_length = datap%n_ptypes(4,j) + datap%pos_ptypes(4,j)
    !                end if
    !                do loop_a11=datap%pos_ptypes(4,j),data_length-1
    !                    num_derivs_mem(loop_a11,num_deriv_loop,j) = A11_CALC(gl,datap%tt(loop_a11),datap%dd(loop_a11),1)
    !                end do
    !            end if
    !            !-------------- A20 - normal
    !            if(datap%pos_ptypes(5,j) .ne. 0) then
    !                if(datap%pos_ptypes(6,j).ne.0) then
    !                    data_length = datap%pos_ptypes(6,j)
    !                else
    !                    data_length = datap%n_ptypes(5,j) + datap%pos_ptypes(5,j)
    !                end if
    !                do loop_a20=datap%pos_ptypes(5,j),data_length-1
    !                    num_derivs_mem(loop_a20,num_deriv_loop,j) = A20_CALC(gl,datap%tt(loop_a20),datap%dd(loop_a20),1)
    !                end do
    !            end if
    !
    !
    !
    !            if(fit_proc%fit_crit.eq.1) then ! fitting of critical datapoints
    !
    !                !-------------- A00 - crit
    !                if(datap%pos_ptypes_crit(1,j) .ne. 0) then
    !                    if(datap%pos_ptypes_crit(2,j).ne.0) then
    !                        data_length = datap%pos_ptypes_crit(2,j)
    !                    else
    !                        data_length = datap%n_ptypes_crit(1,j) + datap%pos_ptypes_crit(1,j)
    !                    end if
    !                    do loop_a00=datap%pos_ptypes_crit(1,j),data_length-1
    !                        num_derivs_mem(loop_a00,num_deriv_loop,j) = A00_CALC(gl,datap%tt(loop_a00),datap%dd(loop_a00),1)
    !                    end do
    !                end if
    !
    !                !-------------- A10 - crit
    !                if(datap%pos_ptypes_crit(2,j) .ne. 0) then
    !                    if(datap%pos_ptypes_crit(3,j).ne.0) then
    !                        data_length = datap%pos_ptypes_crit(3,j)
    !                    else
    !                        data_length = datap%n_ptypes_crit(2,j) + datap%pos_ptypes_crit(2,j)
    !                    end if
    !                    do loop_a10=datap%pos_ptypes_crit(2,j),data_length-1
    !                        num_derivs_mem(loop_a10,num_deriv_loop,j) = A10_CALC(gl,datap%tt(loop_a10),datap%dd(loop_a10),1)
    !                    end do
    !                end if
    !
    !                !-------------- A01 - crit
    !                if(datap%pos_ptypes_crit(3,j) .ne. 0) then
    !                    if(datap%pos_ptypes_crit(4,j).ne.0) then
    !                        data_length = datap%pos_ptypes_crit(4,j)
    !                    else
    !                        data_length = datap%n_ptypes_crit(3,j) + datap%pos_ptypes_crit(3,j)
    !                    end if
    !                    do loop_a01=datap%pos_ptypes_crit(3,j),data_length-1
    !                        num_derivs_mem(loop_a01,num_deriv_loop,j) = A01_CALC(gl,datap%tt(loop_a01),datap%dd(loop_a01),1)
    !                    end do
    !                end if
    !                !-------------- A11 - crit
    !                if(datap%pos_ptypes_crit(4,j) .ne. 0) then
    !                    if(datap%pos_ptypes_crit(5,j).ne.0) then
    !                        data_length = datap%pos_ptypes_crit(5,j)
    !                    else
    !                        data_length = datap%n_ptypes_crit(4,j) + datap%pos_ptypes_crit(4,j)
    !                    end if
    !                    do loop_a11=datap%pos_ptypes_crit(4,j),data_length-1
    !                        num_derivs_mem(loop_a11,num_deriv_loop,j) = A11_CALC(gl,datap%tt(loop_a11),datap%dd(loop_a11),1)
    !                    end do
    !                end if
    !                !-------------- A20 - crit
    !                if(datap%pos_ptypes_crit(5,j) .ne. 0) then
    !                    if(datap%pos_ptypes_crit(6,j).ne.0) then
    !                        data_length = datap%pos_ptypes_crit(6,j)
    !                    else
    !                        data_length = datap%n_ptypes_crit(5,j) + datap%pos_ptypes_crit(5,j)
    !                    end if
    !                    do loop_a20=datap%pos_ptypes_crit(5,j),data_length-1
    !                        num_derivs_mem(loop_a20,num_deriv_loop,j) = A20_CALC(gl,datap%tt(loop_a20),datap%dd(loop_a20),1)
    !                    end do
    !                end if
    !
    !
    !                !-------------- DPD - crit
    !                if(datap%pos_ptypes_crit(6,j) .ne. 0) then
    !                    if(datap%pos_ptypes_crit(7,j).ne.0) then
    !                        data_length = datap%pos_ptypes_crit(7,j)
    !                    else
    !                        data_length = datap%n_ptypes_crit(6,j) + datap%pos_ptypes_crit(6,j)
    !                    end if
    !                    do loop_DPDD=datap%pos_ptypes_crit(6,j),data_length-1
    !                        datap%tt(loop_DPDD) = eqn%t_red(1)
    !                        datap%dd(loop_DPDD) = gl%rhored(1)
    !                        num_derivs_mem(loop_DPDD,num_deriv_loop,j) = DPDD_CALC(gl,datap%tt(loop_DPDD),datap%dd(loop_DPDD),1)
    !                    end do
    !                end if
    !
    !                !-------------- D2PD2 - crit
    !                if(datap%pos_ptypes_crit(7,j) .ne. 0) then
    !                    if(datap%pos_ptypes_crit(8,j).ne.0) then
    !                        data_length = datap%pos_ptypes_crit(8,j)
    !                    else
    !                        data_length = datap%n_ptypes_crit(7,j) + datap%pos_ptypes_crit(7,j)
    !                    end if
    !                    do loop_D2PDD2=datap%pos_ptypes_crit(7,j),data_length-1
    !                        datap%tt(loop_D2PDD2) = eqn%t_red(1)
    !                        datap%dd(loop_D2PDD2) = gl%rhored(1)
    !                        num_derivs_mem(loop_D2PDD2,num_deriv_loop,j) = D2PDD2_CALC(gl,datap%tt(loop_D2PDD2),datap%dd(loop_D2PDD2),1)
    !                    end do
    !                end if
    !
    !
    !            end if
    !
    !
    !            !set back the correct delta for numerical differentiation
    !            if(num_deriv_loop.eq.1) then
    !                gl%tred(1) = gl%tred(1)-num_delta
    !                eqn%t_red(j) = eqn%t_red(j) -num_delta
    !            elseif(num_deriv_loop.eq.2) then
    !                gl%tred(1) = gl%tred(1)+ num_delta
    !                eqn%t_red(j) = eqn%t_red(j) +num_delta
    !            elseif(num_deriv_loop.eq.3) then
    !                gl%rhored(1) = gl%rhored(1) - num_delta
    !                eqn%r_red(j) = eqn%r_red(j) -num_delta
    !            elseif(num_deriv_loop.eq.4) then
    !                gl%rhored(1) = gl%rhored(1) + num_delta
    !                eqn%r_red(j) = eqn%r_red(j) + num_delta
    !
    !            end if
    !
    !        end do
    !
    !        !calc the numerical derivative with all the sub results
    !        eqn%red_num_derivs(1:datap%ndata,1,j) = (num_derivs_mem(1:datap%ndata,1,j)-num_derivs_mem(1:datap%ndata,2,j))/(2.d0*num_delta)!dervis with respect to tc
    !        eqn%red_num_derivs(1:datap%ndata,2,j) = (num_derivs_mem(1:datap%ndata,3,j)-num_derivs_mem(1:datap%ndata,4,j))/(2.d0*num_delta)!dervis with respect to rhocc
    !
    !
    !    end do
    !
    !    !if(controller==3) then
    !    !IF change, change also in general_variations.f90 !!!!!!!!!!!!!!!
    !    eqn%general_coeff_factors(1,:) = 1.d0/eqn%t_red(:)
    !    eqn%general_coeff_factors(2,:) = 1.d0/eqn%t_red(:)*eqn%r_red(:)
    !    eqn%general_coeff_factors(3,:) = 1.d0/eqn%t_red(:)**2
    !    eqn%general_texp_factors(1,:)  = 1.d0/eqn%t_red(:)
    !    eqn%general_texp_factors(2,:)  = 1.d0/eqn%t_red(:)*eqn%r_red(:)
    !    eqn%general_texp_factors(3,:)  = 1.d0/eqn%t_red(:)**2
    !
    !elseif(controller.eq.5) then !calulate cst constraint values differential evolution
    !    !loop over all constraints
    !    ! all grids
    !
    !    do i=1,cst%nr_cst_mem
    !        do j=1,size(cst%grids(i)%grid_1)
    !            do k=1,size(cst%grids(i)%grid_2)
    !                if(trim(cst%calctype(i)) .eq. 'D2PDD2') then
    !                    cst%mem(eqn_id)%calcs(i)%res_row(j)%res_col(k) = D2PDD2_CALC(gl,cst%grids(i)%grid_1(j),cst%grids(i)%grid_2(k),1)
    !                elseif(trim(cst%calctype(i)) .eq. 'DPD') then
    !                    cst%mem(eqn_id)%calcs(i)%res_row(j)%res_col(k) = DPDD_CALC(gl,cst%grids(i)%grid_1(j),cst%grids(i)%grid_2(k),1)
    !                    !write(*,*) cst%mem(eqn_id)%calcs(i)%res_row(j)%res_col(k)
    !                end if
    !            end do
    !        end do
    !    end do
    !
    !elseif(controller.eq.6) then !calulate cst constraint values jacobian matrix (derivatives of constraint functionwith respect to the parameters)
    !
    !    !loop over all constraints
    !    do i=1,cst%nr_cst_mem
    !        j=1
    !        do m=1,size(cst%grids(i)%grid_1)
    !            do n=1,size(cst%grids(i)%grid_2)
    !                cst%mem(eqn_id)%const_dev(m*n*i) = DPDD_CALC(gl,cst%grids(i)%grid_1(m),cst%grids(i)%grid_2(n),1)
    !                !coefficients
    !                do k=1,nterm
    !                    val_prev = gl%eos_coeff%ni(k,1)!save the previous value
    !                    do L=1,2
    !                        gl%eos_coeff%ni(k,1) = gl%eos_coeff%ni(k,1)  + num_del(L)
    !                        val_num_deriv(L) = DPDD_CALC(gl,cst%grids(i)%grid_1(m),cst%grids(i)%grid_2(n),1)
    !                    end do
    !                    cst%mem(eqn_id)%jac_const(m*n*i,k) = (val_num_deriv(1)-val_num_deriv(2))/(2d0*num_del(1))
    !                    gl%eos_coeff%ni(k,1) = val_prev
    !                end do
    !
    !                !texp
    !                do k=1,nterm
    !                    val_prev = gl%eos_coeff%ti(k,1)!save the previous value
    !                    do L=1,2
    !                        gl%eos_coeff%ni(k,1) = gl%eos_coeff%ni(k,1)  + num_del(L)
    !                        val_num_deriv(L) = DPDD_CALC(gl,cst%grids(i)%grid_1(m),cst%grids(i)%grid_2(n),1)
    !                    end do
    !                    cst%mem(eqn_id)%jac_const(m*n*i,nterm+k) = (val_num_deriv(1)-val_num_deriv(2))/(2d0*num_del(1))
    !                    gl%eos_coeff%ti(k,1) = val_prev
    !                end do
    !
    !                !eta
    !                do k=1,eqn%I_G
    !                    val_prev = gl%eos_coeff%eta(k,1)!save the previous value
    !                    do L=1,2
    !                        gl%eos_coeff%eta(k,1) = gl%eos_coeff%eta(k,1)  + num_del(L)
    !                        val_num_deriv(L) = DPDD_CALC(gl,cst%grids(i)%grid_1(m),cst%grids(i)%grid_2(n),1)
    !                    end do
    !                    cst%mem(eqn_id)%jac_const(m*n*i,2*nterm+k) = (val_num_deriv(1)-val_num_deriv(2))/(2d0*num_del(1))
    !                    gl%eos_coeff%eta(k,1)  =val_prev
    !                end do
    !
    !                !beta
    !                do k=1,eqn%I_G
    !                    val_prev = gl%eos_coeff%beta(k,1)!save the previous value
    !                    do L=1,2
    !                        gl%eos_coeff%beta(k,1) = gl%eos_coeff%beta(k,1)  + num_del(L)
    !                        val_num_deriv(L) = DPDD_CALC(gl,cst%grids(i)%grid_1(m),cst%grids(i)%grid_2(n),1)
    !                    end do
    !                    cst%mem(eqn_id)%jac_const(m*n*i,eqn%I_G+2*nterm+k) = (val_num_deriv(1)-val_num_deriv(2))/(2d0*num_del(1))
    !                    gl%eos_coeff%beta(k,1)  =val_prev
    !                end do
    !
    !                !gamma
    !                do k=1,eqn%I_G
    !                    val_prev = gl%eos_coeff%gam(k,1)!save the previous value
    !                    do L=1,2
    !                        gl%eos_coeff%gam(k,1) = gl%eos_coeff%gam(k,1)  + num_del(L)
    !                        val_num_deriv(L) = DPDD_CALC(gl,cst%grids(i)%grid_1(m),cst%grids(i)%grid_2(n),1)
    !                    end do
    !                    cst%mem(eqn_id)%jac_const(m*n*i,2*eqn%I_G+2*nterm+k) = (val_num_deriv(1)-val_num_deriv(2))/(2d0*num_del(1))
    !                    gl%eos_coeff%gam(k,1)  =val_prev
    !                end do
    !
    !                !eps
    !                do k=1,eqn%I_G
    !                    val_prev = gl%eos_coeff%eps(k,1)!save the previous value
    !                    do L=1,2
    !                        gl%eos_coeff%eps(k,1) = gl%eos_coeff%eps(k,1)  + num_del(L)
    !                        val_num_deriv(L) = DPDD_CALC(gl,cst%grids(i)%grid_1(m),cst%grids(i)%grid_2(n),1)
    !                    end do
    !                    cst%mem(eqn_id)%jac_const(m*n*i,3*eqn%I_G+2*nterm+k) = (val_num_deriv(1)-val_num_deriv(2))/(2d0*num_del(1))
    !                    gl%eos_coeff%eps(k,1)  =val_prev
    !                end do
    !
    !            end do
    !        end do
    !    end do
    !
    !
    !end if
    !
    !
    !!------------------------------------------------------------------------------------------------------------------------------------------------------
    !!CONSTRAINTS-------------------------------------------------------------------------------------------------------------------------------------------
    !!------------------------------------------------------------------------------------------------------------------------------------------------------
    !
    !!    if(eqn%constr%constr_on) then  !calculate the constraints
    !
    !!reset the value of the penalty function
    !eqn%constr%penalty_func = 0.d0
    !
    !do loop_constraints=1,eqn%constr%nr_constrains
    !
    !    !check which kind of constraint has to be calculated
    !    !if(eqn%constr%ctype(loop_constraints).eq.2) then !DPD constraint (2-phase region / 1 van der waals loop)
    !    do loop_pots=1,datap%n_pot
    !        do loop_temp =1,eqn%constr%nr_temp(loop_constraints)
    !
    !            !check which constraint
    !            if(eqn%constr%ctype(loop_constraints).eq.1) then
    !                eqn%constr%targ(loop_temp,1,loop_constraints,loop_pots) =  DPDD_CALC(gl,eqn%constr%temp(loop_temp,loop_constraints,loop_pots),eqn%constr%dense(1,loop_constraints,loop_pots),1)
    !            elseif(eqn%constr%ctype(loop_constraints).eq.2) then
    !                eqn%constr%targ(loop_temp,1,loop_constraints,loop_pots) =  D4PD4_CALC(gl,eqn%constr%temp(loop_temp,loop_constraints,loop_pots),eqn%constr%dense(1,loop_constraints,loop_pots),1)
    !            end if
    !
    !            eqn%constr%press(loop_temp,1,loop_constraints,loop_pots) =  P_CALC(gl,eqn%constr%temp(loop_temp,loop_constraints,loop_pots),eqn%constr%dense(1,loop_constraints,loop_pots),1)
    !
    !
    !
    !            eqn%constr%penalty_func =    abs(eqn%constr%targ(loop_temp,1,loop_constraints,loop_pots) - eqn%constr%slack_var(1,loop_constraints,loop_pots))*eqn%constr%psi(loop_temp,1,loop_constraints,loop_pots) + &
    !                & (1.d0/(2.d0*eqn%constr%mue(loop_temp,1,loop_constraints,loop_pots)))*(eqn%constr%targ(loop_temp,1,loop_constraints,loop_pots) - eqn%constr%slack_var(1,loop_constraints,loop_pots))**2 &
    !                & + eqn%constr%penalty_func
    !            calc_val_prev = eqn%constr%targ(loop_temp,1,loop_constraints,loop_pots)
    !            do loop_points=2,eqn%constr%nr_pts(loop_constraints)
    !
    !                !check which constraint
    !                if(eqn%constr%ctype(loop_constraints).eq.1) then
    !                    eqn%constr%targ(loop_temp,loop_points,loop_constraints,loop_pots) =  DPDD_CALC(gl,eqn%constr%temp(loop_temp,loop_constraints,loop_pots),eqn%constr%dense(loop_points,loop_constraints,loop_pots),1)
    !                elseif(eqn%constr%ctype(loop_constraints).eq.2) then
    !                    eqn%constr%targ(loop_temp,loop_points,loop_constraints,loop_pots) =  D4PD4_CALC(gl,eqn%constr%temp(loop_temp,loop_constraints,loop_pots),eqn%constr%dense(loop_points,loop_constraints,loop_pots),1)
    !                end if
    !
    !                eqn%constr%press(loop_temp,loop_points,loop_constraints,loop_pots) =  P_CALC(gl,eqn%constr%temp(loop_temp,loop_constraints,loop_pots),eqn%constr%dense(loop_points,loop_constraints,loop_pots),1)
    !
    !
    !                !eqn%constr%press(loop_temp,loop_points,loop_constraints,loop_pots) =  P_CALC(gl,eqn%constr%temp(loop_temp,loop_constraints,loop_pots),eqn%constr%dense(loop_points,loop_constraints,loop_pots),1)
    !                eqn%constr%penalty_func =    abs(eqn%constr%targ(loop_temp,loop_points,loop_constraints,loop_pots) - eqn%constr%slack_var(loop_points,loop_constraints,loop_pots))*eqn%constr%psi(loop_temp,loop_points,loop_constraints,loop_pots) + &
    !                    & (1.d0/(2.d0*eqn%constr%mue(loop_temp,loop_points,loop_constraints,loop_pots)))*(eqn%constr%targ(loop_temp,loop_points,loop_constraints,loop_pots) -  eqn%constr%slack_var(loop_points,loop_constraints,loop_pots))**2 &
    !                    & + eqn%constr%penalty_func
    !                calc_val_prev = eqn%constr%targ(loop_temp,loop_points,loop_constraints,loop_pots)
    !
    !            end do
    !        end do
    !    end do
    !    !end if
    !end do
    !
    !
    !!write(*,*) 'PENALTY VALUE', eqn%constr%penalty_func
    !! eqn%constr%penalty_func = 0.d0
    !!   end if
    !
    !
    !
    !!Code is not compilated - Output for user information
    !!!DEC$ ELSE
    !!!DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
    !!write(*,*) 'This Subroutine is not available in the Version!'
    !!!DEC$ END IF ! WO_WRITE_TO_CONSOLE
    !!pause
    !
    !end subroutine TREND_FOR_FITTER
    !DEC$ END IF



    !subroutine TEST_CHAR_ARRAYS_STDCALL (save_char_array,simple_character, gl_handle)
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "TEST_CHAR_ARRAYS_STDCALL" :: TEST_CHAR_ARRAYS_STDCALL
    !!DEC$ ATTRIBUTES stdcall, reference :: TEST_CHAR_ARRAYS_STDCALL
    !
    !use, intrinsic :: iso_c_binding



    !
    !implicit none
    !type(type_gl), pointer :: gl
    !type(c_ptr), intent(inout) :: gl_handle
    !integer(int_ptr_kind()), intent(inout) :: save_char_array  !Pointer to a SafeArray structure
    !character(30) :: simple_character
    !character(30), dimension(30) :: fluids
    !character(30), dimension(:), allocatable :: fluids_allo
    !integer, dimension(30) :: EOS_indicator
    !EOS_indicator = 0
    !EOS_indicator(1) = 1
    !EOS_indicator(2) = 1
    !allocate(fluids_allo(30))
    !fluids_allo = fluids
    !call safearray_to_chararray(save_char_array,30, FORTRAN_char_vektor=fluids_allo)
    !!call safearray_to_chararray(save_char_array,30, FORTRAN_char_array=fluids_allo)
    !fluids = fluids_allo
    !
    !call control_fluids(gl,fluids,EOS_indicator,gl_handle)
    !gl%components = fluids
    !gl%eq_type = EOS_indicator
    !fluids(1) = 'methane'
    !fluids(2) = 'thane'
    !fluids(3) = 'ethane'
    !fluids_allo = fluids
    !call chararray_to_safearray(save_char_array,30, FORTRAN_char_vektor=fluids_allo)
    !simple_character = trim(simple_character) // "FORTRAN"
    !!open(unit=197,file="D:\test_char_arrys.txt"
    !deallocate(fluids_allo)
    !end subroutine


    end module interface_routines_module

