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

    ! module for file hdrt_main.f90
    module hdrt_main_module
    !global use inclusion
    use module_all_types
    use module_hdrt_phaselines
    !DEC$ IF DEFINED(HDRT)
    use hdrt_phaselines_calc
    use hdrt_data
    !DEC$ END IF
    use setup_module
    use utility_module
    use initialize_module
    contains




    !DEC$ IF DEFINED(HDRT)

    !integer function exe_mode(gl, mode)
    !

    !implicit none
    !type(type_gl) :: gl
    !integer:: mode
    !
    !if (mode == 1) then
    !    unattendmode = .true.
    !    exe_mode = 1
    !elseif (mode == -1) then
    !    unattendmode = .false.
    !    exe_mode = 0
    !elseif (mode == 0) then
    !    if (unattendmode .eqv. .true.) then
    !        exe_mode = 1
    !    elseif (unattendmode .eqv. .false.) then
    !        exe_mode = 0
    !    endif
    !endif
    !
    !end function

    !integer function DISPLAY_MODE(mode)
    !!DEC$ ATTRIBUTES DLLEXPORT :: DISPLAY_MODE
    !
    !implicit none
    !integer:: mode
    !
    !if (mode == 1) then
    !    gl%show_progress = .true.
    !    display_mode = 1
    !elseif (mode == -1) then
    !    gl%show_progress = .false.
    !    display_mode = 0
    !elseif (mode == 0) then
    !    if (gl%show_progress .eqv. .true.) then
    !        display_mode = 1
    !    elseif (gl%show_progress .eqv. .false.) then
    !        display_mode = 0
    !    endif
    !endif
    !
    !end function

    subroutine messages(gl, pel, flag, str255, str_arr, i, n, k)


    implicit none
    type(type_gl) :: gl
    type(type_hdrt_pheq) :: pel
    integer :: flag, i, ii, n, k
    character (255) :: str255
    character(30), dimension(30) :: str_arr
    !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)

    if (flag == 1) then
        write(*,*)'Error when opening settingsfile.'
    elseif (flag == 2) then

        write(*,*) '*******************************************************'
        write(*,*) 'Components: Composition chosen (sorted!):'

        Do ii = 1, n

            write(*,*) str_arr(ii), ":", gl%molfractions(ii)

        End do
        !write(*,*) 'composition (sorted!):' // str255

        write(*,*) '*******************************************************'

    elseif (flag == 3) then

        write(*,*) 'Default hydrate model (REOS, M4, CONS, BS)'

    elseif (flag == 4) then

        write(*,*) 'Error - equilibrium with DRY ICE => CO2 must be the 2nd component!'

    elseif (flag == 5) then

        write(*,*) trim(str255)
        write(*,*)'Invalid equilibrium type'

    elseif (flag == 6)then

        write(*,'(''+'',A)')'No startvalues for ' // trim(str255) // trim(str_arr(1)) // ' Line for chosen guests.'

    elseif (flag == 7) then

        write(*,*) 'Output file was not found'
        write(*,*) trim(str255)

    elseif (flag == 8) then

        write(*,*)'Calculation of ' // trim(str255) // ' ' // trim(str_arr(1)(2:)) // '          '

    elseif (flag == 9) then

        write(*,*)'Large jump in calculated temperature (> 15%). Values will not be saved.'

    elseif (flag == 10) then

        write(*,*)'Error in threephase_equil_iter. Code: ',i

    elseif (flag == 11) then

        write(*,*)'Error in threephase_equil_iter: maximum iterations reached'

    elseif (flag == 99) then

        write(*,*)''


    elseif (flag == 100) then

        write(*,'(''+'',A,A,I5,A4,I5,A)')trim(str255),': calculating point ',i, ' of ',n+1,' ' // trim(str_arr(1)(2:)) // '          '


    elseif (flag == 200) then

        write(*,*)'Input file does not exist or is used by antoher program.'

    elseif (flag == 201) then

        write(*,*)'Output path does not exist.'

    elseif (flag == 202) then

        write(*,'(''+'',A,I5,A4,I5)')'Calculating Deviation of Point ',i,' of ',n

    elseif (flag == 203) then

        write(*,*)'One EOS is used out of range of validity. Errors MAY occur.'
        write(*,*)

    elseif (flag == 204) then

        write(*,*)'End of data file. Import was successful.'

    elseif (flag == 205) then

        write(*,*)'End of data file. Calculation was successful.'


    elseif (flag == 400) then

        write(*,'(A)')' Init 4phaseline calc'

    elseif (flag == 401) then

        write(*,*)'Densities of Vapor and Guest rich phase very close. Calculation ends.'

    elseif (flag == 402) then

        write(*,'(''+'',A,I5)')trim(str255) // ': calculating point ',i

    elseif (flag == 403) then

        write(*,*)'Calculation of ' // trim(str255)

    elseif (flag == 404) then

        write(*,'(''+'',A,A,I5,A,I2)')trim(str255),': calculating point ', i, ' Finding physical solutions',n

    elseif (flag == 405) then

        write(*,'(''+'',A,A,I5,A,I2,A,I3)')trim(str255),': calculating point ', i, ' Finding physical solutions',n,' Iteration',k

    elseif (flag == 406) then

        write(*,*)'More than two roots found, code needs to be updated!'

    elseif (flag == 407) then

        write(*,'(''+'',A,A,I5,A,I2,A)')trim(str255),': calculating point ', i, ' Found',n ,' physical solutions    '

    elseif ((flag == 408).or.(flag == 409)) then

        write(*,*)
        write(*,*)'Convergence Problem in fourphaselines in ptflash_solid_NC_4P.'
        write(*,*)'It is very likely that there exists no Quadruple point for the binary system of'
        write(*,*)'water and one guest. Calculation of threephase lines starting from ' // pel%phl4(i)%Q_point_mix(k) // 'line'
        write(*,*)'maybe is not possible. gl%Molfractions of Guests from last successfully calculated'
        write(*,*)'point on 4-ph-line:'

        if (flag == 408) then

            write(*,*)pel%phl4(i)%x_Q_ph1_mix(2:3,n-1)*(1.d0-gl%molfractions(1))

        elseif (flag == 409) then

            write(*,*)pel%phl4(i)%x_Q_ph1_mix(2:3,n-1)/(1.d0-pel%phl4(i)%x_Q_ph1_mix(2:3,n-1))

        endif

        write(*,*)

    elseif (flag == 410) then

        write(*,*)'Error in fourphaselines: path does not exist. No File will be written for 4-phase line, but calculation for 3-phase line will continue.'

    elseif (flag == 411) then

        write(*,*)'Quadruple point VLwHIw not present for both hydrate former.'

    elseif (flag == 412) then

        write(*,*)'Quadruple point VLwLxH not present for both hydrate former.'

    elseif (flag == 413) then

        write(*,*)'Quadruple point VLxHIc not present for both hydrate former.'

    elseif (flag == 414) then

        write(*,*)'Quadruple point LxLyHIc not present for both hydrate former.'

    endif
    !DEC$ END IF ! WO_WRITE_TO_CONSOLE

    end subroutine messages

    subroutine HDRT_MAIN(unattendmode, COMPONENTSX, typeofcalc, CompositionTypeIn, moles, threephaselines, threephasepoints, additional_args)
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "HDRT_MAIN" :: HDRT_MAIN


    implicit none
    !type(type_gl), allocatable :: gl
    !For hydrate calculations this is a vector of types with following entries
    !in case of pure hydrates
    !only 1 row : pure hydrate
    !in case of mixed hydrates
    !1st row : mixed hydrate as given
    !2nd row : pure hydrate with 1st given guest
    !3rd row : pure hydrate with 2nd given guest
    !nth row : pure hydrate with (n-1)th given guest
    type(type_gl), dimension(:), allocatable :: gl
    !type for (p)hase (e)equilibrium (l)lines
    type(type_hdrt_pheq), allocatable :: pel

    character(256) :: additional_args
    character(255) :: username               !custom user settings are in User.f95 file
    character(255) :: resultpath         !The directory where the results are stored
    character(255) :: threepl_path       !The directory where the results are stored
    character(255) :: threepd_path       !The directory where the results are stored
    character(255) :: twopd_path         !The directory where the results are stored
    character(255) :: HydNrd_path        !The directory were the results for hydration numbers are stored
    character(255) :: Datenfilepath      !The directory where the daten file is stored
    character(255) :: Fitfilepath        !The directory where the fit file is stored
    character(255) :: Qfilepath          !The directory where the quadruple file with work of formation is stored
    character(16) :: dialogIn          !Variable to store input from console
    character(255) :: Describ_path       !Path for Hydrate model description file
    integer, dimension(30) :: EOS_indicator ! e.g., '1;1;1' for binary mixture Helmholtz EoS + Helmh. mixing rules
    integer :: MIX_indicator
    character (255) :: componentsstrx

    !Variables needed for setup
    character (12) :: input        !Dummy variable for initialization
    !               molesfr.list,fluidlist,molesfr.list for vapor phase, temporary molfraction
    character (255) :: str_temp
    character(8), dimension(20) :: threephaselines

    character(255) :: filepath = ''       !Fluidfile
    double precision, dimension(30):: moles   !molefractions for each component
    double precision, dimension(30):: moles_vap   !molefractions for each component
    double precision, dimension(30):: moles_temp   !molefractions for each component
    double precision, dimension(30):: moles_pure   !molefractions for each component

    double precision:: Temp, press, Dens      !self explaining

    !                   components saved here, components for quaduple point of pure hydrate saved here
    character(30), dimension(30)::COMPONENTSX, componentsx_orig, COMPONENTS_pure
    !                                               molefraction vector of the vapor phase (water-free base)
    double precision, dimension(30)::MOLFRACTIONSX!, molfractions_orig, molfractions_Qpoint, molfraction_vap
    integer, parameter :: settingsunit = 13
    !some integer variables and         Definition of 2-phase equilibrium used at Q-point for work_form
    integer:: i, j, k, n, nn, errorflag,typeofcalc, screen_in, hdrt_model_spec, dT_VHIw_dp_VLwH, threephasepoints, CompositionTypeIn
    integer,dimension(30) :: componentsorder
    !  probably unused, if unattendmode_runonce .eqv. .true. variable will be set to true after tranlasting components from user input to shorter names e.g. methane -> ch4 for short output filenames
    logical:: twophase, translate_components_runonce, unattendmode
    !if .true. then no user inputs in console are required
    !logical:: unattendmode_runonce = .false.
    !variables for translating character inputs to lower case
    character (26):: lower,upper
    character (12):: adj_compname

    !Variables for the uncertainty estimates
    double precision, dimension(2):: dT_ab, dp_ab, work_form_Q, dTdp
    double precision :: chempot_Q, lnfi_Q, cp_WaterIce

    !pure Quadruple point calculation(pure hydrate variables, new vars for mixed hydrates moved to modules)
    character(6), dimension(4)  :: Q_point       !Quadruple point
    double precision, dimension(:,:), allocatable :: x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4
    logical, dimension (4):: Qexist
    integer, dimension(5) :: error_Q, Q_calc
    !double precision, dimension(6,np_4ph) :: T_Q_Hyd_NC, p_Q_Hyd_NC

    if (.not.allocated(x_Q_ph1))allocate(x_Q_ph1(4,30))
    if (.not.allocated(x_Q_ph2))allocate(x_Q_ph2(4,30))
    if (.not.allocated(x_Q_ph3))allocate(x_Q_ph3(4,30))
    if (.not.allocated(x_Q_ph4))allocate(x_Q_ph4(4,30))


    !//////////////////////////////////////////////////////////////////
    if (.not. allocated(gl))  allocate(gl(count(COMPONENTSX(:) /= "")))
    if (.not. allocated(pel)) allocate(pel)
    !Init
    username = ''
    resultpath = ''
    threepl_path = ''
    threepd_path = ''
    twopd_path = ''
    HydNrd_path = ''
    Datenfilepath = ''
    Fitfilepath = ''
    Qfilepath = ''
    !//////////////////////////////////////////////////////////////////
    !Program settings
    !hdrt_model_spec
    !0  :   use settingsfile
    !1  :   default model
    !2  :   custom model
    hdrt_model_spec = 0
    !call user(username)
    open (unit = settingsunit, file='./settings.ini', action='read', iostat=errorflag)
    if (errorflag /= 0) then
        call messages(gl(1), pel,1, trim(''), trim(''),0,0,0)
        return
    else
        do
            read(settingsunit,'(A)', iostat = errorflag)str_temp
            call uppertolower_char(str_temp,len(str_temp))
            if (errorflag < 0) exit
            if (str_temp(1:1) /= '!') then
                i = index(str_temp,'=')
                !username = 'maximilian'
                if (trim(str_temp(1:i-1)) == 'username') then
                    username = adjustl(trim(str_temp(i+1:)))
                elseif (trim(str_temp(1:i-1)) == 'print_head') then
                    if (adjustl(trim(str_temp(i+1:))) == '.false.') then
                        pel%print_head = .false.
                    elseif (adjustl(trim(str_temp(i+1:))) == '.true.') then
                        pel%print_head = .true.
                    endif
                elseif (trim(str_temp(1:i-1)) == 'show_progress') then
                    if (adjustl(trim(str_temp(i+1:))) == '.false.') then
                        pel%show_progress = .false.
                        pel%phl4%show_progress = .false.
                    elseif (adjustl(trim(str_temp(i+1:))) == '.true.') then
                        pel%show_progress = .true.
                        pel%phl4%show_progress = .true.
                    endif
                elseif (trim(str_temp(1:i-1)) == 'fourphase_physical') then
                    if (adjustl(trim(str_temp(i+1:))) == '.false.') then
                        pel%phl4%fourphase_physical = .false.
                    elseif (adjustl(trim(str_temp(i+1:))) == '.true.') then
                        pel%phl4%fourphase_physical = .true.
                    endif
                elseif (trim(str_temp(1:i-1)) == 'resultpath') then
                    resultpath = adjustl(trim(str_temp(i+1:)))
                elseif (trim(str_temp(1:i-1)) == 'datenfilepath') then
                    Datenfilepath = adjustl(trim(str_temp(i+1:)))
                elseif (trim(str_temp(1:i-1)) == 'fitfilepath') then
                    Fitfilepath = adjustl(trim(str_temp(i+1:)))
                elseif (trim(str_temp(1:i-1)) == 'filepath') then
                    filepath = adjustl(trim(str_temp(i+1:)))
                elseif (trim(str_temp(1:i-1)) == 'eos_type') then
                    gl(1)%eos_type = adjustl(trim(str_temp(i+1:)))
                elseif (trim(str_temp(1:i-1)) == 'ra_def') then
                    gl(1)%ra_def = adjustl(trim(str_temp(i+1:)))
                elseif (trim(str_temp(1:i-1)) == 'langmuir') then
                    gl(1)%langmuir = adjustl(trim(str_temp(i+1:)))
                elseif (trim(str_temp(1:i-1)) == 'latt_par') then
                    gl(1)%latt_par = adjustl(trim(str_temp(i+1:)))
                elseif (trim(str_temp(1:i-1)) == 'latt_par_mix') then
                    gl(1)%latt_par_mix = adjustl(trim(str_temp(i+1:)))
                endif
            endif

        enddo
        close(settingsunit)
    endif

    !**************************************************************************************************
    !**************************************************************************************************
    !Several paths have to be set here
    !------
    ! RESULTPATH:   In resultpath the folders ...\"fluidname"\Three_Phase_Lines
    !                                       ...\"fluidname"\Three_Phase_Data
    !                                       ...\"fluidname"\Two_Phase_Data
    !               must exist! e.g. if CO2 hydrate was chose, the folders
    !                                       ...\CO2\Three_Phase_Lines
    !                                       ...\CO2\Three_Phase_Data
    !                                       ...\CO2\Two_Phase_Data
    !               must exist
    !------
    !------
    ! DATENFILEPATH: The path where all the daten files lie.
    !                Daten files have the ending .hdt
    !------
    !------
    ! FITFILEPATH: The path where the fit files will be written to
    !------
    !------
    ! FILEPATH: The path where the folders "fluids" and "binary_mix_files" are
    !------
    If (trim(username) == 'vaclav') then
        !**************************************************************************************************
        !Vaclav
        !--------------------------------------------------------------------------------------------------
        resultpath = 'D:\AV_CR\Gas_hydrates\!Exchange_folder\branches\Files_Results\'
        !resultpath = 'D:\AV_CR\Gas_hydrates\!FORTRAN_projects\Results_mix_branch\'
        Datenfilepath = 'D:\AV_CR\Gas_hydrates\!Exchange_folder\branches\Files_Daten\'
        !Fitfilepath = 'D:\AV_CR\Gas_hydrates\!FORTRAN_projects\Hyd_TR_10_uncer\'
        Fitfilepath = 'D:\AV_CR\Gas_hydrates\!Exchange_folder\branches\Files_CPO\'
        if (gl(1)%EoS_type == 'reos') then
            filepath = 'D:\AV_CR\Gas_hydrates\!Exchange_folder\trunk\Fluid_files\REOS\'
        elseif (gl(1)%EoS_type == 'gerg') then
            filepath = 'D:\AV_CR\Gas_hydrates\!Exchange_folder\trunk\Fluid_files\GERG-2008\'
        end if
        !--------------------------------------------------------------------------------------------------
    ElseIf (trim(username) == 'andreas') then
        !--------------------------------------------------------------------------------------------------
        !Andy
        resultpath = 'D:\Feststoffe\Hydrate\ExchangeSVN\branches\Results\'
        Datenfilepath = 'D:\Feststoffe\Hydrate\ExchangeSVN\trunk\Daten_files\'
        Fitfilepath = 'D:\Feststoffe\Hydrate\ExchangeSVN\branches\Results\'
        if (gl(1)%EoS_type == 'reos') then
            !filepath = 'D:\Feststoffe\Hydrate\HYD_TREND_10\Fluidfiles\REOS\'
            filepath = 'D:\Feststoffe\Hydrate\ExchangeSVN\trunk\Fluid_files\REOS\'
        elseif (gl(1)%EoS_type == 'gerg') then
            filepath = 'D:\Feststoffe\Hydrate\ExchangeSVN\trunk\Fluid_files\REOS\GERG-2008\'
        end if
    ElseIf (trim(username) == 'andreas_ws') then
        !--------------------------------------------------------------------------------------------------
        !Andy
        resultpath = 'C:\Daten\Arbeit\Projekte\Eigene_Projekte\Hydrat_MIX\Results\'
        Datenfilepath = 'C:\Daten\Arbeit\Projekte\Eigene_Projekte\Hydrat_MIX\Hydrate_Development\trunk\Daten_files\'
        Fitfilepath = 'C:\Daten\Arbeit\Projekte\Eigene_Projekte\Hydrat_MIX\Hydrate_Development\trunk\z_CPOs_REOS'
        if (gl(1)%EoS_type == 'reos') then
            !filepath = 'D:\Feststoffe\Hydrate\HYD_TREND_10\Fluidfiles\REOS\'
            filepath = 'C:\Daten\Arbeit\Projekte\Eigene_Projekte\Hydrat_MIX\Hydrate_Development\trunk\Fluid_files\REOS\'
        elseif (gl(1)%EoS_type == 'gerg') then
            filepath = 'C:\Daten\Arbeit\Projekte\Eigene_Projekte\Hydrat_MIX\Hydrate_Development\trunk\Fluid_files\GERG-2008\'
        end if
        !!--------------------------------------------------------------------------------------------------
    ElseIf (trim(username) == 'sebastian') then
        !--------------------------------------------------------------------------------------------------
        !Sebastian
        resultpath = 'D:\OneDrive - ruhr-uni-bochum.de\Projekte\Hydrate\Results\'
        Datenfilepath = 'D:\Repositories\Hydrates\trunk\Daten_files\'
        Fitfilepath = 'D:\Repositories\Hydrates\trunk\z_CPOs_REOS\'
        if (gl(1)%EoS_type == 'reos') then
            !filepath = 'D:\Feststoffe\Hydrate\HYD_TREND_10\Fluidfiles\REOS\'
            filepath = 'D:\Repositories\Hydrates\trunk\Fluid_files\REOS\'
        elseif (gl(1)%EoS_type == 'gerg') then
            filepath = 'D:\Repositories\Hydrates\trunk\Fluid_files\GERG-2008\'
        end if
        !--------------------------------------------------------------------------------------------------
        !!--------------------------------------------------------------------------------------------------
    ElseIf (trim(username) == 'sebastian2') then
        !--------------------------------------------------------------------------------------------------
        !Sebastian
        resultpath = 'D:\OneDrive - ruhr-uni-bochum.de\Projekte\Hydrate\Results2\'
        Datenfilepath = 'D:\Repositories\Hydrates\trunk\Daten_files\'
        Fitfilepath = 'D:\Repositories\Hydrates\trunk\z_CPOs_REOS\'
        if (gl(1)%EoS_type == 'reos') then
            !filepath = 'D:\Feststoffe\Hydrate\HYD_TREND_10\Fluidfiles\REOS\'
            filepath = 'D:\Repositories\Hydrates\trunk\Fluid_files\REOS\'
        elseif (gl(1)%EoS_type == 'gerg') then
            filepath = 'D:\Repositories\Hydrates\trunk\Fluid_files\GERG-2008\'
        end if
        !--------------------------------------------------------------------------------------------------
    ElseIf (trim(username) == 'maximilian') then
        !--------------------------------------------------------------------------------------------------
        !Maxmilian
        resultpath = 'Z:\Hiwis\Bialdyga\Hydrate\Results\'
        Datenfilepath = 'Z:\Hiwis\Bialdyga\Hydrate\Daten_files\'
        Fitfilepath = 'Z:\Programme\Hydrate\Tortoise Repository\trunk\z_CPOs_REOS\'
        if (gl(1)%EoS_type == 'reos') then
            !filepath = 'D:\Feststoffe\Hydrate\HYD_TREND_10\Fluidfiles\REOS\'
            filepath = 'Z:\Hiwis\Bialdyga\Hydrate\fluid_files\REOS\'
        elseif (gl(1)%EoS_type == 'gerg') then
            filepath = 'Z:\Hiwis\Bialdyga\Hydrate\fluid_files\GERG-2008\'
        end if
        !--------------------------------------------------------------------------------------------------
    ElseIf (trim(username) == 'benedikt') then
        !--------------------------------------------------------------------------------------------------
        !Maxmilian
        resultpath = 'Z:\Semrau\Results\'
        Datenfilepath = 'Z:\Semrau\trunk\Daten_files\'
        Fitfilepath = 'Z:\Semrau\trunk\z_CPOs_REOS\'
        if (gl(1)%EoS_type == 'reos') then
            !filepath = 'D:\Feststoffe\Hydrate\HYD_TREND_10\Fluidfiles\REOS\'
            filepath = 'Z:\Semrau\trunk\Fluid_files\REOS\'
        elseif (gl(1)%EoS_type == 'gerg') then
            filepath = 'Z:\Semrau\trunk\Fluid_files\GERG-2008\'
        end if
        !--------------------------------------------------------------------------------------------------
        !default for debugging from python:
        !else
        !    resultpath = 'Z:\Projekte\Hydrate\Results\'
        !    Datenfilepath = 'Z:\Programme\Hydrate\Tortoise Repository\trunk\Daten_files\'
        !    Fitfilepath = 'Z:\Programme\Hydrate\Tortoise Repository\trunk\z_CPOs_REOS\'
        !    if (EoS_type == 'reos') then
        !        !filepath = 'D:\Feststoffe\Hydrate\HYD_TREND_10\Fluidfiles\REOS\'
        !        filepath = 'Z:\Programme\Hydrate\Tortoise Repository\trunk\Fluid_files\REOS\'
        !    elseif (EoS_type == 'gerg') then
        !        filepath = 'Z:\Programme\Hydrate\Tortoise Repository\trunk\Fluid_files\GERG-2008\'
        !    end if
        !**************************************************************************************************
    endif
    !Program settings end
    pel%n_3ph_lines = count(threephaselines(:) /= "")
    if (.not. allocated(pel%phl3)) then
        nn = pel%n_3ph_lines
        if (count(componentsx(:) /= "") >= 3) then
            nn = nn * 2
        endif
        if ((any(componentsx == 'ethanol')).or.(any(componentsx == 'methanol'))) then
            nn = pel%n_3ph_lines
            if (count(componentsx(:) /= "")-1 >= 3) then
                nn = nn * 2
            endif
        endif
        allocate(pel%phl3(nn))
        do i = 1,nn
            pel%phl3(i)%loops = threephasepoints
            pel%phl3(i)%trline_mix = threephaselines(i)
        enddo
    endif

    n = count(componentsx(:) /= "")
    do i = 1,n
        call uppertolower_char(componentsx(i), len(componentsx(i)))
    end do
    !sort guests alphabetically
    componentsx_orig = componentsx
    if (componentsx(1) == "seawater") then
        call Sorter(componentsx(3:n),size(componentsx(3:n),1),len(componentsx(3)))
        n = n - 1
    else
        call Sorter(componentsx(2:n),size(componentsx(2:n),1),len(componentsx(2)))
    endif
    do j = 1,n
        do i = 1,n
            if (trim(componentsx_orig(j)) == trim(componentsx(i))) then
                componentsorder(j) = i
            endif
        enddo
    enddo
    moles(1:n) = moles(componentsorder(1:n))

    ! Generate List of guests for the filenames
    if (componentsx(1) == "seawater") then
        pel%CompsHdrtsAll = trim(COMPONENTSX(2)) // '_' // trim(COMPONENTSX(3))
    else
        pel%CompsHdrtsAll = trim(COMPONENTSX(2))
    endif
    if (n > 2) then
        do i = 3,n
            pel%CompsHdrtsAll = trim(pel%CompsHdrtsAll) // '_' // trim(COMPONENTSX(i))
        end do
    end if
    if (hdrt_model_spec /= 0) then
        call hdrt_custom_model(gl(1), hdrt_model_spec)
    endif
    !/////////////////////////////////////////////////////////////////////////////////////////////////////
    ! Initial calculation of QUADRUPLE POINT for PURE HYDRATE (for water + each guest as binary combination)
    !/////////////////////////////////////////////////////////////////////////////////////////////////////
    translate_components_runonce = .false.
    !fluidl = ''
    !moles = 0.D0
    input = 'tp+'    !Dummy
    Temp = 300.D0   !Dummy
    press = 0.1D0   !Dummy
    twophase = .true.
    EOS_indicator = 0
    EOS_indicator(1:2) = 1
    if (componentsx(1) == "seawater") then
        EOS_indicator(3) = 1
    endif
    MIX_indicator = EOS_indicator(1) !1: default, 110: always linear mixing rule, 120: always lorentz-berthelot
    componentsx_orig = componentsx
    ! ---------------------------
    ! ONLY WATER + ONE GUEST !!!!
    ! ---------------------------
    moles_pure = 0.d0
    moles_pure(1:2) = 0.5d0
    if (componentsx(1) == "seawater") then
        moles_pure(2) = moles(2)
        moles_pure(3) = 0.5d0
    endif
    COMPONENTS_pure = ""
    if (componentsx(1) == "seawater") then
        n = n + 1
    endif

    !Initialize the calculations for (mixed) hydrates
    EOS_indicator(1:n) = 1
    MIX_indicator = 121
    call setup(gl(1),input, Temp, press, COMPONENTSX, moles, filepath, EOS_indicator, MIX_indicator, errorflag)

    if ((any(componentsx == 'ethanol')).or.(any(componentsx == 'methanol'))) then
        n = n - 1
        EOS_indicator(1:n) = 1
        MIX_indicator = 121
        !moles_pure = moles
    endif

    do i = 2,n

        input = gl(1)%inptorig
        EOS_indicator = 0
        EOS_indicator(1:n) = 1
        !Transformation of the fluidlist and moles list to a character
        !This step is only necessary to use setup
        COMPONENTS_pure(1) = COMPONENTSX(1)
        COMPONENTS_pure(2) = COMPONENTSX(i)
        if (componentsx(1) == "seawater") then!.or.((any(componentsx(:) == 'ethanol')).or.(any(componentsx(:) == 'methanol')))
            COMPONENTS_pure(2) = COMPONENTSX(2)
            COMPONENTS_pure(3) = COMPONENTSX(i+1)
            EOS_indicator(3) = 1
        endif

        !Initialize the calculations FOR BINARY MIXTURE only!!
        call setup(gl(i), input, Temp, press, COMPONENTS_pure, moles_pure, filepath, EOS_indicator, MIX_indicator, errorflag)

        !Store the quadruple point information for binary mixtures forming hydrate
        Q_calc = 0
        Q_calc(5) = 1   ! all existing quadruple points are calculated
        call quadruplepoints (gl(i), Fitfilepath, Q_calc, Q_point, Qexist, gl(i)%T_Q_Hyd_2C, gl(i)%p_Q_Hyd_2C, &
            & moles_pure, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)

        if (i == 2) then    !pure hydrate or mixed hydrate guest 1
            pel%phl4(:)%x_Q_ph1_mix(1,1) = x_Q_ph1(:,1) !molfrac of water in ph1
            pel%phl4(:)%x_Q_ph2_mix(1,1) = x_Q_ph2(:,1) !molfrac of water in ph2
            pel%phl4(:)%x_Q_ph3_mix(1,1) = x_Q_ph3(:,1) !molfrac of water in ph3
            pel%phl4(:)%x_Q_ph4_mix(1,1) = x_Q_ph4(:,1) !molfrac of water in ph4

            pel%phl4(:)%x_Q_ph1_mix(i,1) = x_Q_ph1(:,2) !molfrac of guest 1 in ph1
            pel%phl4(:)%x_Q_ph2_mix(i,1) = x_Q_ph2(:,2) !molfrac of guest 1 in ph2
            pel%phl4(:)%x_Q_ph3_mix(i,1) = x_Q_ph3(:,2) !molfrac of guest 1 in ph3
            pel%phl4(:)%x_Q_ph4_mix(i,1) = x_Q_ph4(:,2) !molfrac of guest 1 in ph4

            pel%phl4(:)%Q_point_mix(1) = Q_point       !types of q-points for water + guest 1
            pel%phl4(:)%Temp_Q_mix(1) = gl(i)%T_Q_Hyd_2C(1:4)         !temperature of quarduple point for water +  guest 1
            pel%phl4(:)%press_Q_mix(1) = gl(i)%p_Q_Hyd_2C(1:4)        !pressure of quarduple point for water +  guest 1
            pel%phl4(:)%Qexist_mix(1) = Qexist         !there are up to 4 different quarduple points (e.g. VLwHIw, VLcLwH, LcLwHIc, VLcHIc) if a q point exists for current overall comp the relevant flag is set to true
        elseif (i == 3) then !mixed hydrate guest 2
            pel%phl4(:)%x_Q_ph1_mix(1,np_4ph) = x_Q_ph1(:,1) !molfrac of water in ph1
            pel%phl4(:)%x_Q_ph2_mix(1,np_4ph) = x_Q_ph2(:,1) !molfrac of water in ph2
            pel%phl4(:)%x_Q_ph3_mix(1,np_4ph) = x_Q_ph3(:,1) !molfrac of water in ph3
            pel%phl4(:)%x_Q_ph4_mix(1,np_4ph) = x_Q_ph4(:,1) !molfrac of water in ph4

            pel%phl4(:)%x_Q_ph1_mix(i,np_4ph) = x_Q_ph1(:,2) !molfrac of guest 2 in ph1
            pel%phl4(:)%x_Q_ph2_mix(i,np_4ph) = x_Q_ph2(:,2) !molfrac of guest 2 in ph2
            pel%phl4(:)%x_Q_ph3_mix(i,np_4ph) = x_Q_ph3(:,2) !molfrac of guest 2 in ph3
            pel%phl4(:)%x_Q_ph4_mix(i,np_4ph) = x_Q_ph4(:,2) !molfrac of guest 2 in ph4

            pel%phl4(:)%Q_point_mix(np_4ph) = Q_point       !types of q-points for water + guest 2
            pel%phl4(:)%Temp_Q_mix(np_4ph) = gl(i)%T_Q_Hyd_2C(1:4)         !temperature of quarduple point for water +  guest 2
            pel%phl4(:)%press_Q_mix(np_4ph) = gl(i)%p_Q_Hyd_2C(1:4)        !pressure of quarduple point for water +  guest 2
            pel%phl4(:)%Qexist_mix(np_4ph) = Qexist         !there are up to 4 different quarduple points (e.g. VLwHIw, VLcLwH, LcLwHIc, VLcHIc) if a q point exists for current overall comp the relevant flag is set to true
        endif
    enddo

    !/////////////////////////////////////////////////////////////////////////////////////////////////////
    !/////////////////////////////////////////////////////////////////////////////////////////////////////
    ! From here multicomponent hydrates including pure hydrates are considered.
    !---------------------------------------
    ! Definition of the calculation
    if (typeofcalc /= 0) then

        moles_vap = 0.d0
        if (n == 2) then
            If  (typeofcalc /= 3) then
                moles(1) = 0.5d0
                moles(2) = 0.5d0
            endif
        elseif (n > 2) then
            If  ((typeofcalc /= 3).and.(typeofcalc /= 4)) then
                if (n == 3) then
                    moles(1) = 0.5d0
                    moles(2:3) = 0.25d0
                elseif (n == 4) then
                    moles(1) = 0.25d0
                    moles(2:4) = 0.25d0
                elseif (n == 5) then
                    moles(1:5) = 0.2d0
                end if
                dialogIn = '1'! overall composition is set in this case

            endif




        end if



        !open(unit=36, file='D:\cp-T_ice.txt', status='unknown', action='write', iostat=errorflag)
        !write(36,*)'temp,press,cp'
        !press = 1.d-6
        !do i = 1,3500
        !    temp = i/10.d0
        !    dens = cp_WaterIce(gl(1),Temp, press)
        !    write(36,*)temp,press,dens
        !enddo
        !close(36)
        if (errorflag /= 0) return

        call messages(gl(1), pel,2, str_temp, componentsx,0,n,0)
        call messages(gl(1), pel,3, trim(''), trim(''),0,0,0)

        !"Old" File Format for Pure Hydrates
        !Generate path to daten and fit file
        !The daten file type is .hdt
        if (gl(1)%N_guests == 1) then
            Datenfilepath = trim(Datenfilepath) // trim(pel%CompsHdrtsAll) // '.hdt'
        elseif ((gl(1)%N_guests > 1).and.(translate_components_runonce .eqv. .false.)) then
            !New File Format for mixed hydrates - translate long compnames to short (e.g. methane -> ch4)
            call translate_components(Datenfilepath,pel%CompsHdrtsAll)
            translate_components_runonce = .true.
        endif

        ! Generate text file with description of the Hydrate model
        if (typeofcalc /= 0) then   ! quit criterium

            Describ_path = trim(resultpath) // trim(pel%CompsHdrtsAll) // '\z_Model_description.txt'

            open(unit=36, file=Describ_path, status='unknown', action='write', iostat=errorflag)
            if (errorflag /= 0) then
                write(36,*) ' '
                write(36,*) '******************************************'
                write(36,*) 'Components:'
                Do i = 1, n
                    write(36,*) componentsx(i)
                End do
                write(36,*) ' '
                write(36,*) 'Hydrate model settings'
                write(36,*) '-----------------------------------------'
                write(36,*) 'Equation of state:      ', gl(1)%EoS_type
                write(36,*) 'Compressibility:        ', gl(1)%Latt_par
                write(36,*) 'Shell radii definition: ', gl(1)%Ra_def
            endif
        end if

        ! Storing of Q-points file in component specific folder
        open(unit=12, file=trim(resultpath) // trim(pel%CompsHdrtsAll) // '\Q-points.txt', status='unknown', action='write', iostat=errorflag)
        if (errorflag == 0) then
            j = 1
            do i = 2,n
                !Write Quadruple Points of pure hydrate into file
                !j = 1:   first guest
                !j = np_4ph: second guest
                write(12,*)trim(COMPONENTSX(i))
                do k = 1,4
                    if (pel%phl4(k)%Q_point_mix(j) /= '') then
                        write(12,3024) pel%phl4(k)%Q_point_mix(j), pel%phl4(k)%Temp_Q_mix(j), pel%phl4(k)%press_Q_mix(j), pel%phl4(k)%x_Q_ph1_mix(1,j), pel%phl4(k)%x_Q_ph1_mix(i,j), pel%phl4(k)%x_Q_ph2_mix(1,j), pel%phl4(k)%x_Q_ph2_mix(i,j), &
                            & pel%phl4(k)%x_Q_ph3_mix(1,j), pel%phl4(k)%x_Q_ph3_mix(i,j), pel%phl4(k)%x_Q_ph4_mix(1,j), pel%phl4(k)%x_Q_ph4_mix(i,j)
                    end if
                end do
                j = np_4ph
            end do

            close(12)
        end if


    end if



    Select Case (typeofcalc)
    Case(1)
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*)'Maybe in the future'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        !call threephasepoint()
    Case(2)
        call hydratedatadev(gl, pel, Datenfilepath,Resultpath, filepath, errorflag)
    Case(3)
        call phaselines (gl, pel, Resultpath, moles_vap, errorflag)
    Case(4)
        call phaselines_expert(gl, pel, Resultpath, additional_args, errorflag)
    Case(11)
        !The z_ChemPot_out file type is .cpo
        call GenerateFitFile(gl(1), pel, Fitfilepath,Datenfilepath,typeofcalc,dT_ab,dp_ab, &
            &   Q_point, Qexist,gl(1)%T_Q_Hyd_2C, gl(1)%p_Q_Hyd_2C, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
    Case(12)
        ! Calculation at decreased T
        call GenerateFitFile(gl(1), pel, Fitfilepath,Datenfilepath,typeofcalc,dT_ab,dp_ab, &
            &   Q_point, Qexist,gl(1)%T_Q_Hyd_2C, gl(1)%p_Q_Hyd_2C, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
    Case(13)
        ! Calculation at increased p
        call GenerateFitFile(gl(1), pel, Fitfilepath,Datenfilepath,typeofcalc,dT_ab,dp_ab, &
            &   Q_point, Qexist,gl(1)%T_Q_Hyd_2C, gl(1)%p_Q_Hyd_2C, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
    Case(14)
        ! Optimization of dT and dp from the work of formation
        call GenerateFitFile(gl(1), pel, Fitfilepath,Datenfilepath,typeofcalc,dT_ab,dp_ab, &
            &   Q_point, Qexist,gl(1)%T_Q_Hyd_2C, gl(1)%p_Q_Hyd_2C, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
    Case(15)
        ! Evaluation of the work of formation at the Q-point
        Qfilepath = trim(Fitfilepath) // 'work_form_Q.txt'
        !The path, were the calculated values will be written to
        open(unit=57, file=trim(Qfilepath), status='unknown', action='write', iostat=errorflag)
        if (errorflag /= 0) then
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,*) 'An error occured when trying to generate the fit file'
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        end if
        write(57,5002) 'Q-point [T,p] = ',gl(1)%T_Q_Hyd_2C(1),gl(1)%p_Q_Hyd_2C(1)
        write(57,*) 'dT, dp, work_form(dT,dp)'
        dT_ab = (/0.d0, 0.d0/)
        dp_ab = (/0.d0, 0.d0/)
        !call WorkForm_Qpoint(Tp_Q_Hyd_2C(1,1),Tp_Q_Hyd_2C(1,2),dT_ab,dp_ab,work_form_Q)
        !write(57,5001) dT_ab(2), dp_ab(2), work_form_Q(1)
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'Please define max. TEMPERATURE scatter at the Q-point in K. [0 or 0.1]'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) dT_ab(2)
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'Please define max. PRESSURE scatter at the Q-point in MPa. [0.02]'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) dp_ab(2)
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) ' '

        write(*,*) 'Definition of 2-phase equilibrium used at Q-point for work_form:'
        write(*,*) ' 1 - dT and dp steps calculated from stable VIw or VLw,'
        write(*,*) ' 2 - dT-step and dp-step calculated ONLY from VLw,'
        write(*,*) ' 3 - dT-step and dp-step calculated ONLY from VIw.'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) dT_VHIw_dp_VLwH
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) ' '

        write(*,*) ' Calculate work of formation along an interval? [Yes-1]'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) screen_in

        if (screen_in == 1) then
            NN = 30     ! number of T,p points for the interval calculation
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,*) ' '
            write(*,*) ' Calculate work of formation along quadruple T and p ONLY? [Yes-1]'
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
            read(*,*) k
        else
            NN = 1      ! only one data point calculated
            k = 1
            !call WorkForm_Qpoint(dT_VHIw_dp_VLwH,T_Q_Hyd_2C(1),p_Q_Hyd_2C(1),dT_ab,dp_ab,work_form_Q)
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,*) ' '
            write(*,5000) ' Work form = ',work_form_Q(1),dT_ab(2),work_form_Q(2),dp_ab(2)
            write(*,*) ' '
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE

        end if

        if (k == 1) then
            do j = 1,NN
                ! Work of formation at the lower Q-point VLwHIw
                !call WorkForm_Qpoint(dT_VHIw_dp_VLwH,T_Q_Hyd_2C(1),p_Q_Hyd_2C(1),dT_ab*(j)/NN,dp_ab*(j)/NN,work_form_Q)

                write(57,5001) dT_ab(2)*j/NN, 0.d0, work_form_Q(1)
                write(57,5001) 0.d0, dp_ab(2)*j/NN, work_form_Q(2)
            end do

        else

            dTdp = (/0.d0, 0.d0/)
            NN = 30

            do i = 1,NN+1
                do j = 1,NN+1
                    ! Work of formation at the lower Q-point VLwHIw
                    Temp = gl(1)%T_Q_Hyd_2C(1)-dT_ab(2)*(i-1)/NN
                    !call WorkForm_Qpoint(dT_VHIw_dp_VLwH,Temp,p_Q_Hyd_2C(1),dTdp,dp_ab*(j-1)/NN,work_form_Q)

                    write(57,5001) gl(1)%T_Q_Hyd_2C(1)-Temp,dp_ab(2)*(j-1)/NN,work_form_Q(1)

                end do
            end do

        end if ! if-k work form only along quadruple point T and p

        close(57)

    Case(0)
        !exit
        continue

    End select



    !---------------------------------------

    write(36,*) ''
    write(36,5003) 'Type of calculation: ', typeofcalc
    write(36,5004) 'dT_a_b = ', dT_ab
    write(36,5004) 'dp_a_b = ', dp_ab

    write(36,*) ' '
    write(36,*) 'Langmuir constant:      ', gl(1)%Langmuir

    if (gl(1)%Langmuir=='tr') then
        write(36,*) 'Central-well potential params.:'
        write(36,5007) '  m_const   = ', gl(1)%m_const
        write(36,5007) '  rs_const  = ', gl(1)%rs_const
    else
        write(36,*) 'Potential params.: '
        write(36,*) 'guest       : a_hc [A], sigma [A], eps/kB [K] '
        do i = 2,gl(1)%nrofhydrateformers
            write(36,5050) COMPONENTSX(i), gl(1)%a_hc(i-1), gl(1)%sigma(i-1), gl(1)%eps_k(i-1)
        end do


        write(36,*) ' '
    end if

    write(36,*) 'guest       : a0_hdrt [A]'
    do i = 2,gl(1)%nrofhydrateformers
        write(36,5051) COMPONENTSX(i), gl(1)%a0_hdrt(i-1)
    end do
    write(36,*) ' '
    write(36,5051) 'mixed a0_m  ',gl(1)%a0_m
    write(36,*) ' '


    !write(36,5006) '  a0_hdrt = ', gl(1)%a0_hdrt
    if (gl(1)%Latt_par=='m4') then
        write(36,5007) '  B0_m = ', gl(1)%B0_m(1), gl(1)%B0_m(2)
    else
        write(36,5007) '  kapa = ', gl(1)%kapa(1), gl(1)%kapa(2)
    end if
    write(36,5006) ' gw_B0 = ', gl(1)%gw_B0
    write(36,5006) ' hw_B0 = ', gl(1)%hw_B0
    write(36,*) ' '
    write(36,*) ' '
    write(36,*) 'Notation'
    write(36,*) '-----------------------------------------'
    write(36,*) 'Hydrate compressibility:'
    write(36,*) ' vj - pressure dependent compress. - kapa1,2'
    write(36,*) ' bs - constant compress. - kapa1'
    write(36,*) ' m4 - compress. from Murnaghan EoS (with derB0 = 4)'
    write(36,*) 'Shell radii dependence:'
    write(36,*) ' r(a) = R(T,p) - T,p-dependent Ra'
    write(36,*) ' t0p0 = R(T0,p0) - constant Ra dependent on a0_hdrt(T0,p0)'
    write(36,*) ' cons = R_S, R_L - constant Ra independent of guest-type'
    write(36,*) ' '
    write(36,*) ' Langmuir const. - bs, jh, vj, ks, tr - see modules file'
    write(36,*) ' gw_B0 - Gibbs energy of water in empty lattice [J/mol]'
    write(36,*) ' hw_B0 - Enthalpy of empty lattice [J.mol^-1]'
    write(36,*) ' kapa  - isothermal compressibilty - coefs.'
    write(36,*) ' B0 - bulk modulus and its derivative at 0 Pa'

    close(36)


    !output for call from python
    !py_npoints(1) = np_4ph
    !py_npoints(2) = pel%phl3%nflashs*pel%n_3ph_lines
    !py_npoints(3) = pel%phl3%np_3ph
    !py_press_Q_mix = pel%phl4%press_Q_mix
    !py_Temp_Q_mix = pel%phl4%Temp_Q_mix
    !do i = 1,py_npoints(2)
    !    do j = 1,py_npoints(3)
    !        py_press_tr_mix(i,j) = pel%phl3%press_tr_mix(i,j)
    !        py_temp_tr_mix(i,j) = pel%phl3%temp_tr_mix(i,j)
    !    enddo
    !enddo
    !py_trline_mix = ''
    !do j = 1, py_npoints(2)
    !    py_trline_mix(j) = pel%phl3%trline_mix(j)
    !enddo
    !py_errorflag = errorflag
    !end output for call from python

5000 format(a,e15.8,' at dT[K] = ',f9.4,' and ',e15.8,' at dp[MPa] = ',f9.4)
5001 format(f9.5,f9.5,e15.7)
5002 format(a,f9.4,f9.4)
5003 format(a21,i2)
5004 format(a9,f10.4,' ',f10.4)
5005 format(a,e15.8)
5006 format(a,e20.12)
5007 format(a,e15.8,' ',e15.8)

3024 format (a6,' ',f8.4,' ',f10.6,' ',f16.8,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10,' ',f12.10)
3025 format (a6,' ',f8.4,' ',f10.6,' ',f16.8,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10 &
        &,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10,' ',f12.10)
5050 format(' ',a12,': ',f8.6,', ',f8.6,', ', f8.4)
5051 format(' ',a12,': ',f10.7)

    end subroutine HDRT_MAIN
    !******************************************************************************
    !******************************************************************************

    !>The Routine needs Mainfolder of Data Directory and all Hydrate formers as Input.
    !!
    !!
    !!
    !!The names of thehydrate formers are translated to fit to the name of the datafile:
    !! methane   -> CH4
    !! ethane    -> C2H6
    !! propane   -> C3H8
    !! butane    -> C4H12
    !! argon     -> AR
    !! nitrogen  -> N2
    !! oxygen    -> O2
    !!
    !!These Arguements are merged to one string according to:
    !! path + "\csv\" + "H2O-" guest 1 - guest 2 - ... - guest i
    subroutine translate_components(path, components )
    !******************************************************************************
    implicit none
    character*255 :: path, components, componentslwc
    character*255,dimension(30) :: complwcvec
    integer :: i, nrofcomps


    componentslwc = components
    nrofcomps = 1
    do i = 1,len(componentslwc)
        if (componentslwc(i:i) == '_') then
            componentslwc(i:i) = ' '
            nrofcomps = nrofcomps + 1
        endif
    enddo
    do i = 1,30-nrofcomps
        componentslwc = trim(componentslwc) // " empty"
    enddo
    read(componentslwc,*)complwcvec(1),complwcvec(2),complwcvec(3),complwcvec(4),complwcvec(5),complwcvec(6),complwcvec(7),complwcvec(8),complwcvec(9) &
        & ,complwcvec(10),complwcvec(11),complwcvec(12),complwcvec(13),complwcvec(14),complwcvec(15),complwcvec(16),complwcvec(17),complwcvec(18),complwcvec(19) &
        & ,complwcvec(20),complwcvec(21),complwcvec(22),complwcvec(23),complwcvec(24),complwcvec(25),complwcvec(26),complwcvec(27),complwcvec(28),complwcvec(29),complwcvec(30)
    do i = 1,30
        if (trim(complwcvec(i)) == 'methane') then
            complwcvec(i) = 'ch4'
        elseif (trim(complwcvec(i)) == 'ethane') then
            complwcvec(i) = 'c2h6'
        elseif (trim(complwcvec(i)) == 'propane') then
            complwcvec(i) = 'c3h8'
        elseif (trim(complwcvec(i)) == 'butane') then
            complwcvec(i) = 'c4h12'
        elseif (trim(complwcvec(i)) == 'argon') then
            complwcvec(i) = 'ar'
        elseif (trim(complwcvec(i)) == 'nitrogen') then
            complwcvec(i) = 'n2'
        elseif (trim(complwcvec(i)) == 'oxygen') then
            complwcvec(i) = 'o2'
        elseif (trim(complwcvec(i)) == 'empty') then
            complwcvec(i) = ''
        endif
    enddo
    path = trim(path) // 'csv\h2o'
    do i = 1,nrofcomps
        path = trim(path) // '-' // trim(complwcvec(i))
    enddo

    path = trim(path) // '.csv'


    end subroutine translate_components
    !******************************************************************************
    !******************************************************************************
    subroutine hdrt_custom_model(gl,hdrt_model_spec)

    implicit none
    type(type_gl) :: gl
    integer :: hdrt_model_spec
    !####################################################################################################

    ! Hydrate model definition
    ! Default hydrate model
    if (hdrt_model_spec == 1) then
        gl%EoS_type = 'reos'
        gl%Latt_par = 'm4'
        !Ra_def = 'tref'
        gl%Langmuir = 'bs'

        ! User specified hydrate model
    elseif (hdrt_model_spec == 2) then
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)

        write(*,*) ' '
        write(*,*) 'Please choose fluid EoS (1 = REOS, 2 = GERG, 3 = SRK):'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) hdrt_model_spec
        if (hdrt_model_spec == 1) then
            gl%EoS_type = 'reos'
        elseif (hdrt_model_spec == 2) then
            gl%EoS_type = 'gerg'
        elseif (hdrt_model_spec == 3) then
            gl%EoS_type = 'srk'
        end if
        !EoS_type = 'reos'
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) ' '
        write(*,*) 'Please choose compressibility definition:'
        write(*,*) '(1 = M4 - p-dependent compress. from B0 - Murnaghan EoS)'
        write(*,*) '(2 = BS - constant compress. - kapa1)'
        write(*,*) '(3 = VJ - p-dependent compress. - kapa1,2 )'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) hdrt_model_spec
        if (hdrt_model_spec == 1) then
            gl%Latt_par = 'm4'
        elseif (hdrt_model_spec == 2) then
            gl%Latt_par = 'bs'
        elseif (hdrt_model_spec == 3) then
            gl%Latt_par = 'vj'
        end if
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) ' '
        write(*,*) 'Please choose shell radii dependence:'
        write(*,*) '(1 = R(T,p) - T,p-dependent Ra)'
        write(*,*) '(2 = R(Tref) - constant Ra dependent on a(T_ref))'
        write(*,*) '(3 = R_S, R_L - constant Ra independent of gas-type)'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) hdrt_model_spec
        if (hdrt_model_spec == 1) then
            gl%Ra_def = 'r(a)'
        elseif (hdrt_model_spec == 2) then
            gl%Ra_def = 'tref'
        elseif (hdrt_model_spec == 3) then
            gl%Ra_def = 'cons'
        end if
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) ' '
        write(*,*) 'Please choose potential to be used:'
        write(*,*) '(1 = BS - Multilayered cage - Kihara Potential'
        write(*,*) '(2 = JH - 3 Water shells - Kihara Potential'
        write(*,*) '(3 = VJ - Multilayered + additional layers - Kihara Potential'
        write(*,*) '(4 = KS - Temperature dependent Langmuir constant by Klauda and Sandler'
        write(*,*) '(5 = TR - Central well potential as definited by Trouts group'
        write(*,*) '(6 = MS - Kihara Potential for rod like molecules'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        read(*,*) hdrt_model_spec
        if (hdrt_model_spec == 1) then
            gl%Langmuir = 'bs'
        elseif (hdrt_model_spec == 2) then
            gl%Langmuir = 'jh'
        elseif (hdrt_model_spec == 3) then
            gl%Langmuir = 'vj'
        elseif (hdrt_model_spec == 4) then
            gl%Langmuir = 'ks'
        elseif (hdrt_model_spec == 5) then
            gl%Langmuir = 'tr'
        elseif (hdrt_model_spec == 6) then
            gl%Langmuir = 'MS'
        end if

    end if
    ! Hydrate model definition

    !write(*,*) ' '
    !write(*,*) '*******************************************************'
    !####################################################################################################
    end subroutine hdrt_custom_model
    !DEC$ ELSE
    subroutine HDRT_MAIN()

    end subroutine
    !DEC$ END IF


    end module hdrt_main_module
