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

    ! module for file hdrt_properties_modules.f90
    module hdrt_properties_modules_module
    !global use inclusion
    use module_all_types
    use utility_module
    use waterice_module
    use dryice_module
    use hdrt_params_module
    use hdrt_chem_pot_module

    contains




    !  Properties of the hydrate model
    !  MODULES (property definition) + Hydrate structures definition (sI and sII)

    !******************************************************************************
    !module hdrt_property_definition
    !!******************************************************************************
    !!
    !!  Module with property constants
    !! October 2012; modified May-October 2015
    !! November 2016:
    !    ! universal reference state conditions T0 = T0a = 273.15 K, p0 = p0a = 1 Pa
    !    ! new lattice parameter correlation - conference article EFM 2016
    !    ! REOS - Reference EoS (not GERG)
    !    ! BS - Langmuir from - multilayered water shell by Ballard & Sloan
    !    ! M4 - compressibility from Murnaghan EoS
    !    ! CONS - constant shell radii for all hydrate formers
    !    ! constant gw0 and hw0 for both sI and sII hydrates
    !
    !!/////////////////////////////////////////////////////////////////////////
    !implicit none
    !save
    !logical :: numderiv = .false.
    !
    !!>toggle for turning head row in outputdata files on or off
    !logical :: print_head, show_progress = .false.
    !
    !!>Set all inputs at the beginning in the programcode an let program run "unattended"
    !logical :: unattendmode = .false.
    !!logical :: pythoncalling = .false.
    !
    !integer :: CompositionType !1: Overall Comp   2: Composition of the VaporPhase at Equilibrium (water-free base)
    !
    !!>Calculate values for four phase line only for physical reasonable values of phasefraction (0 <= beta_ph <= 1)
    !logical :: fourphase_physical = .false.
    !
    !
    !!>Variables for Quadruple Lines
    !!!
    !!!Q_map:
    !!!in rows start points of quadruple line "types" are saved:
    !!!VLHI, VLLH, etc.
    !!!in columns start points for a three phase line from a quadruple line are saved
    !!!1: VLwHIw: VLH, LHI
    !!!2: VLwLcH: VLL, LwLcH
    !!!3: not defined for mixtures now
    !!!4: not defined for mixtures now
    !integer, dimension(4, 6) :: Q_map
    !!>Type of the quadruple point. Dimensions:
    !!!(1: Eqtype, 2: steps between pure components 1 is pure guest#1, np_4ph is pure guest#2)
    !character(len=6), dimension(4,np_4ph) :: Q_point_mix
    !!>pressure and Temperature of the corresponding quadruple line
    !!!(1: Eqtype, 2: steps between pure components 1 is pure guest#1, np_4ph is pure guest#2)
    !double precision, dimension(6,np_4ph) :: press_Q_mix, Temp_Q_mix
    !!>Phasefractions are stored in this variable. Dimensions are:
    !!!(1: Eqtype, 2: phase#, 3: steps between pure components 1 is pure guest#1, np_4ph is pure guest#2)
    !double precision, dimension(4,4,np_4ph) :: phasefrac_mix
    !!>Positions of emerging phases at roots are saved in this variable
    !!!(1: Eqtype, 2: root#, 3: emerging phase#)
    !integer, dimension(4,6,2) :: phaseemerge_mix
    !!>Compositions of each phase and each Quadruple Point are saved in these variables. (:,:,1) belongs to the Quadruple Point of pure hydrate: water - guest#1.
    !!!(:,:,100) belongs to the Quadruple Point of pure hydrate: water - guest#2.
    !!!Dimensions are: (1: Eqtype, 2: molfractions of components, 3: steps between pure components 1 is pure guest#1, 100 ist pure guest#2)
    !double precision, dimension(4,30,np_4ph) :: x_Q_ph1_mix = 0.d0, x_Q_ph2_mix = 0.d0, x_Q_ph3_mix = 0.d0, x_Q_ph4_mix = 0.d0
    !logical, dimension (4,np_4ph) :: Qexist_mix
    !
    !!Variables for threephase lines
    !integer :: np_3ph, n_3ph_lines, nflashs, loops
    !character(len=25), dimension(:), allocatable :: trline_mix
    !double precision, dimension(:,:), allocatable :: press_tr_mix, temp_tr_mix
    !
    !! Public (global) arguments
    !character(len=255) :: CompsHdrtsAll = ''  !Guests forming hydrates - variable for the naming of files
    !character(len=255) :: fourphaselinestr, threephaselinestr
    !! Equations of state for fluid phases:
    !character(4) :: EoS_type = 'reos'  ! 'reos' = reference EoSs + EoS-CG mix rules or GERG mix rules if EoS-CG unavailable
    !                                              ! 'gerg' = GERG 2004 with GERG orig. reference state
    !                                              ! 'srk' = SRK with quadratic mixing rules for a and linear mixing rules for b used
    !                                              !         If available, binary interaction parameter kij used
    !! Shell radius def.
    !character(4) :: Ra_def = 'cons'    ! 'r(a)' = T,p dependent shell radii Ra = R0*a0_hdrt/a(T,p)
    !                                              ! 'tref' = constant shell radii Ra = R0(T_ref) * a(T_ref)/a_ref
    !                                              ! 'cons' = constant shell radii Ra, Ra = R_S, R_L ... used since 10.2016
    !
    !!character(2) :: hdrt_structure
    !character(2) :: Langmuir = 'bs'    ! 'jh' = John & Holder (1982)    ... three water shells
    !                                              ! 'ks' = Klauda & Sandler (2002) ... T-fit
    !                                              ! 'bs' = Ballard & Sloan (2002)  ... multilayered cages
    !                                              ! 'vj' = Vins & Jaeger :-) ... multilayered cages + two additional layers
    !                                              ! 'tr' = Trout's group from MIT ... "Experimental" Langmuir constants
    !                                              !                                    --> fit two potential parameters m and r0 to potential derived from Langmuir constant for small and large cavity
    !                                              !                                    Parameters available for: Methane, Ethane, Propane, Argon
    !
    !character(2) :: Latt_par = 'm4' ! 'm4' = correlation based on Murnaghan EoS with dB0 = 4 = const.
    !                                           ! 'bs' = Ballard & Sloan = constant compressibility (kapa1)
    !                                           ! 'vj' = Vins & Jaeger = pressure dependent compress. (kapa1,2)
    !
    !! Reference conditions @ T0, p0
    !double precision :: gw_B0                               ! Reference state Gibbs energy of water in empty lattice [J/mol]
    !double precision :: hw_B0                               ! Reference state Enthalpy of empty lattice [J.mol^-1]
    !double precision :: gw_B0_a0, hw_B0_a0
    !double precision, dimension(2) :: kb_gw                 ! slope and offset of the linear correlation for gw_B0
    !double precision, dimension(2) :: kb_hw                 ! slope and offset of the linear correlation for hw_B0
    !
    !double precision, dimension(3) :: alfa ! thermal expansion parameters
    !double precision :: T0_a, p0_a ! T0_a: [K] reference temperature for a [A]
    !                               ! p0_a: [Pa] reference pressure for a [A] ... p0_a = 1 Pa for 'm4' since October 2016
    !double precision, dimension(30) :: a0_hdrt, a_hc, sigma, eps_k, sigmad, epsd_k       !sigmad and epsd_k set "0" for all single occupand hydrate formers! physically correct: Kihara parameters "JJ"
    !double precision, dimension(30) :: r_vdw    !half Moleculecular distance (vdW Radii) for double occupancy, if =0, then not known for double occ., fitting parameter?!
    !double precision :: a_ref                   ! [A] reference lattice parameter @ T_ref
    !double precision :: T_ref                   ! [K] reference temperature for the shell radii
    !double precision :: p_ref                   ! [Pa] reference pressure for the shell radii
    !double precision :: l_core
    !!double precision, dimension(6) :: R_S, R_L, R_M
    !!integer, dimension(6) :: z_S, z_L, z_M
    !double precision, dimension(3,6) :: R_shell
    !integer, dimension(3,6) :: z_shell
    !
    !integer :: N_hdrts                              !The number of how many components in the mixture form hydrates (with water)
    !integer :: N_guests                             !The number of guest components forming hydrates (WITHOUT water)
    !double precision, dimension(30) :: moles_hdrt, Y_yj_tf   !Sorted molfractions equal to moleslist_hydrate from module_solids
    !double precision, dimension(2,30) :: B0_hdrt         !Bulk modulus of pure hydrate formers
    !double precision, dimension(2) :: B0_m = (/10.d9, 4.d0/)          !Bulk modulus of mixed hydrate
    !                                                                  ! [Pa, -] Bulk modulus and its derivative from Murnaghan EoS - mixed hydrate
    !double precision :: a0_m                        !lattice parameter of mixed hydrate at T0_a & p0_a = 0 Pa
    !double precision :: a_T_ref                     !lattice parameter of mixed hydrate at T_ref & p_ref
    !                                                ! [A] lattice parameter of the gas at T_ref and low p of 1 Mpa
    !
    !integer :: occupmax
    !
    !double precision, dimension(3,30) :: KS_s, KS_l         ! Klauda&Sandler T-fit for the Langmuir constant - only sI & sII
    !
    !double precision, dimension(3) :: kapa
    !double precision, dimension(2,30) :: m_const, rs_const
    !character(6) :: Latt_par_mix = 'vegard'                  !mixing rule for lattice parameter of mixed hydrates
    !                                        !Latt_par_mix = 'volume'/   ! 'volume' = ideal mixing rule for volume of unit cell
    !                                        !Latt_par_mix = 'vegard'/   ! 'vegard' = Vegard's mixing rule for lattice parameter
    !character(3) :: mixdecision                     !mixing rule for Langmuir Constants of double occupied cages
    !!///////////////////////////////////////////////////////////////////////////////////////////////////
    !
    !
    !! -------------------------
    !! Structure of the hydrate:
    !! -------------------------
    !! data hdrt_structure/'s1'/  ! 's1' = hydrate structure sI
    !                             ! 's2' = hydrate structure sII
    !                             ! 'sH' = hydrate structure sH (not programmed yet)
    !! value for 'hdrt_structure' needs to be defined outside of the 'hdrt_property_definition' module
    !! ... needed in searching for stable hydrate structure at given conditions
    !
    !! --------------------------------------------------------
    !! Molar heat capacity of Hydrate (initial values - ice Ih)
    !! --------------------------------------------------------
    !! heat capacity of the empty hydrate lattice can be well approximated
    !! by the cubic, respectively hexagonal, ice Ih ... fit from IAPWS ice Ih EoS
    !! for p0 = 2.0 MPa ... T = 140 - 300 K
    !!data cpH_a/0.12814d0/     ! [J.mol^-1] constant in cp(T): cp = cp_a*T + cp_b
    !!data cpH_b/2.74566d0/     ! [J.mol^-1.K^-1] constant in cp(T) expression
    !! for p0 = 0 MPa ... T = 120 - 300 K
    !double precision ::  cpH_a = 0.12739d0     ! [J.mol^-1] constant in cp(T): cp = cp_a*T + cp_b
    !double precision ::  cpH_b = 2.93401d0     ! [J.mol^-1.K^-1] constant in cp(T) expression
    !! linear correlation for cp(ice Ih) stays the SAME for p0 = 1 Pa!!! - checked November 2016
    !
    !! -----------------------------
    !! Reference conditions @ T0, p0
    !! -----------------------------
    !double precision ::  T0 = 273.15d0        ! [K] reference temperature for Gibbs energy
    !    !data p0/2.0d6/           ! [Pa] reference pressure for Gibbs energy
    !    !........................................................................
    !    !THE REFERENCE POINT OF THE EMPTY BETA LATTICE WAS CHANGED TO p0 = 0 PA
    !    !THE REASON FOR THIS IS, THAT AT 0 PRESSURE ALL HYDRATES OF ONE TYPE MAY
    !    !HAVE COMPARABLE REFERENCE CONDITIONS, SINCE THE PRESSURE TERM VANISHES,
    !    !THE THERMAL EXPANSION TERM IS THE SAME FOR ALL HYDRATES OF THE SAME
    !    !STRUCTURE. THUS ONLY a0_hdrt SHALL MAKE A DIFFERENCE IN THE REFERENCE
    !    !PARAMETERS. Andreas October 2013
    !    !data p0/0.D0/           ! [Pa] reference pressure for Gibbs energy
    !    !........................................................................
    !    !The reference pressure both for the lattice parameter a(T0_a, p0_a) and
    !    !the empty beta-lattice was changed to 1 Pa ... October 2016
    !    !Reason: calculation of entropy of an ideal gas at the ref. pressure 1 Pa
    !double precision ::  p0 = 1.d0           ! [Pa] reference pressure for Gibbs energy = 1 Pa
    !
    !
    !!      INITIATION of other variables
    !
    !!       |           |            |
    !!       |           |            |
    !!       |           |            |
    !!       V           V            V
    !
    !! -----------------------------------
    !! # of H2O molecules in one unit cell
    !! -----------------------------------
    !double precision :: Nw = 0.d0               ! [-] number of water molecules in unit cell of given hydrate structure
    !
    !! ------------------------------
    !! # of cavities per H2O molecule
    !! ------------------------------
    !double precision, dimension(3) :: v_cavi (/0.d0, 0.d0, 0.d0/) ! first small then large cavity (third is for medium cavity in sH structure)
    !integer :: N_cavi = 2                                         !Number of cavities in the given hydrate structure sI & sII = 2, sH = 3
    !                                                              ! number of cavities ... 2 for sI and sII (3 for sH)
    !
    !! --------------------
    !! Lattice parameter a:
    !! --------------------
    !    !data a0_hdrt/0.d0/               ! [A] lattice param @ T0, p0
    !    !data kapa/0.d0, 0.d0, 0.d0/ ! isothermal compressibilty - coefs.
    !
    !! ------------------------------------------------------------
    !! Kihara potential parameters & Reference conditions @ T0, p0:
    !! ------------------------------------------------------------
    !    !data a_hc/0.d0/     ! [A] hard core radius (considered as constant)
    !    !data sigma/0.d0/    ! [A] collision diameter
    !    !data eps_k/0.d0/    ! [K] energy parameter - potential depth
    !
    !
    !    !!New definition if the reference conditions
    !    !!The reference condition according to Ballard and Sloan are devided by a0_hdrt, such that
    !    !!they are universal for each hydrate former of one structure
    !    !data gw_B0_a0/0.d0/        ! Gibbs energy of water in empty lattice [J/mol], devided by the lattice parameter [A] @T0 and p0
    !    !data hw_B0_a0/0.d0/        ! Enthalpy of empty lattice [J.mol^-1], devided by the lattice parameter [A] @T0 and p0
    !
    !! -------------------------------------------
    !! Water shell radii and coordination numbers:
    !! -------------------------------------------
    !    !data R_S/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/      ! small cavity
    !    !data z_S/0, 0, 0, 0, 0, 0/
    !    !data R_L/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/ ! large cavity
    !    !data z_L/0, 0, 0, 0, 0, 0/
    !    !data R_M/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/ ! medium cavity - only for sH hydrate
    !    !data z_M/0, 0, 0, 0, 0, 0/
    !
    !! ----------------------------------
    !! Klauda and Sandler T-fit for C_iJ:
    !! ----------------------------------
    !    !data KS_s/0.d0, 0.d0, 0.d0/
    !    !data KS_l/0.d0, 0.d0, 0.d0/
    !
    !! --------------------------------
    !! Trout's group from MIT for C_iJ:
    !! --------------------------------
    !    !data m_const/0.D0, 0.D0/     !Small & large cavity
    !    !data rs_const/0.D0, 0.D0/    !Small & large cavity
    !
    !end module hdrt_property_definition
    !******************************************************************************
    !******************************************************************************


    !******************************************************************************
    subroutine hdrt_structure_definition(gl,hdrt_structure_flag)
    !******************************************************************************
    !  Definition of the hydrate structure
    !  This subroutine shall be called everytime for the new overall composition!!!!
    !  May 2015 Vaclav
    !---
    !  UPDATE:
    !  With the actual version of the hydrate model, this routine does not need to be called
    !  when the overall composition changes, but the structure of the hydrate is an input now!
    !  Andreas Jäger, November 2017





    implicit none

    type(type_gl) :: gl


    character(30), dimension(30):: COMPONENTSX
    ! Local variables
    double precision :: vol_B
    double precision, dimension(6) :: R_S, R_L, R_M
    double precision, dimension(30) :: fug_g_ref
    integer, dimension(6) :: z_S, z_L, z_M
    integer :: J
    integer :: hdrt_structure_flag !1: sI, 2: sII, 3: sH
    call uppertolower_char(gl%EoS_type, len(gl%EoS_type))
    if ((gl%EoS_type /= 'reos').and.(gl%EoS_type /= 'gerg').and.(gl%EoS_type /= 'srk')) then
        gl%EoS_type = 'reos'
    endif
    ! EoS decision step:
    if (gl%EoS_type=='reos') then
        call Set_coef_CO2(gl,1)        ! Dry ice ref. state is connected to reference equations
        call Set_coef_H20(gl,1)        ! water ice ref. state is connected to reference equations
    elseif (gl%EoS_type=='gerg') then
        call Set_coef_CO2(gl,2)        ! Dry ice ref. state is connected to GERG-2004
        call Set_coef_H20(gl,2)        ! Water ice ref. state is connected to GERG-2004
    elseif (gl%EoS_type=='srk') then
        call Set_coef_CO2(gl,3)        ! Dry ice ref. state is connected to SRK cubic EoS
        call Set_coef_H20(gl,3)        ! Water ice ref. state is connected to SRK cubic EoS
    end if
    if (hdrt_structure_flag /= 0) then

        COMPONENTSX = gl%Fluidlist_hydrate


        call uppertolower_char(gl%ra_def, len(gl%ra_def))
        if ((gl%ra_def /= 'r(a)').and.(gl%ra_def /= 'tref').and.(gl%ra_def /= 'cons')) then
            gl%ra_def = 'cons'
        endif
        call uppertolower_char(gl%Langmuir, len(gl%Langmuir))
        if ((gl%Langmuir /= 'jh').and.(gl%Langmuir /= 'ks').and.(gl%Langmuir /= 'bs').and.(gl%Langmuir /= 'vj').and.(gl%Langmuir /= 'tr')) then
            gl%Langmuir = 'bs'
        endif
        call uppertolower_char(gl%Latt_par, len(gl%Latt_par))
        if ((gl%Latt_par /= 'm4').and.(gl%Latt_par /= 'bs').and.(gl%Latt_par /= 'vj')) then
            gl%Latt_par = 'm4'
        endif
        call uppertolower_char(gl%Latt_par_mix, len(gl%Latt_par_mix))
        if ((gl%Latt_par_mix /= 'vegard').and.(gl%Latt_par_mix /= 'volume')) then
            gl%Latt_par_mix = 'vegard'
        endif

        gl%N_hdrts = gl%nrofhydrateformers    ! number of hydrate formers = water + guests
        ! Note: this variable avoids calling of module_solids in most hdrt subroutines
        ! nrofhydrateformers = water(1st) + guests

        gl%N_guests = gl%nrofhydrateformers - 1   !  number of guest components (only gases - no water)

        gl%moles_hdrt =  gl%moleslist_hydrate !Sorted molfractions according to fluidlist_hydrate (taken from module_solids) - water shall be 1st - CHECK!


        ! Shell radii (R_i) and coordination numbers (z_i)
        z_S = 0
        z_L = 0
        z_M = 0
        R_S = 0.d0
        R_L = 0.d0
        R_M = 0.d0

        ! Default Bulk modulus and reference pressure for lattice parameter
        gl%B0_m(1) = 10.d9         ! [Pa, -] Bulk modulus at zero pressure ... for MIXED hydrate
        gl%B0_m(2) = 4.d0          ! [-] Bulk modulus derivative dB0 = 4 = const. ... for MIXED hydrate
        !p0_a = 0.d0            ! [Pa] reference pressure for a [A]
        gl%p0_a = 1.d0            ! [Pa] reference pressure for a [A] - October 2016 = 1.0 Pa

        ! Bulk moduli for pure hydrate formers
        gl%B0_hdrt(1,1:30) = 10.d9      ! [Pa, -] Bulk modulus at zero pressure (10 GPa)
        gl%B0_hdrt(2,1:30) = 4.d0       ! [-] Bulk modulus derivative dB0 = 4 = const.

        ! Langmuir constant from Temp-fit by Klauda & Sandler (2003)
        gl%KS_s = 0.d0
        gl%KS_l = 0.d0




        ! --------------------------------------------------------
        ! Definition of the parameters for particular hydrate type

        do J = 1,gl%N_guests
            if (COMPONENTSX(J+1)=='co2') then
                call hdrt_params_co2(gl,J,hdrt_structure_flag)
            elseif (COMPONENTSX(J+1)=='methane') then
                call hdrt_params_methane(gl,J,hdrt_structure_flag)
            elseif (COMPONENTSX(J+1)=='ethane') then
                call hdrt_params_ethane(gl,J,hdrt_structure_flag)
            elseif (COMPONENTSX(J+1)=='co') then
                call hdrt_params_co(gl,J,hdrt_structure_flag)
            elseif (COMPONENTSX(J+1)=='nitrogen') then
                gl%numderiv = .true.
                call hdrt_params_nitrogen(gl,J,hdrt_structure_flag)
            elseif (COMPONENTSX(J+1)=='argon') then
                gl%numderiv = .true.
                call hdrt_params_argon(gl,J,hdrt_structure_flag)
            elseif (COMPONENTSX(J+1)=='propane') then
                call hdrt_params_propane(gl,J,hdrt_structure_flag)
            elseif (COMPONENTSX(J+1)=='oxygen') then
                gl%numderiv = .true.
                call hdrt_params_oxygen(gl,J,hdrt_structure_flag)
            end if
            !SH: Added for calculation of pure hydrates in mix hydrate program code
            !if (N_guests == 1) call hdrt_params_refstate_gh(J)
        end do

        ! -----------------------------
        if (hdrt_structure_flag==1) then ! sI hydrate

            ! Lattice parameter a:
            ! --------------------
            !T0_a = 275.15d0        ! [K] reference temperature for a [A] ... sI ... 275.4017d0 = old value for CO2
            !if (Latt_par /= 'm4') then
            !    p0_a = 2.91d6          ! [Pa] reference pressure for a [A] ... for M4 p0_a = 0 Pa
            !end if
            !alfa = (/1.161150840803288d-4, 2.569958160924152d-7, 1.500983807039881d-10/) ! thermal expansion parameters ... our for CO2
            !alfa = (/1.161151d-4, 2.569958d-7, 1.500983666666666d-10/) ! thermal expansion parameters ... our for CO2, rounded (last value *3 = 4.502951 10^-10)

            ! New alfa correlated to both a(T) and a(p) experimental data for T0_a = 273.15 K & p0_a = 1 Pa
            gl%T0_a = 273.15d0        ! [K] reference temperature for a [A]
            gl%alfa = (/8.416128d-5, 7.141596d-8, -2.106472d-10/) ! thermal expansion parameters ... October 2016 (Conference EFM 2016)

            gl%Nw = 46.d0           ! [-] number of water molecules in an unit cell
            gl%v_cavi = (/(1.d0/23.d0), (3.d0/23.d0), 0.d0/) ! small, large cavities in sI

            ! First shell-radius R(1) must be the SMALLEST - numerical reasons

            if (gl%Langmuir=='bs'.or.gl%Langmuir=='MS') then
                ! Water shell radii and coordination numbers (for 'bs'):
                ! ------------------------------------------------------
                R_S = (/3.83d0, 3.96d0, 0.d0, 0.d0, 0.d0, 0.d0/)     ! small cavity
                z_S = (/8, 12, 0, 0, 0, 0/)
                R_L = (/4.06d0, 4.47d0, 4.645d0, 4.25d0, 0.d0, 0.d0/) ! large cavity
                z_L = (/8, 8, 4, 4, 0, 0/)

            elseif (gl%Langmuir=='jh') then
                ! Water shell radii and coordination numbers (for 'jh'):
                ! ------------------------------------------------------
                R_S = (/3.875d0, 6.593d0, 8.056d0, 0.d0, 0.d0, 0.d0/)      ! small cavity
                z_S = (/20, 20, 50, 0, 0, 0/)
                R_L = (/4.152d0, 7.078d0, 8.285d0, 0.d0, 0.d0, 0.d0/) ! large cavity
                z_L = (/21, 24, 50, 0, 0, 0/)

            elseif (gl%Langmuir=='vj') then
                ! Water shell radii and coordination numbers (for 'vj'):
                ! ------------------------------------------------------
                R_S = (/3.83d0, 3.96d0, 6.593d0, 8.056d0, 0.d0, 0.d0/)      ! small cavity
                z_S = (/8, 12, 20, 50, 0, 0/)
                R_L = (/4.06d0, 4.47d0, 4.645d0, 4.25d0, 7.078d0, 8.285d0/) ! large cavity
                z_L = (/8, 8, 4, 4, 24, 50/)

            end if

            gl%a_ref = 12.03d0         ! [A] reference lattice parameter @ T_ref for sI (McMullan & Jeffrey - 1965)
            gl%T_ref = 248.15d0        ! [K] reference temperature for the shell radii for sI (McMullan & Jeffrey - 1965)
            ! p_ref shall be set to 0 MPa for mixed hydrates
            gl%p_ref = 0.d0            ! [Pa] reference pressure for the shell radii for sI (McMullan & Jeffrey - 1965)

            ! Reference state parameters for pure hydrate model - 2014
            ! --------------------------
            !kb_gw = (/-0.1624d0, 1272.4d0/)     ! slope and offset of the linear correlation for gw_B0
            !kb_hw = (/5.588d0, -14476.0d0/)     ! slope and offset of the linear correlation for hw_B0

            ! Reference state parameters for mixed hydrates model - 2016
            ! --------------------------
            gl%kb_gw = (/0.d0, 987.030d0/)      ! CONSTANT gw_B0
            gl%kb_hw = (/0.d0, -4790.157d0/)    ! CONSTANT hw_B0


        elseif (hdrt_structure_flag==2) then ! sII hydrate

            ! Lattice parameter a:
            ! --------------------
            !T0_a = 298.15d0        ! [K] reference temperature for a [A] ... sII - from a(T) by Hester et al.(2007)
            !if (Latt_par /= 'm4') then
            !    p0_a = 2.91d6          ! [Pa] reference pressure for a [A] ... for M4 p0_a = 0 Pa
            !end if
            !alfa = (/6.7659200d-5, 6.170560d-8, -6.2648500d-11/) ! thermal expansion parameters ... Hester et al.(2007) for sII

            ! New alfa correlated to both a(T) and a(p) experimental data for T0_a = 273.15 K & p0_a = 1 Pa
            gl%T0_a = 273.15d0        ! [K] reference temperature for a [A]
            gl%alfa = (/7.988807d-5, 1.686500d-7,  1.363783d-10/) ! thermal expansion parameters ... October 2016 (Conference EFM 2016)

            gl%Nw = 136.d0           ! [-] number of water molecules in an unit cell
            gl%v_cavi = (/(2.d0/17.d0), (1.d0/17.d0), 0.d0/) ! small, large cavities in sII

            ! First shell-radius R(1) must be the SMALLEST - numerical reasons

            if (gl%Langmuir=='bs'.or.gl%Langmuir=='MS') then
                ! Water shell radii and coordination numbers (for 'bs'):
                ! ------------------------------------------------------
                R_S = (/3.748d0, 3.845d0, 3.956d0, 0.d0, 0.d0, 0.d0/)     ! small cavity
                z_S = (/2, 6, 12, 0, 0, 0/)
                R_L = (/4.635d0, 4.715d0, 4.729d0, 0.d0, 0.d0, 0.d0/) ! large cavity
                z_L = (/12, 12, 4, 0, 0, 0/)
                !a_ref = 17.10d0                       ! [A] reference hard core radius ... given in Ballard and Sloan model

            elseif (gl%Langmuir=='jh') then
                ! Water shell radii and coordination numbers (for 'jh'):
                ! ------------------------------------------------------
                R_S = (/3.870d0, 6.667d0, 8.079d0, 0.d0, 0.d0, 0.d0/)      ! small cavity
                z_S = (/20, 20, 50, 0, 0, 0/)
                R_L = (/4.703d0, 7.464d0, 8.782d0, 0.d0, 0.d0, 0.d0/) ! large cavity
                z_L = (/28, 28, 50, 0, 0, 0/)


            elseif (gl%Langmuir=='vj') then
                ! Water shell radii and coordination numbers (for 'vj'):
                ! ------------------------------------------------------
                R_S = (/3.748d0, 3.845d0, 3.956d0, 6.667d0, 8.079d0, 0.d0/)     ! small cavity
                z_S = (/2, 6, 12, 20, 50, 0/)
                R_L = (/4.715d0, 4.729d0, 4.635d0, 7.464d0, 8.782d0, 0.d0/) ! large cavity
                z_L = (/12, 4, 12, 28, 50, 0/)

            end if

            gl%a_ref = 17.31d0         ! [A] reference lattice parameter @ T_ref for sII (Mak & McMullan - 1965)
            gl%T_ref = 253.15d0        ! [K] reference temperature for the shell radii for sII (Mak & McMullan - 1965)
            ! p_ref shall be set to 0 MPa for mixed hydrates
            gl%p_ref = 0.d0            ! [Pa] reference pressure for the shell radii for sII (Mak & McMullan - 1965)

            ! Reference state parameters for pure hydrate model - 2014
            ! --------------------------
            !kb_gw = (/0.d0, 917.816d0/)      ! slope and offset of the linear correlation for gw_B0
            !kb_hw = (/0.d0, -5094.306d0/)    ! slope and offset of the linear correlation for hw_B0

            ! Reference state parameters for mixed hydrates model - 2016
            ! --------------------------
            gl%kb_gw = (/0.d0, 904.887d0/)      ! CONSTANT gw_B0
            gl%kb_hw = (/0.d0, -5044.369d0/)    ! CONSTANT hw_B0

        elseif (hdrt_structure_flag==3) then ! sH hydrate

            gl%N_cavi = 3  ! number of cavities
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)

            write(*,*) ' '
            write(*,*) 'Go AWAY! ... Please let me be!'
            write(*,*) ' '
            write(*,*) 'I do not want to calculate sH hydrate, yet :-('
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE

            ! Reference state parameters
            ! --------------------------
            gl%kb_gw = (/0.d0, 0.d0/)     ! slope and offset of the linear correlation for gw_B0
            gl%kb_hw = (/0.d0, 0.d0 /)    ! slope and offset of the linear correlation for hw_B0

        end if

        gl%R_shell(1,1:6) = R_S
        gl%z_shell(1,1:6) = z_S
        gl%R_shell(2,1:6) = R_L
        gl%z_shell(2,1:6) = z_L
        gl%R_shell(3,1:6) = R_M
        gl%z_shell(3,1:6) = z_M

        !mixing rules for mixed lanmuir constant of double occupied cages
        !gl%mixdecision='add'
        !gl%mixdecision='mul'
        !gl%mixdecision='mol'
        gl%mixdecision='non'

        ! CALL THIS SUBROUTINE EVERYTIME COC IS CHANGED!!! ... COC = constant overall composition
        ! i.e. also outside this subroutine
        call hdrt_ref_a_gw0_hw0_COC(gl)
    endif

    end subroutine hdrt_structure_definition
    !******************************************************************************
    !******************************************************************************




    end module hdrt_properties_modules_module
