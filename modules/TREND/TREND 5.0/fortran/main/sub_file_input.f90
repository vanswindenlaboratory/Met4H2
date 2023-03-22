!!    !*************************************************************************************
!!    !				TREND Version 4.0
!!    !		   Thermodynamic Reference & Engineering Data
!!    !
!!    !- software for the calculation of thermodynamic and other properties -
!!    !
!!    !Copyright (C) 2019,  Prof. Dr.-Ing. R.Span
!!    !                     Lehrstuhl fuer Thermodynamik
!!    !                     Ruhr-Universitaet Bochum
!!    !                     Universitaetsstr. 150
!!    !                     D-44892 Bochum
!!    !
!!    !Cite as: Span, R.; Beckmuller, R.; Eckermann, T.; Herrig, S.; Hielscher, S.;
!!    !          Jager, A.; Mickoleit, E.; Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (2019):
!!    !          TREND. Thermodynamic Reference and Engineering Data 4.0.
!!    !          Lehrstuhl fur Thermodynamik, Ruhr-Universitat Bochum.
!!
!!    !
!!    !This program is free software: you can redistribute it and/or modify
!!    !it under the terms of the GNU General Public License as published by
!!    !the Free Software Foundation, either version 3 of the License, or
!!    !(at your option) any later version.
!!    !
!!    !This program is distributed in the hope that it will be useful,
!!    !but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    !MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    !GNU General Public License for more details.
!!    !
!!    !You should have received a copy of the GNU General Public License
!!    !along with this program.  If not, see < http://www.gnu.org/licenses/ >.
!!    !
!!    !*************************************************************************************
!!
!!    ! module for file sub_file_input.f90
!!    module sub_file_input_module
!!    !global use inclusion
!!    use module_all_types
!!    use calc_functions
!!    use rhomix_pt_module
!!    use utility_module
!!    use initialize_module
!!    use crit_tpd_module
!!    use reduced_parameters_calc_module
!!    use cubic_eos_module
!!
!!    contains
!!
!!
!!
!!
!!
!!    !  --------------------------------------------------
!!    !  Sub for reading Fluid-Files
!!    !  J. Moeller, Denmark, 09-2009
!!    !  --------------------------------------------------
!!
!!
!!    !Input:
!!    !  a) vector of components like
!!    !       character(len=12), dimension(2) :: compos        !has to be transfered to this subroutine
!!    !       compos(1) = 'Propane'
!!    !       compos(2) = 'Argon'
!!    !  b) vector of molfractions like
!!    !       double precision, dimension(2) :: molfracs      !has to be transfered to this subroutine
!!    !       molfracs(1) = 0.5d0
!!    !       molfracs(2) = 0.5d0
!!    !  c) complete directory where to find the folder containing the folder for the fluid files and
!!    !     the folder for the binay_mix_files
!!    !             if your are using the REFPROP folder don't forget to put the binay_mix_files in this folder
!!    !       character(len=28) :: path = 'C:\Programme\REFPROP\
!!
!!    !Output: errornumber (no error = 0)
!!
!!    subroutine read_from_fluidfiles (gl,path, errorflag)
!!
!!
!!    !Used modules:
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!
!!    !Declarations: -----------------------------------------------------------------------------------------
!!
!!
!!
!!    implicit none
!!
!!    type(type_gl) :: gl
!!
!!
!!    !COSTALD equation of state
!!    character*33 :: fluidtait_COSTALD
!!    character*15 :: CAS_COSTALD
!!
!!    double precision :: vstar_COSTALD, omega_COSTALD
!!
!!    !Define Input and Output of the Subroutine:
!!    character(len=255) :: path      !has to be transfered to this subroutine
!!    integer :: errorflag
!!
!!    !different paths for singles or binarys   T. Wiens
!!    character(255) :: pathforsingles      !created from path + fluids
!!    character(255) :: pathforbinarys      !created from path + binary_mix_files
!!
!!    !for reading files:
!!    character(len=255) :: filename! = ''                   !later fitted to component file !Sebastian, Monika, Sven !Keine save zuweisungen �ber "="
!!
!!    !Dummy-Variable:
!!    character(len=255) :: dummy !used to jump over the not needed lines of the inputfile
!!    character(len=255) :: dummygen, dummy_GEOS
!!    integer :: i,j,k,ln,m,n        !to count
!!    integer :: nrsubst        !current substance (number)
!!    integer :: error
!!    !!    integer :: ref_error      !for errors occuring while calculating the reference state
!!    integer :: line                    !for ufc/psrk calculation
!!    logical :: new_group            !for ufc/psrk calculation
!!    logical :: new_subgroup        !for ufc/psrk calculation
!!    double precision :: groupnr       !for ufc/psrk calculation
!!    double precision :: subgroupnr    !for ufc/psrk calculation
!!    double precision :: amount        !for ufc/psrk calculation
!!    double precision :: Qkcoeff             !for ufc/psrk calculation
!!    double precision, dimension(100) :: Qkcoeffs !for ufc/psrk calculation
!!    double precision, dimension(30) :: tccur    !current temperature at critical point [K]
!!    double precision, dimension(30) :: pccur    !current pressure at critical point [MPa]
!!    double precision, dimension(30) :: rhoccur  !current density at critical point [mol/m設
!!    double precision, dimension(30) :: ttpcur   !current temperature at triple point [K]
!!    double precision, dimension(30) :: tnbpcur  !current temperature at normal boiling point [K]
!!    double precision, dimension(30) :: wmcur    !current molecular weight [kg/mol]
!!    double precision, dimension(30) :: accencur !current acentric factor [-]
!!    double precision, dimension(30) :: qc_inv, qd_inv   !input parameters viscosity
!!    integer:: ide,ist,loop,ibwr,lbwr
!!    double precision :: dummy1,dummy2,dummy3,dummy4,dummy5             !dummy for reading and sorting regular pol/exp-terms
!!    integer:: jeta, ieta, keta, leta, meta, neta
!!
!!    ! double precision :: P_CALC  !needed to calculate critical pressure from EOS
!!
!!    !Control-Boolean:
!!    Logical :: cpp_read   !cpp has been read
!!    Logical :: phk_read   !phk has been read
!!    Logical :: pho_read   !pho has been read
!!    !Control-Boolean:
!!    Logical :: feq_read   !eq has been read
!!    Logical :: cp0_read   !cp0 has been read
!!    !    Logical :: eta_read   !viscosity model has been read
!!    ! T. Wiens 07.02.2012
!!    logical :: exists
!!
!!    logical :: gl_aga8
!!    !CAS ID as input
!!    integer :: e
!!    character (30), dimension(:) , allocatable :: casnr, fluidname
!!    character (30), dimension(:,:), allocatable :: dummy_arr
!!    character (255) :: CAS_list
!!    logical :: CAS_list_exists
!!
!!    !Variables needed for reading SRK /PR/LKP/ PCSAFT file
!!    character(255):: fluidnames_SRK, fluidnames_PR, fluidnames_LKP, fluidnames_PCSAFT, fluidnames_PCSAFT1, fluidnames_PCSAFT2, fluidnames_GEOS
!!    character(255):: flname_SRK_alter, flname_PR_alter, flname_LKP_alter, flname_PCSAFT_alter, flname_PCSAFT_alter1, flname_PCSAFT_alter2, flname_GEOS_alter
!!    character(255), dimension(:), allocatable :: CASNR_SRK, CASNR_PR, CASNR_LKP, CASNR_PCSAFT, CASNR_GEOS
!!    character(255) :: CASNR_1, CASNR_2  !In TREND the CASNR of the substances is used for the identification of mixing rules
!!    integer, dimension(30) :: Kompnr_SRK !Variables needed for messer kij indentification
!!    integer :: Kompnr1, Kompnr2          !Variables needed for messer kij indentification
!!    double precision :: kij_read, lij_read
!!    double precision :: accen_SRK, accen_PR, accen_LKP, accen_PCSAFT, accen_GEOS
!!    double precision :: Tc_SRK, Tc_PR, Tc_LKP, Tc_PCSAFT, Tc_GEOS
!!    double precision :: pc_SRK, pc_PR, pc_LKP, dc_PCSAFT, pc_GEOS
!!    double precision :: ttp_SRK, ttp_PR, ttp_LKP, ttp_PCSAFT, ttp_GEOS
!!    double precision :: ptp_SRK, ptp_PR, ptp_LKP, ptp_PCSAFT, ptp_GEOS
!!    double precision :: MW_SRK, MW_PR, MW_LKP, MW_PCSAFT, MW_GEOS
!!    double precision :: A_cv0_SRK, A_cv0_PR, A_cv0_LKP, A_cv0_PCSAFT
!!    double precision :: B_cv0_SRK, B_cv0_PR, B_cv0_LKP, B_cv0_PCSAFT
!!    double precision :: C_cv0_SRK, C_cv0_PR, C_cv0_LKP, C_cv0_PCSAFT
!!    double precision :: D_cv0_SRK, D_cv0_PR, D_cv0_LKP, D_cv0_PCSAFT
!!    double precision :: E_cv0_SRK, E_cv0_PR, E_cv0_LKP, E_cv0_PCSAFT
!!    double precision :: F_cv0_SRK, F_cv0_PR, F_cv0_LKP, F_cv0_PCSAFT
!!    double precision :: G_cv0_SRK, G_cv0_PR, G_cv0_LKP, G_cv0_PCSAFT
!!    character(len=512) :: dummy_SRK, dummy_PR, dummy_LKP, dummy_COSTALD, dummy_PCSAFT !used to jump over the not needed lines of the inputfile
!!    double precision :: m_PCSAFT, sig_PCSAFT, epsk_PCSAFT, epsabk_PCSAFT, kap_PCSAFT, Q_PCSAFTQ, D_PCSAFTD
!!    integer:: io_error, lengname, lengname_ext, nD_PCSAFTD, nQ_PCSAFTQ
!!    character(2):: model_PCSAFT
!!    character(25):: form, lengn, lengntr
!!    double precision :: rhoc_est, tc_eos, rhoc_eos, pc_eos
!!    integer :: iter, number, riter
!!    character(len=4) cpdummy
!!
!!    !Variables needed for reading cp0-model of Joback
!!    double precision, dimension (4) :: dA, dB, dC, dD   !dummy variables
!!    integer,dimension(30) :: N_m     !number of molecules of the type
!!
!!    !Variables needed for calculating General equation
!!    integer ::  g_eq, ncomp_orig
!!    double precision:: rhoc_GEOS, polfac_eq_GEOS
!!
!!    double precision, dimension(:,:), allocatable :: a_sun1   !factor for calculation of the generalized EOS of Sun & Ely
!!    double precision, dimension(:,:), allocatable :: a_sun2   !factor for calculation of the generalized EOS of Sun & Ely
!!
!!    double precision :: factor_inv
!!
!!    !New variables for gE-models
!!    integer :: nr_main_group, nr_sub_group
!!    double precision:: vdW_vol_Rk, vdW_ar_Qk
!!    integer:: nr_of_groups_i_read, errval
!!    integer, dimension(100) :: subgroups_ik_list_read
!!    double precision, dimension(3):: Cji_cubic_read
!!    integer:: line_nr, n_gr, m_gr, n_gr_read, m_gr_read, leng
!!    double precision :: a_nm_read, b_nm_read, c_nm_read, a_mn_read, b_mn_read, c_mn_read
!!    double precision :: tc_est_inp, pc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_
!!
!!    !New variables for COSMO_SAC model
!!    integer, dimension(30) :: index_nr
!!    character(len=255) :: CASNR_COSMO
!!    character(4), dimension(30) :: index_char
!!    character(len=255) :: dummy_cosmo   !used to jump over the not needed lines of the inputfile
!!    character(1) :: spt     !sigma-profile type
!!    integer :: interval, profile_count, err_read
!!    integer, dimension(30,3) :: counter_min, counter_max
!!
!!    !Electrostatic interaction
!!    double precision :: sigmaacc, sigmadon              !Acceptor and Donator (hydrogen bonding interactions) in calculation of electrostatic interaction
!!    double precision, dimension(:,:),allocatable :: delta_w                         !Electrostatic interactions
!!    double precision :: eo, eps, sigma_hb, aeff
!!    double precision :: coord, c_hb, fpol, alpha, alpha_prime
!!    double precision :: c_es, c_OH_OH, c_OT_OT, c_OH_OT
!!    integer :: l
!!    !genralized equations
!!    integer:: gen_i
!!    double precision :: p_i_double
!!
!!    !open units
!!    integer:: cas_unit, fld_unit, srk_unit, cost_unit, rkm_unit, rkminfo_unit, pcsaft_unit, pr_unit, lkp_unit, gen_unit, srkmix_unit, pcsaftmix_unit, ge_unit, geuni_unit, geinter_unit, prmix_unit, lkpmix_unit
!!
!!    integer:: dummycounter
!!
!!    factor_inv = 1.d0/gl%factor
!!
!!    if (.not.allocated(CASNR_SRK)) allocate(CASNR_SRK(30))
!!    if (.not.allocated(CASNR_PR)) allocate(CASNR_PR(30))
!!    if (.not.allocated(CASNR_LKP)) allocate(CASNR_LKP(30))
!!    if (.not.allocated(CASNR_PCSAFT)) allocate(CASNR_PCSAFT(30))
!!    if (.not.allocated(CASNR_GEOS)) allocate(CASNR_GEOS(30))
!!    if (.not.allocated(a_sun1)) allocate(a_sun1(100,30))
!!    if (.not.allocated(a_sun2)) allocate(a_sun2(100,30))
!!    if (.not.allocated(delta_w)) allocate(delta_w(51,51))
!!
!!    CASNR_SRK = ''
!!
!!    !allocate ideal gas parameters
!!    !DEC$ IF DEFINED(HCgen)
!!    !DEC$ ELSE
!!    if (.not.allocated(gl%cp0coeff)) then
!!        allocate(gl%cp0coeff(20,gl%ncomp))
!!        allocate(gl%cp0exp,gl%cp0hyp1,gl%cp0hyp2,gl%cp0hyp3,mold=gl%cp0coeff)
!!    end if
!!    !allocate divers fluid names
!!    if (.not.allocated(gl%substshortname)) then
!!        allocate(gl%substshortname(gl%ncomp))
!!        allocate(gl%substfullname(gl%ncomp))
!!        allocate(gl%substchemname,gl%substsynonym,mold=gl%substfullname)
!!        allocate(gl%substcasnr,mold=gl%substshortname)
!!    end if
!!    !DEC$ ENDIF
!!
!!
!!
!!
!!    if(.not.allocated(gl%vpcoeff)) then
!!        allocate(gl%vpcoeff(20,gl%ncomp))
!!        allocate(gl%vpexp,gl%dlcoeff,gl%dlexp,gl%dvcoeff,gl%dvexp,gl%pmeltcoeff,gl%pmeltexp,gl%psubcoeff,gl%psubexp, mold=gl%vpcoeff)
!!    end if
!!
!!
!!
!!
!!
!!    !set default data:
!!    error = 0
!!    tccur = 0.d0
!!    pccur = 0.d0
!!    rhoccur = 0.d0
!!    ttpcur = 0.d0
!!    tnbpcur = 0.d0
!!    wmcur = 0.d0
!!    accencur = 0.d0
!!    line = 0
!!    new_group = .true.
!!    new_subgroup = .true.
!!    gl_aga8 = .false.
!!    groupnr = 0.d0
!!    subgroupnr = 0.d0
!!    amount = 0.d0
!!    Qkcoeffs = 0.d0
!!    lbwr = 0
!!    !nr_of_groups_mix = 0    !Reset total number of groups in mixture to 0
!!
!!    !casnr = ""
!!    !fluidname = ""
!!
!!    cas_unit =0
!!    fld_unit =0
!!    srk_unit =0
!!    cost_unit =0
!!    rkm_unit =0
!!    rkminfo_unit =0
!!    pcsaft_unit =0
!!    pr_unit =0
!!    lkp_unit =0
!!    gen_unit =0
!!    srkmix_unit =0
!!    pcsaftmix_unit =0
!!    ge_unit =0
!!    geuni_unit =0
!!    geinter_unit =0
!!    prmix_unit =0
!!    lkpmix_unit =0
!!
!!    !Format-Definitions:
!!2003 format (a3)
!!2012 format (a12)
!!2080 format (a255)
!!2100 format (a12, a20)
!!
!!    !Declarations End---------------------------------------------------------------------------------------
!!
!!    !Create pathforsingles and pathforbinarys --------------------------------------------------------------
!!    !changed for gfortran
!!    !pathforsingles = trim(path)//"fluids\"
!!    !pathforbinarys = trim(path)//"BINARY_MIX_FILES\"
!!    pathforsingles = trim(path)//"fluids/"
!!    pathforbinarys = trim(path)//"binary_mix_files/"
!!    !-------------------------------------------------------------------------------------------------------
!!
!!
!!    ! transformation from CAS to fluidname
!!    if (gl%Eq_type(1) == 1 .or. gl%EQ_type(1) == 7) then     !only neccessary for Helmholtz EOS to find the correct fld file.
!!        !check whether CAS-list can be found
!!        
!!        Check if gerg folder is specified
!!        if(gl%EQ_type(1) == 7) then
!!            if(INDEX(pathforbinarys,'gerg').eq.0) then
!!                errorflag = -7900 !no gerg folder specification found - return with error
!!            end if
!!        end if
!!
!!        CAS_list = trim (pathforsingles) // 'cas.txt'
!!        inquire (file = CAS_list, exist = CAS_list_exists)
!!
!!        ! convert all CAS IDs to fluid names
!!        do nrsubst = 1, gl%ncomp
!!            e = index (gl%components(nrsubst), '-')
!!            if (e /= 0) then    ! "-" found in compoment - =>  CAS ID - =>  convert to fluid name
!!                if (.not.CAS_list_exists) then   ! no list with CAS IDs available - >  error, abort setup
!!                    errorflag = -7884
!!                    return
!!                else    ! CAS IDs are available - =>  try to convert to fld
!!                    open(newunit = cas_unit, file = CAS_list, status='old', action='read', iostat=errorflag)
!!                    if (errorflag /= 0) return
!!
!!                    e=0
!!                    do
!!                        read(cas_unit,*,end=998) dummy
!!                        e = e + 1
!!                    enddo
!!998                 rewind(cas_unit)
!!                    if(.not. allocated(dummy_arr)) allocate(dummy_arr(2,e))
!!
!!                    if(.not. allocated(casnr)) allocate(casnr(e))
!!                    if(.not. allocated(fluidname)) allocate(fluidname(e))
!!                    read (cas_unit, 2100) dummy_arr
!!                    casnr = dummy_arr(1,:)
!!                    fluidname = dummy_arr(2,:)
!!                    if (any(casnr(:) == gl%components (nrsubst))) then
!!                        !WORKS ONLY IF IVF 2018 or newer is installed
!!                        !gl%components(nrsubst) = fluidname(findloc(casnr(:), gl%components (nrsubst),1))
!!                        !THIS IS THE OLD FORTRAN-ish STUPID WAY
!!                        do iter  = 1,e
!!                            if (casnr(iter) == gl%components(nrsubst)) then
!!                                gl%components(nrsubst) = fluidname(iter)
!!                                exit
!!                            endif
!!                        enddo
!!                    else
!!                        errorflag = -7890  ! CAS ID not found until end of list
!!                    endif
!!                    !                        if (casnr == gl%components (nrsubst)) then   ! CAS ID is found
!!                    !                            gl%components (nrsubst) = fluidname
!!                    !                            filename = trim (pathforsingles) // trim (fluidname) // '.fld'
!!                    !                            exit
!!                    !                        end if
!!                    !                        number = 0
!!                    !                        read (gl%components(nrsubst)(e+1:e+1),*,iostat=error) number
!!                    !                        if (error /= 0) then
!!                    !                            errorflag = -7879  ! CAS ID not found until end of list
!!                    !                        else
!!                    !999                         errorflag = -7890  ! CAS ID not found until end of list
!!                    !                        end if
!!                    !                        close (cas_unit)
!!                    close (cas_unit)
!!                end if
!!            end if
!!        end do
!!    end if
!!
!!    loopsubst: do nrsubst = 1, gl%ncomp
!!
!!        if (gl%Eq_type(nrsubst) == 1 .or. gl%Eq_type(nrsubst) == 7) then     !Helmholtz EOS. Continue reading the fluid file
!!
!!            if (gl%components(nrsubst) == '') then
!!                exit loopsubst
!!            else
!!
!!                feq_read = .false.   !eq has not been read, yet
!!                cp0_read = .false.   !cp0 has not been read, yet
!!                cpp_read = .false.
!!                phk_read = .false.
!!                pho_read = .false.
!!                gl%bwr_read = .false.
!!                gl%ben_read = .false.
!!                gl%sta_read = .false.
!!                gl%stn_read = .false.
!!
!!
!!                filename = trim (pathforsingles) // trim (gl%components(nrsubst)) // '.fld'  ! check for standard case: fluid file (fld)
!!                inquire (file = filename, exist = exists)
!!                if (.not.exists) then   ! no suitable fld found
!!                    filename = trim (pathforsingles) // trim (gl%components(nrsubst)) // '.ppf' ! check for pseudo pure fluid (ppf)
!!                    inquire (file = filename, exist = exists)
!!                end if
!!                if (.not.exists) then   ! no suitable ppf found
!!                    errorflag = -7878
!!                    return
!!                end if
!!
!!
!!2908            open(newunit=fld_unit, file=trim(filename), status='old', action='read', iostat=error, share='denynone')
!!
!!                !If no error occurs read file
!!                if (error == 0) then
!!
!!                    !!!if(.not.feq_read) then
!!
!!                    !first paragraph till UN-number:
!!                    !1. Name of Substance:
!!                    read(fld_unit,2012) gl%substshortname(nrsubst)
!!                    read(fld_unit,2012) gl%substcasnr(nrsubst)
!!                    read(fld_unit,2080) gl%substfullname(nrsubst)
!!                    read(fld_unit,2080) gl%substchemname(nrsubst)
!!                    read(fld_unit,2080) gl%substsynonym(nrsubst)
!!
!!                    !2. Current Properties of Substance:
!!                    read(fld_unit,*) wmcur(nrsubst)             !current molecular weight [g/mol]
!!                    read(fld_unit,*) ttpcur(nrsubst)            !current temperature at triple point [K]
!!                    read(fld_unit,*) tnbpcur(nrsubst)           !current temperature at normal boiling point [K]
!!                    read(fld_unit,*) tccur(nrsubst)             !current temperature at critical point [K]
!!                    read(fld_unit,*) pccur(nrsubst)             !current pressure at critical point [kPa]
!!                    pccur(nrsubst) = pccur(nrsubst)*factor_inv   ![MPa]
!!                    read(fld_unit,*) rhoccur(nrsubst)           !current density at critical point [mol/L]
!!                    rhoccur(nrsubst) = rhoccur(nrsubst)*gl%Factor ![mol/m設
!!                    read(fld_unit,*) accencur(nrsubst)          !current acentric factor [-]
!!                    read(fld_unit,*) gl%dipole(nrsubst)            !dipole moment [debye]
!!
!!                    !3. Reference State:
!!                    read(fld_unit,2003) gl%refstate(nrsubst)           !default reference state
!!                    if (gl%refstate(nrsubst)(1:2) == 'OT') then !if 'OTher' reference state then read in reference T,p,h,s
!!                        read(fld_unit,*)  gl%tref(nrsubst),gl%pref(nrsubst),gl%href(nrsubst),gl%sref(nrsubst)
!!                        gl%pref(nrsubst) = gl%pref(nrsubst)*factor_inv        ![MPa]
!!                    elseif (gl%refstate(nrsubst)(1:3) == 'NBP') then
!!                        gl%tref(nrsubst) = 0.d0 !later in this subroutine set to normal boiling point temperature
!!                        gl%pref(nrsubst) = 0.101325d0 !always 1 atm = 0.101325 MPa
!!                        gl%href(nrsubst) = 0.d0 ![J/mol]
!!                        gl%sref(nrsubst) = 0.d0 ![J/(mol K)]
!!                    elseif (gl%refstate(nrsubst)(1:3) == 'ASH') then
!!                        gl%tref(nrsubst) = 233.15d0
!!                        gl%pref(nrsubst) = -1.d0 !"-1" shows that pref must be calculated later in another subroutine
!!                        gl%href(nrsubst) = 0.d0 ![J/mol]
!!                        gl%sref(nrsubst) = 0.d0 ![J/(mol K)]
!!                    elseif (gl%refstate(nrsubst)(1:3) == 'IIR') then !Pay attention when using Lennard Jones Units
!!                        gl%tref(nrsubst) = 273.15d0
!!                        gl%pref(nrsubst) = -1.d0  !"-1" shows that pref must be calculated later in another subroutine
!!                        gl%href(nrsubst) = 200.d3 ![J/kg]     -later changed with wm into [J/mol]
!!                        gl%sref(nrsubst) = 1.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
!!                    else !set to default values
!!                        gl%tref(nrsubst) = 0.d0
!!                        gl%pref(nrsubst) = 0.d0
!!                        gl%href(nrsubst) = 0.d0 ![J/mol]
!!                        gl%sref(nrsubst) = 0.d0 ![J/(mol K)]
!!                    endif
!!
!!                    !In the fluids files there is additional information given (e.g., version number, family, etc.).
!!                    !These information is not needed here so that it is not read in.
!!
!!                    do
!!                        read(fld_unit,*) dummy
!!                        if ((dummy == '#EOS').OR.(dummy == '@EOS').OR.(dummy == '@FEK').OR.(dummy == '#FEK').OR.(dummy == '#AUX').OR.(dummy == '#ASO')  &
!!                            .OR.(dummy == '#DE').OR.(dummy == '#STN').OR.(dummy == '#FE1').OR.(dummy == '@FE1')) then
!!                            backspace(fld_unit)
!!                            exit
!!                        end if
!!                    end do
!!
!!
!!
!!                    !Loop to read whole Fluid-file
!!1337                loopread:    do
!!                        read(fld_unit,*) dummygen
!!                        if (dummygen(1:4) == '#EOS') then
!!                            !second paragraph with comments can be ignored (every line starts with "!")
!!                            goover1:    do !until dummy starts with "#" or "@"
!!                                !if next paragraph starts with "#EOS" then check next line for "FEQ"
!!                                read(fld_unit,*) dummy !=Equation-Type - can be one of the following:
!!                                !         'NBS' = NIST-recommended model
!!                                !         'BWR' = modified Benedict-Webb-Rubin equation of state
!!                                !         'FEQ' = fundamental (Helmholtz) equation of state
!!                                !         'FEH' = fundamental (Helmholtz) equation of state with hard shpere term
!!                                !         'ECS' = extended corresponding states model
!!                                !         'HMX' = mixture Helmholtz model
!!                                !         'CPx' = polynomial fit of ideal gas heat capacity
!!                                !         'VSi' = pure fluid viscosity model #i
!!                                !         'TCi' = pure fluid thermal conductivity model #i
!!                                exit goover1
!!                            enddo goover1
!!
!!                            !if Equation-Type is FEQ then there is another comment paragraph:
!!                            if ((dummy(1:2) == 'FE')) then
!!                                if (dummy(1:3) == 'FEH') gl%eos_coeff%hard_sphere(nrsubst) = .true.
!!                                !third paragraph with comments can be ignored (every line starts with "?")
!!                                goover2:        do !until dummy starts with "!"
!!                                    read(fld_unit,*) dummy
!!                                    call uppertolower_char(dummy, len(dummy))
!!                                    if (dummy(1:1) == '!') then
!!                                        exit goover2
!!                                    elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                                        read(fld_unit,'(A)')dummy
!!                                        do while (trim(dummy(1:5)) /= '!````')
!!                                            gl%litref%lit_ref_res(nrsubst) = trim(gl%litref%lit_ref_res(nrsubst)) // dummy(2:)
!!                                            read(fld_unit,'(A)')dummy
!!                                        end do
!!                                        backspace(fld_unit)
!!                                    endif
!!                                enddo goover2
!!
!!                                !data paragraph
!!                                read(fld_unit,*) gl%tminfluid(nrsubst)          !lower temperature limit [K]
!!                                read(fld_unit,*) gl%tmaxfluid(nrsubst)          !upper temperature limit [K]
!!                                read(fld_unit,*) gl%pmaxfluid(nrsubst)          !upper pressure limit [kPa]
!!                                gl%pmaxfluid(nrsubst) = gl%pmaxfluid(nrsubst)*factor_inv    ![MPa]
!!                                read(fld_unit,*) gl%rhomaxfluid(nrsubst)        !upper density limit [mol/L]
!!                                gl%rhomaxfluid(nrsubst) = gl%rhomaxfluid(nrsubst)*gl%Factor ![mol/m設
!!                                read(fld_unit,2003) gl%eos_coeff%cppcheck(nrsubst)        !pointer to Cp0 model
!!                                read(fld_unit,*) gl%wm(nrsubst)                 !molecular weight [g/mol]
!!                                if (gl%refstate(nrsubst)(1:3) == 'IIR') then
!!                                    gl%href(nrsubst) = gl%href(nrsubst)*gl%wm(nrsubst)*factor_inv ![J/mol]
!!                                    gl%sref(nrsubst) = gl%sref(nrsubst)*gl%wm(nrsubst)*factor_inv ![J/(mol*K)]
!!                                endif
!!                                read(fld_unit,*) gl%ttp(nrsubst)          !triple point temperature
!!                                if ((.not. gl%hold_limits).and.(.not. gl%check_solid)) gl%ttp(nrsubst) = gl%ttp(nrsubst)*0.1d0    !small triple point temperature for the extrapolation of the VLE curves
!!                                !     (if "&" occurs in the input variable)
!!                                read(fld_unit,*) gl%ptp(nrsubst)          !pressure at triple point
!!                                gl%ptp(nrsubst) = gl%ptp(nrsubst)*factor_inv  ![MPa]
!!                                read(fld_unit,*) gl%rhotp(nrsubst)        !density at triple point
!!                                gl%rhotp(nrsubst) = gl%rhotp(nrsubst)*gl%Factor ![mol/m設
!!                                read(fld_unit,*) gl%tnbp(nrsubst)         !temperature at normal boiling point
!!                                if (gl%refstate(nrsubst)(1:3) == 'NBP') then
!!                                    gl%tref(nrsubst) = gl%tnbp(nrsubst) !could not be set to normal boiling point temperature before
!!                                endif
!!                                read(fld_unit,*) gl%accen(nrsubst)        !acentric factor
!!                                read(fld_unit,*) gl%tc(nrsubst),gl%pc(nrsubst),gl%rhoc(nrsubst)   !critical parameters
!!                                gl%pc(nrsubst) = gl%pc(nrsubst)*factor_inv       ![MPa]
!!                                gl%rhoc(nrsubst) = gl%rhoc(nrsubst)*gl%Factor   ![mol/m設
!!                                read(fld_unit,*) gl%tred(nrsubst),gl%rhored(nrsubst)  !reducing parameters
!!                                gl%rhored(nrsubst) = gl%rhored(nrsubst)*gl%Factor ![mol/m設
!!                                read(fld_unit,*) gl%Req(nrsubst)          !gas constant used in fit
!!                                !data columns
!!                                err_read = 0
!!                                read(fld_unit,*,iostat=err_read) gl%eos_coeff%nreg(nrsubst),gl%eos_coeff%nlreg(nrsubst), gl%eos_coeff%ncrt(nrsubst),gl%eos_coeff%nlcrt(nrsubst), gl%eos_coeff%nna(nrsubst), dummycounter, gl%eos_coeff%nlna(nrsubst)
!!                                if (err_read /= 0) then
!!                                    backspace(fld_unit)
!!                                    read(fld_unit,*,iostat=err_read) gl%eos_coeff%nreg(nrsubst),gl%eos_coeff%nlreg(nrsubst), gl%eos_coeff%ncrt(nrsubst),gl%eos_coeff%nlcrt(nrsubst), gl%eos_coeff%nna(nrsubst), gl%eos_coeff%nlna(nrsubst)
!!                                endif
!!                                !read(fld_unit,*) nreg(nrsubst),nlreg(nrsubst), ncrt(nrsubst),nlcrt(nrsubst), assomodel(nrsubst),nlna(nrsubst)
!!
!!                                !  the gamma term is a multiplier for the (rho or del) in only the exponential
!!                                !  terms; it is needed for e.g. the Bender EOS; set to 1.0 if not present
!!
!!                                do j=1,gl%eos_coeff%nreg(nrsubst) !for normal terms
!!                                    p_i_double=0d0
!!                                    if (gl%eos_coeff%nlreg(nrsubst) == 5) then ! if 5 columns are stated for regular terms, then gama = 1 (except for BWR equations)
!!                                        read(fld_unit,*) gl%eos_coeff%ni(j,nrsubst),gl%eos_coeff%ti(j,nrsubst),gl%eos_coeff%di(j,nrsubst),p_i_double,gl%eos_coeff%gama(j,nrsubst)
!!                                        gl%eos_coeff%p_i(j,nrsubst) = int(p_i_double)
!!                                    else
!!                                        read(fld_unit,*) gl%eos_coeff%ni(j,nrsubst),gl%eos_coeff%ti(j,nrsubst),gl%eos_coeff%di(j,nrsubst),p_i_double
!!                                        gl%eos_coeff%p_i(j,nrsubst) = int(p_i_double)
!!                                        if (gl%eos_coeff%p_i(j,nrsubst) /= 0) gl%eos_coeff%gama(j,nrsubst) = 1.d0
!!                                    endif
!!                                enddo
!!
!!                                if (gl%eos_coeff%ncrt(nrsubst) > 0) then
!!                                    do j=1,gl%eos_coeff%ncrt(nrsubst) !for critical terms
!!                                        read(fld_unit,*) gl%eos_coeff%ni(j+gl%eos_coeff%nreg(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%ti(j+gl%eos_coeff%nreg(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%di(j+gl%eos_coeff%nreg(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%pli(j,nrsubst), &
!!                                            gl%eos_coeff%tli(j,nrsubst), &
!!                                            gl%eos_coeff%eta(j,nrsubst), &
!!                                            gl%eos_coeff%beta(j,nrsubst), &
!!                                            gl%eos_coeff%gam(j,nrsubst), &
!!                                            gl%eos_coeff%eps(j,nrsubst), &
!!                                            gl%eos_coeff%etana(j,nrsubst), &!!
!!                                            gl%eos_coeff%eidna(j,nrsubst), &!!
!!                                            gl%eos_coeff%eitna(j,nrsubst)!!
!!                                    enddo
!!                                endif
!!
!!                                if (gl%eos_coeff%nna(nrsubst) > 0) then
!!                                    do j=1,gl%eos_coeff%nna(nrsubst) !for special gaussian terms
!!                                        read(fld_unit,*) gl%eos_coeff%ni(j+gl%eos_coeff%nreg(nrsubst)+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%ti(j+gl%eos_coeff%nreg(nrsubst)+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%di(j+gl%eos_coeff%nreg(nrsubst)+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%pli(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%tli(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%eta(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%beta(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%gam(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%eps(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &
!!                                            gl%eos_coeff%etana(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &!!
!!                                            gl%eos_coeff%eidna(j+gl%eos_coeff%ncrt(nrsubst),nrsubst), &!!
!!                                            gl%eos_coeff%eitna(j+gl%eos_coeff%ncrt(nrsubst),nrsubst)!!
!!                                    enddo
!!                                endif
!!
!!
!!                                tc_est_inp = gl%tc(nrsubst)
!!                                rhoc_est_inp = gl%rhoc(nrsubst)
!!                                gl%rhoredmix = 1.d0
!!
!!                                !! find the number of polynomial terms in the equation
!!                                !do ln = 1, nreg(nrsubst)
!!                                !    if (p_i(nrsubst, ln) == 0.d0) I_pol(nrsubst) = I_pol(nrsubst) + 1
!!                                !end do
!!                                ! find the number of gaussian bell-shaped terms in the equation
!!                                gl%eos_coeff%I_GBS(nrsubst) = 0
!!                                do m = 1, gl%eos_coeff%ncrt(nrsubst)
!!                                    if (gl%eos_coeff%etana(m,nrsubst) == 0.d0) gl%eos_coeff%I_GBS(nrsubst) = gl%eos_coeff%I_GBS(nrsubst) + 1
!!                                end do
!!                                !I_exp(nrsubst) = nreg(nrsubst) - I_pol(nrsubst)    ! number of exponential terms
!!                                !I_NA(nrsubst) = ncrt(nrsubst) - I_GBS(nrsubst)     ! number of nonanalytic terms
!!
!!
!!                                !packing fraction parameters for hard sphere terms
!!                                if (gl%eos_coeff%hard_sphere(nrsubst)) then
!!                                    read (fld_unit,*) gl%eos_coeff%h1(nrsubst), gl%eos_coeff%h2(nrsubst), gl%eos_coeff%h3(nrsubst), gl%eos_coeff%h4(nrsubst), gl%eos_coeff%dof2(nrsubst)
!!                                end if
!!
!!                                Goto 1337
!!
!!
!!                            elseif (dummy(1:3) == 'ECS') then
!!                                errorflag = -7999	!chosen file has wrong format
!!                                close(fld_unit)
!!                                return
!!                            endif
!!
!!                        elseif ((dummy(1:3) == 'BWR') .or. (dummy(1:3) == 'BEN') .or. (dummy(1:3) == 'STA'))  then
!!                            if (dummy(1:3) == 'BWR') then
!!                                gl%bwr_read = .true.
!!                                loop = 11
!!                            elseif (dummy(1:3) == 'BEN') then
!!                                gl%ben_read = .true.
!!                                loop = 7
!!                            elseif (dummy(1:3) == 'STA') then
!!                                gl%sta_read = .true.
!!                                loop = 5
!!                            endif
!!                            loopBWR:        do
!!
!!                                gooverBWR:         do !until dummy starts with "!"
!!                                    read(fld_unit,*) dummy
!!                                    if (dummy(1:1) == '!') then
!!                                        exit gooverBWR
!!                                    endif
!!                                end do gooverBWR
!!                                !reading Parameters
!!                                !data paragraph
!!                                read(fld_unit,*) gl%tminfluid(nrsubst)                            !lower temperature limit [K]
!!                                read(fld_unit,*) gl%tmaxfluid(nrsubst)                            !upper temperature limit [K]
!!                                read(fld_unit,*) gl%pmaxfluid(nrsubst)                            !upper pressure limit [kPa]
!!                                gl%pmaxfluid(nrsubst) = gl%pmaxfluid(nrsubst)*factor_inv      ![MPa]
!!                                read(fld_unit,*) gl%rhomaxfluid(nrsubst)                          !upper density limit [mol/L]
!!                                gl%rhomaxfluid(nrsubst) = gl%rhomaxfluid(nrsubst)*gl%Factor  ![mol/m設
!!                                read(fld_unit,*) gl%eos_coeff%cppcheck(nrsubst)                             !pointer to Cp0 model
!!                                read(fld_unit,*) gl%wm(nrsubst)                                   !molecular weight [g/mol]
!!                                read(fld_unit,*) gl%ttp(nrsubst)                                  !triple point temperature
!!                                if (.not. gl%hold_limits) gl%ttp(nrsubst) = gl%ttp(nrsubst)*0.1d0    !small triple point temperature for the extrapolation of the VLE curves
!!                                !     (if "&" occurs in the input variable)
!!                                read(fld_unit,*) gl%ptp(nrsubst)                                  !pressure at triple point
!!                                gl%ptp(nrsubst) = gl%ptp(nrsubst)*factor_inv                  ![MPa]
!!                                read(fld_unit,*) gl%rhotp(nrsubst)                                !density at triple point
!!                                gl%rhotp(nrsubst) = gl%rhotp(nrsubst)*gl%Factor              ![mol/m設
!!                                read(fld_unit,*) gl%tnbp(nrsubst)                                 !temperature at normal boiling point
!!                                read(fld_unit,*) gl%accen(nrsubst)                                !acentric factor
!!                                read(fld_unit,*) gl%tc(nrsubst),gl%pc(nrsubst),gl%rhoc(nrsubst)         !critical parameters
!!                                gl%pc(nrsubst) = gl%pc(nrsubst)*factor_inv                    ![MPa]
!!                                gl%rhoc(nrsubst) = gl%rhoc(nrsubst)*gl%Factor                ![mol/m設
!!                                read(fld_unit,*) gl%tred(nrsubst),gl%rhored(nrsubst)                 !reducing parameters
!!                                gl%rhored(nrsubst) = gl%rhored(nrsubst)*gl%Factor            ![mol/m設
!!                                read(fld_unit,*) gl%gama_bwr(nrsubst)
!!                                read(fld_unit,*) gl%Req(nrsubst)                                  !gas constant used in fit [L-bar/mol-K]
!!                                gl%Req(nrsubst) = gl%Req(nrsubst)*gl%factorrbwr                  ![J/mol-K]
!!                                read(fld_unit,*) gl%eos_coeff%nreg(nrsubst)                                 !number of terms in BWR equation
!!
!!                                do ibwr=1, loop
!!                                    if (lbwr+ibwr < (gl%eos_coeff%nreg(nrsubst)-1.d0)) then
!!                                        read(fld_unit,*) gl%ncoeff(nrsubst,ibwr+lbwr), gl%ncoeff(nrsubst, (ibwr+lbwr+1)), gl%ncoeff(nrsubst, (ibwr+lbwr+2))
!!                                    else
!!                                        if (gl%bwr_read) then
!!                                            read(fld_unit,*)  gl%ncoeff(nrsubst,ibwr+lbwr), gl%ncoeff(nrsubst, (ibwr+lbwr+1))
!!                                        elseif (gl%ben_read) then
!!                                            read(fld_unit,*)  gl%ncoeff(nrsubst,ibwr+lbwr)
!!                                        endif
!!                                    endif
!!                                    lbwr = lbwr+2
!!                                enddo
!!                                call bwr_recalc(gl,nrsubst)
!!
!!                                !do j=1,nreg(nrsubst)
!!                                !    if (p_i(nrsubst,j) == 0.d0) then
!!                                !        I_pol(nrsubst) = I_pol(nrsubst) + 1               ! number of polynomial terms
!!                                !    end if
!!                                !end do
!!                                !I_exp(nrsubst) = nreg(nrsubst) - I_pol(nrsubst)    ! number of exponential terms
!!                                exit loopBWR
!!                            enddo loopBWR
!!
!!
!!
!!
!!                            !elseif (dummy(1:3) == 'ECS') then
!!                            !        errorflag = -7999	!chosen file has wrong format
!!                            !     close(unit=8)
!!                            !     return
!!                            !endif
!!
!!                            Goto 1337
!!
!!                        elseif ((dummygen(1:4) == '#AUX')) then !only the first #AUX, so don't mix with next do loop
!!
!!                            loopaux:    do
!!                                !if (dummy(1:4) == '#AUX') then !only the first #AUX, so don't mix with next do loop
!!                                read(fld_unit,*) dummy
!!                                if ((dummy(1:2) == 'CJ') .and. (.not.cpp_read)) then
!!                                    gl%cpmodel(nrsubst) = .true.
!!                                    do !until dummy starts with "!"
!!                                        read(fld_unit,*) dummy
!!                                        if (dummy(1:1) == '!') then
!!                                            exit
!!                                        endif
!!                                    enddo
!!                                    !paragraph with doubled parameters can be ignored (always 4 lines)
!!                                    do j=1,4
!!                                        read(fld_unit,*) dummy
!!                                    enddo
!!                                    !read reducing parameters for cp0:
!!                                    read(fld_unit,*) gl%tcp0red(nrsubst), gl%cp0red(nrsubst)
!!                                    !read coeffs ande exps:
!!                                    read(fld_unit,*) gl%ncp0poly(nrsubst),gl%ncp0pl(nrsubst), gl%ncp0cosh(nrsubst),gl%ncp0sinh(nrsubst)
!!                                    do j=1,(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst))
!!                                        read(fld_unit,*) gl%cp0coeff(j,nrsubst),gl%cp0exp(j,nrsubst)
!!                                    enddo
!!
!!                                    cpdummy=gl%eos_coeff%cppcheck(nrsubst)
!!                                    cpdummy(1:2) = 'CP'
!!                                    gl%cp0red(nrsubst) = 8.3144621d0
!!                                    gl%CP0COEFF(1,nrsubst) = (gl%cp0coeff(1,nrsubst)-37.93d0)/gl%cp0red(nrsubst)
!!                                    gl%CP0COEFF(2,nrsubst) = (gl%cp0coeff(2,nrsubst)+0.21d0)/gl%cp0red(nrsubst)
!!                                    gl%CP0COEFF(3,nrsubst) = (gl%cp0coeff(3,nrsubst)-3.91d-4)/gl%cp0red(nrsubst)
!!                                    gl%CP0COEFF(4,nrsubst) = (gl%cp0coeff(4,nrsubst)+2.06d-7)/gl%cp0red(nrsubst)
!!
!!
!!                                end if
!!
!!                                if ((dummy(1:2) == 'CP') .and. (.not.cpp_read)) then
!!                                    gl%cpmodel(nrsubst) = .true.
!!                                    !paragraph with comments can be ignored (every line starts with "?")
!!                                    goover3:                do !until dummy starts with "!"
!!                                        read(fld_unit,*) dummy
!!                                        if (dummy(1:1) == '!') then
!!                                            exit goover3
!!                                        endif
!!                                    enddo goover3
!!                                    !paragraph with doubled parameters can be ignored (always 4 lines)
!!                                    do j=1,4
!!                                        read(fld_unit,*) dummy
!!                                    enddo
!!                                    !read reducing parameters for cp0:
!!                                    read(fld_unit,*) gl%tcp0red(nrsubst), gl%cp0red(nrsubst)
!!                                    !read coeffs ande exps:
!!                                    read(fld_unit,*) gl%ncp0poly(nrsubst),gl%ncp0pl(nrsubst), gl%ncp0cosh(nrsubst),gl%ncp0sinh(nrsubst)
!!                                    do j=1,(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst))
!!                                        read(fld_unit,*) gl%cp0coeff(j,nrsubst),gl%cp0exp(j,nrsubst)
!!                                    enddo
!!                                    if ((gl%ncp0cosh(nrsubst)+gl%ncp0sinh(nrsubst)) > 0) then
!!                                        do j=1,(gl%ncp0cosh(nrsubst)+gl%ncp0sinh(nrsubst)) !for hyperbolic
!!                                            read(fld_unit,*) gl%cp0coeff(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst), &
!!                                                gl%cp0exp(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst), &
!!                                                gl%cp0hyp1(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst), &
!!                                                gl%cp0hyp2(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst), &
!!                                                gl%cp0hyp3(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst)
!!                                        enddo
!!
!!                                    endif
!!                                    cpp_read = .true.
!!                                    exit loopaux
!!                                elseif ((dummy(1:3) == 'PHK') .AND. (gl%eos_coeff%cppcheck(nrsubst) == 'PHK') .and. (.not.phk_read)) then !Modified by Andreas July 2011 , Modified by Judith August 2011
!!                                    gl%phkmodel(nrsubst) = .true.
!!                                    gl%cpmodel(nrsubst) = .false.
!!                                    goover25:               do !until dummy starts with "!"
!!                                        read(fld_unit,*) dummy
!!                                        if (dummy(1:1) == '!') then
!!                                            exit goover25
!!                                        endif
!!                                    end do goover25
!!                                    !paragraph with doubled parameters can be ignored (always 4 lines)
!!                                    do j=1,4
!!                                        read(fld_unit,*) dummy
!!                                    end do
!!                                    !set reference state and reducing parameters for cp0:
!!                                    gl%refstate(nrsubst)='OT0'
!!                                    gl%tref(nrsubst) = 298.15d0
!!                                    gl%pref(nrsubst) = 0.101325d0  !"-1" shows that pref must be calculated later in another subroutine
!!                                    gl%href(nrsubst) = 0.d0 ![J/kg]     -later changed with wm into [J/mol]
!!                                    gl%sref(nrsubst) = 0.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
!!                                    gl%tcp0red(nrsubst) = gl%tc(nrsubst)
!!                                    gl%cp0red(nrsubst) = 8.31451D0
!!                                    !read coeffs and exps:
!!                                    read(fld_unit,*) gl%ncp0log(nrsubst),gl%ncp0(nrsubst),gl%ncp0logexp(nrsubst), gl%ncp0cosh(nrsubst),gl%ncp0sinh(nrsubst)
!!                                    !Read polynomials and Constants CI and CII
!!                                    do j=1,(gl%ncp0log(nrsubst)+gl%ncp0(nrsubst)+gl%ncp0logexp(nrsubst))
!!                                        read(fld_unit,*) gl%cp0coeff(j,nrsubst),gl%cp0exp(j,nrsubst)
!!                                    end do
!!                                    if ((gl%ncp0cosh(nrsubst)+gl%ncp0sinh(nrsubst)) > 0) then
!!                                        do j=1,(gl%ncp0cosh(nrsubst)+gl%ncp0sinh(nrsubst)) !for hyperbolic terms
!!                                            read(fld_unit,*) gl%cp0coeff(j+(gl%ncp0log(nrsubst)+gl%ncp0(nrsubst)+gl%ncp0logexp(nrsubst)),nrsubst), &
!!                                                gl%cp0exp(j+(gl%ncp0log(nrsubst)+gl%ncp0(nrsubst)+gl%ncp0logexp(nrsubst)),nrsubst)
!!                                        end do
!!                                    end if
!!                                    phk_read = .true.
!!                                    exit loopaux
!!                                elseif ((dummy(1:3) == 'PH0') .AND. (gl%eos_coeff%cppcheck(nrsubst) == 'PH0').and. (.not.pho_read)) then
!!                                    gl%phmodel(nrsubst) = .true.
!!                                    goover26:               do !until dummy starts with "!"
!!                                        read(fld_unit,*) dummy
!!                                        if (dummy(1:1) == '!') then
!!                                            exit goover26
!!                                        endif
!!                                    end do goover26
!!                                    !paragraph with doubled parameters can be ignored (always 4 lines)
!!                                    do j=1,4
!!                                        read(fld_unit,*) dummy
!!                                    end do
!!                                    !set reducing parameters for cp0:
!!                                    gl%tcp0red(nrsubst) = gl%tc(nrsubst)
!!                                    gl%cp0red(nrsubst) = gl%REQ(nrsubst)
!!                                    !read coeffs and exps:
!!                                    read(fld_unit,*) gl%ncp0c(nrsubst),gl%ncp0poly(nrsubst),gl%ncp0pl(nrsubst), gl%ncp0cosh(nrsubst),gl%ncp0sinh(nrsubst)
!!                                    gl%ncp0poly(nrsubst) = gl%ncp0poly(nrsubst)+gl%ncp0c(nrsubst)
!!                                    !Read polynomials and Constants CI and CII
!!                                    do j=1,(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst))
!!                                        read(fld_unit,*) gl%cp0coeff(j,nrsubst),gl%cp0exp(j,nrsubst)
!!                                    end do
!!                                    if ((gl%ncp0cosh(nrsubst)+gl%ncp0sinh(nrsubst)) > 0) then
!!                                        do j=1,(gl%ncp0cosh(nrsubst)+gl%ncp0sinh(nrsubst)) !for hyperbolic terms
!!                                            read(fld_unit,*) gl%cp0coeff(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst), &
!!                                                gl%cp0exp(j+(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst)),nrsubst)
!!                                        end do
!!                                    end if
!!                                    pho_read = .true.
!!                                    exit loopaux
!!                                elseif ((dummy(1:3) == 'H2O').or.(dummy(1:3) == 'I08').or.(dummy(1:3) == 'CI1').or.(dummy(1:3) == 'CI2')) then
!!                                    if (.not.allocated(gl%visco)) then
!!                                        exit loopaux
!!                                    endif
!!                                    if ((dummy(1:3) == 'H2O') .AND. (gl%visco%pointer_lambda(nrsubst) == 'H2O')) then
!!
!!                                        ! elseif ((dummy(1:3) == 'H2O') .AND. (gl%visco%pointer_lambda(nrsubst) == 'H2O')) then
!!                                        goover27:              do !until dummy starts with "!"
!!                                            read(fld_unit,*) dummy
!!                                            if (dummy(1:1) == '!') then
!!                                                exit goover27
!!                                            endif
!!                                        end do goover27
!!                                        read(fld_unit,*) gl%visco%lowbo_temp_mue2(nrsubst)       !lower boundary [K] of region around the crit. Temp.
!!                                        read(fld_unit,*) gl%visco%uppbo_temp_mue2(nrsubst)       !upper boundary [K] of region around the crit. Temp.
!!                                        read(fld_unit,*) gl%visco%lowbo_dens_mue2(nrsubst)       !lower boundary	[kg/m^3] of region around the crit. Dens.
!!                                        read(fld_unit,*) gl%visco%uppbo_dens_mue2(nrsubst)       !upper boundary	[kg/m^3] of region around the crit. Dens.
!!                                        read(fld_unit,*) gl%visco%corr_ref_mue2(nrsubst)         !reference correlation length
!!                                        read(fld_unit,*) gl%visco%X_mue2(nrsubst)                !x_mue
!!                                        read(fld_unit,*) qc_inv(nrsubst)               !q_c^-1  [1/nm]
!!                                        gl%visco%qc_mue2(nrsubst)=1.d0/(qc_inv(nrsubst)*1.d-9)    ![m]
!!                                        read(fld_unit,*) qd_inv(nrsubst)               !q_d^-1  [1/nm]
!!                                        gl%visco%qd_mue2(nrsubst)=1.d0/(qd_inv(nrsubst)*1.d-9)    ![m]
!!                                        read(fld_unit,*) gl%visco%nue_mue2(nrsubst)              !nue
!!                                        read(fld_unit,*) gl%visco%gamma_mue2(nrsubst)            !gamma
!!                                        read(fld_unit,*) gl%visco%Xi0_mue2(nrsubst)              !Xi [nm]
!!                                        gl%visco%Xi0_mue2(nrsubst)=gl%visco%Xi0_mue2(nrsubst)*1.d-9  ![m]
!!                                        read(fld_unit,*) gl%visco%gamma0_mue2(nrsubst)           !GAMMA_0
!!                                        read(fld_unit,*) gl%visco%T_mue2(nrsubst)                !T_quer_R
!!                                        exit loopaux
!!
!!
!!                                        ! elseif ((dummy(1:3) == 'I08') .AND. allocated(gl%visco)) then
!!                                    elseif ((dummy(1:3) == 'I08') .AND. (gl%visco%pointerceaux(nrsubst) == 'I08')) then
!!                                        goover31:               do !until dummy starts with "!"
!!                                            read(fld_unit,*) dummy
!!                                            if (dummy(1:1) == '!') then
!!                                                exit goover31
!!                                            endif
!!                                        end do goover31
!!                                        read(fld_unit,*) gl%visco%lowbo_temp_mue2(nrsubst)       !lower boundary [K] of region around the crit. Temp.
!!                                        read(fld_unit,*) gl%visco%uppbo_temp_mue2(nrsubst)       !upper boundary [K] of region around the crit. Temp.
!!                                        read(fld_unit,*) gl%visco%lowbo_dens_mue2(nrsubst)       !lower boundary	[kg/m^3] of region around the crit. Dens.
!!                                        read(fld_unit,*) gl%visco%uppbo_dens_mue2(nrsubst)       !upper boundary	[kg/m^3] of region around the crit. Dens.
!!                                        read(fld_unit,*) gl%visco%corr_ref_mue2(nrsubst)         !reference correlation length
!!                                        read(fld_unit,*) gl%visco%X_mue2(nrsubst)                !x_mue
!!                                        read(fld_unit,*) gl%visco%qc_mue2(nrsubst)               !q_c^-1
!!                                        read(fld_unit,*) gl%visco%qd_mue2(nrsubst)               !q_d^-1
!!                                        read(fld_unit,*) gl%visco%nue_mue2(nrsubst)              !nue
!!                                        read(fld_unit,*) gl%visco%gamma_mue2(nrsubst)            !gamma
!!                                        read(fld_unit,*) gl%visco%Xi0_mue2(nrsubst)              !Xi
!!                                        read(fld_unit,*) gl%visco%gamma0_mue2(nrsubst)           !GAMMA_0
!!                                        read(fld_unit,*) gl%visco%T_mue2(nrsubst)                !T_quer_R
!!                                        gl%visco%lowbo_dens_mue2(nrsubst)=gl%visco%lowbo_dens_mue2(nrsubst)*gl%Factor/gl%wm(nrsubst)
!!                                        gl%visco%uppbo_dens_mue2(nrsubst)=gl%visco%uppbo_dens_mue2(nrsubst)*gl%Factor/gl%wm(nrsubst)
!!                                        gl%visco%h2o_read= .true.
!!                                        exit loopaux
!!
!!                                        !elseif ((dummy(1:3) == 'I08') .AND. allocated(gl%visco)) then
!!                                    elseif ((dummy(1:3) == 'CI1') .AND. (gl%visco%pointereta(nrsubst) == 'CI1')) then
!!                                        goover32:               do !until dummy starts with "!"
!!                                            read(fld_unit,*) dummy
!!                                            if (dummy(1:1) == '!') then
!!                                                exit goover32
!!                                            endif
!!                                        end do goover32
!!                                        read(fld_unit,*) gl%visco%low_tempcoll(nrsubst)               !lower temperature limit [K]
!!                                        read(fld_unit,*) gl%visco%upp_tempcoll(nrsubst)               !upper temperature limit [K]
!!                                        read(fld_unit,*) gl%visco%low_prescoll(nrsubst)               !(dummy) upper pressure limit
!!                                        read(fld_unit,*) gl%visco%upp_denscoll(nrsubst)               !(dummy) maximum density
!!                                        read(fld_unit,*) gl%visco%colltermnr(nrsubst)               !number of terms
!!                                        do k = 1, gl%visco%colltermnr(nrsubst), 1            !coeff, power of Tstar
!!                                            read(fld_unit,*) gl%visco%coeffcoll(nrsubst,k),gl%visco%powercoll(nrsubst,k)
!!                                        end do
!!                                        gl%visco%coll_read=.true.
!!                                        exit loopaux
!!
!!                                    elseif ((dummy(1:3) == 'CI2') .AND. (gl%visco%pointereta(nrsubst) == 'CI2')) then
!!                                        goover33:               do !until dummy starts with "!"
!!                                            read(fld_unit,*) dummy
!!                                            if (dummy(1:1) == '!') then
!!                                                exit goover33
!!                                            endif
!!                                        end do goover33
!!                                        read(fld_unit,*) gl%visco%low_tempcoll(nrsubst)               !lower temperature limit [K]
!!                                        read(fld_unit,*) gl%visco%upp_tempcoll(nrsubst)               !upper temperature limit [K]
!!                                        read(fld_unit,*) gl%visco%low_prescoll(nrsubst)               !(dummy) upper pressure limit
!!                                        read(fld_unit,*) gl%visco%upp_denscoll(nrsubst)               !(dummy) maximum density
!!                                        read(fld_unit,*) gl%visco%colltermnr(nrsubst)                 !number of terms
!!                                        do k = 1, gl%visco%colltermnr(nrsubst)                 !coeff, power of Tstar
!!                                            read(fld_unit,*) gl%visco%coeffcoll(nrsubst,k)
!!                                        end do
!!                                        gl%visco%coll_read=.true.
!!                                        exit loopaux
!!                                    else
!!                                        goto 1337
!!                                    end if
!!                                elseif((dummy(1:3) == 'TK1').or.(dummy(1:3) == 'TK3').or.(dummy(1:3) == 'TK6')) then
!!                                    if (.not.allocated(gl%tcx)) then
!!                                        exit loopaux
!!                                    endif
!!
!!                                    if ((dummy(1:3) == 'TK1') .AND. (gl%tcx%tkmodel(nrsubst) == 'TK1')) then
!!                                        goover34:               do !until dummy starts with "!"
!!                                            read(fld_unit,*) dummy
!!                                            if (dummy(1:1) == '!') then
!!                                                exit goover34
!!                                            endif
!!                                        end do goover34
!!                                        read(fld_unit,*) gl%tcx%tmin_tk(nrsubst)               !lower temperature limit [K]
!!                                        read(fld_unit,*) gl%tcx%tmax_tk(nrsubst)               !upper temperature limit [K]
!!                                        read(fld_unit,*) gl%tcx%pmax_tk(nrsubst)               !(dummy) upper pressure limit
!!                                        read(fld_unit,*) gl%tcx%rhomax_tk(nrsubst)               !(dummy) maximum density
!!                                        read(fld_unit,*) gl%tcx%npnum_tk1(nrsubst),gl%tcx%npdenom_tk1(nrsubst),gl%tcx%nexp_tk1(nrsubst),gl%tcx%nspare_tk1(nrsubst)        !number of terms
!!                                        read(fld_unit,*) gl%tcx%tred_tk(nrsubst), gl%tcx%rhored_tk(nrsubst), gl%tcx%tcx_tk(nrsubst)
!!                                        do i=1,gl%tcx%npnum_tk1(nrsubst)
!!                                            read(fld_unit,*) gl%tcx%a_tk1(nrsubst,i), gl%tcx%tsum_tk1(nrsubst,i), gl%tcx%texp_tk1(nrsubst,i), gl%tcx%dsum_tk1(nrsubst,i), gl%tcx%dexp_tk1(nrsubst,i)
!!                                        end do
!!                                        do i=gl%tcx%npnum_tk1(nrsubst)+1,gl%tcx%npnum_tk1(nrsubst)+gl%tcx%npdenom_tk1(nrsubst)
!!                                            read(fld_unit,*) gl%tcx%a_tk1(nrsubst,i), gl%tcx%tsum_tk1(nrsubst,i), gl%tcx%texp_tk1(nrsubst,i), gl%tcx%dsum_tk1(nrsubst,i), gl%tcx%dexp_tk1(nrsubst,i)
!!                                        end do
!!                                        read(fld_unit,*) gl%tcx%tred_tk(nrsubst), gl%tcx%rhored_tk(nrsubst), gl%tcx%tcx_tk(nrsubst)
!!                                        gl%tcx%rhored_tk(nrsubst)=gl%tcx%rhored_tk(nrsubst)*gl%factor
!!                                        do i=gl%tcx%npnum_tk1(nrsubst)+gl%tcx%npdenom_tk1(nrsubst)+1,gl%tcx%npnum_tk1(nrsubst)+gl%tcx%npdenom_tk1(nrsubst)+gl%tcx%nexp_tk1(nrsubst)
!!                                            read(fld_unit,*) gl%tcx%a_tk1(nrsubst,i), gl%tcx%tsum_tk1(nrsubst,i), gl%tcx%texp_tk1(nrsubst,i), gl%tcx%dsum_tk1(nrsubst,i), gl%tcx%dexp_tk1(nrsubst,i)
!!                                        end do
!!                                        gl%tcx%tcx_crit_read=.true.
!!                                        exit loopaux
!!
!!                                    elseif ((dummy(1:3) == 'TK3') .AND. (gl%tcx%tkmodel(nrsubst) == 'TK3')) then
!!                                        goover35:               do !until dummy starts with "!"
!!                                            read(fld_unit,*) dummy
!!                                            if (dummy(1:1) == '!') then
!!                                                exit goover35
!!                                            endif
!!                                        end do goover35
!!                                        read(fld_unit,*) gl%tcx%tmin_tk(nrsubst)               !lower temperature limit [K]
!!                                        read(fld_unit,*) gl%tcx%tmax_tk(nrsubst)               !upper temperature limit [K]
!!                                        read(fld_unit,*) gl%tcx%pmax_tk(nrsubst)               !(dummy) upper pressure limit
!!                                        read(fld_unit,*) gl%tcx%rhomax_tk(nrsubst)               !(dummy) maximum density
!!                                        read(fld_unit,*) gl%tcx%term_tk(nrsubst)                 !number of terms
!!                                        read(fld_unit,*) gl%tcx%tred_tk(nrsubst), gl%tcx%rhored_tk(nrsubst), gl%tcx%tcx_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%gnu_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%gamma_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%R0_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%z_visonly_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%c_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%xi0_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%gam0_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%qd_inverse_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%tref_tk(nrsubst)
!!                                        gl%tcx%tcx_crit_read=.true.
!!                                        exit loopaux
!!
!!                                    elseif ((dummy(1:3) == 'TK6') .AND. (gl%tcx%tkmodel(nrsubst) == 'TK6')) then
!!                                        goover36:               do !until dummy starts with "!"
!!                                            read(fld_unit,*) dummy
!!                                            if (dummy(1:1) == '!') then
!!                                                exit goover36
!!                                            endif
!!                                        end do goover36
!!                                        read(fld_unit,*) gl%tcx%tmin_tk(nrsubst)               !lower temperature limit [K]
!!                                        read(fld_unit,*) gl%tcx%tmax_tk(nrsubst)               !upper temperature limit [K]
!!                                        read(fld_unit,*) gl%tcx%pmax_tk(nrsubst)               !(dummy) upper pressure limit
!!                                        read(fld_unit,*) gl%tcx%rhomax_tk(nrsubst)               !(dummy) maximum density
!!                                        read(fld_unit,*) gl%tcx%term_tk(nrsubst)                 !number of terms
!!                                        read(fld_unit,*) gl%tcx%tred_tk(nrsubst), gl%tcx%rhored_tk(nrsubst), gl%tcx%tcx_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%gnu_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%gamma_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%R0_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%z_visonly_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%c_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%xi0_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%gam0_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%qd_inverse_tk(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%tref_tk(nrsubst)
!!                                        gl%tcx%tcx_crit_read=.true.
!!                                        exit loopaux
!!                                    else
!!                                        goto 1337
!!                                    endif
!!                                else
!!                                    goto 1337
!!                                end if
!!                            enddo loopaux
!!                            Goto 1337
!!
!!                        elseif (dummygen(1:4) == '#ASO') then
!!                            read(fld_unit,*) gl%assoc(nrsubst)                 ! multiplicative coefficient
!!                            read(fld_unit,*) gl%assovolred(nrsubst)            ! reduced normalization volume
!!                            read(fld_unit,*) gl%assovolint(nrsubst)            ! dimensionless volume of interaction
!!                            read(fld_unit,*) gl%assoenergy(nrsubst)            ! dimensionless association energy
!!                            gl%assoexist(nrsubst) = .true.
!!                            Goto 1337
!!
!!                        elseif ((dummygen(1:3) == '#DE') .and. (gl%transport == 3) .and. (allocated(gl%de)))   then
!!                            gl%de%de_read = .true.
!!                            loopDE:    do
!!                                read(fld_unit,*) gl%de%demodel(nrsubst)
!!                                gooverDE:               do !until dummy starts with "!"
!!                                    read(fld_unit,*) dummy
!!                                    call uppertolower_char(dummy, len(dummy))
!!                                    if (dummy(1:1) == '!') then
!!                                        exit gooverDE
!!                                    elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                                        read(fld_unit,'(A)')dummy
!!                                        do while (trim(dummy(1:5)) /= '!````')
!!                                            gl%litref%lit_ref_de(nrsubst) = trim(gl%litref%lit_ref_de(nrsubst)) // dummy(2:)
!!                                            read(fld_unit,'(A)')dummy
!!                                        end do
!!                                        backspace(fld_unit)
!!                                    endif
!!                                end do gooverDE
!!                                !reading Parameters for DE3 and DE4 together, spereated from DE2
!!                                if ((gl%de%demodel(nrsubst) == 'DE3') .OR. (gl%de%demodel(nrsubst) == 'DE4'))  then
!!                                    read(fld_unit,*) gl%de%tmin_de(nrsubst)
!!                                    read(fld_unit,*) gl%de%tmax_de(nrsubst)
!!                                    read(fld_unit,*) gl%de%pmax_de(nrsubst)
!!                                    read(fld_unit,*) gl%de%rhomax_de(nrsubst)
!!                                    read(fld_unit,*) gl%de%ref_temp_de(nrsubst)
!!                                    read(fld_unit,*) gl%de%term_num1_de(nrsubst), gl%de%term_num2_de(nrsubst), gl%de%term_num3_de(nrsubst), &
!!                                        gl%de%term_num4_de(nrsubst), gl%de%term_num5_de(nrsubst), gl%de%term_num6_de(nrsubst)
!!
!!                                    !reading coefficients
!!                                    if ((gl%de%term_num1_de(nrsubst) /= 0).or.(gl%de%term_num2_de(nrsubst) /= 0).or.(gl%de%term_num3_de(nrsubst) /= 0)) then
!!                                        do ide=1,(gl%de%term_num1_de(nrsubst)+gl%de%term_num2_de(nrsubst)+gl%de%term_num3_de(nrsubst))
!!                                            read(fld_unit,*) gl%de%coeffde(nrsubst,ide),gl%de%texpde(nrsubst,ide),gl%de%dexpde(nrsubst,ide)
!!                                        end do
!!                                    end if
!!                                    !reading Parameters for DE2
!!                                elseif (gl%de%demodel(nrsubst) == 'DE2') then
!!                                    read(fld_unit,*) gl%de%tmin_de(nrsubst)
!!                                    read(fld_unit,*) gl%de%tmax_de(nrsubst)
!!                                    read(fld_unit,*) gl%de%pmax_de(nrsubst)
!!                                    read(fld_unit,*) gl%de%rhomax_de(nrsubst)
!!                                    read(fld_unit,*) gl%de%ref_temp_de(nrsubst) , gl%de%ref_dens_de(nrsubst), gl%de%ref_press_de(nrsubst)
!!                                    read(fld_unit,*) gl%de%term_num1_de(nrsubst), gl%de%term_num2_de(nrsubst), gl%de%term_num3_de(nrsubst), &
!!                                        gl%de%term_num4_de(nrsubst), gl%de%term_num5_de(nrsubst), gl%de%term_num6_de(nrsubst)
!!
!!                                    if ((gl%de%term_num1_de(nrsubst) /= 0).or.(gl%de%term_num2_de(nrsubst) /= 0).or.(gl%de%term_num3_de(nrsubst) /= 0)) then
!!                                        do ide=1,(gl%de%term_num1_de(nrsubst)+gl%de%term_num2_de(nrsubst)+gl%de%term_num3_de(nrsubst))
!!                                            read(fld_unit,*) gl%de%coeffde(nrsubst,ide),gl%de%texpde(nrsubst,ide),gl%de%dexpde(nrsubst,ide),gl%de%pexpde(nrsubst,ide)
!!                                        end do
!!                                    end if
!!
!!                                endif
!!                                exit loopDE
!!                            enddo loopDE
!!                            Goto 1337
!!
!!
!!                        elseif (dummygen(1:4) == '#STN') then
!!                            gl%stn_read = .true.
!!                            !allocate the stn type to the size of components
!!                            if(.not.allocated(gl%stn)) allocate(surface_tension_t::gl%stn(gl%ncomp))
!!
!!                            loopST:    do
!!                                read(fld_unit,*) gl%stn(nrsubst)%stmodel
!!                                gooverST:               do !until dummy starts with "!"
!!                                    read(fld_unit,*) dummy
!!                                    call uppertolower_char(dummy, len(dummy))
!!                                    if (dummy(1:1) == '!') then
!!                                        exit gooverST
!!                                    elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                                        do while (trim(dummy(1:5)) /= '!````')
!!                                            gl%litref%lit_ref_stn(nrsubst) = trim(gl%litref%lit_ref_stn(nrsubst)) // dummy(2:)
!!                                            read(fld_unit,'(A)')dummy
!!                                        end do
!!                                        backspace(fld_unit)
!!                                    endif
!!                                end do gooverST
!!                                !reading Parameters for ST1
!!
!!                                if (gl%stn(nrsubst)%stmodel == 'ST1') then
!!                                    read(fld_unit,*) gl%stn(nrsubst)%low_temp_st
!!                                    read(fld_unit,*) gl%stn(nrsubst)%upp_temp_st
!!                                    read(fld_unit,*) gl%stn(nrsubst)%upp_press_st
!!                                    read(fld_unit,*) gl%stn(nrsubst)%max_dens_st
!!                                    read(fld_unit,*) gl%stn(nrsubst)%term_num_st
!!                                    read(fld_unit,*) gl%stn(nrsubst)%Temp_crit_st
!!
!!
!!
!!                                    do ist=1,(gl%stn(nrsubst)%term_num_st)
!!                                        read(fld_unit,*) gl%stn(nrsubst)%sigma0_st(ist),gl%stn(nrsubst)%n_exp_st(ist)
!!                                    end do
!!                                endif
!!                                exit loopST
!!                            enddo loopST
!!                            Goto 1337
!!
!!
!!                            ! Lars / Andreas (written by Monika Thol??), March 2014
!!                            !Implemented for ECS models (transport properties)
!!                            !-----------
!!                        elseif ((dummygen(1:4) == '#TCX') .and. (gl%transport == 2) .and. (allocated(gl%tcx)) ) then
!!                            gl%tcx%tcx_read = .true.
!!                            loopTC:          do
!!                                read(fld_unit,*) gl%tcx%tcmodel(nrsubst)
!!                                gooverTC:                do
!!                                    read(fld_unit,*) dummy
!!                                    call uppertolower_char(dummy, len(dummy))
!!                                    if (dummy(1:1) == '!') then
!!                                        exit gooverTC
!!                                    elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                                        read(fld_unit,'(A)')dummy
!!                                        do while (trim(dummy(1:5)) /= '!````')
!!                                            gl%litref%lit_ref_tcx(nrsubst) = trim(gl%litref%lit_ref_tcx(nrsubst)) // dummy(2:)
!!                                            read(fld_unit,'(A)')dummy
!!                                        end do
!!                                        backspace(fld_unit)
!!                                    endif
!!                                end do gooverTC
!!                                if (gl%tcx%tcmodel(nrsubst) == 'TC1') then
!!                                    read(fld_unit,*) gl%tcx%tmin_tc(nrsubst)                                  ! K
!!                                    read(fld_unit,*) gl%tcx%tmax_tc(nrsubst)                                  ! K
!!                                    read(fld_unit,*) gl%tcx%pmax_tc(nrsubst)                                  ! kPa
!!                                    gl%tcx%pmax_tc(nrsubst) = gl%tcx%pmax_tc(nrsubst)*factor_inv              !MPa
!!                                    read(fld_unit,*) gl%tcx%rhomax_tc(nrsubst)                                ! mol/L
!!                                    gl%tcx%rhomax_tc(nrsubst) = gl%tcx%rhomax_tc(nrsubst)*gl%Factor          ! mol/m?
!!                                    read(fld_unit,*) gl%tcx%num_dil_tc(nrsubst), gl%tcx%den_dil_tc(nrsubst)
!!                                    read(fld_unit,*) gl%tcx%tred_dil(nrsubst), gl%tcx%tcx_dil(nrsubst)
!!                                    do i=1, gl%tcx%num_dil_tc(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%num_tc_coeff_dil(nrsubst,i), gl%tcx%num_tc_powerT_dil(nrsubst,i)
!!                                    end do
!!                                    do i=1, gl%tcx%den_dil_tc(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%den_tc_coeff_dil(nrsubst,i), gl%tcx%den_tc_exp_dil(nrsubst,i)
!!                                    enddo
!!                                    read(fld_unit,*) gl%tcx%num_bgrd_tc(nrsubst), gl%tcx%den_bgrd_tc(nrsubst)
!!                                    read(fld_unit,*) gl%tcx%tred_bgrd(nrsubst), gl%tcx%rhored_bgrd(nrsubst), gl%tcx%tcx_bgrd(nrsubst)
!!                                    do i=1, gl%tcx%num_bgrd_tc(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%num_tc_coeff_bgrd(nrsubst,i), gl%tcx%num_tc_powerT_bgrd(nrsubst,i), gl%tcx%num_tc_rho_bgrd(nrsubst,i), gl%tcx%rho_exp_bgrd(nrsubst,i)
!!                                    end do
!!                                    read(fld_unit,*) gl%tcx%tkmodel(nrsubst)
!!                                elseif (gl%tcx%tcmodel(nrsubst) == 'TC0') then
!!                                    read(fld_unit,*) gl%tcx%tmin_tc(nrsubst)                                  ! K
!!                                    read(fld_unit,*) gl%tcx%tmax_tc(nrsubst)                                  ! K
!!                                    read(fld_unit,*) gl%tcx%pmax_tc(nrsubst)                                  ! kPa
!!                                    gl%tcx%pmax_tc(nrsubst) = gl%tcx%pmax_tc(nrsubst)*factor_inv              ! MPa
!!                                    read(fld_unit,*) gl%tcx%rhomax_tc(nrsubst)                                ! mol/L
!!                                    gl%tcx%rhomax_tc(nrsubst) = gl%tcx%rhomax_tc(nrsubst)*gl%Factor
!!                                    read(fld_unit,*) gl%tcx%pointer_hardtc(nrsubst)
!!                                    read(fld_unit,*) gl%tcx%tcxlam0_num(nrsubst), gl%tcx%tcxlam0_den(nrsubst), gl%tcx%dellam(nrsubst), gl%tcx%dellamcr(nrsubst)
!!
!!                                    If(gl%tcx%pointer_hardtc(nrsubst) == 'R23') then
!!                                        read(fld_unit,*) gl%tcx%tred_tc(nrsubst), gl%tcx%rhored_tc(nrsubst), gl%tcx%etared_tc(nrsubst)
!!                                        read(fld_unit,*) gl%tcx%rholl_r23
!!                                        gl%tcx%rholl_r23=gl%tcx%rholl_r23*gl%factor
!!                                        read(fld_unit,*) gl%tcx%b1_r23
!!                                        read(fld_unit,*) gl%tcx%b2_r23
!!                                        read(fld_unit,*) gl%tcx%c1l_r23
!!                                        read(fld_unit,*) gl%tcx%c2l_r23
!!                                        read(fld_unit,*) gl%tcx%delgl_r23
!!                                        read(fld_unit,*) gl%tcx%lmax_r23
!!                                        exit loopTC
!!                                    end if
!!
!!                                    if (gl%tcx%pointer_hardtc(nrsubst) /= 'ETY') then
!!                                        read(fld_unit,*) gl%tcx%tred_tc(nrsubst), gl%tcx%rhored_tc(nrsubst), gl%tcx%etared_tc(nrsubst)      ! use tred_bgrd/rhored_tc instead of tred_tc/rhored_tc?!
!!                                        do i=1, gl%tcx%tcxlam0_num(nrsubst)
!!                                            read(fld_unit,*) gl%tcx%lam0_coeff(nrsubst,i), gl%tcx%lam0_exp(nrsubst,i)
!!                                        enddo
!!                                        ! here lamo0_den_coeff...
!!                                        do i=1, gl%tcx%dellam(nrsubst)
!!                                            read(fld_unit,*) gl%tcx%dellam_coeff(nrsubst,i), gl%tcx%dellam_exp(nrsubst,i), gl%tcx%dellam_exp2(nrsubst,i)
!!                                        enddo
!!                                        do i=1, gl%tcx%dellamcr(nrsubst)
!!                                            read(fld_unit,*) gl%tcx%dellamcr_coeff(nrsubst,i), gl%tcx%dellamcr_exp(nrsubst,i)
!!                                        enddo
!!                                    endif
!!                                    read(fld_unit,*) gl%tcx%tkmodel(nrsubst)
!!
!!                                elseif (gl%tcx%tcmodel(nrsubst) == 'TC3') then                 ! xenon
!!                                    read(fld_unit,*) gl%tcx%tmin_tc(nrsubst)                          ! K
!!                                    read(fld_unit,*) gl%tcx%tmax_tc(nrsubst)                          ! K
!!                                    read(fld_unit,*) gl%tcx%pmax_tc(nrsubst)                          ! kPa
!!                                    gl%tcx%pmax_tc(nrsubst) = gl%tcx%pmax_tc(nrsubst)*factor_inv      ! MPa
!!                                    read(fld_unit,*) gl%tcx%rhomax_tc(nrsubst)                        ! mol/L
!!                                    gl%tcx%rhomax_tc(nrsubst) = gl%tcx%rhomax_tc(nrsubst)*gl%Factor  ! mol/m?
!!                                    read(fld_unit,*) gl%tcx%lj_sigma(nrsubst)                         ! nm
!!                                    read(fld_unit,*) gl%tcx%lj_epskap(nrsubst)                        ! K
!!                                    read(fld_unit,*) gl%tcx%const_eq20(nrsubst)
!!                                    read(fld_unit,*) gl%tcx%Texp_eq20(nrsubst)
!!                                    do i=1, 9   ! manually
!!                                        read(fld_unit,*) gl%tcx%eta0(nrsubst,i)
!!                                    enddo
!!                                    do i=1, 4   ! manually
!!                                        read(fld_unit,*) gl%tcx%Fvi(nrsubst,i)
!!                                    enddo
!!                                    do i=1, 8   ! manually
!!                                        read(fld_unit,*) gl%tcx%Evi(nrsubst,i)
!!                                    enddo
!!                                    read(fld_unit,*) gl%tcx%F_tc3(nrsubst)
!!                                    read(fld_unit,*) gl%tcx%rm(nrsubst)
!!                                    read(fld_unit,*) gl%tcx%tkmodel(nrsubst)
!!                                endif
!!                                exit loopTC
!!                            end do loopTC
!!                            Goto 1337
!!                            !-----------
!!                            !-----------
!!
!!
!!
!!                        elseif((dummygen(1:4) == '#ETA').and.((gl%transport == 1) .or. (gl%transport == 2)).and.(allocated(gl%visco))) then  !read parameters for viscosity equations
!!                            gl%visco%eta_read = .true.
!!                            loopETA:    do
!!                                read(fld_unit,*) gl%visco%etamodel(nrsubst)
!!                                gooverETA:               do!until dummy starts with "!"
!!                                    read(fld_unit,*) dummy
!!                                    call uppertolower_char(dummy, len(dummy))
!!                                    if (dummy(1:1) == '!') then
!!                                        exit gooverETA
!!                                    elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                                        read(fld_unit,'(A)')dummy
!!                                        do while (trim(dummy(1:5)) /= '!````')
!!                                            gl%litref%lit_ref_eta(nrsubst) = trim(gl%litref%lit_ref_eta(nrsubst)) // dummy(2:)
!!                                            read(fld_unit,'(A)')dummy
!!                                        end do
!!                                        backspace(fld_unit)
!!                                    endif
!!                                end do gooverETA
!!
!!                                !reading parameters for VS0 model specification
!!                                if (gl%visco%etamodel(nrsubst) == 'VS0') then
!!                                    read(fld_unit,*) gl%visco%tmineta(nrsubst)    ![K]
!!                                    read(fld_unit,*) gl%visco%tmaxeta(nrsubst)    ![K]
!!                                    read(fld_unit,*) gl%visco%pmaxeta(nrsubst)   ![kPa]
!!                                    gl%visco%pmaxeta(nrsubst)=gl%visco%pmaxeta(nrsubst)*factor_inv
!!                                    read(fld_unit,*) gl%visco%rhomaxeta(nrsubst)    ![mol/L]
!!                                    gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%factor
!!                                    read(fld_unit,*) gl%visco%pointer_lambda(nrsubst)     !Pointer to hardcoded model
!!                                    if (trim(gl%visco%pointer_lambda(nrsubst))  == 'LJF' .or. trim(gl%visco%pointer_lambda(nrsubst))  == 'LJ1' ) then
!!                                        return
!!                                    end if
!!                                    if((trim(gl%visco%pointer_lambda(nrsubst)) == 'HE') &
!!                                        & .or.(trim(gl%visco%pointer_lambda(nrsubst)) == 'ETY').or.(trim(gl%visco%pointer_lambda(nrsubst)) == 'NEO')) exit
!!                                        
!!                                        read(fld_unit,*) gl%visco%nterm1_eta(nrsubst), gl%visco%nterm2_eta(nrsubst), gl%visco%nterm3_eta(nrsubst), gl%visco%nterm4_eta(nrsubst), &
!!                                            & gl%visco%nterm5_eta(nrsubst), gl%visco%nterm6_eta(nrsubst), gl%visco%nterm7_eta(nrsubst), gl%visco%nterm8_eta(nrsubst)
!!
!!                                    if (trim(gl%visco%pointer_lambda(nrsubst)) == 'R23') then
!!                                        read(fld_unit,*) gl%visco%pointereta(nrsubst)
!!                                        read(fld_unit,*) gl%visco%lejo_sigma(nrsubst)
!!                                        read(fld_unit,*) gl%visco%lejo_epka(nrsubst)
!!
!!                                        read(fld_unit,*) gl%visco%red_temp_eta(nrsubst), gl%visco%red_vis_eta(nrsubst)
!!                                        read(fld_unit,*) gl%visco%Chapman_Enskog1(nrsubst), gl%visco%Chapman_Enskog2(nrsubst)
!!
!!                                        read(fld_unit,*) gl%visco%red_temp_eta(nrsubst) , gl%visco%red_dens_eta(nrsubst), gl%visco%red_vis_eta(nrsubst)    !reducing parameters
!!                                        gl%visco%red_dens_eta(nrsubst)=gl%visco%red_dens_eta(nrsubst)*gl%Factor
!!                                        read(fld_unit,*) gl%visco%rhol_R23
!!                                        gl%visco%rhol_R23=gl%visco%rhol_R23*gl%factor
!!                                        read(fld_unit,*) gl%visco%c1_R23
!!                                        read(fld_unit,*) gl%visco%c2_R23
!!                                        read(fld_unit,*) gl%visco%delg_R23
!!                                        read(fld_unit,*) gl%visco%etamax_R23
!!                                        read(fld_unit,*) gl%visco%pointercrit_eta(nrsubst)
!!                                        exit loopETA
!!                                    end if
!!
!!                                    read(fld_unit,*) gl%visco%red_temp_eta(nrsubst) , gl%visco%red_dens_eta(nrsubst), gl%visco%red_vis_eta(nrsubst)    !reducing parameters
!!                                    gl%visco%red_dens_eta(nrsubst)=gl%visco%red_dens_eta(nrsubst)*gl%Factor
!!
!!                                    do i = 1,gl%visco%nterm1_eta(nrsubst)   ! Loop for reduced viscosity at high density according to the Enskog theory
!!                                        read(fld_unit,*) gl%visco%coeff_ei(nrsubst,i),gl%visco%exp_ei(nrsubst,i)
!!                                    end do
!!
!!                                    do i = 1,gl%visco%nterm2_eta(nrsubst)   ! Loop for dilute-gas viscosity Mue0
!!                                        read(fld_unit,*) gl%visco%coeff_dil(nrsubst,i),gl%visco%exp_dil(nrsubst,i)
!!                                    end do
!!
!!                                    do i = 1,(gl%visco%nterm3_eta(nrsubst)+gl%visco%nterm4_eta(nrsubst)) ! Loop for water and ethylene
!!                                        read(fld_unit,*) gl%visco%coeff_hieta(nrsubst,i),dummy1,dummy2, gl%visco%exp1eta(nrsubst,i), gl%visco%exp2eta(nrsubst,i)
!!                                        gl%visco%iexpeta(nrsubst,i) = int(dummy1)
!!                                        gl%visco%jexpeta(nrsubst,i) = int(dummy2)
!!                                        if (trim(gl%components(nrsubst)) == 'water') then
!!                                            gl%visco%hieta(nrsubst, gl%visco%iexpeta(nrsubst,i)+1, gl%visco%jexpeta(nrsubst,i)+1)=gl%visco%coeff_hieta(nrsubst,i)
!!                                        end if
!!                                    end do
!!
!!                                    do i = 1,gl%visco%nterm5_eta(nrsubst)   ! Loop for viscosity in the zero density limit
!!                                        read(fld_unit,*) gl%visco%coeff_ai(nrsubst,i),gl%visco%exp_ai(nrsubst,i)
!!                                    end do
!!
!!                                    do i = 1,gl%visco%nterm6_eta(nrsubst)   ! Loop for second viscosity virial coefficient
!!                                        read(fld_unit,*) gl%visco%coeff_bi(nrsubst,i),gl%visco%exp_bi(nrsubst,i)
!!                                    end do
!!
!!                                    do i = 1,gl%visco%nterm7_eta(nrsubst)   ! Loop for third viscosity virial coefficient
!!                                        read(fld_unit,*) gl%visco%coeff_ci(nrsubst,i),gl%visco%exp_ci(nrsubst,i)
!!                                    end do
!!
!!                                    do i = 1,gl%visco%nterm8_eta(nrsubst)   ! Loop for reduced viscosity at high density according to the Enskog theory
!!                                        read(fld_unit,*) gl%visco%coeff_di(nrsubst,i),gl%visco%exp_di(nrsubst,i)
!!                                    end do
!!
!!                                    read(fld_unit,*) gl%visco%pointercrit_eta(nrsubst)
!!
!!                                    !reading parameters for VS1 model specification
!!                                elseif (gl%visco%etamodel(nrsubst) == 'VS1') then
!!                                    read(fld_unit,*) gl%visco%tmineta(nrsubst)    ![K]
!!                                    read(fld_unit,*) gl%visco%tmaxeta(nrsubst)    ![K]
!!                                    read(fld_unit,*) gl%visco%pmaxeta(nrsubst)   ![kPa]
!!                                    gl%visco%pmaxeta(nrsubst) = gl%visco%pmaxeta(nrsubst)*Factor_inv
!!                                    read(fld_unit,*) gl%visco%rhomaxeta(nrsubst)    ![mol/l]
!!                                    gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%Factor
!!                                    read(fld_unit,*) gl%visco%dilgas_term(nrsubst)     !number of terms associated with dilute-gas function
!!                                    read(fld_unit,*) gl%visco%pointereta(nrsubst)     !Points out which collision integral model is used
!!                                    read(fld_unit,*) gl%visco%slj(nrsubst)      !Lennard-Jones coefficient sigma [nm]
!!                                    read(fld_unit,*) gl%visco%eklj(nrsubst)       !Lennard-Jones coefficient epsilon/kappa [K]
!!                                    read(fld_unit,*) gl%visco%tredetadg(nrsubst), gl%visco%visredetadg(nrsubst)
!!                                    read(fld_unit,*) gl%visco%ce(nrsubst), gl%visco%cep(nrsubst)   !ce (=0.0266958*sqrt(M)), cep (=Power of T)
!!                                    if((gl%substcasnr(nrsubst) == '64-17-5').or.(gl%substcasnr(nrsubst) == '306-83-2')) then !ethanol and R123
!!                                        read(fld_unit,*) gl%visco%a0_vs1(nrsubst)
!!                                        read(fld_unit,*) gl%visco%a1_vs1(nrsubst)
!!                                        read(fld_unit,*) gl%visco%a2_vs1(nrsubst)
!!                                    end if
!!                                    read(fld_unit,*) gl%visco%denstermnr(nrsubst)
!!                                    if (gl%visco%denstermnr(nrsubst) /= 0.d0) THEN
!!                                        read(fld_unit,*) gl%visco%treddens(nrsubst), gl%visco%etaB2(nrsubst) !reducing parameters for T (= eps/k), etaB2 (= 0.6022137*sigma**3)
!!                                        DO jeta=1,gl%visco%denstermnr(nrsubst),1
!!                                            read(fld_unit,*) gl%visco%coeffbetastar(nrsubst,jeta), gl%visco%powerbetastar(nrsubst,jeta)
!!                                        END DO
!!                                    end if
!!                                    read(fld_unit,*) gl%visco%term_num1_eta(nrsubst), gl%visco%term_num2_eta(nrsubst), gl%visco%term_num3_eta(nrsubst), &
!!                                        & gl%visco%term_num4_eta(nrsubst), gl%visco%term_num5_eta(nrsubst), gl%visco%term_num6_eta(nrsubst)
!!                                    read(fld_unit,*) gl%visco%tredeta(nrsubst) , gl%visco%rhoredeta(nrsubst), gl%visco%visredeta(nrsubst)    !reducing parameters
!!                                    gl%visco%rhoredeta(nrsubst)=gl%visco%rhoredeta(nrsubst)*gl%Factor
!!                                    if (gl%visco%term_num1_eta(nrsubst) /= 0) then
!!                                        do ieta = 1, abs(gl%visco%term_num1_eta(nrsubst)), 1
!!                                            read(fld_unit,*) gl%visco%a_vs1(nrsubst,ieta), gl%visco%powertau(nrsubst,ieta)
!!                                        end do
!!                                    end if
!!                                    if (gl%visco%term_num2_eta(nrsubst) /= 0) then
!!                                        do jeta = 1,gl%visco%term_num2_eta(nrsubst),1
!!                                            read(fld_unit,*) gl%visco%b_vs1(nrsubst,jeta), gl%visco%ptausp(nrsubst,jeta),gl%visco%pdelsp(nrsubst,jeta),gl%visco%pdel0sp(nrsubst,jeta),gl%visco%expdel(nrsubst,jeta)
!!                                        enddo
!!                                    end if
!!                                    if (gl%visco%term_num3_eta(nrsubst) /= 0) then
!!                                        do keta = gl%visco%term_num2_eta(nrsubst)+1,gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst),1
!!                                            read(fld_unit,*) gl%visco%b_vs1(nrsubst,keta), gl%visco%ptausp(nrsubst,keta),gl%visco%pdelsp(nrsubst,keta),gl%visco%pdel0sp(nrsubst,keta),gl%visco%expdel(nrsubst,keta)
!!                                        enddo
!!                                    end if
!!                                    if (gl%visco%term_num4_eta(nrsubst) /= 0) then
!!                                        do leta = gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+1,&
!!                                            & gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst),1
!!                                            read(fld_unit,*) gl%visco%b_vs1(nrsubst,leta), gl%visco%ptausp(nrsubst,leta),gl%visco%pdelsp(nrsubst,leta),gl%visco%pdel0sp(nrsubst,leta),gl%visco%expdel(nrsubst,leta)
!!                                        enddo
!!                                    end if
!!                                    if (gl%visco%term_num5_eta(nrsubst) /= 0) then
!!                                        do meta = gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst)+1, &
!!                                            & gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst)+gl%visco%term_num5_eta(nrsubst),1
!!                                            read(fld_unit,*) gl%visco%b_vs1(nrsubst,meta), gl%visco%ptausp(nrsubst,meta),gl%visco%pdelsp(nrsubst,meta),gl%visco%pdel0sp(nrsubst,meta),gl%visco%expdel(nrsubst,meta)
!!                                        enddo
!!                                    end if
!!                                    if (gl%visco%term_num6_eta(nrsubst) /= 0) then
!!                                        do neta = gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst)+gl%visco%term_num5_eta(nrsubst)+1,&
!!                                            & gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst)+gl%visco%term_num5_eta(nrsubst)+gl%visco%term_num6_eta(nrsubst),1
!!                                            read(fld_unit,*) gl%visco%b_vs1(nrsubst,neta), gl%visco%ptausp(nrsubst,neta),gl%visco%pdelsp(nrsubst,neta),gl%visco%pdel0sp(nrsubst,neta),gl%visco%expdel(nrsubst,neta)
!!                                        enddo
!!                                    end if
!!
!!                                    !reading parameters for VS2 model specification
!!                                    !B. A. Younglove, J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982
!!                                elseif (gl%visco%etamodel(nrsubst) == 'VS2') then
!!                                    read(fld_unit,*) gl%visco%tmineta(nrsubst)          ![K]
!!                                    read(fld_unit,*) gl%visco%tmaxeta(nrsubst)          ![K]
!!                                    read(fld_unit,*) gl%visco%pmaxeta(nrsubst)          ![kPa]
!!                                    read(fld_unit,*) gl%visco%rhomaxeta(nrsubst)        ![mol/l]
!!                                    gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%Factor
!!                                    read(fld_unit,*) gl%visco%pointereta(nrsubst)       !Points out which collision integral model is used
!!                                    read(fld_unit,*) gl%visco%slj(nrsubst)              !Lennard-Jones coefficient sigma [nm]
!!                                    read(fld_unit,*) gl%visco%eklj(nrsubst)             !Lennard-Jones coefficient epsilon/kappa [K]
!!                                    read(fld_unit,*) gl%visco%const19_VS2(nrsubst)      !const in Eq 19
!!                                    read(fld_unit,*) gl%visco%exp19_VS2(nrsubst)        !exponent in Eq 19 for T
!!                                    read(fld_unit,*) gl%visco%Fv_VS2(nrsubst,1)         !coefficients for initial density dependence of viscosity for eq. (21)
!!                                    do i=2,4
!!                                        read(fld_unit,*) gl%visco%Fv_VS2(nrsubst,i)
!!                                    enddo
!!                                    read(fld_unit,*) gl%visco%Ev_VS2(nrsubst,1)         !coefficients for residual viscosity for eqs. (22) - (25)
!!                                    do i=2,8
!!                                        read(fld_unit,*) gl%visco%Ev_VS2(nrsubst,i)
!!                                    enddo
!!                                    read(fld_unit,*) gl%visco%pointercrit_eta(nrsubst)  !pointer to critical enhancement auxiliary function (none used)
!!                                    !! unter #AUX CI2 k霵nten die Koeffizienten f�r Omega eingelesen werden
!!
!!                                    !reading parameters for VS4 model specification
!!                                elseif (gl%visco%etamodel(nrsubst) == 'VS4') then
!!                                    read(fld_unit,*) gl%visco%tmineta(nrsubst)    ![K]
!!                                    read(fld_unit,*) gl%visco%tmaxeta(nrsubst)    ![K]
!!                                    read(fld_unit,*) gl%visco%pmaxeta(nrsubst)   ![kPa]
!!                                    gl%visco%pmaxeta(nrsubst) = gl%visco%pmaxeta(nrsubst)*Factor_inv
!!                                    read(fld_unit,*) gl%visco%rhomaxeta(nrsubst)    ![mol/l]
!!                                    gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%Factor
!!                                    read(fld_unit,*) gl%visco%dilgas_term(nrsubst)     !number of terms associated with dilute-gas function
!!                                    read(fld_unit,*) gl%visco%pointereta(nrsubst)     !Points out which collision integral model is used
!!                                    read(fld_unit,*) gl%visco%slj(nrsubst)      !Lennard-Jones coefficient sigma [nm]
!!                                    read(fld_unit,*) gl%visco%eklj(nrsubst)       !Lennard-Jones coefficient epsilon/kappa [K]
!!                                    read(fld_unit,*) gl%visco%tredetadg(nrsubst), gl%visco%visredetadg(nrsubst)
!!                                    read(fld_unit,*) gl%visco%ce(nrsubst), gl%visco%cep(nrsubst)   !ce (=0.0266958*sqrt(M)), cep (=Power of T)
!!                                    do i = 1, gl%visco%dilgas_term(nrsubst)-1
!!                                        read(fld_unit,*) gl%visco%d0_vs4(nrsubst,i), gl%visco%d0exp(nrsubst,i)
!!                                    end do
!!                                    read(fld_unit,*) gl%visco%denstermnr(nrsubst)
!!                                    if (gl%visco%denstermnr(nrsubst) /= 0.d0) THEN !?????
!!                                        read(fld_unit,*) gl%visco%treddens(nrsubst), gl%visco%etaB2(nrsubst) !reducing parameters for T (= eps/k), etaB2 (= 0.6022137*sigma**3)
!!                                        DO jeta=1,gl%visco%denstermnr(nrsubst),1
!!                                            read(fld_unit,*) gl%visco%coeffbetastar(nrsubst,jeta), gl%visco%powerbetastar(nrsubst,jeta)
!!                                        END DO
!!                                    end if
!!                                    read(fld_unit,'(A255)')dummy
!!                                    if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%a_0(nrsubst), gl%visco%a_1(nrsubst), gl%visco%a_2(nrsubst)
!!                                    read(fld_unit,'(A255)')dummy
!!                                    if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%b0_vs4(nrsubst), gl%visco%b1_vs4(nrsubst), gl%visco%b2_vs4(nrsubst)
!!                                    read(fld_unit,'(A255)')dummy
!!                                    if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%c0(nrsubst), gl%visco%c_1(nrsubst), gl%visco%c_2(nrsubst)
!!                                    read(fld_unit,'(A255)')dummy
!!                                    if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%aa0(nrsubst), gl%visco%aa1(nrsubst), gl%visco%aa2(nrsubst)
!!                                    read(fld_unit,'(A255)')dummy
!!                                    if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%bb0(nrsubst), gl%visco%bb1(nrsubst), gl%visco%bb2(nrsubst)
!!                                    read(fld_unit,'(A255)')dummy
!!                                    if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%cc0(nrsubst), gl%visco%cc1(nrsubst), gl%visco%cc2(nrsubst)
!!                                    read(fld_unit,'(A255)')dummy
!!                                    if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%dd_0(nrsubst), gl%visco%dd_1(nrsubst), gl%visco%dd_2(nrsubst)
!!                                    read(fld_unit,'(A255)')dummy
!!                                    if (dummy(1:3) /= 'NUL') read(dummy,*) gl%visco%e_0(nrsubst), gl%visco%e_1(nrsubst), gl%visco%e_2(nrsubst)
!!
!!                                    !reading parameters for VS5 model specification
!!                                    !T-H. Chung, M. Ajlan, L.L. Lee and K.E. Starling, Ind. Eng. Chem. Res. 1998, 27, 671-679.
!!                                elseif (gl%visco%etamodel(nrsubst) == 'VS5') then
!!                                    read(fld_unit,*) gl%visco%tmineta(nrsubst)                     ![K]
!!                                    read(fld_unit,*) gl%visco%tmaxeta(nrsubst)                     ![K]
!!                                    read(fld_unit,*) gl%visco%pmaxeta(nrsubst)                     ![kPa]
!!                                    read(fld_unit,*) gl%visco%rhomaxeta(nrsubst)                   ![mol/l]
!!                                    gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%Factor
!!                                    read(fld_unit,*) gl%visco%dilgas_term(nrsubst)                 !number of terms associated with dilute-gas function
!!                                    read(fld_unit,*) gl%visco%pointereta(nrsubst)                  !Points out which collision integral model is used
!!                                    read(fld_unit,*) gl%visco%slj(nrsubst)                         !Lennard-Jones coefficient sigma [nm]
!!                                    gl%visco%slj_VS5(nrsubst) = gl%visco%slj(nrsubst)*gl%factor_VS5slj  !unit conversion to Angstroem
!!                                    read(fld_unit,*) gl%visco%eklj(nrsubst)                        !Lennard-Jones coefficient epsilon/kappa [K] =Tc/1.2593
!!                                    gl%visco%tc_VS5(nrsubst) = gl%visco%eklj(nrsubst)*1.2593d0       !Chung uses a different tc value
!!                                    read(fld_unit,*) gl%visco%tredetadg(nrsubst), gl%visco%visredetadg(nrsubst)
!!                                    read(fld_unit,*) gl%visco%ce(nrsubst), gl%visco%cep(nrsubst)            !ce (=0.021357*SQRT(MW)), cep (=Power of T)
!!                                    read(fld_unit,*) gl%visco%denstermnr(nrsubst)
!!                                    if (gl%visco%denstermnr(nrsubst) /= 0.d0) THEN
!!                                        !not yet used for VS5
!!                                    end if
!!                                    read(fld_unit,*) gl%visco%accen_VS5(nrsubst), &
!!                                        & gl%visco%dipolered_VS5(nrsubst), gl%visco%kappa_VS5(nrsubst)       !w, mur, kappa for Chung, fit
!!                                    read(fld_unit,*) gl%visco%addchung_VS5(nrsubst)                !additional parameters for Chung used if residual parameters are not the same as dilute gas
!!                                    read(fld_unit,*) gl%visco%pointercrit_eta(nrsubst)             !pointer to critical enhancement auxiliary function (none used)
!!
!!                                    !reading parameters for VS7 model specification
!!                                    !cf. Vogel, Span, Herrmann, J. Phys. Chem. Ref. Data 44, 043101 (2015)
!!                                elseif (gl%visco%etamodel(nrsubst) == 'VS7') then
!!                                    read(fld_unit,*) gl%visco%tmineta(nrsubst)    ![K]
!!                                    read(fld_unit,*) gl%visco%tmaxeta(nrsubst)    ![K]
!!                                    read(fld_unit,*) gl%visco%pmaxeta(nrsubst)   ![kPa]
!!                                    gl%visco%pmaxeta(nrsubst) = gl%visco%pmaxeta(nrsubst)*Factor_inv
!!                                    read(fld_unit,*) gl%visco%rhomaxeta(nrsubst)    ![mol/l]
!!                                    gl%visco%rhomaxeta(nrsubst)=gl%visco%rhomaxeta(nrsubst)*gl%Factor
!!                                    read(fld_unit,*) gl%visco%omega_model(nrsubst)
!!                                    read(fld_unit,*) gl%visco%ndens_initial(nrsubst)
!!                                    read(fld_unit,*) gl%visco%tredetadg(nrsubst), gl%visco%rhoredeta(nrsubst), gl%visco%visredetadg(nrsubst) !reducing parameters for T, rho, eta
!!
!!                                    !read diluted gas terms
!!                                    ieta = 1
!!                                    do
!!                                        !err_read = 0
!!                                        read(fld_unit,*,iostat = err_read) gl%visco%coeff_hi(nrsubst,ieta)
!!                                        if (err_read .ne. 0) exit
!!                                        ieta = ieta + 1
!!                                    end do
!!                                    gl%visco%dilgas_term(nrsubst) = ieta - 1 !number of terms associated with dilute-gas function
!!
!!                                    !read(fld_unit,*) gl%visco%tinit_red(nrsubst), gl%visco%rhoredeta(nrsubst), gl%visco%etainit_red(nrsubst) !reducing parameters for T, rho, eta
!!                                    !gl%visco%rhoredeta(nrsubst)=gl%visco%rhoredeta(nrsubst)*gl%Factor
!!
!!                                    ieta = 1
!!                                    !read initial density terms
!!                                    if (gl%visco%ndens_initial(nrsubst) .ne. 0) then
!!
!!                                        read(fld_unit,*) gl%visco%tinit_red(nrsubst), gl%visco%rhoredeta(nrsubst), gl%visco%etainit_red(nrsubst) !reducing parameters for T, rho, eta
!!                                        gl%visco%rhoredeta(nrsubst)=gl%visco%rhoredeta(nrsubst)*gl%Factor
!!
!!                                        do ieta = 1, gl%visco%ndens_initial(nrsubst)
!!                                            read(fld_unit,*,iostat = err_read) gl%visco%cieta(nrsubst,ieta), gl%visco%tieta(nrsubst,ieta)
!!                                        end do
!!                                        read(fld_unit,*) dummy
!!                                    end if
!!
!!                                    read(fld_unit,*) gl%visco%tredeta(nrsubst), gl%visco%rhoredeta(nrsubst), gl%visco%visredeta(nrsubst) !reducing parameters for T, rho, eta
!!                                    gl%visco%rhoredeta(nrsubst) = gl%visco%rhoredeta(nrsubst) * gl%factor
!!
!!                                    !read polynomial terms
!!                                    ieta = 1
!!                                    do
!!                                        read(fld_unit,*,iostat = err_read) gl%visco%coeff(nrsubst,ieta), gl%visco%texpo(nrsubst,ieta), gl%visco%dexpo(nrsubst,ieta), gl%visco%pi_VS7(nrsubst,ieta), gl%visco%li_VS7(nrsubst,ieta)
!!                                        if (err_read .ne. 0) exit
!!                                        ieta = ieta + 1
!!                                    end do
!!                                    gl%visco%term_num1_eta(nrsubst) = ieta - 1
!!
!!                                    if (gl%substcasnr(nrsubst)=="106-97-8") then !butane
!!                                        !read complex exponential terms
!!                                        do
!!                                            read(fld_unit,*,iostat = err_read) gl%visco%coeff(nrsubst,ieta), gl%visco%texpo(nrsubst,ieta), gl%visco%dexpo(nrsubst,ieta), gl%visco%pi_VS7(nrsubst,ieta), gl%visco%li_VS7(nrsubst,ieta)
!!                                            if (err_read .ne. 0) exit
!!                                            ieta = ieta + 1
!!                                        end do
!!
!!                                        gl%visco%term_com_expo_eta(nrsubst) = ieta - 1 - gl%visco%term_num1_eta(nrsubst)   !number of complex exponential terms
!!                                        gl%visco%term_num1_eta(nrsubst) = ieta - 1  !number of simple terms (polynomial + exponential terms)
!!                                    else
!!                                        !read simple exponential terms
!!                                        do
!!                                            read(fld_unit,*,iostat = err_read) gl%visco%coeff(nrsubst,ieta), gl%visco%texpo(nrsubst,ieta), gl%visco%dexpo(nrsubst,ieta), gl%visco%pi_VS7(nrsubst,ieta), gl%visco%li_VS7(nrsubst,ieta)
!!                                            if (err_read .ne. 0) exit
!!                                            ieta = ieta + 1
!!                                        end do
!!
!!                                        gl%visco%term_expo_eta(nrsubst) = ieta - 1 - gl%visco%term_num1_eta(nrsubst)   !number of exponential terms
!!                                        gl%visco%term_num1_eta(nrsubst) = ieta - 1  !number of simple terms (polynomial + exponential terms)
!!
!!                                    end if
!!                                    !read complex terms
!!                                    do
!!                                        read(fld_unit,*,iostat = err_read) gl%visco%coeff(nrsubst,ieta), gl%visco%texpo(nrsubst,ieta), gl%visco%dexpo(nrsubst,ieta), gl%visco%pi_VS7(nrsubst,ieta), gl%visco%li_VS7(nrsubst,ieta), gl%visco%eta_VS7(nrsubst,ieta), gl%visco%beta_VS7(nrsubst,ieta), gl%visco%par1(nrsubst,ieta), gl%visco%par2(nrsubst,ieta)
!!                                        if (err_read .ne. 0) exit
!!                                        ieta = ieta + 1
!!                                    end do
!!                                    gl%visco%term_num2_eta(nrsubst) = ieta - 1 - gl%visco%term_num1_eta(nrsubst)
!!
!!                                END IF
!!
!!                                exit loopETA
!!                            enddo loopETA
!!                            Goto 1337
!!
!!
!!                        elseif (dummygen(1:3) == '#PS') then
!!                            read(fld_unit,*) dummy
!!                            gl%vptype(nrsubst) = ichar(dummy(3:3)) - 48 !get the number at 3. place
!!                            !paragraph with comments can be ignored (every line starts with "?")
!!                            goover4:            do !until dummy starts with "!"
!!                                read(fld_unit,*) dummy
!!                                if (dummy(1:1) == '!') then
!!                                    exit goover4
!!                                endif
!!                            enddo goover4
!!                            !paragraph with doubled parameters can be ignored (always 4 lines)
!!                            do j=1,4
!!                                read(fld_unit,*) dummy
!!                            enddo
!!                            !read reducing parameters for vapour pressure:
!!                            read(fld_unit,*) gl%tvpred(nrsubst), gl%vpred(nrsubst)
!!                            gl%vpred(nrsubst) = gl%vpred(nrsubst)*factor_inv ![MPa]
!!                            !read coeffs ande exps:
!!                            read(fld_unit,*) gl%nvpcoeff(nrsubst)
!!
!!                            do j=1,gl%nvpcoeff(nrsubst)
!!                                read(fld_unit,*) gl%vpcoeff(j,nrsubst),gl%vpexp(j,nrsubst)
!!                            enddo
!!                            Goto 1337
!!                            !saturated liquid density:
!!                        elseif (dummygen(1:3) == '#DL') then
!!                            read(fld_unit,*) dummy
!!                            gl%dltype(nrsubst) = ichar(dummy(3:3)) - 48 !get the number at 3. place
!!                            !paragraph with comments can be ignored (every line starts with "?")
!!                            goover5:            do !until dummy starts with "!"
!!                                read(fld_unit,*) dummy
!!                                if (dummy(1:1) == '!') then
!!                                    exit goover5
!!                                endif
!!                            enddo goover5
!!                            !paragraph with doubled parameters can be ignored (always 4 lines)
!!                            do j=1,4
!!                                read(fld_unit,*) dummy
!!                            enddo
!!                            !read reducing parameters for saturated liquid density:
!!                            read(fld_unit,*) gl%tdlred(nrsubst), gl%dlred(nrsubst)
!!                            gl%dlred(nrsubst) = gl%dlred(nrsubst)*gl%Factor ![mol/m設
!!                            !read coeffs ande exps:
!!                            read(fld_unit,*) gl%ndlcoeff(nrsubst)
!!                            do j=1,gl%ndlcoeff(nrsubst)
!!                                read(fld_unit,*) gl%dlcoeff(j,nrsubst),gl%dlexp(j,nrsubst)
!!                            enddo
!!                            Goto 1337
!!                            !saturated vapour density:
!!                        elseif (dummygen(1:3) == '#DV') then
!!                            read(fld_unit,*) dummy
!!                            gl%dvtype(nrsubst) = ichar(dummy(3:3)) - 48 !get the number at 3. place
!!                            !paragraph with comments can be ignored (every line starts with "?")
!!                            goover6:            do !until dummy starts with "!"
!!                                read(fld_unit,*) dummy
!!                                if (dummy(1:1) == '!') then
!!                                    exit goover6
!!                                endif
!!                            enddo goover6
!!                            !paragraph with doubled parameters can be ignored (always 4 lines)
!!                            do j=1,4
!!                                read(fld_unit,*) dummy
!!                            enddo
!!                            !read reducing parameters for saturated vapour density:
!!                            read(fld_unit,*) gl%tdvred(nrsubst), gl%dvred(nrsubst)
!!                            gl%dvred(nrsubst) = gl%dvred(nrsubst)*gl%Factor ![mol/m設
!!                            !read coeffs ande exps:
!!                            read(fld_unit,*) gl%ndvcoeff(nrsubst)
!!                            do j=1,gl%ndvcoeff(nrsubst)
!!                                read(fld_unit,*) gl%dvcoeff(j,nrsubst),gl%dvexp(j,nrsubst)
!!                            enddo
!!                            Goto 1337
!!                            !melting pressure:
!!                        elseif (dummygen(1:4) == '#MLT') then
!!                            read(fld_unit,*) dummy
!!                            gl%pmelttype(nrsubst) = ichar(dummy(3:3)) - 48 !get the number at 3. place
!!                            !paragraph with comments can be ignored (every line starts with "?")
!!                            if (gl%pmelttype(nrsubst) == 39) then
!!
!!                                goover7:              do !until dummy starts with "!"
!!                                    gl%meltexist(nrsubst) = .true.
!!                                    read(fld_unit,*) dummy
!!                                    call uppertolower_char(dummy, len(dummy))
!!                                    if (dummy(1:1) == '!') then
!!                                        exit goover7
!!                                    elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                                        read(fld_unit,'(A)')dummy
!!                                        do while (trim(dummy(1:5)) /= '!````')
!!                                            gl%litref%lit_ref_mlt(nrsubst) = trim(gl%litref%lit_ref_mlt(nrsubst)) // dummy(2:)
!!                                            read(fld_unit,'(A)')dummy
!!                                        end do
!!                                        backspace(fld_unit)
!!                                    endif
!!                                enddo goover7
!!                                !paragraph with doubled parameters can be ignored (always 4 lines)
!!                                do j=1,4
!!                                    read(fld_unit,*) dummy
!!                                enddo
!!                                !read reducing parameters for melting pressure:
!!                                read(fld_unit,*) gl%tpmeltred(nrsubst), gl%pmeltred(nrsubst)
!!                                gl%pmeltred(nrsubst) = gl%pmeltred(nrsubst)*1.d-3 ![MPa]
!!                                !read coeffs and exponents:
!!                                read(fld_unit,*) gl%npmeltcoeff1(nrsubst), gl%npmeltcoeff2(nrsubst), gl%npmeltcoeff3(nrsubst)
!!                                do j=1,(gl%npmeltcoeff1(nrsubst) + gl%npmeltcoeff2(nrsubst) + gl%npmeltcoeff3(nrsubst))
!!                                    read(fld_unit,*) gl%pmeltcoeff(j,nrsubst),gl%pmeltexp(j,nrsubst)
!!                                enddo
!!                            elseif (gl%pmelttype(nrsubst) == 24) then                                        !fluid is water
!!
!!                                !  paragraph with comments can be ignored (every line starts with "?")
!!                                goover77:              do !until dummy starts with "!"
!!                                    read(fld_unit,*) dummy
!!                                    call uppertolower_char(dummy, len(dummy))
!!                                    if (dummy(1:1) == '!') then
!!                                        exit goover77
!!                                    elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                                        read(fld_unit,'(A)')dummy
!!                                        do while (trim(dummy(1:5)) /= '!````')
!!                                            gl%litref%lit_ref_mlt(nrsubst) = trim(gl%litref%lit_ref_mlt(nrsubst)) // dummy(2:)
!!                                            read(fld_unit,'(A)')dummy
!!                                        end do
!!                                        backspace(fld_unit)
!!                                    endif
!!                                enddo goover77
!!                                !paragraph with doubled parameters can be ignored (always 4 lines)
!!                                read(fld_unit,*) gl%pmeltmintemp(nrsubst)
!!                                read(fld_unit,*) gl%pmeltmaxtemp(nrsubst)
!!                                read(fld_unit,*) dummy
!!                                read(fld_unit,*) dummy
!!                                !read reducing parameters for melting pressure:
!!                                read(fld_unit,*) gl%tpmeltred(nrsubst), gl%pmeltred(nrsubst)
!!                                gl%pmeltred(nrsubst) = gl%pmeltred(nrsubst)*1.d-3 ![MPa]
!!
!!                                !read coeffs ande exps:
!!                                read(fld_unit,*) gl%npmeltcoeff1(nrsubst), gl%npmeltcoeff2(nrsubst), gl%npmeltcoeff3(nrsubst),&
!!                                    gl%npmeltcoeff4(nrsubst), gl%npmeltcoeff5(nrsubst)
!!                                do j=1,(gl%npmeltcoeff1(nrsubst) + gl%npmeltcoeff2(nrsubst) + gl%npmeltcoeff3(nrsubst)&
!!                                    +gl%npmeltcoeff4(nrsubst)+ gl%npmeltcoeff5(nrsubst))
!!                                    read(fld_unit,*) gl%pmeltcoeff(j,nrsubst),gl%pmeltexp(j,nrsubst)
!!                                enddo
!!                            elseif (gl%pmelttype(nrsubst) == 20) then                                        !fluid is heavy water
!!                                goover111:              do !until dummy starts with "!"
!!                                    read(fld_unit,*) dummy
!!                                    call uppertolower_char(dummy, len(dummy))
!!                                    if (dummy(1:1) == '!') then
!!                                        exit goover111
!!                                    elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                                        read(fld_unit,'(A)')dummy
!!                                        do while (trim(dummy(1:5)) /= '!````')
!!                                            gl%litref%lit_ref_mlt(nrsubst) = trim(gl%litref%lit_ref_mlt(nrsubst)) // dummy(2:)
!!                                            read(fld_unit,'(A)')dummy
!!                                        end do
!!                                        backspace(fld_unit)
!!                                    endif
!!                                enddo goover111
!!                                !paragraph with doubled parameters can be ignored (always 4 lines)
!!                                read(fld_unit,*) gl%pmeltmintemp(nrsubst)
!!                                read(fld_unit,*) gl%pmeltmaxtemp(nrsubst)
!!                                read(fld_unit,*) dummy
!!                                read(fld_unit,*) dummy
!!                                read(fld_unit,*) gl%tpmeltred(nrsubst), gl%pmeltred(nrsubst)
!!
!!                            else
!!                                goover88:              do !until dummy starts with "!"
!!                                    read(fld_unit,*) dummy
!!                                    call uppertolower_char(dummy, len(dummy))
!!                                    if (dummy(1:1) == '!') then
!!                                        exit goover88
!!                                    elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                                        read(fld_unit,'(A)')dummy
!!                                        do while (trim(dummy(1:5)) /= '!````')
!!                                            gl%litref%lit_ref_mlt(nrsubst) = trim(gl%litref%lit_ref_mlt(nrsubst)) // dummy(2:)
!!                                            read(fld_unit,'(A)')dummy
!!                                        end do
!!                                        backspace(fld_unit)
!!                                    endif
!!                                enddo goover88
!!                                !paragraph with doubled parameters can be ignored (always 4 lines)
!!                                read(fld_unit,*) gl%pmeltmintemp(nrsubst)
!!                                read(fld_unit,*) gl%pmeltmaxtemp(nrsubst)
!!                                read(fld_unit,*) dummy
!!                                read(fld_unit,*) dummy
!!                                read(fld_unit,*) gl%tpmeltred(nrsubst), gl%pmeltred(nrsubst)
!!                                !Andreas, September 2014, I guess the 1.D-30 cant be right? Replaced by 1.d-3
!!                                !pmeltred(nrsubst) = pmeltred(nrsubst)*1.d-30 ![MPa]
!!                                gl%pmeltred(nrsubst) = gl%pmeltred(nrsubst)*1.d-3 ![MPa]
!!                                read(fld_unit,*) gl%npmeltcoeff1(nrsubst), gl%npmeltcoeff2(nrsubst), gl%npmeltcoeff3(nrsubst),&
!!                                    gl%npmeltcoeff4(nrsubst), gl%npmeltcoeff5(nrsubst)
!!                                do j=1,(gl%npmeltcoeff1(nrsubst) + gl%npmeltcoeff2(nrsubst) + gl%npmeltcoeff3(nrsubst)&
!!                                    +gl%npmeltcoeff4(nrsubst)+ gl%npmeltcoeff5(nrsubst))
!!                                    read(fld_unit,*) gl%pmeltcoeff(j,nrsubst),gl%pmeltexp(j,nrsubst)
!!                                enddo
!!                            end if
!!                            Goto 1337
!!                            !sublimation pressure:
!!                        elseif (dummygen(1:4) == '#SBL') then
!!                            read(fld_unit,*) dummy
!!                            gl%psubtype(nrsubst) = ichar(dummy(3:3)) - 48 !get the number at 3. place
!!                            !paragraph with comments can be ignored (every line starts with "?")
!!                            goover8:            do !until dummy starts with "!"
!!                                read(fld_unit,*) dummy
!!                                call uppertolower_char(dummy, len(dummy))
!!                                if (dummy(1:1) == '!') then
!!                                    exit goover8
!!                                elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                                    read(fld_unit,'(A)')dummy
!!                                    do while (trim(dummy(1:5)) /= '!````')
!!                                        gl%litref%lit_ref_sbl(nrsubst) = trim(gl%litref%lit_ref_sbl(nrsubst)) // dummy(2:)
!!                                        read(fld_unit,'(A)')dummy
!!                                    end do
!!                                    backspace(fld_unit)
!!                                endif
!!                            enddo goover8
!!                            !paragraph with doubled parameters can be ignored (always 4 lines)
!!                            do j=1,4
!!                                read(fld_unit,*) dummy
!!                            enddo
!!                            !read reducing parameters for sublimation pressure:
!!                            read(fld_unit,*) gl%tpsubred(nrsubst), gl%psubred(nrsubst)
!!                            gl%psubred(nrsubst) = gl%psubred(nrsubst)*factor_inv ![MPa]
!!                            !read coeffs ande exps:
!!                            read(fld_unit,*) gl%npsubcoeff1(nrsubst), gl%npsubcoeff2(nrsubst), gl%npsubcoeff3(nrsubst)
!!                            do j=1,(gl%npsubcoeff1(nrsubst) + gl%npsubcoeff2(nrsubst) + gl%npsubcoeff3(nrsubst))
!!                                read(fld_unit,*) gl%psubcoeff(j,nrsubst),gl%psubexp(j,nrsubst)
!!                            enddo
!!                            Goto 1337
!!                            !for PSRK calculation (Theresa):
!!                            !elseif (dummygen(1:3) == 'UFC') then
!!                            !!goover9:            do !until dummy starts with "!"
!!                            !!                        read(fld_unit,*) dummy
!!                            !!                        if (dummy(1:1) == '!') then
!!                            !!                            exit goover9
!!                            !!                        endif
!!                            !!                    enddo goover9
!!                            !read(fld_unit,*) gl%ccoeff(nrsubst,1), gl%ccoeff(nrsubst,2), gl%ccoeff(nrsubst,3)
!!                            !ufcloop1:           do
!!                            !    read(fld_unit,*) groupnr, subgroupnr, amount, Qkcoeff !if no different subgroups exist then subgroupnr should be 1
!!                            !    if (int(groupnr) == 0) then !if end of ufc part in file then exit loop
!!                            !        exit ufcloop1
!!                            !    else
!!                            !        new_group = .true.
!!                            !        new_subgroup = .true.
!!                            !        goover10:                   do j=1, gl%ngroup_nsubgroup !check, if the same group and subgroup is already in ufc matrix
!!                            !            if (int(gl%ufc(j,1)) == int(groupnr)) then
!!                            !                new_group = .false.
!!                            !                if (int(gl%ufc(j,2)) == int(subgroupnr)) then
!!                            !                    new_subgroup = .false.
!!                            !                    line = j !line in matrix (groupnr + subgroupnr-1 build one line)
!!                            !                    exit goover10
!!                            !                endif
!!                            !            endif
!!                            !        enddo goover10
!!                            !        if ((new_group).or.(new_subgroup)) then !add new group and subgroup to ufc matrix
!!                            !            gl%ngroup_nsubgroup = gl%ngroup_nsubgroup+1
!!                            !            line = gl%ngroup_nsubgroup !new groupnr or new subgroupnr so add new line to ufc matrix
!!                            !            gl%ufc(line,1) = groupnr
!!                            !            gl%ufc(line,2) = subgroupnr
!!                            !            gl%ufc(line, 3) = Qkcoeff
!!                            !            Qkcoeffs(line) = Qkcoeff
!!                            !        endif
!!                            !        gl%ufc(line, 3 + nrsubst) = gl%ufc(line, 3 + nrsubst) + amount
!!                            !    endif
!!                            !enddo ufcloop1
!!                            !Goto 1337
!!                            !End of file
!!                        elseif (dummygen(1:4) == '@END') then
!!                            exit loopread
!!                        endif
!!                    enddo loopread    !belongs to the "do" at "loopread"
!!
!!                    !if (trim(components(nrsubst)) == 'pec5') rhoc_est_inp = 550.d0
!!                    !if (trim(components(nrsubst)) == 'pec7') rhoc_est_inp = 410.d0
!!                    !call find_crit_tpd (tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errval, nrsubst)
!!
!!                    !Close file
!!                    close(fld_unit)
!!
!!                    !change coefficients of the ideal gas if cp0red = 1.d0 (instead of the molar gas constant)
!!                    if((gl%cp0red(nrsubst)-1.d0) < 1.d-8) then
!!                        do j=1,(gl%ncp0poly(nrsubst)+gl%ncp0pl(nrsubst))
!!                            gl%cp0coeff(j,nrsubst)=gl%cp0coeff(j,nrsubst)/gl%req(nrsubst)
!!                            gl%cp0red(nrsubst)=gl%req(nrsubst)
!!                        end do
!!                    end if
!!
!!                    !change unit of molar mass from g/mol into kg/mol:
!!                    gl%wm(nrsubst) = gl%wm(nrsubst) *factor_inv           ![kg/mol]
!!                    wmcur(nrsubst) = wmcur(nrsubst) *factor_inv    ![kg/mol]
!!
!!
!!
!!                else
!!                    errorflag = -7879
!!                    return
!!                    !write(*,*) 'Error: occured when opening the fluid-file: ', trim(filename)
!!                endif
!!
!!
!!
!!            end if !substance name was not ''
!!
!!
!!            gl%rhoredmix = gl%rhored(nrsubst) ! needed, because the module variable hasn't been set, yet
!!            gl%Tredmix = gl%Tred(nrsubst)     ! needed, because the module variable hasn't been set, yet
!!
!!For AGA8, the critical pressure cannot be calculated. Therefore, we use the GERG-2008 as workaround
!!            !Trieu, Moni: 09/2019
!!            if(gl%eq_type(nrsubst) == 7) then
!!                gl_aga8 = .true.
!!                gl%eq_type(nrsubst) = 1
!!	end if
!!
!!            gl%pc(nrsubst) = P_CALC(gl, gl%tc(nrsubst), gl%rhoc(nrsubst),nrsubst) * gl%factor  !critical pressure calculated from EOS [kPa]
!!            gl%pc(nrsubst) = gl%pc(nrsubst)*factor_inv       ![MPa]
!!
!!            !Change back to AGA8 if that was the original input
!!            !Trieu, Moni: 09/2019
!!            if(gl_aga8) then
!!                gl%eq_type(nrsubst) = 7
!!                gl_aga8 = .false.
!!	end if
!!
!!            !calculate acentric factor, if error it takes it from fluid file
!!
!!            !gl%accen(nrsubst) = OMEGA_CALC(gl,nrsubst)
!!
!!        elseif ((gl%Eq_type(nrsubst) == 2) .or. (gl%Eq_type(nrsubst) == 8) .or. (gl%Eq_type(nrsubst) == 81) .or. (gl%Eq_type(nrsubst) == 9)) then     !SRK is used for that fluid. Andreas Oktober 2012; SRK is needed for Costald equation. Monika 12/2014
!!
!!
!!            if (gl%components(nrsubst) == '') then
!!                exit loopsubst
!!            end if
!!
!!
!!            !filename = trim (pathforsingles) // trim (components(nrsubst)) // '.fld'  ! check for standard case: fluid file (fld)
!!            !changed for compatibility with gfortran
!!            !filename = trim(path) // 'SRK\FLUIDS_SRK\SRKFLUIDS.fld'  ! check for standard case: fluid file (fld)
!!            filename = trim(path) // 'srk/fluids_srk/srkfluids.fld'  ! check for standard case: fluid file (fld)
!!            inquire (file = filename, exist = exists)
!!            if (.not.exists) then   ! no suitable ppf found
!!                errorflag = -8878
!!                return
!!            end if
!!
!!            open(newunit=srk_unit, file=trim(filename), status='old', action='read', iostat=error)
!!            if (error /= 0) then
!!                errorflag = -7885
!!                return
!!            end if
!!
!!            !Skip the first line in the file
!!            read(srk_unit,*,end=900) dummy_SRK
!!            !CASNR_SRK = ''
!!            Do i = 1, 1000
!!                read(srk_unit,'(A)',end=900) dummy_SRK
!!                if (dummy_SRK(1:4) == '@END') exit
!!                !**************************************************************************************************************
!!                !MESSER SRK FLUIDFILE FORMAT IS NOW THE SAME AS OUR FLUID FILE FORMAT EXCEPT A - G ideal part parameters
!!                !Andreas, July 2014
!!                !**************************************************************************************************************
!!                !New information in srkfluids.fld: Mathias Copeman Parameters C1, C2, C3 and number of groups and group numbers for UNIFAC groups (needed for PSRK)
!!                !Andreas J輍er, January 2017
!!                !Try to read in additional information first.
!!                !SH, Init of Variables for new model was missing
!!                Cji_cubic_read = 0.D0
!!                nr_of_groups_i_read = 0
!!                subgroups_ik_list_read = 0
!!
!!                !WARNING: AT THE MOMENT ONLY 30 GROUPS ARE READ IN. IF MORE GROUPS ARE NEEDED, THE FOLLOWING READ STATEMENT AS WELL AS THE FILE srkfluids.fld NEED TO BE ADJUSTED
!!                read(dummy_SRK,*,iostat= error) fluidnames_SRK, flname_SRK_alter, CASNR_SRK(nrsubst), MW_SRK, accen_SRK, pc_SRK, Tc_SRK, ptp_SRK, ttp_SRK, A_cv0_SRK,&
!!                    & B_cv0_SRK, C_cv0_SRK, D_cv0_SRK,E_cv0_SRK,F_cv0_SRK,G_cv0_SRK, Cji_cubic_read(1), Cji_cubic_read(2), &
!!                    & Cji_cubic_read(3), nr_of_groups_i_read, subgroups_ik_list_read(1), subgroups_ik_list_read(2), &
!!                    & subgroups_ik_list_read(3), subgroups_ik_list_read(4), subgroups_ik_list_read(5), subgroups_ik_list_read(6), &
!!                    & subgroups_ik_list_read(7), subgroups_ik_list_read(8), subgroups_ik_list_read(9), subgroups_ik_list_read(10), &
!!                    & subgroups_ik_list_read(11), subgroups_ik_list_read(12), subgroups_ik_list_read(13), subgroups_ik_list_read(14), &
!!                    & subgroups_ik_list_read(15), subgroups_ik_list_read(16), subgroups_ik_list_read(17), subgroups_ik_list_read(18), &
!!                    & subgroups_ik_list_read(19), subgroups_ik_list_read(20), subgroups_ik_list_read(21), subgroups_ik_list_read(22), &
!!                    & subgroups_ik_list_read(23), subgroups_ik_list_read(24), subgroups_ik_list_read(25), subgroups_ik_list_read(26), &
!!                    & subgroups_ik_list_read(27), subgroups_ik_list_read(28), subgroups_ik_list_read(29), subgroups_ik_list_read(30)
!!
!!                !if an error occured, proceed with the usual reading routine. Maybe this step is not needed and can be deleted as it seems that even if the read routine before fails, the existing values are nevertheless read. However, for safety reasons keep this at the moment.
!!                !SH:08/17:        !Messer format is different we definitely need this fallback code
!!                    gl%TMINFLUID(nrsubst) = ttp_PCSAFT
!!                    gl%TMAXFLUID(nrsubst) = 2000.d0  !dummy, see SRK
!!                if (error /= 0) then
!!                    error = 0
!!                    fluidnames_SRK = ''
!!                    read(dummy_SRK,*,iostat= error) fluidnames_SRK, flname_SRK_alter, CASNR_SRK(nrsubst), MW_SRK, accen_SRK, pc_SRK, Tc_SRK, ptp_SRK, ttp_SRK, A_cv0_SRK,&
!!                        B_cv0_SRK, C_cv0_SRK, D_cv0_SRK,E_cv0_SRK,F_cv0_SRK,G_cv0_SRK
!!                    !DEC$ IF DEFINED(MESSER)
!!                    !****************************************************************************************
!!                    if (error /= 0) then
!!                        !Messerformat begin
!!                        error = 0
!!                        fluidnames_SRK = ''
!!                        ptp_SRK = 0.d0
!!                        ttp_SRK = 0.d0
!!                        A_cv0_SRK = 0.d0
!!                        B_cv0_SRK = 0.d0
!!                        C_cv0_SRK = 0.d0
!!                        D_cv0_SRK = 0.d0
!!                        E_cv0_SRK = 0.d0
!!                        F_cv0_SRK = 0.d0
!!                        G_cv0_SRK = 0.d0
!!
!!                        !read(dummy_SRK,*,iostat= error) fluidnames_SRK, flname_SRK_alter, CASNR_SRK(nrsubst), Kompnr_SRK(nrsubst), accen_SRK, pc_SRK, Tc_SRK
!!                        return
!!
!!                    End if
!!                    !DEC$ ENDIF
!!                    !****************************************************************************************
!!
!!                end if
!!                continue
!!                if ((trim(fluidnames_SRK(1:30)) == trim(gl%components(nrsubst))) .OR. (trim(flname_SRK_alter(1:30)) == trim(gl%components(nrsubst)))  &
!!                    & .OR. (CASNR_SRK(nrsubst) == gl%components(nrsubst))) then
!!                    gl%substfullname(nrsubst) = fluidnames_SRK
!!                    gl%substshortname(nrsubst) = flname_SRK_alter
!!                    gl%substcasnr(nrsubst) = CASNR_SRK(nrsubst)
!!                    gl%tc(nrsubst) = tc_SRK
!!                    gl%pc(nrsubst) = pc_SRK
!!                    gl%tred(nrsubst) = gl%tc(nrsubst)
!!                    gl%accen(nrsubst) = accen_SRK
!!                    !Triple point information !Andreas, July 2014
!!                    gl%ptp(nrsubst) = ptp_SRK
!!                    gl%ttp(nrsubst) = ttp_SRK
!!                    !New variables for Mathias Copeman              !Andreas J輍er, January 2017
!!                    !Only use Mathias Copeman parameters if the PSRK is used (mix_type == 22)
!!                    if (gl%mix_type == 22) then
!!                        gl%Cji_cubic(:,nrsubst) = Cji_cubic_read
!!                    else
!!                        gl%Cji_cubic(:,nrsubst) = 0.D0
!!                    end if
!!                    if (gl%VLE_needed) then
!!                        !Andreas, November 2015: PURE SRK equation needs to be evaluated here!!
!!                        !Set tredmix to tc(nrsubst)
!!                        !THE NEXT LINE IS EXTREMELY IMPORTANT if the SRK is used in mixtures, because without this line the pure fluid is evaluated
!!                        !at the wrong temperature and thus the critical density and all other calculated densities are wrong if tredmix = tc(nrsubst) is missing!
!!                        !Andreas J輍er, September 2017
!!                        gl%tredmix = gl%tc(nrsubst)
!!                        gl%rhored(nrsubst) = rho_SRK (gl, gl%tc(nrsubst), gl%pc(nrsubst), 0, nrsubst)!tc(nrsubst), pc(nrsubst), 1,
!!                        gl%rhoc(nrsubst) = gl%rhored(nrsubst)
!!                        if ((gl%ttp(nrsubst) > 1.d-14) .and. (gl%ptp(nrsubst) > 1.d-14)) then !Sebastian: Theresa neu funktioniert nicht wenn ttp und ptp = 0
!!                            gl%rhotp(nrsubst) = rho_SRK (gl, gl%ttp(nrsubst), gl%ptp(nrsubst), 0, 1) !Theresa neu
!!                        end if
!!                        !Set back tredmix, Andreas September 2017
!!                        gl%tredmix = 1.D0
!!                    end if
!!                    gl%pmaxfluid(nrsubst) = 1000.d0    !Needed for PhaseDet_TD (mixtures) (Dummy value)
!!                    gl%tminfluid(nrsubst) = ttp_SRK       !Needed for Flash calculations, arbitrary value!! Andreas, May 2016
!!                    !Changed from 1.d0 to ttp because for PH and PS Flash calculations ran into problems. Sebastian, October 2016
!!                    gl%tmaxfluid(nrsubst) = 2000.d0    !Needed for Flash calculations, arbitrary value!! Andreas, May 2016
!!                    gl%wm(nrsubst) = MW_SRK
!!                    gl%REQ(nrsubst) = 8.3144598D0  !J / mol K
!!                    !TEST_AJ - For paper Bell & J輍er 2016 values for SRK similar to Stradi et al. have been chosen August 2016
!!                    !REQ(nrsubst) = 8.31451D0
!!                    gl%A_cv0(nrsubst)= A_cv0_SRK
!!                    gl%B_cv0(nrsubst)= B_cv0_SRK
!!                    gl%C_cv0(nrsubst)= C_cv0_SRK
!!                    gl%D_cv0(nrsubst)= D_cv0_SRK
!!                    gl%E_cv0(nrsubst)= E_cv0_SRK
!!                    gl%F_cv0(nrsubst)= F_cv0_SRK
!!                    gl%G_cv0(nrsubst)= G_cv0_SRK
!!                    if (((dabs(gl%B_cv0(nrsubst)) > 1.D-8) .and. (dabs(gl%C_cv0(nrsubst)) > 1.D-8))) gl%cpmodel(nrsubst) = .true.
!!                    gl%refstate(nrsubst) = "OT0"
!!                    gl%tref(nrsubst) = 298.15d0
!!                    gl%pref(nrsubst) = 0.101325d0  !"-1" shows that pref must be calculated later in another subroutine
!!                    gl%href(nrsubst) = 0.d0 ![J/kg]     -later changed with wm into [J/mol]
!!                    gl%sref(nrsubst) = 0.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
!!                    gl%tcp0red(nrsubst) = gl%tc(nrsubst)
!!                    gl%cp0red(nrsubst) = 8.3144598D0
!!                    gl%bi_SRK(nrsubst) = 0.08664D0 * gl%REQ(nrsubst) * gl%tc(nrsubst) / (gl%pc(nrsubst)*1.D6)
!!                    gl%rhomaxfluid(nrsubst) = 1.D0/gl%bi_SRK(nrsubst)
!!                    if (.not. gl%hold_limits) gl%ttp(nrsubst) = gl%ttp(nrsubst)*0.1d0    !small triple point temperature for the extrapolation of the VLE curves
!!                    !     (if "&" occurs in the input variable)
!!                    !New variables for group contribution methods   !Andreas J輍er, January 2017
!!                    gl%nr_of_groups_i(nrsubst) = nr_of_groups_i_read
!!                    gl%subgroups_ik_list(nrsubst,:) = subgroups_ik_list_read
!!                    if (gl%mix_type /= 12) then        !Necessary for using cubics in new mixture model 12, because elsewise the groups are counted twice
!!                        gl%nr_of_groups_mix = gl%nr_of_groups_mix + gl%nr_of_groups_i(nrsubst)
!!                    end if
!!                    exit
!!                end if
!!            End do
!!
!!            close(srk_unit)
!!
!!            !Moni, 2017-11-09: that might be dangerous for the Costald, RKM or eRKM
!!            if (dabs(gl%tc(nrsubst)) < 1.D-14) then
!!                errorflag = -8878
!!                return
!!            end if
!!
!!
!!            if ((gl%eq_type(nrsubst) == 8) .or. (gl%eq_type(nrsubst) == 81)) then
!!                !**************************************************************************************************************
!!                !in case of the revised Klosek-McKinley equation of state:
!!                !reads the k1 and k2 correction factors as well as values for the molar volumes Vi and different typical lng
!!                !compositions
!!                !Taken from: R. D. McCarty, J. Chem. Thermodynamics, 14:837-854 (1982)
!!                !**************************************************************************************************************
!!                !Christopher / January 2016
!!
!!
!!                !##############################################################################################
!!                !Saves the correction factor k1 and k2 in matrix k1(i,j) and k2(i,j), respectively
!!                !##############################################################################################
!!                filename = trim(path) // 'RKM\RKM_correction_factors.txt'
!!                open(newunit=rkm_unit, file=trim(filename), status='old', action='read', iostat=error)
!!                if (error /= 0) then
!!                    errorflag = -7893
!!                    return
!!                end if
!!
!!                read(rkm_unit,*) gl%dummy_RKM
!!                read(rkm_unit,*) gl%dummy_RKM
!!
!!                Do i = 0, 10
!!                    read(rkm_unit,*) gl%k1(i,0), gl%k1(i,1), gl%k1(i,2), gl%k1(i,3), gl%k1(i,4), gl%k1(i,5), gl%k1(i,6), gl%k1(i,7), gl%k1(i,8), gl%k1(i,9), gl%k1(i,10)
!!                end do
!!
!!                read(rkm_unit,*) gl%dummy_RKM
!!                read(rkm_unit,*) gl%dummy_RKM
!!
!!                Do i = 0, 10
!!                    read(rkm_unit,*) gl%k2(i,0), gl%k2(i,1), gl%k2(i,2), gl%k2(i,3), gl%k2(i,4), gl%k2(i,5), gl%k2(i,6), gl%k2(i,7), gl%k2(i,8), gl%k2(i,9), gl%k2(i,10)
!!                end do
!!
!!                close(rkm_unit)
!!
!!                gl%k1(1:10,1:10) = gl%k1(1:10,1:10) *factor_inv
!!                gl%k2(1:10,1:10) = gl%k2(1:10,1:10) *factor_inv
!!
!!                !##############################################################################################
!!                !Saves the given molar volumes in matrix V(i,j), the molar mass in RKM_M(i) and the list of
!!                !components in RKM_fluids(i)
!!                !##############################################################################################
!!
!!                filename = trim(path) // 'rkm\RKM_fluid_information.txt'
!!                open(newunit=rkminfo_unit, file=trim(filename), status='old', action='read', iostat=error)
!!                if (error /= 0) then
!!                    errorflag = -7893
!!                    return
!!                end if
!!
!!                !Skip the first line in the file
!!                read(rkminfo_unit,*) gl%dummy_RKM
!!                read(rkminfo_unit,*) gl%RKM_fluids(1,1),gl%RKM_fluids(2,1),gl%RKM_fluids(3,1),gl%RKM_fluids(4,1),gl%RKM_fluids(5,1),gl%RKM_fluids(6,1),gl%RKM_fluids(7,1),gl%RKM_fluids(8,1)
!!                read(rkminfo_unit,*) gl%RKM_fluids(1,2),gl%RKM_fluids(2,2),gl%RKM_fluids(3,2),gl%RKM_fluids(4,2),gl%RKM_fluids(5,2),gl%RKM_fluids(6,2),gl%RKM_fluids(7,2),gl%RKM_fluids(8,2)
!!                read(rkminfo_unit,*) gl%RKM_fluids(1,3),gl%RKM_fluids(2,3),gl%RKM_fluids(3,3),gl%RKM_fluids(4,3),gl%RKM_fluids(5,3),gl%RKM_fluids(6,3),gl%RKM_fluids(7,3),gl%RKM_fluids(8,3)
!!                read(rkminfo_unit,*) gl%dummy_RKM
!!                read(rkminfo_unit,*) gl%RKM_Tc(1),gl%RKM_Tc(2),gl%RKM_Tc(3),gl%RKM_Tc(4),gl%RKM_Tc(5),gl%RKM_Tc(6),gl%RKM_Tc(7),gl%RKM_Tc(8)
!!                read(rkminfo_unit,*) gl%dummy_RKM
!!                read(rkminfo_unit,*) gl%RKM_M(1),gl%RKM_M(2),gl%RKM_M(3),gl%RKM_M(4),gl%RKM_M(5),gl%RKM_M(6),gl%RKM_M(7),gl%RKM_M(8)
!!                read(rkminfo_unit,*) gl%dummy_RKM
!!                Do i = 0, 100
!!                    read(rkminfo_unit,*) gl%RKM_ps_me(i,0), gl%RKM_ps_me(i,1)
!!                    if (gl%RKM_ps_me(i,0) == 140.D0) then
!!                        exit
!!                    end if
!!                end do
!!                read(rkminfo_unit,*) gl%dummy_RKM
!!                Do i = 0, 100
!!                    read(rkminfo_unit,*) gl%RKM_V(i,0), gl%RKM_V(i,1), gl%RKM_V(i,2), gl%RKM_V(i,3), gl%RKM_V(i,4), gl%RKM_V(i,5), gl%RKM_V(i,6), gl%RKM_V(i,7), gl%RKM_V(i,8)
!!                    if (gl%RKM_V(i,0) == 140.D0) then
!!                        exit
!!                    end if
!!                end do
!!                close(rkminfo_unit)
!!
!!
!!                !##############################################################################################
!!                !Sort Mole Fractions
!!                !##############################################################################################
!!
!!                do i = 1, 8
!!                    do j = 1, 8
!!                        if ((gl%components(j) == gl%RKM_fluids(i,1)) .or. (gl%components(j) == gl%RKM_fluids(i,2)) .or. (gl%components(j) == gl%RKM_fluids(i,3))) then
!!                            gl%wm(j) = gl%RKM_M(i)
!!                            exit
!!                        end if
!!                    end do
!!                end do
!!
!!                gl%wm(:) = gl%wm(:) *factor_inv           ![kg/mol]
!!
!!                exit
!!
!!            end if
!!
!!            if (gl%eq_type(nrsubst) == 9) then
!!                !**************************************************************************************************************
!!                !in case of the COSTALD equation of state:
!!                !reads the Costald parameters omega and vstar
!!                !Taken from: R. W. Hankinson and G. H. Thomson, AIChE Journal, 25(4):653-663 (1979)
!!                !**************************************************************************************************************
!!                dummy_COSTALD = ''
!!
!!                filename = trim(path) // 'costald/costald.txt'
!!                !filename = trim(path) // 'COSTALD/COSTALD_API.txt'
!!                !filename = trim(path) // 'COSTALD/COSTALD_REF.txt'
!!                inquire (file = filename, exist = exists)
!!                if (.not.exists) then   ! no suitable text file found
!!                    errorflag = -8881
!!                    return
!!                end if
!!
!!                open(newunit=cost_unit, file=trim(filename), status='old', action='read', iostat=error)
!!                if (error /= 0) then
!!                    errorflag = -7891
!!                    return
!!                end if
!!
!!                !Skip the first line in the file
!!                read(cost_unit,*,end=900) dummy_COSTALD
!!
!!                Do i = 1, 1000
!!                    read(cost_unit,*,end=900) dummy_COSTALD
!!                    if (dummy_COSTALD(1:4) == '@END') then
!!                        errorflag =-8880
!!                        close(cost_unit)
!!                        return
!!                    end if
!!                    backspace(cost_unit)
!!
!!                    read(cost_unit,'(A33,A15,2F11.5)') fluidtait_COSTALD, CAS_COSTALD, omega_COSTALD, vstar_COSTALD
!!                    !read(cost_unit,*) fluidtait_COSTALD, CAS_COSTALD, omega_COSTALD, vstar_COSTALD
!!                    !if ((trim(fluidtait_COSTALD)) == trim(components(nrsubst))) then
!!                    if ((trim(CAS_COSTALD)) == trim(CASNR_SRK(nrsubst))) then
!!                        gl%omega(nrsubst) = omega_COSTALD
!!                        gl%vstar(nrsubst) = vstar_COSTALD
!!                        exit
!!                    end if
!!                End do
!!
!!                close(cost_unit)
!!            end if
!!
!!            !Fluid not found, return
!!            if (gl%tc(nrsubst) == 0.D0) then
!!                errorflag = -8878
!!                close(cost_unit)
!!                return
!!            end if
!!
!!
!!            if (gl%components(nrsubst) == 'oil') then
!!                !Thu 2014/05
!!                !cp0-model of Joback (1986)
!!                !derivation 1...2%, increase with more complexe molecule structures
!!                !source: Proofs W酺meatlas, Teil D1, p.21
!!
!!                gl%cpmodel = .true.
!!
!!                !set numbers of moleculegroup of substance
!!                N_m(1)=4         !CH3 (1)
!!                N_m(2)=26        !CH2 (2)
!!                N_m(3)=1         !CH  (3)
!!                N_m(4)=2         !C   (4)
!!
!!                dA(1) = 19.5d0
!!                dB(1) = -8.08d-3
!!                dC(1) = 1.53d-4
!!                dD(1) = -9.67d-8
!!                dA(2) = -0.909d0
!!                dB(2) = 9.5d-2
!!                dC(2) = -5.44d-5
!!                dD(2) = 1.19d-8
!!                dA(3) = -23.d0
!!                dB(3) = 2.04d-1
!!                dC(3) = -2.65d-4
!!                dD(3) = 1.2d-7
!!                dA(4) = -66.2d0
!!                dB(4) = 4.27d-01
!!                dC(4) = -6.41d-04
!!                dD(4) = 3.01d-7
!!
!!                do i=1,4
!!                    gl%cp0_A=gl%cp0_A+N_m(i)*dA(i)
!!                    gl%cp0_B=gl%cp0_B+N_m(i)*dB(i)
!!                    gl%cp0_C=gl%cp0_C+N_m(i)*dC(i)
!!                    gl%cp0_D=gl%cp0_D+N_m(i)*dD(i)
!!                end do
!!
!!                ! setting up common CPP module variables
!!                cpdummy=gl%eos_coeff%cppcheck(nrsubst)
!!                cpdummy(1:2) = 'CP'
!!                gl%ncp0poly(nrsubst) = 4
!!                gl%tcp0red(nrsubst) = 1.d0
!!                gl%cp0red(nrsubst) = 8.3144621d0
!!                gl%cp0exp(1,nrsubst) = 0.d0
!!                gl%cp0exp(2,nrsubst) = 1.d0
!!                gl%cp0exp(3,nrsubst) = 2.d0
!!                gl%cp0exp(4,nrsubst) = 3.d0
!!                gl%CP0COEFF(1,nrsubst) = (gl%cp0_A-37.93d0)/gl%cp0red(nrsubst)
!!                gl%CP0COEFF(2,nrsubst) = (gl%cp0_B+0.21d0)/gl%cp0red(nrsubst)
!!                gl%CP0COEFF(3,nrsubst) = (gl%cp0_C-3.91d-4)/gl%cp0red(nrsubst)
!!                gl%CP0COEFF(4,nrsubst) = (gl%cp0_D+2.06d-7)/gl%cp0red(nrsubst)
!!
!!            end if
!!
!!            !change unit of molar mass from g/mol into kg/mol:
!!            gl%wm(nrsubst) = gl%wm(nrsubst) *factor_inv           ![kg/mol]
!!            wmcur(nrsubst) = wmcur(nrsubst) *factor_inv    ![kg/mol]
!!
!!
!!        elseif (gl%Eq_type(nrsubst) == 6) then     !PCSAFT is used for that fluid. Henning Markgraf Dezember 2015
!!
!!            if (gl%components(nrsubst) == '') then
!!                exit loopsubst
!!            end if
!!
!!
!!            !filename = trim (pathforsingles) // trim (components(nrsubst)) // '.fld'  ! check for standard case: fluid file (fld)
!!            !changed for compatibility with gfortran
!!            !filename = trim(path) // 'SRK\FLUIDS_SRK\SRKFLUIDS.fld'  ! check for standard case: fluid file (fld)
!!            filename = trim(path) // 'pc-saft/pcsaftfluids.fld'  ! check for standard case: fluid file (fld)
!!            inquire (file = filename, exist = exists)
!!            if (.not.exists) then   ! no suitable ppf found
!!                errorflag = -8878
!!                return
!!            end if
!!
!!            open(newunit=pcsaft_unit, file=trim(filename), status='old', action='read', iostat=error)
!!            if (error /= 0) then
!!                errorflag = -7885
!!                return
!!            end if
!!
!!            !Skip the first line in the file
!!            read(pcsaft_unit,*,end=900) dummy_PCSAFT
!!
!!            Do i = 1, 1000
!!                io_error = 0
!!                read(pcsaft_unit,*,end=900, iostat=io_error) fluidnames_PCSAFT
!!                if (fluidnames_PCSAFT(1:4) == '@END') exit
!!                if (fluidnames_PCSAFT(1:1) .ne. '#') goto 5454
!!                backspace(pcsaft_unit)
!!
!!                io_error = 0
!!                read(pcsaft_unit,*,iostat=io_error) fluidnames_PCSAFT, flname_PCSAFT_alter, model_PCSAFT, CASNR_PCSAFT(nrsubst), MW_PCSAFT, m_PCSAFT, sig_PCSAFT, epsk_PCSAFT, epsabk_PCSAFT, kap_PCSAFT, Q_PCSAFTQ, nQ_PCSAFTQ, D_PCSAFTD, nD_PCSAFTD, accen_PCSAFT, dc_PCSAFT, Tc_PCSAFT, ptp_PCSAFT, ttp_PCSAFT, A_cv0_PCSAFT,&
!!                    B_cv0_PCSAFT, C_cv0_PCSAFT, D_cv0_PCSAFT,E_cv0_PCSAFT,F_cv0_PCSAFT,G_cv0_PCSAFT
!!                fluidnames_PCSAFT = fluidnames_PCSAFT(2:)
!!                if (io_error == 0) then
!!                    continue
!!                else
!!                    lengname = len(trim(fluidnames_PCSAFT))
!!                    lengname_ext = lengname + 25
!!                    backspace(pcsaft_unit)
!!                    read(pcsaft_unit,'(A<lengname_ext>)',end=900) fluidnames_PCSAFT
!!                    fluidnames_PCSAFT = fluidnames_PCSAFT(2:)
!!                    if(fluidnames_PCSAFT(lengname+1:lengname+1) == ' ') then    !fluid name consists of two words
!!                        backspace(pcsaft_unit)
!!                        read(pcsaft_unit,*,end=900) fluidnames_PCSAFT1, fluidnames_PCSAFT2
!!                        fluidnames_PCSAFT = trim(fluidnames_PCSAFT1) // ' ' // trim(fluidnames_PCSAFT2)
!!                        fluidnames_PCSAFT = fluidnames_PCSAFT(2:)
!!                        lengname = len(trim(fluidnames_PCSAFT))
!!                        backspace(pcsaft_unit)
!!                        io_error = 0
!!                        read(pcsaft_unit,*,end=900,iostat=io_error) fluidnames_PCSAFT1, fluidnames_PCSAFT2, flname_PCSAFT_alter, model_PCSAFT, CASNR_PCSAFT(nrsubst), MW_PCSAFT, m_PCSAFT, sig_PCSAFT, epsk_PCSAFT, epsabk_PCSAFT, kap_PCSAFT, Q_PCSAFTQ, nQ_PCSAFTQ, D_PCSAFTD, nD_PCSAFTD, accen_PCSAFT, dc_PCSAFT, Tc_PCSAFT, ptp_PCSAFT, ttp_PCSAFT, A_cv0_PCSAFT,&
!!                            B_cv0_PCSAFT, C_cv0_PCSAFT, D_cv0_PCSAFT,E_cv0_PCSAFT,F_cv0_PCSAFT,G_cv0_PCSAFT
!!                        if (io_error == 0) then
!!                            continue
!!                        else    !fluid name and alternative name consist of two words
!!                            backspace(pcsaft_unit)
!!                            read(pcsaft_unit,*,end=900) fluidnames_PCSAFT1, fluidnames_PCSAFT2, flname_PCSAFT_alter1, flname_PCSAFT_alter2, model_PCSAFT, CASNR_PCSAFT(nrsubst), MW_PCSAFT, m_PCSAFT, sig_PCSAFT, epsk_PCSAFT, epsabk_PCSAFT, kap_PCSAFT, Q_PCSAFTQ, nQ_PCSAFTQ, D_PCSAFTD, nD_PCSAFTD, accen_PCSAFT, dc_PCSAFT, Tc_PCSAFT, ptp_PCSAFT, ttp_PCSAFT, A_cv0_PCSAFT,&
!!                                B_cv0_PCSAFT, C_cv0_PCSAFT, D_cv0_PCSAFT,E_cv0_PCSAFT,F_cv0_PCSAFT,G_cv0_PCSAFT
!!                            flname_PCSAFT_alter = trim(flname_PCSAFT_alter1) // ' ' // trim(flname_PCSAFT_alter2)
!!                            continue
!!                        end if
!!                    else    !fluid name consists of one word, but alternative name consists of two words
!!                        backspace(pcsaft_unit)
!!                        read(pcsaft_unit,*,end=900) fluidnames_PCSAFT, flname_PCSAFT_alter1, flname_PCSAFT_alter2, model_PCSAFT, CASNR_PCSAFT(nrsubst), MW_PCSAFT, m_PCSAFT, sig_PCSAFT, epsk_PCSAFT, epsabk_PCSAFT, kap_PCSAFT, Q_PCSAFTQ, nQ_PCSAFTQ, D_PCSAFTD, nD_PCSAFTD, accen_PCSAFT, dc_PCSAFT, Tc_PCSAFT, ptp_PCSAFT, ttp_PCSAFT, A_cv0_PCSAFT,&
!!                            B_cv0_PCSAFT, C_cv0_PCSAFT, D_cv0_PCSAFT,E_cv0_PCSAFT,F_cv0_PCSAFT,G_cv0_PCSAFT
!!                        flname_PCSAFT_alter = trim(flname_PCSAFT_alter1) // ' ' // trim(flname_PCSAFT_alter2)
!!                        fluidnames_PCSAFT = fluidnames_PCSAFT(2:)
!!                        continue
!!                    end if
!!                end if
!!
!!
!!                if ((trim(fluidnames_PCSAFT(1:30)) == trim(gl%components(nrsubst))) .OR. (trim(flname_PCSAFT_alter(1:25)) == trim(gl%components(nrsubst)))  &
!!                    & .OR. (CASNR_PCSAFT(nrsubst) == gl%components(nrsubst))) then
!!
!!
!!                    if(CASNR_PCSAFT(nrsubst) == '999-999-999') then
!!                        gl%Factor = 1.D0
!!                        gl%factorpress=1.D6
!!                        gl%factortrans=1.D0
!!                        gl%factorrbwr=1.D0
!!                        gl%factor_VS5eta=1.D0
!!                        gl%factor_VS5slj=1.D0
!!                    else                         !for calculating real fluids
!!                        gl%Factor = 1.D3
!!                        gl%factorpress=1.D0
!!                        gl%factortrans=1.D6
!!                        gl%factorrbwr=1.D2
!!                        gl%factor_VS5eta=1.D5
!!                        gl%factor_VS5slj=1.D1
!!                    end if
!!
!!                    gl%SAFTmodel_assoc(nrsubst) = model_PCSAFT(2:2)
!!                    read (model_PCSAFT(1:1),*) gl%n_sites(nrsubst)
!!
!!                    !leng = len(trim(components(2)))
!!                    !n_sites(2) = components(2)(leng-1:leng-1)
!!
!!                    !gl%tc(nrsubst) = tc_PCSAFT
!!                    !gl%rhoc(nrsubst) = dc_PCSAFT !rho_SRK (gl, gl%tc(nrsubst), gl%pc(nrsubst), 0, nrsubst)
!!                    gl%accen(nrsubst) = accen_PCSAFT
!!                    gl%tred(nrsubst) = 1.d0 ! derivatives with respect to temperature can be used but need to be transformed
!!                    gl%rhored(nrsubst) = 1.d0 ! derivatives with respect to density can be used
!!                    gl%Eq_type(nrsubst) = 6
!!                    gl%Mix_type = 6
!!                    gl%pmaxfluid(nrsubst) = 1000.d0    !Needed for PhaseDet_TD (mixtures) (Dummy value)
!!                    gl%wm(nrsubst) = MW_PCSAFT/gl%factor
!!                    if (tc_PCSAFT < 5.d0) then
!!                        gl%REQ(nrsubst) = 1.D0
!!                    else
!!                        gl%REQ(nrsubst) = 8.31451D0
!!                    end if
!!                    gl%A_cv0(nrsubst)= A_cv0_PCSAFT
!!                    gl%B_cv0(nrsubst)= B_cv0_PCSAFT
!!                    gl%C_cv0(nrsubst)= C_cv0_PCSAFT
!!                    gl%D_cv0(nrsubst)= D_cv0_PCSAFT
!!                    gl%E_cv0(nrsubst)= E_cv0_PCSAFT
!!                    gl%F_cv0(nrsubst)= F_cv0_PCSAFT
!!                    gl%G_cv0(nrsubst)= G_cv0_PCSAFT
!!                    if (((dabs(gl%B_cv0(nrsubst)) > 1.D-8) .and. (dabs(gl%C_cv0(nrsubst)) > 1.D-8))) gl%cpmodel(nrsubst) = .true.
!!                    gl%refstate(nrsubst) = "OT0"
!!                    gl%tref(nrsubst) = 298.15d0
!!                    gl%pref(nrsubst) = 0.101325d0  !"-1" shows that pref must be calculated later in another subroutine
!!                    gl%href(nrsubst) = 0.d0 ![J/kg]     -later changed with wm into [J/mol]
!!                    gl%sref(nrsubst) = 0.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
!!                    gl%tcp0red(nrsubst) = gl%tc(nrsubst)
!!                    gl%cp0red(nrsubst) = 8.31451D0
!!                    !Triple point information
!!                    gl%ptp(nrsubst) = ptp_PCSAFT
!!                    gl%ttp(nrsubst) = ttp_PCSAFT
!!                    !rhotp(nrsubst) = rho_pt (ttp(nrsubst), ptp(nrsubst), 0, 1)
!!                    !PC-SAFT parameters m, sigma, epsilon/k
!!                    gl%mPCSAFT(nrsubst) = m_PCSAFT
!!                    gl%sigPCSAFT(nrsubst) = sig_PCSAFT
!!                    gl%epskPCSAFT(nrsubst) = epsk_PCSAFT
!!                    gl%epsabkPCSAFT(nrsubst) = epsabk_PCSAFT
!!                    gl%kabPCSAFT(nrsubst)= kap_PCSAFT
!!                    gl%QPCSAFTQ(nrsubst) = Q_PCSAFTQ
!!                    gl%MyPCSAFTD(nrsubst) = D_PCSAFTD
!!                    gl%nPCSAFTQ(nrsubst) = nQ_PCSAFTQ !dummy, unklar woher der Wert kommt
!!                    gl%nPCSAFTD(nrsubst) = nD_PCSAFTD !dummy, unklar woher der Wert kommt
!!                    gl%rhoredmix = 1.d0
!!
!!                    !special case elongation given
!!                    if (m_PCSAFT < 0.d0 .and. (Q_PCSAFTQ > 0.d0 .or. (Q_PCSAFTQ == 0.d0 .and. D_PCSAFTD == 0.d0))) then
!!                        gl%mPCSAFT(nrsubst) = 1.d0 + 0.2177d0*dabs(m_PCSAFT/sig_PCSAFT) + 3.1498d0*dabs(m_PCSAFT/sig_PCSAFT)**2 &
!!                            & - 3.6738d0*dabs(m_PCSAFT/sig_PCSAFT)**3 + 1.3063d0*dabs(m_PCSAFT/sig_PCSAFT)**4
!!                    else if (m_PCSAFT < 0.d0 .and. D_PCSAFTD > 0.d0) then
!!                        gl%mPCSAFT(nrsubst) = 1.d0 + 0.1795d0*dabs(m_PCSAFT/sig_PCSAFT) + 3.3283d0*(m_PCSAFT/sig_PCSAFT)**2 &
!!                            & - 3.8855d0*(dabs(m_PCSAFT/sig_PCSAFT))**3 + 1.3777d0*(m_PCSAFT/sig_PCSAFT)**4
!!                    end if
!!
!!
!!                    !determine true critical point of EOS
!!                    tc_est_inp = tc_PCSAFT
!!                    rhoc_est_inp = dc_PCSAFT*gl%factor
!!                    !if (trim(gl%components(nrsubst)) == 'pec5') rhoc_est_inp = 550.d0
!!                    !if (trim(gl%components(nrsubst)) == 'pec7') rhoc_est_inp = 410.d0
!!
!!                    call find_crit_tpd (gl,tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errval, nrsubst)
!!
!!
!!                    riter = 1
!!                    rhoc_est_inp = dc_PCSAFT*gl%factor
!!                    do while (errval .ne. 0)
!!                        errval = 0
!!                        rhoc_est_inp = rhoc_est_inp*1.02d0
!!                        call find_crit_tpd (gl,tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errval, nrsubst)
!!                        riter = riter + 1
!!                        if (errval == 0 .or. riter > 105) exit
!!                    end do
!!
!!                    riter = 1
!!                    rhoc_est_inp = dc_PCSAFT*gl%factor
!!                    do while (errval  .ne. 0)
!!                        errval = 0
!!                        rhoc_est_inp = rhoc_est_inp*0.98d0
!!                        call find_crit_tpd (gl,tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errval, nrsubst)
!!                        riter = riter + 1
!!                        if (errval == 0 .or. riter > 105) exit
!!                    end do
!!
!!                    if (errval .ne. 0) then  !calculate the critical parameters from the Lennard-Jones fluid if the given ones are too far off
!!                        errval = 0
!!                        !pc_est_inp = 0.13d0 * 1.3806504d-23 * epsk_PCSAFT/(1.d-30*sig_PCSAFT**3)*1.d-6
!!                        tc_est_inp = 1.32d0 * epsk_PCSAFT
!!                        rhoc_est_inp = 0.31d0 / (6.02214d-7*sig_PCSAFT**3) /2.d0    !TE /2.d0
!!                        call find_crit_tpd (gl,tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errval, nrsubst)
!!                        if (errval .ne. 0) then
!!                            errorflag = errval
!!                            return
!!                        end if
!!                    end if
!!
!!                    gl%tc(nrsubst) = tc_eos_
!!                    gl%rhoc(nrsubst) = rhoc_eos_
!!                    gl%pc(nrsubst) = pc_eos_
!!
!!                    if (gl%ttp(nrsubst) < 1.d-13) then   !calculate the triple point temperaturefrom the Lennard-Jones fluid if there is nothing given
!!                        gl%ttp(nrsubst) = 0.66d0 * epsk_PCSAFT
!!                        gl%rhotp(nrsubst) = 0.865d0 / (6.02214d-7*sig_PCSAFT**3)
!!                        gl%ptp(nrsubst) = P_CALC(gl, gl%ttp(nrsubst), gl%rhotp(nrsubst), nrsubst)
!!                        gl%rhomaxfluid(nrsubst) = gl%rhotp(nrsubst) * 1.5d0
!!                    else
!!                        gl%rhomaxfluid(nrsubst) = 4.d0 * gl%rhoc(nrsubst)
!!                        gl%rhotp(nrsubst) = rhomix_calc (gl, gl%ttp(nrsubst), gl%ptp(nrsubst), 0.d0, 1, nrsubst)
!!                        if (gl%rhotp(nrsubst) < gl%rhoc(nrsubst)) then
!!                            gl%rhotp(nrsubst) = 3.d0 * gl%rhoc(nrsubst)
!!                            gl%rhomaxfluid(nrsubst) = 4.d0 * gl%rhoc(nrsubst)
!!                        else
!!                            gl%rhomaxfluid(nrsubst) = 1.5d0 * gl%rhotp(nrsubst)
!!                        end if
!!                    end if
!!
!!                    !rhomaxfluid(nrsubst) = 3.9d0*rhoc(nrsubst)
!!                    exit
!!                end if
!!5454            continue
!!
!!            End do
!!
!!            close(pcsaft_unit)
!!
!!            if (dabs(gl%tc(nrsubst)) < 1.D-14) then
!!                errorflag = -8878
!!                return
!!            end if
!!
!!
!!        elseif (gl%Eq_type(nrsubst) == 3) then     !PR is used for that fluid. Stefan Feb 2014
!!
!!            if (gl%components(nrsubst) == '') then
!!                exit loopsubst
!!            else
!!
!!                !filename = trim (pathforsingles) // trim (components(nrsubst)) // '.fld'  ! check for standard case: fluid file (fld)
!!                !changed for compatibility with gfortran
!!                !filename = trim(path) // 'PR\FLUIDS_PR\PRFLUIDS.fld'  ! check for standard case: fluid file (fld)
!!                filename = trim(path) // 'pr/fluids_pr/prfluids.fld'  ! check for standard case: fluid file (fld)
!!                inquire (file = filename, exist = exists)
!!                if (.not.exists) then   ! no suitable ppf found
!!                    errorflag = -8878
!!                    return
!!                end if
!!
!!                open(newunit=pr_unit, file=trim(filename), status='old', action='read', iostat=error)
!!                if (error /= 0) then
!!                    errorflag = -7886
!!                    return
!!                end if
!!
!!
!!                !Skip the first line in the file
!!                read(pr_unit,*,end=900) dummy_PR
!!
!!                Do i = 1, 1000
!!                    read(pr_unit,*,end=900) dummy_PR
!!                    if (dummy_PR(1:4) == '@END') exit
!!                    backspace(pr_unit)
!!                    read(pr_unit,*) fluidnames_PR, flname_PR_alter, CASNR_PR(nrsubst), MW_PR, accen_PR, pc_PR, Tc_PR, ptp_PR, ttp_PR, A_cv0_PR,&
!!                        B_cv0_PR, C_cv0_PR, D_cv0_PR,E_cv0_PR,F_cv0_PR,G_cv0_PR
!!                    if ((trim(fluidnames_PR(1:30)) == trim(gl%components(nrsubst))) .OR. (trim(flname_PR_alter) == trim(gl%components(nrsubst)))  &
!!                        & .OR. (CASNR_PR(nrsubst) == gl%components(nrsubst))) then
!!                        gl%substfullname(nrsubst) = fluidnames_PR
!!                        gl%substshortname(nrsubst) = flname_PR_alter
!!                        gl%substcasnr(nrsubst) = CASNR_PR(nrsubst)
!!                        gl%tc(nrsubst) = tc_PR
!!                        gl%pc(nrsubst) = pc_PR
!!                        gl%tred(nrsubst) = gl%tc(nrsubst)
!!                        gl%accen(nrsubst) = accen_PR
!!                        !Andreas, November 2015: PURE PR equation needs to be evaluated here!!
!!                        !Set tredmix to tc(nrsubst)
!!                        gl%tredmix = gl%tc(nrsubst)
!!                        gl%rhored(nrsubst) = rho_PR (gl, gl%tc(nrsubst), gl%pc(nrsubst), 0, nrsubst)!tc(nrsubst), pc(nrsubst), 1,
!!                        !Set back tredmix
!!                        gl%tredmix = 1.D0
!!                        gl%rhoc(nrsubst) = gl%rhored(nrsubst)
!!                        gl%pmaxfluid(nrsubst) = 1000.d0    !Needed for PhaseDet_TD (mixtures) (Dummy value)
!!                        gl%tminfluid(nrsubst) = ttp_PR       !Needed for Flash calculations, arbitrary value!! Andreas, May 2016
!!                        !Changed from 1.d0 to ttp because for PH and PS Flash calculations ran into problems. Sebastian, October 2016
!!                        gl%tmaxfluid(nrsubst) = 2000.d0    !Needed for Flash calculations, arbitrary value!! Andreas, May 2016
!!                        gl%wm(nrsubst) = MW_PR/gl%factor
!!                        gl%REQ(nrsubst) = 8.3144598D0
!!                        gl%A_cv0(nrsubst)= A_cv0_PR
!!                        gl%B_cv0(nrsubst)= B_cv0_PR
!!                        gl%C_cv0(nrsubst)= C_cv0_PR
!!                        gl%D_cv0(nrsubst)= D_cv0_PR
!!                        gl%E_cv0(nrsubst)= E_cv0_PR
!!                        gl%F_cv0(nrsubst)= F_cv0_PR
!!                        gl%G_cv0(nrsubst)= G_cv0_PR
!!                        if (((dabs(gl%B_cv0(nrsubst)) > 1.D-8) .and. (dabs(gl%C_cv0(nrsubst)) > 1.D-8))) gl%cpmodel(nrsubst) = .true.
!!                        gl%refstate(nrsubst) = "OT0"
!!                        gl%tref(nrsubst) = 298.15d0
!!                        gl%pref(nrsubst) = 0.101325d0  !"-1" shows that pref must be calculated later in another subroutine
!!                        gl%href(nrsubst) = 0.d0 ![J/kg]     -later changed with wm into [J/mol]
!!                        gl%sref(nrsubst) = 0.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
!!                        gl%tcp0red(nrsubst) = gl%tc(nrsubst)
!!                        gl%cp0red(nrsubst) = 8.3144598D0
!!                        gl%bi_PR(nrsubst) = 0.0778D0*gl%R_PR * gl%tc(nrsubst) / (gl%pc(nrsubst)*1.D6)
!!                        gl%rhomaxfluid(nrsubst) = 1/gl%bi_PR(nrsubst)
!!                        !Triple point information !Andreas, July 2014
!!                        gl%ptp(nrsubst) = ptp_PR
!!                        gl%ttp(nrsubst) = ttp_PR
!!                        if (.not. gl%hold_limits) gl%ttp(nrsubst) = gl%ttp(nrsubst)*0.1d0    !small triple point temperature for the extrapolation of the VLE curves
!!                        !     (if "&" occurs in the input variable)
!!                        !Andreas, May 2018: PURE PR equation needs to be evaluated here!! (important for Peng-Robinson in Multi-Fluid mixture model)
!!                        !Set tredmix to tc(nrsubst)
!!                        gl%tredmix = gl%tc(nrsubst)
!!                        gl%rhoredmix = gl%rhoc(nrsubst)
!!                        if ((gl%ttp(nrsubst) > 1.d-14) .and. (gl%ptp(nrsubst) > 1.d-14)) then !Sebastian: Theresa neu funktioniert nicht wenn ttp und ptp = 0
!!                            gl%rhotp(nrsubst) = rho_PR (gl, gl%ttp(nrsubst), gl%ptp(nrsubst), 0, 1) !Theresa neu
!!                        end if
!!                        !Set back tredmix and rhoredmix
!!                        gl%tredmix = 1.D0
!!                        gl%rhoredmix = 1.D0
!!                        exit
!!                    end if
!!                End do
!!
!!            End if
!!
!!            close(pr_unit)
!!
!!            if (dabs(gl%tc(nrsubst)) < 1.D-14) then
!!                errorflag = -8878
!!                return
!!            end if
!!
!!        elseif (gl%Eq_type(nrsubst) == 4) then     !LKP is used for that fluid. Stefan May 2014
!!
!!            if (gl%components(nrsubst) == '') then
!!                exit loopsubst
!!            else
!!
!!                !filename = trim (pathforsingles) // trim (components(nrsubst)) // '.fld'  ! check for standard case: fluid file (fld)
!!                !changed for compatibility with gfortran
!!                !filename = trim(path) // 'LKP\FLUIDS_LKP\LKPFLUIDS.fld'  ! check for standard case: fluid file (fld)
!!                filename = trim(path) // 'lkp/fluids_lkp/lkpfluids.fld'  ! check for standard case: fluid file (fld)
!!                inquire (file = filename, exist = exists)
!!                if (.not.exists) then   ! no suitable ppf found
!!                    errorflag = -8878
!!                    return
!!                end if
!!
!!                open(newunit=lkp_unit, file=trim(filename), status='old', action='read', iostat=error)
!!                if (error /= 0) then
!!                    errorflag = -7887
!!                    return
!!                end if
!!
!!
!!                !Skip the first line in the file
!!                read(lkp_unit,*,end=900) dummy_LKP
!!
!!                Do i = 1, 1000
!!                    read(lkp_unit,*,end=900) dummy_LKP
!!                    if (dummy_LKP(1:4) == '@END') exit
!!                    backspace(lkp_unit)
!!                    read(lkp_unit,*) fluidnames_LKP, flname_LKP_alter, CASNR_LKP(nrsubst), MW_LKP, accen_LKP, pc_LKP, Tc_LKP, ptp_LKP, ttp_LKP, A_cv0_LKP,&
!!                        B_cv0_LKP, C_cv0_LKP, D_cv0_LKP,E_cv0_LKP,F_cv0_LKP,G_cv0_LKP
!!                    if ((trim(fluidnames_LKP(1:30)) == trim(gl%components(nrsubst))) .OR. (trim(flname_LKP_alter(1:30)) == trim(gl%components(nrsubst)))  &
!!                        & .OR. (CASNR_LKP(nrsubst) == gl%components(nrsubst))) then
!!                        gl%substfullname(nrsubst) = fluidnames_LKP
!!                        gl%substshortname(nrsubst) = flname_LKP_alter
!!                        gl%substcasnr(nrsubst) = CASNR_LKP(nrsubst)
!!                        gl%tc(nrsubst) = tc_LKP
!!                        gl%pc(nrsubst) = pc_LKP
!!                        gl%tred(nrsubst) = gl%tc(nrsubst)
!!                        gl%accen(nrsubst) = accen_LKP
!!                        gl%REQ(nrsubst) = 8.3144598D0
!!                        gl%rhored(nrsubst) = 1.d0/((0.2905d0 - 0.085d0 * gl%accen(nrsubst)) * gl%tc(nrsubst) * gl%REQ(nrsubst) / (gl%pc(nrsubst)*1.d6)) !rho_LKP (tc(nrsubst), pc(nrsubst), 0, 1)!tc(nrsubst), pc(nrsubst), 1,
!!                        gl%rhoc(nrsubst) = gl%rhored(nrsubst)
!!                        gl%vc_LKP(nrsubst) = 1.d0/gl%rhoc(nrsubst)
!!                        gl%pmaxfluid(nrsubst) = 1000.d0    !Needed for PhaseDet_TD (mixtures) (Dummy value)
!!                        gl%tminfluid(nrsubst) = ttp_LKP!1.D0       !Needed for Flash calculations, arbitrary value!! Andreas, May 2016
!!                        !Changed from 1.d0 to ttp because for PH and PS Flash calculations ran into problems. Sebastian, October 2016
!!                        gl%tmaxfluid(nrsubst) = 2000.D0    !Needed for Flash calculations, arbitrary value!! Andreas, May 2016
!!                        gl%wm(nrsubst) = MW_LKP
!!                        gl%A_cv0(nrsubst)= A_cv0_LKP
!!                        gl%B_cv0(nrsubst)= B_cv0_LKP
!!                        gl%C_cv0(nrsubst)= C_cv0_LKP
!!                        gl%D_cv0(nrsubst)= D_cv0_LKP
!!                        gl%E_cv0(nrsubst)= E_cv0_LKP
!!                        gl%F_cv0(nrsubst)= F_cv0_LKP
!!                        gl%G_cv0(nrsubst)= G_cv0_LKP
!!                        if (((dabs(gl%B_cv0(nrsubst)) > 1.D-8) .and. (dabs(gl%C_cv0(nrsubst)) > 1.D-8))) gl%cpmodel(nrsubst) = .true.
!!                        gl%refstate(nrsubst) = "OT0"
!!                        gl%tref(nrsubst) = 298.15d0
!!                        gl%pref(nrsubst) = 0.101325d0  !"-1" shows that pref must be calculated later in another subroutine
!!                        gl%href(nrsubst) = 0.d0 ![J/kg]     -later changed with wm into [J/mol]
!!                        gl%sref(nrsubst) = 0.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
!!                        gl%tcp0red(nrsubst) = gl%tc(nrsubst)
!!                        gl%cp0red(nrsubst) = 8.3144598D0
!!                        gl%rhomaxfluid(nrsubst) = 80000.d0
!!                        !Triple point information !Andreas, July 2014
!!                        gl%ptp(nrsubst) = ptp_LKP
!!                        gl%ttp(nrsubst) = ttp_LKP
!!                        if (.not. gl%hold_limits) gl%ttp(nrsubst) = gl%ttp(nrsubst)*0.1d0    !small triple point temperature for the extrapolation of the VLE curves
!!                        !     (if "&" occurs in the input variable)
!!                        if (gl%ncomp  == 1) Then
!!                            gl%zcLKP=(gl%pc(1)*1.d6)/(gl%tc(1)*gl%REQ(1)*gl%rhoc(1))
!!                        end if
!!
!!                        exit
!!                    end if
!!                End do
!!
!!            end if
!!
!!            close(lkp_unit)
!!
!!            !Fluid not found, return
!!            if (dabs(gl%tc(nrsubst)) < 1.D-14) then
!!                errorflag = -8878
!!                return
!!            end if
!!
!!            !change unit of molar mass from g/mol into kg/mol:
!!            gl%wm(nrsubst) = gl%wm(nrsubst) *factor_inv           ![kg/mol]
!!            wmcur(nrsubst) = wmcur(nrsubst) *factor_inv    ![kg/mol]
!!
!!
!!        elseif ((gl%Eq_type(nrsubst) == 51) .or. (gl%Eq_type(nrsubst) == 52) .or. (gl%Eq_type(nrsubst) == 53)) then     !Generalized EOS is used for that fluid; M. Thol, M. Kluge 04/2015
!!
!!            if (gl%components(nrsubst) == '') then
!!                exit loopsubst
!!            end if
!!
!!            filename = trim(path) // 'gen_eq\fluids_gen_eq\gen_eq_fluids.txt'
!!
!!            open(newunit=gen_unit, file=trim(filename), status='old', action='read', iostat=error)
!!            if (error /= 0) then
!!                errorflag = -9000
!!                return
!!            end if
!!
!!            !Skip the first line in the file
!!            read(gen_unit,*,end=900) dummy_GEOS
!!
!!            Do i = 1, 1000
!!                read(gen_unit,*,end=900) dummy_GEOS
!!                if (dummy_GEOS(1:4) == '@END') exit
!!                backspace(gen_unit)
!!                read(gen_unit,*) fluidnames_GEOS, flname_GEOS_alter, CASNR_GEOS(nrsubst), MW_GEOS, accen_GEOS, rhoc_GEOS, Tc_GEOS, ptp_GEOS, ttp_GEOS, polfac_eq_GEOS
!!                if ((trim(fluidnames_GEOS) == trim(gl%components(nrsubst))) .OR. (trim(flname_GEOS_alter) == trim(gl%components(nrsubst)))  &
!!                    & .OR. (CASNR_GEOS(nrsubst) == gl%components(nrsubst))) then
!!                    gl%tc(nrsubst) = tc_GEOS
!!                    gl%rhoc(nrsubst) = rhoc_GEOS*gl%factor   !change unit to mol/m?
!!                    gl%rhored(nrsubst) = gl%rhoc(nrsubst)
!!                    gl%tred(nrsubst) = gl%tc(nrsubst)
!!                    gl%accen(nrsubst) = accen_GEOS
!!                    gl%pmaxfluid(nrsubst) = 1.d12    !Needed for PhaseDet_TD (mixtures) (Dummy value)
!!                    gl%wm(nrsubst) = MW_GEOS
!!                    gl%REQ(nrsubst) = 8.31451D0
!!                    if (CASNR_GEOS(nrsubst) == '999-999-999') gl%REQ(nrsubst) = 1.d0
!!                    !A_cv0(nrsubst)= A_cv0_SRK
!!                    !B_cv0(nrsubst)= B_cv0_SRK
!!                    !C_cv0(nrsubst)= C_cv0_SRK
!!                    !D_cv0(nrsubst)= D_cv0_SRK
!!                    !E_cv0(nrsubst)= E_cv0_SRK
!!                    !F_cv0(nrsubst)= F_cv0_SRK
!!                    !G_cv0(nrsubst)= G_cv0_SRK
!!                    gl%refstate(nrsubst) = "OT0"
!!                    gl%tref(nrsubst) = 298.15d0
!!                    gl%pref(nrsubst) = 0.101325d0
!!                    gl%href(nrsubst) = 0.d0 ![J/kg]     -later changed with wm into [J/mol]
!!                    gl%sref(nrsubst) = 0.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
!!                    gl%tcp0red(nrsubst) = gl%tc(nrsubst)
!!                    gl%cp0red(nrsubst) = 8.31451D0
!!                    gl%rhomaxfluid(nrsubst) = 6.5D0*gl%rhoc(nrsubst)  !dummy value
!!                    !Triple point information
!!                    gl%ptp(nrsubst) = ptp_GEOS
!!                    gl%ttp(nrsubst) = ttp_GEOS
!!                    gl%rhotp(nrsubst) = 1.d0 !rho_SRK(ttp(nrsubst), ptp(nrsubst), 0, 1)  Theresa: there is no pc defined so far!
!!                    gl%polfac_eq = polfac_eq_GEOS
!!
!!2209                format (D14.7,F12.3,F6.1,F6.0)
!!
!!
!!                    if (gl%Eq_type(nrsubst) == 51) then    !calculation based on Alexandrov
!!                        if (gl%accen(nrsubst) < 0.d0) then
!!                            errorflag = -6001
!!                            return
!!                        end if
!!
!!
!!
!!                        call initialize_Gen_Eq_Igor(gl)
!!                        gl%eos_coeff%nreg(nrsubst) = gl%anz_term_igor
!!                        do gen_i = 1, gl%anz_term_igor
!!                            gl%eos_coeff%ni(gen_i, nrsubst) = gl%coefficientsi(gen_i, 1) + (gl%coefficientsi(gen_i, 2)*gl%accen(nrsubst)) + (gl%coefficientsi(gen_i, 3)*gl%accen(nrsubst)**(gl%coefficientsi(gen_i, 4)))
!!                            gl%eos_coeff%di(gen_i, nrsubst) = gl%exponentsi(gen_i, 1)
!!                            gl%eos_coeff%ti(gen_i, nrsubst) = gl%exponentsi(gen_i, 2)
!!                            gl%eos_coeff%p_i(gen_i, nrsubst) = gl%exponentsi(gen_i, 3)
!!                            if (gl%eos_coeff%p_i(gen_i, nrsubst) /= 0) gl%eos_coeff%gama(gen_i, nrsubst) = 1.d0
!!                            !write (71,2209) gl%eos_coeff%ni(gen_i, nrsubst), gl%eos_coeff%ti(gen_i, nrsubst), gl%eos_coeff%di(gen_i, nrsubst), gl%eos_coeff%p_i(gen_i, nrsubst)
!!                        end do
!!
!!                    elseif (gl%Eq_type(nrsubst) == 52) then  !calculation based on Span/Wagner
!!
!!                        call initialize_Gen_Eq_Span(gl)
!!                        gl%eos_coeff%nreg(nrsubst) = gl%anz_term_span
!!
!!                        do gen_i = 1, gl%anz_term_span
!!                            gl%eos_coeff%ni(gen_i, nrsubst) = gl%coefficientss(gen_i, 1) + (gl%coefficientss(gen_i, 2)*gl%accen(nrsubst)) + (gl%coefficientss(gen_i, 3)*(gl%accen(nrsubst)**4))
!!                            gl%eos_coeff%ti(gen_i, nrsubst) = gl%exponentss(gen_i, 1)
!!                            gl%eos_coeff%di(gen_i, nrsubst) = gl%exponentss(gen_i, 2)
!!                            gl%eos_coeff%gama(gen_i, nrsubst) = gl%exponentss(gen_i, 3)
!!                            gl%eos_coeff%p_i(gen_i, nrsubst) = gl%exponentss(gen_i, 4)
!!                            !write (71,2209) gl%eos_coeff%ni(gen_i, nrsubst), gl%eos_coeff%ti(gen_i, nrsubst), gl%eos_coeff%di(gen_i, nrsubst), gl%eos_coeff%p_i(gen_i, nrsubst)
!!                        end do
!!
!!                    elseif (gl%Eq_type(nrsubst) == 53) then  !calculation based on Sun/Ely
!!
!!                        call initialize_Gen_Eq_Sun(gl)
!!                        gl%eos_coeff%nreg(nrsubst) = gl%anz_term_sun
!!
!!                        do gen_i = 1, gl%anz_term_sun
!!                            a_sun1(gen_i, nrsubst) = (gl%coefficientssun(gen_i,2)-gl%coefficientssun(gen_i,1))/(gl%acenfactor_2-gl%acenfactor_1)
!!                            a_sun2(gen_i, nrsubst) = 1/(gl%polarfactor_3-gl%polarfactor_1)*((gl%coefficientssun(gen_i,3)-gl%coefficientssun(gen_i,1))-((gl%acenfactor_3-gl%acenfactor_1)/(gl%acenfactor_2-gl%acenfactor_1)*(gl%coefficientssun(gen_i,2)-gl%coefficientssun(gen_i,1))))
!!
!!                            gl%eos_coeff%ni(gen_i, nrsubst) = gl%coefficientssun(gen_i,1)+((gl%accen(nrsubst)-gl%acenfactor_1)*a_sun1(gen_i, nrsubst))+((gl%polfac_eq-gl%polarfactor_1)*a_sun2(gen_i, nrsubst))
!!                            gl%eos_coeff%ti(gen_i, nrsubst) = gl%exponentssun(gen_i,2)
!!                            gl%eos_coeff%di(gen_i, nrsubst) = gl%exponentssun(gen_i,1)
!!                            gl%eos_coeff%p_i(gen_i, nrsubst) = gl%exponentssun(gen_i,3)
!!                            if (gl%eos_coeff%p_i(gen_i, nrsubst) /= 0) gl%eos_coeff%gama(gen_i, nrsubst) = 1.d0
!!                        end do
!!
!!                    end if
!!
!!                    tc_est_inp = gl%tc(nrsubst)
!!                    rhoc_est_inp = gl%rhoc(nrsubst)
!!                    call find_crit_tpd (gl,tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errval, nrsubst)
!!                    if (errval == 0) then
!!                        gl%tc(nrsubst) = tc_eos_
!!                        gl%rhoc(nrsubst) = rhoc_eos_
!!                        gl%pc(nrsubst) = pc_eos_
!!                    end if
!!
!!
!!                    !if (ncomp > 1) call reduced_parameters_calc(300.D0) !Dummy temperature 300 K for the SRK
!!                    if (gl%VLE_needed) then
!!                        ncomp_orig = gl%ncomp
!!                        !gl%ncomp = 1
!!                        gl%rhoredmix = gl%rhored(nrsubst)
!!                        gl%Tredmix = gl%Tred(nrsubst)
!!                        gl%pc(nrsubst) = p_calc(gl,gl%tc(nrsubst),gl%rhoc(nrsubst),nrsubst)
!!                        gl%rhotp(nrsubst) = rho_SRK(gl, gl%ttp(nrsubst), gl%ptp(nrsubst), 0, 1)
!!                        !gl%ncomp = ncomp_orig
!!                    end if
!!                    exit
!!                end if
!!            End do
!!            !Ideal part oil
!!            if (gl%components(nrsubst) == 'oil') then
!!                gl%cpmodel(nrsubst) = .true.
!!                ! setting up common CPP module variables
!!                cpdummy=gl%eos_coeff%cppcheck(nrsubst)
!!                cpdummy(1:2) = 'CP'
!!                gl%tcp0red(nrsubst) = 1.d0
!!                gl%cp0red(nrsubst) = 8.3144621d0
!!                gl%ncp0poly(nrsubst) = 1
!!                gl%CP0COEFF(1,nrsubst) = 0.414793476829D+03
!!                gl%cp0exp(1,nrsubst) = 0.d0
!!                gl%ncp0pl(nrsubst) = 1
!!                gl%CP0COEFF(2,nrsubst) = 0.878393464574D+03
!!                gl%cp0exp(2,nrsubst) = 2281.40745904d0
!!            end if
!!
!!            close(gen_unit)
!!
!!            !Fluid not found, return
!!            if (dabs(gl%tc(nrsubst)) < 1.D-14) then
!!                errorflag = -8878
!!                return
!!            end if
!!
!!            !change unit of molar mass from g/mol into kg/mol:
!!            gl%wm(nrsubst) = gl%wm(nrsubst) *factor_inv           ![kg/mol]
!!            wmcur(nrsubst) = wmcur(nrsubst) *factor_inv    ![kg/mol]
!!
!!        end if
!!
!!    enddo loopsubst !read files end
!!
!!
!!    !for mixtures check if psrk calculation is possible:
!!    if (gl%ncomp == 1) then !no mixture, no psrk...
!!        gl%psrk_ok = .false.
!!    else !mixture, psrk possible...
!!        do j=1,gl%ngroup_nsubgroup
!!            if (int(Qkcoeffs(j)) == 0) then
!!                gl%psrk_ok = .false.
!!            else
!!                gl%psrk_ok = .true.
!!            end if
!!        enddo
!!        if (gl%psrk_ok) then
!!            !changed for compatibility with gfortran
!!            !call read_atcoeff_file(pathforsingles // '\UFC\',error)
!!            !call read_atcoeff_file(gl,pathforsingles // '/gl%ufc/',error)
!!            gl%psrk_ok = .false.
!!            !if (error /= 0) then
!!            !    !error = error + 4000 !high numbers indicate errors in inner subroutines
!!            !    error = 0
!!            !    gl%psrk_ok = .false.
!!            !else
!!            !    !check if all needed at-coefficients are given:
!!            !    atloop:         do j=1,gl%ngroup_nsubgroup-1
!!            !        if ((int(gl%atcoeff(int(gl%ufc(j,1)),int(gl%ufc(j+1,1)))) == 0) &
!!            !            .and.(int(gl%atcoeff(int(gl%ufc(j+1,1)),int(gl%ufc(j,1)))) == 0)) then
!!            !            gl%psrk_ok = .false. !no coefficient given, but is needed
!!            !            exit atloop
!!            !        endif
!!            !    enddo atloop
!!            !endif
!!            !possibility to read btcoeff file and ctcoeff file here, if psrk_ok
!!        endif
!!    endif
!!
!!    !This call does not make sense since the mixing parameters have not been read in yet --> Commented this line
!!    !Andreas J輍er March 2017
!!    !for mixtures recalculate the reduced parameters
!!    !call reduced_parameters_calc(300.D0) !Dummy temperature 300 K for the SRK
!!
!!    !pressure at triple point temperature calculated from vapor pressure equation
!!    ptrloop: do nrsubst = 1, gl%ncomp
!!
!!        !if (pvexist(nrsubst)) then
!!        !    ptp(nrsubst) = vp_eq(ttp(nrsubst), nrsubst)
!!        !end if
!!        !    ptpmod(nrsubst) = ptp(nrsubst)*0.1d0
!!
!!        ! Andreas October 2012
!!        ! If no information about the triple point temperature is given, set it to an arbitrary value, since it is needed for the phase equilibrium routines
!!        ! The same is true for the triple point pressure
!!        ! This is especially important when using different equations than Helmholtz equations!!!
!!
!!        !if (ttp(nrsubst) == 0.D0) then
!!        !    ttp(nrsubst) = 50.D0   !Arbitrary Temp of 50 K
!!        !end if
!!        !if (ptp(nrsubst) == 0.D0) then
!!        !    ptp(nrsubst) = 0.01D0   !Abitrary press of 0,01 MPa
!!        !end if
!!        ! Changed the estimation of the triple point pressure and temperature by more reasonable estimates depending on the critical values
!!        ! Drastically wrong estimates can cause the phase envelope to fail. (Example: 50% Hydrogen, 50% Methane)
!!        ! Andreas Jul 2013
!!        if (dabs(gl%ttp(nrsubst)) < 1.D-14) then
!!            !ttp(nrsubst) = 150.D0   !Arbitrary Temp of 150 K
!!            gl%ttp(nrsubst) = gl%tc(nrsubst) * 0.5d0
!!            if (gl%components(nrsubst) == "oil") gl%ttp(nrsubst) = 200.d0
!!        end if
!!        if ((dabs(gl%ptp(nrsubst)) < 1.D-14) .and. (gl%components(nrsubst) /= "oil"))then
!!            gl%ptp(nrsubst) = 0.01D0   !Abitrary press of 0,01 MPa
!!            !ptp(nrsubst) = pc(nrsubst) / 2.D0
!!        end if
!!
!!    enddo ptrloop
!!
!!
!!    !Read mixing rules
!!    !Andreas, Stefan / June 2014
!!    If (gl%ncomp > 1) then !read mixingcoefficients if more than one component:
!!
!!        !allocation of variables used in mixtures_i
!!        if(.not.allocated(gl%dfn)) then
!!            allocate(gl%dfn(25,gl%ncomp,gl%ncomp))
!!            allocate(gl%dfd,gl%dft,gl%dfl,gl%dfp,gl%dfeta,gl%dfeps,gl%dfbeta,gl%dfgamma,gl%dfgeta,gl%dfgbeta,gl%dfggam,gl%dfgeps,mold=gl%dfn)
!!            gl%dfn     = 0d0
!!            gl%dfd     = 0d0
!!            gl%dft     = 0d0
!!            gl%dfl     = 0d0
!!            gl%dfp     = 0d0
!!            gl%dfeta   = 0d0
!!            gl%dfeps   = 0d0
!!            gl%dfbeta  = 0d0
!!            gl%dfgamma = 0d0
!!            gl%dfgeta  = 0d0
!!            gl%dfgbeta = 0d0
!!            gl%dfggam  = 0d0
!!            gl%dfgeps  = 0d0
!!
!!
!!            allocate(gl%dfd_coeff_structure(gl%ncomp,gl%ncomp))
!!            allocate(gl%dfl_coeff_structure(gl%ncomp,gl%ncomp))
!!            gl%dfd_coeff_structure = coeff_structure_all_zeros
!!            gl%dfl_coeff_structure = coeff_structure_all_zeros
!!        end if
!!
!!        If ((gl%Mix_type == 1) .or. (gl%Mix_type == 7) .or. (gl%Mix_type == 110) .or. (gl%Mix_type == 111) .or. (gl%Mix_type == 120) .or. (gl%Mix_type == 121)) then
!!            do nrsubst = 1, gl%ncomp
!!                if (gl%components(nrsubst+1) /= '') then !'' means the current substance is also the last one
!!                    do j = nrsubst+1,gl%ncomp !always compare current substance to next ones
!!                        !if (Eq_type(nrsubst) == 1) then         !Do this only if for both fluids Helmholtz EoS are used! Andreas Okt. 2012. Commented, because Helmholtz and other models may be combinated now, Andreas J輍er, March 2017.
!!                        !110, 120: special cases => all components are mixed according to linear or LB mixing rule
!!                        if (gl%Mix_type == 110) then
!!                            !All components are mixed according to the linear mixing rule
!!                            gl%rfbetat(nrsubst,j) = 1.d0
!!                            gl%rfbetarho(nrsubst,j) = 1.d0
!!                            gl%rfbetat(j,nrsubst) = 1.d0
!!                            gl%rfbetarho(j,nrsubst) = 1.d0
!!
!!                            gl%rfgammat(nrsubst,j) = 0.5d0*(gl%tc(nrsubst)+gl%tc(j))/dsqrt(gl%tc(nrsubst)*gl%tc(j))
!!                            gl%rfgammat(j,nrsubst) = gl%rfgammat(nrsubst,j)
!!                            gl%rfgammarho(nrsubst,j) = 4.d0*(1.d0/gl%rhoc(nrsubst)+1.d0/gl%rhoc(j))/(gl%rhoc(nrsubst)**(-1.d0/3.d0)+gl%rhoc(j)**(-1.d0/3.d0))**3
!!                            gl%rfgammarho(j,nrsubst) = gl%rfgammarho(nrsubst,j)
!!                            if(nrsubst == gl%ncomp) gl%mix_type = 1
!!                            !---------------------------------------------------------------------------------------------------------------
!!                        elseif (gl%Mix_type == 120) then
!!                            !All components are mixed according to the Lorentz-Berthelot mixing rule.
!!                            gl%rfbetat(nrsubst,j) = 1.d0
!!                            gl%rfbetarho(nrsubst,j) = 1.d0
!!                            gl%rfgammat(nrsubst,j) = 1.d0
!!                            gl%rfgammarho(nrsubst,j) = 1.d0
!!                            gl%rfbetat(j,nrsubst) = 1.d0
!!                            gl%rfbetarho(j,nrsubst) = 1.d0
!!                            gl%rfgammat(j,nrsubst) = 1.d0
!!                            gl%rfgammarho(j,nrsubst) = 1.d0
!!                            if(nrsubst == gl%ncomp) gl%mix_type = 1
!!                        else
!!                            !check if binary mix files are available
!!                            call read_mix_file(gl,nrsubst,j,pathforbinarys, error)
!!                            if (error /= 0) then
!!                                !Added error handling if binary mix file was not read correctly, Andreas August 2013
!!                                !error = error + 2000 !high numbers show errors in inner subroutines
!!
!!                                if (gl%Mix_type == 111) then
!!                                    gl%rfbetat(nrsubst,j) = 1.d0
!!                                    gl%rfbetarho(nrsubst,j) = 1.d0
!!                                    gl%rfbetat(j,nrsubst) = 1.d0
!!                                    gl%rfbetarho(j,nrsubst) = 1.d0
!!
!!                                    gl%rfgammat(nrsubst,j) = 0.5d0*(gl%tc(nrsubst)+gl%tc(j))/dsqrt(gl%tc(nrsubst)*gl%tc(j))
!!                                    gl%rfgammat(j,nrsubst) = gl%rfgammat(nrsubst,j)
!!                                    gl%rfgammarho(nrsubst,j) = 4.d0*(1.d0/gl%rhoc(nrsubst)+1.d0/gl%rhoc(j))/(gl%rhoc(nrsubst)**(-1.d0/3.d0)+gl%rhoc(j)**(-1.d0/3.d0))**3
!!                                    gl%rfgammarho(j,nrsubst) = gl%rfgammarho(nrsubst,j)
!!                                    error = 0
!!                                    if(nrsubst == gl%ncomp) gl%mix_type = 1
!!                                    !---------------------------------------------------------------------------------------------------------------
!!                                elseif (gl%Mix_type == 121) then
!!                                    gl%rfbetat(nrsubst,j) = 1.d0
!!                                    gl%rfbetarho(nrsubst,j) = 1.d0
!!                                    gl%rfgammat(nrsubst,j) = 1.d0
!!                                    gl%rfgammarho(nrsubst,j) = 1.d0
!!                                    gl%rfbetat(j,nrsubst) = 1.d0
!!                                    gl%rfbetarho(j,nrsubst) = 1.d0
!!                                    gl%rfgammat(j,nrsubst) = 1.d0
!!                                    gl%rfgammarho(j,nrsubst) = 1.d0
!!                                    error = 0
!!                                    if(nrsubst == gl%ncomp) gl%mix_type = 1
!!                                else
!!                                    !exit
!!                                    errorflag = error
!!                                    return
!!                                end if
!!                            endif
!!                        end if
!!                        !end if
!!                    enddo
!!                endif
!!            end do
!!        end if
!!
!!        !if mix_types 110, 111, 120, 121 are used, the same mixing rules are used as in the case of mix_type = 1.
!!        !In order to avoid adding the four special mix_types in every single routine, it is changed to mix_type = 1 here.
!!        if ((gl%mix_type .eq. 110) .or. (gl%mix_type .eq. 111) .or. (gl%mix_type .eq. 120) .or. (gl%mix_type .eq. 121)) gl%mix_type = 1
!!
!!        !!Andreas, November 2015. If Multiparameter EOS are mixed with other EOS, read Helmholtz mixing rules here;
!!        !!Monika, June 2016: 110 - 121 special cases if binary mix files not available
!!        !If (Mix_type == 1) then
!!        !    !read mixingcoefficients if more than one component:
!!        !    do i = 1, ncomp
!!        !        do j=i+i, ncomp
!!        !            call read_mix_file(i,j,pathforbinarys, error)
!!        !            if (error /= 0) then
!!        !                !Added error handling if binary mix file was not read correctly, Andreas August 2013
!!        !                !error = error + 2000 !high numbers show errors in inner subroutines
!!        !                !exit
!!        !                errorflag = error
!!        !                return
!!        !            endif
!!        !        end do
!!        !    end do
!!        !end if
!!        !---------------------------------------------------------------------------------------------------------------
!!
!!        !---------------------------------------------------------------------------------------------------------------
!!
!!
!!        !Read mixing rules for the SRK
!!        If ((gl%Mix_type == 2) .or. (gl%Mix_type == 21)) then
!!            !DEC$ IF DEFINED(HCgen)
!!            !DEC$ ELSE
!!            if(.not.allocated(gl%kij_SRK)) then
!!                allocate(gl%kij_SRK(gl%ncomp,gl%ncomp))
!!                allocate(gl%aij_SRK,gl%kij_PR,gl%aij_PR,gl%lij_SRK,gl%bij_SRK,gl%lij_PR,gl%bij_PR,mold=gl%kij_SRK)
!!            end if
!!            !DEC$ ENDIF
!!
!!            gl%kij_SRK = 0
!!            gl%lij_SRK = 0
!!            !Andreas, November 2015: Why is REQ set to this value here?? R for SRK is defined differently...
!!            !REQ = 8.31451D0
!!
!!            !Sebastian, Oct 2016 Killed the line below
!!            !REQ = R_SRK
!!
!!            !changed for compatibility with gfortran
!!            !filename = trim(path) // 'SRK\BINARY_SRK\BIN_SRK.MIX'  ! check for binary mixing parameters file for the SRK
!!            filename = trim(path) // 'srk/binary_srk/bin_srk.mix'  ! check for binary mixing parameters file for the SRK
!!            inquire (file = filename, exist = exists)
!!            if (.not.exists) then
!!                errorflag = -8878
!!                return
!!            end if
!!            open(newunit=srkmix_unit, file=trim(filename), status='old', action='read', iostat=error)
!!            if (error /= 0) then
!!                errorflag = -7888
!!                return
!!            end if
!!
!!
!!            !Skip the first three line in the file
!!            read(srkmix_unit,*,end=900) dummy_SRK
!!            read(srkmix_unit,*,end=900) dummy_SRK
!!            read(srkmix_unit,*,end=900) dummy_SRK
!!
!!            Do i = 1, 1000
!!                read(srkmix_unit,*,end=900) dummy_SRK
!!                if (dummy_SRK(1:4) == '@END') exit
!!                backspace(srkmix_unit)
!!                !!DEC$ IF DEFINED(MESSER)
!!                !!****************************************************************************************
!!                !!MESSER IDENTIFICATION OF MIXING RULES
!!                !read(srkmix_unit,*,iostat = error)  kompnr1, kompnr2, kij_read
!!                !if (error == 0) then
!!                !    do nrsubst = 1, gl%ncomp
!!                !        do j = nrsubst + 1,  gl%ncomp
!!                !            if ((kompnr_srk(j) == kompnr1) .and. (kompnr_srk(nrsubst) == kompnr2)) then
!!                !                gl%kij_srk(j,nrsubst) = kij_read
!!                !                gl%kij_srk(nrsubst,j) = kij_read
!!                !            end if
!!                !            if ((kompnr_srk(j) == kompnr2) .and. (kompnr_srk(nrsubst) == kompnr1)) then
!!                !                gl%kij_srk(j,nrsubst) = kij_read
!!                !                gl%kij_srk(nrsubst,j) = kij_read
!!                !            end if
!!                !        end do
!!                !    end do
!!                !end if
!!                !!DEC$ ELSE
!!                !****************************************************************************************
!!                !New identification of mixing rules using the CAS numbers of the fluids
!!                !Adjusted read statement because the covolume parameter "b" can now also be mixed with a quadratic mixing rule. For this, the parameter lij was introduced in the file bin_srk.mix
!!                !Andreas J輍er, January 2017
!!                !Try new reading routine first
!!                read(srkmix_unit,*,iostat = error)  CASNR_1, CASNR_2, kij_read, lij_read
!!                !Use old routine if error occurs. Safety net, can maybe be deleted in the future.
!!                if (error /= 0) then
!!                    backspace(srkmix_unit)
!!                    read(srkmix_unit,*,iostat = error)  CASNR_1, CASNR_2, kij_read
!!                    !If no 4th parameter is given, set lij to 0
!!                    lij_read = 0.D0
!!                    if (error/=0) then
!!                        errorflag = -7895
!!                    end if
!!                end if
!!
!!
!!
!!
!!                do nrsubst = 1, gl%ncomp
!!                    Do j = nrsubst + 1, gl%ncomp
!!                        if ((CASNR_SRK(j) == CASNR_1) .and. (CASNR_SRK(nrsubst) == CASNR_2)) then
!!                            gl%kij_SRK(j,nrsubst) = kij_read
!!                            gl%kij_SRK(nrsubst,j) = kij_read
!!                            gl%lij_SRK(j,nrsubst) = lij_read
!!                            gl%lij_SRK(nrsubst,j) = lij_read
!!                            backspace(srkmix_unit)
!!                            read(srkmix_unit,'(A)')dummy
!!                            k = 1
!!                            do while (dummy(k:k) /= '!')
!!                                k = k + 1
!!                            enddo
!!                            if (allocated(gl%litref)) then
!!                                gl%litref%lit_ref_mix = dummy(k+1:)
!!                            end if
!!                        end if
!!                        if ((CASNR_SRK(j) == CASNR_2) .and. (CASNR_SRK(nrsubst) == CASNR_1)) then
!!                            gl%kij_SRK(j,nrsubst) = kij_read
!!                            gl%kij_SRK(nrsubst,j) = kij_read
!!                            gl%lij_SRK(j,nrsubst) = lij_read
!!                            gl%lij_SRK(nrsubst,j) = lij_read
!!                            backspace(srkmix_unit)
!!                            read(srkmix_unit,'(A)')dummy
!!                            k = 1
!!                            do while (dummy(k:k) /= '!')
!!                                k = k + 1
!!                            enddo
!!                            if (allocated(gl%litref)) then
!!                                gl%litref%lit_ref_mix = dummy(k+1:)
!!                            end if
!!                        end if
!!                    End do
!!                end do
!!                !!DEC$ END IF !END DEFINE MESSER
!!            End do
!!
!!            close(srkmix_unit)
!!
!!
!!            !Read mixing rules for PC-Saft
!!        else If (gl%Mix_type == 6) then
!!
!!
!!            !changed for compatibility with gfortran
!!            !filename = trim(path) // 'SRK\BINARY_SRK\BIN_SRK.MIX'  ! check for binary mixing parameters file for the PC-Saft
!!            filename = trim(path) // 'pc-saft/binary_saft/bin_saft.mix'  ! check for binary mixing parameters file for PC-Saft
!!            inquire (file = filename, exist = exists)
!!            if (.not.exists) then
!!                errorflag = -8878
!!                return
!!            end if
!!            open(newunit=pcsaftmix_unit, file=trim(filename), status='old', action='read', iostat=error)
!!            if (error /= 0) then
!!                errorflag = -7888
!!                return
!!            end if
!!
!!
!!            !Skip the first three line in the file
!!            read(pcsaftmix_unit,*,end=900) dummy_PCSAFT
!!            read(pcsaftmix_unit,*,end=900) dummy_PCSAFT
!!            read(pcsaftmix_unit,*,end=900) dummy_PCSAFT
!!
!!            Do i = 1, 1000
!!                read(pcsaftmix_unit,*,end=900) dummy_PCSAFT
!!                if (dummy_PCSAFT(1:4) == '@END') exit
!!                backspace(pcsaftmix_unit)
!!                !****************************************************************************************
!!                !MESSER IDENTIFICATION OF MIXING RULES
!!                !read(10,*)  Kompnr1, Kompnr2, kij_read
!!                !do nrsubst = 1, ncomp
!!                !    Do j = nrsubst + 1, ncomp
!!                !            if ((Kompnr_SRK(j) == Kompnr1) .and. (Kompnr_SRK(nrsubst) == Kompnr2)) then
!!                !                kij_SRK(j,nrsubst) = kij_read
!!                !                kij_SRK(nrsubst,j) = kij_read
!!                !            end if
!!                !            if ((Kompnr_SRK(j) == Kompnr2) .and. (Kompnr_SRK(nrsubst) == Kompnr1)) then
!!                !                kij_SRK(j,nrsubst) = kij_read
!!                !                kij_SRK(nrsubst,j) = kij_read
!!                !            end if
!!                !    End do
!!                !end do
!!                !****************************************************************************************
!!                !New identification of mixing rules using the CAS numbers of the fluids
!!                !Adjusted read statement because the covolume parameter "b" can now also be mixed with a quadratic mixing rule. For this, the parameter lij was introduced in the file bin_srk.mix
!!                !Andreas J輍er, January 2017
!!                !Try new reading routine first
!!                read(pcsaftmix_unit,*,iostat = error)  CASNR_1, CASNR_2, kij_read
!!                do nrsubst = 1, gl%ncomp
!!                    Do j = nrsubst + 1, gl%ncomp
!!                        if ((CASNR_PCSAFT(j) == CASNR_1) .and. (CASNR_PCSAFT(nrsubst) == CASNR_2)) then
!!                            gl%kij_PCSAFT(j,nrsubst) = kij_read
!!                            gl%kij_PCSAFT(nrsubst,j) = kij_read
!!                            backspace(pcsaftmix_unit)
!!                            read(pcsaftmix_unit,'(A)')dummy
!!                            k = 1
!!                            do while (dummy(k:k) /= '!')
!!                                k = k + 1
!!                            enddo
!!                            if (allocated(gl%litref)) then
!!                                gl%litref%lit_ref_mix = dummy(k+1:)
!!                            end if
!!                        end if
!!                        if ((CASNR_PCSAFT(j) == CASNR_2) .and. (CASNR_PCSAFT(nrsubst) == CASNR_1)) then
!!                            gl%kij_PCSAFT(j,nrsubst) = kij_read
!!                            gl%kij_PCSAFT(nrsubst,j) = kij_read
!!                            backspace(pcsaftmix_unit)
!!                            read(pcsaftmix_unit,'(A)')dummy
!!                            k = 1
!!                            do while (dummy(k:k) /= '!')
!!                                k = k + 1
!!                            enddo
!!                            if (allocated(gl%litref)) then
!!                                gl%litref%lit_ref_mix = dummy(k+1:)
!!                            end if
!!                        end if
!!                    End do
!!                end do
!!            End do
!!
!!            close(pcsaftmix_unit)
!!
!!            !PSRK mixing rules for parameter "a". For PSRK the UNIFAC parameters are needed which are stored in the files "Interac.par" and "rk_qk.par"
!!            !Andreas J輍er, January 2017
!!            !Use same routine for reading Helmholtz+gE parameters (Mix_type == 12). Andreas J輍er, August 2017
!!        Elseif ((gl%Mix_Type == 22) .or. (gl%Mix_type == 12)) then
!!
!!            !For Helm+gE model, read subgroups of fluids from file
!!            if (gl%mix_type == 12) then
!!
!!                do nrsubst = 1, gl%ncomp   !Loop over all components in the mixture
!!
!!                    !Check if groups file exist at given path
!!                    filename = trim(path) // 'ge_models/unifac/helm_ge/groups.par'
!!                    inquire (file = filename, exist = exists)
!!                    if (.not.exists) then
!!                        errorflag = -8879
!!                        return
!!                    end if
!!                    !Try to open groups file for reading
!!                    open(newunit=ge_unit, file=trim(filename), status='old', action='read', iostat=error)
!!                    if (error /= 0) then
!!                        errorflag = -8879
!!                        return
!!                    end if
!!
!!                    !Skip the first line in the file
!!                    read(ge_unit,*,iostat= error) dummy_SRK
!!                    if (error /= 0) then
!!                        errorflag = -7898
!!                        close(ge_unit)
!!                        return
!!                    end if
!!
!!                    Do i = 1, 1000
!!
!!                        !Read next line of file
!!                        read(ge_unit,*,iostat= error) fluidnames_SRK, flname_SRK_alter, CASNR_SRK(nrsubst), &
!!                            & nr_of_groups_i_read, subgroups_ik_list_read(1), subgroups_ik_list_read(2), &
!!                            & subgroups_ik_list_read(3), subgroups_ik_list_read(4), subgroups_ik_list_read(5), subgroups_ik_list_read(6), &
!!                            & subgroups_ik_list_read(7), subgroups_ik_list_read(8), subgroups_ik_list_read(9), subgroups_ik_list_read(10), &
!!                            & subgroups_ik_list_read(11), subgroups_ik_list_read(12), subgroups_ik_list_read(13), subgroups_ik_list_read(14), &
!!                            & subgroups_ik_list_read(15), subgroups_ik_list_read(16), subgroups_ik_list_read(17), subgroups_ik_list_read(18), &
!!                            & subgroups_ik_list_read(19), subgroups_ik_list_read(20), subgroups_ik_list_read(21), subgroups_ik_list_read(22), &
!!                            & subgroups_ik_list_read(23), subgroups_ik_list_read(24), subgroups_ik_list_read(25), subgroups_ik_list_read(26), &
!!                            & subgroups_ik_list_read(27), subgroups_ik_list_read(28), subgroups_ik_list_read(29), subgroups_ik_list_read(30)
!!
!!                        if (error /= 0) then
!!                            errorflag = -7898
!!                            close(ge_unit)
!!                            return
!!                        end if
!!
!!                        if ((trim(fluidnames_SRK(1:30)) == trim(gl%components(nrsubst))) .OR. (trim(flname_SRK_alter(1:30)) == trim(gl%components(nrsubst)))  &
!!                            & .OR. (CASNR_SRK(nrsubst) == gl%components(nrsubst))) then
!!                            !Save group information for substance found
!!                            gl%nr_of_groups_i(nrsubst) = nr_of_groups_i_read
!!                            gl%subgroups_ik_list(nrsubst,:) = subgroups_ik_list_read
!!                            gl%nr_of_groups_mix = gl%nr_of_groups_mix + gl%nr_of_groups_i(nrsubst)
!!                            close(ge_unit)
!!                            exit
!!                        end if
!!                    end do
!!
!!                end do
!!
!!
!!
!!
!!                !Andreas J輍er, April 2018. Use linear mixing rules as default
!!                gl%LIN_or_LB = 2
!!                if (gl%LIN_or_LB == 1) then
!!                    !LB-mixing rules
!!                    gl%rfbetat(1,2) = 1.d0
!!                    gl%rfbetarho(1,2) = 1.d0
!!                    gl%rfgammat(1,2) = 1.d0
!!                    gl%rfgammarho(1,2) = 1.d0
!!                    gl%rfbetat(2,1) = 1.d0
!!                    gl%rfbetarho(2,1) = 1.d0
!!                    gl%rfgammat(2,1) = 1.d0
!!                    gl%rfgammarho(2,1) = 1.d0
!!                elseif (gl%LIN_or_LB == 2) then
!!                    !Linear mixing rules
!!                    gl%rfbetat(1,2) = 1.d0
!!                    gl%rfbetarho(1,2) = 1.d0
!!                    gl%rfbetat(2,1) = 1.d0
!!                    gl%rfbetarho(2,1) = 1.d0
!!                    gl%rfgammat(1,2) = 0.5d0*(gl%tc(1)+gl%tc(2))/dsqrt(gl%tc(1)*gl%tc(2))
!!                    gl%rfgammat(2,1) = gl%rfgammat(1,2)
!!                    gl%rfgammarho(1,2) = 4.d0*(1.d0/gl%rhoc(1)+1.d0/gl%rhoc(2))/(gl%rhoc(1)**(-1.d0/3.d0)+gl%rhoc(2)**(-1.d0/3.d0))**3
!!                    gl%rfgammarho(2,1) = gl%rfgammarho(1,2)
!!                else
!!                    errorflag = -11223344   !Dummy
!!                    return
!!                End if
!!            end if
!!
!!            !Determine the start indices for the groups of each pure fluid when they are written in one common vector
!!            gl%group_start_index(1) = 1
!!            Do i = 2, gl%ncomp
!!                gl%group_start_index(i) = gl%group_start_index(i-1) + gl%nr_of_groups_i(i-1)
!!            End do
!!            m = 0
!!
!!            !Read the van der Waals volume RK and the van der Waals area Qk of the groups in the mixture
!!            if (gl%mix_type == 22) then
!!                filename = trim(path) // 'ge_models/unifac/psrk/rk_qk.par'
!!            elseif (gl%mix_type == 12) then
!!                filename = trim(path) // 'ge_models/unifac/helm_ge/rk_qk.par'
!!            end if
!!            inquire (file = filename, exist = exists)
!!            if (.not.exists) then
!!                errorflag = -8879
!    integer :: firstsubst, secondsubst, j, readerr
!    integer :: dfnreg, dfnlreg, dfnspec, dfnlspec, dummynr1,dummynr2,dummynr3, dfngauss, dfnlgauss,dfpli,dftli, ndummy
!    ndummy = 0
!    readerr = 0
!!                return
!!            end if
!!            open(newunit=geuni_unit, file=trim(filename), status='old', action='read', iostat=error)
!!            if (error /= 0) then
!!                errorflag = -8879
!!                return
!        read(dummy,*,iostat = readerr) modeltype,betat,gammat,betav,gammav,fij_dat,ndummy !(3)
!        if (readerr .ne. 0) then
!            read(dummy,*) modeltype,betat,gammat,betav,gammav,fij_dat!,ndummy !(3)
!            readerr = 0
!        end if
!!            end if
!!
!!            read(geuni_unit,*,end=900) dummy_SRK
!!            Do j = 1, 2000
!!                read(geuni_unit,*,end=900) dummy_SRK
!!                if (dummy_SRK(1:4) == '@END') exit
!!                backspace(geuni_unit)
!!                !Read the next line of the file and check if the fluids in the mixture have this group
!!                read(geuni_unit,*,iostat = error) nr_main_group, nr_sub_group, vdW_vol_Rk, vdW_ar_Qk
!!                if (error /= 0) then
!!                    errorflag = -7896
!!                    close(geuni_unit)
!!                    return
!!                end if
!!
!!                if(.not.allocated(gl%R_ik)) then
!!                    allocate(gl%R_ik(gl%ncomp,maxval(gl%nr_of_groups_i)))
!!                    allocate(gl%Q_ik,gl%v_ik,mold=gl%R_ik)
!!                end if
!!
!!                do i=1,gl%ncomp
!!                    do k = 1, gl%nr_of_groups_i(i)
!!                        if (gl%subgroups_ik_list(i,k) == nr_sub_group) then
!!                            gl%R_ik(i,k) = vdW_vol_Rk
!!                            gl%Q_ik(i,k) = vdW_ar_Qk
!!                            gl%v_ik(i,k) = 1   !ALL GROUPS ARE TREATED SEPERATELY AT THE MOMENT, EVEN IF THEY EXIST MULTIPLE TIMES IN A SUBSTANCE. THIS COULD BE CHANGED IN THE FUTURE BY COUNTING HOW OFTEN CERTAIN GROUPS OCCUR, BUT IT IS NOT CLEAR AT THE MOMENT IF THIS IS ADVANTAGEOUS OR NOT
!!                            !sort the main group into the correct position of the vector maingroups_mix_list
!!                            gl%maingroups_mix_list(gl%group_start_index(i) + k - 1) = nr_main_group
!!                            !Count how many groups have already been assigned. If all groups have been found, quit reading
!!                            m = m + 1
!!                            if (m == gl%nr_of_groups_mix) exit
!!                        end if
!!                    end do
!!                    if (m == gl%nr_of_groups_mix) exit
!!                end do
!!                if (m == gl%nr_of_groups_mix) exit
!!            end do
!!
!!            if (m < gl%nr_of_groups_mix) then
!!                errorflag = -8879
!!                close(geuni_unit)
!!                return
!!            end if
!!
!!            close(geuni_unit)
!!
!!
!!            !Read the binary interactions parameters "anm", "bnm", and "cnm" for all groups in the mixture
!!            if (gl%mix_type == 22) then
!!                filename = trim(path) // 'ge_models/unifac/psrk/Interac.par'
!!            elseif (gl%mix_type == 12) then
!!                filename = trim(path) // 'ge_models/unifac/helm_ge/Interac.par'
!!            end if
!!            inquire (file = filename, exist = exists)
!!            if (.not.exists) then
!!                errorflag = -8879
!!                return
!!            end if
!!            open(newunit=geinter_unit, file=trim(filename), status='old', action='read', iostat=error)
!!            if (error /= 0) then
!!                errorflag = -8879
!!                return
!!            end if
!!
!!            !The binary parameters file has a well defined structure. In the first column the main groups are listed in numerical order (1...85, Index "n"), in the second column
!!            !the main groups interacting with the main groups of the first column are given, starting at (maingroup_column1 + 1), index m. After that, the up to six interaction parameters
!!            !are given: anm, bnm, cnm, amn, bmn, cmn
!!            !Example: maingroup n   maingroup m     anm     bnm     cnm     amn     bmn     cmn
!!            !                   1           2       ....
!!            !                   1           3       ....
!!            !                   1           4       ....
!!            !                   ...
!!            !                   1           85      ....
!!            !                   2           3       ....
!!            !                   2           4       ....
!!            !                   ...
!!            !Thus, if the groups n and m are given for which the interaction parameters are needed, the line where this information stands can be calculated with the formula
!!            ! line = (n-1)*85-n*(n-1)/2+(m-n)
!!            !allocate a_nm var
!!            if(.not.allocated(gl%a_nm)) then
!!                allocate(gl%a_nm(gl%nr_of_groups_mix,gl%nr_of_groups_mix))
!!                allocate(gl%b_nm,gl%c_nm,mold=gl%a_nm)
!!            end if
!!
!!            Do i = 1, gl%nr_of_groups_mix
!!                n = gl%maingroups_mix_list(i)
!!                Do j = i, gl%nr_of_groups_mix
!!                    m = gl%maingroups_mix_list(j)
!!                    if (n < m) then     !This is the way that the groups are stored in the file, no switching necessary
!!                        n_gr = n
!!                        m_gr = m
!!                        line_nr = (n_gr-1)*85-n_gr*(n_gr-1)/2+(m_gr-n_gr)
!!                        !call fseek(10, line_nr, 0, error)      !Unfortunately does not work, because fseek goes through the file entry by entry and not line by line
!!                        Do k = 1, line_nr-1
!!                            Read(geinter_unit,*,iostat = error) n_gr_read     !Dummy read statement
!!                            if (error /= 0) then
!!                                errorflag = -7897
!!                                close(geinter_unit)
!!                                return
!!                            end if
!!                        end do
!!                        if (error == 0) then
!!                            read(geinter_unit,*,iostat = error) n_gr_read, m_gr_read, a_nm_read, b_nm_read, c_nm_read, a_mn_read, b_mn_read, c_mn_read
!!                            if ((n_gr /= n_gr_read) .or. (m_gr /= m_gr_read) .or. (error /= 0)) then
!!                                errorflag = -7897
!!                                close(geinter_unit)
!!                                return
!!                            end if
!!
!!                            gl%a_nm(i,j) = a_nm_read
!!                            gl%b_nm(i,j) = b_nm_read
!!                            gl%c_nm(i,j) = c_nm_read
!!                            gl%a_nm(j,i) = a_mn_read
!!                            gl%b_nm(j,i) = b_mn_read
!!                            gl%c_nm(j,i) = c_mn_read
!!                            rewind(geinter_unit)  !Go back to the beginning of the file
!!                        else
!!                            errorflag = -7897
!!                            close(geinter_unit)
!!                            return
!!                        end if
!!                    elseif (n > m) then              !Switch groups, such that n_gr < m_gr
!!                        n_gr = m
!!                        m_gr = n
!!                        line_nr = (n_gr-1)*85-n_gr*(n_gr-1)/2+(m_gr-n_gr)
!!                        !call fseek(10, line_nr, 0, error)      !Unfortunately does not work, because fseek goes through the file entry by entry and not line by line
!!                        Do k = 1, line_nr-1
!!                            Read(geinter_unit,*,iostat = error) n_gr_read     !Dummy read statement
!!                            if (error /= 0) then
!!                                errorflag = -7897
!!                                close(geinter_unit)
!!                                return
!!                            end if
!!                        end do
!!                        if (error == 0) then
!!                            read(geinter_unit,*,iostat = error) n_gr_read, m_gr_read, a_nm_read, b_nm_read, c_nm_read, a_mn_read, b_mn_read, c_mn_read
!!                            if ((n_gr /= n_gr_read) .or. (m_gr /= m_gr_read) .or. (error /= 0)) then
!!                                errorflag = -7897
!!                                close(geinter_unit)
!!                                return
!!                            end if
!!                            gl%a_nm(j,i) = a_nm_read
!!                            gl%b_nm(j,i) = b_nm_read
!!                            gl%c_nm(j,i) = c_nm_read
!!                            gl%a_nm(i,j) = a_mn_read
!!                            gl%b_nm(i,j) = b_mn_read
!!                            gl%c_nm(i,j) = c_mn_read
!!                            rewind(geinter_unit)  !Go back to the beginning of the file
!!                        else
!!                            errorflag = -7897
!!                            close(geinter_unit)
!!                            return
!!                        end if
!!                    elseif (n == m) then       !Interaction of same main group -> all interaction parameters 0
!!                        gl%a_nm(j,i) = 0.D0
!!                        gl%b_nm(j,i) = 0.D0
!!                        gl%c_nm(j,i) = 0.D0
!!                        gl%a_nm(i,j) = 0.D0
!!                        gl%b_nm(i,j) = 0.D0
!!                        gl%c_nm(i,j) = 0.D0
!!                    end if
!!                end do
!!            end do
!!
!!            !Erik, April 2018
!!            !Read COSMO-SAC files and generate sigma-profiles
!!            !Read index number
!!        elseif (gl%Mix_type == 13) then
!!
!!            !------------------------------------------------
!!            !Set COSMO-SAC Version
!!            gl%cosmo%COSMO_ver =  1
!!            !analytical derivations for COSMO-SAC Version 1?
!!            gl%cosmo%analytical = .false.
!!            !successive substitution or Newton-Raphson methode for solving segment activity coefficients; Newton-Raphson = true; successive substitution = false
!!            gl%cosmo%solv_method = .false.
!!            !------------------------------------------------
!!
!!            interval=51
!!            filename = trim(path) // 'ge_models/cosmo-sac/VT-2005_Sigma_Profile_Database_Index_v2.txt'
!!
!!            inquire (file = filename, exist = exists)
!!            if (.not.exists) then
!!                errorflag = -7898
!!                return
!!            end if
!!
!!            do nrsubst = 1, gl%ncomp
!!                open(unit=10, file=trim(filename), status='old', action='read', iostat=error)
!!                if (error /= 0) then
!!                    errorflag = -7898
!!                    return
!!                end if
!!
!!                index_nr(nrsubst) = 0
!!                index_char(nrsubst) = ""
!!
!!                read(10,*,iostat = error) dummy_cosmo !Dummy read statement
!!                do k = 1, 1600
!!                    Read(10,*,iostat = error) index_nr(nrsubst), CASNR_COSMO, dummy_cosmo, gl%cosmo%Vcosmo(nrsubst)
!!                    !Search Indes number for component via CAS # instead of name or formula
!!                    if (CASNR_COSMO == gl%substcasnr(nrsubst)) then
!!                        close(10)
!!                        write(index_char(nrsubst),'(i4)') index_nr(nrsubst)
!!                        if (index_nr(nrsubst) < 10) then
!!                            index_char(nrsubst) = "000" // index_char(nrsubst)(4:4)
!!                        end if
!!                        if ((index_nr(nrsubst) < 100) .and. (index_nr(nrsubst) > 9)) then
!!                            index_char(nrsubst) = "00" // index_char(nrsubst)(3:4)
!!                        end if
!!                        if ((index_nr(nrsubst) < 1000) .and. (index_nr(nrsubst) > 99)) then
!!                            index_char(nrsubst) = "0" // index_char(nrsubst)(2:4)
!!                        end if
!!                        exit
!!                    end if
!!                    if (error /= 0) then
!!                        errorflag = -7898
!!                        close(10)
!!                        return
!!                    end if
!!                end do
!!                if (index_char(nrsubst) == "") then
!!                    errorflag = -7898
!!                    return
!!                end if
!!            end do
!!
!!            !---------------------------------------------------------------------------------------------------------------------------------
!!            !Read Sigma-Profiles
!!            !---------------------------------------------------------------------------------------------------------------------------------
!!            if (gl%cosmo%COSMO_ver == 1) then
!!                !Read Sigma-Profile Files for components and calculate COSMO surface area Acosmo (COSMO-SAC Version 1)
!!                do nrsubst = 1, gl%ncomp
!!                    !filename = trim(path) // 'ge_models/cosmo-sac/sigma-profiles/VT2005-' // index_char(nrsubst) // '-PROF.txt'
!!                    filename = trim(path) // 'ge_models/cosmo-sac/sigma-profiles/VT2005-' // index_char(nrsubst) // '.txt'
!!                    inquire (file = filename, exist = exists)
!!                    if (.not.exists) then
!!                        errorflag = -7898
!!                        return
!!                    end if
!!
!!                    open(unit=10, file=trim(filename), status='old', action='read', iostat=error)
!!                    if (error /= 0) then
!!                        errorflag = -7898
!!                        return
!!                    end if
!!                    gl%cosmo%Acosmo(nrsubst) = 0
!!                    Do k = 1, interval
!!                        Read(10,*,iostat = error) gl%cosmo%counter(k), gl%cosmo%sigma(k, nrsubst)
!!                        gl%cosmo%Acosmo(nrsubst) = gl%cosmo%Acosmo(nrsubst) + gl%cosmo%sigma(k, nrsubst)
!!                        if (error /= 0) then
!!                            errorflag = -7898
!!                            close(10)
!!                            return
!!                        end if
!!                    end do
!!                    !November 2018, Erik
!!                    !read min and max values for occupied intervals
!!                    Read(10,*,iostat = error) counter_min(nrsubst,1)
!!                    Read(10,*,iostat = error) counter_max(nrsubst,1)
!!                    profile_count = 1
!!
!!                    !Read total Volume of cavity
!!                    read(10,*,iostat = error) gl%cosmo%Vcosmo(nrsubst)
!!                    close(10)
!!                end do
!!            elseif ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3) .or. (gl%cosmo%COSMO_ver == 4)) then
!!                !Read Sigma-Profile Files for components and calculate COSMO surface area Acosmo (COSMO-SAC Version 2 and 3)
!!                !update version 4 (Xiong), Feb. 7th 2018
!!                do nrsubst = 1, gl%ncomp
!!                    gl%cosmo%Acosmo(nrsubst) = 0
!!                    if (gl%cosmo%COSMO_ver == 4) then
!!                        filename = trim(path) // 'ge_models/cosmo-sac/sigma-profiles_ver_4/VT2005-' // index_char(nrsubst) // '.txt'
!!                    else
!!                        filename = trim(path) // 'ge_models/cosmo-sac/sigma-profiles3/VT2005-' // index_char(nrsubst) // '.txt'
!!                    end if
!!                    inquire (file = filename, exist = exists)
!!                    if (.not.exists) then
!!                        errorflag = -7898
!!                        return
!!                    end if
!!
!!                    open(unit=10, file=trim(filename), status='old', action='read', iostat=error)
!!                    if (error /= 0) then
!!                        errorflag = -7898
!!                        return
!!                    end if
!!                    do i = 1, 3
!!                        do k = 1, interval
!!                            Read(10,*,iostat = error) gl%cosmo%counter_v23(k, i), gl%cosmo%sigma_v23(k, nrsubst, i)
!!                            gl%cosmo%Acosmo(nrsubst) = gl%cosmo%Acosmo(nrsubst) + gl%cosmo%sigma_v23(k, nrsubst, i)
!!                            if (error /= 0) then
!!                                errorflag = -7898
!!                                close(10)
!!                                return
!!                            end if
!!                        end do
!!                    end do
!!                    !Read dispersive interaction parameters
!!                    if ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3)) then
!!                        read(10,*,iostat = error) gl%cosmo%molecule_type(nrsubst)
!!                        read(10,*,iostat = error) gl%cosmo%eps_molecule(nrsubst)
!!                        if (error /= 0) then
!!                            errorflag = -7898
!!                            close(10)
!!                            return
!!                        end if
!!                    elseif (gl%cosmo%COSMO_ver == 4) then
!!                        read(10,*,iostat = error) gl%cosmo%molecule_type(nrsubst)
!!                        do j = 1, 17
!!                            read(10,*,iostat = error) dummy_cosmo, gl%cosmo%m_tau(nrsubst,j)
!!                            if (error /= 0) then
!!                                errorflag = -7898
!!                                close(10)
!!                                return
!!                            end if
!!                        end do
!!                    end if
!!                    !November 2018, Erik
!!                    !read min and max values for occupied intervals
!!                    do i = 1, 3
!!                        Read(10,*,iostat = error) counter_min(nrsubst,i)
!!                        Read(10,*,iostat = error) counter_max(nrsubst,i)
!!                    end do
!!                    profile_count = 3
!!
!!                    !Read total Volume of cavity
!!                    read(10,*,iostat = error) gl%cosmo%Vcosmo(nrsubst)
!!                    close(10)
!!                end do
!!            end if
!!
!!            !---------------------------------------------------------------------------------------------------------------------------------
!!            !November 2018, Erik
!!            !min and max value for occupied intervals for all compounds
!!            !Min
!!            do i = 1, profile_count
!!                gl%cosmo%interval_min(i) = 51
!!                do nrsubst = 1, gl.ncomp
!!                    if ((counter_min(nrsubst,i) < gl%cosmo%interval_min(i)) .and. (counter_min(nrsubst,i) > 0.D0))  then
!!                        gl%cosmo%interval_min(i) = counter_min(nrsubst,i)
!!                    end if
!!                end do
!!                !set value for selected profile to 0, if not occupied at all
!!                if (gl%cosmo%interval_min(i) == 51) then
!!                    gl%cosmo%interval_min(i) = 0
!!                end if
!!            end do
!!            !Max
!!            do i = 1, profile_count
!!                gl%cosmo%interval_max(i) = 0
!!                do nrsubst = 1, gl.ncomp
!!                    if (counter_max(nrsubst,i) > gl%cosmo%interval_max(i)) then
!!                        gl%cosmo%interval_max(i) = counter_max(nrsubst,i)
!!                    end if
!!                end do
!!            end do
!!
!!            !Erik, December 2018
!!            !to increase calculation speed allocate large matrices with the lowest size
!!            !very extensive to change, because everything would have to be rearranged
!!            !gl%cosmo%int = gl%cosmo%interval_max(1) - gl%cosmo%interval_min(1) + 1
!!
!!            !---------------------------------------------------------------------------------------------------------------------------------
!!            !Calculate electrostatic interaction
!!            !---------------------------------------------------------------------------------------------------------------------------------
!!
!!            if ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3)) then    !Hsieh
!!                aeff = 7.25    !surface area of a standard segment in Angstrom^2
!!            elseif (gl%cosmo%COSMO_ver == 1) then    !Mullins
!!                aeff = 7.5    !surface area of a standard segment in Angstrom^2
!!            elseif (gl%cosmo%COSMO_ver == 4) then    !Xiong
!!                aeff = 7.9    !surface area of a standard segment in Angstrom^2
!!            end if
!!
!!            !Version 1
!!
!!            if (gl%cosmo%COSMO_ver == 1) then
!!                !Calculate electrostatic interaction COSMO-SAC Version 1
!!                eo=2.395D-4
!!                eps=3.667   !(LIN AND SANDLER USE A CONSTANT FPOL WHICH YEILDS EPS=3.68)
!!                sigma_hb=0.0084
!!                c_hb=85580.0
!!                fpol = (eps-1.0)/(eps+0.5)
!!                alpha = (0.3*aeff**(1.5))/(eo)
!!                alpha_prime = fpol*alpha
!!                !do o = 1, 1000
!!                !do j = begin_sigma_mix(1), end_sigma_mix(1)
!!                ! do k = begin_sigma_mix(1), end_sigma_mix(1)
!!                do j = gl%cosmo%interval_min(1), gl%cosmo%interval_max(1)
!!                    do k = gl%cosmo%interval_min(1), gl%cosmo%interval_max(1)
!!                        if (gl%cosmo%counter(j)>=gl%cosmo%counter(k)) then
!!                            sigmaacc= gl%cosmo%counter(j)
!!                            sigmadon = gl%cosmo%counter(k)
!!                        end if
!!                        if (gl%cosmo%counter(j)<gl%cosmo%counter(k)) then
!!                            sigmadon = gl%cosmo%counter(j)
!!                            sigmaacc = gl%cosmo%counter(k)
!!                        end if
!!
!!                        gl%cosmo%delta_w_gl(j,k) = (alpha_prime/2.0)*(gl%cosmo%counter(j)+gl%cosmo%counter(k))**2.0 + c_hb *   &
!!                            max(0.D0,(sigmaacc - sigma_hb)) * min(0.D0,(sigmadon + sigma_hb))
!!                    end do
!!                end do
!!            end if
!!
!!            !Version 2 to 4
!!            if ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3) .or. (gl%cosmo%COSMO_ver == 4)) then
!!                c_OH_OH = 4013.78 * 4184.D0     !hydrogen bonding interaction coefficient
!!                c_OT_OT = 932.31 * 4184.D0      !hydrogen bonding interaction coefficient
!!                c_OH_OT = 3016.43 * 4184.D0     !hydrogen bonding interaction coefficient
!!                !in Versions 2 and 3 of COSMO-SAC electrostatic interactions are temperature dependend. This is formulated in c_ES and has to be calculated
!!                !in COSMO_SAC_CALC. The rest can be calculated here
!!
!!                do k = 1, 3
!!                    if (gl%cosmo%interval_min(k) == 0) then
!!                        cycle
!!                    end if
!!                    do l = 1, 3
!!                        if (gl%cosmo%interval_min(l) == 0) then
!!                            cycle
!!                        end if
!!                        do i = gl%cosmo%interval_min(k), gl%cosmo%interval_max(k)
!!                            do j = gl%cosmo%interval_min(l), gl%cosmo%interval_max(l)
!!                                if ((k == 2) .and. (l == 2) .and. (gl%cosmo%counter_v23(i,k) * gl%cosmo%counter_v23(j,l) < 0)) then
!!                                    c_hb = c_OH_OH
!!                                elseif ((((k == 2) .and. (l == 3)) .or. ((k == 3) .and. (l == 2))) .and. (gl%cosmo%counter_v23(i,k) * gl%cosmo%counter_v23(j,l) < 0)) then
!!                                    c_hb = c_OH_OT
!!                                elseif ((k == 3) .and. (l == 3) .and. (gl%cosmo%counter_v23(i,k) * gl%cosmo%counter_v23(j,l) < 0)) then
!!                                    c_hb = c_OT_OT
!!                                else
!!                                    c_hb = 0.D0
!!                                end if
!!                                gl%cosmo%sum_square_counter(i,j,k,l) = (gl%cosmo%counter_v23(i,k) + gl%cosmo%counter_v23(j,l)) ** 2.D0
!!                                gl%cosmo%hb_inter(i,j,k,l) = c_hb * (gl%cosmo%counter_v23(i,k) - gl%cosmo%counter_v23(j,l)) ** 2.D0
!!
!!                                !delta_w_v23(i,j,k,l) = c_es * (gl%cosmo%counter_v23(i,k) + gl%cosmo%counter_v23(j,l)) ** 2 - c_hb * (gl%cosmo%counter_v23(i,k) - gl%cosmo%counter_v23(j,l)) ** 2
!!                            end do
!!                        end do
!!                    end do
!!                end do
!!
!!            end if
!!
!!            !---------------------------------------------------------------------------------------------------------------------------------
!!
!!            !---------------------------------------------------------------------------------------------------------------------------------
!!            !Read NBT experimental Data for pure fluids (Version 4 (Xiong)
!!            if (gl%cosmo%COSMO_ver == 4) then
!!                filename = trim(path) // 'ge_models/cosmo-sac/NBT_pure.txt'
!!                do nrsubst = 1, gl%ncomp
!!                    open(unit=12, file=trim(filename), status='old', action='read', iostat=error)
!!                    if (error /= 0) then
!!                        errorflag = -7898
!!                        return
!!                    end if
!!                    read(12,*,iostat = error) dummy_cosmo
!!                    read(12,*,iostat = error) dummy_cosmo
!!                    do i = 1, 2000
!!                        read(12,*,iostat = error) dummy_cosmo, gl%cosmo%NBT_cosmo(nrsubst), dummy_cosmo, dummy_cosmo, dummy_cosmo, dummy_cosmo, CASNR_COSMO
!!                        if (trim(CASNR_COSMO) == trim(gl%substcasnr(nrsubst))) then
!!                            exit
!!                        end if
!!                        if (error /= 0) then
!!                            errorflag = -7898
!!                            close(12)
!!                            return
!!                        end if
!!                    end do
!!                    close(12)
!!                end do
!!
!!            end if
!!
!!            !if (LIN_or_LB == 1) then
!!            !    !LB-mixing rules
!!            !    rfbetat(1,2) = 1.d0
!!            !    rfbetarho(1,2) = 1.d0
!!            !    rfgammat(1,2) = 1.d0
!!            !    rfgammarho(1,2) = 1.d0
!!            !    rfbetat(2,1) = 1.d0
!!            !    rfbetarho(2,1) = 1.d0
!!            !    rfgammat(2,1) = 1.d0
!!            !    rfgammarho(2,1) = 1.d0
!!            !elseif (LIN_or_LB == 2) then
!!            !Linear mixing rules
!!            gl%rfbetat(1,2) = 1.d0
!!            gl%rfbetarho(1,2) = 1.d0
!!            gl%rfbetat(2,1) = 1.d0
!!            gl%rfbetarho(2,1) = 1.d0
!!            gl%rfgammat(1,2) = 0.5d0*(gl%tc(1)+gl%tc(2))/dsqrt(gl%tc(1)*gl%tc(2))
!!            gl%rfgammat(2,1) = gl%rfgammat(1,2)
!!            gl%rfgammarho(1,2) = 4.d0*(1.d0/gl%rhoc(1)+1.d0/gl%rhoc(2))/(gl%rhoc(1)**(-1.d0/3.d0)+gl%rhoc(2)**(-1.d0/3.d0))**3
!!            gl%rfgammarho(2,1) = gl%rfgammarho(1,2)
!!            !else
!!            !    errorflag = -11223344   !Dummy
!!            !    return
!!            !End if
!!
!!
!!            !_________________________________
!!
!!            !Read mixing rules for the PR
!!        Elseif ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then
!!
!!            !Andreas, November 2015
!!            !Sebastian, Oct 2016 Killed the line below
!!            !REQ = R_PR
!!
!!            filename = trim(path) // 'PR/BINARY_PR/BIN_PR.MIX'  ! check for binary mixing parameters file for the SRK
!!            inquire (file = filename, exist = exists)
!!            if (.not.exists) then
!!                errorflag = -8878
!!                return
!!            end if
!!            open(newunit=prmix_unit, file=trim(filename), status='old', action='read', iostat=error)
!!            if (error /= 0) then
!!                errorflag = -7888
!!                return
!!            end if
!!
!!
!!            !Skip the first three line in the file
!!            read(prmix_unit,*,end=900) dummy_PR
!!            read(prmix_unit,*,end=900) dummy_PR
!!            read(prmix_unit,*,end=900) dummy_PR
!!
!!            Do i = 1, 1000
!!                read(prmix_unit,*,end=900) dummy_PR
!!                if (dummy_PR(1:4) == '@END') exit
!!                backspace(prmix_unit)
!!
!!                !Adjusted read statement because the covolume parameter "b" can now also be mixed with a quadratic mixing rule. For this, the parameter lij was introduced in the file bin_pr.mix
!!                !Andreas J輍er, March 2017
!!
!!                !Try new reading routine first
!!                read(prmix_unit,*,iostat = error)  CASNR_1, CASNR_2, kij_read, lij_read
!!                !Use old routine if error occurs. Safety net, can maybe be deleted in the future.
!!                if (error /= 0) then
!!                    backspace(prmix_unit)
!!                    read(prmix_unit,*,iostat = error)  CASNR_1, CASNR_2, kij_read
!!                    !If no 4th parameter is given, set lij to 0
!!                    lij_read = 0.D0
!!                    if (error/=0) then
!!                        errorflag = -7895
!!                    end if
!!                end if
!!
!!                do nrsubst = 1, gl%ncomp
!!                    Do j = nrsubst + 1, gl%ncomp
!!                        if (((CASNR_PR(j) == CASNR_1) .and. (CASNR_PR(nrsubst) == CASNR_2)) .or. ((CASNR_PR(j) == CASNR_2) .and. (CASNR_PR(nrsubst) == CASNR_1))) then
!!                            gl%kij_PR(j,nrsubst) = kij_read
!!                            gl%kij_PR(nrsubst,j) = kij_read
!!                            gl%lij_PR(j,nrsubst) = lij_read
!!                            gl%lij_PR(nrsubst,j) = lij_read
!!                            backspace(prmix_unit)
!!                            read(prmix_unit,'(A)')dummy
!!                            k = 1
!!                            do while (dummy(k:k) /= '!')
!!                                k = k + 1
!!                                if (k > len(dummy)) exit
!!                            enddo
!!                            if (k > len(dummy)) then
!!                                if (allocated(gl%litref)) then
!!                                    gl%litref%lit_ref_mix = 'no literature reference documented'
!!                                end if
!!                            else
!!                                if (allocated(gl%litref)) then
!!                                    gl%litref%lit_ref_mix = dummy(k+1:)
!!                                end if
!!                            endif
!!                        end if
!!                    End do
!!                end do
!!            End do
!!
!!            close(prmix_unit)
!!
!!
!!            !Read mixing rules for LKP
!!        Elseif (gl%Mix_type == 4) then
!!
!!            !Sebastian, Oct 2016 Killed the line below
!!            !REQ = 8.31451D0
!!
!!            !changed for compatibility with gfortran
!!            !filename = trim(path) // 'LKP\BINARY_LKP\BIN_LKP.MIX'  ! check for binary mixing parameters file for the LKP
!!            filename = trim(path) // 'lkp/binary_lkp/bin_lkp.mix'  ! check for binary mixing parameters file for the LKP
!!            inquire (file = filename, exist = exists)
!!            if (.not.exists) then
!!                errorflag = -8878
!!                return
!!            end if
!!            open(newunit=lkpmix_unit, file=trim(filename), status='old', action='read', iostat=error)
!!            if (error /= 0) then
!!                errorflag = -7889
!!                return
!!            end if
!!
!!            !Skip the first three line in the file
!!            read(lkpmix_unit,*,end=900) dummy_LKP
!!            read(lkpmix_unit,*,end=900) dummy_LKP
!!            read(lkpmix_unit,*,end=900) dummy_LKP
!!
!!            Do i = 1, 1000
!!                read(lkpmix_unit,*,end=900) dummy_LKP
!!                if (dummy_LKP(1:4) == '@END') exit
!!                backspace(lkpmix_unit)
!!
!!                !New identification of mixing rules using the CAS numbers of the fluids
!!                read(lkpmix_unit,*)  CASNR_1, CASNR_2, kij_read
!!                do nrsubst = 1, gl%ncomp
!!                    Do j = nrsubst + 1, gl%ncomp
!!                        if ((CASNR_LKP(j) == CASNR_1) .and. (CASNR_LKP(nrsubst) == CASNR_2)) then
!!                            gl%kij_LKP(j,nrsubst) = kij_read
!!                            gl%kij_LKP(nrsubst,j) = kij_read
!!                            backspace(lkpmix_unit)
!!                            read(lkpmix_unit,'(A)')dummy
!!                            k = 1
!!                            do while (dummy(k:k) /= '!')
!!                                k = k + 1
!!                            enddo
!!                            if (allocated(gl%litref)) then
!!                                gl%litref%lit_ref_mix = dummy(k+1:)
!!                            endif
!!                        end if
!!                        if ((CASNR_LKP(j) == CASNR_2) .and. (CASNR_LKP(nrsubst) == CASNR_1)) then
!!                            gl%kij_LKP(j,nrsubst) = kij_read
!!                            gl%kij_LKP(nrsubst,j) = kij_read
!!                            backspace(lkpmix_unit)
!!                            read(lkpmix_unit,'(A)')dummy
!!                            k = 1
!!                            do while (dummy(k:k) /= '!')
!!                                k = k + 1
!!                            enddo
!!                            if (allocated(gl%litref)) then
!!                                gl%litref%lit_ref_mix = dummy(k+1:)
!!                            endif
!!                        end if
!!                    End do
!!                end do
!!            End do
!!
!!            close(lkpmix_unit)
!!        End if
!!
!!    End if
!!
!!
!!    errorflag = error
!!    return
!!    !only enter if end of file is reached during read:
!!900 errorflag = -7880
!!    close(rkm_unit)
!!    close(srk_unit)
!!    !write(*,*) 'Error: end of fluid-file is reached during read'
!!
!!
!!
!!    end subroutine read_from_fluidfiles
!!
!!
!!
!!    !called by sub "read_from_fluidfiles", need not to be called from somewhere else:
!!    subroutine read_mix_file(gl,substance1,substance2,filepath,internalioerror)
!!
!!
!!
!!
!!
!!
!!    implicit none
!!
!!    type(type_gl) :: gl
!!
!!
!!    character(len=*) :: filepath          !has to be transfered to this subroutine
!!    integer :: substance1,substance2    !has to be transfered to this subroutine
!!    integer :: internalioerror                 !no error as startvalue
!!
!!    character(len=255) :: filename1
!!    character(len=255) :: filename2
!!    character(len=255) :: dummy
!!    character(len=3) :: modeltype
!!    integer :: firstsubst, secondsubst, j
!!    integer :: dfnreg, dfnlreg, dfnspec, dfnlspec, dummynr1,dummynr2,dummynr3, dfngauss, dfnlgauss,dfpli,dftli
!!    double precision :: betat,gammat,betav,gammav,fij_dat
!!    double precision :: dfni,dfti,dfdi,dfli,dfetai,dfepsi,dfbetai,dfgammai,dfgetai,dfgepsi,dfgbetai,dfggami,dfgdummy1,dfgdummy2,dfgdummy3, dfpi
!!    integer:: mix_unit, pos
!!    filename1=''
!!    filename2=''
!!    dummy=''
!!    modeltype=''
!!
!!    firstsubst = 0
!!    secondsubst = 0
!!    dfnreg = 0
!!    dfnlreg = 0
!!    dfnspec = 0
!!    dfnlspec = 0
!!    dummynr1 = 0
!!    dummynr2 = 0
!!    dummynr3 = 0
!!    betat=0.d0
!!    gammat=0.d0
!!    betav=0.d0
!!    gammav=0.d0
!!    fij_dat=0.d0
!!    dfni=0.d0
!!    dfti=0.d0
!!    dfdi=0.d0
!!    dfli=0.D0
!!    dfetai=0.d0
!!    dfepsi=0.d0
!!    dfbetai=0.d0
!!    dfgammai=0.d0
!!    mix_unit = 0
!!    dfpi = 0.d0
!!    pos = 0
!!
!!    filename1 = trim(filepath) // trim(gl%components(substance1)) //'-'// &
!!        trim(gl%components(substance2)) //'.mix'
!!    filename2 = trim(filepath) // trim(gl%components(substance2)) //'-'// &
!!        trim(gl%components(substance1)) //'.mix'
!!
!!
!!    open(newunit=mix_unit, file=trim(filename1), status='old', action='read', iostat=internalioerror)
!!    firstsubst = substance1
!!    secondsubst = substance2
!!    if (internalioerror /= 0) then !file not found, so try other way round
!!        open(newunit=mix_unit, file=trim(filename2), status='old', action='read', iostat=internalioerror)
!!        firstsubst = substance2
!!        secondsubst = substance1
!!        if(internalioerror /= 0) then !next file also not found so give default values and exit the subroutine
!!            !Andreas Nov 2013
!!            !If no file for the binary mixture exists, exit with error
!!            internalioerror = -7881
!!
!!            !write(*,*) 'Error: mix-files not found: ', trim(filename1), ' or ', trim(filename2)
!!            return
!!        endif
!!    endif
!!
!!    if (internalioerror == 0) then
!!
!!        read(mix_unit,*,end=999) dummy !information not used (1)
!!        !if (dummy(1:1) /= "?") then
!!        !    internalioerror = -7882
!!        !    return
!!        !end if
!!        if (allocated(gl%litref)) then
!!            read(mix_unit,'(a255)',end=999)gl%litref%lit_ref_mix
!!        else
!!            read(mix_unit,*,end=999) dummy !information not used (2)
!!        endif
!!        read(mix_unit,'(a255)',end=999) dummy
!!        pos = index(dummy,' ')
!!        dummy = dummy(pos:)
!!        read(dummy,*) modeltype,betat,gammat,betav,gammav,fij_dat,dummy !(3)
!!        !if not modeltype not "KW" then set default values to 1:
!!        If ((modeltype(1:2) == 'KW').or.(modeltype(1:1) == "B")) then
!!            gl%rfbetat(firstsubst,secondsubst) = betat
!!            gl%rfbetarho(firstsubst,secondsubst) = betav !gammat
!!            gl%rfgammat(firstsubst,secondsubst) = gammat !betav
!!            gl%rfgammarho(firstsubst,secondsubst) = gammav
!!            gl%rfbetat(secondsubst,firstsubst) = 1/betat
!!            gl%rfbetarho(secondsubst,firstsubst) = 1/betav !gammat
!!            gl%rfgammat(secondsubst,firstsubst) = gammat !1/betav
!!            gl%rfgammarho(secondsubst,firstsubst) = gammav
!!        else
!!            gl%rfbetat(firstsubst,secondsubst) = 1.d0
!!            gl%rfbetarho(firstsubst,secondsubst) = 1.d0
!!            gl%rfgammat(firstsubst,secondsubst) = 1.d0
!!            gl%rfgammarho(firstsubst,secondsubst) = 1.d0
!!            gl%rfbetat(secondsubst,firstsubst) = 1.d0
!!            gl%rfbetarho(secondsubst,firstsubst) = 1.d0
!!            gl%rfgammat(secondsubst,firstsubst) = 1.d0
!!            gl%rfgammarho(secondsubst,firstsubst) = 1.d0
!!        endif
!!
!!        !read the last part with departure function coefficients only if not modeltype is "KW0" and Fij=0:
!!        !if ((modeltype == 'KW0').and.(int(fij_dat) == 0))then    !Problem: if Fij < 0.5 e.g. 0.13 its rounded down to 0 due to the int statement Andreas Aug2010
!!        if ((((modeltype(1:2) == 'KW').or.(modeltype(1:1) == "B")).and.(fij_dat == 0.0d0)))then
!!            do
!!                read(mix_unit,*,end=998) dummy !information not used
!!                if ((trim(dummy(1:5)) == '?LITE').and.(allocated(gl%litref))) then
!!                    read(mix_unit,'(A)')dummy
!!                    do while (trim(dummy(1:5)) /= '!end ')
!!                        gl%litref%lit_ref_mix = trim(gl%litref%lit_ref_mix) // ' ' // dummy(2:)
!!                        read(mix_unit,'(A)')dummy
!!                    end do
!!                endif
!!            end do
!!998         close(mix_unit)
!!            return
!!        endif
!!
!!        goover1: do !jump over comments until "!end" or "! end!
!!            read(mix_unit,*,end=999) dummy !information not used
!!            call uppertolower_char(dummy, len(dummy))
!!            if ((dummy(1:4) == '!end').or.(dummy(1:5) == '! end')) then
!!                exit goover1
!!                !endif
!!                !if (dummy(1:1) == '!') then
!!                !                exit goover2
!!            elseif ((trim(dummy(1:5)) == '?````').and.(allocated(gl%litref))) then
!!                read(mix_unit,'(A)')dummy
!!                do while (trim(dummy(1:5)) /= '!````')
!!                    gl%litref%lit_ref_mix = trim(gl%litref%lit_ref_mix) // dummy(2:)
!!                    read(mix_unit,'(A)')dummy
!!                end do
!!                backspace(mix_unit)
!!            endif
!!        enddo goover1
!!
!!        read(mix_unit,*,end=999) dummy !information not used
!!        if ((modeltype(1:2) /= 'KW').and.(modeltype(1:1) /= "B")) then !use the default values:
!!            read(mix_unit,*,end=999) betat,gammat,betav,gammav,gl%fij,dummy
!!            gl%rfbetat(firstsubst,secondsubst) = betat
!!            gl%rfbetarho(firstsubst,secondsubst) = gammat
!!            gl%rfgammat(firstsubst,secondsubst) = betav
!!            gl%rfgammarho(firstsubst,secondsubst) = gammav
!!            gl%rfbetat(secondsubst,firstsubst) = 1/betat
!!            gl%rfbetarho(secondsubst,firstsubst) = gammat
!!            gl%rfgammat(secondsubst,firstsubst) = 1/betav
!!            gl%rfgammarho(secondsubst,firstsubst) = gammav
!!        else
!!            read(mix_unit,*,end=999) dummy !information not used
!!        endif
!!
!!        read(mix_unit,*,end=999) dfnreg, dfnlreg, dummynr1, dfnspec, dfnlspec, dfngauss, dfnlgauss
!!        if ((dfnreg+dfnspec+dfngauss) > 0) then !departure function coefficients exist
!!            !Andreas Aug 2010
!!            gl%dfpol(firstsubst,secondsubst) = dfnreg
!!            gl%dfpol(secondsubst,firstsubst) = dfnreg
!!            gl%dfexp(firstsubst,secondsubst) = dfnspec
!!            gl%dfexp(secondsubst,firstsubst) = dfnspec
!!            gl%dfgau(firstsubst,secondsubst) = dfngauss
!!            gl%dfgau(secondsubst,firstsubst) = dfngauss
!!            gl%Fij(firstsubst,secondsubst) = fij_dat
!!            gl%Fij(secondsubst,firstsubst) = fij_dat
!!            !----------------------------------------
!!        endif
!!
!!
!!        do j=1,dfnreg !for normal terms
!!            if (dfnlreg == 4) then  !checks, if a value for dfli exists for this mixture - J.Gernert, Sept. 2010
!!                read(mix_unit,*) dfni,dfti,dfdi, dfli
!!            elseif (dfnlreg == 5) then
!!                read(mix_unit,*) dfni,dfti,dfdi, dfli, dfpi
!!            else
!!                read(mix_unit,*) dfni,dfti,dfdi
!!            end if
!!            !---------------------------------------
!!            !Andreas Aug 2010
!!            gl%dfn(j,firstsubst,secondsubst) = dfni
!!            gl%dfn(j,secondsubst,firstsubst) = dfni
!!            gl%dfd(j,firstsubst,secondsubst) = dfdi
!!            gl%dfd(j,secondsubst,firstsubst) = dfdi
!!            gl%dft(j,firstsubst,secondsubst) = dfti
!!            gl%dft(j,secondsubst,firstsubst) = dfti
!!            gl%dfl(j,firstsubst,secondsubst) = dfli
!!            gl%dfl(j,secondsubst,firstsubst) = dfli
!!            gl%dfp(j,firstsubst,secondsubst) = dfpi
!!            gl%dfp(j,secondsubst,firstsubst) = dfpi
!!
!!            select case(gl%dfd_coeff_structure(firstsubst,secondsubst))
!!
!!            case (coeff_structure_all_zeros)
!!                if (dfdi /= 0d0) then
!!                    if (DBLE(INT(dfdi)) /= dfdi) then
!!                        gl%dfd_coeff_structure(firstsubst,secondsubst) = coeff_structure_real
!!                    else
!!                        gl%dfd_coeff_structure(firstsubst,secondsubst) = coeff_structure_int
!!                    endif
!!                endif
!!
!!            case (coeff_structure_int)
!!                if (DBLE(INT(dfdi)) /= dfdi) then
!!                    gl%dfd_coeff_structure(firstsubst,secondsubst) = coeff_structure_real
!!                endif
!!
!!            end select
!!
!!
!!            select case(gl%dfl_coeff_structure(firstsubst,secondsubst))
!!
!!            case (coeff_structure_all_zeros)
!!                if (dfli /= 0d0) then
!!                    if (DBLE(INT(dfli)) /= dfli) then
!!                        gl%dfl_coeff_structure(firstsubst,secondsubst) = coeff_structure_real
!!                    else
!!                        gl%dfl_coeff_structure(firstsubst,secondsubst) = coeff_structure_int
!!                    endif
!!                endif
!!
!!            case (coeff_structure_int)
!!                if (DBLE(INT(dfli)) /= dfli) then
!!                    gl%dfl_coeff_structure(firstsubst,secondsubst) = coeff_structure_real
!!                endif
!!
!!            end select
!!
!!            !---------------------------------------
!!        enddo
!!
!!        gl%dfd_coeff_structure(secondsubst,firstsubst) = gl%dfd_coeff_structure(firstsubst,secondsubst)
!!        !if (gl%dfd_coeff_structure(firstsubst,secondsubst) == coeff_structure_int) then
!!        !    if (.not.allocated(gl%dfd_int)) allocate(gl%dfd_int(25,gl%ncomp,gl%ncomp))
!!        !    gl%dfd_int(:,firstsubst,secondsubst) = int(gl%dfd(:,firstsubst,secondsubst))
!!        !    gl%dfd_int(:,secondsubst,firstsubst) = int(gl%dfd(:,firstsubst,secondsubst))
!!        !endif
!!
!!
!!        gl%dfl_coeff_structure(secondsubst,firstsubst) = gl%dfl_coeff_structure(firstsubst,secondsubst)
!!        !if (gl%dfl_coeff_structure(firstsubst,secondsubst) == coeff_structure_int) then
!!        !    if (.not.allocated(gl%dfl_int)) allocate(gl%dfl_int(25,gl%ncomp,gl%ncomp))
!!        !    gl%dfl_int(:,firstsubst,secondsubst) = int(gl%dfl(:,firstsubst,secondsubst))
!!        !    gl%dfl_int(:,secondsubst,firstsubst) = int(gl%dfl(:,firstsubst,secondsubst))
!!        !endif
!!
!!
!!
!!
!!        if (dfnspec > 0) then !for special terms
!!            do j=1,dfnspec
!!                read(mix_unit,*) dfni,dfti,dfdi,dfetai,dfepsi,dfbetai,dfgammai
!!
!!                !---------------------------------------
!!                !Andreas Aug 2010
!!                gl%dfn(j+dfnreg,firstsubst,secondsubst) = dfni
!!                gl%dfn(j+dfnreg,secondsubst,firstsubst) = dfni
!!                gl%dfd(j+dfnreg,firstsubst,secondsubst) = dfdi
!!                gl%dfd(j+dfnreg,secondsubst,firstsubst) = dfdi
!!                gl%dft(j+dfnreg,firstsubst,secondsubst) = dfti
!!                gl%dft(j+dfnreg,secondsubst,firstsubst) = dfti
!!
!!                gl%dfeta(j,firstsubst,secondsubst) = dfetai
!!                gl%dfeta(j,secondsubst,firstsubst) = dfetai
!!                gl%dfeps(j,firstsubst,secondsubst) = dfepsi
!!                gl%dfeps(j,secondsubst,firstsubst) = dfepsi
!!                gl%dfbeta(j,firstsubst,secondsubst) = dfbetai
!!                gl%dfbeta(j,secondsubst,firstsubst) = dfbetai
!!                gl%dfgamma(j,firstsubst,secondsubst) = dfgammai
!!                gl%dfgamma(j,secondsubst,firstsubst) = dfgammai
!!                !----------------------------------------
!!
!!            enddo
!!        endif
!!
!!        if (dfngauss > 0) then !for gaussian terms
!!            do j=1,dfngauss
!!                read(mix_unit,*) dfni,dfti,dfdi,dfpli,dftli,dfgetai,dfgbetai,dfggami,dfgepsi,dfgdummy1,dfgdummy2,dfgdummy3
!!
!!                !---------------------------------------
!!                !Andreas Aug 2010
!!                gl%dfn(j+dfnreg+dfnspec,firstsubst,secondsubst) = dfni
!!                gl%dfn(j+dfnreg+dfnspec,secondsubst,firstsubst) = dfni
!!                gl%dfd(j+dfnreg+dfnspec,firstsubst,secondsubst) = dfdi
!!                gl%dfd(j+dfnreg+dfnspec,secondsubst,firstsubst) = dfdi
!!                gl%dft(j+dfnreg+dfnspec,firstsubst,secondsubst) = dfti
!!                gl%dft(j+dfnreg+dfnspec,secondsubst,firstsubst) = dfti
!!
!!                gl%dfgeta(j,firstsubst,secondsubst) = dfgetai
!!                gl%dfgeta(j,secondsubst,firstsubst) = dfgetai
!!                gl%dfgeps(j,firstsubst,secondsubst) = dfgepsi
!!                gl%dfgeps(j,secondsubst,firstsubst) = dfgepsi
!!                gl%dfgbeta(j,firstsubst,secondsubst) = dfgbetai
!!                gl%dfgbeta(j,secondsubst,firstsubst) = dfgbetai
!!                gl%dfggam(j,firstsubst,secondsubst) = dfggami
!!                gl%dfggam(j,secondsubst,firstsubst) = dfggami
!!                !----------------------------------------
!!
!!            enddo
!!        endif
!!    else
!!        internalioerror = -7881
!!        return
!!        !        if (firstsubst == substance1) then
!!        !            write(*,*) 'Error: occured when opening the mix-file: ', trim(filename1)
!!        !        else
!!        !            write(*,*) 'Error: occured when opening the mix-file: ', trim(filename2)
!!        !        endif
!!    endif
!!
!!    close(mix_unit)
!!    return
!!
!!    !only enter if end of file is reached during read:
!!999 internalioerror=-7882
!!    !write(*,*) 'Error: end of mix-file is reached during read'
!!
!!    end subroutine read_mix_file
!!
!!    !    !called by sub "read_from_fluidfiles", may be called from somewhere else too:
!!    !    subroutine read_atcoeff_file(gl,ufcfilepath,internalioerror)
!!    !
!!    !
!!    !
!!    !    !Declarations: -----------------------------------------------------------------------------------------
!!    !
!!    !    implicit none
!!    !
!!    !    type(type_gl) :: gl
!!    !
!!    !
!!    !    !Define Input and Output of the Subroutine:
!!    !    character(len=*) :: ufcfilepath      !has to be transfered to this subroutine
!!    !    integer :: internalioerror                 !no error as startvalue
!!    !
!!    !    !to count:
!!    !    integer :: k,l,nrcoeffs
!!    !    character(len=1) :: dummy
!!    !
!!    !    !unit integers
!!    !    integer:: atcoeff_unit
!!    !    dummy = ''
!!    !
!!    !    !pre-set:
!!    !    internalioerror = 0
!!    !    atcoeff_unit = 0
!!    !
!!    !    !open file and read coefficients:
!!    !    open(newunit=atcoeff_unit, file= ufcfilepath // 'atcoeff.txt', status='old', action='read', iostat=internalioerror)
!!    !    if (internalioerror == 0) then
!!    !        read(atcoeff_unit,*,end=9000) nrcoeffs
!!    !        read(atcoeff_unit,*,end=9000) dummy
!!    !        kloop:  do k = 1, nrcoeffs
!!    !            !backspace(atcoeff_unit) !if less than 100 elements in txt-line, then last read command ends in jump into next line, so backspace command necessary
!!    !            read(atcoeff_unit,*,end=9000) (gl%atcoeff(k,l),l=1,nrcoeffs)
!!    !        enddo kloop
!!    !    else
!!    !        internalioerror = -7883
!!    !        !write(*,*) 'Error: occured when opening the atcoeff-file: ', ufcfilepath , 'atcoeff.txt'
!!    !    endif
!!    !
!!    !9000 continue
!!    !
!!    !
!!    !    end subroutine read_atcoeff_file
!!    !
!!    !    !called for psrk callclations:
!!    !    subroutine read_btcoeff_file(gl,ufcfilepath,internalioerror)
!!    !
!!    !
!!    !
!!    !    !Declarations: -----------------------------------------------------------------------------------------
!!    !
!!    !    implicit none
!!    !
!!    !    type(type_gl) :: gl
!!    !
!!    !
!!    !    !Define Input and Output of the Subroutine:
!!    !    character(len=*) :: ufcfilepath      !has to be transfered to this subroutine
!!    !    integer :: internalioerror                 !no error as startvalue
!!    !
!!    !    !to count:
!!    !    integer :: k,l,nrcoeffs
!!    !    character(len=1) :: dummy
!!    !
!!    !    integer:: btcoeff_unit
!!    !
!!    !    !pre-set:
!!    !    dummy = ''
!!    !    internalioerror = 0
!!    !    btcoeff_unit = 0
!!    !
!!    !    !open file and read coefficients:
!!    !    open(newunit=btcoeff_unit, file= ufcfilepath // 'btcoeff.txt', status='old', action='read', iostat=internalioerror)
!!    !    if (internalioerror == 0) then
!!    !        read(btcoeff_unit,*,end=9000) nrcoeffs
!!    !        read(btcoeff_unit,*,end=9000) dummy
!!    !        kloop:  do k = 1, nrcoeffs
!!    !            read(btcoeff_unit,*,end=9000) (gl%btcoeff(k,l),l=1,nrcoeffs)
!!    !        enddo kloop
!!    !        !else
!!    !        !write(*,*) 'Error: occured when opening the btcoeff-file: ', ufcfilepath , 'btcoeff.txt'
!!    !    endif
!!    !
!!    !9000 continue
!!    !
!!    !
!!    !    end subroutine read_btcoeff_file
!!    !
!!    !    subroutine read_ctcoeff_file(gl,ufcfilepath,internalioerror)
!!    !
!!    !
!!    !
!!    !    !Declarations: -----------------------------------------------------------------------------------------
!!    !
!!    !    implicit none
!!    !
!!    !    type(type_gl) :: gl
!!    !
!!    !
!!    !    !Define Input and Output of the Subroutine:
!!    !    character(len=*) :: ufcfilepath      !has to be transfered to this subroutine
!!    !    integer :: internalioerror                 !no error as startvalue
!!    !
!!    !    !to count:
!!    !    integer :: k,l,nrcoeffs
!!    !    character(len=1) :: dummy
!!    !
!!    !    integer:: ctcoeff_unit
!!    !
!!    !    !pre-set:
!!    !    dummy = ''
!!    !    internalioerror = 0
!!    !    ctcoeff_unit = 0
!!    !
!!    !    !open file and read coefficients:
!!    !    open(newunit=ctcoeff_unit, file= ufcfilepath // 'ctcoeff.txt', status='old', action='read', iostat=internalioerror)
!!    !    if (internalioerror == 0) then
!!    !        read(ctcoeff_unit,*,end=9000) nrcoeffs
!!    !        read(ctcoeff_unit,*,end=9000) dummy
!!    !        kloop:  do k = 1, nrcoeffs
!!    !            read(ctcoeff_unit,*,end=9000) (gl%ctcoeff(k,l),l=1,nrcoeffs)
!!    !        enddo kloop
!!    !        !else
!!    !        !write(*,*) 'Error: occured when opening the ctcoeff-file: ', ufcfilepath , 'ctcoeff.txt'
!!    !    endif
!!    !
!!    !9000 continue
!!    !
!!    !
!!    !    end subroutine read_ctcoeff_file
!!
!!
!!
!!
!!    subroutine bwr_recalc(gl,nrsubst)
!!
!!
!!    implicit none
!!
!!    type(type_gl) :: gl
!!
!!
!!
!!    integer :: i, l, count, nrsubst, pol, pol_add,loop
!!
!!
!!    count = 19           !for MBWR exp-coeffs 20-40
!!    pol = 0              !number of polynomial terms
!!    pol_add = 0          !number of additional polynomial terms
!!
!!
!!    !only pol di
!!    if (gl%eos_coeff%nreg(nrsubst) == 32) then            !MBWR
!!        do l=1, 5
!!            gl%eos_coeff%di(l,nrsubst) = 1.d0
!!        enddo
!!        do l= 6, 9
!!            gl%eos_coeff%di(l,nrsubst) = 2.d0
!!        enddo
!!        do l=10, 12
!!            gl%eos_coeff%di(l,nrsubst) = 3.d0
!!        enddo
!!        gl%eos_coeff%di(13,nrsubst) = 4.d0
!!        do l=14, 15
!!            gl%eos_coeff%di(l,nrsubst) = 5.d0
!!        enddo
!!        gl%eos_coeff%di(16,nrsubst) = 6.d0
!!        do l=17, 18
!!            gl%eos_coeff%di(l,nrsubst) = 7.d0
!!        enddo
!!        gl%eos_coeff%di(19,nrsubst) = 8.d0
!!
!!        do l=1, 6
!!            gl%eos_coeff%di(l+count,nrsubst) = 2.d0*l
!!            gl%eos_coeff%di(l+count+1,nrsubst) = 2.d0*l
!!            count=count+1
!!        enddo
!!        gl%eos_coeff%di(32,nrsubst) = 12.d0
!!
!!
!!        !only pol ti
!!        gl%eos_coeff%ti(1,nrsubst) = 0.d0
!!        gl%eos_coeff%ti(2,nrsubst) = 0.5d0
!!        gl%eos_coeff%ti(3,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(4,nrsubst) = 2.d0
!!        gl%eos_coeff%ti(5,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(6,nrsubst) = 0.d0
!!        gl%eos_coeff%ti(7,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(8,nrsubst) = 2.d0
!!        gl%eos_coeff%ti(9,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(10,nrsubst) = 0.d0
!!        gl%eos_coeff%ti(11,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(12,nrsubst) = 2.d0
!!        gl%eos_coeff%ti(13,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(14,nrsubst) = 2.d0
!!        gl%eos_coeff%ti(15,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(16,nrsubst) = 2.d0
!!        gl%eos_coeff%ti(17,nrsubst) = 2.d0
!!        gl%eos_coeff%ti(18,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(19,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(20,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(21,nrsubst) = 4.d0
!!        gl%eos_coeff%ti(22,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(23,nrsubst) = 5.d0
!!        gl%eos_coeff%ti(24,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(25,nrsubst) = 4.d0
!!        gl%eos_coeff%ti(26,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(27,nrsubst) = 5.d0
!!        gl%eos_coeff%ti(28,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(29,nrsubst) = 4.d0
!!        gl%eos_coeff%ti(30,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(31,nrsubst) = 4.d0
!!        gl%eos_coeff%ti(32,nrsubst) = 5.d0
!!
!!    elseif (gl%eos_coeff%nreg(nrsubst) == 19) then        !Bender
!!        do l=1, 5
!!            gl%eos_coeff%di(l,nrsubst) = 1.d0
!!        enddo
!!        do l=6, 8
!!            gl%eos_coeff%di(l,nrsubst) = 2.d0
!!        enddo
!!        do l=9, 10
!!            gl%eos_coeff%di(l,nrsubst) = 3.d0
!!        enddo
!!        do l=11, 12
!!            gl%eos_coeff%di(l,nrsubst) = 4.d0
!!        enddo
!!        gl%eos_coeff%di(13,nrsubst) = 5.d0
!!        do l=14, 16
!!            gl%eos_coeff%di(l,nrsubst) = 2.d0
!!        enddo
!!        do l=17,19
!!            gl%eos_coeff%di(l,nrsubst) = 4.d0
!!        enddo
!!
!!
!!        gl%eos_coeff%ti(1,nrsubst) = 0.d0
!!        gl%eos_coeff%ti(2,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(3,nrsubst) = 2.d0
!!        gl%eos_coeff%ti(4,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(5,nrsubst) = 4.d0
!!        gl%eos_coeff%ti(6,nrsubst) = 0.d0
!!        gl%eos_coeff%ti(7,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(8,nrsubst) = 2.d0
!!        gl%eos_coeff%ti(9,nrsubst) = 0.d0
!!        gl%eos_coeff%ti(10,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(11,nrsubst) = 0.d0
!!        gl%eos_coeff%ti(12,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(13,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(14,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(15,nrsubst) = 4.d0
!!        gl%eos_coeff%ti(16,nrsubst) = 5.d0
!!        gl%eos_coeff%ti(17,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(18,nrsubst) = 4.d0
!!        gl%eos_coeff%ti(19,nrsubst) = 5.d0
!!
!!    elseif (gl%eos_coeff%nreg(nrsubst) == 12) then        !Starling
!!        do l=1, 5
!!            gl%eos_coeff%di(l,nrsubst) = 1.d0
!!        enddo
!!        do l=6, 8
!!            gl%eos_coeff%di(l,nrsubst) = 2.d0
!!        enddo
!!        do l=9, 10
!!            gl%eos_coeff%di(l,nrsubst) = 5.d0
!!        enddo
!!        gl%eos_coeff%di(11,nrsubst) = 2.d0
!!        gl%eos_coeff%di(12,nrsubst) = 4.d0
!!
!!
!!        gl%eos_coeff%ti(1,nrsubst) = 0.d0
!!        gl%eos_coeff%ti(2,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(3,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(4,nrsubst) = 4.d0
!!        gl%eos_coeff%ti(5,nrsubst) = 5.d0
!!        gl%eos_coeff%ti(6,nrsubst) = 0.d0
!!        gl%eos_coeff%ti(7,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(8,nrsubst) = 2.d0
!!        gl%eos_coeff%ti(9,nrsubst) = 1.d0
!!        gl%eos_coeff%ti(10,nrsubst) = 2.d0
!!        gl%eos_coeff%ti(11,nrsubst) = 3.d0
!!        gl%eos_coeff%ti(12,nrsubst) = 3.d0
!!
!!    endif
!!
!!
!!    !n2i to ni  (Multiparameter Equations of State (R.Span, 2000) page 26 Eq[3.28] to Eq[3.26])
!!    if (((gl%Req(nrsubst) /= 1.d0) .or. ((gl%Req(nrsubst) == 1.d0) .and. gl%vir)) .and. (.not.gl%sta_read) .and. (.not.gl%ben_read)) then !(gl%Req(nrsubst) == 1.d0) .and. gl%vir) Exception for mbwr equations in reduced units
!!        do i=1, gl%eos_coeff%nreg(nrsubst)                                      !
!!            gl%ncoeff(nrsubst,i) = gl%ncoeff(nrsubst,i) * ((gl%rhored(nrsubst)/gl%Factor)**gl%eos_coeff%di(i,nrsubst)) / &
!!                (gl%tred(nrsubst)**(gl%eos_coeff%ti(i,nrsubst))) /(gl%Req(nrsubst)/gl%factorrbwr)
!!        enddo
!!    elseif ((gl%sta_read) .or. (gl%ben_read)) then
!!        do i=1, gl%eos_coeff%nreg(nrsubst)
!!            gl%ncoeff(nrsubst,i) = gl%ncoeff(nrsubst,i) * ((gl%rhored(nrsubst))**gl%eos_coeff%di(i,nrsubst)) / &
!!                (gl%tred(nrsubst)**(gl%eos_coeff%ti(i,nrsubst))) /gl%Req(nrsubst)
!!        enddo
!!    endif
!!
!!    !term number after integration
!!    if (gl%bwr_read) then
!!        gl%eos_coeff%nreg(nrsubst) = 40
!!        pol = 22
!!        pol_add = 3
!!    elseif (gl%ben_read) then
!!        gl%eos_coeff%nreg(nrsubst) = 22
!!        pol = 16
!!        pol_add = 3
!!    elseif (gl%sta_read) then
!!        gl%eos_coeff%nreg(nrsubst) = 13
!!        pol = 11
!!        pol_add = 1
!!    endif
!!
!!
!!    do l=1, pol
!!        gl%eos_coeff%p_i(l,nrsubst) = 0
!!    enddo
!!    do l=(pol+1), gl%eos_coeff%nreg(nrsubst)
!!        gl%eos_coeff%p_i(l,nrsubst) = 2
!!    enddo
!!
!!
!!    !poly-term coefficients transformation
!!    Do i=1, (pol-pol_add)
!!        gl%eos_coeff%ni(i,nrsubst) = gl%ncoeff(nrsubst,i)/gl%eos_coeff%di(i,nrsubst)
!!    enddo
!!
!!
!!    !gamma calculation
!!    do i=1, (pol-pol_add)
!!        gl%eos_coeff%gama(i,nrsubst)=  0.d0
!!    enddo
!!    do i=(pol-pol_add+1), gl%eos_coeff%nreg(nrsubst)
!!        gl%eos_coeff%gama(i,nrsubst)= ((gl%rhoc(nrsubst)/gl%Factor)/gl%gama_bwr(nrsubst))**2.d0
!!    enddo
!!
!!
!!
!!    !integrated ni-terms      from Multiparameter Equations of State(R.Span) page 28 / Table 3.4 and 3.5
!!    !MBWR terms
!!    if (gl%bwr_read) then
!!        gl%eos_coeff%ni(20,nrsubst) = gl%ncoeff(nrsubst,20)/(2.d0*gl%eos_coeff%gama(20,nrsubst)) + gl%ncoeff(nrsubst,22)/(2.d0*gl%eos_coeff%gama(20,nrsubst)**2) &
!!            + gl%ncoeff(nrsubst,24)/gl%eos_coeff%gama(20,nrsubst)**3 + 3.d0*gl%ncoeff(nrsubst,26)/gl%eos_coeff%gama(20,nrsubst)**4 &
!!            + 12.d0*gl%ncoeff(nrsubst,28)/gl%eos_coeff%gama(20,nrsubst)**5 + 60.d0*gl%ncoeff(nrsubst,30)/gl%eos_coeff%gama(20,nrsubst)**6
!!        gl%eos_coeff%ni(21,nrsubst) = gl%ncoeff(nrsubst,21)/(2.d0*gl%eos_coeff%gama(21,nrsubst)) + gl%ncoeff(nrsubst,25)/(gl%eos_coeff%gama(21,nrsubst)**3) &
!!            + 12.d0*gl%ncoeff(nrsubst,29)/(gl%eos_coeff%gama(21,nrsubst)**5) + 60.d0*gl%ncoeff(nrsubst,31)/gl%eos_coeff%gama(21,nrsubst)**6
!!        gl%eos_coeff%ni(22,nrsubst) = (gl%ncoeff(nrsubst,23)/(2.d0*gl%eos_coeff%gama(22,nrsubst)**2)) + (3.d0*gl%ncoeff(nrsubst,27)/(gl%eos_coeff%gama(22,nrsubst)**4)) &
!!            + (60.d0*gl%ncoeff(nrsubst,32)/(gl%eos_coeff%gama(22,nrsubst)**6))
!!        gl%eos_coeff%ni(23,nrsubst) = -gl%eos_coeff%ni(20,nrsubst)
!!        gl%eos_coeff%ni(24,nrsubst) = -gl%eos_coeff%ni(21,nrsubst)
!!        gl%eos_coeff%ni(25,nrsubst) = -gl%eos_coeff%ni(22,nrsubst)
!!        gl%eos_coeff%ni(26,nrsubst) = -gl%ncoeff(nrsubst,22)/(2.d0*gl%eos_coeff%gama(26,nrsubst)) - gl%ncoeff(nrsubst,24)/(gl%eos_coeff%gama(26,nrsubst)**2) &
!!            - 3.d0*gl%ncoeff(nrsubst,26)/(gl%eos_coeff%gama(26,nrsubst)**3) - 12.d0*gl%ncoeff(nrsubst,28)/(gl%eos_coeff%gama(26,nrsubst)**4) &
!!            - 60.d0*gl%ncoeff(nrsubst,30)/(gl%eos_coeff%gama(26,nrsubst)**5)
!!        gl%eos_coeff%ni(27,nrsubst) = -gl%ncoeff(nrsubst,25)/(gl%eos_coeff%gama(27,nrsubst)**2) - 12.d0*gl%ncoeff(nrsubst,29)/(gl%eos_coeff%gama(27,nrsubst)**4) &
!!            - 60.d0*gl%ncoeff(nrsubst,31)/(gl%eos_coeff%gama(27,nrsubst)**5)
!!        gl%eos_coeff%ni(28,nrsubst) = -gl%ncoeff(nrsubst,23)/(2.d0*gl%eos_coeff%gama(28,nrsubst)) - 3.d0*gl%ncoeff(nrsubst,27)/(gl%eos_coeff%gama(28,nrsubst)**3) &
!!            - 60.d0*gl%ncoeff(nrsubst,32)/(gl%eos_coeff%gama(28,nrsubst)**5)
!!        gl%eos_coeff%ni(29,nrsubst) = -gl%ncoeff(nrsubst,24)/(2.d0*gl%eos_coeff%gama(29,nrsubst)) - 3.d0*gl%ncoeff(nrsubst,26)/(2.d0*gl%eos_coeff%gama(29,nrsubst)**2) &
!!            - 6.d0*gl%ncoeff(nrsubst,28)/(gl%eos_coeff%gama(29,nrsubst)**3) - 30.d0*gl%ncoeff(nrsubst,30)/(gl%eos_coeff%gama(29,nrsubst)**4)
!!        gl%eos_coeff%ni(30,nrsubst) = -gl%ncoeff(nrsubst,25)/(2.d0*gl%eos_coeff%gama(30,nrsubst)) - 6.d0*gl%ncoeff(nrsubst,29)/(gl%eos_coeff%gama(30,nrsubst)**3) &
!!            - 30.d0*gl%ncoeff(nrsubst,31)/(gl%eos_coeff%gama(30,nrsubst)**4)
!!        gl%eos_coeff%ni(31,nrsubst) = -3.d0*gl%ncoeff(nrsubst,27)/(2.d0*gl%eos_coeff%gama(31,nrsubst)**2) - 30.d0*gl%ncoeff(nrsubst,32)/(gl%eos_coeff%gama(31,nrsubst)**4)
!!        gl%eos_coeff%ni(32,nrsubst) = -gl%ncoeff(nrsubst,26)/(2.d0*gl%eos_coeff%gama(32,nrsubst)) - 2.d0*gl%ncoeff(nrsubst,28)/(gl%eos_coeff%gama(32,nrsubst)**2) &
!!            - 10.d0*gl%ncoeff(nrsubst,30)/(gl%eos_coeff%gama(32,nrsubst)**3)
!!        gl%eos_coeff%ni(33,nrsubst) = -2.d0*gl%ncoeff(nrsubst,29)/(gl%eos_coeff%gama(33,nrsubst)**2) - 10.d0*gl%ncoeff(nrsubst,31)/(gl%eos_coeff%gama(33,nrsubst)**3)
!!        gl%eos_coeff%ni(34,nrsubst) = -gl%ncoeff(nrsubst,27)/(2.d0*gl%eos_coeff%gama(34,nrsubst)) - 10.d0*gl%ncoeff(nrsubst,32)/(gl%eos_coeff%gama(34,nrsubst)**3)
!!        gl%eos_coeff%ni(35,nrsubst) = -gl%ncoeff(nrsubst,28)/(2.d0*gl%eos_coeff%gama(35,nrsubst)) - 5.d0*gl%ncoeff(nrsubst,30)/(2.d0*(gl%eos_coeff%gama(35,nrsubst)**2))
!!        gl%eos_coeff%ni(36,nrsubst) = -gl%ncoeff(nrsubst,29)/(2.d0*gl%eos_coeff%gama(36,nrsubst)) - 5.d0*gl%ncoeff(nrsubst,31)/(2.d0*(gl%eos_coeff%gama(36,nrsubst)**2))
!!        gl%eos_coeff%ni(37,nrsubst) = -5.d0*gl%ncoeff(nrsubst,32)/(2.d0*(gl%eos_coeff%gama(37,nrsubst)**2))
!!        gl%eos_coeff%ni(38,nrsubst) = -gl%ncoeff(nrsubst,30)/(2.d0*gl%eos_coeff%gama(38,nrsubst))
!!        gl%eos_coeff%ni(39,nrsubst) = -gl%ncoeff(nrsubst,31)/(2.d0*gl%eos_coeff%gama(39,nrsubst))
!!        gl%eos_coeff%ni(40,nrsubst) = -gl%ncoeff(nrsubst,32)/(2.d0*gl%eos_coeff%gama(40,nrsubst))
!!
!!        !Bender terms
!!    elseif (gl%ben_read) then
!!        gl%eos_coeff%ni(14,nrsubst) = gl%ncoeff(nrsubst,14)/(2.d0*gl%eos_coeff%gama(14,nrsubst)) + gl%ncoeff(nrsubst,17)/(2.d0*(gl%eos_coeff%gama(14,nrsubst)**2))
!!        gl%eos_coeff%ni(15,nrsubst) = gl%ncoeff(nrsubst,15)/(2.d0*gl%eos_coeff%gama(15,nrsubst)) + gl%ncoeff(nrsubst,18)/(2.d0*(gl%eos_coeff%gama(15,nrsubst)**2))
!!        gl%eos_coeff%ni(16,nrsubst) = gl%ncoeff(nrsubst,16)/(2.d0*gl%eos_coeff%gama(16,nrsubst)) + gl%ncoeff(nrsubst,19)/(2.d0*(gl%eos_coeff%gama(16,nrsubst)**2))
!!        gl%eos_coeff%ni(17,nrsubst) = -gl%eos_coeff%ni(14,nrsubst)
!!        gl%eos_coeff%ni(18,nrsubst) = -gl%eos_coeff%ni(15,nrsubst)
!!        gl%eos_coeff%ni(19,nrsubst) = -gl%eos_coeff%ni(16,nrsubst)
!!        gl%eos_coeff%ni(20,nrsubst) = -gl%ncoeff(nrsubst,17)/(2.d0*gl%eos_coeff%gama(20,nrsubst))
!!        gl%eos_coeff%ni(21,nrsubst) = -gl%ncoeff(nrsubst,18)/(2.d0*gl%eos_coeff%gama(21,nrsubst))
!!        gl%eos_coeff%ni(22,nrsubst) = -gl%ncoeff(nrsubst,19)/(2.d0*gl%eos_coeff%gama(22,nrsubst))
!!
!!        !Starling terms
!!    elseif (gl%sta_read) then
!!        gl%eos_coeff%ni(11,nrsubst) = gl%ncoeff(nrsubst,11)/(2.d0*gl%eos_coeff%gama(11,nrsubst)) + gl%ncoeff(nrsubst,12)/(2.d0*(gl%eos_coeff%gama(11,nrsubst)**2))
!!        gl%eos_coeff%ni(12,nrsubst) = -gl%eos_coeff%ni(11,nrsubst)
!!        gl%eos_coeff%ni(13,nrsubst) = -gl%ncoeff(nrsubst,12)/(2.d0*gl%eos_coeff%gama(13,nrsubst))
!!    endif
!!
!!
!!    !integrated di coeff
!!    if (.not.gl%sta_read) then
!!        do l=(pol-2), (pol+3)
!!            gl%eos_coeff%di(l,nrsubst) = 0.d0
!!        enddo
!!        do l=(pol+4), (pol+6)
!!            gl%eos_coeff%di(l,nrsubst) = 2.d0
!!        enddo
!!        do l=29, 31
!!            gl%eos_coeff%di(l,nrsubst) = 4.d0
!!        enddo
!!
!!        do l=32, 34
!!            gl%eos_coeff%di(l,nrsubst) = 6.d0
!!        enddo
!!
!!        do l=35, 37
!!            gl%eos_coeff%di(l,nrsubst) = 8.d0
!!        enddo
!!        do l=38, 40
!!            gl%eos_coeff%di(l,nrsubst) = 10.d0
!!        enddo
!!    else
!!        gl%eos_coeff%di(11,nrsubst) = 0.d0
!!        gl%eos_coeff%di(12,nrsubst) = 0.d0
!!        gl%eos_coeff%di(13,nrsubst) = 2.d0
!!    endif
!!
!!
!!    !integrated ti coeff
!!    if (.not.gl%sta_read) then
!!        if (gl%bwr_read) then
!!            loop = 7
!!            count = pol-pol_add
!!        elseif (gl%ben_read) then
!!            loop = 3
!!            count = pol
!!        endif
!!        do l=1, loop
!!            gl%eos_coeff%ti((l+count),nrsubst) = 3.d0
!!            gl%eos_coeff%ti((l+count+1),nrsubst) = 4.d0
!!            gl%eos_coeff%ti((l+count+2),nrsubst) = 5.d0
!!            count = count+2
!!        enddo
!!    else
!!        do l=11, 13
!!            gl%eos_coeff%ti(l,nrsubst) = 3.d0
!!        enddo
!!    endif
!!
!!
!!
!!    end subroutine
!!
!!
!!
!!
!!
!!
!!
!!    end module sub_file_input_module
