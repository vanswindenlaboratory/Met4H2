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

    ! module for file initialize.f90
    module initialize_module
    !global use inclusion
    use module_all_types


    contains


    subroutine initialize(gl)
    !!DEC$ ATTRIBUTES DLLEXPORT :: initialize
    !Set default values for the arrays



























    implicit none

    type(type_gl) :: gl


    !module_ref
    gl%calc_ref = .true.
    gl%ref_set = .false.
    gl%VLE_needed = .true.
    !gl%transport = 0

    !module_fluid_parameters

    !gl%components_old = ""
    !gl%ncomp_old = 0
    gl%same_components = .false.
    !gl%path_old = ""

    !gl%components = ""
    gl%molfractions = 0.d0
    !gl%ncomp = 0

    !gl%substfullname = ''
    !gl%substshortname = ''
    !gl%substchemname = ''
    !gl%substsynonym = ''
    !gl%substcasnr = ''

    gl%tc = 0.d0
    gl%pc = 0.d0
    gl%rhoc = 0.d0

    gl%ttp = 0.d0
    gl%ptp = 0.d0

    gl%ptpmod = 0.d0
    gl%rhotp = 0.d0

    gl%tnbp = 0.d0

    gl%wm = 0.d0
    gl%accen = 0.d0
    gl%dipole = 0.d0

    gl%tminfluid = 0.d0
    gl%tmaxfluid = 0.d0
    gl%pmaxfluid = 0.d0
    gl%rhomaxfluid = 0.d0

    gl%Req = 0.d0
    !SH: 01/16 commented -> LJF was broken
    !Factor = 1.D3
    !factorpress=1.D0
    !factortrans=1.D6
    !factorrbwr=1.D2

    !inptorig='' SH: 11/2015 because the input is now checked before the initialization precedure, there is "theoretecially" now need to init this var

    gl%already_converted = 0
    
    gl%pvexist = .false.
    gl%meltexist = .false.

    !gl%Eq_type = 1
    !gl%Mix_type = 1

    gl%twophasecalc = .false.
    gl%bad_mixmodel = .false.

    gl%factor_decide4vap = 0.4d0

    !hold_limits = .true.

    gl%phase_id = 0

    gl%A_cv0 = 0.d0
    gl%B_cv0 = 0.d0
    gl%C_cv0 = 0.d0
    gl%D_cv0 = 0.d0
    gl%E_cv0 = 0.d0
    gl%F_cv0 = 0.d0
    gl%G_cv0 = 0.d0

    gl%cp0_A = 0.d0
    gl%cp0_B = 0.d0
    gl%cp0_C = 0.d0
    gl%cp0_D = 0.d0

    !module_ancillary_parameters
    gl%vpred = 0.d0
    gl%tvpred = 0.d0
    gl%dlred = 0.d0
    gl%tdlred = 0.d0
    gl%dvred = 0.d0
    gl%tdvred = 0.d0
    gl%pmeltred = 0.d0
    gl%tpmeltred = 0.d0
    gl%psubred = 0.d0
    gl%tpsubred = 0.d0
    gl%pmeltmintemp = 0.d0
    gl%pmeltmaxtemp = 0.d0
    gl%tsatmelt_old = 0.d0
    gl%pmelt_high = 0.d0
    gl%melt_it=.false.
    gl%p_fake = .false.
    gl%melt_p = 0.d0

    gl%nvpcoeff = 0
    gl%ndlcoeff = 0
    gl%ndvcoeff = 0
    gl%npmeltcoeff1 = 0
    gl%npmeltcoeff2 = 0
    gl%npmeltcoeff3 = 0
    gl%npsubcoeff1 = 0
    gl%npsubcoeff2 = 0
    gl%npsubcoeff3 = 0
    gl%vptype = 0
    gl%dltype = 0
    gl%dvtype = 0
    gl%pmelttype = 0
    gl%psubtype = 0

    gl%vpcoeff = 0.d0
    gl%vpexp = 0.d0
    gl%dlcoeff = 0.d0
    gl%dlexp = 0.d0
    gl%dvcoeff = 0.d0
    gl%dvexp = 0.d0
    gl%pmeltcoeff = 0.d0
    gl%pmeltexp = 0.d0
    gl%psubcoeff = 0.d0
    gl%psubexp = 0.d0

    call initialize_general_eos_parameters(gl)

    if(allocated(gl%eos_coeff)) call initialize_eos_coefficients(gl)

    call initialize_ideal_gas_coefficients(gl)

    !module_mixture_parameters

    gl%rfbetat = 0.d0
    gl%rfbetarho = 0.d0
    gl%rfgammat = 0.d0
    gl%rfgammarho = 0.d0

    gl%dfn = 0.d0
    gl%dfd = 0.d0
    gl%dft = 0.d0
    gl%dfl = 0.d0
    gl%dfp = 0.d0
    gl%dfeta = 0.d0
    gl%dfeps = 0.d0
    gl%dfbeta = 0.d0
    gl%dfgamma = 0.d0
    gl%dfgeps = 0.d0
    gl%dfgbeta = 0.d0
    gl%dfggam = 0.d0
    gl%dfgeta = 0.d0

    gl%Fij = 0.d0
    gl%dfpol = 0
    gl%dfexp = 0
    gl%dfgau = 0
    !module_psrk
    gl%psrk_ok = .false.
    !gl%ufc = 0.d0
    gl%ccoeff = 0.d0
    !gl%atcoeff = 0.d0
    !gl%btcoeff = 0.d0
    !gl%ctcoeff = 0.d0
    gl%ngroup = 0
    gl%ngroup_nsubgroup = 0

    !module_asso
    gl%assoc = 0.d0
    gl%assovolred = 0.d0
    gl%assovolint = 0.d0
    gl%assoenergy = 0.d0
    gl%assoexist = .false.

    !module_VLE
    call initialize_vle(gl)

    !module_phasedet_pure
    gl%rho_it = 0.d0
    gl%p_it = 0.d0

    if (allocated(gl%uncty))call initialize_uncty(gl)

    !Andreas Oct 2012
    !module_cubic
    gl%kij_SRK = 0.D0
    gl%lij_SRK = 0.D0
    gl%rhored_SRK = 1.D0   !mol/m³
    gl%pred_SRK = 1.D0     !Pa
    gl%Tred_SRK = 1.D0     !K
    gl%Cji_cubic = 0.D0
    gl%R_SRK = 8.3144598D0  !J / mol K     !Andreas Nov 2015, new value of gas constant
    !TEST_AJ - For paper Bell & Jäger 2016 values for SRK similar to Stradi et al. have been chosen August 2016
    !R_SRK = 8.31451D0

    !Stefan Feb 2014
    !module_cubic
    gl%kij_PR = 0.D0
    gl%lij_PR = 0.D0
    gl%rhored_PR = 1.D0   !mol/m³
    gl%pred_PR = 1.D0     !Pa
    gl%Tred_PR = 1.D0     !K
    gl%R_PR = 8.3144598D0  !J / mol K     !Andreas Nov 2015, new value of gas constant

    !Stefan Feb 2014
    !module LKP
    gl%kij_LKP = 1.d0
    gl%tc_ij_LKP = 0.d0
    gl%vc_ij_LKP = 0.d0
    gl%accenLKPMix = 0.d0
    gl%zcLKP = 0.d0
    gl%wLKP = 0.d0

    gl%lkp_b1_0 = 0.1181193d0
    gl%lkp_b2_0 = 0.265728d0
    gl%lkp_b3_0 = 0.15479d0
    gl%lkp_b4_0 = 0.030323d0
    gl%lkp_c1_0 = 0.0236744d0
    gl%lkp_c2_0 = 0.0186984d0
    gl%lkp_c3_0 = 0.d0
    gl%lkp_c4_0 = 0.042724d0
    gl%lkp_d1_0 = 0.155428d-4
    gl%lkp_d2_0 = 0.623689d-4
    gl%lkp_beta_0 = 0.65392d0
    gl%lkp_gamma_0 = 0.060167d0
    gl%lkp_w_ref = 0.3978d0
    gl%lkp_b1_ref = 0.2026579d0
    gl%lkp_b2_ref = 0.331511d0
    gl%lkp_b3_ref = 0.027655d0
    gl%lkp_b4_ref = 0.203488d0
    gl%lkp_c1_ref = 0.0313385d0
    gl%lkp_c2_ref = 0.0503618d0
    gl%lkp_c3_ref = 0.016901d0
    gl%lkp_c4_ref = 0.041577d0
    gl%lkp_d1_ref = 0.48736d-4
    gl%lkp_d2_ref = 0.0740336d-4
    gl%lkp_beta_ref = 1.226d0
    gl%lkp_gamma_ref = 0.03754d0
    if (gl%check_solid .eqv. .true.) then
        call initialize_solids(gl)
    endif

    if(allocated(gl%de)) call initialize_de(gl)

    !module_surface_tension
    gl%stn_read = .false.
    !gl%stmodel = ' '
    !gl%sigma0_st = 0.d0
    !gl%n_exp_st = 0.d0
    !gl%low_temp_st = 0.d0
    !gl%upp_temp_st = 0.d0
    !gl%upp_press_st = 0.d0
    !gl%max_dens_st = 0.d0
    !gl%term_num_st = 0
    !gl%Temp_crit_st = 0.d0

    if(allocated(gl%visco)) call initialize_viscosity(gl)

    if(allocated(gl%tcx)) call initialize_thermal_conductivity(gl)

    call initialize_transport_properties(gl)

    !module BWR
    gl%ncoeff=0.d0
    gl%gama_bwr=0.d0



    !module FNR
    !gl%TEMP_FNR_OLD = 0.d0
    !gl%DENS_FNR_OLD  = 0.d0
    !gl%tredmix_OLD = 0.d0
    !gl%rhoredmix_OLD = 0.d0
    !gl%nrsubst_old = 0
    !gl%FNR_OLD = 0.d0
    !gl%FNRDER_OLD = 0.d0
    !gl%delpi_old = 0.d0
    !gl%reg_term_old = 0.d0
    !gl%gauss_term_old = 0.d0
    !gl%nacalc_old = 0.d0
    !gl%deleps_old = 0.d0
    !gl%taugam_old = 0.d0
    !gl%zaehler = 0
    !gl%fnr_c = 0
    !gl%eq_type_old = 0
    !gl%mix_type_old = 0



    !module FNR_MIX
    !gl%TEMP_FNR_MIX_OLD = 0.d0
    !gl%DENS_FNR_MIX_OLD  = 0.d0
    !gl%tredmix_mix_OLD = 0.d0
    !gl%rhoredmix_mix_OLD = 0.d0
    !gl%FNR_MIX_OLD = 0.d0
    !gl%FNRDER_MIX_OLD = 0.d0
    !gl%zaehler_mix = 0
    !gl%cn_fnm = 0
    !gl%za_fnm = 0
    !gl%uncomp = 0

    !not necessary anymore    
    !Andreas Jäger, March 2017
    !module_NRTL
    !gl%alpha_ij_0_NRTL = 0.D0
    !gl%alpha_ij_1_NRTL = 0.D0
    !gl%tau_Aij_NRTL = 0.D0
    !gl%tau_Bij_NRTL = 0.D0
    !gl%tau_Cij_NRTL = 0.D0
    !gl%tau_Dij_NRTL = 0.D0
    !gl%tau_Eij_NRTL = 0.D0
    !gl%tau_Fij_NRTL = 0.D0
    !gl%alpha_ij_NRTL = 0.D0
    !gl%tau_ij_NRTL = 0.D0
    !gl%G_ij_NRTL = 0.D0

    !not needed anymore
    !Andreas Jäger, March 2017
    !module_UNIFAC
    !gl%R_ik = 0.D0    !Group volume parameters for group k of molecule i (gl,i=30, k=100)
    !gl%Q_ik = 0.D0    !Group area parameters for group k of molecule i (i=30, k=100)
    !gl%v_ik = 0.D0    !Number of groups k in molecules of type i (i=30, k=100)
    !gl%a_nm = 0.D0    !Interaction parameter a for group n and m
    !gl%b_nm = 0.D0    !Interaction parameter b for group n and m
    !gl%c_nm = 0.D0    !Interaction parameter c for group n and m

    gl%nr_of_groups_i = 0       !Number of different groups of which molecule i is composed of (example: propane 2x1 + 1x2 -> nr_of_groups(propane) = 2)
    gl%group_start_index = 0    !This vector contains the indices where the groups of each component start, when all groups are written in one array (in maingroups_mix_list)
    gl%nr_of_groups_mix = 0                     !Total number of groups in the mixture (with double groups)
    gl%subgroups_ik_list = 0     !List that contains the subgroups k (max 100) of every fluid i (max 30) in the mixture
    gl%maingroups_mix_list= 0     !List that contains all maingroups that appear in the mixture


    !module_HelmgE
    !PACKING FRACTION
    !!!!!!!!!!!!!!!!!!!!!!
    gl%u_pack = 1.17D0
    !!!!!!!!!!!!!!!!!!!!!!
    gl%LIN_or_LB = 2       !1: LB mixing rules, 2: linear mixing rules for Helm+gE

    !Erik ,April 2018
    !SH June 18: Check if gl%cosmo was already allocated.
    if (allocated(gl%cosmo)) call initialize_COSMO_SAC(gl)

    if (allocated(gl%litref)) call initialize_literature(gl)

    call initialize_flash(gl)

    call initialize_pcSAFT(gl)






    !module_costald_eq

    gl%omega = 0.d0
    gl%vstar = 0.d0

    gl%omegaCOS = 0.d0
    gl%vstarCOS = 0.d0
    gl%TcCOS = 0.d0

    gl%k1 = 0.D0
    gl%k2 = 0.D0
    gl%RKM_V = 0.D0
    gl%RKM_ps_me = 0.D0
    gl%RKM_x = 0.D0
    gl%RKM_M = 0.D0
    gl%RKM_Tc = 0.D0

    gl%RKM_fluids = ""

    !end type

    call initialize_hdrt_property_definition(gl)

    call initialize_ecs(gl)

    end subroutine initialize

    
    
    subroutine initialize_viscosity(gl)
    
    implicit none

    type(type_gl) :: gl

    !Module_viscosity

    gl%visco%eta_read =.false.
    gl%visco%h2o_read=.false.
    gl%visco%coll_read=.false.
    gl%visco%etamodel=' '
    gl%visco%pointercrit_eta=' '
    gl%visco%pointerceaux=' '
    gl%visco%pointer_lambda=' '
    gl%visco%pointer_collmod=' '
    gl%visco%low_temp_eta=0.d0
    gl%visco%upp_temp_eta=0.d0
    gl%visco%max_dens_eta=0.d0
    gl%visco%red_temp_eta=0.d0
    gl%visco%red_dens_eta=0.d0
    gl%visco%red_vis_eta=0.d0
    gl%visco%Chapman_Enskog1=0.d0
    gl%visco%Chapman_Enskog2=0.d0
    gl%visco%nterm1_eta=0
    gl%visco%nterm2_eta=0
    gl%visco%nterm3_eta=0
    gl%visco%nterm4_eta=0
    gl%visco%nterm5_eta=0
    gl%visco%nterm6_eta=0
    gl%visco%nterm7_eta=0
    gl%visco%nterm8_eta=0
    gl%visco%nterm9_eta=0
    gl%visco%nterm_eta=0
    gl%visco%coeff_dil=0.d0
    gl%visco%exp_dil=0.d0
    gl%visco%coeff_hieta=0.d0
    gl%visco%hieta=0.d0
    gl%visco%iexpeta=0
    gl%visco%jexpeta=0
    gl%visco%exp1eta=0.d0
    gl%visco%exp2eta=0.d0
    gl%visco%coeff_ai=0.d0
    gl%visco%exp_ai=0.d0
    gl%visco%coeff_bi=0.d0
    gl%visco%exp_bi=0.d0
    gl%visco%coeff_ci=0.d0
    gl%visco%exp_ci=0.d0
    gl%visco%coeff_di=0.d0
    gl%visco%exp_di=0.d0
    gl%visco%coeff_ei=0.d0
    gl%visco%exp_ei=0.d0
    gl%visco%bieta=0.d0
    gl%visco%tieta=0.d0
    gl%visco%cieta=0.d0
    gl%visco%ti_eta=0.d0
    gl%visco%alpha_eta=0.d0
    gl%visco%exp_tau_eta=0.d0
    gl%visco%exp_del_eta=0.d0
    gl%visco%ndilgas_term=0
    gl%visco%ndens_initial=0
    gl%visco%coeff_in=0.d0
    gl%visco%texp_in=0.d0
    gl%visco%ndilgas_term=0
    gl%visco%tinit_red=0.d0
    gl%visco%etainit_red=0.d0
    gl%visco%pointer_eta=' '
    gl%visco%lejo_sigma=0.d0
    gl%visco%lejo_epka=0.d0
    gl%visco%c1_eta=0.d0
    !gl%visco%tred_eta=0.d0
    !gl%visco%rhored_eta=0.d0
    !gl%visco%etared=0.d0

    gl%visco%Ninuer=0.d0
    gl%visco%tinuer=0.d0
    gl%visco%dinuer=0.d0
    gl%visco%di0nuer=0.d0
    gl%visco%linuer=0.d0

    gl%visco%lowbo_temp_mue2=0.d0
    gl%visco%uppbo_temp_mue2=0.d0
    gl%visco%lowbo_dens_mue2=0.d0
    gl%visco%uppbo_dens_mue2=0.d0
    gl%visco%rhoc_eta=0.d0
    gl%visco%corr_ref_mue2=0.d0
    gl%visco%qc_mue2=0.d0
    gl%visco%qd_mue2=0.d0
    gl%visco%nue_mue2=0.d0
    gl%visco%gamma_mue2=0.d0
    gl%visco%Xi0_mue2=0.d0
    gl%visco%gamma0_mue2=0.d0
    gl%visco%T_mue2=0.d0
    gl%visco%X_mue2=0.d0

    gl%visco%rhoc_eta=0.d0
    gl%visco%Tc_eta=0.d0

    gl%visco%low_tempcoll=0.d0
    gl%visco%upp_tempcoll=0.d0
    gl%visco%low_prescoll=0.d0
    gl%visco%upp_denscoll=0.d0
    gl%visco%nterm_coll=0
    gl%visco%coeffcoll=0.d0
    gl%visco%ticoll=0.d0

    gl%visco%d0exp=0.d0
    gl%visco%a_0=0.d0
    gl%visco%a_1=0.d0
    gl%visco%a_2=0.d0
    gl%visco%aa0=0.d0
    gl%visco%aa1=0.d0
    gl%visco%aa2=0.d0
    gl%visco%bb0=0.d0
    gl%visco%bb1=0.d0
    gl%visco%bb2=0.d0
    gl%visco%c0=0.d0
    gl%visco%c_1=0.d0
    gl%visco%c_2=0.d0
    gl%visco%cc0=0.d0
    gl%visco%cc1=0.d0
    gl%visco%cc2=0.d0

    !variables specifically for VS2
    gl%visco%const19_VS2 = 0.d0
    gl%visco%exp19_VS2 = 0.d0
    gl%visco%Fv_VS2 = 0.d0
    gl%visco%Ev_VS2 = 0.d0
    
    !variables specifically for VS4
    gl%visco%d0_vs4 = 0.d0
    gl%visco%d0exp = 0.d0
    gl%visco%a_0 = 0.d0
    gl%visco%a_1 = 0.d0
    gl%visco%a_2 = 0.d0
    gl%visco%aa0 = 0.d0
    gl%visco%aa1 = 0.d0
    gl%visco%aa2 = 0.d0
    gl%visco%b0_vs4 = 0.d0
    gl%visco%b1_vs4 = 0.d0
    gl%visco%b2_vs4 = 0.d0
    gl%visco%bb0 = 0.d0
    gl%visco%bb1 = 0.d0
    gl%visco%bb2 = 0.d0
    gl%visco%c0 = 0.d0
    gl%visco%c_1 = 0.d0
    gl%visco%c_2 = 0.d0
    gl%visco%cc0 = 0.d0
    gl%visco%cc1 = 0.d0
    gl%visco%cc2 = 0.d0
    gl%visco%dd_0 = 0.d0
    gl%visco%dd_1 = 0.d0
    gl%visco%dd_2 = 0.d0
    gl%visco%e_0 = 0.d0
    gl%visco%e_1 = 0.d0
    gl%visco%e_2 = 0.d0
    gl%visco%d0_vs4 = 0.d0
    gl%visco%denstermnr = 0
    gl%visco%coeffbetastar = 0.d0
    gl%visco%powerbetastar = 0.d0

    !variables specifically for VS5
    gl%visco%slj_VS5 = 0.d0
    gl%visco%tc_VS5 = 0.d0
    gl%visco%accen_VS5 = 0.d0
    gl%visco%dipolered_VS5 = 0.d0
    gl%visco%kappa_VS5 = 0.d0
    gl%visco%addchung_VS5 = 0.d0
    gl%visco%pointereta = ''
    gl%visco%slj = 0.d0
    gl%visco%eklj = 0.d0
    gl%visco%ce = 0.d0
    gl%visco%cep = 0.d0

    !variables specifically for VS9
    gl%visco%coeff = 0.d0
    gl%visco%texpo = 0.d0
    gl%visco%dexpo = 0.d0
    gl%visco%beta_VS9 = 0.d0
    gl%visco%eta_VS9 = 0.d0
    gl%visco%pi_VS9 = 0.d0
    gl%visco%li_VS9 = 0.d0
    gl%visco%par1 = 0.d0
    gl%visco%par2 = 0.d0
    
    gl%visco%dilgas_term = 0
    gl%visco%tmineta = 0.d0
    gl%visco%tmaxeta = 0.d0
    gl%visco%pmaxeta = 0.d0
    gl%visco%rhomaxeta = 0.d0
    gl%visco%omega_model = ''
    
    gl%visco%tredetadg = 0.d0
    gl%visco%tredeta = 0.d0
    gl%visco%rhoredeta = 0.d0
    gl%visco%visredetadg = 0.d0
    gl%visco%visredeta = 0.d0
    
    gl%visco%term_num1_eta = 0
    gl%visco%term_num2_eta = 0

    gl%visco%term_expo_eta = 0
    gl%visco%term_com_expo_eta = 0
    gl%visco%coeff_hi = 0.d0 

    end subroutine initialize_viscosity

    subroutine initialize_general_eos_parameters(gl)


    implicit none

    type(type_gl) :: gl



    !module_general_eos_parameters

    gl%tred = 0.d0
    gl%rhored = 0.d0
    gl%tredmix = 0.d0
    gl%rhoredmix = 0.d0
    gl%tref = 0.d0
    gl%pref = 0.d0
    gl%rhoref = 0.d0
    gl%href = 0.d0
    gl%sref = 0.d0

    gl%refstate = ''

    end subroutine initialize_general_eos_parameters

    subroutine initialize_vle(gl)


    implicit none

    type(type_gl) :: gl


    !module_VLE
    gl%rho_vap = 0.D0
    gl%rho_liq = 0.D0
    gl%rho_liq1 = 0.D0
    gl%rho_liq2 = 0.D0
    gl%T_pts = 0.D0
    gl%p_pts = 0.D0
    gl%rholiq_pts = 0.D0
    gl%rhovap_pts = 0.D0
    gl%x_pts = 0.D0
    gl%pointID = 0
    gl%phasenv_pts = 0
    !init SHOULD not be neccesary because logical elements in types are set to .false. during allocation
    !gl%savebounds_p = .false.   !MT 2018/03 new definition for types
    gl%min_factor = 0.5d0       !MT 2018/03 new definition for types

    end subroutine initialize_vle

    subroutine initialize_ideal_gas_coefficients(gl)



    implicit none

    type(type_gl) :: gl




    !module_ideal_gas_coefficients

    gl%c1 = 0.d0
    gl%c2 = 0.d0
    gl%tcp0red = 0.d0
    gl%cp0red = 0.d0

    gl%ncp0c = 0
    gl%ncp0poly = 0
    gl%ncp0pl = 0
    gl%ncp0log = 0
    gl%ncp0 = 0
    gl%ncp0logexp = 0
    gl%ncp0cosh = 0
    gl%ncp0sinh = 0
    gl%ncp0el = 0

    !gl%cp0coeff = 0.d0
    !gl%cp0exp = 0.d0
    !gl%cp0hyp1 = 0.d0
    !gl%cp0hyp2 = 0.d0
    !gl%cp0hyp3 = 0.d0

    gl%cpmodel = .false.
    gl%phkmodel = .false.
    gl%phmodel = .false.

    end subroutine initialize_ideal_gas_coefficients

    subroutine  initialize_thermal_conductivity(gl)


    implicit none

    type(type_gl) :: gl



    gl%tcx%tcmodel = ''
    gl%tcx%tkmodel = ''
    gl%tcx%pointer_hardtc = ''

    gl%tcx%num_dil_tc = 0
    gl%tcx%den_dil_tc = 0
    gl%tcx%num_bgrd_tc = 0
    gl%tcx%den_bgrd_tc = 0


    gl%tcx%tcx_read = .false.

    gl%tcx%tcx_dil = 0.d0
    gl%tcx%tcx_bgrd = 0.d0
    gl%tcx%tmin_tc = 0.d0
    gl%tcx%tmax_tc = 0.d0
    gl%tcx%pmax_tc = 0.d0
    gl%tcx%rhomax_tc = 0.d0
    gl%tcx%tred_dil  = 0.d0
    gl%tcx%tred_bgrd = 0.d0
    gl%tcx%rhored_bgrd = 0.d0

    gl%tcx%lj_sigma = 0.d0
    gl%tcx%lj_epskap = 0.d0
    gl%tcx%const_eq20 = 0.d0
    gl%tcx%Texp_eq20 = 0.d0
    gl%tcx%F_tc3 = 0.d0
    gl%tcx%rm = 0.d0
    gl%tcx%eta0 = 0.d0
    gl%tcx%Fvi = 0.d0
    gl%tcx%Evi = 0.d0
    gl%tcx%tcxlam0_num = 0
    gl%tcx%tcxlam0_den = 0
    gl%tcx%dellam = 0
    gl%tcx%dellamcr = 0

    gl%tcx%tred_tc = 0.d0
    gl%tcx%rhored_tc = 0.d0
    gl%tcx%etared_tc = 0.d0

    gl%tcx%lam0_coeff = 0.d0
    gl%tcx%dellam_coeff = 0.d0
    gl%tcx%dellamcr_coeff = 0.d0
    gl%tcx%lam0_exp = 0.d0
    gl%tcx%dellam_exp = 0.d0
    gl%tcx%dellam_exp2 = 0.d0
    gl%tcx%dellamcr_exp = 0.d0
    gl%tcx%num_tc_coeff_dil = 0.d0
    gl%tcx%den_tc_coeff_dil = 0.d0
    gl%tcx%den_tc_exp_dil = 0.d0
    gl%tcx%num_tc_powerT_dil = 0.d0
    gl%tcx%tmin_tk = 0.d0
    gl%tcx%tmax_tk = 0.d0
    gl%tcx%pmax_tk = 0.d0
    gl%tcx%rhomax_tk = 0.d0
    gl%tcx%term_tk = 0
    gl%tcx%tred_tk = 0.d0
    gl%tcx%rhored_tk = 0.d0
    gl%tcx%tcx_tk = 0.d0
    gl%tcx%gnu_tk = 0.d0
    gl%tcx%gamma_tk = 0.d0
    gl%tcx%R0_tk = 0.d0
    gl%tcx%z_visonly_tk = 0.d0
    gl%tcx%c_tk = 0.d0
    gl%tcx%xi0_tk = 0.d0
    gl%tcx%gam0_tk = 0.d0
    gl%tcx%qd_inverse_tk = 0.d0
    gl%tcx%tref_tk = 0.d0
    gl%tcx%tcx_crit_read = .false.
    gl%tcx%rholl_r23 = 0.d0
    gl%tcx%b1_r23 = 0.d0
    gl%tcx%b2_r23 = 0.d0
    gl%tcx%c1l_r23 = 0.d0
    gl%tcx%c2l_r23 = 0.d0
    gl%tcx%delgl_r23 = 0.d0
    gl%tcx%lmax_r23 = 0.d0

    gl%tcx%npnum_tk1 = 0
    gl%tcx%npdenom_tk1 = 0
    gl%tcx%nexp_tk1 = 0
    gl%tcx%nspare_tk1 = 0

    gl%tcx%a_tk1 = 0.d0
    gl%tcx%tsum_tk1 = 0.d0
    gl%tcx%texp_tk1 = 0.d0
    gl%tcx%dsum_tk1 = 0.d0
    gl%tcx%dexp_tk1 = 0.d0


    end subroutine  initialize_thermal_conductivity

    subroutine initialize_transport_properties(gl)


    implicit none

    type(type_gl) :: gl



    gl%tmintrans = 0.d0
    gl%tmaxtrans = 0.d0
    gl%pmaxtrans = 0.d0
    gl%rhomaxtrans = 0.d0

    end subroutine initialize_transport_properties

    subroutine initialize_flash(gl)


    implicit none

    type(type_gl) :: gl



    gl%startvaluespin = .false.

    end subroutine initialize_flash



    subroutine initialize_hdrt_property_definition(gl)

    implicit none

    type(type_gl) :: gl

    !gl%numderiv = .false.
    !gl%print_head = .false.
    !gl%show_progress = .false.
    !gl%CompositionType = 0
    !gl%fourphase_physical = .false.
    !
    !gl%Q_map = 0
    !gl%Q_point_mix = ""
    !gl%press_Q_mix = 0.d0
    !gl%Temp_Q_mix = 0.d0
    !gl%phasefrac_mix = 0.d0
    !gl%phaseemerge_mix = 0
    !
    !gl%x_Q_ph1_mix = 0.d0
    !gl%x_Q_ph2_mix = 0.d0
    !gl%x_Q_ph3_mix = 0.d0
    !gl%x_Q_ph4_mix = 0.d0
    !gl%Qexist_mix = .false.
    !
    !gl%np_3ph = 0
    !gl%n_3ph_lines = 0
    !gl%nflashs = 0
    !gl%loops = 0
    !gl%trline_mix = ""
    !gl%press_tr_mix = 0.d0
    !gl%temp_tr_mix = 0.d0


    gl%gw_B0 = 0.d0
    gl%hw_B0 = 0.d0
    gl%gw_B0_a0 = 0.d0
    gl%hw_B0_a0 = 0.d0
    gl%kb_gw = 0.d0
    gl%kb_hw = 0.d0

    gl%alfa = 0.d0
    gl%T0_a = 0.d0
    gl%p0_a = 0.d0

    gl%a0_hdrt = 0.d0
    gl%a_hc = 0.d0
    gl%sigma = 0.d0
    gl%eps_k = 0.d0
    gl%sigmad = 0.d0
    gl%epsd_k = 0.d0
    gl%r_vdw = 0.d0
    gl%a_ref = 0.d0
    gl%T_ref = 0.d0
    gl%p_ref = 0.d0
    gl%l_core = 0.d0
    gl%R_shell = 0.d0
    gl%z_shell = 0

    gl%N_hdrts = 0
    gl%N_guests = 0
    gl%moles_hdrt = 0.d0
    gl%Y_yj_tf = 0.d0
    gl%B0_hdrt = 0.d0
    gl%B0_m = (/10.d9, 4.d0/)
    gl%a0_m = 0.d0
    gl%a_T_ref = 0.d0

    gl%occupmax = 0

    gl%KS_s = 0.d0
    gl%KS_l = 0.d0

    gl%kapa = 0.d0
    gl% m_const = 0.d0
    gl%rs_const = 0.d0

    gl%mixdecision = ""

    gl%cpH_a = 0.12739d0     ! [J.mol^-1] constant in cp(T): cp = cp_a*T + cp_b
    gl%cpH_b = 2.93401d0     ! [J.mol^-1.K^-1] constant in cp(T) expression
    gl%T0 = 273.15d0        ! [K] reference temperature for Gibbs energy
    gl%p0 = 1.d0           ! [Pa] reference pressure for Gibbs energy = 1 Pa
    gl%Nw = 0.d0               ! [-] number of water molecules in unit cell of given hydrate structure
    gl%v_cavi = (/0.d0, 0.d0, 0.d0/) ! first small then large cavity (third is for medium cavity in sH structure)
    gl%N_cavi = 2                                         !Number of cavities in the given hydrate structure sI & sII = 2, sH = 3



    !!Constants g0 - g10, n
    !double precision :: g0, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10
    !double precision :: exp_n
    !!constants g0_alpha to g8_alpha
    !double precision :: g0_a, g1_a, g2_a, g3_a, g4_a, g5_a, g6_a, g7_a, g8_a
    !!constants g0_kappa to g2_kappa
    !double precision :: g0_k, g1_k, g2_k
    !!Reference point
    !double precision :: T0_dryice, p0_dryice, R_CO2
    !!constant pi
    !double precision :: const_pi

    end subroutine


    subroutine initialize_ecs(gl)

    implicit none

    type(type_gl) :: gl


    !type module_ecs


    gl%input_0 = 'td'
    gl%t_0 = 300.d0
    gl%d_0 = 10.d0
    gl%fluids_0 = ""
    gl%moles_0 = "1"
    gl%molev_0 = 0.d0
    gl%eos_indicator_0 = ""
    gl%eos_indicator_0_call = "1"
    gl%tc_0 = 0.d0
    gl%dc_0 = 0.d0
    gl%m_0 = 0.d0
    gl%sigma_0 = 0.d0
    gl%epsilon_0 = 0.d0

    !fluid variables
    gl%tc_fluid = 0.d0
    gl%dc_fluid = 0.d0
    gl%w_fluid = 0.d0
    gl%m_fluid = 0.d0
    gl%ncomp_fluid = 0
    gl%alpha_temp = 0.d0
    gl%dalpha_ddelta_temp = 0.d0
    gl%d2alpha_ddelta2_temp = 0.d0
    gl%dalpha_dtau_temp = 0.d0
    gl%d2alpha_ddelta_dtau_temp = 0.d0
    gl%z_temp = 0.d0
    gl%residual_ecs = 0.d0

    !end type module_ecs
    end subroutine


    !Erik, April 2018
    subroutine initialize_COSMO_SAC(gl)

    implicit none

    type(type_gl) :: gl
    !allocate(gl%cosmo)
    gl%cosmo%molefraction_cosmo_prev = 0.D0
    gl%cosmo%temp_cosmo_prev = 0.D0
    gl%cosmo%seggamma_pure_prev = 0.D0
    gl%cosmo%seggamma_prev = 0.D0
    gl%cosmo%ln_gamma_C_cosmo_prev = 0.D0
    gl%cosmo%ln_gamma_R_cosmo_prev = 0.D0
    gl%cosmo%gE_C_cosmo_prev = 0.D0
    gl%cosmo%gE_R_cosmo_prev = 0.D0
    end subroutine


    subroutine initialize_pcSAFT(gl)

    ! Henning Markgraf, June 2016

    ! Initializations and dimension allocation of the module_PCSAFT variables


    implicit none

    type(type_gl) :: gl





    gl%SAFTmodel_assoc = ''

    gl%molfractions_save = 0.d0
    gl%mol_save = .false.


    ! binary interaction parameter
    gl%kij_PCSAFT = 0.d0

    !parameters for PC SAFT m, sigma, epsilon/k
    gl%mPCSAFT = 0.d0
    gl%sigPCSAFT = 0.d0
    gl%epskPCSAFT = 0.d0
    !----------------------------------------------------
    !function part variables to avoid double calculation

    ! in the beginning, the arrays in the module variables are not allocated yet
    if (gl%already_allocated_PCSAFT) then

        deallocate(gl%mmeanQij_PCSAFT)
        deallocate(gl%mmeanQijk_PCSAFT)
        deallocate(gl%ab_PCSAFTQ)
        deallocate(gl%c_PCSAFTQ)
        deallocate(gl%mmeanDij_PCSAFT)
        deallocate(gl%mmeanDijk_PCSAFT)
        deallocate(gl%ab_PCSAFTD)
        deallocate(gl%c_PCSAFTD)

        deallocate(gl%ar_calculated_PCSAFT)

        ! the dxi (x1) derivatives:
        deallocate(gl%arx1_calculated_PCSAFT)

        ! the dxi dxj (x2) derivatives
        deallocate(gl%arx2_calculated_PCSAFT)

        ! the dxi dxj dxk (x3) derivatives
        deallocate(gl%arx3_calculated_PCSAFT)

        !--------------------------------------------------------------
        ! arrays that store the results of the calculated derivatives:
        ! T,rho-derivative variables
        deallocate(gl%ar_PCSAFT) ! note: is the option to catch errors allocate(ar_PCSAFT,stat=allocation_status) necessary?
        deallocate(gl%adisp_PCSAFT)
        deallocate(gl%c_PCSAFT)
        deallocate(gl%i1_PCSAFT)
        deallocate(gl%i2_PCSAFT)
        deallocate(gl%meo1_PCSAFT)
        deallocate(gl%meo2_PCSAFT)
        deallocate(gl%ahc_PCSAFT)
        deallocate(gl%ahs_PCSAFT)
        deallocate(gl%z0_PCSAFT)
        deallocate(gl%z1_PCSAFT)
        deallocate(gl%z2_PCSAFT)
        deallocate(gl%z3_PCSAFT)
        deallocate(gl%gii_PCSAFT)
        deallocate(gl%di_PCSAFT)
        deallocate(gl%gij_PCSAFT)

        ! dxi (x1) derivative variables
        deallocate(gl%arx1_PCSAFT)
        deallocate(gl%adispx1_PCSAFT)
        deallocate(gl%cx1_PCSAFT)
        deallocate(gl%i1x1_PCSAFT)
        deallocate(gl%i2x1_PCSAFT)
        deallocate(gl%meo1x1_PCSAFT)
        deallocate(gl%meo2x1_PCSAFT)
        deallocate(gl%ahcx1_PCSAFT)
        deallocate(gl%ahsx1_PCSAFT)
        deallocate(gl%z0x1_PCSAFT)
        deallocate(gl%z1x1_PCSAFT)
        deallocate(gl%z2x1_PCSAFT)
        deallocate(gl%z3x1_PCSAFT)
        deallocate(gl%giix1_PCSAFT)
        deallocate(gl%abx1_PCSAFT)

        ! dxi dxj (x2) derivative variables
        deallocate(gl%arx2_PCSAFT)
        deallocate(gl%adispx2_PCSAFT)
        deallocate(gl%cx2_PCSAFT)
        deallocate(gl%i1x2_PCSAFT)
        deallocate(gl%i2x2_PCSAFT)
        deallocate(gl%meo1x2_PCSAFT)
        deallocate(gl%meo2x2_PCSAFT)
        deallocate(gl%ahcx2_PCSAFT)
        deallocate(gl%ahsx2_PCSAFT)
        deallocate(gl%giix2_PCSAFT)
        deallocate(gl%abx2_PCSAFT)

        ! dxi dxj dxk (x3) derivative variables
        deallocate(gl%arx3_PCSAFT)
        deallocate(gl%adispx3_PCSAFT)
        deallocate(gl%cx3_PCSAFT)
        deallocate(gl%i1x3_PCSAFT)
        deallocate(gl%i2x3_PCSAFT)
        deallocate(gl%ahcx3_PCSAFT)
        deallocate(gl%ahsx3_PCSAFT)
        deallocate(gl%giix3_PCSAFT)
        deallocate(gl%abx3_PCSAFT)

        deallocate(gl%J2_PCSAFTQ)
        deallocate(gl%J3_PCSAFTQ)
        deallocate(gl%A2_PCSAFTQ)
        deallocate(gl%A3_PCSAFTQ)
        deallocate(gl%AQQ_PCSAFTQ)
        deallocate(gl%A2_PCSAFTD)
        deallocate(gl%A3_PCSAFTD)
        deallocate(gl%ADD_PCSAFTD)
        deallocate(gl%J2_PCSAFTD)
        deallocate(gl%J3_PCSAFTD)
        deallocate(gl%AASSOC_PCSAFT)
        deallocate(gl%xA_PCSAFT)
        deallocate(gl%delta_AB)

        ! dxi (x1) derivative variables
        deallocate(gl%A2X1_PCSAFTQ)
        deallocate(gl%A3X1_PCSAFTQ)
        deallocate(gl%AQQX1_PCSAFTQ)
        deallocate(gl%A2X1_PCSAFTD)
        deallocate(gl%A3X1_PCSAFTD)
        deallocate(gl%ADDX1_PCSAFTD)
        deallocate(gl%AASSOCX1_PCSAFT)
        deallocate(gl%J2X1_PCSAFTQ)
        deallocate(gl%J2X1_PCSAFTD)
        deallocate(gl%J3X1_PCSAFTQ)
        deallocate(gl%J3X1_PCSAFTD)

        ! dxi dxj (x2) derivative variables
        deallocate(gl%A2X2_PCSAFTQ)
        deallocate(gl%A3X2_PCSAFTQ)
        deallocate(gl%AQQX2_PCSAFTQ)
        deallocate(gl%A2X2_PCSAFTD)
        deallocate(gl%A3X2_PCSAFTD)
        deallocate(gl%ADDX2_PCSAFTD)
        deallocate(gl%AASSOCX2_PCSAFT)
        deallocate(gl%J2X2_PCSAFTQ)
        deallocate(gl%J2X2_PCSAFTD)
        deallocate(gl%J3X2_PCSAFTQ)
        deallocate(gl%J3X2_PCSAFTD)


    endif
    gl%already_allocated_PCSAFT = .FALSE.

    ! Initialize all variables
    !        call init_derivs_PCSAFT ! NOT ANY MORE, now everything in the main pc saft routines, on first run allocate_arrays_PCSAFT is called and initializes!

    !---------------------------------------

    end subroutine initialize_pcSAFT


    subroutine initialize_Gen_Eq_Igor(gl)


    implicit none

    type(type_gl) :: gl



    gl%coefficientsi(1,1) = 5.3410734d0
    gl%coefficientsi(2,1) = -2.2778189d0
    gl%coefficientsi(3,1) = -3.8785499d0
    gl%coefficientsi(4,1) = -0.012190959d0
    gl%coefficientsi(5,1) = 0.00092942159d0
    gl%coefficientsi(6,1) = -0.016631229d0
    gl%coefficientsi(7,1) = -1.6572887d0
    gl%coefficientsi(8,1) = 1.2642606d0
    gl%coefficientsi(9,1) = 0.096008662d0
    gl%coefficientsi(10,1) = 0.09295083d0
    gl%coefficientsi(11,1) = -0.38271299d0
    gl%coefficientsi(12,1) = 0.34936066d0
    gl%coefficientsi(13,1) = 0.041718709d0
    gl%coefficientsi(14,1) = -0.012149915d0

    gl%coefficientsi(1,2) = 6.6819473d0
    gl%coefficientsi(2,2) = -1.2846893d0
    gl%coefficientsi(3,2) = -8.6095696d0
    gl%coefficientsi(4,2) = 0.36869492d0
    gl%coefficientsi(5,2) = 0.080731074d0
    gl%coefficientsi(6,2) = -0.080314182d0
    gl%coefficientsi(7,2) = 21.646346d0
    gl%coefficientsi(8,2) = 2.1645843d0
    gl%coefficientsi(9,2) = 4.4221976d0
    gl%coefficientsi(10,2) = -0.057463893d0
    gl%coefficientsi(11,2) = -2.0429713d0
    gl%coefficientsi(12,2) = 6.4055642d0
    gl%coefficientsi(13,2) = -0.90287649d0
    gl%coefficientsi(14,2) = -0.15474203d0

    gl%coefficientsi(1,3) = 1.6692414d0
    gl%coefficientsi(2,3) = 1.3795302d0
    gl%coefficientsi(3,3) = -2.6707821d0
    gl%coefficientsi(4,3) = -0.20627285d0
    gl%coefficientsi(5,3) = -0.081358186d0
    gl%coefficientsi(6,3) = -0.35343719d0
    gl%coefficientsi(7,3) = -16.018967d0
    gl%coefficientsi(8,3) = -2.5726222d0
    gl%coefficientsi(9,3) = 1.1591367d0
    gl%coefficientsi(10,3) = 0.44419682d0
    gl%coefficientsi(11,3) = 1.1751452d0
    gl%coefficientsi(12,3) = -8.3598749d0
    gl%coefficientsi(13,3) = 0.23069811d0
    gl%coefficientsi(14,3) = 0.23233099d0

    gl%coefficientsi(1,4) = 2.9446922d0
    gl%coefficientsi(2,4) = 2.3284396d0
    gl%coefficientsi(3,4) = 2.7960114d0
    gl%coefficientsi(4,4) = 0.6373147d0
    gl%coefficientsi(5,4) = 0.99619992d0
    gl%coefficientsi(6,4) = 1.1870929d0
    gl%coefficientsi(7,4) = 1.0375103d0
    gl%coefficientsi(8,4) = 1.3733437d0
    gl%coefficientsi(9,4) = 1.1168557d0
    gl%coefficientsi(10,4) = 0.7639042d0
    gl%coefficientsi(11,4) = 1.4829049d0
    gl%coefficientsi(12,4) = 1.0080516d0
    gl%coefficientsi(13,4) = 1.3320474d0
    gl%coefficientsi(14,4) = 1.2062411d0

    gl%exponentsi(1,1) = 1d0
    gl%exponentsi(2,1) = 1d0
    gl%exponentsi(3,1) = 1d0
    gl%exponentsi(4,1) = 3d0
    gl%exponentsi(5,1) = 7d0
    gl%exponentsi(6,1) = 2d0
    gl%exponentsi(7,1) = 1d0
    gl%exponentsi(8,1) = 1d0
    gl%exponentsi(9,1) = 2d0
    gl%exponentsi(10,1) = 5d0
    gl%exponentsi(11,1) = 1d0
    gl%exponentsi(12,1) = 1d0
    gl%exponentsi(13,1) = 4d0
    gl%exponentsi(14,1) = 2d0

    gl%exponentsi(1,2) = 0.686d0
    gl%exponentsi(2,2) = 1.118d0
    gl%exponentsi(3,2) = 0.857d0
    gl%exponentsi(4,2) = 0.559d0
    gl%exponentsi(5,2) = 0.442d0
    gl%exponentsi(6,2) = 0.831d0
    gl%exponentsi(7,2) = 0.484d0
    gl%exponentsi(8,2) = 2.527d0
    gl%exponentsi(9,2) = 1.549d0
    gl%exponentsi(10,2) = 0.757d0
    gl%exponentsi(11,2) = 3.355d0
    gl%exponentsi(12,2) = 1.905d0
    gl%exponentsi(13,2) = 4.941d0
    gl%exponentsi(14,2) = 12.805d0

    gl%exponentsi(1,3) = 0.d0
    gl%exponentsi(2,3) = 0.d0
    gl%exponentsi(3,3) = 0.d0
    gl%exponentsi(4,3) = 0.d0
    gl%exponentsi(5,3) = 0.d0
    gl%exponentsi(6,3) = 0.d0
    gl%exponentsi(7,3) = 1.d0
    gl%exponentsi(8,3) = 1.d0
    gl%exponentsi(9,3) = 1.d0
    gl%exponentsi(10,3) = 1.d0
    gl%exponentsi(11,3) = 2.d0
    gl%exponentsi(12,3) = 2.d0
    gl%exponentsi(13,3) = 2.d0
    gl%exponentsi(14,3) = 3.d0

    gl%anz_term_igor = 14

    end subroutine initialize_Gen_Eq_Igor

    subroutine initialize_Gen_Eq_Span(gl)


    implicit none

    type(type_gl) :: gl



    gl%coefficientss(1,1) = 0.636479524d0
    gl%coefficientss(2,1) = -1.74667493d0
    gl%coefficientss(3,1) = -0.0144442644d0
    gl%coefficientss(4,1) = 0.06799731d0
    gl%coefficientss(5,1) = 0.0000767320032d0
    gl%coefficientss(6,1) = 0.218194143d0
    gl%coefficientss(7,1) = 0.0810318494d0
    gl%coefficientss(8,1) = -0.0907368899d0
    gl%coefficientss(9,1) = 0.025312225d0
    gl%coefficientss(10,1) = -0.0209937023d0

    gl%coefficientss(1,2) = 0.82247342d0
    gl%coefficientss(2,2) = -0.954932692d0
    gl%coefficientss(3,2) = -0.745462328d0
    gl%coefficientss(4,2) = 0.182685593d0
    gl%coefficientss(5,2) = 0.0000547120142d0
    gl%coefficientss(6,2) = 0.761697913d0
    gl%coefficientss(7,2) = 0.415691324d0
    gl%coefficientss(8,2) = -0.825206373d0
    gl%coefficientss(9,2) = -0.240558288d0
    gl%coefficientss(10,2) = -0.0643818403d0

    gl%coefficientss(1,3) = -1.86193063d0
    gl%coefficientss(2,3) = 10.5083555d0
    gl%coefficientss(3,3) = 1.6403233d0
    gl%coefficientss(4,3) = -0.613747797d0
    gl%coefficientss(5,3) = -0.00069318829d0
    gl%coefficientss(6,3) = -7.05727791d0
    gl%coefficientss(7,3) = -2.90006245d0
    gl%coefficientss(8,3) = -0.232497527d0
    gl%coefficientss(9,3) = -0.282346515d0
    gl%coefficientss(10,3) = 2.54250643d0

    gl%exponentss(1,1) = 0.125d0
    gl%exponentss(2,1) = 1.125d0
    gl%exponentss(3,1) = 1.25d0
    gl%exponentss(4,1) = 0.25d0
    gl%exponentss(5,1) = 0.75d0
    gl%exponentss(6,1) = 0.625d0
    gl%exponentss(7,1) = 2.d0
    gl%exponentss(8,1) = 4.125d0
    gl%exponentss(9,1) = 4.125d0
    gl%exponentss(10,1) = 17.d0

    gl%exponentss(1,2) = 1.d0
    gl%exponentss(2,2) = 1.d0
    gl%exponentss(3,2) = 2.d0
    gl%exponentss(4,2) = 3.d0
    gl%exponentss(5,2) = 8.d0
    gl%exponentss(6,2) = 2.d0
    gl%exponentss(7,2) = 3.d0
    gl%exponentss(8,2) = 1.d0
    gl%exponentss(9,2) = 4.d0
    gl%exponentss(10,2) = 3.d0

    gl%exponentss(1,3) = 0.d0
    gl%exponentss(2,3) = 0.d0
    gl%exponentss(3,3) = 0.d0
    gl%exponentss(4,3) = 0.d0
    gl%exponentss(5,3) = 0.d0
    gl%exponentss(6,3) = 1.d0
    gl%exponentss(7,3) = 1.d0
    gl%exponentss(8,3) = 1.d0
    gl%exponentss(9,3) = 1.d0
    gl%exponentss(10,3) = 1.d0

    gl%exponentss(1,4) = 0.d0
    gl%exponentss(2,4) = 0.d0
    gl%exponentss(3,4) = 0.d0
    gl%exponentss(4,4) = 0.d0
    gl%exponentss(5,4) = 0.d0
    gl%exponentss(6,4) = 1.d0
    gl%exponentss(7,4) = 1.d0
    gl%exponentss(8,4) = 2.d0
    gl%exponentss(9,4) = 2.d0
    gl%exponentss(10,4) = 3.d0

    gl%anz_term_span = 10

    end subroutine initialize_Gen_Eq_Span



    subroutine initialize_Gen_Eq_Sun(gl)


    implicit none

    type(type_gl) :: gl



    gl%coefficientssun(1,1) = 0.970439249d0
    gl%coefficientssun(2,1) = 0.973671323d0
    gl%coefficientssun(3,1) = -2.96661981d0
    gl%coefficientssun(4,1) = 0.0784340496d0
    gl%coefficientssun(5,1) = 0.000278440866d0
    gl%coefficientssun(6,1) = -0.0677622221d0
    gl%coefficientssun(7,1) = -0.0856371936d0
    gl%coefficientssun(8,1) = 0.177467443d0
    gl%coefficientssun(9,1) = 0.391636018d0
    gl%coefficientssun(10,1) = -0.00803312946d0
    gl%coefficientssun(11,1) = -0.260385851d0
    gl%coefficientssun(12,1) = -0.0191104746d0
    gl%coefficientssun(13,1) = -0.063133147d0
    gl%coefficientssun(14,1) = -0.0227769095d0

    gl%coefficientssun(1,2) = 1.57750154d0
    gl%coefficientssun(2,2) = 1.15745614d0
    gl%coefficientssun(3,2) = -3.54867092d0
    gl%coefficientssun(4,2) = 0.118030671d0
    gl%coefficientssun(5,2) = 0.000302753897d0
    gl%coefficientssun(6,2) = -0.263074957d0
    gl%coefficientssun(7,2) = 0.0255299486d0
    gl%coefficientssun(8,2) = -0.126632996d0
    gl%coefficientssun(9,2) = 0.448343319d0
    gl%coefficientssun(10,2) = -0.00946702997d0
    gl%coefficientssun(11,2) = -0.443927529d0
    gl%coefficientssun(12,2) = -0.0168224827d0
    gl%coefficientssun(13,2) = -0.11586464d0
    gl%coefficientssun(14,2) = -0.0132417591d0

    gl%coefficientssun(1,3) = 0.34682192d0
    gl%coefficientssun(2,3) = 0.503423025d0
    gl%coefficientssun(3,3) = -0.35105957d0
    gl%coefficientssun(4,3) = 0.0507004866d0
    gl%coefficientssun(5,3) = 0.000199939129d0
    gl%coefficientssun(6,3) = -0.569888763d0
    gl%coefficientssun(7,3) = -0.196198912d0
    gl%coefficientssun(8,3) = -2.02509554d0
    gl%coefficientssun(9,3) = -1.09353609d0
    gl%coefficientssun(10,3) = 0.0725785202d0
    gl%coefficientssun(11,3) = 0.216072642d0
    gl%coefficientssun(12,3) = -0.10154263d0
    gl%coefficientssun(13,3) = 0.0746926106d0
    gl%coefficientssun(14,3) = 0.00218830463d0

    gl%exponentssun(1,1) = 1.d0
    gl%exponentssun(2,1) = 1.d0
    gl%exponentssun(3,1) = 1.d0
    gl%exponentssun(4,1) = 3.d0
    gl%exponentssun(5,1) = 7.d0
    gl%exponentssun(6,1) = 2.d0
    gl%exponentssun(7,1) = 1.d0
    gl%exponentssun(8,1) = 1.d0
    gl%exponentssun(9,1) = 2.d0
    gl%exponentssun(10,1) = 5.d0
    gl%exponentssun(11,1) = 1.d0
    gl%exponentssun(12,1) = 1.d0
    gl%exponentssun(13,1) = 4.d0
    gl%exponentssun(14,1) = 2.d0

    gl%exponentssun(1,2) = 1.5d0
    gl%exponentssun(2,2) = 0.25d0
    gl%exponentssun(3,2) = 1.25d0
    gl%exponentssun(4,2) = 0.25d0
    gl%exponentssun(5,2) = 0.875d0
    gl%exponentssun(6,2) = 1.375d0
    gl%exponentssun(7,2) = 0.d0
    gl%exponentssun(8,2) = 2.375d0
    gl%exponentssun(9,2) = 2.d0
    gl%exponentssun(10,2) = 2.125d0
    gl%exponentssun(11,2) = 3.5d0
    gl%exponentssun(12,2) = 6.5d0
    gl%exponentssun(13,2) = 4.75d0
    gl%exponentssun(14,2) = 12.5d0

    gl%exponentssun(1,3) = 0.d0
    gl%exponentssun(2,3) = 0.d0
    gl%exponentssun(3,3) = 0.d0
    gl%exponentssun(4,3) = 0.d0
    gl%exponentssun(5,3) = 0.d0
    gl%exponentssun(6,3) = 0.d0
    gl%exponentssun(7,3) = 1.d0
    gl%exponentssun(8,3) = 1.d0
    gl%exponentssun(9,3) = 1.d0
    gl%exponentssun(10,3) = 1.d0
    gl%exponentssun(11,3) = 2.d0
    gl%exponentssun(12,3) = 2.d0
    gl%exponentssun(13,3) = 2.d0
    gl%exponentssun(14,3) = 3.d0

    !Propane
    gl%acenfactor_1 = 0.153d0
    gl%polarfactor_1 = 0.d0

    !Octane
    gl%acenfactor_2 = 0.391d0
    gl%polarfactor_2 = 0.d0

    !Water
    gl%acenfactor_3 = 0.348d0
    gl%polarfactor_3 = 1.d0

    gl%anz_term_sun = 14

    end subroutine initialize_Gen_Eq_Sun

    subroutine initialize_literature(gl)

    implicit none
    type(type_gl) :: gl
    !module_literature_ref
    gl%litref%lit_ref_res = ''
    gl%litref%lit_ref_mix = ''

    gl%litref%lit_ref_de = ''
    gl%litref%lit_ref_tcx = ''
    gl%litref%lit_ref_stn = ''
    gl%litref%lit_ref_eta = ''
    gl%litref%lit_ref_eta_aux = ''
    gl%litref%lit_ref_mlt = ''
    gl%litref%lit_ref_sbl = ''


    end subroutine initialize_literature


    subroutine initialize_eos_coefficients(gl)

    implicit none
    type(type_gl) :: gl
    !module_eos_coefficients

    gl%eos_coeff%nreg = 0
    gl%eos_coeff%ncrt = 0
    gl%eos_coeff%nna = 0
    gl%eos_coeff%nlreg = 0
    gl%eos_coeff%nlcrt = 0
    gl%eos_coeff%nlna = 0

    gl%eos_coeff%I_pol = 0
    gl%eos_coeff%I_exp = 0
    gl%eos_coeff%I_GBS = 0
    gl%eos_coeff%I_NA = 0

    gl%eos_coeff%cppcheck = ''

    gl%eos_coeff%ni = 0.d0
    gl%eos_coeff%ti = 0.d0
    gl%eos_coeff%tli = 0.d0
    gl%eos_coeff%pli = 0.d0
    gl%eos_coeff%di = 0.d0
    gl%eos_coeff%p_i = 0
    gl%eos_coeff%gama = 0.d0
    gl%eos_coeff%eps = 0.d0
    gl%eos_coeff%beta = 0.d0
    gl%eos_coeff%gam = 0.d0
    gl%eos_coeff%eta = 0.d0
    gl%eos_coeff%etana = 0.d0
    gl%eos_coeff%eidna = 0.d0
    gl%eos_coeff%eitna = 0.d0
    if (allocated(gl%eos_coeff%h1)) then
        gl%eos_coeff%h1 = 0.d0
        gl%eos_coeff%h2 = 0.d0
        gl%eos_coeff%h3 = 0.d0
        gl%eos_coeff%h4 = 0.d0
        gl%eos_coeff%dof2 = 0.d0
        gl%eos_coeff%hard_sphere=.false.
    end if

    end subroutine initialize_eos_coefficients

    subroutine initialize_de(gl)

    implicit none
    type(type_gl) :: gl
    !module_dielectric_constant
    gl%de%de_read = .true.
    gl%de%demodel = ' '
    gl%de%tmin_de = 0.d0
    gl%de%tmax_de = 0.d0
    gl%de%pmax_de = 0.d0
    gl%de%rhomax_de = 0.d0
    gl%de%ref_temp_de = 0.d0
    gl%de%ref_dens_de = 0.d0
    gl%de%ref_press_de = 0.d0
    gl%de%term_num1_de = 0
    gl%de%term_num2_de = 0
    gl%de%term_num3_de = 0
    gl%de%term_num4_de = 0
    gl%de%term_num5_de = 0
    gl%de%term_num6_de = 0
    gl%de%coeffde = 0.d0
    gl%de%texpde = 0.d0
    gl%de%dexpde = 0.d0
    gl%de%pexpde = 0.d0

    end subroutine

    subroutine initialize_uncty(gl)

    implicit none
    type(type_gl) :: gl
    !module_uncertainty_parameters
    gl%uncty%blnCheck2Phase = .true.
    !UNCTY_FILE = ''
    gl%uncty%dblUncty_CP = 0.D0
    gl%uncty%dblUncty_D = 0.D0
    gl%uncty%dblUncty_WS = 0.D0
    gl%uncty%dblUncty_P = 0.D0
    gl%uncty%dblUncty_CV = 0.D0
    gl%uncty%dblUncty_CP_vle = 0.D0
    gl%uncty%dblUncty_D_vle = 0.D0
    gl%uncty%dblUncty_WS_vle = 0.D0
    gl%uncty%dblUncty_P_vle = 0.D0
    gl%uncty%dblUncty_CV_vle = 0.D0
    gl%uncty%estUncty_out = 0.D0
    gl%uncty%ni_varied_pos = 0.D0
    gl%uncty%ni_varied_neg = 0.D0
    gl%uncty%ni_varied_pos_ideal = 0.D0
    gl%uncty%ni_varied_neg_ideal = 0.D0
    !Andreas April 2013
    gl%uncty%Uncty_cp = 0.D0
    gl%uncty%Uncty_d = 0.D0
    gl%uncty%Uncty_w = 0.D0
    gl%uncty%Uncty_p = 0.D0
    gl%uncty%Uncty_cv = 0.D0
    gl%uncty%Uncty_h = 0.D0
    gl%uncty%Uncty_s = 0.D0
    gl%uncty%Uncty_u = 0.D0
    gl%uncty%temp_last = 0.D0
    gl%uncty%dens_last = 0.D0

    gl%uncty%rho_limit_uncty = 0.D0   !MT 2018/03 new definition for types
    gl%uncty%cp_limit_uncty = 0.D0    !MT 2018/03 new definition for types
    gl%uncty%cv_limit_uncty = 0.D0    !MT 2018/03 new definition for types
    gl%uncty%ws_limit_uncty = 0.D0    !MT 2018/03 new definition for types


    end subroutine initialize_uncty

    subroutine initialize_solids(gl)

    implicit none
    type(type_gl) :: gl

    !module_solids
    !Andreas May 2012
    !Hydrate_list contains all components for which hydrate equations are available
    gl%Hydrate_list = ''
    gl%Hydrate_list(1) = 'water'
    !old order
    !Hydrate_list(2) = 'co2'
    !Hydrate_list(3) = 'methane'
    !Hydrate_list(4) = 'ethane'
    !Hydrate_list(5) = 'co'
    !Hydrate_list(6) = 'nitrogen'
    !Hydrate_list(7) = 'propane'
    !Hydrate_list(8) = 'argon'
    !Hydrate_list(9) = 'oxygen'
    !changed old order to alphabetical order
    gl%Hydrate_list(2) = 'argon'
    gl%Hydrate_list(3) = 'co'
    gl%Hydrate_list(4) = 'co2'
    gl%Hydrate_list(5) = 'ethane'
    gl%Hydrate_list(6) = 'methane'
    gl%Hydrate_list(7) = 'nitrogen'
    gl%Hydrate_list(8) = 'oxygen'
    gl%Hydrate_list(9) = 'propane'

    !solid_list contains all components for which solid equations are available
    gl%Solid_list = ''
    gl%Solid_list(1) = 'water'
    gl%Solid_list(2) = 'co2'

    gl%Hydrate_formers  = ''
    !   Hydrate_pos = 0
    gl%nrofhydrateformers = 0
    gl%Fluidlist_hydrate = ''
    gl%moleslist_hydrate = 0.D0
    gl%mapping = 0

    gl%solid_subst = ''
    gl%solid_indicator = 0
    gl%nrofsolids = 0

    gl%p_Q_Hyd_2C = 0.D0
    gl%T_Q_Hyd_2C = 0.D0

    gl%pmax_Dryice = 500.D0
    gl%Tmin_Dryice = 80.D0
    gl%pmax_Waterice = 211.D0
    gl%Tmin_Waterice = 0.0D0

    gl%solidtype = 0

    !If the code should check for formation of solid phases (and three phase equilibria), this variable must be assigned the value .true., elsewise it must be .false.
    !Checking for solids might slow down the code
    !In the standard version, the code does not check for solids or three phase equilibria. Activate this check by putting a "+" at the end of the input combination (e.g. "TP+", "PH+")
    !Commented, because setinput always sets check_solid = .false.. Andreas, Feb 2016
    !check_solid = .false.

    !Save information of the equilibrium found
    gl%solidtype_akt_phase = 0
    gl%solidpos_akt_phase = 0

    end subroutine


    end module initialize_module
