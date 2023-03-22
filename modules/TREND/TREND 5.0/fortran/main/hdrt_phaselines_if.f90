    !DEC$ IF DEFINED(HDRT)
    module hdrt_phaselines_calc

    use module_all_types
    !use module_regula_falsi
    !use rhomix_pt_module
    !use flash_module
    !use phasedet_mix_module
    !use hdrt_chem_pot_module
    !use ptflash_solids_mix_module
    !use hdrt_properties_modules_module
    !use reduced_parameters_calc_module
    !use hdrt_main_module

    use module_hdrt_phaselines
    !use omp_lib

    interface
    !******************************************************************************
    module subroutine quadruplepoints (gl, path, Q_calc, Q_point, Qexist, Temp_Q, press_Q, &
        & x, x_ph1, x_ph2, x_ph3, x_ph4)
    !******************************************************************************
    ! Calculation of all existing quadruple points for pure hydrate



    implicit none
    type(type_gl) :: gl
    character(len=255) :: path

    double precision, dimension(30) :: x, lnfi !<-lnfi is currently used in the commented part ONLY
    double precision:: Temp, Dens
    double precision:: beta, rhoph1_est, rhoph2_est, rhoph3_est, rhophX_est  !, fug_gas
    integer :: error, i, solidnr, solid !<- solidnr and solid are used in the commented part ONLY
    integer :: iPhase, iFlash, iter, iPhase_sol

    !Pressures and Temperatures at each Quadrupelpoint
    character(len=6), dimension(4)  :: Q_point       !Quadruple point
    double precision, dimension (6) :: press_Q, Temp_Q
    double precision, dimension(4,30) :: x_ph1, x_ph2, x_ph3, x_ph4, x_phX
    logical, dimension (4):: Qexist
    integer, dimension(5) :: error_Q, Q_calc

    ! Mixed hydrates 2015
    double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double !<- CiJ and occup are used in the commented part ONLY
    double precision, dimension(30):: xH, fug_gas!, occup_ls, occup_ld, occup_sms, occup_smd !<- xH and fug_gas are used in the commented part ONLY

    !New quadruple point routine
    double precision, dimension(4):: rho_qpts, Phasefrac
    double precision, dimension(30):: x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd
    double precision:: rhofluid1_est, rhofluid2_est, rhofluid3_est, rhofluid4_est
    integer:: Phasefrac_0

    !//////////////////////////////////////////////////////////////////

    end subroutine quadruplepoints
    !******************************************************************************
    !******************************************************************************



    !******************************************************************************
    module subroutine threephase_equil(gl, pel, Eqtype, cur_3ph,iFlash,Temp, press, x , x_ph1, x_ph2, x_ph3, &
        &   rho_ph1, rho_ph2, rho_ph3, Chempot, lnfi, occup, error)!, &
    !&   Q_point, Qexist,Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)
    !******************************************************************************
    ! Calculation of a single three-phase equilibrium point by calling 'ptflash_solid_NC_3P'
    ! iTp = 1 ... temperature given, iTp = 2 ... pressure given




    implicit none
    type(type_gl) :: gl
    type(type_hdrt_pheq) :: pel
    character(6), intent(in) :: Eqtype        !Type of Equilibrium, which should be calculated
    integer, intent(in) :: iflash

    !Variables at the three-phase equilibrium
    double precision :: Temp, press
    double precision, dimension(30) :: x, x_ph1, x_ph2, x_ph3, x_phX
    double precision :: rho_ph1, rho_ph2, rho_ph3
    double precision, dimension(30) :: Chempot, lnfi!, occup_ls, occup_ld, occup_sms, occup_smd
    double precision:: Dens
    double precision:: beta, rhoph1_est, rhoph2_est, rhoph3_est, rhophX_est
    double precision, dimension(5):: rho
    double precision, dimension(3):: rho_sol
    integer :: error
    integer :: iPhase, iter, iPhase_sol!, iFlash
    integer :: cur_3ph

    !Variables at the Quadruple point
    !character(len=6), dimension(4)  :: Q_point
    !double precision, dimension (4) :: press_Q, Temp_Q
    !double precision, dimension(4,30) :: x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4
    !logical, dimension (4):: Qexist

    double precision, dimension(3) :: Phasefrac

    ! Mixed hydrates 2015
    double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double
    double precision, dimension(30):: xH, fug_gas

    !//////////////////////////////////////////////////////////////////

    end subroutine threephase_equil
    !******************************************************************************
    !******************************************************************************

    module subroutine phaselines_init(phl3)

    !use module_hdrt_phaselines

    implicit none

    type(type_hdrt_3ph_lines) :: phl3


    end subroutine phaselines_init

    module subroutine phaselines_out(gl, pel, phaseline, pathIn, errorflag)

    !use module_all_types
    !use module_hdrt_phaselines

    implicit none

    type(type_gl) :: gl
    type(type_hdrt_pheq) :: pel

    integer :: threephlineunit, fourphlineunit, eqtype_ind, phaseline, cur_3ph
    integer :: i, k, errorflag, row
    character(255) :: path, pathIn


    end subroutine phaselines_out


    module subroutine phaselines_expert(gl, pel, pathIn, additional_args, errorflag)
    !use module_all_types
    !use module_hdrt_phaselines

    implicit none

    type(type_gl), dimension(:) :: gl
    type(type_hdrt_pheq) :: pel

    double precision, dimension(30) :: xIn
    double precision, dimension(3,30) :: x_ph
    double precision, dimension(3) :: rho
    double precision :: p_start, p_end, press, press_step
    double precision :: T_start, T_end, Temp, Temp_step

    ! iTp = 1 ... temperature given, iTp = 2 ... pressure given
    integer :: iFlash, iTp
    integer :: errorflag

    character(255) :: pathIn
    character(100) :: additional_args
    character(30), dimension(30):: args_split


    integer :: cur_3ph, k, j, i

    end subroutine phaselines_expert



    !******************************************************************************
    module subroutine phaselines (gl, pel, pathIn, x_vap_given, errorflag)
    !******************************************************************************
    ! Subroutine calculates three- and four-phase equilibrium lines and calls 'threephase_equil'




    implicit none
    type(type_gl), dimension(:) :: gl
    type(type_gl), dimension(:), allocatable :: gl_parallel
    type(type_hdrt_pheq) :: pel
    type(type_hdrt_pheq), dimension(:), allocatable :: pel_parallel

    !> current three phase line local variable
    integer :: cur_3ph
    !> unit for outputfiles
    integer :: threephlineunit, fourphlineunit

    character(255) :: pathIn, path, filename, dummy
    character(6) :: EqtypeIn, Eqtype        !Type of Equilibrium, which should be calculated
    logical :: Eqtype_found, VLwH_done


    !Variables at/for the three-phase equilibrium
    double precision :: Temp, press, presstemp
    double precision, dimension(30) :: xIn, x, x_ph1, x_ph2, x_ph3, x_phX, x_vap_given
    double precision :: rho_ph1, rho_ph2, rho_ph3
    double precision, dimension(30) :: Chempot, lnfi
    double precision:: Dens, press_step, max_press, min_press
    double precision:: beta, rhoph1_est, rhoph2_est, rhoph3_est, rhophX_est

    integer :: trloop, i, k, iTp, iflash

    integer :: row
    integer :: iPhase, iPhase_sol
    logical :: write_3ph
    !>Indicate which Qpoint is used for start values
    !!EqType is for quadruple line, iQ is for startpoint with desired phasefractions (plural!) = 0
    integer :: EqTypeStart, EqTypeEnd, iQ

    ! Mixed hydrates 2015
    double precision, dimension(30) :: xH
    double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double
    double precision, dimension(3) :: Phasefrac

    ! Properties for enthalpy of formation and hydration number
    integer :: ii, J, refin_VHIw
    double precision :: h_vap, h_hyd, s_hyd, h_liq, h_melt, h_phase, h_sol, h_WaterIce
    double precision, dimension(30) :: fug_g, hyd_nr!, occup_ls, occup_ld, occup_sms, occup_smd

    integer :: errorflag, eqtype_ind, avail_cores

    !New variables for tpd calculation
    double precision:: tpd_sol, ChemPot_hyd, a0_other_struct
    integer:: hdrt_structure_flag

    !//////////////////////////////////////////////////////////////////

    end subroutine phaselines

    module subroutine phaselines_3 (gl, pel, xIn, x_vap_given, cur_3ph, errorflag)
    !******************************************************************************
    ! Subroutine calculates three-phase equilibrium lines and calls 'threephase_equil'




    implicit none
    type(type_gl) :: gl
    type(type_gl), dimension(:), allocatable :: gl_parallel
    type(type_hdrt_pheq) :: pel
    type(type_hdrt_pheq), dimension(:), allocatable :: pel_parallel

    !> current three phase line local variable
    integer :: cur_3ph

    character(255) :: path, filename, dummy
    character(6) :: EqtypeIn, Eqtype        !Type of Equilibrium, which should be calculated
    logical :: Eqtype_found, VLwH_done


    !Variables at/for the three-phase equilibrium
    double precision :: Temp, press, presstemp
    double precision, dimension(30) :: xIn, x, x_ph1, x_ph2, x_ph3, x_phX, x_vap_given
    double precision :: rho_ph1, rho_ph2, rho_ph3
    double precision, dimension(30) :: Chempot, lnfi
    double precision:: Dens, press_step, max_press, min_press
    double precision:: beta, rhoph1_est, rhoph2_est, rhoph3_est, rhophX_est

    integer :: trloop, i, k, iTp, iflash

    integer, dimension(1) :: row
    integer :: iPhase, iPhase_sol
    logical :: write_3ph
    !>Indicate which Qpoint is used for start values
    !!EqType is for quadruple line, iQ is for startpoint with desired phasefractions (plural!) = 0
    integer :: EqTypeStart, EqTypeEnd, iQ

    ! Mixed hydrates 2015
    double precision, dimension(30) :: xH
    double precision, dimension(3,30):: CiJ, occup, occup_single, occup_double
    double precision, dimension(3) :: Phasefrac

    ! Properties for enthalpy of formation and hydration number
    integer :: ii, J, refin_VHIw
    double precision :: h_vap, h_hyd, s_hyd, h_liq, h_melt, h_phase, h_sol, h_WaterIce
    double precision, dimension(30) :: fug_g, hyd_nr!, occup_ls, occup_ld, occup_sms, occup_smd

    integer :: errorflag, eqtype_ind

    !New variables for tpd calculation
    double precision:: tpd_sol, ChemPot_hyd, a0_other_struct
    integer:: hdrt_structure_flag
    !//////////////////////////////////////////////////////////////////

    end subroutine phaselines_3
    !******************************************************************************

    !******************************************************************************
    module subroutine phaseemerge_detection(pel, EqType, phases, iQ)
    !******************************************************************************
    !
    !


    implicit none
    type(type_hdrt_pheq) :: pel
    integer, intent(in) :: phases(2), EqType
    integer, intent(out) :: iQ

    end subroutine
    !******************************************************************************


    !******************************************************************************
    module subroutine threephase_equil_iter(gl, pel, Eqtype,cur_3ph, iTp,Temp, press, x, x_ph1, x_ph2, x_ph3, rho_ph1, rho_ph2, rho_ph3,  &
        & Chempot, lnfi, occup, error, &!Q_point, Qexist, Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4, &
        x_vap_given, bisect_mode)
    !******************************************************************************
    ! This subroutine iterates the overall composition to match the given ratio of
    ! molefractions for (atm) two components of the gasphase



    implicit none
    type(type_gl) :: gl
    type(type_hdrt_pheq) :: pel
    character(len=6), intent(in) :: Eqtype        !Type of Equilibrium, which should be calculated
    integer, intent(in) :: iTp

    !Variables at the three-phase equilibrium
    double precision :: Temp, press
    double precision, dimension(30) :: x, x_ph1, x_ph2, x_ph3, x_vap_given, x_lower, x_upper
    double precision :: rho_ph1, rho_ph2, rho_ph3, x_diff_exit, x_change, x_diff_prev
    double precision, dimension(30) :: Chempot, lnfi


    integer :: cur_3ph
    integer :: error, itermax, iterations
    integer :: iPhase, iFlash, bisect_mode

    !Variables at the Quadruple point
    !character(len=6), dimension(4), intent(in)  :: Q_point
    !double precision, dimension (4), intent(in) :: press_Q, Temp_Q
    !double precision, dimension(4,30), intent(in) :: x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4
    logical, dimension (4):: Qexist

    ! Mixed hydrates 2015
    double precision, dimension(3,30):: occup


    !//////////////////////////////////////////////////////////////////

    end subroutine threephase_equil_iter

    module subroutine phaselines_4(gl, phl4, x, x_vap_given, eqtype_ind, error)

    implicit none
    type(type_gl) :: gl
    type(type_hdrt_4ph_lines) :: phl4
    double precision, dimension(30) :: x, x_vap_given
    integer :: eqtype_ind
    integer :: error


    end subroutine phaselines_4

    module subroutine phaselines_4_move_q_point(phl4, i_in, pure_Q_flag)

    implicit none
    type(type_hdrt_4ph_lines) :: phl4
    integer :: pure_Q_flag,i,i_in

    end subroutine phaselines_4_move_q_point


    module subroutine massbalance_4ph(gl, press, press_prev, x, x_fluid1, x_fluid2, x_fluid3, x_sol, x_hyd, vec_compphafra, vec_compphafra_prev, eqtype_ind, inphysicalx, root_dirx, root_between, root_sysofeqs)

    implicit none
    type(type_gl) :: gl
    double precision :: press, press_prev, press_sol(4,3,2)
    double precision, dimension(60,60) :: mat_matbal
    double precision, dimension(60,30) :: vec_compphafra, vec_compphafra_prev
    double precision, dimension(30):: x, x_fluid1, x_fluid2, x_fluid3, x_sol, x_hyd
    integer :: eqtype_ind, inphysicalx, inphysical, nsysofeqs, betas_good
    character(255):: herr
    !integer variables for loops etc
    integer :: ierr, neq, l, ll, i, ii, j
    !sign_changed is a vector that indicates if if the sign of beta_ph_i of sysofeqs l has changed from previous to current step
    integer :: root_between, root_dirx, root_dir, root_dir_changed, sign_changed(4), root_sysofeqs

    end subroutine massbalance_4ph


    module subroutine GenerateFitFile(gl, pel, Filepath,Datenfilepath,typeofcalc,dT_ab_Q,dp_ab_Q, &
        &          Q_point, Qexist,Temp_Q, press_Q, x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4)

    implicit none

    type(type_gl) :: gl
    type(type_hdrt_pheq) :: pel
    character(len=255) :: Datenfilepath
    double precision, dimension(30) :: x_ph1, x_ph2, x_ph3
    integer :: typeofcalc
    character(len=255) :: Filepath
    double precision, dimension(2) ::dT_ab_Q, dp_ab_Q
    character(len=6), dimension(4), intent(in)  :: Q_point
    double precision, dimension (6), intent(in) :: press_Q, Temp_Q
    double precision, dimension(4,30), intent(in) :: x_Q_ph1, x_Q_ph2, x_Q_ph3, x_Q_ph4
    logical, dimension (4):: Qexist


    end subroutine GenerateFitFile
    !******************************************************************************
    !******************************************************************************






    module subroutine cpo_2phase_at_3phase (ph_eq_type, Temp, press, x, &
        &           Chempot_ph1, &
        &           lnfi_ph1, x_ph1, x_ph2, &
        &           x_ph3, rho_out, error)

    implicit none

    integer :: ph_eq_type

    double precision, dimension(30) :: x, x_ph1, x_ph2, x_ph3
    double precision, dimension(30) :: Chempot_ph1
    double precision, dimension(30) :: lnfi_ph1

    double precision:: Temp, press
    integer :: error

    double precision, dimension(5):: rho_out

    end subroutine cpo_2phase_at_3phase
    !******************************************************************************

    end interface

    end module hdrt_phaselines_calc

    !DEC$ ELSE
    module hdrt_phaselines_calc


    end module hdrt_phaselines_calc
    !DEC$ END IF