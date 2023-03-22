
    module plot_info
    type plot_i
        real(4) :: plot_x,plot_y
        character(50):: font_t
        double precision :: font_size
        integer:: resolution
    end type

    type(plot_i):: plt_i

    end module

    module curves

    !DEC$ IF DEFINED(IDEAL_CURVES)
    use module_all_types
    use calc_functions

    implicit none

    !class for calculation of the different ideal curves in parallel
    type curve
        !arrays with the results
        double precision, allocatable :: t_red(:),p_red(:)
        double precision :: d_last
        double precision :: t_last
        INTEGER :: Max_Iterations, Iterations
        DOUBLE PRECISION :: dstart_min,dstart_max
        DOUBLE PRECISION :: dmin_allowed,dmax_allowed
        DOUBLE PRECISION :: Delta_allowed
        INTEGER :: error
        type(type_additional_parameters) :: t_params
        double precision :: t_trip,t_step
    end type

    !**********************************************************************************************************************
    !interface for standard constructor
    interface curve
    module procedure init_curve
    end interface
    !**********************************************************************************************************************

    contains

    !----------------------------------------------------------------------------------------------------------------------
    !Constructor for all variables and checking the inputs
    !----------------------------------------------------------------------------------------------------------------------
    type(curve) function init_curve(gl,nr_pts,crv_type,nrsubst)
    integer:: nr_pts,nrsubst,crv_type,i
    double precision :: t_step
    type(type_gl) :: gl
    allocate(init_curve%t_red(nr_pts))
    allocate(init_curve%p_red(nr_pts))
    init_curve%d_last = 0d0
    init_curve%t_last = 0d0
    init_curve%Max_Iterations = 50
    init_curve%Iterations = 0

    !specific for each fluid
    init_curve%dstart_min = 1.d-4
    init_curve%dstart_max = gl%rhored(nrsubst)*3d0
    init_curve%dmin_allowed =  1.d-9
    init_curve%dmax_allowed = gl%rhored(nrsubst)*3d0
    init_curve%Delta_allowed = 1.d-9
    init_curve%error = 0
    !init_curve%param = 0d0


    if((crv_type.eq.1).or.(crv_type.eq.2)) then
        init_curve%t_red(1) = TBoyle_calc(gl,nrsubst)
    elseif(crv_type.eq.3) then
        init_curve%t_red(1) = TJTINV_calc(gl,nrsubst)
        if(init_curve%t_red(1) .eq. 0d0) then
            init_curve%t_red(1) = 4d0*gl%tc(nrsubst)
        end if
    elseif(crv_type.eq.4) then
        init_curve%t_red(1) = TJT_calc(gl,nrsubst)
        if(init_curve%t_red(1) .eq. 0d0) then
            init_curve%t_red(1) = 15d0*gl%tc(nrsubst)
        end if
    end if

    !triple point temperature
    if(crv_type.eq.1 .or. crv_type.eq.3 .or. crv_type.eq.4) then
        init_curve%t_trip = gl%ttp(nrsubst)
    else
        init_curve%t_trip = gl%tc(nrsubst)
    end if

    t_step = (init_curve%t_red(1) - init_curve%t_trip)/(nr_pts-1)
    do i = 2,nr_pts
        init_curve%t_red(i) = init_curve%t_red(1) - i*t_step
    end do

    !set the parameter
    init_curve%t_params%a_p(1) = 0.d0
    init_curve%t_params%a_p(2) = 0.d0
    init_curve%t_params%a_p(3) = nrsubst

    end function

    end module
    !-

    !**********************************************************************************************************************
    !**********************************************************************************************************************
    !**********************************************************************************************************************
    ! Module for the calculation of the ideal curves
    ! By Sven Pohl (2019)
    !**********************************************************************************************************************
    !**********************************************************************************************************************
    !**********************************************************************************************************************
    module ideal_curves

    use module_all_types
    use controlling
    use calc_functions
    use curves
    use plot_info
    use flash_pure_module
    use setup_module
    use calc_functions
    implicit none

    !**********************************************************************************************************************

    !interface for residual function
    abstract interface
    double precision function custom_name(gl,d_akt,parameters)
    use module_all_types
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION :: d_akt,parameters(65)
    end function custom_name
    end interface
    type container
        procedure(custom_name), pointer, nopass :: ptr !(gl,d_akt,parameters)
    end type
    type func_ptrs
        type(container), dimension(4):: fcn
    end type

    !class type definition
    type ideal
        type(c_ptr) :: gl_handle
        character(255) :: path
        integer, dimension (30) :: EOS_indicator
        character(30), dimension(30):: fluids
        integer:: mixtype
        double precision , dimension(30):: moles
        character (12) :: input
        character (20) :: unit
        type(curve), allocatable:: crv(:)
        type(func_ptrs):: res_func
        !for p_vap calculations
        double precision :: t_start,t_end,t_step
        integer:: t_stp_nr
        double precision , allocatable:: temps(:),press_v(:)
    contains

    !**********************************************************************************************************************
    !member functions of ideal curves
    procedure :: get_classic => ideal_curve

    end type

    !**********************************************************************************************************************
    !interface for standard constructor
    interface ideal
    module procedure  init_ideal
    end interface
    !**********************************************************************************************************************

    contains

    !----------------------------------------------------------------------------------------------------------------------
    !Constructor for all variables and checking the inputs
    !----------------------------------------------------------------------------------------------------------------------
    type(ideal) function init_ideal(gl,input_in,fluids_in,moles_in,EOS_indicator_in,path_in,unit_in,gl_handle_in,nrsubst)

    implicit none

    type(type_gl),pointer :: gl
    type(c_ptr) :: gl_handle_in
    character(255) :: path_in
    integer, dimension (30) :: EOS_indicator_in
    character(30), dimension(30):: fluids_in
    character (12) :: input_in
    integer:: errorflag,nrsubst,i
    double precision :: prop1,prop2
    character (12) :: unit
    double precision , dimension(30):: moles_in
    character (20) :: unit_in

    !----------------------------------------------------------------------------------------------------------------------
    !set structure values to in values
    init_ideal%input = input_in
    init_ideal%fluids = fluids_in
    init_ideal%EOS_indicator = EOS_indicator_in
    init_ideal%gl_handle = gl_handle_in
    init_ideal%moles = moles_in
    init_ideal%path = path_in

    !check the control fluids input
    call control_fluids(gl,init_ideal%input,init_ideal%fluids,init_ideal%EOS_indicator,init_ideal%gl_handle)

    !default settings for all calculations:
    gl%calc_ref = .false.
    gl%transport = 0
    gl%VLE_needed = .false.
    gl%converttype = 0
    gl%vir = .false.

    !call the setup
    init_ideal%mixtype = 1
    errorflag = 0

    !Check characters of unit and set intern unit parameter
    call uppertolower_char(unit_in,len(unit_in))

    if(trim(unit_in) == "specific") then
        gl%unitin=2
    elseif(trim(unit_in) == "molar") then
        gl%unitin=1
    elseif(trim(unit_in) == "reduced") then
        gl%unitin=3
    end if

    !write the unit in class structure
    init_ideal%unit = unit_in

    !set some dummy input
    init_ideal%input = input_in
    prop1= 200d0
    prop2= 1d0

    !initalize the fluid with the use of the setup
    call setup(gl,init_ideal%input, prop1, prop2, init_ideal%fluids, init_ideal%moles, init_ideal%path, init_ideal%EOS_indicator, init_ideal%mixtype, errorflag)

    !init ptrs for all functions
    init_ideal%res_func%fcn(1)%ptr => Res_ideal
    init_ideal%res_func%fcn(2)%ptr => Res_boyle
    init_ideal%res_func%fcn(3)%ptr => Res_JTIC
    init_ideal%res_func%fcn(4)%ptr => Res_JIC

    !init values for p_vap calculations
    if(plt_i%resolution.eq.1) then
        init_ideal%t_stp_nr = 100
    elseif(plt_i%resolution.eq.2) then
        init_ideal%t_stp_nr = 200
    elseif(plt_i%resolution.eq.2) then
        init_ideal%t_stp_nr = 500
    end if

    allocate(init_ideal%temps(init_ideal%t_stp_nr))
    allocate(init_ideal%press_v(init_ideal%t_stp_nr))

    init_ideal%t_end =  gl%ttp(nrsubst)
    init_ideal%t_start = gl%tc(nrsubst)
    init_ideal%t_step = (init_ideal%t_start-init_ideal%t_end)/init_ideal%t_stp_nr
    init_ideal%temps(1) = init_ideal%t_start
    do i=2,init_ideal%t_stp_nr-1
        init_ideal%temps(i) =init_ideal%t_start  - i*init_ideal%t_step
    end do

    end function

    !-------------------------------------------------------------------------------------------------
    !Span - Multiparameter equation of state -  page 168
    !Ideal curve: Condition: (Z=1) dalphar/ddelta = 0
    !Boyle curve: Condition: (dZ/drho) = 0
    !Joule-Thomson inv. curve: Condition: (dZ/dT)p = 0
    !Joulve inversion curve: Condition: (dZ/dT)rho = 0
    !-------------------------------------------------------------------------------------------------
    subroutine ideal_curve(id,gl,nrsubst)

    implicit none
    class(ideal):: id
    type(type_gl) :: gl
    double precision :: T,D
    integer:: nrsubst
    integer:: lc,loop_t
    integer:: steps

    !number of steps
    steps = 200

    !Create four curve objects for each ideal type
    if(.not.allocated(id%crv)) allocate(id%crv(4))
    !Initalize the starting values ofr each curve
    id%crv(1) = curve(gl,steps,1,nrsubst)
    id%crv(2) = curve(gl,steps,2,nrsubst)
    id%crv(3) = curve(gl,steps,3,nrsubst)
    id%crv(4) = curve(gl,steps,4,nrsubst)


    !loop over all ideal curves
    do lc=1,4
        do loop_t = 1,steps

            if(loop_t.ne.1) then !at the first point density is zero

                id%crv(lc)%t_params%a_P(2) = id%crv(lc)%t_red(loop_t)

                !D = de_root(gl,id%res_func%fcn(lc)%ptr,id%crv(lc)%dmin_allowed,id%crv(lc)%dmax_allowed,id%crv(lc)%param)
                !D = improved_regulafalsi(gl,id%res_func%fcn(lc)%ptr,id%crv(lc)%dmin_allowed,id%crv(lc)%dmax_allowed,id%crv(lc)%param)

                !D = 0d0
                !Iterate the density
                CALL Regula_Falsi(gl,id%res_func%fcn(lc)%ptr,D,id%crv(lc)%dstart_min,id%crv(lc)%dstart_max,id%crv(lc)%Delta_allowed,&
                    id%crv(lc)%dmin_allowed,id%crv(lc)%dmax_allowed,id%crv(lc)%Max_iterations,id%crv(lc)%Iterations,id%crv(lc)%error,id%crv(lc)%t_params)

                !calculate the pressure
                id%crv(lc)%p_red(loop_t) = P_CALC(gl,id%crv(lc)%t_red(loop_t),D,nrsubst)/P_CALC(gl,gl%tred(nrsubst),gl%rhored(nrsubst),nrsubst)

                id%crv(lc)%t_red(loop_t) = id%crv(lc)%t_red(loop_t)/gl%tred(nrsubst)

                if(D.ne.0d0) then
                    id%crv(lc)%dstart_min = D
                end if
            else
                !at starting point
                D = 1d-12
                id%crv(lc)%p_red(loop_t) = P_CALC(gl,id%crv(lc)%t_red(loop_t),D,nrsubst)/P_CALC(gl,gl%tred(nrsubst),gl%rhored(nrsubst),nrsubst)
                !id%crv(lc)%t_red(loop_t+1) = id%crv(lc)%t_red(loop_t)-2d0
                id%crv(lc)%t_params%a_p(2) = id%crv(lc)%t_red(loop_t)
                id%crv(lc)%t_red(loop_t) = id%crv(lc)%t_red(loop_t)/gl%tred(nrsubst)
            end if
        end do
    end do

    !calculate p_vap
    call calc_p_vap(gl,id,nrsubst)

    !clean the results - convert all negative values to zero - then write file to plot the curves
    call write_data_for_plot(id,steps)

    !ideal_curve = 1

    end subroutine

    !****************************************************************************************************************************************
    !****************************************************************************************************************************************
    !RESIDUAL FUNCTIONS FOR ALL TYPE OF IDEAL CURVES
    !****************************************************************************************************************************************
    !****************************************************************************************************************************************
    !--------------------------------------------------------------------------
    !Ideal curve
    DOUBLE PRECISION FUNCTION Res_ideal(gl,d_akt,parameters)
    !--------------------------------------------------------------------------
    !Parameter(1) not used because a root is searched: goal is zero!
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION :: d_akt,parameters(65)
    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER
    GETDERR = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/) !Define derivatives needed
    CALL FNRDERIVS(gl,parameters(2),d_akt,GETDERR,FNRDER,int(parameters(3)))
    Res_ideal= FNRDER(2)
    END function
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !Boyle curve
    DOUBLE PRECISION FUNCTION Res_boyle(gl,d_akt,parameters)
    !--------------------------------------------------------------------------
    !Parameter(1) not used because a root is searched: goal is zero!
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION :: d_akt,parameters(65)
    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER
    GETDERR = (/0,1,1,0,0,0,0,0,0,0,0,0,0,0,0/)
    CALL FNRDERIVS(gl,parameters(2),d_akt,GETDERR,FNRDER,int(parameters(3)))
    Res_boyle = FNRDER(2)+FNRDER(3)
    END function
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    !Joule-Thomson inv. curve
    DOUBLE PRECISION FUNCTION Res_JTIC(gl,d_akt,parameters)
    !--------------------------------------------------------------------------
    !Parameter(1) not used because a root is searched: goal is zero!
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION :: d_akt,parameters(65)
    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER
    GETDERR = (/0,1,1,0,0,1,0,0,0,0,0,0,0,0,0/)
    CALL FNRDERIVS(gl,parameters(2),d_akt,GETDERR,FNRDER,int(parameters(3)))
    Res_JTIC = FNRDER(2)+FNRDER(3)+FNRDER(6)
    END function
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !Joulve inversion curve
    DOUBLE PRECISION FUNCTION Res_JIC(gl,d_akt,parameters)
    !--------------------------------------------------------------------------
    !Parameter(1) not used because a root is searched: goal is zero!
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION :: d_akt,parameters(65)
    integer,DIMENSION(nderivs) :: GETDERR
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER
    GETDERR = (/0,1,1,0,0,1,0,0,0,0,0,0,0,0,0/)
    CALL FNRDERIVS(gl,parameters(2),d_akt,GETDERR,FNRDER,int(parameters(3)))
    Res_JIC = FNRDER(6)
    END function
    !--------------------------------------------------------------------------

    !calculation of p_vap for boyle curve
    subroutine calc_p_vap(gl,id,nrsubst)
    implicit none
    type(type_gl) :: gl
    class(ideal):: id
    integer:: i,nrsubst,iter,iFlash,errval
    double precision :: rhovap_est,rholiq_est
    iFlash = 1 ! pvap
    rhovap_est = 0d0
    rholiq_est = 0d0
    errval = 0

    id%press_v(1) = P_CALC(gl,gl%tred(nrsubst),gl%rhored(nrsubst),nrsubst)
    id%temps(1) = gl%tred(nrsubst)
    do i=2,size(id%temps)
        call Flash_Pure_PhaseBoundary (gl,id%press_v(i),id%temps(i), rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst)
    end do

    id%press_v = id%press_v /P_CALC(gl,gl%tred(nrsubst),gl%rhored(nrsubst),nrsubst)
    id%temps = id%temps/gl%tred(nrsubst)

    end subroutine calc_p_vap





    !write data to file - (gnuplot format and plot it)
    subroutine write_data_for_plot(id,steps)
    implicit none
    class(ideal):: id
    integer:: steps,i,j
    integer:: loc_classic_ideal,loc_boyle,loc_jti,loc_ji
    double precision :: x_ci_l,y_ci_l !label coordinates
    double precision :: x_b_l,y_b_l !label coordinates
    double precision :: x_ji_l,y_ji_l !label coordinates
    double precision :: x_j_l,y_j_l !label coordinates

    ! Body of write_data_for_plot
    open(unit=60,file='ic_vals.dat')
    write(60,*) "##########################################################"
    write(60,*) "################# IDEAL CURVES ###########################"
    write(60,*) "##########################################################"
    write(60,*) "# t/tc - p/pc"
    write(60,*) "##########################################################"
    write(60,*) ""

    do j=1,4

        if(j.eq.1) then
            write(60,*) '# CLASSIC IDEAL CURVE '
        end if

        do i=1,steps
            if(id%crv(j)%p_red(i).gt.0d0) then
                write(60,*) id%crv(j)%t_red(i), id%crv(j)%p_red(i) , j
            end if
        end do

        write(60,*) ""
        if(j.eq.1) then
            write(60,*) '# BOYLE CURVE '
        elseif(j.eq.2) then
            write(60,*) '# JOULE-THOMSON INVERSION CURVE '
        elseif(j.eq.3) then
            write(60,*) '# JOULE INVERSION CURVE '
        end if
        write(60,*) ""
    end do

    !write pv data to file
    write(60,*) ''
    do i=1,size(id%press_v)-1
        write(60,*) id%temps(i), id%press_v(i)
    end do

    close(60)

    !find positions for the labels: depends on shape of ideal curves
    !first classic ideal curve
    !Max loc of p/pc
    loc_classic_ideal = steps/2!maxloc(id%crv(1)%p_red,1)
    loc_boyle = maxloc(id%crv(2)%p_red,1)
    loc_jti= maxloc(id%crv(3)%p_red,1)
    loc_ji= maxloc(id%crv(4)%p_red,1)

    x_ci_l = id%crv(1)%t_red(1)
    y_ci_l = id%crv(1)%p_red(loc_classic_ideal)
    x_b_l = id%crv(2)%t_red(loc_boyle)
    y_b_l = id%crv(2)%p_red(loc_boyle)
    x_ji_l = id%crv(3)%t_red(loc_jti)
    y_ji_l = id%crv(3)%p_red(loc_jti)
    x_j_l = id%crv(4)%t_red(loc_ji)
    y_j_l= id%crv(4)%p_red(loc_ji)


    open(unit=60,file='plot_ideal.plt')
    write(60,*) 'set term pdfcairo font', '"', trim(plt_i%font_t) , ",",int(plt_i%font_size),'" \'
    write(60,*) " size",  plt_i%plot_x ,"cm,",plt_i%plot_y ,"cm"
    write(60,*) 'set output "ideal_c.pdf"'
    write(60,*) 'set nokey'

    write(60,*) 'set label 1  "CLASSIC IDEAL" at', int(x_ci_l),",",int(y_ci_l), 'font ",8"'
    write(60,*) 'set label 2  "BOYLE CURVE" at', 1.02,",",1,'font ",8"'
    write(60,*) 'set label 3  "JOULE-THOMSON INVERSION" at', int(x_ji_l),",",int(y_ji_l)+5,'font ",8"'
    write(60,*) 'set label 4  "JOULE INVERSION" at', int(x_j_l)+5,",",int(y_j_l)+20,'font ",8"'

    write(60,'(A)') 'set object circle at first  1,1 radius char 0.25 fillcolor rgb "black" fillstyle solid noborder'
    write(60,*) 'set linestyle  1 lt 1 lw 1 lc black'
    write(60,*) 'set logscale xy'
    write(60,*) 'set yrange [0.1:400]'
    write(60,*) 'set xrange [0.5:20]'
    write(60,*) 'set xlabel("{/:Italic T} / {/:Italic T}_{c}")'
    write(60,*) 'set ylabel("{/:Italic p} / {/:Italic p}_{c}")'
    write(60,*) 'plot "ic_vals.dat" u 1:2 w lines ls 1'

    close(60)

    Call SYSTEM('gnuplot plot_ideal.plt')

    end subroutine write_data_for_plot
    !DEC$ ELSE
    contains
    subroutine dummy()
    end subroutine
    !DEC$ ENDIF
    end module


    !  !differential evolution to find the root
    !double precision function improved_regulafalsi(gl,func,val_start,val_end,parameters)
    !! Variables
    !use module_all_types
    !implicit none
    !type(type_gl) :: gl
    !Double Precision :: func
    !DOUBLE PRECISION :: parameters(65)
    !double precision :: val_start,val_end,x1,x2
    !double precision :: fz,f1,f2
    !double precision :: z,s,m
    !integer:: n_iter
    !fz=1d3
    !m=0.5d0
    !!get the starting values
    !x1 = val_start
    !x2 = val_end
    !n_iter = 0
    !do while(dabs(x2-x1).gt.1D-7.AND. dabs(fz).gt.1.d-7)
    !
    !    f1 = func(gl,x1,parameters)
    !    f2 = func(gl,x2,parameters)
    !    if((f1*f2).gt.0.d0) then
    !        do while((f1*f2).gt.0d0 .and. (x1 .lt. x2))
    !            x1= x1 + 1D-1
    !            f1 = func(gl,x1,parameters)
    !            f2 = func(gl,x2,parameters)
    !        end do
    !        improved_regulafalsi = 0d0
    !    end if
    !
    !    s = (f2-f1)/(x2-x1)
    !    z= x1-func(gl,x1,parameters)/s
    !    fz = func(gl,z,parameters)
    !
    !    !decision block
    !    if(fz*f2 .lt. 0d0) then
    !        x1 = x2
    !        f1 = f2
    !        x2 = z
    !        f2 = fz
    !    elseif(fz*f2 .gt. 0d0) then
    !        f1 = 0.5d0*f1
    !        x2 = z
    !        f2 = fz
    !    else
    !        x1 = z
    !        f1 = fz
    !        x2 = z
    !        f2 = fz
    !    end if
    !
    !    n_iter = n_iter +1
    !
    !    if(n_iter.gt.100) then
    !        exit
    !    end if
    !
    !end do
    !
    !improved_regulafalsi = z
    !n_iter  =0
    !end function
    !
    !
    !double precision function de_root(gl,func,val_start,val_end,parameters)
    !! Variables
    !use module_all_types
    !implicit none
    !type(type_gl) :: gl
    !
    !Double Precision :: func
    !DOUBLE PRECISION :: parameters(65)
    !double precision :: val_start,val_end
    !double precision, allocatable:: val_sols(:),val(:),val_old(:),dev(:),dev_old(:),R(:)
    !double precision, allocatable:: r_pos(:,:)
    !double precision :: min_dev,F,C
    !integer:: nr_sols,i,m,rd_pos
    !double precision , dimension(3):: random_dbl
    !integer, dimension(3):: random_ints
    !nr_sols = 500
    !F = 0.01d0
    !C = 0.1d0
    !
    !allocate(val_sols(nr_sols),val(nr_sols),val_old(nr_sols),dev(nr_sols),dev_old(nr_sols),r(nr_sols))
    !allocate(r_pos(nr_sols,3))
    !call RANDOM_SEED()
    !call RANDOM_NUMBER(R)
    !
    !
    !!create start solution
    !do i=1,nr_sols
    !    val(i) = val_start+R(i)*(val_end-val_start)
    !    dev(i) = func(gl,val(i),parameters)
    !end do
    !
    !val_old = val
    !dev_old = dev
    !
    !min_dev = (minval(dev))
    !
    !!loop generations
    !do while((min_dev).gt.1d-8)
    !
    !    do i=1,nr_sols
    !        call RANDOM_NUMBER(random_dbl)
    !        random_ints = ceiling(random_dbl* nr_sols)
    !        val(i) = val(random_ints(1))*C + F*(val(random_ints(2))-val(random_ints(3)))
    !        !if(val(i).lt.0d0) val(i) = val(i)*-1d0
    !    end do
    !
    !    !calc new residual
    !    do i=1,nr_sols
    !        dev(i) = func(gl,val(i),parameters)
    !        if(dev(i).lt.dev_old(i)) then
    !            val_old(i) = val(i)
    !            dev_old(i) = dev(i)
    !        else
    !            val(i) = val_old(i)
    !        end if
    !    end do
    !
    !    min_dev =  (minval(dev_old))
    !    write(*,*) min_dev
    !end do
    !
    !de_root = val_old(minloc(dev_old,1))
    !
    !end function de_root