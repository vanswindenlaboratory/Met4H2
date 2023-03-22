    module interface_fitter

    use module_all_types
    use calc_functions
    use controlling
    use rhomix_pt_module
    use flash_module
    use flash_pure_module

    implicit none

    contains

    ! Set the helmholtz parameters from extern
    subroutine set_helmholtz_param_from_x(gl,x,n_i,d_i,t_i,p_i,gam_exp,gauss,fit_mode,nterm)
    type(type_gl) :: gl
    integer :: fit_mode,nterm
    double precision , dimension(:) :: x
    integer, dimension(:) ::  d_i,p_i
    double precision, dimension(:) :: t_i,n_i
    double precision, dimension(:,:) :: gauss
    double precision , dimension(:) :: gam_exp
    integer :: i,nr_gam_exp_fitted,n_exp,n_gauss,nr_gam_exp

    gl%eos_coeff%p_i(:,1) =  0
    gl%eos_coeff%gama = 0d0
    n_exp = gl%eos_coeff%nreg(1)-gl%eos_coeff%I_pol(1)
    n_gauss = gl%eos_coeff%I_gbs(1)

    if( any(gam_exp .ne. 1d0)  ) then
        nr_gam_exp_fitted = size(gam_exp,1)
    else
        nr_gam_exp_fitted = 0
    end if

    if(fit_mode .eq. 1) then      ! all parameters of the function are fitted

        ! coefficients position: 1:nterm
        gl%eos_coeff%ni(1:nterm,1) = x(1:nterm)

        ! temperature exponents position: nterm+1:2*nterm
        gl%eos_coeff%ti(1:nterm,1) = x(1+nterm:2*nterm)

        ! density exponents position: 2*nterm+1:3nterm
        gl%eos_coeff%di(1:nterm,1) = x(2*nterm+1:3*nterm)

        ! p_i from x vector 3*nterm+1:3*nterm+n_exp
        gl%eos_coeff%p_i(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) =  x(3*nterm+1:3*nterm+n_exp)

        ! gama_i parameter position: 3*nterm+1+n_exp:3*nterm+2*n_exp
        gl%eos_coeff%gama(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) =  x(3*nterm+1+n_exp:3*nterm+2*n_exp)

        !open(90, file='testdatei.txt')
        !write(90,*) gl%eos_coeff%p_i(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1)
        !close(90)

        ! gaussian parameters position: 3*nterm+2*n_exp+a*n_gauss:3*nterm+2*n_exp+(a+1)*n_gauss
        if(gl%eos_coeff%I_gbs(1) .gt. 0) then
            gl%eos_coeff%pli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
            gl%eos_coeff%tli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
            gl%eos_coeff%eta(1:gl%eos_coeff%I_gbs(1),1) =  x(3*nterm+2*n_exp+1:3*nterm+2*n_exp+n_gauss)
            gl%eos_coeff%beta(1:gl%eos_coeff%I_gbs(1),1) = x(3*nterm+2*n_exp+1+n_gauss:3*nterm+2*n_exp+2*n_gauss)
            gl%eos_coeff%gam(1:gl%eos_coeff%I_gbs(1),1) =  x(3*nterm+2*n_exp+1+2*n_gauss:3*nterm+2*n_exp+3*n_gauss)
            gl%eos_coeff%eps(1:gl%eos_coeff%I_gbs(1),1) =  x(3*nterm+2*n_exp+1+3*n_gauss:3*nterm+2*n_exp+4*n_gauss)
        end if


    elseif(fit_mode == 2) then ! all but no density exponents
        gl%eos_coeff%ni(1:nterm,1) = x(1:nterm)
        gl%eos_coeff%ti(1:nterm,1) = x(1+nterm:2*nterm)

        if(gl%eos_coeff%I_gbs(1) .gt. 0) then
            gl%eos_coeff%pli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
            gl%eos_coeff%tli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
            gl%eos_coeff%eta(1:gl%eos_coeff%I_gbs(1),1) =  x(2*nterm+1+nr_gam_exp_fitted:2*nterm+gl%eos_coeff%I_gbs(1)+nr_gam_exp_fitted)
            gl%eos_coeff%beta(1:gl%eos_coeff%I_gbs(1),1) = x(2*nterm+1+gl%eos_coeff%I_gbs(1)+nr_gam_exp_fitted:2*nterm+2*gl%eos_coeff%I_gbs(1)+nr_gam_exp_fitted)
            gl%eos_coeff%gam(1:gl%eos_coeff%I_gbs(1),1) =  x(2*nterm+1+2*gl%eos_coeff%I_gbs(1)+nr_gam_exp_fitted:2*nterm+3*gl%eos_coeff%I_gbs(1)+nr_gam_exp_fitted)
            gl%eos_coeff%eps(1:gl%eos_coeff%I_gbs(1),1) =  x(2*nterm+1+3*gl%eos_coeff%I_gbs(1)+nr_gam_exp_fitted:2*nterm+4*gl%eos_coeff%I_gbs(1)+nr_gam_exp_fitted)
        end if
        gl%eos_coeff%di(1:nterm,1) = d_i(1:nterm)


        if( nr_gam_exp_fitted .ne. 0)  then
            gl%eos_coeff%gama(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) =  x(2*nterm+1:2*nterm+nr_gam_exp_fitted)
        else
            gl%eos_coeff%gama(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) = gam_exp(1:nr_gam_exp)
        end if

    elseif(fit_mode .eq. 3) then !only n are fitted
        gl%eos_coeff%ni(1:nterm,1) = x(1:nterm)
        gl%eos_coeff%ti(1:nterm,1) = t_i
        gl%eos_coeff%pli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
        gl%eos_coeff%tli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
        gl%eos_coeff%eta(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,1)
        gl%eos_coeff%beta(1:gl%eos_coeff%I_gbs(1),1) = gauss(:,2)
        gl%eos_coeff%gam(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,3)
        gl%eos_coeff%eps(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,4)
        gl%eos_coeff%di(1:nterm,1) = d_i(1:nterm)
        gl%eos_coeff%gama(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) = gam_exp(1:nr_gam_exp)
    elseif(fit_mode .eq. 4) then !only t are fitted
        gl%eos_coeff%ni(1:nterm,1) =  n_i
        gl%eos_coeff%ti(1:nterm,1) =  x(1:nterm)
        gl%eos_coeff%pli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
        gl%eos_coeff%tli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
        gl%eos_coeff%eta(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,1)
        gl%eos_coeff%beta(1:gl%eos_coeff%I_gbs(1),1) = gauss(:,2)
        gl%eos_coeff%gam(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,3)
        gl%eos_coeff%eps(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,4)
        gl%eos_coeff%di(1:nterm,1) = d_i(1:nterm)
        gl%eos_coeff%gama(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) = gam_exp(1:nr_gam_exp)
    elseif(fit_mode .eq. 5) then ! only t and n are fitted
        gl%eos_coeff%ni(1:nterm,1) =  x(1:nterm)
        gl%eos_coeff%ti(1:nterm,1) =  x(nterm+1:2*nterm)
        gl%eos_coeff%pli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
        gl%eos_coeff%tli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
        gl%eos_coeff%eta(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,1)
        gl%eos_coeff%beta(1:gl%eos_coeff%I_gbs(1),1) = gauss(:,2)
        gl%eos_coeff%gam(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,3)
        gl%eos_coeff%eps(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,4)
        gl%eos_coeff%di(1:nterm,1) = d_i(1:nterm)
        gl%eos_coeff%gama(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) = gam_exp(1:nr_gam_exp)
    elseif( fit_mode .eq. 6 ) then ! only gaussians are fitted
        gl%eos_coeff%ni(1:nterm,1) =  n_i
        gl%eos_coeff%ti(1:nterm,1) =  t_i
        gl%eos_coeff%pli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
        gl%eos_coeff%tli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
        gl%eos_coeff%eta(1:gl%eos_coeff%I_gbs(1),1) =  x(1:gl%eos_coeff%I_gbs(1))
        gl%eos_coeff%beta(1:gl%eos_coeff%I_gbs(1),1) = x(gl%eos_coeff%I_gbs(1)+1:2*gl%eos_coeff%I_gbs(1))
        gl%eos_coeff%gam(1:gl%eos_coeff%I_gbs(1),1) =  x(2*gl%eos_coeff%I_gbs(1)+1:3*gl%eos_coeff%I_gbs(1))
        gl%eos_coeff%eps(1:gl%eos_coeff%I_gbs(1),1) =  x(3*gl%eos_coeff%I_gbs(1)+1:4*gl%eos_coeff%I_gbs(1))
        gl%eos_coeff%di(1:nterm,1) = d_i(1:nterm)
        gl%eos_coeff%gama(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) = gam_exp(1:nr_gam_exp)
    end if



    !gl%eos_coeff%p_i(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) =  p_i(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1))


    !gl%tred(1) = x(size(x,1)-1)
    !gl%rhored(1) = x(size(x,1))
    !gl%tc(1) = gl%tred(1)
    !gl%accen = 0d0
    !gl%rhoredmix = gl%rhored(1)
    !gl%tredmix =  gl%tred(1)


    end subroutine set_helmholtz_param_from_x

    !************************************************************************************************************************
    !Set TREND parameters from extern (MATLAB interface)
    subroutine  set_parameter_from_extern_x(gl,x,n_i,t_i,d_i,p_i,gam_exp,gauss,add_vals,nterms,nrsubst,fit_mode,fit_red)
    ! Variables
    use module_all_types
    implicit none
    type(type_gl):: gl
    double precision, dimension(:) :: x
    double precision, dimension(2):: add_vals
    integer, dimension(:):: nterms,d_i,p_i
    double precision,dimension(:):: n_i,t_i
    double precision,dimension(:,:)::gauss
    double precision , dimension(:) :: gam_exp
    integer:: i,nrsubst,nterm,fit_mode,fit_red,n_red_param
    integer*4 :: Int_4



    if(fit_red .eq. 1) then
        n_red_param = 2
    else
        n_red_param = 0
    end if

    ! Body of set_parameter_from_extern_x
    gl%EQ_type = 0
    gl%components(1) = 'ljf'
    gl%EQ_type(1) =  1 !has to be changed later!!!
    gl%mix_Type = 1
    gl%Factor = 1.D0
    gl%factorpress = 1.D6
    gl%factortrans=1.D0
    gl%factorrbwr=1.D0
    gl%factor_VS5eta=1.D0
    gl%factor_VS5slj=1.D0
    gl%ncomp = 1
    gl%REQ(1)  = 1
    gl%ttp(1) = 0.1d0
    gl%tred(1) = add_vals(1)
    gl%rhored(1) = add_vals(2)
    gl%tc(1) = gl%tred(1)
    gl%rhoc(1) =  gl%rhored(1)
    gl%accen = 0d0
    gl%rhoredmix = gl%rhored(1)
    gl%tredmix =  gl%tred(1)
    if(.not. allocated(gl%eos_coeff)) allocate(gl%eos_coeff)
    gl%eos_coeff%I_pol(1) = nterms(1)
    gl%eos_coeff%I_exp(1) = nterms(2)
    gl%eos_coeff%I_gbs(1) = nterms(3)
    gl%eos_coeff%nreg(1) =  sum(nterms(1:2))
    nterm  = sum(nterms(1:3))

    if(.not.allocated(gl%eos_coeff%ni)) then
        allocate(gl%eos_coeff%ni(nterm,1))
        allocate(gl%eos_coeff%ti(nterm,1))
        allocate(gl%eos_coeff%di(nterm,1))
        allocate(gl%eos_coeff%p_i(nterm,1))
        allocate(gl%eos_coeff%gama(nterm,1))
        allocate(gl%eos_coeff%pli(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%tli(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%eta(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%beta(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%gam(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%eps(gl%eos_coeff%I_gbs(1),1))
    end if


    ! set the helmholtz parameters
    call set_helmholtz_param_from_x(gl,x,n_i,d_i,t_i,p_i,gam_exp,gauss,fit_mode,nterm)


    end subroutine set_parameter_from_extern_x


    subroutine set_parameter_from_arrays(gl,n_i,t_i,d_i,p_i,gam_exp,gauss,add_vals,nterms)
    implicit none
    double precision,dimension(:):: n_i,t_i
    double precision, dimension(:) :: d_i,p_i
    double precision,dimension(:,:)::gauss
    double precision , dimension(:) :: gam_exp
    integer, dimension(:):: nterms
    double precision, dimension(2):: add_vals
    integer :: nterm,nr_gam_exp
    type(type_gl) :: gl

    nterm = sum(nterms)
    nr_gam_exp = size(gam_exp,1)


    ! Body of set_parameter_from_extern_x
    gl%EQ_type = 0
    gl%components(1) = 'ljf'
    gl%EQ_type(1) =  1 !has to be changed later!!!
    gl%mix_Type = 1
    gl%Factor = 1.D0
    gl%factorpress = 1.D6
    gl%factortrans=1.D0
    gl%factorrbwr=1.D0
    gl%factor_VS5eta=1.D0
    gl%factor_VS5slj=1.D0
    gl%ncomp = 1
    gl%REQ(1)  = 1
    gl%ttp(1) = 0.1d0
    gl%tred(1) = add_vals(1)
    gl%rhored(1) = add_vals(2)
    gl%tc(1) = gl%tred(1)
    gl%rhoc(1) =  gl%rhored(1)
    gl%accen = 0d0
    gl%rhoredmix = gl%rhored(1)
    gl%tredmix =  gl%tred(1)
    if(.not. allocated(gl%eos_coeff)) allocate(gl%eos_coeff)
    gl%eos_coeff%I_pol(1) = nterms(1)
    gl%eos_coeff%I_exp(1) = nterms(2)
    gl%eos_coeff%I_gbs(1) = nterms(3)
    gl%eos_coeff%nreg(1) =  sum(nterms(1:2))
    nterm  = sum(nterms(1:3))

    if(.not.allocated(gl%eos_coeff%ni)) then
        allocate(gl%eos_coeff%ni(nterm,1))
        allocate(gl%eos_coeff%ti(nterm,1))
        allocate(gl%eos_coeff%di(nterm,1))
        allocate(gl%eos_coeff%p_i(nterm,1))
        allocate(gl%eos_coeff%gama(nterm,1))
        allocate(gl%eos_coeff%pli(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%tli(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%eta(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%beta(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%gam(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%eps(gl%eos_coeff%I_gbs(1),1))
    end if

    gl%eos_coeff%p_i(:,1) =  0
    gl%eos_coeff%gama = 0d0
    gl%eos_coeff%ni(1:nterm,1) = n_i
    gl%eos_coeff%ti(1:nterm,1) = t_i
    gl%eos_coeff%di(1:nterm,1) = d_i

    if(gl%eos_coeff%I_gbs(1) > 0) then
        gl%eos_coeff%pli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
        gl%eos_coeff%tli(1:gl%eos_coeff%I_gbs(1),1) = 2d0
        gl%eos_coeff%eta(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,1)
        gl%eos_coeff%beta(1:gl%eos_coeff%I_gbs(1),1) = gauss(:,2)
        gl%eos_coeff%gam(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,3)
        gl%eos_coeff%eps(1:gl%eos_coeff%I_gbs(1),1) =  gauss(:,4)
    end if

    gl%eos_coeff%p_i(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) =  p_i(1:nr_gam_exp)
    gl%eos_coeff%gama(gl%eos_coeff%I_pol(1)+1:gl%eos_coeff%nreg(1),1) = gam_exp(1:nr_gam_exp)

    end subroutine


    !************************************************************************************************************************
    ! Function to calculate the gradient of the sum of square function (objective function)
    function calc_gradient(gl,datas,datas_id,val_calc,sizes,num_virial_mode,nrsubst,x,n_i,d_i,t_i,p_i,gam_exp,gauss,fit_mode,nterm)

    type(type_gl) :: gl
    double precision , dimension(:) :: x
    integer, dimension(:) ::  d_i,p_i
    double precision, dimension(:) :: t_i,n_i
    double precision, dimension(:,:) :: gauss
    double precision , dimension(size(x)) :: x_tmp,calc_gradient
    double precision, dimension(:,:) :: datas
    integer , dimensioN(:) :: datas_id
    integer :: fit_mode
    double precision :: val_calc
    double precision , dimension(size(datas,1)) :: val_calc_dummy
    double precision, dimension(size(x,1)) :: val_calc_forward
    double precision , dimension(:) :: gam_exp
    integer, dimension(:) :: sizes
    integer :: i,nrsubst,nterm
    logical:: num_virial_mode ! if = true: numerical approximation of virial coeffs because the d_i are no integers (fit process)
    double precision , parameter :: d_num = 1D-8
    ! First stencial was already calculated during sum of square calculation (reuse it!!!!)

    ! save x
    x_tmp = x

    !loop over all x_s
    do i=1,size(x)

        x(i) = x(i) + d_num

        call set_helmholtz_param_from_x(gl,x,n_i,d_i,t_i,p_i,gam_exp,gauss,fit_mode,nterm)

        val_calc_forward(i) = calc_ssq(gl,x,datas,datas_id,val_calc_dummy,sizes,num_virial_mode,fit_mode,nterm,nrsubst)

        x = x_tmp
    end do

    do i=1,size(x)
        calc_gradient(i) = (val_calc_forward(i) - val_calc)/d_num
    end do


    end function calc_gradient
    !************************************************************************************************************************


    !************************************************************************************************************************
    ! Function to calculate the actual sum of squares
    double precision function calc_ssq(gl,x,datas,datas_id,val_calc,sizes,num_virial_mode,fit_mode,nterm,nrsubst)

    type(type_gl) :: gl
    double precision, dimension(:,:) :: datas
    double precision , dimension(:) :: x
    integer , dimensioN(:) :: datas_id
    double precision, dimension(:) :: val_calc
    integer, dimension(:) :: sizes
    integer :: i,nrsubst,fit_mode
    logical:: num_virial_mode ! if = true: numerical approximation of virial coeffs because the d_i are no integers (fit process)
    double precision :: delta_loss = 1d-2
    double precision :: d_tmp(nterm)
    integer :: nterm


    do i=1,sizes(3)
        ! if density exponents are fitted and the data id  ==  12,41 Virials
        ! then the exponents will be rounded
        if(datas(i,5) .ne. 0d0) then

            if( (datas_id(i) .eq. 12 .or. datas_id(i) .eq. 41) .and. fit_mode   == 1) then
                d_tmp  = x(2*nterm+1:3*nterm)
                gl%eos_coeff%di(1:nterm,1) =  floor (d_tmp)
            end if


            if( datas_id(i) == 19) then
                val_calc(i) = (datas(i,3)-A00_calc(gl,datas(i,1),datas(i,2),nrsubst))**2/datas(i,4)*datas(i,5)
            elseif(datas_id(i) == 20) then
                val_calc(i) = (datas(i,3)-A01_calc(gl,datas(i,1),datas(i,2),nrsubst))**2/datas(i,4)*datas(i,5)
            elseif(datas_id(i) == 21) then
                val_calc(i) = (datas(i,3)-A02_calc(gl,datas(i,1),datas(i,2),nrsubst))**2/datas(i,4)*datas(i,5)
            elseif(datas_id(i) == 22) then
                val_calc(i) = (datas(i,3)-A03_calc(gl,datas(i,1),datas(i,2),nrsubst))**2/datas(i,4)*datas(i,5)
            elseif( datas_id(i) == 23) then
                val_calc(i) = (datas(i,3)-A10_calc(gl,datas(i,1),datas(i,2),nrsubst))**2/datas(i,4)*datas(i,5)
            elseif( datas_id(i) == 24) then
                val_calc(i) = (datas(i,3)-A11_calc(gl,datas(i,1),datas(i,2),nrsubst))**2/datas(i,4)*datas(i,5)
            elseif( datas_id(i) == 25) then
                val_calc(i) = (datas(i,3)-A11_calc(gl,datas(i,1),datas(i,2),nrsubst))**2/datas(i,4)*datas(i,5)
            elseif( datas_id(i) == 26) then
                val_calc(i) = (datas(i,3)-A20_calc(gl,datas(i,1),datas(i,2),nrsubst))**2/datas(i,4)*datas(i,5)
            elseif( datas_id(i) == 27) then
                val_calc(i) = (datas(i,3)-A21_calc(gl,datas(i,1),datas(i,2),nrsubst))**2/datas(i,4)*datas(i,5)
            elseif( datas_id(i) == 28) then
                val_calc(i) = (datas(i,3)-A30_calc(gl,datas(i,1),datas(i,2),nrsubst))**2/datas(i,4)*datas(i,5)
            elseif( datas_id(i) == 12) then
                if(.not.num_virial_mode) then
                    val_calc(i) = (datas(i,3)-B_CALC(gl,datas(i,1),nrsubst))**2/datas(i,4)*datas(i,5)
                else
                    val_calc(i) = (datas(i,3)-B_CALC_num(gl,datas(i,1),nrsubst))**2/datas(i,4)*datas(i,5)
                end if
            elseif( datas_id(i) == 13) then
                if(.not.num_virial_mode) then
                    val_calc(i) = (datas(i,3)-C_CALC(gl,datas(i,1),nrsubst))**2/datas(i,4)*datas(i,5)
                else
                    val_calc(i) = (datas(i,3)-C_CALC_num(gl,datas(i,1),nrsubst))**2/datas(i,4)*datas(i,5)
                end if

            elseif( datas_id(i) == 41) then
                if(.not.num_virial_mode) then
                    val_calc(i) = (datas(i,3)-D_CALC(gl,datas(i,1),nrsubst))**2/datas(i,4)*datas(i,5)
                else
                    val_calc(i) = (datas(i,3)-D_CALC_num(gl,datas(i,1),nrsubst))**2/datas(i,4)*datas(i,5)
                end if
            end if

            ! set d_i back to float numbers
            if( (datas_id(i) .eq. 12 .or. datas_id(i) .eq. 41) .and. fit_mode   == 1) then
                gl%eos_coeff%di(1:nterm,1) = d_tmp
            end if
            !val_calc(i) = delta_loss**2  * (sqrt(1d0 +  (val_calc(i)/delta_loss)**2) - 1d0)
        else
            val_calc(i) = 0d0
        end if

    end do

    ! sum up the ssq
    calc_ssq = sum(val_calc)


    end function



    ! Hessian of trend_ssq with respect to x_i
    subroutine hessian_ssq(x,n_i,t_i,d_i,p_i,gam_exp,gauss,nterms,add_vals,datas,datas_id,sizes,mode,input,fit_mode,fit_red,hessian)
    use module_all_types
    implicit none
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "hessian_ssq" :: hessian_ssq
    ! Variables
    integer, dimension(3):: sizes
    integer:: mode,fit_mode ! 1: ssq calc, 2: single prop
    double precision, dimension(sizes(1)) :: x,x_tmp
    double precision, dimension(sizes(1),sizes(1)) :: hessian
    double precision, dimension(2):: add_vals
    integer, dimension(sizes(2)) :: nterms
    double  precision, dimension(nterms(3),4)::gauss
    double precision, dimension(sizes(3),5):: datas
    integer, dimension(sizes(3)):: datas_id
    character(2):: input
    double precision, dimension(sizes(3)):: val_calc
    integer, dimension(sum(nterms)) :: d_i,p_i
    double precision , dimension(nterms(2)) :: gam_exp
    double precision, dimension(sum(nterms))::n_i,t_i !these are necessary if the parameters are not fitted
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    double precision , parameter :: num_step = 1d-4
    double precision , parameter ,dimension(3) :: d_num = (/-num_step,0d0,num_step/)
    double precision , parameter ,dimension(4) :: d_num_x1 = (/num_step,num_step,-num_step,-num_step/)
    double precision , parameter ,dimension(4) :: d_num_x2 = (/num_step,-num_step,num_step,-num_step/)
    double precision , parameter ,dimension(4) :: num_factors_x1_x2 = (/1d0,-1d0,-1d0,1d0/)
    double precision , parameter ,dimension(3) :: num_factors = (/1d0,-2d0,1d0/)
    double precision, dimension(3) :: stencil_ii
    double precision, dimension(4) :: stencil_ij
    type(c_ptr):: handle
    type(type_gl):: gl
    integer:: nterm,i,j,k,iter,nrsubst,fit_red
    logical :: num_virial_mode
    double precision , dimension(size(datas,1)) :: val_calc_dummy


    !init phase
    handle = c_null_ptr
    call  CREATE_FLUID(handle)
    CALL C_F_POINTER(handle,gl)
    nrsubst = 1
    hessian = 0d0


    ! Set the TREND parameters from x
    call set_parameter_from_extern_x(gl,x,n_i,t_i,d_i,p_i,gam_exp,gauss,add_vals,nterms,nrsubst,fit_mode,fit_red)


    ! calculate the critical pressure
    gl%pc(1) = p_calc(gl,gl%tred(1),gl%rhored(1),1)

    nterm = sum(nterms,1)

    !all numbers are real integers
    if(all(int(gl%eos_coeff%di(1:nterm,1)) .eq. gl%eos_coeff%di(1:nterm,1))) then
        num_virial_mode = .false.
    else
        num_virial_mode = .true.
    end if

    ! save x
    x_tmp = x

    !loop over all x_s
    do i=1,size(x)
        do j=i,size(x)
            ! ---------------------------------------------------------------------------
            stencil_ii = 0d0
            stencil_ij = 0d0
            if(i == j) then !diagonal element
                do k = 1,3
                    x(i) = x(i) + d_num(k)
                    call set_helmholtz_param_from_x(gl,x,n_i,d_i,t_i,p_i,gam_exp,gauss,fit_mode,nterm)
                    stencil_ii(k) = num_factors(k) * calc_ssq(gl,x,datas,datas_id,val_calc_dummy,sizes,num_virial_mode,fit_mode,nterm,nrsubst)
                    x = x_tmp
                end do
                hessian(i,j) = sum(stencil_ii)/num_step**2
                ! ---------------------------------------------------------------------------
            else
                ! ---------------------------------------------------------------------------
                do k = 1,4
                    x(i) = x(i) + d_num_x1(k)
                    x(j) = x(j) + d_num_x2(k)
                    call set_helmholtz_param_from_x(gl,x,n_i,d_i,t_i,p_i,gam_exp,gauss,fit_mode,nterm)
                    stencil_ij(k) = num_factors_x1_x2(k) * calc_ssq(gl,x,datas,datas_id,val_calc_dummy,sizes,num_virial_mode,fit_mode,nterm,nrsubst)
                    x = x_tmp
                end do
                hessian(i,j) = sum(stencil_ij)/(4d0*num_step**2)
            end if
            ! ---------------------------------------------------------------------------
            hessian(j,i) = hessian(i,j) ! fill upper triangle matrix
        end do
    end do

    



    end subroutine
    ! ------------------------------------------------------------------------------------------------------------------------------------------

    !Main interface function
    subroutine trend_ssq(x,n_i,t_i,d_i,p_i,gam_exp,gauss,nterms,add_vals,datas,datas_id,sizes,mode,input,fit_mode,fit_red,output,gradient_needed)
    use module_all_types
    implicit none
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "trend_ssq" :: trend_ssq
    ! Variables

    integer, dimension(3):: sizes
    integer:: mode,fit_mode ! 1: ssq calc, 2: single prop
    double precision, dimension(sizes(1)) :: x,gradient
    double precision, dimension(sizes(1) + 1) :: output !first position ssq, then gradient
    double precision, dimension(2):: add_vals
    integer, dimension(sizes(2)) :: nterms
    double  precision, dimension(nterms(3),4)::gauss
    double precision, dimension(sizes(3),5):: datas
    integer, dimension(sizes(3)):: datas_id
    character(2):: input
    double precision, dimension(sizes(3)):: val_calc
    integer, dimension(sum(nterms)) :: d_i,p_i
    double precision , dimension(nterms(2)) :: gam_exp
    double precision, dimension(sum(nterms))::n_i,t_i !these are necessary if the parameters are not fitted
    ! Body of trend_ssq
    integer:: nterm,i,errval, iter, nrsubst,IPHASE,fit_red,gradient_needed
    double  precision:: T,D,rho_vap,rho_liq,press,RHO_EST_GIVEN
    logical:: num_virial_mode ! if = true: numerical approximation of virial coeffs because the d_i are no integers (fit process)
    type(c_ptr):: handle
    type(type_gl):: gl
    double precision:: test1,test2

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !init phase
    handle = c_null_ptr
    call  CREATE_FLUID(handle)
    CALL C_F_POINTER(handle,gl)
    nrsubst = 1
    output = 0d0


    ! Set the TREND parameters from x
    call set_parameter_from_extern_x(gl,x,n_i,t_i,d_i,p_i,gam_exp,gauss,add_vals,nterms,nrsubst,fit_mode,fit_red)


    ! calculate the critical pressure
    gl%pc(1) = p_calc(gl,gl%tred(1),gl%rhored(1),1)

    nterm = sum(nterms,1)

    !all numbers are real integers
    if(all(int(gl%eos_coeff%di(1:nterm,1)) .eq. gl%eos_coeff%di(1:nterm,1))) then
        num_virial_mode = .false.
    else
        num_virial_mode = .true.
    end if

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !sum of squares calculation
    if(mode  ==  1) then
        output(1) = calc_ssq(gl,x,datas,datas_id,val_calc,sizes,num_virial_mode,fit_mode,nterm,nrsubst)

        if(gradient_needed == 1) then
            ! calculate the jacobian matrix val_calc(n) can be used as first stencil (forward difference)
            output(2:) = calc_gradient(gl,datas,datas_id,output(1),sizes,num_virial_mode,nrsubst,x,n_i,d_i,t_i,p_i,gam_exp,gauss,fit_mode,nterm)
        end if
    else

        if( datas_id(1) == 19) then
            output(1) = A00_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif(datas_id(1) == 20) then
            output(1) = A01_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif(datas_id(1) == 21) then
            output(1) = A02_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif(datas_id(1) == 22) then
            output(1) = A03_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif( datas_id(1) == 23) then
            output(1) = A10_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif( datas_id(1) == 24) then
            output(1) = A11_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif( datas_id(1) == 25) then
            output(1) = A11_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif( datas_id(1) == 26) then
            output(1) = A20_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif( datas_id(1) == 27) then
            output(1) = A21_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif( datas_id(1) == 28) then
            output(1) = A30_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif(datas_id(1) == 59) then
            output(1) = DPDD_calc(gl,datas(1,1),datas(1,2),nrsubst)
        elseif(datas_id(1) == 60) then
            output(1) = D2PDD2_CALC(gl,datas(1,1),datas(1,2),nrsubst)
        elseif(datas_id(1) == 35) then
            output(1) = PIP_CALC(gl,datas(1,1),datas(1,2),nrsubst)
        elseif(datas_id(1) == 12) then
            output(1) = B_CALC(gl,datas(1,1),nrsubst)
        elseif(datas_id(1) == 13) then
            output(1) = C_CALC(gl,datas(1,1),nrsubst)
        elseif(datas_id(1) == 41) then
            output(1) = D_CALC(gl,datas(1,1),nrsubst)
        elseif(datas_id(1) == 58) then
            output(1) = DCDT_CALC(gl,datas(1,1),nrsubst)
        elseif(datas_id(1) == 74) then
            output(1) = D2CDT2_CALC(gl,datas(1,1),nrsubst)
        elseif(datas_id(1) == 75) then
            !output(1) = D_GRUEN_DT_CALC(gl,datas(1,1),datas(1,2),nrsubst)
        end if

    endif


    call DESTROY_FLUID(handle)

    end subroutine

    !#####################################################################################################################################################
    !subroutine trend_constraint(x,n_i,t_i,d_i,p_i,gauss,nterms,add_vals,datas,datas_id,sizes,fit_mode,fit_red,output)
    !!(double*,double*,double*,int*,int*,double*,int*,double*,double*,int*,int*,int*,int*);
    !!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "trend_constraint" :: trend_constraint
    !integer, dimension(3):: sizes
    !integer:: fit_red,fit_mode,nrsubst,nterm
    !double precision, dimension(sizes(1)) :: x,gradient
    !double precision, dimension(sizes(1) + 1) :: output !first position ssq, then gradient
    !double precision, dimension(2):: add_vals
    !integer, dimension(sizes(2)) :: nterms
    !double  precision, dimension(nterms(3),4)::gauss
    !double precision, dimension(sizes(3),5):: datas
    !integer, dimension(sizes(3)):: datas_id
    !character(2):: input
    !logical:: num_virial_mode ! if = true: numerical approximation of virial coeffs because the d_i are no integers (fit process)
    !integer, dimension(sum(nterms)) :: d_i,p_i
    !double precision, dimension(sum(nterms))::n_i,t_i !these are necessary if the parameters are not fitted
    !type(c_ptr):: handle
    !type(type_gl):: gl
    !double precision :: d_tmp(sum(nterms))
    !
    !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !!init phase
    !handle = c_null_ptr
    !call  CREATE_FLUID(handle)
    !CALL C_F_POINTER(handle,gl)
    !nrsubst = 1
    !call set_parameter_from_extern(gl,x,n_i,t_i,d_i,p_i,gauss,add_vals,nterms,nrsubst,fit_mode,fit_red)
    !gl%pc(1) = p_calc(gl,gl%tred(1),gl%rhored(1),1)
    !output = 0d0
    !nterm = sum(nterms,1)
    !!all numbers are real integers
    !if(all(int(gl%eos_coeff%di(1:nterm,1)) .eq. gl%eos_coeff%di(1:nterm,1))) then
    !    num_virial_mode = .false.
    !else
    !    num_virial_mode = .true.
    !end if
    !
    !
    !!! if density exponents are fitted and the data id  ==  12,41 Virials
    !!! then the exponents will be rounded
    !!if( (datas_id(1) .eq. 12 .or. datas_id(1) .eq. 41) .and. fit_mode   == 1) then
    !!    d_tmp  = x(2*nterm+1:3*nterm)
    !!    gl%eos_coeff%di(1:nterm,1) =  floor (d_tmp)
    !!end if
    !
    !output(1) = 0d0
    !if( datas_id(1) == 19) then
    !    output(1) = A00_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif(datas_id(1) == 20) then
    !    output(1) = A01_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif(datas_id(1) == 21) then
    !    output(1) = A02_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif(datas_id(1) == 22) then
    !    output(1) = A03_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif( datas_id(1) == 23) then
    !    output(1) = A10_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif( datas_id(1) == 24) then
    !    output(1) = A11_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif( datas_id(1) == 25) then
    !    output(1) = A11_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif( datas_id(1) == 26) then
    !    output(1) = A20_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif( datas_id(1) == 27) then
    !    output(1) = A21_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif( datas_id(1) == 28) then
    !    output(1) = A30_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif(datas_id(1) == 59) then
    !    output(1) = DPDD_calc(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif(datas_id(1) == 60) then
    !    output(1) = D2PDD2_CALC(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif(datas_id(1) == 35) then
    !    output(1) = PIP_CALC(gl,datas(1,1),datas(1,2),nrsubst)
    !elseif(datas_id(1) == 12) then
    !    output(1) = B_CALC(gl,datas(1,1),nrsubst)
    !elseif(datas_id(1) == 13) then
    !    output(1) = C_CALC(gl,datas(1,1),nrsubst)
    !elseif(datas_id(1) == 41) then
    !    output(1) = D_CALC(gl,datas(1,1),nrsubst)
    !elseif(datas_id(1) == 58) then
    !    output(1) = DCDT_CALC(gl,datas(1,1),nrsubst)
    !elseif(datas_id(1) == 74) then
    !    output(1) = D2CDT2_CALC(gl,datas(1,1),nrsubst)
    !elseif(datas_id(1) == 75) then
    !    !output(1) = D_GRUEN_DT_CALC(gl,datas(1,1),datas(1,2),nrsubst)
    !end if
    !
    !
    !!! set d_i back to float numbers
    !!if( (datas_id(1) .eq. 12 .or. datas_id(1) .eq. 41) .and. fit_mode   == 1) then
    !!    gl%eos_coeff%di(1:nterm,1) = d_tmp
    !!end if
    !
    !
    !call DESTROY_FLUID(handle)
    !
    !end subroutine
    !#####################################################################################################################################################



    !#####################################################################################################################################################
    !#####################################################################################################################################################
    double precision function trend_constraint(x,n_i,t_i,d_i,p_i,gam_exp,gauss,nterms,add_vals,datas,datas_id,sizes,fit_mode,fit_red,val_calc)
    use module_all_types
    implicit none
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "TRENDCONSTR" :: trend_constraint
    ! Variables

    integer, dimension(3):: sizes
    integer:: mode,fit_mode ! 1: ssq calc, 2: single prop
    double precision, dimension(sizes(1)) :: x
    double precision, dimension(2):: add_vals
    integer, dimension(sizes(2)) :: nterms
    double  precision, dimension(nterms(3),4)::gauss
    double precision, dimension(sizes(3),5):: datas
    integer, dimension(sizes(3)):: datas_id
    character(2):: input
    double precision, dimension(sizes(3)):: val_calc
    integer, dimension(sum(nterms)) :: d_i,p_i
    double precision, dimension(sum(nterms))::n_i,t_i
    double precision , dimension(nterms(2)) :: gam_exp
    ! Body of trend_ssq
    integer:: nterm,i,errval, iter, nrsubst,IPHASE,fit_red  ,iFlash
    double  precision:: T,D,rho_vap,rho_liq,press,RHO_EST_GIVEN
    logical:: num_virial_mode ! if = true: numerical approximation of virial coeffs because the d_i are no integers (fit process)
    type(c_ptr):: handle
    type(type_gl):: gl
    double precision:: test1,test2

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !init phase
    handle = c_null_ptr
    call  CREATE_FLUID(handle)
    CALL C_F_POINTER(handle,gl)
    nrsubst = 1
    call set_parameter_from_extern_x(gl,x,n_i,t_i,d_i,p_i,gam_exp,gauss,add_vals,nterms,nrsubst,fit_mode,fit_red)
    gl%pc(1) = p_calc(gl,gl%tred(1),gl%rhored(1),1)
    trend_constraint = 0d0
    nterm = sum(nterms,1)
    !all numbers are real integers
    if(all(int(gl%eos_coeff%di(1:nterm,1)) .eq. gl%eos_coeff%di(1:nterm,1))) then
        num_virial_mode = .false.
    else
        num_virial_mode = .true.
    end if

    do i=1,sizes(3)
        if(datas(i,1) == 0d0 .and. (datas_id(i) .ne. (13 .or. 14 .or.41) )) then

            press = 0d0
            rho_vap = 0d0
            rho_liq = 0d0
            iflash = 1
            iter = 0

            call Flash_Pure_PhaseBoundary(gl ,press, datas(i,2), rho_vap, rho_liq, iFlash, errval, iter, nrsubst)
            datas(i,1) =   rho_liq

        end if

        if( datas_id(i) == 19) then
            val_calc(i) = A00_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 23) then
            val_calc(i) = A10_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif(datas_id(i) == 20) then
            val_calc(i) = A01_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 26) then
            val_calc(i) = A20_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 27) then
            val_calc(i) = A21_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 28) then
            val_calc(i) = A30_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 86) then
            val_calc(i) = A40_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 35) then
            val_calc(i) = pip_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 12) then
            val_calc(i) = B_CALC(gl,datas(i,1),nrsubst)
        elseif( datas_id(i) == 13) then
            val_calc(i) = C_CALC(gl,datas(i,1),nrsubst)
        elseif( datas_id(i) == 41) then
            val_calc(i) = D_CALC(gl,datas(i,1),nrsubst)
        elseif(datas_id(i) == 59) then
            val_calc(i) = DPDD_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif(datas_id(i) == 60) then
            val_calc(i) = D2PDD2_CALC(gl,datas(i,1),datas(i,2),nrsubst)
        elseif(datas_id(i) == 34) then
            val_calc(i) = GRUEN_WO_CV0_CALC(gl,datas(i,1),datas(i,2),nrsubst)
        elseif(datas_id(i) == 75) then
            val_calc(i) = GRUEN_WO_CV0_CALC_DT(gl,datas(i,1),datas(i,2),nrsubst)
        elseif(datas_id(i) == 87) then
            val_calc(i) = D_PIP_DT_CALC(gl,datas(i,1),datas(i,2),nrsubst)
        end if
    end do
    !!$omp end parallel do
    trend_constraint = 1234d0

    call DESTROY_FLUID(handle)

    end function trend_constraint

    !#####################################################################################################################################################
    double precision function trendvals(x,n_i,t_i,d_i,p_i,gam_exp,gauss,nterms,add_vals,datas,datas_id,sizes,fit_mode,fit_red,val_calc)
    use module_all_types
    implicit none
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "TRENDVAL" :: trendvals
    ! Variables

    integer, dimension(3):: sizes
    integer:: mode,fit_mode ! 1: ssq calc, 2: single prop
    double precision, dimension(sizes(1)) :: x
    double precision, dimension(2):: add_vals
    integer, dimension(sizes(2)) :: nterms
    double  precision, dimension(nterms(3),4)::gauss
    double precision, dimension(sizes(3),5):: datas
    integer, dimension(sizes(3)):: datas_id
    double precision , dimension(nterms(2)) :: gam_exp
    character(2):: input
    double precision, dimension(sizes(3)):: val_calc
    double precision , dimension(sum(nterms)) :: d_i,p_i
    double precision, dimension(sum(nterms))::n_i,t_i
    ! Body of trend_ssq
    integer:: nterm,i,errval, iter, nrsubst,IPHASE,fit_red
    double  precision:: T,D,rho_vap,rho_liq,press,RHO_EST_GIVEN
    logical:: num_virial_mode ! if = true: numerical approximation of virial coeffs because the d_i are no integers (fit process)
    type(c_ptr):: handle
    type(type_gl):: gl
    double precision:: test1,test2

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !init phase
    handle = c_null_ptr

    call  CREATE_FLUID(handle)

    CALL C_F_POINTER(handle,gl)

    nrsubst = 1

    call set_parameter_from_arrays(gl,n_i,t_i,d_i,p_i,gam_exp,gauss,add_vals,nterms)

    gl%pc(1) = p_calc(gl,gl%tred(1),gl%rhored(1),1)

    trendvals = 0d0

    nterm = sum(nterms,1)

    !all numbers are real integers
    if(all(int(gl%eos_coeff%di(1:nterm,1)) .eq. gl%eos_coeff%di(1:nterm,1))) then
        num_virial_mode = .false.
    else
        num_virial_mode = .true.
    end if

    do i=1,sizes(3)
        if( datas_id(i) == 19) then
            val_calc(i) = A00_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 23) then
            val_calc(i) = A10_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif(datas_id(i) == 20) then
            val_calc(i) = A01_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 26) then
            val_calc(i) = A20_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 27) then
            val_calc(i) = A21_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 28) then
            val_calc(i) = A30_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 86) then
            val_calc(i) = A40_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 35) then
            val_calc(i) = pip_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif( datas_id(i) == 12) then
            val_calc(i) = B_CALC(gl,datas(i,1),nrsubst)
        elseif( datas_id(i) == 13) then
            val_calc(i) = C_CALC(gl,datas(i,1),nrsubst)
        elseif( datas_id(i) == 41) then
            val_calc(i) = D_CALC(gl,datas(i,1),nrsubst)
        elseif(datas_id(i) == 59) then
            val_calc(i) = DPDD_calc(gl,datas(i,1),datas(i,2),nrsubst)
        elseif(datas_id(i) == 60) then
            val_calc(i) = D2PDD2_CALC(gl,datas(i,1),datas(i,2),nrsubst)
        end if
    end do
    !!$omp end parallel do
    trendvals = 1234d0

    call DESTROY_FLUID(handle)

    end function trendvals

    ! end intern trend_Sqq


    ! calculation of real properties .e.g. argon
    double precision function trendreals(sizes,x,n,t,d,p,gauss,datas,datas_id,nterms,add_vals)
    !DEC$ ATTRIBUTES DLLEXPORT, decorate, alias: "trendreals" :: trendreals
    integer, dimension(4) :: sizes
    integer, dimension(3) :: nterms
    double precision, dimension(2) :: add_vals
    double precision, dimension(sizes(1)) :: x
    double precision, dimension(sizes(2)) :: n,t,d,p
    double precision, dimension(sizes(3),4) :: gauss
    double precision, dimension(sizes(4),4) :: datas
    integer, dimension(sizes(4)) :: datas_id
    double precision, allocatable,dimension(:) :: val_calc
    integer :: i
    double precision :: val
    type(c_ptr):: handle
    type(type_gl) :: gl
    double precision :: press, Temp
    double precision :: rhovap_est, rholiq_est
    integer :: iFlash, nrsubst,IPHASE
    integer :: errval, iter,nterm

    handle = c_null_ptr
    call  CREATE_FLUID(handle)
    CALL C_F_POINTER(handle,gl)
    trendreals = 0d0
    ! Body of set_parameter_from_extern_x
    gl%EQ_type = 0
    gl%components(1) = 'argon'
    gl%EQ_type(1) =  1
    gl%mix_Type = 1
    gl%Factor = 1.D3
    gl%factorpress=1.D0
    gl%factortrans=1.D6
    gl%factorrbwr=1.D2
    gl%factor_VS5eta=1.D5
    gl%factor_VS5slj=1.D1
    gl%ncomp = 1
    gl%REQ(1)  = 8.314472d0
    gl%ttp(1) = 0.66d0
    gl%tred(1) = add_vals(1)
    gl%rhored(1) = add_vals(2)
    gl%tc(1) = gl%tred(1)
    gl%rhoc(1) = gl%rhored(1)
    gl%accen =0d0
    gl%rhoredmix = gl%rhored(1)
    gl%tredmix =  gl%tred(1)
    nrsubst = 1
    gl%molfractions = 0d0
    gl%molfractions(1) = 1d0
    gl%rhomaxfluid = 0d0
    gl%rhomaxfluid(1) = 1d6
    allocate(val_calc(sizes(4)))

    if(.not. allocated(gl%eos_coeff)) allocate(gl%eos_coeff)
    gl%eos_coeff%I_pol(1) = nterms(1)
    gl%eos_coeff%I_exp(1) = nterms(2)
    gl%eos_coeff%I_gbs(1) = nterms(3)
    gl%eos_coeff%nreg(1) =  sum(nterms(1:2))
    nterm  = sum(nterms(1:3))

    if(.not.allocated(gl%eos_coeff%ni)) then
        allocate(gl%eos_coeff%ni(nterm,1))
        allocate(gl%eos_coeff%ti(nterm,1))
        allocate(gl%eos_coeff%di(nterm,1))
        allocate(gl%eos_coeff%p_i(nterm,1))
        allocate(gl%eos_coeff%gama(nterm,1))
        allocate(gl%eos_coeff%pli(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%tli(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%eta(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%beta(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%gam(gl%eos_coeff%I_gbs(1),1))
        allocate(gl%eos_coeff%eps(gl%eos_coeff%I_gbs(1),1))
    end if

    ! set the parameters
    gl%eos_coeff%ni(:,1)   = n
    gl%eos_coeff%ti(:,1)   = t
    gl%eos_coeff%di(:,1)   = d
    gl%eos_coeff%p_i(:,1)  = p
    gl%eos_coeff%gama(:,1) = 1d0
    gl%eos_coeff%pli(:,1)  = 2d0
    gl%eos_coeff%tli(:,1)  = 2d0
    gl%eos_coeff%eta(:,1)  = gauss(:,1)
    gl%eos_coeff%beta(:,1) = gauss(:,2)
    gl%eos_coeff%gam(:,1)  =  gauss(:,3)
    gl%eos_coeff%eps(:,1)  =  gauss(:,4)

    gl%pc(1) = p_calc(gl,gl%tred(1),gl%rhored(1),1)
    !calculation of properties
    ! datas(:;1) -> p
    ! datas(:;2) -> t
    ! datas(:,3) -> d
    ! datas(:;1) -> prop

    do i=1,sizes(4)
        val = 0d0
        Temp = 0d0
        press = 0d0
        rhovap_est = 0d0
        rholiq_est = 0d0
        if(datas(i,4) .ne. 0d0) then
            if(datas_id(i) .eq. 111) then ! DL calculation
                iflash = 1
                iter = 0
                call Flash_Pure_PhaseBoundary(gl ,press, datas(i,2), rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst)
                val_calc(i) = rholiq_est/1d3
                !if(isnan(val_calc(i))) val_calc(i) = 0d0
                trendreals = datas(i,4)*(datas(i,3) - val_calc(i))**2 + trendreals
            elseif(datas_id(i) .eq. 112) then ! DV calculation
                iflash = 1
                iter = 0
                call Flash_Pure_PhaseBoundary (gl,press, datas(i,2), rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst)
                val_calc(i) = rhovap_est/1d3
                !if(isnan(val_calc(i))) val_calc(i) = 0d0
                trendreals = datas(i,4)*(datas(i,3) - val_calc(i))**2 + trendreals
            elseif(datas_id(i)  .eq. 113) then !PVT calculation
                IPHASE = 0
                val = rhomix_calc(gl,datas(i,2),datas(i,1), rhovap_est, IPHASE, nrsubst)
                val_calc(i) = val/1d3
                !if(isnan(val_calc(i))) val_calc(i) = 0d0
                trendreals = datas(i,4)*(datas(i,3) - val_calc(i))**2 + trendreals
            elseif(datas_id(i) .eq. 114) then ! PV calculation
                iflash = 1
                iter = 0
                call Flash_Pure_PhaseBoundary (gl, press, datas(i,2), rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst)
                val_calc(i) = press
                !if(isnan(val_calc(i))) val_calc(i) = 0d0
                trendreals = datas(i,4)*(datas(i,1) - val_calc(i))**2 + trendreals
            end if
        end if
        if(isnan(trendreals)) then
            write(*,*)  trendreals
        end if
    end do
    call DESTROY_FLUID(handle)
    end function

    end module interface_fitter