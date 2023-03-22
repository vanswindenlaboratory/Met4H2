    !   AGA Report No. 8 Part 1
    !   Thermodynamic Properties of Natural Gas and Related Gases
    !   DETAIL Equation of State
    !
    !   Calculation of the residual Part with AGA8-DC92

    module mixtures_AGA8_module

    ! necessary modules from modules.f90
    use module_all_types
    use module_parameters
    use variables_transformation_module
    implicit none                   ! all variables need to be declared explicitly

    ! declare global variables
    type aga_8_type

        double precision, dimension(21) :: x_old
        double precision, dimension(:), allocatable :: x
        double precision :: T, D
        double precision :: B, D_red
        double precision, dimension(18) :: B_s
        double precision, dimension(13:58) :: C_sn
        double precision, dimension(58) :: zeta_n, psi_n, beta_n, omega_n

    end type

    type(aga_8_type) :: ag8

    ! iteration constants
    integer, parameter :: i_elements = 21
    ! integer, parameter :: N_coeff = 58


    contains


    subroutine compare_fluids(gl, errorflag)


    type(type_gl) :: gl
    integer :: i, errorflag

    !=========================================================================================================================!
    ! assign mole fractions of matching fluids() to x(i) in following order
    ! Element Compounds i => 1: Methane,     2: Nitrogen,         3: Carbon dioxide, 4: Ethane,            5: Propane,
    !                        6: i-butane,    7: n-Butane,         8: i-pentane,      9: n-Pentane,        10: n-Hexane,
    !                       11: n-Heptane,  12: n-Octane,        13: n-Nonane,      14: n-Decane,         15: Hydrogen,
    !                       16: Oxygen,     17: Carbon monoxide, 18: Water,         19: Hydrogen sulfide, 20: Helium,
    !                       21: Argon
    !=========================================================================================================================!

    do i = 1, gl%ncomp
        if      (gl%components(i) == 'methane') then
            ag8%x_old(1) = gl%molfractions(i)
        else if (gl%components(i) == 'nitrogen') then
            ag8%x_old(2) = gl%molfractions(i)
        else if (gl%components(i) == 'co2') then
            ag8%x_old(3) = gl%molfractions(i)
        else if (gl%components(i) == 'ethane') then
            ag8%x_old(4) = gl%molfractions(i)
        else if (gl%components(i) == 'propane') then
            ag8%x_old(5) = gl%molfractions(i)
        else if (gl%components(i) == 'isobutan') then
            ag8%x_old(6) = gl%molfractions(i)
        else if (gl%components(i) == 'butane') then
            ag8%x_old(7) = gl%molfractions(i)
        else if (gl%components(i) == 'ipentane') then
            ag8%x_old(8) = gl%molfractions(i)
        else if (gl%components(i) == 'pentane') then
            ag8%x_old(9) = gl%molfractions(i)
        else if (gl%components(i) == 'hexane') then
            ag8%x_old(10) = gl%molfractions(i)
        else if (gl%components(i) == 'heptane') then
            ag8%x_old(11) = gl%molfractions(i)
        else if (gl%components(i) == 'octane') then
            ag8%x_old(12) = gl%molfractions(i)
        else if (gl%components(i) == 'nonane') then
            ag8%x_old(13) = gl%molfractions(i)
        else if (gl%components(i) == 'decane') then
            ag8%x_old(14) = gl%molfractions(i)
        else if (gl%components(i) == 'hydrogen') then
            ag8%x_old(15) = gl%molfractions(i)
        else if (gl%components(i) == 'oxygen') then
            ag8%x_old(16) = gl%molfractions(i)
        else if (gl%components(i) == 'co') then
            ag8%x_old(17) = gl%molfractions(i)
        else if (gl%components(i) == 'water') then
            ag8%x_old(18) = gl%molfractions(i)
        else if (gl%components(i) == 'h2s') then
            ag8%x_old(19) = gl%molfractions(i)
        else if (gl%components(i) == 'helium') then
            ag8%x_old(20) = gl%molfractions(i)
        else if (gl%components(i) == 'argon') then
            ag8%x_old(21) = gl%molfractions(i)
        else
            errorflag = -9996
            return
        end if
    end do

    end subroutine  compare_fluids

    !=========================================================================================================================!


    !=========================================================================================================================!

    subroutine mixture_parameters(gl,T,D,B)
    ! necessary parameters to calculate residual part

    use module_parameters

    type(type_gl) :: gl

    integer :: i, j, n, N_components

    !=========================================================================================================================!
    ! adjusted parameters from the module_parameters
    double precision, dimension(:), allocatable :: M_i, E_i, K_i, G_i, Q_i, F_i, S_i, W_i
    double precision, dimension(:,:), allocatable :: E_ij, E_ij_temp2, U_ij, U_ij_temp2, K_ij, K_ij_temp2, G_ij, G_ij_temp2
    !double precision, dimension(:,:), allocatable :: n0_ik, theta0_ik
    !=========================================================================================================================!

    double precision :: T, D
    double precision :: U_5, U_sum1, U_sum2, U
    double precision :: G, G_sum1, G_sum2
    double precision :: Q
    double precision :: F
    double precision :: K_5, K_sum1, K_sum2, K
    double precision :: D_red
    double precision :: B
    double precision, dimension (18) ::  B_s, B_s_sum, B_s_sum1, B_s_sum2, B_snii, B_snij
    double precision, dimension (:,:,:), allocatable :: B_s_mat
    double precision, dimension (13:58) :: C_sn
    double precision, dimension(58) :: zeta_n, psi_n, omega_n, beta_n
    double precision, dimension(21) :: x_D_iteration
    integer, dimension(21) :: components_number
    integer, dimension(:), allocatable :: components_index

    !init vars
    ag8%D = D/1000.d0
    U_sum1= 0.d0
    U_sum2= 0.d0
    G_sum1= 0.d0
    G_sum2= 0.d0
    Q = 0.d0
    F = 0.d0
    K_sum1= 0.d0
    K_sum2= 0.d0
    B_s_sum = 0d0
    B_s_sum1= 0.d0
    B_s_sum2 = 0d0
    B_snii = 0d0
    B_snij =0d0
    B_s    =0d0
    C_sn   =0d0
    zeta_n =0d0
    psi_n  =0d0
    beta_n=0d0
    omega_n=0d0
    U_5= 0d0
    U = 0d0
    G = 0d0
    K_5= 0d0
    K = 0d0
    D_red = 0d0
    B = 0d0

    x_D_iteration = 0
    components_number = 0


    do i = 1, 21
        ! workaround for density iteration if it is not an input property
        if (T == gl%tc(i)) then
            x_D_iteration(i) = 1.d0
        end if
    end do

    if (all(x_D_iteration == 0)) then
        ! reduce the module parameters to the used components
        ag8%x = pack(ag8%x_old, ag8%x_old /= 0)                                 ! new x only consisiting of the components with a mole fraction above 0
        do i = 1, size(ag8%x_old)                                               ! replace moles with the indexes of the components used
            if (ag8%x_old(i) /= 0) then
                components_number(i) = i
            end if
        end do
    else
        ag8%x = pack(x_D_iteration, x_D_iteration /=0)
        do i = 1, size(x_D_iteration)                                           ! replace moles with the indexes of the components used
            if (x_D_iteration(i) /= 0) then
                components_number(i) = i
            end if
        end do
    end if

    N_components = size(ag8%x)                                                  ! the number of components

    if(.not.allocated(components_index)) then
        allocate(components_index(i_elements))
        components_index = 0
    end if
    components_index = pack(components_number, components_number /=0)           ! indexes of components to define columns and rows of the parameters


    if(.not.allocated(M_i)) then
        allocate(M_i(size(components_index)))
        allocate(E_i,K_i,G_i,Q_i,F_i,S_i,W_i,mold=M_i)
    end if

    do j = 1, size(components_index)
        i = components_index(j)
        M_i(j) = M_i_temp(i)
        E_i(j) = E_i_temp(i)
        K_i(j) = K_i_temp(i)
        G_i(j) = G_i_temp(i)
        Q_i(j) = Q_i_temp(i)
        F_i(j) = F_i_temp(i)
        S_i(j) = S_i_temp(i)
        W_i(j) = W_i_temp(i)



        if (.not. allocated(E_ij_temp2)) allocate(E_ij_temp2(21, size(components_index)))
        E_ij_temp2(:,j) = E_ij_temp(:,i)

        if (.not. allocated(U_ij_temp2)) allocate(U_ij_temp2(21, size(components_index)))
        U_ij_temp2(:,j) = U_ij_temp(:,i)

        if (.not. allocated(K_ij_temp2)) allocate(K_ij_temp2(21, size(components_index)))
        K_ij_temp2(:,j) = K_ij_temp(:,i)

        if (.not. allocated(G_ij_temp2)) allocate(G_ij_temp2(21, size(components_index)))
        G_ij_temp2(:,j) = G_ij_temp(:,i)


        !        if (.not. allocated(n0_ik)) allocate(n0_ik(7, size(components_index)))
        !        n0_ik(:,j) = n0_ik_temp(:,i)

        !        if (.not. allocated(theta0_ik)) allocate(theta0_ik(4, size(components_index)))
        !        theta0_ik(:,j) = theta0_ik_temp(:,i)
    end do

    do j = 1, size(components_index)
        i = components_index(j)
        if (.not. allocated(E_ij)) allocate(E_ij(size(components_index), size(components_index)))
        E_ij(j,:) = E_ij_temp2(i,:)

        if (.not. allocated(U_ij)) allocate(U_ij(size(components_index), size(components_index)))
        U_ij(j,:) = U_ij_temp2(i,:)

        if (.not. allocated(K_ij)) allocate(K_ij(size(components_index), size(components_index)))
        K_ij(j,:) = K_ij_temp2(i,:)

        if (.not. allocated(G_ij)) allocate(G_ij(size(components_index), size(components_index)))
        G_ij(j,:) = G_ij_temp2(i,:)
    end do



    U_sum1 = sum(ag8%x(1:N_components)*(E_i(1:N_components))**2.5)

    do i = 1, N_components-1
        do j = i+1, N_components
            U_sum2 = U_sum2+ag8%x(i)*ag8%x(j)*((U_ij(i,j))**5-1)*(E_i(i)*E_i(j))**2.5
        end do
    end do


    U_5 = U_sum1**2+2*U_sum2            ! summands combined to U^5
    U = U_5**0.2                        ! U without exponent



    G_sum1 = sum(ag8%x(1:N_components)*G_i(1:N_components))

    do i = 1, N_components-1
        do j = i+1, N_components
            G_sum2 = G_sum2+ag8%x(i)*ag8%x(j)*(G_ij(i,j)-1)*(G_i(i)+G_i(j))
        end do
    end do

    G = G_sum1+G_sum2                   ! summands combined to G



    Q = sum(ag8%x(1:N_components)*Q_i(1:N_components))



    F = sum(ag8%x(1:N_components)**2*F_i(1:N_components))



    K_sum1 = sum(ag8%x(1:N_components)*K_i(1:N_components)**2.5)

    do i = 1, N_components-1
        do j = i+1, N_components
            K_sum2 = K_sum2+ag8%x(i)*ag8%x(j)*(K_ij(i,j)**5-1)*(K_i(i)*K_i(j))**2.5
        end do
    end do

    K_5 = K_sum1**2+2*K_sum2            ! summands combined to K^5
    K = K_5**0.2                        ! K without exponent


    ! reduced D
    D_red = (K**3)*ag8%D


    ! 2nd virial coefficient B

    !do n = 1, 18
    !
    !    ! pure fluid contribution
    !    do i = 1, N_components
    !        B_snii(n) = 1
    !        ! different cases for B_snii
    !        if  (g_n(n) == 1) then
    !            B_snii(n) = G_i(i)
    !        end if
    !        if (q_n(n) == 1) then
    !            B_snii(n) = Q_i(i)**2
    !        end if
    !        if (f_n(n) == 1) then
    !            B_snii(n) = F_i(i)**2
    !        end if
    !        if (s_n(n) == 1) then
    !            B_snii(n) = S_i(i)**2
    !        end if
    !        if (w_n(n) == 1) then
    !            B_snii(n) = W_i(i)**2
    !        end if
    !
    !        !B_s_sum1(n) = B_s_sum1(n)+ag8%x(i)**2*E_i(i)**(u_n(n))*K_i(i)**3*B_snii(n)
    !        B_s_sum1(n) = B_s_sum1(n)+a_n(n)*E_i(i)**(u_n(n))*K_i(i)**3*B_snii(n)
    !    end do
    !
    !    do i = 1, N_components-1
    !        do j = i+1, N_components
    !            ! different cases for B_snij
    !            B_snij(n) = 1
    !            if  (g_n(n) == 1) then
    !                B_snij(n) = G_ij(i,j)*((G_i(i)+G_i(j)))/2
    !            end if
    !            if (q_n(n) == 1) then
    !                B_snij(n) = Q_i(i)*Q_i(j)
    !            end if
    !            if (f_n(n) == 1) then
    !                B_snij(n) = F_i(i)*F_i(j)
    !            end if
    !            if (s_n(n) == 1) then
    !                B_snij(n) = S_i(i)*S_i(j)
    !            end if
    !            if (w_n(n) == 1) then
    !                B_snij(n) = W_i(i)*W_i(j)
    !            end if
    !
    !            !B_s_sum2(n) = ag8%x(i)*ag8%x(j)*((E_ij(i,j)*(E_i(i)*E_i(j))**0.5)**u_n(n))*(K_i(i)*K_i(j))**1.5*B_snij(n)
    !            B_s_sum2(n) = a_n(n)*((E_ij(i,j)*(E_i(i)*E_i(j))**0.5)**u_n(n))*(K_i(i)*K_i(j))**1.5*B_snij(n)
    !        end do
    !    end do
    !
    !    B_s(n) = B_s_sum1(n)+2*B_s_sum2(n)
    !end do

    if (.not. allocated(B_s_mat)) allocate(B_s_mat(N_components, N_components, 18))
    B_s_mat = 0.d0

    do n = 1, 18
        do i = 1, N_components
            do j = i, N_components
                B_snij(n) = 1d0
                if (g_n(n) == 1) then
                    B_snij(n) = G_ij(i,j)*(G_i(i)+G_i(j))/2
                end if
                if (q_n(n) == 1) then
                    B_snij(n) = Q_i(i)*Q_i(j)
                end if
                if (f_n(n) == 1) then
                    B_snij(n) = F_i(i)*F_i(j)
                end if
                if (s_n(n) == 1) then
                    B_snij(n) = S_i(i)*S_i(j)
                end if
                if (w_n(n) == 1) then
                    B_snij(n) = W_i(i)*W_i(j)
                end if

                B_s_mat(i,j,n) = a_n(n)*((E_ij(i,j)*(E_i(i)*E_i(j))**0.5)**u_n(n))*(K_i(i)*K_i(j))**1.5*B_snij(n)
            end do
        end do
    end do

    ! pure fluid contributions to Bvir
    do i = 1, N_components
        do n = 1, 18
            B_s(n) = B_s(n)+ag8%x(i)**2*B_s_mat(i,i,n)
        end do
    end do

    ! binary pair contributions to Bvir
    do i = 1, N_components - 1
        do j = i + 1, N_components
            do n = 1, 18
                B_s(n) = B_s(n)+2*ag8%x(i)*ag8%x(j)*B_s_mat(i,j,n)
            end do
        end do
    end do

    B = sum(T**(-u_n(1:18))*B_s(1:18))




    do n = 13, 58
        ! different cases for C_sn
        if      (g_n(n) == 1) then
            C_sn(n) = a_n(n)*G*U**u_n(n)
        else if (q_n(n) == 1) then
            C_sn(n) = a_n(n)*(Q**2)*U**u_n(n)
        else if (f_n(n) == 1) then
            C_sn(n) = a_n(n)*F*U**u_n(n)
        else
            C_sn(n) = a_n(n)*U**u_n(n)
        end if

        ! nomenclature to simplify the apperance of the derivative equations
        zeta_n(n)   = C_sn(n)*T**(-u_n(n))
        psi_n(n)    = C_sn(n)*T**(-u_n(n))*(D_red**b_n(n))*exp(-c_n(n)*D_red**k_n(n))
        omega_n(n)  = c_n(n)*k_n(n)*D_red**k_n(n)
        beta_n(n)   = b_n(n)-omega_n(n)
    end do


    ! pass parameter which need to be global
    ag8%B = B
    ag8%D_red = D_red
    ag8%B_s =     B_s
    ag8%C_sn =    C_sn
    ag8%zeta_n =  zeta_n
    ag8%psi_n =   psi_n
    ag8%beta_n =   beta_n
    ag8%omega_n = omega_n

    end subroutine mixture_parameters

    !=========================================================================================================================!

    !=========================================================================================================================!

    subroutine residualpart_AGA8 (T, D, GETDER, FNRDER)

    use module_parameters
    ! brackets [] = fnrder() in TREND
    double precision :: ar_RT,          ar                              ! [1] residual part
    double precision :: dardD_DRT,      DdardD                          ! [2] 1st derivative with respect to density
    double precision :: d2ardD2_D2RT,   D2d2ardD2                       ! [3] 2nd derivative with respect to density
    double precision :: d3ardD3_D3RT,   D3d3ardD3                       ! [8] 3rd derivative with respect to density
    double precision :: dardT_R,        dardT                           ! [4] 1st derivative with respect to temperature
    double precision :: d2ardT2_TR,     d2ardT2_T                       ! [5] 2nd derivative with respect to temperature
    double precision :: d2ardDdT_DR,    d2ardDdT_D                      ! [6] 1st derivative with respect to density and temperature
    double precision :: TdBdT, T2d2BdT2                                 ! 1st and 2nd derivate of 2nd virial coefficient B with respect to T
    double precision :: sum1, sum2, sum3, sum4, sum5                    ! parts of the derivatives
    double precision :: sum6, sum7, sum8, sum9, sum10
    double precision :: dB, d2B                                         ! derivatives of B with respect to T
    double precision :: darRTdT, darRTdT_T                             ! derivatives for the converter
    double precision :: d2arRTdT2, d2arRTdT2_T2
    double precision :: d2arRTdTdD, d2arRTdTdD_T

    integer, dimension(nderivs) :: GETDER                               ! array specifier to indicate which derivative is needed
    double precision, dimension(nderivs) :: AG8DER
    double precision, dimension(nderivs) :: FNRDER                      ! array with the computed values for the derivatives

    type(type_gl) :: gl


    ! parameters
    integer :: n, i, errorflag
    double precision :: T, D
    double precision :: B, D_red
    double precision, dimension(18) :: B_s
    double precision, dimension(13:58) :: C_sn
    double precision, dimension(58) :: zeta_n, psi_n, beta_n, omega_n, dzeta_n, dpsi_n, d2zeta_n, d2psi_n

    FNRDER = 0d0
    AG8DER = 0d0


    B       = ag8%B
    D_red   = ag8%D_red
    B_s     = ag8%B_s
    C_sn    = ag8%C_sn
    zeta_n  = ag8%zeta_n
    psi_n   = ag8%psi_n
    beta_n  = ag8%beta_n
    omega_n = ag8%omega_n



    sum1  = sum(zeta_n(13:18))                                                              !sum1 necessary for ar_RT and dardD_DRT
    sum2  = sum((u_n(13:18)-1)*zeta_n(13:18))												!sum2 necessary for dardT_R and d2ardDdT_DR
    sum3  = sum(u_n(13:18)*(u_n(13:18)-1)*zeta_n(13:18))                                    !sum3 necessary for d2ardT2_TR
    sum4  = sum(psi_n(13:58))                                                               !sum4 necessary for ar_RT
    sum5  = sum(beta_n(13:58)*psi_n(13:58))                                                 !sum5 necessary for dardD_DRT
    sum6  = sum((beta_n(13:58)*(beta_n(13:58)-1)-k_n(13:58)*omega_n(13:58))*psi_n(13:58))   !sum6 necessary for d2ardD2_D2RT
    sum7  = sum(((beta_n(13:58)-2)*(beta_n(13:58)*(beta_n(13:58)-1)-k_n(13:58)*omega_n(13:58))+k_n(13:58)*omega_n(13:58)*(1-k_n(13:58)-2*beta_n(13:58)))*psi_n(13:58))  !sum7 necessary for d3ardD3_D3RT
    sum8  = sum((u_n(13:58)-1)*psi_n(13:58))                                                !sum8 necessary for dardT_R and d2ardDdT_DR
    sum9  = sum(u_n(13:58)*(u_n(13:58)-1)*psi_n(13:58))                                     !sum9 necessary for d2ardT2_TR
    sum10 = sum((u_n(13:58)-1)*(b_n(13:58)-omega_n(13:58))*psi_n(13:58))                    !sum10 necessary for d2ardDdT_DR



    TdBdT    = sum(-u_n(1:18)*(T**(-u_n(1:18)))*B_s(1:18))
    T2d2BdT2 = sum(u_n(1:18)*(u_n(1:18)+1)*(T**(-u_n(1:18)))*B_s(1:18))

    ! residual part and its derivatives
    ar_RT        = B*ag8%D-D_red*sum1+sum4                      ! [1] residual part
    dardD_DRT    = B*ag8%D-D_red*sum1+sum5                      ! [2] 1st derivative with respect to density
    d2ardD2_D2RT = sum6                                         ! [3] 2nd derivative with respect to density
    d3ardD3_D3RT = sum7                                         ! [8] 3rd derivative with respect to density
    !dardT_R      = B*ag8%D+ag8%D*TdBdT+D_red*sum2-sum8                  ! [4] 1st derivative with respect to temperature
    d2ardT2_TR   = 2*ag8%D*TdBdT+ag8%D*T2d2BdT2-D_red*sum3-sum9         ! [5] 2nd derivative with respect to temperature
    !d2ardDdT_DR  = B*ag8%D+ag8%D*TdBdT+D_red*sum2-sum10                 ! [6] 1st derivative with respect to density and temperature



    ! derivatives of B, zeta, psi with respect to T
    dB = 1/T * sum(-u_n(1:18)*T**(-u_n(1:18))*B_s(1:18))
    !dB = sum(-u_n(1:18)*T**(-u_n(1:18)-1)*B_s(1:18))
    d2B = sum(u_n(1:18)*(u_n(1:18)+1)*T**(-u_n(1:18)-2)*B_s(1:18))

    
    
    do n = 13, 18
        dzeta_n(n)  = -u_n(n)*C_sn(n)*T**(-u_n(n)-1)
        d2zeta_n(n) =  u_n(n)*(u_n(n)+1)*C_sn(n)*T**(-u_n(n)-2)
    end do

    do n = 13, 58
        dpsi_n(n)  = -u_n(n)*C_sn(n)*T**(-u_n(n)-1)*(D_red**b_n(n))*exp(-c_n(n)*D_red**k_n(n))
        d2psi_n(n) =  u_n(n)*(u_n(n)+1)*C_sn(n)*T**(-u_n(n)-2)*(D_red**b_n(n))*exp(-c_n(n)*D_red**k_n(n))
    end do

    darRTdT    = ag8%D*dB-D_red*sum(dzeta_n(13:18))+sum(dpsi_n(13:58))
    d2arRTdT2  = ag8%D*d2B-D_red*sum(d2zeta_n(13:18))+sum(d2psi_n(13:58))
    d2arRTdTdD = ag8%D*dB-D_red*sum(dzeta_n(13:18))+sum(beta_n(13:58)*dpsi_n(13:58))
    
    darRTdT_T    = T*darRTdT                                ! [4] 1st derivative with respect to temperature
    d2arRTdT2_T2 = T**2*d2arRTdT2                           ! [5] 2nd derivative with respect to temperature
    d2arRTdTdD_T = T*d2arRTdTdD                             ! [6] 1st derivative with respect to density and temperature

    !    dardT       = dardT_R*T                            ! [4] 1st derivative with respect to temperature
    !    d2ardT2_T   = d2ardT2_TR*Rag8                      ! [5] 2nd derivative with respect to temperature
    !    d2ardDdT_D  = d2ardDdT_DR*Rag8                     ! [6] 1st derivative with respect to density and temperature


    ! fill in AG8DER
    AG8DER(1) = ar_RT
    AG8DER(2) = dardD_DRT
    AG8DER(3) = d2ardD2_D2RT
    AG8DER(4) = darRTdT_T
    AG8DER(5) = d2arRTdT2_T2
    AG8DER(6) = d2arRTdTdD_T
    AG8DER(8) = d3ardD3_D3RT


    ! check if this routine can be used: check functions needed
    ! GETDER has the dimension of 15, but AGA8 only provides 7 functions, so there will be an error if other functions than those 7 are needed for the calculations

    if (GETDER(7) == 1) then
        errorflag = -9995
        return
    end if

    do i = 9, 15
        if (GETDER(i) == 1) then
            errorflag = -9995
            return
        end if
    end do


    ! Transform T to tau
    do i = 4, 6
        if (GETDER(i) == 1) then
            call convert_T_derivs_Trho(gl, GETDER, AG8DER, FNRDER)
        end if
    end do


    ! fill functions in FNRDER
    if (GETDER(1) == 1) then
        FNRDER(1) = ar_RT
    end if

    if (GETDER(2) == 1) then
        FNRDER(2) = dardD_DRT
    end if

    if (GETDER(3) == 1) then
        FNRDER(3) = d2ardD2_D2RT
    end if

    if (GETDER(8) == 1) then
        FNRDER(8) = d3ardD3_D3RT
    end if


    end subroutine residualpart_AGA8


    end module mixtures_AGA8_module