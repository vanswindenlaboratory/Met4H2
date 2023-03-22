subroutine px_moni(gl,Temp, p_points, T_points, x_points, rhovap_points, rholiq_points, points, fileout, switchfluids, errval)
    
    !--------------------------------------------------------------------------------
    ! THIS SUBROUTINE CALCULATES THE PHASE BOUNDARIES OF A BINARY MIXTURE AT GIVEN
    ! TEMPERATURE. THE RESULTS ARE RETURNED IN VECTORS AND CAN BE PRINTED INTO A FILE
    !--------------------------------------------------------------------------------
    ! M. Thol, May 2019

    use module_all_types
    use ancillary_equations_mix_module
    use flash_module
    use flash_pure_module
    implicit none

    type(type_gl) :: gl
    !--------------------------------------------------------------------------------
    integer, parameter:: maxi = 2000 ! maximum number of calculated points, lengt of the return vectors!
    double precision, intent(in):: Temp     ! specified temperature for the p-x diagram
    double precision, dimension(maxi), intent(out):: T_points, p_points, rhovap_points, rholiq_points    ! return vectors of the calculated points
    double precision, intent(out):: x_points(maxi,4)   ! this array is filled according to:
    ! x_points(i,1) = x_vap(1), x_points(i,2) = x_vap(2)
    ! x_points(i,3) = x_liq(1), x_points(i,4) = x_liq(2)
    character(255), intent(in):: fileout    ! optional path for output file
    !switchfluids = 1: fluid order was switched -> this has to be undone by swapping point order and molefractions x1 -> x2, x2 -> x1
    !switchfluids = 2: former azeotropic mixture, splited in two vle regions (co2-ethane T > 290 K)
    !swtichfluids = 3: Mostly relevant for open VLE regions -> start with component with smaller vapor pressure, because this component will most likely be attached to y axis while the other component is detached from y axis
    !switchsupercrit: indicates if the point of minimal x(1) on the dewline for a binary mixtures where one component is supercritical was passed
    integer:: points, switchfluids, errval                ! number of calculated points, error code
    !--------------------------------------------------------------------------------
    
    double precision:: psat_1, psat_2, press
    double precision::  increment, press_old, x_old
    integer:: iter, errpsat, i
    double precision:: T, pstart, d_vap, d_liq, vapfrac
    
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, Nr_x_given, iPhase_try

    
    !initalize variables
    p_points = 0.d0; T_points = 0.d0; x_points = 0.d0; rhovap_points = 0.d0; rholiq_points = 0.d0
    points = 1; errval = 0;
    T = Temp
    iter = 0; errpsat = 0
    d_vap = 0.d0; d_liq = 0.d0
    press_old = 0.d0
    x_old = 0.d0
   
    !p_pts_right = 0.d0; T_pts_right = 0.d0; x_pts_right = 0.d0; rhovap_pts_right = 0.d0; rholiq_pts_right = 0.d0;
    !x_vap = 0.d0; x_liq = 0.d0;  errMSG = ''
    !d_vap = 0.d0; d_liq = 0.d0
    !azeo = .false.
    !
    !deltap_max = 0.5d0  !maximum pressure step
    !psat_high = 0.d0
    !switchsupercrit = 0
    !force_spline = 0
    !save_calc = 0
    !dpdx = 0.d0
    !dpdx_old = 0.d0
    !open_count = 0

    !--------------------------------------------------------------------------------
    ! STEP 1: START VALUES
    !--------------------------------------------------------------------------------
    ! get the vapor pressures of the pure fluids
    psat_1 = 0.d0
    psat_2 = 0.d0

    psat_1 = vp_eq(gl,T, 1)
    if (T < gl%tc(2)) then
        psat_2 = vp_eq(gl,T, 2)
    end if 
    
    ! calculate and save the pure fluid as the first point
    call Flash_Pure_PhaseBoundary(gl,psat_1, t, d_vap, d_liq, 1, errPsat, iter, 1)
    
    T_points(points) = T
    p_points(points) = psat_1
    pstart = psat_1
    x_points(points, 1) = 1.d0 !xvap1
    x_points(points, 2) = 0.d0 !xvap2
    x_points(points, 3) = 1.d0 !xliq1
    x_points(points, 4) = 0.d0 !xliq2
    rhovap_points(points) = d_vap
    rholiq_points(points) = d_liq
    points = points + 1

    open(6932, file = 'D:\Zustandsgleichungen\Fortranfiles SVN\trunk\TREND 4.0\px-diag-out.txt') 
    write(6932, '(2F20.10)') x_points(1, 3), psat_1
    
    
    !   - BUBBLE POINT: T AND x' VEXTOR GIVEN   --  iFlash = 3
    !   - DEW POINT:    T AND x" VECTOR GIVEN   --  iFlash = 4
    
    rhovap_est = d_vap
    rholiq_est = d_liq
    
    x_known = 0.d0
    x_known(1) = 1.d0
    x_known(2) = 1.d0 - x_known(1)
    
    press = 0.d0
    vapfrac = 0.d0
    Nr_x_given = 0
    iter = 0
    x_vap = 0.d0; x_liq = 0.d0
    i = 1
    
    increment = 0.01d0
    
    !calculate bubble line
    iflash = 3
    iPhase_try = 1
    
    do while (x_known(1) .gt. 1.d-10)
        errval = 0
        i = i + 1
        x_known(1) = x_known(1) - increment
        x_known(2) = 1.d0 - x_known(1)
        x_liq = x_known
        call Flash_PhaseBoundary_calc(gl,press, T, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                   & iPhase_try, Nr_x_given, errval, iter)  
        if (errval == 0) then
            p_points(i) = press
            x_points(i,3) = x_liq(1)
            x_points(i,4) = x_liq(2)
        
            write(6932, '(2F20.10)') x_liq(1), press
        end if
    end do
    
    increment = 0.0001d0 
    x_known(1) = 0.005d0

    do while (x_known(1) .gt. 1.d-10)
        errval = 0
        i = i + 1
        x_known(1) = x_known(1) - increment
        x_known(2) = 1.d0 - x_known(1)
        x_liq = x_known
        call Flash_PhaseBoundary_calc(gl,press, T, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                   & iPhase_try, Nr_x_given, errval, iter)  
        if (errval == 0) then
            p_points(i) = press
            x_points(i,3) = x_liq(1)
            x_points(i,4) = x_liq(2)
        
            write(6932, '(2F20.10)') x_liq(1), press
        end if
    end do
    
    
    if (psat_2 .gt. 1.d-12) write(6932, '(2F20.10)') x_points(1, 4), psat_2
    
    !calculate dew line
    i = 1
    rhovap_est = d_vap
    rholiq_est = d_liq
    press = 0.d0
    vapfrac = 0.d0
    Nr_x_given = 0
    iter = 0
    x_vap = 0.d0; x_liq = 0.d0
    x_known(1) = 1.d0
    x_known(2) = 1.d0 - x_known(1)

    iflash = 4
    iPhase_try = 5

    write(6932, '(2F20.10)') x_points(1, 1), psat_1
    
    increment = 0.005d0
    
    do while (x_known(1) .gt. 1.d-3)
        i = i + 1
        x_known(1) = x_known(1) - increment
        x_known(2) = 1.d0 - x_known(1)
        x_vap = x_known
        call Flash_PhaseBoundary_calc(gl,press, T, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                   & iPhase_try, Nr_x_given, errval, iter)   
        if (errval == 0) then
            p_points(i) = press
            x_points(i,1) = x_vap(1)
            x_points(i,2) = x_vap(2)
            write(6932, '(2F20.10)') x_vap(1), press
        end if
    end do

    increment = 0.00001d0 
    
    x_known(1) = 0.003d0
    
    do while (x_known(1) .gt. 1.d-5)
        i = i + 1
        x_known(1) = x_known(1) - increment
        x_known(2) = 1.d0 - x_known(1)
        x_vap = x_known
        call Flash_PhaseBoundary_calc(gl,press, T, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                   & iPhase_try, Nr_x_given, errval, iter)   
        if (errval == 0) then
            p_points(i) = press
            x_points(i,1) = x_vap(1)
            x_points(i,2) = x_vap(2)
        
            write(6932, '(2F20.10)') x_vap(1), press
            press_old = press
            x_old = x_vap(1)
        end if
    end do

   
    x_known(1) = x_old + increment
    press = press_old

    increment = 0.00001d0 
    
    do while (x_known(1) .gt. 1.d-10)
        i = i + 1
        x_known(1) = x_known(1) - increment
        x_known(2) = 1.d0 - x_known(1)
        x_vap = x_known
        call Flash_PhaseBoundary_calc(gl,press, T, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                   & iPhase_try, Nr_x_given, errval, iter)   
        if (errval == 0) then
            p_points(i) = press
            x_points(i,1) = x_vap(1)
            x_points(i,2) = x_vap(2)
        
            write(6932, '(2F20.10)') x_vap(1), press
            press_old = press
            x_old = x_vap(1)
        end if
    end do

    
    x_known(1) = x_old + increment
    press = press_old

    increment = 0.000001d0 
    
    do while (x_known(1) .gt. 1.d-10)
        i = i + 1
        write(*,*) i
        x_known(1) = x_known(1) - increment
        x_known(2) = 1.d0 - x_known(1)
        x_vap = x_known
        call Flash_PhaseBoundary_calc(gl,press, T, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                   & iPhase_try, Nr_x_given, errval, iter)   
        if (errval == 0) then
            p_points(i) = press
            x_points(i,1) = x_vap(1)
            x_points(i,2) = x_vap(2)
        
            write(6932, '(2F20.10)') x_vap(1), press
            press_old = press
            x_old = x_vap(1)
        else
            exit
        end if
    end do

    x_known(1) = x_old + increment
    press = press_old *1.05d0

    increment = 0.000001d0 
    
    do while (x_known(1) .gt. 1.d-10)
        i = i + 1
        write(*,*) i
        x_known(1) = x_known(1) + increment
        x_known(2) = 1.d0 - x_known(1)
        x_vap = x_known
        call Flash_PhaseBoundary_calc(gl,press, T, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                   & iPhase_try, Nr_x_given, errval, iter)   
        if (errval == 0) then
            p_points(i) = press
            x_points(i,1) = x_vap(1)
            x_points(i,2) = x_vap(2)
        
            write(6932, '(2F20.10)') x_vap(1), press
            press_old = press
            x_old = x_vap(1)
        else
            exit
        end if
    end do

    x_known(1) = x_old + increment
    press = press_old *1.05d0

    increment = 0.000001d0 
    
    do while (x_known(1) .gt. 1.d-10)
        i = i + 1
        write(*,*) i
        x_known(1) = x_known(1) - increment
        x_known(2) = 1.d0 - x_known(1)
        x_vap = x_known
        call Flash_PhaseBoundary_calc(gl,press, T, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
                   & iPhase_try, Nr_x_given, errval, iter)   
        if (errval == 0) then
            p_points(i) = press
            x_points(i,1) = x_vap(1)
            x_points(i,2) = x_vap(2)
        
            write(6932, '(2F20.10)') x_vap(1), press
            press_old = press
            x_old = x_vap(1)
        else
            exit
        end if
    end do

    
    continue
    
end subroutine