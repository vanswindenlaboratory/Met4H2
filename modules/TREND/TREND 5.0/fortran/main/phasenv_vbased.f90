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

    ! module for file phasenv_vbased.f90
    module phasenv_vbased_module
    !global use inclusion
    use module_all_types
    use phasenv_pbased_module
    use phasenv_pure_module
    use setup_module
    use cubic_eos_module

    contains


    
    
!**************************************************************************
subroutine PTDIAG_OUT_VBASED(gl,fluids, moles, path, p_spec, T_spec,t_pts_out, p_pts_out, rholiq_pts_out, rhovap_pts_out, &
                        & pointID_pts_out, x_pts_out, fileout, EOS_indicator, errval)
!**************************************************************************
! Interface function for plotting a pt-diagram for mixtures
!DEC$ ATTRIBUTES DLLEXPORT :: PTDIAG_OUT_VBASED





implicit none

    type(type_gl) :: gl


!!! if no pressure / temperature is specified it has to be set to zero !!!
double precision:: p_spec, T_spec                           ! specified pressure / temperature
!integer, intent(inout):: errval                             ! error handling
!character(255), intent(in):: fileout    ! optional path for output file
integer:: errval                             ! error handling
character(255):: fileout    ! optional path for output file
integer, dimension(30):: EOS_indicator
integer :: MIX_indicator

!Output, maximum: 1000 pts
!double precision, dimension(200),intent(out):: t_pts_out, p_pts_out, rholiq_pts_out, rhovap_pts_out
!integer, dimension(200):: pointID_pts_out
!double precision, dimension(200, 30),intent(out):: x_pts_out
double precision, dimension(400),intent(out):: t_pts_out, p_pts_out, rholiq_pts_out, rhovap_pts_out
integer, dimension(400):: pointID_pts_out
double precision, dimension(400, 30),intent(out):: x_pts_out

double precision, dimension(6):: pt_return                  ! return value(s) (temp. or press) at the specified pressure/temperature
double precision, dimension(6,2):: rho_return               ! rho_liq (i,1) and rho_vap (i,2) at the specified pressure/temperature
double precision, dimension(30,6):: x_return                ! return values for the compositions of the emerging phase at specified pressure/temperature
integer:: points_found, trunc_fac, i, k


!-------------------------------------------------------------------
!Initialize
character (12) :: input
double precision :: prop1, prop2
character (30), dimension(30) :: fluids
character (255) :: path
double precision, dimension(30):: moles              ! composition of the system
!double precision, dimension(200):: p_points, T_points, rhovap_points, rholiq_points
!Andreas August 2013
double precision, dimension(400):: p_points, T_points, rhovap_points, rholiq_points
integer:: points
character (255):: filepure
!Dummy values for input, Temp and press
input = 'tp'
prop1 = 300.D0
prop2 = 0.1D0

!call setup
errval = 0
call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator, errval)
!------------------------------------------------------------------
if (errval == 0) then
    t_pts_out = 0.D0
    p_pts_out = 0.D0
    rholiq_pts_out = 0.D0
    rhovap_pts_out = 0.D0
    pointID_pts_out = 0
    x_pts_out = 0.D0
    filepure=''
else
    t_pts_out = errval
    p_pts_out = errval
    rholiq_pts_out = errval
    rhovap_pts_out = errval
    pointID_pts_out = errval
    x_pts_out = errval
    return
end if

if (gl%ncomp == 1) then
    call SATPLOT(gl,p_spec, T_spec, p_points, T_points, rhovap_points, rholiq_points, points, filepure, errval)
    t_pts_out = T_points
    p_pts_out = p_points
    rholiq_pts_out = rholiq_points
    rhovap_pts_out = rhovap_points
    pointID_pts_out = gl%pointID(1:400)
    return
else
    call ptdiag_Vbased(gl,moles, p_spec, T_spec, pt_return, rho_return, x_return, points_found, errval)
end if

!Leave some space in the return array for specified and special points like critical pts, cricondenbar, cricondentherm etc.
if (gl%phasenv_pts > 0) then
    if (gl%phasenv_pts > 385) then
        trunc_fac = ceiling(gl%phasenv_pts / 385.D0)
    !if (phasenv_pts > 185) then
    !   trunc_fac = ceiling(phasenv_pts / 185.D0)
        k = 1
        t_pts_out(k) = gl%T_pts(1)
        p_pts_out(k) = gl%p_pts(1)
        rholiq_pts_out(k) = gl%rholiq_pts(1)
        rhovap_pts_out(k) = gl%rhovap_pts(1)
        pointID_pts_out(k) = gl%pointID(1)
        x_pts_out(k,:) = gl%x_pts(1,:)  
        k = k + 1      
        Do i = 1, gl%phasenv_pts
            if(mod(i,trunc_fac) == 0) then
                t_pts_out(k) = gl%T_pts(i)
                p_pts_out(k) = gl%p_pts(i)
                rholiq_pts_out(k) = gl%rholiq_pts(i)
                rhovap_pts_out(k) = gl%rhovap_pts(i)
                pointID_pts_out(k) = gl%pointID(i)
                x_pts_out(k,:) = gl%x_pts(i,:)  
                k = k + 1       
            End if
            if (k == 400) return
            !if (k == 200) return
            if(gl%pointID(i) /= 0) then
                t_pts_out(k) = gl%T_pts(i)
                p_pts_out(k) = gl%p_pts(i)
                rholiq_pts_out(k) = gl%rholiq_pts(i)
                rhovap_pts_out(k) = gl%rhovap_pts(i)
                pointID_pts_out(k) = gl%pointID(i)
                x_pts_out(k,:) = gl%x_pts(i,:)  
                k = k + 1       
            End if
            if (k == 400) return   
           ! if (k == 200) return 
        End do
    else
        t_pts_out = gl%T_pts(1:400)
        p_pts_out = gl%p_pts(1:400)
        rholiq_pts_out = gl%rholiq_pts(1:400)
        rhovap_pts_out = gl%rhovap_pts(1:400)
        pointID_pts_out = gl%pointID(1:400)
        x_pts_out = gl%x_pts(1:400,:)  
        !t_pts_out = T_pts(1:200)
        !p_pts_out = p_pts(1:200)
        !rholiq_pts_out = rholiq_pts(1:200)
        !rhovap_pts_out = rhovap_pts(1:200)
        !pointID_pts_out = pointID(1:200)
        !x_pts_out = x_pts(1:200,:)  
    End if
else 
    t_pts_out = errval
    p_pts_out = errval
    rholiq_pts_out = errval
    rhovap_pts_out = errval
    pointID_pts_out = errval
    x_pts_out = errval 
End if

! output file format
1000 format(f10.5, f10.4, f10.3, f10.3)

if (fileout /= '') then
    open(unit = 13,file=fileout, status='unknown', action='write', iostat=errval)
    if (errval /= 0) return
    write(13,*)'!---------------------------------------------------------------------------------'
    write(13,*)'! P-T DIAGRAM DATA FILE'
    write(13,*)'!'
    write(13,*)'!---------------------------------------------------------------------------------'
    write(13,*)'!'
    write(13,*)'! p/MPa    ','T/K       ','   dl/mol/l ', 'dv/mol/l'
    do i = 1, 400
        write(13,1000) p_pts_out(i), t_pts_out(i), rholiq_pts_out(i)/1000.d0, rhovap_pts_out(i)/1000.d0
    end do
    close(13)
end if

end subroutine PTDIAG_OUT_Vbased
!**************************************************************************


!**************************************************************************
subroutine ptdiag_Vbased(gl,z, p_spec, T_spec, pt_return, rho_return, x_return, points_found, errval)
!**************************************************************************
! This routine attempts to construct the whole phase envelope for a given mixture.
!   -----------------------------------------------------------------------
! The points on the phase envelope are saved in the module following variables of the module module_VLE, 
! starting with the dew line at low pressure and following the phase envelope past the crit. point until the end of the bubble line is reached.
!    double precision, dimension(imax)::T_pts, p_pts, rhovap_pts, rholiq_pts    --> temperatures, pressures and densities of the exisiting phase (rhovap_pts) and the emerging phase (rholiq_pts)
!    double precision, dimension(imax, 30):: x_pts  --> the composition of the emerging phase
!    integer, dimension(imax)::pointID --> IDs of special points: 0: normal points, 1: crit. point(s), 2: temp. maximum, 3: press. maximum, 4: specified temp, 5: specified press.
!    imax is the maximum number of points that can be stored in the variables and is set in the module module_VLE
! with p_spec and T_spec the user can specify points on the phase envelope, for which the unknown variables are calculated and returned.
!   -----------------------------------------------------------------------
! The construction of the phase envelope starts with the bubble line at low pressures.
! If the calculation of the bubble line succeeds until the vicinity of the crit. point, it continues with the dew line calculation
! if the bubble line calculation fails, the system might have an open phase envelope. The algorithm starts with the dew line 
! and attempts to "jump" over the crit. point to continue with the dew line (or LLE line).
! The crit. point(s) is/are interpolated from adjoining points on the phase envelope.
! They are thus not exact, but usually very accurate estimates, since the interpolation is done with points very close to the crit. point.
! The max. pressure /max. temperature points are constructed by estimating p_max or T_max from cubic interpolations,
! then calculating T(p_max) or p(T_max) also from cubic interpolation, and using those values for a T,x" flash or p,x" flash, respectively.
! Thus, those points are results of an actual EOS flash calculation, but only (again very accurate) approximations of the real maxima.
!   -----------------------------------------------------------------------
! J. Gernert and A. Jäger, August 2012




implicit none

    type(type_gl) :: gl


double precision, dimension(30):: z              ! composition of the system
!!! if no pressure / temperature is specified it has to be set to zero !!!
double precision:: p_spec, T_spec                           ! specified pressure / temperature
double precision, dimension(6):: pt_return                  ! return value(s) (temp. or press) at the specified pressure/temperature
double precision, dimension(6,2):: rho_return               ! rho_liq (i,1) and rho_vap (i,2) at the specified pressure/temperature
double precision, dimension(30,6):: x_return                ! return values for the compositions of the emerging phase at specified pressure/temperature
integer:: points_found                         ! number of points found at specified temperature / pressure
integer:: errval                             ! error handling

integer:: FlashSpec                                         ! indicates if the bubble line (1,3) or the dew line (2, 4) needs to be calculated
double precision, dimension(imax):: t_save, p_save, rholiq_save, rhovap_save 
integer, dimension(imax):: pointID_save
double precision, dimension(imax, 30):: x_save
double precision:: p_crit, T_crit, rho_crit, x_crit(30), para(3), root(3), sum
integer:: bubblepoints, dewpoints, i, k
integer::a(4)
double precision, dimension(4):: p_coeffs, T_coeffs, rholiq_coeffs, rhovap_coeffs
double precision, dimension(30, 4):: x_coeffs

t_save=0.d0; p_save=0.d0; rhovap_save=0.d0; rholiq_save=0.d0; x_save = 0.d0; bubblepoints=0
! first try to calculate the bubble line
Flashspec = -1
call phasenv_Vbased(gl,z, p_spec, T_spec, pt_return, rho_return, x_return, points_found, FlashSpec, errval)
! save the phase envelope points for later use
if (gl%phasenv_pts /= 0) then
    bubblepoints = gl%phasenv_pts
    do i = 1, gl%phasenv_pts
        pointID_save(i) = gl%pointID(gl%phasenv_pts+1-i)
        t_save(i) = gl%T_pts(gl%phasenv_pts+1-i)
        p_save(i) = gl%p_pts(gl%phasenv_pts+1-i)
        rholiq_save(i) = gl%rholiq_pts(gl%phasenv_pts+1-i)
        rhovap_save(i) = gl%rhovap_pts(gl%phasenv_pts+1-i)
        x_save(i,1:gl%ncomp) = gl%x_pts(gl%phasenv_pts+1-i,1:gl%ncomp)
    end do
    FlashSpec = 0
    if (errval == 0) FlashSpec = -2 ! only calculate the dew line until the crit. point
else
    FlashSpec = 0   ! try to jump over the crit. point
end if

errval = 0
gl%T_pts= 0.D0
gl%p_pts = 0.D0
gl%rholiq_pts = 0.D0
gl%rhovap_pts = 0.D0
gl%x_pts = 0.D0
gl%pointID = 0
gl%phasenv_pts = 0
! try to calculate the dew line. Jump over the crit. point, if FlashSpec = 0
call phasenv_Vbased(gl,z, p_spec, T_spec, pt_return, rho_return, x_return, points_found, FlashSpec, errval)
dewpoints = gl%phasenv_pts
if ((dewpoints + bubblepoints) > imax) dewpoints = imax - bubblepoints
! write all the points of the dew and bubble line into the module variables
gl%phasenv_pts = dewpoints+bubblepoints
do i = dewpoints+1, gl%phasenv_pts
    gl%pointID(i) = pointID_save(i-dewpoints)
    gl%T_pts(i) = t_save(i-dewpoints)
    gl%p_pts(i) = p_save(i-dewpoints)
    gl%rholiq_pts(i) = rholiq_save(i-dewpoints)
    gl%rhovap_pts(i) = rhovap_save(i-dewpoints)
    gl%x_pts(i,1:gl%ncomp) = x_save(i-dewpoints,1:gl%ncomp)
end do

! interpolate the critical point in case the phase envelope was constructed with the bubble- and the dew line separately
if ((FlashSpec == -2) .and. (errval == 0)) then
   do i = 6, gl%phasenv_pts
       if ((gl%rhovap_pts(i) - gl%rholiq_pts(i))*(gl%rhovap_pts(i-1) - gl%rholiq_pts(i-1)) < 0.d0) then
            ! interpolate the crit. point using the cubic equation x1_liq(rhovap) from the initial estimates generation
            ! and the relation: x1_liq(rhovap) - z1 = rhovap^3*A + rhovap^2*B + rhovap*C + D - z1 = 0
                !---------------------------------------------------------------
            ! make a polynomial expansion for each variable in the form: 
            ! X(rho_vap) = A*rho_vap^3 + B*rho_vap^2 + C*rho_vap + D
            !---------------------------------------------------------------
            ! This expansion is used to calculate the estimated critical point
            ! set up a linear set of equations from the last four points of each variable    
            a = (/i+5, i, i-1, i-5/)
            call Cubic_Coeffs(gl,34, a, p_coeffs, T_coeffs, x_coeffs, rholiq_coeffs, rhovap_coeffs, errval)

            !-----------------------------------------------------------
            ! generate the three coefficients for the form of the cubic equation: y(x) =x^3 + a*x^2 + b*x + c   
            rho_crit = 0.d0
            para(1) = x_coeffs(1,2)/x_coeffs(1,1)
            para(2) = x_coeffs(1,3)/x_coeffs(1,1)
            para(3) = (x_coeffs(1,4) - z(1))/x_coeffs(1,1)
            x_crit = 0.d0
            ! call the routine that calculates the roots of the cubic equation
            call cubic_nt(gl,para,root)
            ! find the correct root
            do k = 1, 3
                if ((root(k) > 0.d0) .and. ((gl%rhovap_pts(i) - root(k))*(gl%rhovap_pts(i-1) - root(k)) < 0.d0)) rho_crit = root(k)
            end do
            if (rho_crit /= 0.d0) then
                ! calculate the estimated critical parameters from the cubic interpolation and rho_crit
                p_crit = rho_crit**3*p_coeffs(1) + rho_crit**2*p_coeffs(2) + rho_crit*p_coeffs(3) + p_coeffs(4)
                T_crit = rho_crit**3*T_coeffs(1) + rho_crit**2*T_coeffs(2) + rho_crit*T_coeffs(3) + T_coeffs(4)
                sum = 0.d0
                do k = 1, gl%ncomp 
                    x_crit(k) = rho_crit**3*x_coeffs(k,1) + rho_crit**2*x_coeffs(k,2) + rho_crit*x_coeffs(k,3) + x_coeffs(k,4)
                    if (x_crit(k) < 0.d0) x_crit(k) = 1.d-10
                    sum = sum + x_crit(k)
                end do
                x_crit = x_crit/sum
                ! move all points after the critical point one up

                ! overwrite the i point with the crit. point
                gl%x_pts(i, 1:gl%ncomp) = x_crit(1:gl%ncomp)
                gl%T_pts(i) = T_crit
                gl%P_pts(i) = p_crit
                gl%rhovap_pts(i) = rho_crit
                gl%rholiq_pts(i) = rho_crit
                gl%pointID(i) = 1
            end if
            exit
        end if
    end do
end if  

end subroutine ptdiag_Vbased
!**************************************************************************

!**************************************************************************
subroutine phasenv_vbased(gl,z, p_spec, T_spec, pt_return, rho_return, x_return, points_found, FlashSpec, errval)
!**************************************************************************
! This routine calculates the phase envelope of a mixture at given composition.
! The calculated points are saved in the following module variables of the module module_VLE:
!   double precision, dimension(imax):: T_pts, p_pts, rholiq_pts, rhovap_pts
!   double precision, dimension(imax, 31):: x_pts
! For PhaseSpec = 1 or 3:
!   The bubble line is calculated until the specified pressure/temperature is found twice 
!   (retrograde) or the crit. point is reached
! For PhaseSpec = 2 or 4:
!   The dew line is calculated until the specified pressure/temperature is found twice 
!   (retrograde) or the crit. point is reached
! For PhaseSpec = -1 or -2:
!   The bubble line (-1) or dew line (-2) is calculated until the critical point, incl. cricondenbar/cricondentherm 
! For PhaseSpec = 0:
!   Starting on the dew line, either, the maximum number of imax points on the phase boundaries 
!   is calculated or the routine stops when the minimum or maximum pressure on the bubble line is reached
!---------------------------------------------------------------------------
! J. Gernert & A. Jäger, 2010-2012
!---------------------------------------------------------------------------
! Input / Output Variables:
! INPUT:
!   z               feed composition of the mixture
!   p_spec          specified pressure up to which the phase envelope shall be calculated (has to be set to zero if not used!)
!   T_spec          specified temperature up to which the phase envelope shall be calculated (has to be set to zero if not used!)
! OUTPUT
!   pt_return       values for temperature / pressure for the specified value of p_spec / T_spec (vector, length 3)
!   rho_return      values for rho_liq(i,1) and rho_vap(i,2) for the specified value of p_spec / T_spec (matrix, length 3,2)
!   x_return        values for the composition of the emerging phase at the specified value of p_spec / T_spec (matrix, length 30,3)
!   points_found    number of points found at the specified value of p_spec / T_spec (max. value: 3)
!   errval          error value
!---------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl


double precision, dimension(30):: z              ! composition of the system
!!! if no pressure / temperature is specified it has to be set to zero !!!
double precision:: p_spec, T_spec            ! specified pressure / temperature
double precision, dimension(6):: pt_return      ! return value(s) (temp. or press) at the specified pressure/temperature
double precision, dimension(6,2):: rho_return   ! rho_liq (i,1) and rho_vap (i,2) at the specified pressure/temperature
double precision, dimension(30,6):: x_return   ! return values for the compositions of the emerging phase at specified pressure/temperature
integer:: points_found                         ! number of points found at specified temperature / pressure
integer:: errval                             ! error handling
integer:: FlashSpec                             ! indicates if the bubble line (1,3) or the dew line (2, 4) needs to be calculated

double precision, dimension(30)::x_vap, x_liq, x_calc
double precision, dimension(30, 4):: x_coeffs
double precision, dimension(4):: p_coeffs, T_coeffs, rholiq_coeffs, rhovap_coeffs
double precision:: p_start, press, Temp, vapfrac, delta, sum
integer:: i, k, m, iFlash, iter, Nr_x_given, maxVar
logical:: step_ok, step_critical, jump_over_crit, firstcrit_found, cricondentherm, cricondenbar
character(255):: filename
double precision, dimension(30):: dx_drho
double precision:: dT_drho, dp_drho, delta_x_max, C, B, drho_drho
integer, dimension(1):: max_dx_drho
double precision:: step_reduction, step_crit, step_sign, rhocrit_est
double precision:: rhovap_est, rholiq_est, delta_max, rho_vap_new, p_test
logical:: converged
double precision:: p_min, p_max, dpmax, step_factor(5)
integer, dimension(4):: a
double precision:: rho(5), x_Phase(30,5)
integer:: phasetype(5), nrofphases, count, itercheck(5)
double precision::p, q
double precision :: para(3)                ! a,b,c input parameters for cubic equation root solver
double precision :: root(3)              ! output parameters of cubic equation root solver
double precision:: rho_crit, T_crit, p_crit, x_crit(30)

integer:: NumVer        !Definiert das numerische Verfahren

double precision, dimension(4):: sensitiv
integer, dimension(1):: sens_pos

NumVer = 3            !NumVer = 1: Newton
                        !NumVer = 2: Levenberg-Marquardt
                        !NumVer = 3: Dog Leg

!initialization of variables
x_vap = 0.D0; x_liq = 0.d0; x_calc = 0.d0
Temp = 0.D0; vapfrac = 0.D0
p_coeffs = 0.D0; T_coeffs = 0.D0; x_coeffs = 0.D0

Nr_x_given = 1
step_ok = .false.
step_critical = .false.
firstcrit_found = .false.
cricondentherm = .false.
cricondenbar = .false.
rhovap_est = 0.D0; rholiq_est = 0.D0
jump_over_crit = .false.
p_min = 0.1D0   !Lower pressure limit on the bubble line
p_max = 300.d0 ! upper pressure limit on the bubble line (open phase envelope)
dpmax = 2.d0    ! maximum pressure step
points_found = 0
errval = 0

pt_return = 0.d0; rho_return = 0.d0; x_return = 0.d0
rho_crit = 0.d0; T_crit = 0.d0; p_crit = 0.d0; root = 0.d0; x_crit = 0.d0
itercheck = (/2, 4, 7, 15, 20/) ! step checkpoints for the number of iterations from the previous calculation
step_factor = (/0.5d0, 2.d0, 2.d0, 2.d0, 2.d0/) ! factors for the step reduction, depending on the number of iterations used for the previous step
dx_drho = 0.d0 
 
! if no pressure is specified set p_spec to a very large value to evade wrong if-queries
If (p_spec == 0.d0) p_spec = 1.d6   
If (T_spec == 0.d0) T_spec = 1.d6 
!---------------------------------------------------------------
! Step 1: define a start pressure for the phase envelope
!---------------------------------------------------------------
! normal start pressure is 1 bar
p_start = 0.1D0

! check if start pressure is lower than a ficticious mixture triple point pressure
p_test = 0.d0
do k = 1, gl%ncomp
    p_test = p_test + gl%ptp(k)*z(k)
end do
if (p_start < p_test) p_start = p_test

! check if start pressure is larger than the specified pressure
if (p_start > p_spec) p_start = p_spec*0.9

press = p_start

iFlash = FlashSpec
if ((iFlash == 4) .OR. (iFlash == 0) .OR. (iFlash == -2)) iFlash = 2
if ((iFlash == 3) .OR. (iFlash == -1)) iFlash = 1
if ((iFlash < 1) .OR. (iFlash > 2)) then
    errval = -5500
    return
end if

!---------------------------------------------------------------
! Step 2: create four points at low pressures
!---------------------------------------------------------------
! start with pressure-controlled steps on the dew line
i = 1
! try to generate a starting point

call Flash_PhaseBoundary_calc(gl,press, Temp, z, x_vap, x_liq, 0.d0, 0.d0, vapfrac, iFlash, 5, &
    & Nr_x_given, errval,iter)  !changed the IPhase_try from 0 to 5 (assumption: VLE) - J.G. 11/2012
if ((((minval(x_liq(1:gl%ncomp)) < 1.d-12) .OR. (minval(x_vap(1:gl%ncomp)) < 1.d-12)) .and. (errval /= 0)) &
            & .OR. (dabs(gl%rho_vap - gl%rho_liq) < 1.d-5)) then
    errval = 0; Temp = 0.d0; x_vap = 0.d0; x_liq = 0.d0
    call Succ_Sub(gl,press, Temp, z, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, errval, 100, converged)
    if (.not. converged .OR. ((gl%rho_vap - gl%rho_liq)/gl%rho_liq < 1.d-2)) errval = -2222
end if
if (iFlash == 1) x_calc = x_vap
if (iFlash == 2) x_calc = x_liq
! check if the flash calculation succeeded
if ((errval /= 0) .OR. (dabs(gl%rho_vap - gl%rho_liq) < 1.d-5)) then !generation of starting point failed. Try with Stability analysis, using a temperature from Rachford Rice Eqn. as input
    Temp = 0.d0
    x_calc = 0.d0; x_vap = 0.d0; x_liq = 0.d0
    nrofphases = 1
    count = 0
    call PTX_startvals_PhaseBoundary (gl,press, Temp, z, x_vap, x_liq, vapfrac, iFlash, errval)
    if (errval /= 0) then
        errval = -5501
        return
    end if
    do while (nrofphases < 2)
        count = count + 1
        call PhaseDet2(gl,press, Temp, z, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
        if (errval /= 0) nrofphases = 1
        if ((nrofphases < 2) .and. (abs(iFlash) == 2)) Temp = Temp*0.98d0
        if ((nrofphases < 2) .and. (abs(iFlash) == 1)) Temp = Temp/0.98d0
        if (count > 30) then ! Generation of start value failed!   !changed criterion from count > 3 to count > 30 since for very hydrogen rich mixtures ptx_startvals generates a nonsense temperature estimate. 
            errval = -5502
            return
        end if
    end do

    ! Use the calculated phase composition and the Temperature as initial values for the p,x flash
    if (iFlash == 1) x_vap = x_Phase(:,phasetype(iFlash))
    if (iFlash == 2) x_liq = x_Phase(:,phasetype(iFlash))
    rhovap_est = gl%rho_vap
    rholiq_est = gl%rho_liq
    call Flash_PhaseBoundary_calc(gl,press, Temp, z, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, 5, &
        & Nr_x_given, errval,iter)  !changed the IPhase_try from 0 to 5 (assumption: VLE) - J.G. 11/2012
    if ((errval /= 0) .OR. (dabs(gl%rho_vap - gl%rho_liq) < 1.d-5)) then ! Generation of start value failed!
        if (iFlash == 1) x_vap = x_Phase(:,phasetype(iFlash))
        if (iFlash == 2) x_liq = x_Phase(:,phasetype(iFlash))
        call Succ_Sub(gl,press, Temp, z, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, errval, 100, converged)
        if (.not. converged .OR. ((gl%rho_vap - gl%rho_liq)/gl%rho_liq < 1.d-2)) then
            errval =  -5503
            return
        else
            errval = 0
        end if
    end if
end if

!write the results in the module variables
gl%T_pts(i) = Temp
gl%P_pts(i) = press
gl%phasenv_pts = i
if (iFlash == 1) then
    gl%x_pts(i, 1:gl%ncomp) = x_vap(1:gl%ncomp)
    gl%rholiq_pts(i) = gl%rho_vap
    gl%rhovap_pts(i) = gl%rho_liq
else
    gl%x_pts(i, 1:gl%ncomp) = x_liq(1:gl%ncomp)
    gl%rholiq_pts(i) = gl%rho_liq
    gl%rhovap_pts(i) = gl%rho_vap
end if

press = press*1.1D0 !increase the pressure for the next step
errval = 0

! generate three more points at low pressures, using the previuos results as initial guesses
do i = 2, 4
    call Flash_PhaseBoundary_calc(gl,press, Temp, z, x_vap, x_liq, gl%rho_vap, gl%rho_liq, vapfrac, iFlash, 5, &
        & Nr_x_given, errval, iter)  !changed the IPhase_try from 0 to 5 (assumption: VLE) - J.G. 11/2012
    if (((minval(x_liq(1:gl%ncomp)) < 1.d-12) .and. (errval /= 0)) .OR. (dabs(gl%rho_vap - gl%rho_liq) < 1.d-5)) then
        !Andreas August 2012
        !errval = 0; Temp = 0.d0; x_vap = 0.d0; x_liq = 0.d0
        !Do not generate new starting values but take the ones from the last iteration
        errval = 0
        Temp = gl%T_pts(i-1)
        x_liq(:) = gl%x_pts(i-1,:)
        call Succ_Sub(gl,press, Temp, z, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, errval, 100, converged)
        if (.NOT. converged) errval = -2222
        iter = 3
    end if
    if (errval /= 0.D0) then
        errval = -5504
        return
    end if
    !write the results in the module variables
    gl%T_pts(i) = Temp
    gl%P_pts(i) = press
    gl%phasenv_pts = i
    if (iFlash == 1) then
        gl%x_pts(i, 1:gl%ncomp) = x_vap(1:gl%ncomp)
        gl%rholiq_pts(i) = gl%rho_vap
        gl%rhovap_pts(i) = gl%rho_liq
    else
        gl%x_pts(i, 1:gl%ncomp) = x_liq(1:gl%ncomp)
        gl%rholiq_pts(i) = gl%rho_liq
        gl%rhovap_pts(i) = gl%rho_vap
    end if
    
    ! check if the last step went over the specified pressure
    if ((gl%p_pts(i) - p_spec)*(gl%p_pts(i-1) - p_spec) < 0.d0) then
        Temp = (gl%T_pts(i) + gl%T_pts(i-1))/2.d0
        x_calc(1:gl%ncomp) = (gl%x_pts(i, 1:gl%ncomp) + gl%x_pts(i-1, 1:gl%ncomp))/2.d0
        if (iFlash == 1) x_vap = x_calc
        if (iFlash == 2) x_liq = x_calc
        call Flash_PhaseBoundary_calc(gl,p_spec, Temp, z, x_vap, x_liq, gl%rho_vap, gl%rho_liq, vapfrac, iFlash, 5, &
            & Nr_x_given, errval, iter)  !changed the IPhase_try from 0 to 5 (assumption: VLE) - J.G. 11/2012
        if (errval /= 0) then
            errval = -5505
            return
        end if
        points_found = points_found + 1
        if (iFlash == 1) then
            gl%x_pts(i, 1:gl%ncomp) = x_vap(1:gl%ncomp)
            gl%rholiq_pts(i) = gl%rho_vap
            gl%rhovap_pts(i) = gl%rho_liq
            rho_return(points_found,1) = gl%rho_vap
            rho_return(points_found,2) = gl%rho_liq
            x_return(1:gl%ncomp, points_found) = x_vap(1:gl%ncomp)
        else
            gl%x_pts(i, 1:gl%ncomp) = x_liq(1:gl%ncomp)
            gl%rholiq_pts(i) = gl%rho_liq
            gl%rhovap_pts(i) = gl%rho_vap
            rho_return(points_found,1) = gl%rho_liq
            rho_return(points_found,2) = gl%rho_vap
            x_return(1:gl%ncomp, points_found) = x_liq(1:gl%ncomp)
        end if
        
        pt_return(points_found) = Temp
        !write the results in the module variables
        gl%T_pts(i) = Temp
        gl%P_pts(i) = press
        gl%phasenv_pts = i
        gl%pointID(i) = 5
    end if
    ! check if the last step went over the specified temperature
    if ((gl%T_pts(i) - T_spec)*(gl%T_pts(i-1) - T_spec) < 0.d0) then
        press = (gl%p_pts(i) + gl%p_pts(i-1))/2.d0
        x_calc(1:gl%ncomp) = (gl%x_pts(i, 1:gl%ncomp) + gl%x_pts(i-1, 1:gl%ncomp))/2.d0
        if (iFlash == 1) x_vap = x_calc
        if (iFlash == 2) x_liq = x_calc
        call Flash_PhaseBoundary_calc(gl,press, T_spec, z, x_vap, x_liq, gl%rho_vap, gl%rho_liq, vapfrac, iFlash+2, 5,&
                & Nr_x_given, errval, iter)  !changed the IPhase_try from 0 to 5 (assumption: VLE) - J.G. 11/2012
        if (errval /= 0) then
            errval = -5506
            return
        end if
        points_found = points_found + 1
        if (iFlash == 1) then
            gl%x_pts(i, 1:gl%ncomp) = x_vap(1:gl%ncomp)
            gl%rholiq_pts(i) = gl%rho_vap
            gl%rhovap_pts(i) = gl%rho_liq
            rho_return(points_found,1) = gl%rho_vap
            rho_return(points_found,2) = gl%rho_liq
            x_return(1:gl%ncomp, points_found) = x_vap(1:gl%ncomp)
        else
            gl%x_pts(i, 1:gl%ncomp) = x_liq(1:gl%ncomp)
            gl%rholiq_pts(i) = gl%rho_liq
            gl%rhovap_pts(i) = gl%rho_vap
            rho_return(points_found,1) = gl%rho_liq
            rho_return(points_found,2) = gl%rho_vap
            x_return(1:gl%ncomp, points_found) = x_liq(1:gl%ncomp)
        end if
        pt_return(points_found) = press
        !write the results in the module variables
        gl%T_pts(i) = Temp
        gl%P_pts(i) = press
        gl%phasenv_pts = i
        gl%pointID(i) = 4
    end if

    press = press*1.1D0 !increase the pressure for the next step
end do

if (iFlash == 1) then
    x_calc(1:gl%ncomp) = x_vap(1:gl%ncomp)
    gl%rho_vap = gl%rhovap_pts(i-1)
    gl%rho_liq = gl%rholiq_pts(i-1)
else
    x_calc(1:gl%ncomp) = x_liq(1:gl%ncomp)
end if
x_vap = z

do while (i < imax)        
    !---------------------------------------------------------------
    ! Step 3: make a polynomial expansion for each variable in the form: 
    ! X(rho_vap) = A*rho_vap^3 + B*rho_vap^2 + C*rho_vap + D
    !---------------------------------------------------------------
    ! This expansion is used to calculate start values for the next point on the phase envelope
    ! set up a linear set of equations from the last four points of each variable    
    a = (/i-1, i-2, i-3, i-4/)
    call Cubic_Coeffs(gl,34, a, p_coeffs, T_coeffs, x_coeffs, rholiq_coeffs, rhovap_coeffs, errval)
    
    !Derivatives of all molefractions with respect to rho_vap
    do k = 1, gl%ncomp 
        dx_drho(k) = 3.d0 * gl%rho_vap**2 * x_coeffs(k,1) + 2.d0 * gl%rho_vap * x_coeffs(k,2) + x_coeffs(k,3)
        dx_drho(k) = dabs(dx_drho(k)) !/ x_pts(i-1,k))
    end do 
    
    !Derivatives of the temperature with respect to rho_vap
    dT_drho = 3.d0 * gl%rho_vap**2 * T_coeffs(1) + 2.d0 * gl%rho_vap * T_coeffs(2) + T_coeffs(3) 
    dT_drho = dabs(dT_drho / Temp)

    !Derivatives of the pressure with respect to rho_vap
    dp_drho = 3.d0 * gl%rho_vap**2 * p_coeffs(1) + 2.d0 * gl%rho_vap * p_coeffs(2) + p_coeffs(3) 
    dp_drho = dabs(dp_drho / press)
    
    !Ableitung der Dichte nach der Dichte immer 1
    drho_drho = 1
    drho_drho = 1 / gl%rho_vap    
        
    !Get the position of the maximum sensitivity of the compositions w.r.t. density
    max_dx_drho = maxloc(abs(dx_drho)) 
    
    !If the last component is the most sensitive, the component with the next largest sensivity is chosen
    if (max_dx_drho(1) == gl%ncomp) then
        dx_drho(max_dx_drho(1)) = 0.D0
        max_dx_drho = maxloc(abs(dx_drho))
    end if
    
    !get the maximum sensitivity from max(dx_drho), dp_drho, dT_drho and drho_drho
        sensitiv(1) = dp_drho
        sensitiv(2) = dT_drho
        sensitiv(3) = dx_drho(max_dx_drho(1))
        sensitiv(4) = drho_drho
    
        !get the position of the highest value of the vector sensitiv
        sens_pos = 0
        sens_pos = maxloc(abs(sensitiv))
        
        !set iFlash depending on the most sensitive variable
        maxVar = 0
        if (sens_pos(1) == 1) then
            maxVar = 31 !the pressure is the most sensitve variable!
            iFlash = 2
        else if (sens_pos(1) == 2) then
            maxVar = 32 !the temperature is the most sensitve variable!
            iFlash = 4
        else if (sens_pos(1) == 3) then
            maxVar = max_dx_drho(1) !one of the compositions is the most sensitvie variable!
            iFlash = 6
        else if (sens_pos(1) == 4) then
            maxVar = 33 !the density is the most sensitive variable!
            iFlash = 8
        end if
        
                !maxVar = 0
                !if (dx_drho(max_dx_drho(1)) > dp_drho) then
                !    if (dx_drho(max_dx_drho(1)) > dT_drho) then
                !        if (dx_drho(max_dx_drho(1)) > ) then
                !        maxVar = max_dx_drho(1) ! one of the compositions is the most sensitvie variable!
                !        iFlash = 6 !set iFlash depending on the most sensitive variable
                !    else
                !        maxVar = 32 ! the temperature is the most sensitve variable!
                !        iFlash = 4 !set iFlash depending on the most sensitive variable
                !    end if
                !else
                !    if (dp_drho > dT_drho) then 
                !        maxVar = 31 ! the pressure is the most sensitve variable!
                !        iFlash = 2 !set iFlash depending on the most sensitive variable
                !    else
                !        maxVar = 32 ! the temperature is the most sensitve variable!
                !        iFlash = 4 !set iFlash depending on the most sensitive variable
                !    end if
                !end if
    
    !---------------------------------------------------------------
    ! Step 4: select the step size for the next pressure / fraction / density step
    !---------------------------------------------------------------
    !The step depends on the iterations necessary for calculating the last point 
    step_reduction = 4.D0
    if (iter < itercheck(1)) step_reduction = step_reduction*step_factor(1)
    if (iter > itercheck(2)) step_reduction = step_reduction*step_factor(2) 
    if (iter > itercheck(3)) step_reduction = step_reduction*step_factor(3)
    if (iter > itercheck(4)) step_reduction = step_reduction*step_factor(4)
    if (iter > itercheck(5)) step_reduction = step_reduction*step_factor(5)

    ! define the maximum step length for the next density step
    if (maxVar < 31) then
        delta_max = 0.05d0 ! maximum composition step
        delta_max = (delta_max/dx_drho(maxVar)) ! calculate a maximum density step from delta_rho_max = delta_x_max/(dx_max_drho)
    else if (maxVar == 31) then 
        delta_max = .5d0    ! maximum pressure step (in MPa)
        delta_max = (delta_max/dp_drho) ! calculate a maximum density step from delta_rho_max = delta_p/(dp_drho)
    else if (maxVar == 32) then
        delta_max = 5.d0    ! maximum temperature step (in K)
        delta_max = (delta_max/dT_drho) ! calculate a maximum density step from delta_rho_max = delta_T/(dT_drho)
    else if (maxVar == 33) then
        delta_max = 50.D0
    end if
    
   do while (.Not.step_ok) 
        ! calculate the standard density step size
        delta = abs(((gl%rho_liq + gl%rho_vap)/2.d0 - gl%rho_vap)/step_reduction)
        !check if standard step size or maximum step size are used
        if (delta > delta_max) delta = delta_max
        !check if step size is too large at small density values
        if (delta > 0.1d0*gl%rho_vap) delta = 0.1d0*gl%rho_vap
        !set the new variable
        if ((FlashSpec == 1) .or. (FlashSpec == 3).or. (FlashSpec == -1)) delta = -delta
        rho_vap_new = gl%rho_vap + delta  
              
        ! check if the density is very close to the critical density
        if (.NOT. firstcrit_found) then 
            ! estimate for the rel. distance to the crit. point of a vapor-liquid equilibrium (1st crit. point):
            ! use the mean of rho_liq and rho_vap as estimate for the crit. density
            step_crit = abs(0.5d0*(gl%rho_liq  - gl%rho_vap)/gl%rho_vap)
        else
            if (rho_crit /= 0.d0) then  ! if the crit. point was interpolated before, use this value
                rhocrit_est = rho_crit
            else
                ! estimate for the rel. distance to the crit. point of a liquid-liquid equilibrium (> 1st crit. point):
                ! use intersection of linear extrapolations of rho_liq1(p) and rho_liq2(p) as estimate for the crit. density
                C = (gl%rholiq_pts(i-1) - gl%rholiq_pts(i-3))/(gl%p_pts(i-1) - gl%p_pts(i-3))
                B = (gl%rhovap_pts(i-1) - gl%rhovap_pts(i-3))/(gl%p_pts(i-1) - gl%p_pts(i-3))
                rhocrit_est = B*((gl%rhovap_pts(i-1) - gl%rholiq_pts(i-1) - (B - C)*gl%p_pts(i-1))/(C - B) - gl%p_pts(i-1)) + gl%rhovap_pts(i-1)
            end if
            step_crit = dabs((rhocrit_est - gl%rho_vap)/gl%rho_vap)
        end if
        
        !check, if the jump over the crit. point was successful
        if (step_critical) then
            ! if the expression (rho1/rho2 - 1) changed signs in the last two steps
            step_sign = sign(1.d0,((gl%rhovap_pts(i-2)/gl%rholiq_pts(i-2) - 1.D0)*(gl%rhovap_pts(i-1)/gl%rholiq_pts(i-1) - 1.D0)))
            if (step_sign < 0.D0) then 
                jump_over_crit = .true.
                firstcrit_found = .true. ! memorize that at least one crit. point was found
            end if
        end if
        
        if (step_crit < 0.05) then ! first jump criterion:  density is very close to crit. density
            delta_x_max = 0.D0
            do k = 1, gl%ncomp ! get the maximum relative distance between the composition and the crit. composition
                if (dabs((x_calc(k) - z(k))/z(k)) > delta_x_max) delta_x_max = dabs((x_calc(k) - z(k))/z(k))
            end do
            if (delta_x_max < 0.05d0) then ! second jump criterion: composition is very close to crit. composition 
                !----------------------------------------------------------------------------
                ! if only the bubble line needs to be calculated, the algorithm exits here!
                !----------------------------------------------------------------------------
                if (FlashSpec /= 0) return
                !----------------------------------------------------------------------------
                !check if the jump over the crit. point was done in one of the previous steps AND the compositions approach the crit. comp.
                if (.NOT. jump_over_crit .and. (dabs(gl%x_pts(i-2, 1) - z(1)) > dabs(gl%x_pts(i-1, 1)- z(1)))) then 
                    x_calc(max_dx_drho(1)) = 2.d0 * z(max_dx_drho(1)) - x_calc(max_dx_drho(1)) 
                    step_ok = .true.
                    step_critical = .true.
                    iFlash = 6
                else
                    step_critical = .false.     !If successfully jumped once, carry on with density steps
                end if          
            end if
        else ! once outside the crit. region, reset the jump criteria
            step_critical = .false.
            jump_over_crit = .false.
        end if

    !---------------------------------------------------------------
    ! Step 5: calculate new starting values for all unknown variables
    !---------------------------------------------------------------
        sum = 0.D0
        if (.NOT. step_critical) then
            do k = 1, gl%ncomp 
                x_calc(k) = rho_vap_new**3*x_coeffs(k,1) + rho_vap_new**2*x_coeffs(k,2) &
                        & + rho_vap_new*x_coeffs(k,3) + x_coeffs(k,4)
                sum = sum + x_calc(k)
            end do
            press = rho_vap_new**3*p_coeffs(1) + rho_vap_new**2*p_coeffs(2) + rho_vap_new*p_coeffs(3) + p_coeffs(4)
            Temp = rho_vap_new**3*T_coeffs(1) + rho_vap_new**2*T_coeffs(2) + rho_vap_new*T_coeffs(3) + T_coeffs(4)
            rhovap_est = rho_vap_new**3*rhovap_coeffs(1) + rho_vap_new**2*rhovap_coeffs(2) + rho_vap_new*rhovap_coeffs(3)&
                        &  + rhovap_coeffs(4)
            rholiq_est = rho_vap_new**3*rholiq_coeffs(1) + rho_vap_new**2*rholiq_coeffs(2) + rho_vap_new*rholiq_coeffs(3)&
                        &  + rholiq_coeffs(4)
            ! renorming of the composition vector
            if (sum /= 1.D0) x_calc(1:gl%ncomp) = x_calc(1:gl%ncomp) / sum 
        else
            ! calculate the new most sensitive x_k from x_pts(i-1, k) + dx/drho(k)*delta_rho - J. Gernert, Mai 2012
            !x_calc(max_dx_drho(1)) = x_calc(max_dx_drho(1)) + dx_drho(max_dx_drho(1))*delta
            sum = sum + x_calc(max_dx_drho(1))
            !Linear extrapolation of all variables in dependence of the most sensitive x 
            do k = 1, gl%ncomp 
                if (k /= max_dx_drho(1)) then
                    x_calc(k) = (gl%x_pts(i-2, k) - gl%x_pts(i-1, k))/ (gl%x_pts(i-2,  max_dx_drho(1)) - gl%x_pts(i-1,  max_dx_drho(1)))* &
                                & ( x_calc(max_dx_drho(1))  - gl%x_pts(i-1,  max_dx_drho(1))) + gl%x_pts(i-1, k)
                    sum = sum + x_calc(k)
                end if
            end do   
            Temp = (gl%T_pts(i-2) - gl%T_pts(i-1))/ (gl%x_pts(i-2,  max_dx_drho(1)) - gl%x_pts(i-1,  max_dx_drho(1)))* &
                            & ( x_calc(max_dx_drho(1))  - gl%x_pts(i-1,  max_dx_drho(1))) + gl%T_pts(i-1)
            press = (gl%p_pts(i-2) - gl%p_pts(i-1))/(gl%x_pts(i-2,  max_dx_drho(1)) - gl%x_pts(i-1,  max_dx_drho(1)))* &
                            & ( x_calc(max_dx_drho(1))  - gl%x_pts(i-1,  max_dx_drho(1))) + gl%p_pts(i-1)
            rholiq_est = (gl%rholiq_pts(i-2) - gl%rholiq_pts(i-1))/&
                            & (gl%x_pts(i-2,  max_dx_drho(1)) - gl%x_pts(i-1,  max_dx_drho(1)))* &
                            & ( x_calc(max_dx_drho(1))  - gl%x_pts(i-1,  max_dx_drho(1))) + gl%rholiq_pts(i-1) 
            rhovap_est = (gl%rhovap_pts(i-2) - gl%rhovap_pts(i-1))/&
                            & (gl%x_pts(i-2,  max_dx_drho(1)) - gl%x_pts(i-1,  max_dx_drho(1)))* &
                            & ( x_calc(max_dx_drho(1))  - gl%x_pts(i-1,  max_dx_drho(1))) + gl%rhovap_pts(i-1)
                            
            ! renorming of the composition vector
            if (abs(sum - 1.d0) > 1.D-10) then
                do k = 1, gl%ncomp
                    if (k/= max_dx_drho(1)) then
                        x_calc(k) = x_calc(k) / sum *(1.D0 - x_calc(max_dx_drho(1)))
                    end if
                end do
            end if     
        end if

        step_ok = .true.
        ! check if step reduction becomes too large
        if (step_reduction > 1.d+4) then
            errval = -5507
            return
        end if
        ! check if one of the components exceeds the bounds 0 <= x <= 1
        do k = 1, gl%ncomp
            if (x_calc(k) <= 0.d0) then
                step_reduction = step_reduction * 2.D0  ! if so, try again to calculate new estimates with reduced step size 
                step_ok = .false.
                exit
            end if
        end do
        !check if temperature or pressure estimates are negative or the estimated pressure step is too large
        if ((Temp < 0.d0) .OR. (press < 0.d0) .OR. (abs(press - gl%p_pts(i-1)) > dpmax)) then
            step_reduction = step_reduction * 2.D0  ! if so, try again to calculate new estimates with reduced step size 
            step_ok = .false.
        end if
    end do
    
    step_ok = .false.
   
    
    !---------------------------------------------------------------  
    ! Step 6: calculate the new point on the phase envelope
    !---------------------------------------------------------------
    ! Try calculation with the Newton-Raphson method if the last component mol fraction is not too small
    ! if the last mol fraction xN is very small the Newton method is likely to fail, due to xN = 1 - sum(xk), which may be rounded to zero
    if (x_calc(gl%ncomp) > 1.d-14) then
        !call Flash_PhaseBoundary_calc_Vbased(press,Temp,z,x_vap,x_calc,rhovap_est,rholiq_est,vapfrac,iFlash, 0, max_dx_drho(1),errval,iter)
        if(NumVer == 1)then
            call Flash_PhaseBoundary_Vbased_calc(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)         !max_dx_drho(1) = Nrxgiven, sensitivste Variable wird festgehalten
                 
        else if (NumVer == 2)then                                                                                                                           ! Vektor aus syntax gründen, steht immer nur eine komponente drin
            call Flash_PhaseBoundary_Vbased_calc_LevMar(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)  ! z = xknown, x_liq = x_calc
        else
            call Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)
        end if
        gl%rho_vap = rhovap_est
        gl%rho_liq = rholiq_est
    else
        errval = -1 ! set the error to nonzero, so the successive substitution method is used
    end if
    ! if one of the mol fractions is very small and/or the Newton method fails, try the slower successive substitution method
    if ((minval(x_calc(1:gl%ncomp)) <= 1.d-12) .and. (errval /= 0)) then
        errval = 0
        press = rhovap_est**3*p_coeffs(1) + rhovap_est**2*p_coeffs(2) + rhovap_est*p_coeffs(3) + p_coeffs(4)
        Temp = rhovap_est**3*T_coeffs(1) + rhovap_est**2*T_coeffs(2) + rhovap_est*T_coeffs(3) + T_coeffs(4)
        sum = 0.d0
        do k = 1, gl%ncomp 
            x_calc(k) = rhovap_est**3*x_coeffs(k,1) + rhovap_est**2*x_coeffs(k,2) &
                    & + rhovap_est*x_coeffs(k,3) + x_coeffs(k,4)
            sum = sum + x_calc(k)
        end do
        x_calc = x_calc/sum
        x_vap = z
        call Succ_Sub(gl,press, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 4, errval, 100, converged)
        ! the successive substitution method is much faster if the temperature is set instead of the pressure (iFlash = 4)
        ! however, this may fail at or near the temperature maximum (cricondentherm). In this case try again with given pressure (iFlash = 2)
        if (.not. converged) then
            press = rhovap_est**3*p_coeffs(1) + rhovap_est**2*p_coeffs(2) + rhovap_est*p_coeffs(3) + p_coeffs(4)
            sum = 0.d0
            do k = 1, gl%ncomp 
                x_calc(k) = rhovap_est**3*x_coeffs(k,1) + rhovap_est**2*x_coeffs(k,2) &
                        & + rhovap_est*x_coeffs(k,3) + x_coeffs(k,4)
                sum = sum + x_calc(k)
            end do
            x_calc = x_calc/sum
            call Succ_Sub(gl,press, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 2, errval, 100, converged)
        end if
        if (.not. converged) errval = -2222
        iter = 3
    end if
    
    !---------------------------------------------------------------  
    ! Step 7: Errorhandling
    !---------------------------------------------------------------
    if (errval /= 0) then
        errval = -5508
        return
    else    
        ! write the last point into the module variables
        gl%T_pts(i) = Temp
        gl%P_pts(i) = press
        gl%rhovap_pts(i) = gl%rho_vap
        gl%rholiq_pts(i) = gl%rho_liq
        gl%x_pts(i, 1:gl%ncomp) = x_calc(1:gl%ncomp)
        gl%phasenv_pts = i
        ! count one up
        i = i + 1
    end if

    !---------------------------------------------------------------  
    ! Step 8: Check if the specified pressure / temperature is between the last two steps. 
    !---------------------------------------------------------------
    if ((gl%P_pts(i-1) - p_spec)*(gl%P_pts(i-2) - p_spec) < 0.d0) then
        ! If so, calculate the phase boundary at the specified pressure, using mean values of the last two steps as initial guesses
        Temp = (gl%T_pts(i-1) + gl%T_pts(i-2))/2.d0
        rhovap_est = (gl%rhovap_pts(i-1) + gl%rhovap_pts(i-2))/2.d0
        rholiq_est = (gl%rholiq_pts(i-1) + gl%rholiq_pts(i-2))/2.d0
        do k = 1, gl%ncomp
            x_calc(k) = (gl%x_pts(i-1, k) + gl%x_pts(i-2, k))/2.d0
        end do
        !call Flash_PhaseBoundary_calc(p_spec, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 2, 0, Nr_x_given, errval, iter)
        if(NumVer == 1)then
            call Flash_PhaseBoundary_Vbased_calc(gl,p_spec, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)         !max_dx_drho(1) = Nrxgiven, sensitivste Variable wird festgehalten
        else if (NumVer == 2)then                                                                                                                           ! Vektor aus syntax gründen, steht immer nur eine komponente drin
            call Flash_PhaseBoundary_Vbased_calc_LevMar(gl,p_spec, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)  ! z = xknown, x_liq = x_calc
        else
            call Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,p_spec, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)
        end if
        gl%rho_vap = rhovap_est
        gl%rho_liq = rholiq_est
        
        if ((minval(x_calc(1:gl%ncomp)) < 1.d-12) .and. (errval /= 0)) then
            errval = 0
            press = p_spec
            Temp = (gl%T_pts(i-1) + gl%T_pts(i-2))/2.d0
            rhovap_est = (gl%rhovap_pts(i-1) + gl%rhovap_pts(i-2))/2.d0
            rholiq_est = (gl%rholiq_pts(i-1) + gl%rholiq_pts(i-2))/2.d0
            do k = 1, gl%ncomp
                x_calc(k) = (gl%x_pts(i-1, k) + gl%x_pts(i-2, k))/2.d0
            end do
            x_vap = z
            call Succ_Sub(gl,press, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 2, errval, 100, converged)
            if (.not. converged) errval = -2222
            iter = 3
        end if
        if (errval /= 0) then
            errval = -5509
            return
        end if
        points_found = points_found + 1
        pt_return(points_found) = Temp
        rho_return(points_found,1) = gl%rho_liq
        rho_return(points_found,2) = gl%rho_vap
        x_return(1:gl%ncomp, points_found) = x_calc(1:gl%ncomp)
        ! overwrite the last values in the arrays
        press = p_spec
        gl%T_pts(i-1) = Temp
        gl%P_pts(i-1) = press
        gl%rhovap_pts(i-1) = gl%rho_vap
        gl%rholiq_pts(i-1) = gl%rho_liq
        gl%x_pts(i-1, 1:gl%ncomp) = x_calc(1:gl%ncomp)
        gl%pointID(i-1) = 5
    end if  
    if ((gl%T_pts(i-1) - T_spec)*(gl%T_pts(i-2) - T_spec) < 0.d0) then
        ! If so, calculate the phase boundary at the specified pressure, using mean values of the last two steps as initial guesses
        press = (gl%p_pts(i-1) + gl%p_pts(i-2))/2.d0
        rhovap_est = (gl%rhovap_pts(i-1) + gl%rhovap_pts(i-2))/2.d0
        rholiq_est = (gl%rholiq_pts(i-1) + gl%rholiq_pts(i-2))/2.d0
        do k = 1, gl%ncomp
            x_calc(k) = (gl%x_pts(i-1, k) + gl%x_pts(i-2, k))/2.d0
        end do
         if(NumVer == 1)then
            call Flash_PhaseBoundary_Vbased_calc(gl,press, T_spec, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)         !max_dx_drho(1) = Nrxgiven, sensitivste Variable wird festgehalten
        else if (NumVer == 2)then                                                                                                                           ! Vektor aus syntax gründen, steht immer nur eine komponente drin
            call Flash_PhaseBoundary_Vbased_calc_LevMar(gl,press, T_spec, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)  ! z = xknown, x_liq = x_calc
        else
            call Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,press, T_spec, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)
        end if
        gl%rho_vap = rhovap_est
        gl%rho_liq = rholiq_est
        
        if ((minval(x_calc(1:gl%ncomp)) < 1.d-12) .and. (errval /= 0)) then
            errval = 0
            press = (gl%p_pts(i-1) + gl%p_pts(i-2))/2.d0
            Temp = T_spec
            rhovap_est = (gl%rhovap_pts(i-1) + gl%rhovap_pts(i-2))/2.d0
            rholiq_est = (gl%rholiq_pts(i-1) + gl%rholiq_pts(i-2))/2.d0
            do k = 1, gl%ncomp
                x_calc(k) = (gl%x_pts(i-1, k) + gl%x_pts(i-2, k))/2.d0
            end do
            x_vap = z
            call Succ_Sub(gl,press, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 4, errval, 100, converged)
            if (.not. converged) errval = -2222
            iter = 3
        end if
        if (errval /= 0) then
            errval = -5510
            return
        end if
        points_found = points_found + 1
        pt_return(points_found) = press
        rho_return(points_found,1) = gl%rho_liq
        rho_return(points_found,2) = gl%rho_vap
        x_return(1:gl%ncomp, points_found) = x_calc(1:gl%ncomp)
        ! overwrite the last values in the arrays
        Temp = T_spec
        gl%T_pts(i-1) = Temp
        gl%P_pts(i-1) = press
        gl%rhovap_pts(i-1) = gl%rho_vap
        gl%rholiq_pts(i-1) = gl%rho_liq
        gl%x_pts(i-1, 1:gl%ncomp) = x_calc(1:gl%ncomp)
        gl%pointID(i-1) = 4
    end if  
    
    !---------------------------------------------------------------  
    ! Step 9: Check if the cricondentherm / cricondenbar is between the last steps. 
    !---------------------------------------------------------------  
    ! check for cricondentherm
    if ((cricondentherm) .and. (FlashSpec <= 0))then
        !the temperature maximum is in the vicinity of point i-3 on the phase boundary
        ! search for the maximum with the cubic interpolation, using the last 4 points
        !--------------------------------------------------------------- 
        ! Determine the temperature maximum using the polynomial approximation 
        ! T(rho_vap) = rho_vap**3*T_coeffs(1) + rho_vap**2*T_coeffs(2) + rho_vap*T_coeffs(3) + T_coeffs(4)
        ! in order to find the maximum, solve the equation rho^2 + p*rho + q = 0
        p = 2.d0/3.d0*T_coeffs(2)/T_coeffs(1)
        q = T_coeffs(3)/T_coeffs(1)/3.d0
        if (4*q <= p*p) then
            rho_vap_new = 0.5d0*(-p + sqrt(p*p - 4.d0*q)) ! first solution
            ! check first solution for maximum
            if ((6.d0*T_coeffs(1)*rho_vap_new + 2.d0*T_coeffs(2)) >= 0) then ! no maximum
                rho_vap_new = 0.5d0*(-p - sqrt(p*p - 4.d0*q)) ! second solution
            end if
        end if
        ! continue only if the solution is within the density interval of the last 3 points
        if ((rho_vap_new - gl%rhovap_pts(i-4))*(rho_vap_new - gl%rhovap_pts(i-2)) < 0.d0) then
            ! create initial estimates using the density found above
            sum = 0.d0
            do k = 1, gl%ncomp 
                x_calc(k) = rho_vap_new**3*x_coeffs(k,1) + rho_vap_new**2*x_coeffs(k,2) &
                        & + rho_vap_new*x_coeffs(k,3) + x_coeffs(k,4)
                if (x_calc(k) < 0.d0) x_calc(k) = 1.d-10
                sum = sum + x_calc(k)
            end do
            press = rho_vap_new**3*p_coeffs(1) + rho_vap_new**2*p_coeffs(2) + rho_vap_new*p_coeffs(3) + p_coeffs(4)
            Temp = rho_vap_new**3*T_coeffs(1) + rho_vap_new**2*T_coeffs(2) + rho_vap_new*T_coeffs(3) + T_coeffs(4)
            rhovap_est = rho_vap_new**3*rhovap_coeffs(1) + rho_vap_new**2*rhovap_coeffs(2) + rho_vap_new*rhovap_coeffs(3)&
                        &  + rhovap_coeffs(4)
            rholiq_est = rho_vap_new**3*rholiq_coeffs(1) + rho_vap_new**2*rholiq_coeffs(2) + rho_vap_new*rholiq_coeffs(3)&
                        &  + rholiq_coeffs(4)
            ! renorming of the composition vector
            if (sum /= 1.D0) then
                do k = 1, gl%ncomp
                    x_calc(k) = x_calc(k) / sum 
                end do
            end if
            ! calculate the flash on the phase boundary with given pressure according to press = p(Tmax)
            !call Flash_PhaseBoundary_calc(press, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 2, 0, Nr_x_given, errval, iter)
             if(NumVer == 1)then
                call Flash_PhaseBoundary_Vbased_calc(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 2, max_dx_drho(1), errval, iter)         !max_dx_drho(1) = Nrxgiven, sensitivste Variable wird festgehalten
            else if (NumVer == 2)then                                                                                                                           ! Vektor aus syntax gründen, steht immer nur eine komponente drin
                call Flash_PhaseBoundary_Vbased_calc_LevMar(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 2, max_dx_drho(1), errval, iter)  ! z = xknown, x_liq = x_calc
            else
                call Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 2, max_dx_drho(1), errval, iter)
            end if
            gl%rho_vap = rhovap_est
            gl%rho_liq = rholiq_est
        
            if (errval /= 0) then
                ! do nothing, this point is not very important
                gl%pointID(i-3) = -2 ! this indicates: cricondentherm somewhere there, but not found
            else    
                ! overwrite the i-2 point, which was the one with the maximum temp before (but not with dT/dp = 0)
                gl%x_pts(i-3, 1:gl%ncomp) = x_calc(1:gl%ncomp)
                gl%T_pts(i-3) = Temp
                gl%P_pts(i-3) = press
                gl%rhovap_pts(i-3) = gl%rho_vap
                gl%rholiq_pts(i-3) = gl%rho_liq
                gl%pointID(i-3) = 2
                ! write the results of point i-1 back into the working variables for the next calculation step
                x_calc(1:gl%ncomp) = gl%x_pts(i-1, 1:gl%ncomp)
                Temp = gl%T_pts(i-1)
                press = gl%p_pts(i-1)
                gl%rho_vap = gl%rhovap_pts(i-1)
                gl%rho_liq = gl%rholiq_pts(i-1)
            end if
        end if
        cricondentherm = .false.
    end if
    if ((gl%T_pts(i-1) - gl%T_pts(i-2))*(gl%T_pts(i-2) - gl%T_pts(i-3)) < 0.d0) cricondentherm = .true.
    ! check for cricondenbar
    if ((cricondenbar) .and. (FlashSpec <= 0)) then
        !the pressure maximum is in the vicinity of point i-3 on the phase boundary
        !--------------------------------------------------------------- 
        ! Determine the pressure maximum using the polynomial approximation 
        ! p(rho_vap) = rho_vap**3*p_coeffs(1) + rho_vap**2*p_coeffs(2) + rho_vap*p_coeffs(3) + p_coeffs(4)
        ! in order to find the maximum, solve the equation rho^2 + p*rho + q = 0
        p = 2.d0/3.d0*p_coeffs(2)/p_coeffs(1)
        q = p_coeffs(3)/p_coeffs(1)/3.d0
        if (4*q <= p*p) then
            rho_vap_new = 0.5d0*(-p + sqrt(p*p - 4.d0*q)) ! first solution
            ! check first solution for maximum
            if ((6.d0*p_coeffs(1)*rho_vap_new + 2.d0*p_coeffs(2)) >= 0) then ! no maximum
                rho_vap_new = 0.5d0*(-p - sqrt(p*p - 4.d0*q)) ! second solution
            end if
        end if
        ! continue only if the solution is within the density interval of the last 3 points
        if ((rho_vap_new - gl%rhovap_pts(i-4))*(rho_vap_new - gl%rhovap_pts(i-2)) < 0.d0) then
            ! create initial estimates using the density found above
            sum = 0.d0
            do k = 1, gl%ncomp 
                x_calc(k) = rho_vap_new**3*x_coeffs(k,1) + rho_vap_new**2*x_coeffs(k,2) &
                        & + rho_vap_new*x_coeffs(k,3) + x_coeffs(k,4)
                if (x_calc(k) < 0.d0) x_calc(k) = 1.d-10
                sum = sum + x_calc(k)
            end do
            press = rho_vap_new**3*p_coeffs(1) + rho_vap_new**2*p_coeffs(2) + rho_vap_new*p_coeffs(3) + p_coeffs(4)
            Temp = rho_vap_new**3*T_coeffs(1) + rho_vap_new**2*T_coeffs(2) + rho_vap_new*T_coeffs(3) + T_coeffs(4)
            rhovap_est = rho_vap_new**3*rhovap_coeffs(1) + rho_vap_new**2*rhovap_coeffs(2) + rho_vap_new*rhovap_coeffs(3)&
                        &  + rhovap_coeffs(4)
            rholiq_est = rho_vap_new**3*rholiq_coeffs(1) + rho_vap_new**2*rholiq_coeffs(2) + rho_vap_new*rholiq_coeffs(3)&
                        &  + rholiq_coeffs(4)
            ! renorming of the composition vector
            if (sum /= 1.D0) x_calc = x_calc / sum 
            
            ! calculate the flash on the phase boundary with given pressure according to Temp = T(pmax)
             if(NumVer == 1)then
                call Flash_PhaseBoundary_Vbased_calc(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 4, max_dx_drho(1), errval, iter)         !max_dx_drho(1) = Nrxgiven, sensitivste Variable wird festgehalten
            else if (NumVer == 2)then                                                                                                                           ! Vektor aus syntax gründen, steht immer nur eine komponente drin
                call Flash_PhaseBoundary_Vbased_calc_LevMar(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 4, max_dx_drho(1), errval, iter)  ! z = xknown, x_liq = x_calc
            else
                call Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 4, max_dx_drho(1), errval, iter)
            end if
            gl%rho_vap = rhovap_est
            gl%rho_liq = rholiq_est
            if (errval /= 0) then
                ! do nothing, this point is not very important
                gl%pointID(i-3) = -3 ! this indicates: cricondentherm somewhere there, but not found
            else    
                ! overwrite the i-2 point, which was the one with the maximum temp before (but not with dT/dp = 0)
                gl%x_pts(i-3, 1:gl%ncomp) = x_calc(1:gl%ncomp)
                gl%T_pts(i-3) = Temp
                gl%P_pts(i-3) = press
                gl%rhovap_pts(i-3) = gl%rho_vap
                gl%rholiq_pts(i-3) = gl%rho_liq
                gl%pointID(i-3) = 3
                ! write the results of point i-1 back into the working variables for the next calculation step
                x_calc(1:gl%ncomp) = gl%x_pts(i-1, 1:gl%ncomp)
                Temp = gl%T_pts(i-1)
                press = gl%p_pts(i-1)
                gl%rho_vap = gl%rhovap_pts(i-1)
                gl%rho_liq = gl%rholiq_pts(i-1)
            end if
        end if
        cricondenbar = .false.
    end if
    if ((gl%P_pts(i-1) - gl%P_pts(i-2))*(gl%P_pts(i-2) - gl%P_pts(i-3)) < 0.d0) cricondenbar = .true.

    !---------------------------------------------------------------  
    ! Step 10: Check if the crtical point is between the last steps. 
    !---------------------------------------------------------------  
    if (((gl%rhovap_pts(i-2) - gl%rholiq_pts(i-2))*(gl%rhovap_pts(i-3) - gl%rholiq_pts(i-3)) < 0.d0) .and. (FlashSpec == 0)) then
        ! interpolate the crit. point using the cubic equation x1_liq(rhovap) from the initial estimates generation
        ! and the relation: x1_liq(rhovap) - z1 = rhovap^3*A + rhovap^2*B + rhovap*C + D - z1 = 0
        !-----------------------------------------------------------
        ! generate the three coefficients for the form of the cubic equation: y(x) =x^3 + a*x^2 + b*x + c   
        rho_crit = 0.d0
        para(1) = x_coeffs(1,2)/x_coeffs(1,1)
        para(2) = x_coeffs(1,3)/x_coeffs(1,1)
        para(3) = (x_coeffs(1,4) - z(1))/x_coeffs(1,1)
        ! call the routine that calculates the roots of the cubic equation
        call cubic_nt(gl,para,root)
        ! find the correct root
        do k = 1, 3
            if ((root(k) > 0.d0) .and. ((gl%rhovap_pts(i-2) - root(k))*(gl%rhovap_pts(i-3) - root(k)) < 0.d0)) rho_crit = root(k)
        end do
        if (rho_crit /= 0.d0) then
            ! calculate the estimated critical parameters from the cubic interpolation and rho_crit
            p_crit = rho_crit**3*p_coeffs(1) + rho_crit**2*p_coeffs(2) + rho_crit*p_coeffs(3) + p_coeffs(4)
            T_crit = rho_crit**3*T_coeffs(1) + rho_crit**2*T_coeffs(2) + rho_crit*T_coeffs(3) + T_coeffs(4)
            sum = 0.d0
            do k = 1, gl%ncomp 
                x_crit(k) = rho_crit**3*x_coeffs(k,1) + rho_crit**2*x_coeffs(k,2) + rho_crit*x_coeffs(k,3) + x_coeffs(k,4)
                if (x_crit(k) < 0.d0) x_crit(k) = 1.d-10
                sum = sum + x_crit(k)
            end do
            x_crit = x_crit/sum
            ! move the last two points one up
            do k = 1,2
                gl%x_pts(i-k+1, 1:gl%ncomp) = gl%x_pts(i-k, 1:gl%ncomp)
                gl%T_pts(i-k+1) = gl%T_pts(i-k)
                gl%p_pts(i-k+1) = gl%p_pts(i-k)
                gl%rhovap_pts(i-k+1) = gl%rhovap_pts(i-k)
                gl%rholiq_pts(i-k+1) = gl%rholiq_pts(i-k)
            end do
            ! overwrite the i-2 point with the crit. point
            gl%x_pts(i-2, 1:gl%ncomp) = x_crit(1:gl%ncomp)
            gl%T_pts(i-2) = T_crit
            gl%P_pts(i-2) = p_crit
            gl%rhovap_pts(i-2) = rho_crit
            gl%rholiq_pts(i-2) = rho_crit
            gl%pointID(i-2) = 1
            gl%phasenv_pts = i
            i = i + 1
        end if
    end if
    
    !---------------------------------------------------------------  
    ! Step 11: Break criteria: 
    !       No. 1: minimum pressure for "normal" phase envelopes, 
    !       No. 2: maximum pressure for open phase envelopes
    !       No. 3: exit if the sepcified point is found and the routine is called from the phase boundary flash algorithm
    !       No. 4: exit if one of the following variables does not change any more: T, p, rho_vap   --> added 08.2012, J.Gernert
    !---------------------------------------------------------------
    if ((press < p_min) .OR. (press > p_max)) return
    if ((FlashSpec > 0) .and. (points_found > 0)) return
    if ((abs(gl%p_pts(i-1)-gl%p_pts(i-2)) < 1.d-5) .OR. (abs(gl%T_pts(i-1)-gl%T_pts(i-2)) < 1.d-5) &
                                & .OR. (abs(gl%rhovap_pts(i-1)-gl%rhovap_pts(i-2)) < 1.d-5)) then
        errval = -5511
        return
    end if
end do

end subroutine phasenv_Vbased
!**************************************************************************


    
!**************************************************************************
!subroutine PTDIAG_OUT_Vbased_mix(fluids, moles, path, p_spec, T_spec,t_pts_out, p_pts_out, rholiq_pts_out, rhovap_pts_out, &
!                        & pointID_pts_out, x_pts_out, fileout, EOS_indicator, errval)
!!**************************************************************************
!! Interface function for plotting a pt-diagram for mixtures
!!!DEC$ ATTRIBUTES DLLEXPORT :: PTDIAG_OUT
!


!

!implicit none
!
!    type(type_gl) :: gl
!
!
!!!! if no pressure / temperature is specified it has to be set to zero !!!
!double precision:: p_spec, T_spec                           ! specified pressure / temperature
!integer, intent(inout):: errval                             ! error handling
!character(255), intent(in):: fileout    ! optional path for output file
!character (255) :: EOS_indicator
!
!!Output, maximum: 1000 pts
!!double precision, dimension(200),intent(out):: t_pts_out, p_pts_out, rholiq_pts_out, rhovap_pts_out
!!integer, dimension(200):: pointID_pts_out
!!double precision, dimension(200, 30),intent(out):: x_pts_out
!double precision, dimension(400),intent(out):: t_pts_out, p_pts_out, rholiq_pts_out, rhovap_pts_out
!integer, dimension(400):: pointID_pts_out
!double precision, dimension(400, 30),intent(out):: x_pts_out
!
!double precision, dimension(6):: pt_return                  ! return value(s) (temp. or press) at the specified pressure/temperature
!double precision, dimension(6,2):: rho_return               ! rho_liq (i,1) and rho_vap (i,2) at the specified pressure/temperature
!double precision, dimension(30,6):: x_return                ! return values for the compositions of the emerging phase at specified pressure/temperature
!integer:: points_found, trunc_fac, i, k
!
!
!!-------------------------------------------------------------------
!!Initialize
!character (12) :: input
!double precision :: prop1, prop2
!character (255) :: fluids
!character (255) :: moles
!character (255) :: path
!double precision, dimension(30):: molev              ! composition of the system
!!double precision, dimension(200):: p_points, T_points, rhovap_points, rholiq_points
!!Andreas August 2013
!double precision, dimension(400):: p_points, T_points, rhovap_points, rholiq_points
!integer:: points
!character (255):: filepure
!!Dummy values for input, Temp and press
!input = 'tp'
!prop1 = 300.D0
!prop2 = 0.1D0
!
!!call setup
!molev = 0.d0
!errval = 0
!call setup (gl,input, prop1, prop2, fluids, moles, molev, path, EOS_indicator, errval)
!!------------------------------------------------------------------
!if (errval == 0) then
!    t_pts_out = 0.D0
!    p_pts_out = 0.D0
!    rholiq_pts_out = 0.D0
!    rhovap_pts_out = 0.D0
!    pointID_pts_out = 0
!    x_pts_out = 0.D0
!    filepure=''
!else
!    t_pts_out = errval
!    p_pts_out = errval
!    rholiq_pts_out = errval
!    rhovap_pts_out = errval
!    pointID_pts_out = errval
!    x_pts_out = errval
!    return
!end if
!
!if (gl%ncomp == 1) then
!    call SATPLOT(gl,p_spec, T_spec, p_points, T_points, rhovap_points, rholiq_points, points, filepure, errval)
!    t_pts_out = T_points
!    p_pts_out = p_points
!    rholiq_pts_out = rholiq_points
!    rhovap_pts_out = rhovap_points
!    pointID_pts_out = gl%pointID(1:400)
!    return
!else
!    call ptdiag_Vbased_mix(gl,molev, p_spec, T_spec, pt_return, rho_return, x_return, points_found, errval)
!end if
!
!!Leave some space in the return array for specified and special points like critical pts, cricondenbar, cricondentherm etc.
!if (gl%phasenv_pts > 0) then
!    if (gl%phasenv_pts > 385) then
!        trunc_fac = ceiling(gl%phasenv_pts / 385.D0)
!    !if (phasenv_pts > 185) then
!    !   trunc_fac = ceiling(phasenv_pts / 185.D0)
!        k = 1
!        t_pts_out(k) = gl%T_pts(1)
!        p_pts_out(k) = gl%p_pts(1)
!        rholiq_pts_out(k) = gl%rholiq_pts(1)
!        rhovap_pts_out(k) = gl%rhovap_pts(1)
!        pointID_pts_out(k) = gl%pointID(1)
!        x_pts_out(k,:) = gl%x_pts(1,:)  
!        k = k + 1      
!        Do i = 1, gl%phasenv_pts
!            if(mod(i,trunc_fac) == 0) then
!                t_pts_out(k) = gl%T_pts(i)
!                p_pts_out(k) = gl%p_pts(i)
!                rholiq_pts_out(k) = gl%rholiq_pts(i)
!                rhovap_pts_out(k) = gl%rhovap_pts(i)
!                pointID_pts_out(k) = gl%pointID(i)
!                x_pts_out(k,:) = gl%x_pts(i,:)  
!                k = k + 1       
!            End if
!            if (k == 400) return
!            !if (k == 200) return
!            if(gl%pointID(i) /= 0) then
!                t_pts_out(k) = gl%T_pts(i)
!                p_pts_out(k) = gl%p_pts(i)
!                rholiq_pts_out(k) = gl%rholiq_pts(i)
!                rhovap_pts_out(k) = gl%rhovap_pts(i)
!                pointID_pts_out(k) = gl%pointID(i)
!                x_pts_out(k,:) = gl%x_pts(i,:)  
!                k = k + 1       
!            End if
!            if (k == 400) return   
!           ! if (k == 200) return 
!        End do
!    else
!        t_pts_out = gl%T_pts(1:400)
!        p_pts_out = gl%p_pts(1:400)
!        rholiq_pts_out = gl%rholiq_pts(1:400)
!        rhovap_pts_out = gl%rhovap_pts(1:400)
!        pointID_pts_out = gl%pointID(1:400)
!        x_pts_out = gl%x_pts(1:400,:)  
!        !t_pts_out = T_pts(1:200)
!        !p_pts_out = p_pts(1:200)
!        !rholiq_pts_out = rholiq_pts(1:200)
!        !rhovap_pts_out = rhovap_pts(1:200)
!        !pointID_pts_out = pointID(1:200)
!        !x_pts_out = x_pts(1:200,:)  
!    End if
!else 
!    t_pts_out = errval
!    p_pts_out = errval
!    rholiq_pts_out = errval
!    rhovap_pts_out = errval
!    pointID_pts_out = errval
!    x_pts_out = errval 
!End if
!
!! output file format
!1000 format(f10.5, f10.4, f10.3, f10.3)
!
!if (fileout /= '') then
!    open(unit = 13,file=fileout, status='unknown', action='write', iostat=errval)
!    if (errval /= 0) return
!    write(13,*)'!---------------------------------------------------------------------------------'
!    write(13,*)'! P-T DIAGRAM DATA FILE'
!    write(13,*)'!'
!    write(13,*)'!---------------------------------------------------------------------------------'
!    write(13,*)'!'
!    write(13,*)'! p/MPa    ','T/K       ','   dl/mol/l ', 'dv/mol/l'
!    do i = 1, 400
!        write(13,1000) p_pts_out(i), t_pts_out(i), rholiq_pts_out(i)/1000.d0, rhovap_pts_out(i)/1000.d0
!    end do
!    close(13)
!end if
!
!end subroutine PTDIAG_OUT_Vbased_mix
!!**************************************************************************


!**************************************************************************
subroutine ptdiag_Vbased_mix(gl,z, p_spec, T_spec, pt_return, rho_return, x_return, points_found, errval)
!**************************************************************************
! This routine attempts to construct the whole phase envelope for a given mixture.
!   -----------------------------------------------------------------------
! The points on the phase envelope are saved in the module following variables of the module module_VLE, 
! starting with the dew line at low pressure and following the phase envelope past the crit. point until the end of the bubble line is reached.
!    double precision, dimension(imax)::T_pts, p_pts, rhovap_pts, rholiq_pts    --> temperatures, pressures and densities of the exisiting phase (rhovap_pts) and the emerging phase (rholiq_pts)
!    double precision, dimension(imax, 30):: x_pts  --> the composition of the emerging phase
!    integer, dimension(imax)::pointID --> IDs of special points: 0: normal points, 1: crit. point(s), 2: temp. maximum, 3: press. maximum, 4: specified temp, 5: specified press.
!    imax is the maximum number of points that can be stored in the variables and is set in the module module_VLE
! with p_spec and T_spec the user can specify points on the phase envelope, for which the unknown variables are calculated and returned.
!   -----------------------------------------------------------------------
! The construction of the phase envelope starts with the bubble line at low pressures.
! If the calculation of the bubble line succeeds until the vicinity of the crit. point, it continues with the dew line calculation
! if the bubble line calculation fails, the system might have an open phase envelope. The algorithm starts with the dew line 
! and attempts to "jump" over the crit. point to continue with the dew line (or LLE line).
! The crit. point(s) is/are interpolated from adjoining points on the phase envelope.
! They are thus not exact, but usually very accurate estimates, since the interpolation is done with points very close to the crit. point.
! The max. pressure /max. temperature points are constructed by estimating p_max or T_max from cubic interpolations,
! then calculating T(p_max) or p(T_max) also from cubic interpolation, and using those values for a T,x" flash or p,x" flash, respectively.
! Thus, those points are results of an actual EOS flash calculation, but only (again very accurate) approximations of the real maxima.
!   -----------------------------------------------------------------------
! J. Gernert and A. Jäger, August 2012




implicit none

    type(type_gl) :: gl


double precision, dimension(30):: z              ! composition of the system
!!! if no pressure / temperature is specified it has to be set to zero !!!
double precision:: p_spec, T_spec                           ! specified pressure / temperature
double precision, dimension(6):: pt_return                  ! return value(s) (temp. or press) at the specified pressure/temperature
double precision, dimension(6,2):: rho_return               ! rho_liq (i,1) and rho_vap (i,2) at the specified pressure/temperature
double precision, dimension(30,6):: x_return                ! return values for the compositions of the emerging phase at specified pressure/temperature
integer:: points_found                         ! number of points found at specified temperature / pressure
integer:: errval                             ! error handling

integer:: FlashSpec                                         ! indicates if the bubble line (1,3) or the dew line (2, 4) needs to be calculated
double precision, dimension(imax):: t_save, p_save, rholiq_save, rhovap_save 
integer, dimension(imax):: pointID_save
double precision, dimension(imax, 30):: x_save
double precision:: p_crit, T_crit, rho_crit, x_crit(30), para(3), root(3), sum
integer:: bubblepoints, dewpoints, i, k
integer::a(4)
double precision, dimension(4):: p_coeffs, T_coeffs, rholiq_coeffs, rhovap_coeffs
double precision, dimension(30, 4):: x_coeffs

t_save=0.d0; p_save=0.d0; rhovap_save=0.d0; rholiq_save=0.d0; x_save = 0.d0; bubblepoints=0
! first try to calculate the bubble line
Flashspec = -1
call phasenv_Vbased_mix(gl,z, p_spec, T_spec, pt_return, rho_return, x_return, points_found, FlashSpec, errval)
! save the phase envelope points for later use
if (gl%phasenv_pts /= 0) then
    bubblepoints = gl%phasenv_pts
    do i = 1, gl%phasenv_pts
        pointID_save(i) = gl%pointID(gl%phasenv_pts+1-i)
        t_save(i) = gl%T_pts(gl%phasenv_pts+1-i)
        p_save(i) = gl%p_pts(gl%phasenv_pts+1-i)
        rholiq_save(i) = gl%rholiq_pts(gl%phasenv_pts+1-i)
        rhovap_save(i) = gl%rhovap_pts(gl%phasenv_pts+1-i)
        x_save(i,1:gl%ncomp) = gl%x_pts(gl%phasenv_pts+1-i,1:gl%ncomp)
    end do
    FlashSpec = 0
    if (errval == 0) FlashSpec = -2 ! only calculate the dew line until the crit. point
else
    FlashSpec = 0   ! try to jump over the crit. point
end if

errval = 0
gl%T_pts= 0.D0
gl%p_pts = 0.D0
gl%rholiq_pts = 0.D0
gl%rhovap_pts = 0.D0
gl%x_pts = 0.D0
gl%pointID = 0
gl%phasenv_pts = 0
! try to calculate the dew line. Jump over the crit. point, if FlashSpec = 0
call phasenv_Vbased_mix(gl,z, p_spec, T_spec, pt_return, rho_return, x_return, points_found, FlashSpec, errval)
dewpoints = gl%phasenv_pts
if ((dewpoints + bubblepoints) > imax) dewpoints = imax - bubblepoints
! write all the points of the dew and bubble line into the module variables
gl%phasenv_pts = dewpoints+bubblepoints
do i = dewpoints+1, gl%phasenv_pts
    gl%pointID(i) = pointID_save(i-dewpoints)
    gl%T_pts(i) = t_save(i-dewpoints)
    gl%p_pts(i) = p_save(i-dewpoints)
    gl%rholiq_pts(i) = rholiq_save(i-dewpoints)
    gl%rhovap_pts(i) = rhovap_save(i-dewpoints)
    gl%x_pts(i,1:gl%ncomp) = x_save(i-dewpoints,1:gl%ncomp)
end do

! interpolate the critical point in case the phase envelope was constructed with the bubble- and the dew line separately
if ((FlashSpec == -2) .and. (errval == 0)) then
   do i = 6, gl%phasenv_pts
       if ((gl%rhovap_pts(i) - gl%rholiq_pts(i))*(gl%rhovap_pts(i-1) - gl%rholiq_pts(i-1)) < 0.d0) then
            ! interpolate the crit. point using the cubic equation x1_liq(rhovap) from the initial estimates generation
            ! and the relation: x1_liq(rhovap) - z1 = rhovap^3*A + rhovap^2*B + rhovap*C + D - z1 = 0
                !---------------------------------------------------------------
            ! make a polynomial expansion for each variable in the form: 
            ! X(rho_vap) = A*rho_vap^3 + B*rho_vap^2 + C*rho_vap + D
            !---------------------------------------------------------------
            ! This expansion is used to calculate the estimated critical point
            ! set up a linear set of equations from the last four points of each variable    
            a = (/i+5, i, i-1, i-5/)
            call Cubic_Coeffs(gl,34, a, p_coeffs, T_coeffs, x_coeffs, rholiq_coeffs, rhovap_coeffs, errval)

            !-----------------------------------------------------------
            ! generate the three coefficients for the form of the cubic equation: y(x) =x^3 + a*x^2 + b*x + c   
            rho_crit = 0.d0
            para(1) = x_coeffs(1,2)/x_coeffs(1,1)
            para(2) = x_coeffs(1,3)/x_coeffs(1,1)
            para(3) = (x_coeffs(1,4) - z(1))/x_coeffs(1,1)
            x_crit = 0.d0
            ! call the routine that calculates the roots of the cubic equation
            call cubic_nt(gl,para,root)
            ! find the correct root
            do k = 1, 3
                if ((root(k) > 0.d0) .and. ((gl%rhovap_pts(i) - root(k))*(gl%rhovap_pts(i-1) - root(k)) < 0.d0)) rho_crit = root(k)
            end do
            if (rho_crit /= 0.d0) then
                ! calculate the estimated critical parameters from the cubic interpolation and rho_crit
                p_crit = rho_crit**3*p_coeffs(1) + rho_crit**2*p_coeffs(2) + rho_crit*p_coeffs(3) + p_coeffs(4)
                T_crit = rho_crit**3*T_coeffs(1) + rho_crit**2*T_coeffs(2) + rho_crit*T_coeffs(3) + T_coeffs(4)
                sum = 0.d0
                do k = 1, gl%ncomp 
                    x_crit(k) = rho_crit**3*x_coeffs(k,1) + rho_crit**2*x_coeffs(k,2) + rho_crit*x_coeffs(k,3) + x_coeffs(k,4)
                    if (x_crit(k) < 0.d0) x_crit(k) = 1.d-10
                    sum = sum + x_crit(k)
                end do
                x_crit = x_crit/sum
                ! move all points after the critical point one up

                ! overwrite the i point with the crit. point
                gl%x_pts(i, 1:gl%ncomp) = x_crit(1:gl%ncomp)
                gl%T_pts(i) = T_crit
                gl%P_pts(i) = p_crit
                gl%rhovap_pts(i) = rho_crit
                gl%rholiq_pts(i) = rho_crit
                gl%pointID(i) = 1
            end if
            exit
        end if
    end do
end if  

end subroutine ptdiag_Vbased_mix
!**************************************************************************

!**************************************************************************
subroutine phasenv_Vbased_mix(gl,z, p_spec, T_spec, pt_return, rho_return, x_return, points_found, FlashSpec, errval)
!**************************************************************************
! This routine calculates the phase envelope of a mixture at given composition.
! The calculated points are saved in the following module variables of the module module_VLE:
!   double precision, dimension(imax):: T_pts, p_pts, rholiq_pts, rhovap_pts
!   double precision, dimension(imax, 31):: x_pts
! For PhaseSpec = 1 or 3:
!   The bubble line is calculated until the specified pressure/temperature is found twice 
!   (retrograde) or the crit. point is reached
! For PhaseSpec = 2 or 4:
!   The dew line is calculated until the specified pressure/temperature is found twice 
!   (retrograde) or the crit. point is reached
! For PhaseSpec = -1 or -2:
!   The bubble line (-1) or dew line (-2) is calculated until the critical point, incl. cricondenbar/cricondentherm 
! For PhaseSpec = 0:
!   Starting on the dew line, either, the maximum number of imax points on the phase boundaries 
!   is calculated or the routine stops when the minimum or maximum pressure on the bubble line is reached
!---------------------------------------------------------------------------
! J. Gernert & A. Jäger, 2010-2012
!---------------------------------------------------------------------------
! Input / Output Variables:
! INPUT:
!   z               feed composition of the mixture
!   p_spec          specified pressure up to which the phase envelope shall be calculated (has to be set to zero if not used!)
!   T_spec          specified temperature up to which the phase envelope shall be calculated (has to be set to zero if not used!)
! OUTPUT
!   pt_return       values for temperature / pressure for the specified value of p_spec / T_spec (vector, length 3)
!   rho_return      values for rho_liq(i,1) and rho_vap(i,2) for the specified value of p_spec / T_spec (matrix, length 3,2)
!   x_return        values for the composition of the emerging phase at the specified value of p_spec / T_spec (matrix, length 30,3)
!   points_found    number of points found at the specified value of p_spec / T_spec (max. value: 3)
!   errval          error value
!---------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl


double precision, dimension(30):: z              ! composition of the system
!!! if no pressure / temperature is specified it has to be set to zero !!!
double precision:: p_spec, T_spec            ! specified pressure / temperature
double precision, dimension(6):: pt_return      ! return value(s) (temp. or press) at the specified pressure/temperature
double precision, dimension(6,2):: rho_return   ! rho_liq (i,1) and rho_vap (i,2) at the specified pressure/temperature
double precision, dimension(30,6):: x_return   ! return values for the compositions of the emerging phase at specified pressure/temperature
integer:: points_found                         ! number of points found at specified temperature / pressure
integer:: errval                             ! error handling
integer:: FlashSpec                             ! indicates if the bubble line (1,3) or the dew line (2, 4) needs to be calculated

double precision, dimension(30)::x_vap, x_liq, x_calc, x_vap_orig, x_calc_orig
double precision, dimension(30, 4):: x_coeffs
double precision, dimension(4):: p_coeffs, T_coeffs, rholiq_coeffs, rhovap_coeffs
double precision:: p_start, press, Temp, vapfrac, delta, sum, press_orig, Temp_orig
integer:: i, k, m, iFlash, iter, Nr_x_given, maxVar
logical:: step_ok, step_critical, jump_over_crit, firstcrit_found, cricondentherm, cricondenbar
character(255):: filename
double precision, dimension(30):: dx_drho
double precision:: dT_drho, dp_drho, delta_x_max, C, B, drho_drho
integer, dimension(1):: max_dx_drho
double precision:: step_reduction, step_crit, step_sign, rhocrit_est
double precision:: rhovap_est, rholiq_est, delta_max, rho_vap_new, p_test, rhovap_est_orig, rholiq_est_orig, GibbsEQN_b
logical:: converged
double precision:: p_min, p_max, dpmax, step_factor(5)
integer, dimension(4):: a
double precision:: rho(5), x_Phase(30,5)
integer:: phasetype(5), nrofphases, count, itercheck(5)
double precision::p, q
double precision :: para(3)                ! a,b,c input parameters for cubic equation root solver
double precision :: root(3)              ! output parameters of cubic equation root solver
double precision:: rho_crit, T_crit, p_crit, x_crit(30)

double precision, dimension(4):: sensitiv
integer, dimension(1):: sens_pos

integer:: NumVer
NumVer = 1



!initialization of variables
x_vap = 0.D0; x_liq = 0.d0; x_calc = 0.d0
Temp = 0.D0; vapfrac = 0.D0
p_coeffs = 0.D0; T_coeffs = 0.D0; x_coeffs = 0.D0

Nr_x_given = 1
step_ok = .false.
step_critical = .false.
firstcrit_found = .false.
cricondentherm = .false.
cricondenbar = .false.
rhovap_est = 0.D0; rholiq_est = 0.D0
jump_over_crit = .false.
p_min = 0.1D0   !Lower pressure limit on the bubble line
p_max = 300.d0 ! upper pressure limit on the bubble line (open phase envelope)
dpmax = 2.d0    ! maximum pressure step
points_found = 0
errval = 0

pt_return = 0.d0; rho_return = 0.d0; x_return = 0.d0
rho_crit = 0.d0; T_crit = 0.d0; p_crit = 0.d0; root = 0.d0; x_crit = 0.d0
itercheck = (/2, 4, 7, 15, 20/) ! step checkpoints for the number of iterations from the previous calculation
step_factor = (/0.5d0, 2.d0, 2.d0, 2.d0, 2.d0/) ! factors for the step reduction, depending on the number of iterations used for the previous step
dx_drho = 0.d0 
 
! if no pressure is specified set p_spec to a very large value to evade wrong if-queries
If (p_spec == 0.d0) p_spec = 1.d6   
If (T_spec == 0.d0) T_spec = 1.d6 
!---------------------------------------------------------------
! Step 1: define a start pressure for the phase envelope
!---------------------------------------------------------------
! normal start pressure is 1 bar
p_start = 0.1D0

! check if start pressure is lower than a ficticious mixture triple point pressure
p_test = 0.d0
do k = 1, gl%ncomp
    p_test = p_test + gl%ptp(k)*z(k)
end do
if (p_start < p_test) p_start = p_test

! check if start pressure is larger than the specified pressure
if (p_start > p_spec) p_start = p_spec*0.9

press = p_start

iFlash = FlashSpec
if ((iFlash == 4) .OR. (iFlash == 0) .OR. (iFlash == -2)) iFlash = 2
if ((iFlash == 3) .OR. (iFlash == -1)) iFlash = 1
if ((iFlash < 1) .OR. (iFlash > 2)) then
    errval = -5500
    return
end if

!---------------------------------------------------------------
! Step 2: create four points at low pressures
!---------------------------------------------------------------
! start with pressure-controlled steps on the dew line
i = 1
! try to generate a starting point

call Flash_PhaseBoundary_calc(gl,press, Temp, z, x_vap, x_liq, 0.d0, 0.d0, vapfrac, iFlash, 5, &
    & Nr_x_given, errval,iter)  !changed the IPhase_try from 0 to 5 (assumption: VLE) - J.G. 11/2012
if ((((minval(x_liq(1:gl%ncomp)) < 1.d-12) .OR. (minval(x_vap(1:gl%ncomp)) < 1.d-12)) .and. (errval /= 0)) &
            & .OR. (dabs(gl%rho_vap - gl%rho_liq) < 1.d-5)) then
    errval = 0; Temp = 0.d0; x_vap = 0.d0; x_liq = 0.d0
    call Succ_Sub(gl,press, Temp, z, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, errval, 100, converged)
    if (.not. converged .OR. ((gl%rho_vap - gl%rho_liq)/gl%rho_liq < 1.d-2)) errval = -2222
end if
if (iFlash == 1) x_calc = x_vap
if (iFlash == 2) x_calc = x_liq
! check if the flash calculation succeeded
if ((errval /= 0) .OR. (dabs(gl%rho_vap - gl%rho_liq) < 1.d-5)) then !generation of starting point failed. Try with Stability analysis, using a temperature from Rachford Rice Eqn. as input
    Temp = 0.d0
    x_calc = 0.d0; x_vap = 0.d0; x_liq = 0.d0
    nrofphases = 1
    count = 0
    call PTX_startvals_PhaseBoundary (gl,press, Temp, z, x_vap, x_liq, vapfrac, iFlash, errval)
    if (errval /= 0) then
        errval = -5501
        return
    end if
    do while (nrofphases < 2)
        count = count + 1
        call PhaseDet2(gl,press, Temp, z, rho, x_Phase, phasetype, vapfrac, nrofphases, errval)
        if (errval /= 0) nrofphases = 1
        if ((nrofphases < 2) .and. (abs(iFlash) == 2)) Temp = Temp*0.98d0
        if ((nrofphases < 2) .and. (abs(iFlash) == 1)) Temp = Temp/0.98d0
        if (count > 30) then ! Generation of start value failed!   !changed criterion from count > 3 to count > 30 since for very hydrogen rich mixtures ptx_startvals generates a nonsense temperature estimate. 
            errval = -5502
            return
        end if
    end do

    ! Use the calculated phase composition and the Temperature as initial values for the p,x flash
    if (iFlash == 1) x_vap = x_Phase(:,phasetype(iFlash))
    if (iFlash == 2) x_liq = x_Phase(:,phasetype(iFlash))
    rhovap_est = gl%rho_vap
    rholiq_est = gl%rho_liq
    call Flash_PhaseBoundary_calc(gl,press, Temp, z, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, 5, &
        & Nr_x_given, errval,iter)  !changed the IPhase_try from 0 to 5 (assumption: VLE) - J.G. 11/2012
    if ((errval /= 0) .OR. (dabs(gl%rho_vap - gl%rho_liq) < 1.d-5)) then ! Generation of start value failed!
        if (iFlash == 1) x_vap = x_Phase(:,phasetype(iFlash))
        if (iFlash == 2) x_liq = x_Phase(:,phasetype(iFlash))
        call Succ_Sub(gl,press, Temp, z, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, errval, 100, converged)
        if (.not. converged .OR. ((gl%rho_vap - gl%rho_liq)/gl%rho_liq < 1.d-2)) then
            errval =  -5503
            return
        else
            errval = 0
        end if
    end if
end if

!write the results in the module variables
gl%T_pts(i) = Temp
gl%P_pts(i) = press
gl%phasenv_pts = i
if (iFlash == 1) then
    gl%x_pts(i, 1:gl%ncomp) = x_vap(1:gl%ncomp)
    gl%rholiq_pts(i) = gl%rho_vap
    gl%rhovap_pts(i) = gl%rho_liq
else
    gl%x_pts(i, 1:gl%ncomp) = x_liq(1:gl%ncomp)
    gl%rholiq_pts(i) = gl%rho_liq
    gl%rhovap_pts(i) = gl%rho_vap
end if

press = press*1.1D0 !increase the pressure for the next step
errval = 0

! generate three more points at low pressures, using the previuos results as initial guesses
do i = 2, 4
    call Flash_PhaseBoundary_calc(gl,press, Temp, z, x_vap, x_liq, gl%rho_vap, gl%rho_liq, vapfrac, iFlash, 5, &
        & Nr_x_given, errval, iter)  !changed the IPhase_try from 0 to 5 (assumption: VLE) - J.G. 11/2012
    if (((minval(x_liq(1:gl%ncomp)) < 1.d-12) .and. (errval /= 0)) .OR. (dabs(gl%rho_vap - gl%rho_liq) < 1.d-5)) then
        !Andreas August 2012
        !errval = 0; Temp = 0.d0; x_vap = 0.d0; x_liq = 0.d0
        !Do not generate new starting values but take the ones from the last iteration
        errval = 0
        Temp = gl%T_pts(i-1)
        x_liq(:) = gl%x_pts(i-1,:)
        call Succ_Sub(gl,press, Temp, z, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, errval, 100, converged)
        if (.NOT. converged) errval = -2222
        iter = 3
    end if
    if (errval /= 0.D0) then
        errval = -5504
        return
    end if
    !write the results in the module variables
    gl%T_pts(i) = Temp
    gl%P_pts(i) = press
    gl%phasenv_pts = i
    if (iFlash == 1) then
        gl%x_pts(i, 1:gl%ncomp) = x_vap(1:gl%ncomp)
        gl%rholiq_pts(i) = gl%rho_vap
        gl%rhovap_pts(i) = gl%rho_liq
    else
        gl%x_pts(i, 1:gl%ncomp) = x_liq(1:gl%ncomp)
        gl%rholiq_pts(i) = gl%rho_liq
        gl%rhovap_pts(i) = gl%rho_vap
    end if
    
    ! check if the last step went over the specified pressure
    if ((gl%p_pts(i) - p_spec)*(gl%p_pts(i-1) - p_spec) < 0.d0) then
        Temp = (gl%T_pts(i) + gl%T_pts(i-1))/2.d0
        x_calc(1:gl%ncomp) = (gl%x_pts(i, 1:gl%ncomp) + gl%x_pts(i-1, 1:gl%ncomp))/2.d0
        if (iFlash == 1) x_vap = x_calc
        if (iFlash == 2) x_liq = x_calc
        call Flash_PhaseBoundary_calc(gl,p_spec, Temp, z, x_vap, x_liq, gl%rho_vap, gl%rho_liq, vapfrac, iFlash, 5, &
            & Nr_x_given, errval, iter)  !changed the IPhase_try from 0 to 5 (assumption: VLE) - J.G. 11/2012
        if (errval /= 0) then
            errval = -5505
            return
        end if
        points_found = points_found + 1
        if (iFlash == 1) then
            gl%x_pts(i, 1:gl%ncomp) = x_vap(1:gl%ncomp)
            gl%rholiq_pts(i) = gl%rho_vap
            gl%rhovap_pts(i) = gl%rho_liq
            rho_return(points_found,1) = gl%rho_vap
            rho_return(points_found,2) = gl%rho_liq
            x_return(1:gl%ncomp, points_found) = x_vap(1:gl%ncomp)
        else
            gl%x_pts(i, 1:gl%ncomp) = x_liq(1:gl%ncomp)
            gl%rholiq_pts(i) = gl%rho_liq
            gl%rhovap_pts(i) = gl%rho_vap
            rho_return(points_found,1) = gl%rho_liq
            rho_return(points_found,2) = gl%rho_vap
            x_return(1:gl%ncomp, points_found) = x_liq(1:gl%ncomp)
        end if
        
        pt_return(points_found) = Temp
        !write the results in the module variables
        gl%T_pts(i) = Temp
        gl%P_pts(i) = press
        gl%phasenv_pts = i
        gl%pointID(i) = 5
    end if
    ! check if the last step went over the specified temperature
    if ((gl%T_pts(i) - T_spec)*(gl%T_pts(i-1) - T_spec) < 0.d0) then
        press = (gl%p_pts(i) + gl%p_pts(i-1))/2.d0
        x_calc(1:gl%ncomp) = (gl%x_pts(i, 1:gl%ncomp) + gl%x_pts(i-1, 1:gl%ncomp))/2.d0
        if (iFlash == 1) x_vap = x_calc
        if (iFlash == 2) x_liq = x_calc
        call Flash_PhaseBoundary_calc(gl,press, T_spec, z, x_vap, x_liq, gl%rho_vap, gl%rho_liq, vapfrac, iFlash+2, 5,&
                & Nr_x_given, errval, iter)  !changed the IPhase_try from 0 to 5 (assumption: VLE) - J.G. 11/2012
        if (errval /= 0) then
            errval = -5506
            return
        end if
        points_found = points_found + 1
        if (iFlash == 1) then
            gl%x_pts(i, 1:gl%ncomp) = x_vap(1:gl%ncomp)
            gl%rholiq_pts(i) = gl%rho_vap
            gl%rhovap_pts(i) = gl%rho_liq
            rho_return(points_found,1) = gl%rho_vap
            rho_return(points_found,2) = gl%rho_liq
            x_return(1:gl%ncomp, points_found) = x_vap(1:gl%ncomp)
        else
            gl%x_pts(i, 1:gl%ncomp) = x_liq(1:gl%ncomp)
            gl%rholiq_pts(i) = gl%rho_liq
            gl%rhovap_pts(i) = gl%rho_vap
            rho_return(points_found,1) = gl%rho_liq
            rho_return(points_found,2) = gl%rho_vap
            x_return(1:gl%ncomp, points_found) = x_liq(1:gl%ncomp)
        end if
        pt_return(points_found) = press
        !write the results in the module variables
        gl%T_pts(i) = Temp
        gl%P_pts(i) = press
        gl%phasenv_pts = i
        gl%pointID(i) = 4
    end if

    press = press*1.1D0 !increase the pressure for the next step
end do

if (iFlash == 1) then
    x_calc(1:gl%ncomp) = x_vap(1:gl%ncomp)
    gl%rho_vap = gl%rhovap_pts(i-1)
    gl%rho_liq = gl%rholiq_pts(i-1)
else
    x_calc(1:gl%ncomp) = x_liq(1:gl%ncomp)
end if
x_vap = z

do while (i < imax)        
    !---------------------------------------------------------------
    ! Step 3: make a polynomial expansion for each variable in the form: 
    ! X(rho_vap) = A*rho_vap^3 + B*rho_vap^2 + C*rho_vap + D
    !---------------------------------------------------------------
    ! This expansion is used to calculate start values for the next point on the phase envelope
    ! set up a linear set of equations from the last four points of each variable    
    a = (/i-1, i-2, i-3, i-4/)
    call Cubic_Coeffs(gl,34, a, p_coeffs, T_coeffs, x_coeffs, rholiq_coeffs, rhovap_coeffs, errval)
    
    !Derivatives of all molefractions with respect to rho_vap
    do k = 1, gl%ncomp 
        dx_drho(k) = 3.d0 * gl%rho_vap**2 * x_coeffs(k,1) + 2.d0 * gl%rho_vap * x_coeffs(k,2) + x_coeffs(k,3)
        dx_drho(k) = dabs(dx_drho(k)) !/ x_pts(i-1,k))
    end do 
    
    !Derivatives of the temperature with respect to rho_vap
    dT_drho = 3.d0 * gl%rho_vap**2 * T_coeffs(1) + 2.d0 * gl%rho_vap * T_coeffs(2) + T_coeffs(3) 
    dT_drho = dabs(dT_drho / Temp)

    !Derivatives of the pressure with respect to rho_vap
    dp_drho = 3.d0 * gl%rho_vap**2 * p_coeffs(1) + 2.d0 * gl%rho_vap * p_coeffs(2) + p_coeffs(3) 
    dp_drho = dabs(dp_drho / press)
    
    !Ableitung der Dichte nach der Dichte immer 1
    drho_drho = 1
    drho_drho = 1 / gl%rho_vap    
        
    !Get the position of the maximum sensitivity of the compositions w.r.t. density
    max_dx_drho = maxloc(abs(dx_drho)) 
    
    !If the last component is the most sensitive, the component with the next largest sensivity is chosen
    if (max_dx_drho(1) == gl%ncomp) then
        dx_drho(max_dx_drho(1)) = 0.D0
        max_dx_drho = maxloc(abs(dx_drho))
    end if
    
    !get the maximum sensitivity from max(dx_drho), dp_drho, dT_drho and drho_drho
        sensitiv(1) = dp_drho
        sensitiv(2) = dT_drho
        sensitiv(3) = dx_drho(max_dx_drho(1))
        sensitiv(4) = drho_drho
    
        !get the position of the highest value of the vector sensitiv
        sens_pos = 0
        sens_pos = maxloc(abs(sensitiv))
        
        !set iFlash depending on the most sensitive variable
        maxVar = 0
        if (sens_pos(1) == 1) then
            maxVar = 31 !the pressure is the most sensitve variable!
            iFlash = 2
        else if (sens_pos(1) == 2) then
            maxVar = 32 !the temperature is the most sensitve variable!
            iFlash = 4
        else if (sens_pos(1) == 3) then
            maxVar = max_dx_drho(1) !one of the compositions is the most sensitvie variable!
            iFlash = 6
        else if (sens_pos(1) == 4) then
            maxVar = 33 !the density is the most sensitive variable!
            iFlash = 8
        end if
        
                !maxVar = 0
                !if (dx_drho(max_dx_drho(1)) > dp_drho) then
                !    if (dx_drho(max_dx_drho(1)) > dT_drho) then
                !        if (dx_drho(max_dx_drho(1)) > ) then
                !        maxVar = max_dx_drho(1) ! one of the compositions is the most sensitvie variable!
                !        iFlash = 6 !set iFlash depending on the most sensitive variable
                !    else
                !        maxVar = 32 ! the temperature is the most sensitve variable!
                !        iFlash = 4 !set iFlash depending on the most sensitive variable
                !    end if
                !else
                !    if (dp_drho > dT_drho) then 
                !        maxVar = 31 ! the pressure is the most sensitve variable!
                !        iFlash = 2 !set iFlash depending on the most sensitive variable
                !    else
                !        maxVar = 32 ! the temperature is the most sensitve variable!
                !        iFlash = 4 !set iFlash depending on the most sensitive variable
                !    end if
                !end if
    
    !---------------------------------------------------------------
    ! Step 4: select the step size for the next pressure / fraction / density step
    !---------------------------------------------------------------
    !The step depends on the iterations necessary for calculating the last point 
    step_reduction = 4.D0
    if (iter < itercheck(1)) step_reduction = step_reduction*step_factor(1)
    if (iter > itercheck(2)) step_reduction = step_reduction*step_factor(2) 
    if (iter > itercheck(3)) step_reduction = step_reduction*step_factor(3)
    if (iter > itercheck(4)) step_reduction = step_reduction*step_factor(4)
    if (iter > itercheck(5)) step_reduction = step_reduction*step_factor(5)

    ! define the maximum step length for the next density step
    if (maxVar < 31) then
        delta_max = 0.05d0 ! maximum composition step
        delta_max = (delta_max/dx_drho(maxVar)) ! calculate a maximum density step from delta_rho_max = delta_x_max/(dx_max_drho)
    else if (maxVar == 31) then 
        delta_max = .5d0    ! maximum pressure step (in MPa)
        delta_max = (delta_max/dp_drho) ! calculate a maximum density step from delta_rho_max = delta_p/(dp_drho)
    else if (maxVar == 32) then
        delta_max = 5.d0    ! maximum temperature step (in K)
        delta_max = (delta_max/dT_drho) ! calculate a maximum density step from delta_rho_max = delta_T/(dT_drho)
    else if (maxVar == 33) then
        delta_max = 50.D0
    end if
    
   do while (.Not.step_ok) 
        ! calculate the standard density step size
        delta = abs(((gl%rho_liq + gl%rho_vap)/2.d0 - gl%rho_vap)/step_reduction)
        !check if standard step size or maximum step size are used
        if (delta > delta_max) delta = delta_max
        !check if step size is too large at small density values
        if (delta > 0.1d0*gl%rho_vap) delta = 0.1d0*gl%rho_vap
        !set the new variable
        if ((FlashSpec == 1) .or. (FlashSpec == 3).or. (FlashSpec == -1)) delta = -delta
        rho_vap_new = gl%rho_vap + delta  
              
        ! check if the density is very close to the critical density
        if (.NOT. firstcrit_found) then 
            ! estimate for the rel. distance to the crit. point of a vapor-liquid equilibrium (1st crit. point):
            ! use the mean of rho_liq and rho_vap as estimate for the crit. density
            step_crit = abs(0.5d0*(gl%rho_liq  - gl%rho_vap)/gl%rho_vap)
        else
            if (rho_crit /= 0.d0) then  ! if the crit. point was interpolated before, use this value
                rhocrit_est = rho_crit
            else
                ! estimate for the rel. distance to the crit. point of a liquid-liquid equilibrium (> 1st crit. point):
                ! use intersection of linear extrapolations of rho_liq1(p) and rho_liq2(p) as estimate for the crit. density
                C = (gl%rholiq_pts(i-1) - gl%rholiq_pts(i-3))/(gl%p_pts(i-1) - gl%p_pts(i-3))
                B = (gl%rhovap_pts(i-1) - gl%rhovap_pts(i-3))/(gl%p_pts(i-1) - gl%p_pts(i-3))
                rhocrit_est = B*((gl%rhovap_pts(i-1) - gl%rholiq_pts(i-1) - (B - C)*gl%p_pts(i-1))/(C - B) - gl%p_pts(i-1)) + gl%rhovap_pts(i-1)
            end if
            step_crit = dabs((rhocrit_est - gl%rho_vap)/gl%rho_vap)
        end if
        
        !check, if the jump over the crit. point was successful
        if (step_critical) then
            ! if the expression (rho1/rho2 - 1) changed signs in the last two steps
            step_sign = sign(1.d0,((gl%rhovap_pts(i-2)/gl%rholiq_pts(i-2) - 1.D0)*(gl%rhovap_pts(i-1)/gl%rholiq_pts(i-1) - 1.D0)))
            if (step_sign < 0.D0) then 
                jump_over_crit = .true.
                firstcrit_found = .true. ! memorize that at least one crit. point was found
            end if
        end if
        
        if (step_crit < 0.05) then ! first jump criterion:  density is very close to crit. density
            delta_x_max = 0.D0
            do k = 1, gl%ncomp ! get the maximum relative distance between the composition and the crit. composition
                if (dabs((x_calc(k) - z(k))/z(k)) > delta_x_max) delta_x_max = dabs((x_calc(k) - z(k))/z(k))
            end do
            if (delta_x_max < 0.05d0) then ! second jump criterion: composition is very close to crit. composition 
                !----------------------------------------------------------------------------
                ! if only the bubble line needs to be calculated, the algorithm exits here!
                !----------------------------------------------------------------------------
                if (FlashSpec /= 0) return
                !----------------------------------------------------------------------------
                !check if the jump over the crit. point was done in one of the previous steps AND the compositions approach the crit. comp.
                if (.NOT. jump_over_crit .and. (dabs(gl%x_pts(i-2, 1) - z(1)) > dabs(gl%x_pts(i-1, 1)- z(1)))) then 
                    x_calc(max_dx_drho(1)) = 2.d0 * z(max_dx_drho(1)) - x_calc(max_dx_drho(1)) 
                    step_ok = .true.
                    step_critical = .true.
                    iFlash = 6
                else
                    step_critical = .false.     !If successfully jumped once, carry on with density steps
                end if          
            end if
        else ! once outside the crit. region, reset the jump criteria
            step_critical = .false.
            jump_over_crit = .false.
        end if

    !---------------------------------------------------------------
    ! Step 5: calculate new starting values for all unknown variables
    !---------------------------------------------------------------
        sum = 0.D0
        if (.NOT. step_critical) then
            do k = 1, gl%ncomp 
                x_calc(k) = rho_vap_new**3*x_coeffs(k,1) + rho_vap_new**2*x_coeffs(k,2) &
                        & + rho_vap_new*x_coeffs(k,3) + x_coeffs(k,4)
                sum = sum + x_calc(k)
            end do
            press = rho_vap_new**3*p_coeffs(1) + rho_vap_new**2*p_coeffs(2) + rho_vap_new*p_coeffs(3) + p_coeffs(4)
            Temp = rho_vap_new**3*T_coeffs(1) + rho_vap_new**2*T_coeffs(2) + rho_vap_new*T_coeffs(3) + T_coeffs(4)
            rhovap_est = rho_vap_new**3*rhovap_coeffs(1) + rho_vap_new**2*rhovap_coeffs(2) + rho_vap_new*rhovap_coeffs(3)&
                        &  + rhovap_coeffs(4)
            rholiq_est = rho_vap_new**3*rholiq_coeffs(1) + rho_vap_new**2*rholiq_coeffs(2) + rho_vap_new*rholiq_coeffs(3)&
                        &  + rholiq_coeffs(4)
            ! renorming of the composition vector
            if (sum /= 1.D0) x_calc(1:gl%ncomp) = x_calc(1:gl%ncomp) / sum 
        else
            ! calculate the new most sensitive x_k from x_pts(i-1, k) + dx/drho(k)*delta_rho - J. Gernert, Mai 2012
            !x_calc(max_dx_drho(1)) = x_calc(max_dx_drho(1)) + dx_drho(max_dx_drho(1))*delta
            sum = sum + x_calc(max_dx_drho(1))
            !Linear extrapolation of all variables in dependence of the most sensitive x 
            do k = 1, gl%ncomp 
                if (k /= max_dx_drho(1)) then
                    x_calc(k) = (gl%x_pts(i-2, k) - gl%x_pts(i-1, k))/ (gl%x_pts(i-2,  max_dx_drho(1)) - gl%x_pts(i-1,  max_dx_drho(1)))* &
                                & ( x_calc(max_dx_drho(1))  - gl%x_pts(i-1,  max_dx_drho(1))) + gl%x_pts(i-1, k)
                    sum = sum + x_calc(k)
                end if
            end do   
            Temp = (gl%T_pts(i-2) - gl%T_pts(i-1))/ (gl%x_pts(i-2,  max_dx_drho(1)) - gl%x_pts(i-1,  max_dx_drho(1)))* &
                            & ( x_calc(max_dx_drho(1))  - gl%x_pts(i-1,  max_dx_drho(1))) + gl%T_pts(i-1)
            press = (gl%p_pts(i-2) - gl%p_pts(i-1))/(gl%x_pts(i-2,  max_dx_drho(1)) - gl%x_pts(i-1,  max_dx_drho(1)))* &
                            & ( x_calc(max_dx_drho(1))  - gl%x_pts(i-1,  max_dx_drho(1))) + gl%p_pts(i-1)
            rholiq_est = (gl%rholiq_pts(i-2) - gl%rholiq_pts(i-1))/&
                            & (gl%x_pts(i-2,  max_dx_drho(1)) - gl%x_pts(i-1,  max_dx_drho(1)))* &
                            & ( x_calc(max_dx_drho(1))  - gl%x_pts(i-1,  max_dx_drho(1))) + gl%rholiq_pts(i-1) 
            rhovap_est = (gl%rhovap_pts(i-2) - gl%rhovap_pts(i-1))/&
                            & (gl%x_pts(i-2,  max_dx_drho(1)) - gl%x_pts(i-1,  max_dx_drho(1)))* &
                            & ( x_calc(max_dx_drho(1))  - gl%x_pts(i-1,  max_dx_drho(1))) + gl%rhovap_pts(i-1)
                            
            ! renorming of the composition vector
            if (abs(sum - 1.d0) > 1.D-10) then
                do k = 1, gl%ncomp
                    if (k/= max_dx_drho(1)) then
                        x_calc(k) = x_calc(k) / sum *(1.D0 - x_calc(max_dx_drho(1)))
                    end if
                end do
            end if     
        end if

        step_ok = .true.
        ! check if step reduction becomes too large
        if (step_reduction > 1.d+4) then
            errval = -5507
            return
        end if
        ! check if one of the components exceeds the bounds 0 <= x <= 1
        do k = 1, gl%ncomp
            if (x_calc(k) <= 0.d0) then
                step_reduction = step_reduction * 2.D0  ! if so, try again to calculate new estimates with reduced step size 
                step_ok = .false.
                exit
            end if
        end do
        !check if temperature or pressure estimates are negative or the estimated pressure step is too large
        if ((Temp < 0.d0) .OR. (press < 0.d0) .OR. (abs(press - gl%p_pts(i-1)) > dpmax)) then
            step_reduction = step_reduction * 2.D0  ! if so, try again to calculate new estimates with reduced step size 
            step_ok = .false.
        end if
    end do
    
    step_ok = .false.
   
    
    !---------------------------------------------------------------  
    ! Step 6: calculate the new point on the phase envelope
    !---------------------------------------------------------------
    
    !save initial values
    press_orig = press
    Temp_orig = Temp
    rhovap_est_orig = rhovap_est
    rholiq_est_orig = rholiq_est
    x_vap_orig = x_vap
    x_calc_orig = x_calc
    
    ! Try calculation with the Newton-Raphson method if the last component mol fraction is not too small
    ! if the last mol fraction xN is very small the Newton method is likely to fail, due to xN = 1 - sum(xk), which may be rounded to zero
    if (x_calc(gl%ncomp) > 1.d-14) then
        call Flash_PhaseBoundary_Vbased_calc(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)
        gl%rho_vap = rhovap_est
        gl%rho_liq = rholiq_est
    else
        errval = -1 ! set the error to nonzero, so the successive substitution method is used
    end if
    ! if one of the mol fractions is very small and/or the Newton method fails, try the slower successive substitution method
    if ((minval(x_calc(1:gl%ncomp)) <= 1.d-12) .and. (errval /= 0)) then
        errval = 0
        press = rhovap_est**3*p_coeffs(1) + rhovap_est**2*p_coeffs(2) + rhovap_est*p_coeffs(3) + p_coeffs(4)
        Temp = rhovap_est**3*T_coeffs(1) + rhovap_est**2*T_coeffs(2) + rhovap_est*T_coeffs(3) + T_coeffs(4)
        sum = 0.d0
        do k = 1, gl%ncomp 
            x_calc(k) = rhovap_est**3*x_coeffs(k,1) + rhovap_est**2*x_coeffs(k,2) &
                    & + rhovap_est*x_coeffs(k,3) + x_coeffs(k,4)
            sum = sum + x_calc(k)
        end do
        x_calc = x_calc/sum
        x_vap = z
        call Succ_Sub(gl,press, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 4, errval, 100, converged)
        ! the successive substitution method is much faster if the temperature is set instead of the pressure (iFlash = 4)
        ! however, this may fail at or near the temperature maximum (cricondentherm). In this case try again with given pressure (iFlash = 2)
        if (.not. converged) then
            press = rhovap_est**3*p_coeffs(1) + rhovap_est**2*p_coeffs(2) + rhovap_est*p_coeffs(3) + p_coeffs(4)
            sum = 0.d0
            do k = 1, gl%ncomp 
                x_calc(k) = rhovap_est**3*x_coeffs(k,1) + rhovap_est**2*x_coeffs(k,2) &
                        & + rhovap_est*x_coeffs(k,3) + x_coeffs(k,4)
                sum = sum + x_calc(k)
            end do
            x_calc = x_calc/sum
            call Succ_Sub(gl,press, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 2, errval, 100, converged)
        end if
        if (.not. converged) errval = -2222
        iter = 3
    end if
    
    !---------------------------------------------------------------
    !If Newton failed, set up Iteration with Levenberg Marquardt again
    !---------------------------------------------------------------
    
    if (errval /= 0) then
        press = press_orig
        Temp = Temp_orig
        rhovap_est = rhovap_est_orig
        rholiq_est = rholiq_est_orig
        x_vap = x_vap_orig
        x_calc = x_calc_orig
        call Flash_PhaseBoundary_Vbased_calc_LevMar_mix(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter, GibbsEQN_b)
        if ((iter > 20) .OR. (GibbsEQN_b > 1.D-8)) then
            call Flash_PhaseBoundary_Vbased_calc(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)
            gl%rho_vap = rhovap_est
            gl%rho_liq = rholiq_est
        end if
    end if 

    !---------------------------------------------------------------  
    ! Step 7: Errorhandling
    !---------------------------------------------------------------
    if (errval /= 0) then
        errval = -5508
        return
    else    
        ! write the last point into the module variables
        gl%T_pts(i) = Temp
        gl%P_pts(i) = press
        gl%rhovap_pts(i) = gl%rho_vap
        gl%rholiq_pts(i) = gl%rho_liq
        gl%x_pts(i, 1:gl%ncomp) = x_calc(1:gl%ncomp)
        gl%phasenv_pts = i
        ! count one up
        i = i + 1
    end if

    !---------------------------------------------------------------  
    ! Step 8: Check if the specified pressure / temperature is between the last two steps. 
    !---------------------------------------------------------------
    if ((gl%P_pts(i-1) - p_spec)*(gl%P_pts(i-2) - p_spec) < 0.d0) then
        ! If so, calculate the phase boundary at the specified pressure, using mean values of the last two steps as initial guesses
        Temp = (gl%T_pts(i-1) + gl%T_pts(i-2))/2.d0
        rhovap_est = (gl%rhovap_pts(i-1) + gl%rhovap_pts(i-2))/2.d0
        rholiq_est = (gl%rholiq_pts(i-1) + gl%rholiq_pts(i-2))/2.d0
        do k = 1, gl%ncomp
            x_calc(k) = (gl%x_pts(i-1, k) + gl%x_pts(i-2, k))/2.d0
        end do
        !call Flash_PhaseBoundary_calc(p_spec, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 2, 0, Nr_x_given, errval, iter)
        if(NumVer == 1)then
            call Flash_PhaseBoundary_Vbased_calc(gl,p_spec, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)         !max_dx_drho(1) = Nrxgiven, sensitivste Variable wird festgehalten
        else if (NumVer == 2)then                                                                                                                           ! Vektor aus syntax gründen, steht immer nur eine komponente drin
            call Flash_PhaseBoundary_Vbased_calc_LevMar(gl,p_spec, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)  ! z = xknown, x_liq = x_calc
        else
            call Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,p_spec, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)
        end if
        gl%rho_vap = rhovap_est
        gl%rho_liq = rholiq_est
        
        if ((minval(x_calc(1:gl%ncomp)) < 1.d-12) .and. (errval /= 0)) then
            errval = 0
            press = p_spec
            Temp = (gl%T_pts(i-1) + gl%T_pts(i-2))/2.d0
            rhovap_est = (gl%rhovap_pts(i-1) + gl%rhovap_pts(i-2))/2.d0
            rholiq_est = (gl%rholiq_pts(i-1) + gl%rholiq_pts(i-2))/2.d0
            do k = 1, gl%ncomp
                x_calc(k) = (gl%x_pts(i-1, k) + gl%x_pts(i-2, k))/2.d0
            end do
            x_vap = z
            call Succ_Sub(gl,press, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 2, errval, 100, converged)
            if (.not. converged) errval = -2222
            iter = 3
        end if
        if (errval /= 0) then
            errval = -5509
            return
        end if
        points_found = points_found + 1
        pt_return(points_found) = Temp
        rho_return(points_found,1) = gl%rho_liq
        rho_return(points_found,2) = gl%rho_vap
        x_return(1:gl%ncomp, points_found) = x_calc(1:gl%ncomp)
        ! overwrite the last values in the arrays
        press = p_spec
        gl%T_pts(i-1) = Temp
        gl%P_pts(i-1) = press
        gl%rhovap_pts(i-1) = gl%rho_vap
        gl%rholiq_pts(i-1) = gl%rho_liq
        gl%x_pts(i-1, 1:gl%ncomp) = x_calc(1:gl%ncomp)
        gl%pointID(i-1) = 5
    end if  
    if ((gl%T_pts(i-1) - T_spec)*(gl%T_pts(i-2) - T_spec) < 0.d0) then
        ! If so, calculate the phase boundary at the specified pressure, using mean values of the last two steps as initial guesses
        press = (gl%p_pts(i-1) + gl%p_pts(i-2))/2.d0
        rhovap_est = (gl%rhovap_pts(i-1) + gl%rhovap_pts(i-2))/2.d0
        rholiq_est = (gl%rholiq_pts(i-1) + gl%rholiq_pts(i-2))/2.d0
        do k = 1, gl%ncomp
            x_calc(k) = (gl%x_pts(i-1, k) + gl%x_pts(i-2, k))/2.d0
        end do
         if(NumVer == 1)then
            call Flash_PhaseBoundary_Vbased_calc(gl,press, T_spec, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)         !max_dx_drho(1) = Nrxgiven, sensitivste Variable wird festgehalten
        else if (NumVer == 2)then                                                                                                                           ! Vektor aus syntax gründen, steht immer nur eine komponente drin
            call Flash_PhaseBoundary_Vbased_calc_LevMar(gl,press, T_spec, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)  ! z = xknown, x_liq = x_calc
        else
            call Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,press, T_spec, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, iFlash, max_dx_drho(1), errval, iter)
        end if
        gl%rho_vap = rhovap_est
        gl%rho_liq = rholiq_est
        
        if ((minval(x_calc(1:gl%ncomp)) < 1.d-12) .and. (errval /= 0)) then
            errval = 0
            press = (gl%p_pts(i-1) + gl%p_pts(i-2))/2.d0
            Temp = T_spec
            rhovap_est = (gl%rhovap_pts(i-1) + gl%rhovap_pts(i-2))/2.d0
            rholiq_est = (gl%rholiq_pts(i-1) + gl%rholiq_pts(i-2))/2.d0
            do k = 1, gl%ncomp
                x_calc(k) = (gl%x_pts(i-1, k) + gl%x_pts(i-2, k))/2.d0
            end do
            x_vap = z
            call Succ_Sub(gl,press, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 4, errval, 100, converged)
            if (.not. converged) errval = -2222
            iter = 3
        end if
        if (errval /= 0) then
            errval = -5510
            return
        end if
        points_found = points_found + 1
        pt_return(points_found) = press
        rho_return(points_found,1) = gl%rho_liq
        rho_return(points_found,2) = gl%rho_vap
        x_return(1:gl%ncomp, points_found) = x_calc(1:gl%ncomp)
        ! overwrite the last values in the arrays
        Temp = T_spec
        gl%T_pts(i-1) = Temp
        gl%P_pts(i-1) = press
        gl%rhovap_pts(i-1) = gl%rho_vap
        gl%rholiq_pts(i-1) = gl%rho_liq
        gl%x_pts(i-1, 1:gl%ncomp) = x_calc(1:gl%ncomp)
        gl%pointID(i-1) = 4
    end if  
    
    !---------------------------------------------------------------  
    ! Step 9: Check if the cricondentherm / cricondenbar is between the last steps. 
    !---------------------------------------------------------------  
    ! check for cricondentherm
    if ((cricondentherm) .and. (FlashSpec <= 0))then
        !the temperature maximum is in the vicinity of point i-3 on the phase boundary
        ! search for the maximum with the cubic interpolation, using the last 4 points
        !--------------------------------------------------------------- 
        ! Determine the temperature maximum using the polynomial approximation 
        ! T(rho_vap) = rho_vap**3*T_coeffs(1) + rho_vap**2*T_coeffs(2) + rho_vap*T_coeffs(3) + T_coeffs(4)
        ! in order to find the maximum, solve the equation rho^2 + p*rho + q = 0
        p = 2.d0/3.d0*T_coeffs(2)/T_coeffs(1)
        q = T_coeffs(3)/T_coeffs(1)/3.d0
        if (4*q <= p*p) then
            rho_vap_new = 0.5d0*(-p + sqrt(p*p - 4.d0*q)) ! first solution
            ! check first solution for maximum
            if ((6.d0*T_coeffs(1)*rho_vap_new + 2.d0*T_coeffs(2)) >= 0) then ! no maximum
                rho_vap_new = 0.5d0*(-p - sqrt(p*p - 4.d0*q)) ! second solution
            end if
        end if
        ! continue only if the solution is within the density interval of the last 3 points
        if ((rho_vap_new - gl%rhovap_pts(i-4))*(rho_vap_new - gl%rhovap_pts(i-2)) < 0.d0) then
            ! create initial estimates using the density found above
            sum = 0.d0
            do k = 1, gl%ncomp 
                x_calc(k) = rho_vap_new**3*x_coeffs(k,1) + rho_vap_new**2*x_coeffs(k,2) &
                        & + rho_vap_new*x_coeffs(k,3) + x_coeffs(k,4)
                if (x_calc(k) < 0.d0) x_calc(k) = 1.d-10
                sum = sum + x_calc(k)
            end do
            press = rho_vap_new**3*p_coeffs(1) + rho_vap_new**2*p_coeffs(2) + rho_vap_new*p_coeffs(3) + p_coeffs(4)
            Temp = rho_vap_new**3*T_coeffs(1) + rho_vap_new**2*T_coeffs(2) + rho_vap_new*T_coeffs(3) + T_coeffs(4)
            rhovap_est = rho_vap_new**3*rhovap_coeffs(1) + rho_vap_new**2*rhovap_coeffs(2) + rho_vap_new*rhovap_coeffs(3)&
                        &  + rhovap_coeffs(4)
            rholiq_est = rho_vap_new**3*rholiq_coeffs(1) + rho_vap_new**2*rholiq_coeffs(2) + rho_vap_new*rholiq_coeffs(3)&
                        &  + rholiq_coeffs(4)
            ! renorming of the composition vector
            if (sum /= 1.D0) then
                do k = 1, gl%ncomp
                    x_calc(k) = x_calc(k) / sum 
                end do
            end if
            ! calculate the flash on the phase boundary with given pressure according to press = p(Tmax)
            !call Flash_PhaseBoundary_calc(press, Temp, z, x_vap, x_calc, rhovap_est, rholiq_est, vapfrac, 2, 0, Nr_x_given, errval, iter)
             if(NumVer == 1)then
                call Flash_PhaseBoundary_Vbased_calc(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 2, max_dx_drho(1), errval, iter)         !max_dx_drho(1) = Nrxgiven, sensitivste Variable wird festgehalten
            else if (NumVer == 2)then                                                                                                                           ! Vektor aus syntax gründen, steht immer nur eine komponente drin
                call Flash_PhaseBoundary_Vbased_calc_LevMar(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 2, max_dx_drho(1), errval, iter)  ! z = xknown, x_liq = x_calc
            else
                call Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 2, max_dx_drho(1), errval, iter)
            end if
            gl%rho_vap = rhovap_est
            gl%rho_liq = rholiq_est
        
            if (errval /= 0) then
                ! do nothing, this point is not very important
                gl%pointID(i-3) = -2 ! this indicates: cricondentherm somewhere there, but not found
            else    
                ! overwrite the i-2 point, which was the one with the maximum temp before (but not with dT/dp = 0)
                gl%x_pts(i-3, 1:gl%ncomp) = x_calc(1:gl%ncomp)
                gl%T_pts(i-3) = Temp
                gl%P_pts(i-3) = press
                gl%rhovap_pts(i-3) = gl%rho_vap
                gl%rholiq_pts(i-3) = gl%rho_liq
                gl%pointID(i-3) = 2
                ! write the results of point i-1 back into the working variables for the next calculation step
                x_calc(1:gl%ncomp) = gl%x_pts(i-1, 1:gl%ncomp)
                Temp = gl%T_pts(i-1)
                press = gl%p_pts(i-1)
                gl%rho_vap = gl%rhovap_pts(i-1)
                gl%rho_liq = gl%rholiq_pts(i-1)
            end if
        end if
        cricondentherm = .false.
    end if
    if ((gl%T_pts(i-1) - gl%T_pts(i-2))*(gl%T_pts(i-2) - gl%T_pts(i-3)) < 0.d0) cricondentherm = .true.
    ! check for cricondenbar
    if ((cricondenbar) .and. (FlashSpec <= 0)) then
        !the pressure maximum is in the vicinity of point i-3 on the phase boundary
        !--------------------------------------------------------------- 
        ! Determine the pressure maximum using the polynomial approximation 
        ! p(rho_vap) = rho_vap**3*p_coeffs(1) + rho_vap**2*p_coeffs(2) + rho_vap*p_coeffs(3) + p_coeffs(4)
        ! in order to find the maximum, solve the equation rho^2 + p*rho + q = 0
        p = 2.d0/3.d0*p_coeffs(2)/p_coeffs(1)
        q = p_coeffs(3)/p_coeffs(1)/3.d0
        if (4*q <= p*p) then
            rho_vap_new = 0.5d0*(-p + sqrt(p*p - 4.d0*q)) ! first solution
            ! check first solution for maximum
            if ((6.d0*p_coeffs(1)*rho_vap_new + 2.d0*p_coeffs(2)) >= 0) then ! no maximum
                rho_vap_new = 0.5d0*(-p - sqrt(p*p - 4.d0*q)) ! second solution
            end if
        end if
        ! continue only if the solution is within the density interval of the last 3 points
        if ((rho_vap_new - gl%rhovap_pts(i-4))*(rho_vap_new - gl%rhovap_pts(i-2)) < 0.d0) then
            ! create initial estimates using the density found above
            sum = 0.d0
            do k = 1, gl%ncomp 
                x_calc(k) = rho_vap_new**3*x_coeffs(k,1) + rho_vap_new**2*x_coeffs(k,2) &
                        & + rho_vap_new*x_coeffs(k,3) + x_coeffs(k,4)
                if (x_calc(k) < 0.d0) x_calc(k) = 1.d-10
                sum = sum + x_calc(k)
            end do
            press = rho_vap_new**3*p_coeffs(1) + rho_vap_new**2*p_coeffs(2) + rho_vap_new*p_coeffs(3) + p_coeffs(4)
            Temp = rho_vap_new**3*T_coeffs(1) + rho_vap_new**2*T_coeffs(2) + rho_vap_new*T_coeffs(3) + T_coeffs(4)
            rhovap_est = rho_vap_new**3*rhovap_coeffs(1) + rho_vap_new**2*rhovap_coeffs(2) + rho_vap_new*rhovap_coeffs(3)&
                        &  + rhovap_coeffs(4)
            rholiq_est = rho_vap_new**3*rholiq_coeffs(1) + rho_vap_new**2*rholiq_coeffs(2) + rho_vap_new*rholiq_coeffs(3)&
                        &  + rholiq_coeffs(4)
            ! renorming of the composition vector
            if (sum /= 1.D0) x_calc = x_calc / sum 
            
            ! calculate the flash on the phase boundary with given pressure according to Temp = T(pmax)
             if(NumVer == 1)then
                call Flash_PhaseBoundary_Vbased_calc(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 4, max_dx_drho(1), errval, iter)         !max_dx_drho(1) = Nrxgiven, sensitivste Variable wird festgehalten
            else if (NumVer == 2)then                                                                                                                           ! Vektor aus syntax gründen, steht immer nur eine komponente drin
                call Flash_PhaseBoundary_Vbased_calc_LevMar(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 4, max_dx_drho(1), errval, iter)  ! z = xknown, x_liq = x_calc
            else
                call Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,press, Temp, rhovap_est, rholiq_est, z, x_vap, x_calc, vapfrac, 4, max_dx_drho(1), errval, iter)
            end if
            gl%rho_vap = rhovap_est
            gl%rho_liq = rholiq_est
            if (errval /= 0) then
                ! do nothing, this point is not very important
                gl%pointID(i-3) = -3 ! this indicates: cricondentherm somewhere there, but not found
            else    
                ! overwrite the i-2 point, which was the one with the maximum temp before (but not with dT/dp = 0)
                gl%x_pts(i-3, 1:gl%ncomp) = x_calc(1:gl%ncomp)
                gl%T_pts(i-3) = Temp
                gl%P_pts(i-3) = press
                gl%rhovap_pts(i-3) = gl%rho_vap
                gl%rholiq_pts(i-3) = gl%rho_liq
                gl%pointID(i-3) = 3
                ! write the results of point i-1 back into the working variables for the next calculation step
                x_calc(1:gl%ncomp) = gl%x_pts(i-1, 1:gl%ncomp)
                Temp = gl%T_pts(i-1)
                press = gl%p_pts(i-1)
                gl%rho_vap = gl%rhovap_pts(i-1)
                gl%rho_liq = gl%rholiq_pts(i-1)
            end if
        end if
        cricondenbar = .false.
    end if
    if ((gl%P_pts(i-1) - gl%P_pts(i-2))*(gl%P_pts(i-2) - gl%P_pts(i-3)) < 0.d0) cricondenbar = .true.

    !---------------------------------------------------------------  
    ! Step 10: Check if the crtical point is between the last steps. 
    !---------------------------------------------------------------  
    if (((gl%rhovap_pts(i-2) - gl%rholiq_pts(i-2))*(gl%rhovap_pts(i-3) - gl%rholiq_pts(i-3)) < 0.d0) .and. (FlashSpec == 0)) then
        ! interpolate the crit. point using the cubic equation x1_liq(rhovap) from the initial estimates generation
        ! and the relation: x1_liq(rhovap) - z1 = rhovap^3*A + rhovap^2*B + rhovap*C + D - z1 = 0
        !-----------------------------------------------------------
        ! generate the three coefficients for the form of the cubic equation: y(x) =x^3 + a*x^2 + b*x + c   
        rho_crit = 0.d0
        para(1) = x_coeffs(1,2)/x_coeffs(1,1)
        para(2) = x_coeffs(1,3)/x_coeffs(1,1)
        para(3) = (x_coeffs(1,4) - z(1))/x_coeffs(1,1)
        ! call the routine that calculates the roots of the cubic equation
        call cubic_nt(gl,para,root)
        ! find the correct root
        do k = 1, 3
            if ((root(k) > 0.d0) .and. ((gl%rhovap_pts(i-2) - root(k))*(gl%rhovap_pts(i-3) - root(k)) < 0.d0)) rho_crit = root(k)
        end do
        if (rho_crit /= 0.d0) then
            ! calculate the estimated critical parameters from the cubic interpolation and rho_crit
            p_crit = rho_crit**3*p_coeffs(1) + rho_crit**2*p_coeffs(2) + rho_crit*p_coeffs(3) + p_coeffs(4)
            T_crit = rho_crit**3*T_coeffs(1) + rho_crit**2*T_coeffs(2) + rho_crit*T_coeffs(3) + T_coeffs(4)
            sum = 0.d0
            do k = 1, gl%ncomp 
                x_crit(k) = rho_crit**3*x_coeffs(k,1) + rho_crit**2*x_coeffs(k,2) + rho_crit*x_coeffs(k,3) + x_coeffs(k,4)
                if (x_crit(k) < 0.d0) x_crit(k) = 1.d-10
                sum = sum + x_crit(k)
            end do
            x_crit = x_crit/sum
            ! move the last two points one up
            do k = 1,2
                gl%x_pts(i-k+1, 1:gl%ncomp) = gl%x_pts(i-k, 1:gl%ncomp)
                gl%T_pts(i-k+1) = gl%T_pts(i-k)
                gl%p_pts(i-k+1) = gl%p_pts(i-k)
                gl%rhovap_pts(i-k+1) = gl%rhovap_pts(i-k)
                gl%rholiq_pts(i-k+1) = gl%rholiq_pts(i-k)
            end do
            ! overwrite the i-2 point with the crit. point
            gl%x_pts(i-2, 1:gl%ncomp) = x_crit(1:gl%ncomp)
            gl%T_pts(i-2) = T_crit
            gl%P_pts(i-2) = p_crit
            gl%rhovap_pts(i-2) = rho_crit
            gl%rholiq_pts(i-2) = rho_crit
            gl%pointID(i-2) = 1
            gl%phasenv_pts = i
            i = i + 1
        end if
    end if
    
    !---------------------------------------------------------------  
    ! Step 11: Break criteria: 
    !       No. 1: minimum pressure for "normal" phase envelopes, 
    !       No. 2: maximum pressure for open phase envelopes
    !       No. 3: exit if the sepcified point is found and the routine is called from the phase boundary flash algorithm
    !       No. 4: exit if one of the following variables does not change any more: T, p, rho_vap   --> added 08.2012, J.Gernert
    !---------------------------------------------------------------
    if ((press < p_min) .OR. (press > p_max)) return
    if ((FlashSpec > 0) .and. (points_found > 0)) return
    if ((abs(gl%p_pts(i-1)-gl%p_pts(i-2)) < 1.d-5) .OR. (abs(gl%T_pts(i-1)-gl%T_pts(i-2)) < 1.d-5) &
                                & .OR. (abs(gl%rhovap_pts(i-1)-gl%rhovap_pts(i-2)) < 1.d-5)) then
        errval = -5511
        return
    end if
end do

    end subroutine phasenv_Vbased_mix
!**************************************************************************

    
!Subroutines für die numerischen Verfahren: Levenberg Marquardt, Dog Leg    

subroutine GibbsEQN_betrag(gl,eqn,GibbsEQN,GibbsEQN_b)


implicit none

    type(type_gl) :: gl


integer:: i, eqn 
double precision:: GibbsEQN_b
double precision, dimension(60):: GibbsEQN

GibbsEQN_b = 0.D0

do i = 1, eqn
    GibbsEQN_b = GibbsEQN_b + GibbsEQN(i)**2
end do

GibbsEQN_b = GibbsEQN_b**0.5

end subroutine GibbsEQN_betrag 


subroutine t_m_rechner(gl,eqn,t,m,GibbsEQN_new,GibbsEQN)


implicit none

    type(type_gl) :: gl


integer:: i,j,eqn
double precision:: t,m
double precision, dimension(60):: GibbsEQN_new, GibbsEQN

t = 0.D0

do i = 1, eqn
    t = t + GibbsEQN_new(i)**2
end do

m = 0.D0

do j = 1, eqn
    m = m + GibbsEQN(j)**2
end do


end subroutine t_m_rechner


subroutine D_rechner(gl,eqn,mu,A,D)


implicit none

    type(type_gl) :: gl


integer:: i, eqn
double precision:: mu
double precision, dimension(60,60):: A,D

D = 0.D0

D = A

do i = 1, eqn
    D(i,i) = D(i,i) + mu
end do

end subroutine


subroutine rho_rechner(gl,eqn,GibbsEQN_b,GibbsEQN_b_new,mu,h_lm,g,rho)


implicit none

    type(type_gl) :: gl


integer::i,eqn
double precision:: a,b,c,rho,GibbsEQN_b,GibbsEQN_b_new,mu
double precision, dimension(60)::h_lm,g

a = 0.D0
b = 0.D0
c = 0.D0

!a:
a = GibbsEQN_b**2

!b:
b = GibbsEQN_b_new**2

!c:
do i = 1, eqn
    c = c + h_lm(i)*(mu*h_lm(i)-g(i))
end do

!rho:
rho = (0.5*(a - b)) / (0.5*c)


end subroutine


subroutine alpha_rechner(gl,eqn,g,J,alpha)
 

implicit none

    type(type_gl) :: gl


integer:: i, eqn, k
double precision:: alpha, a, b
double precision, dimension(60):: g,c
double precision, dimension(60,60):: J

a = 0.D0

do i = 1, eqn
    a = a + g(i)**2
end do

b = 0.D0
c = 0.D0

c = matmul(J,g)

do k = 1, eqn
    b = b + c(k)**2
end do

alpha = 0.D0

alpha = a / b


end subroutine


subroutine h_gn_betrag(gl,eqn,h_gn,h_gn_b)


implicit none

    type(type_gl) :: gl


integer:: i, eqn 
double precision:: h_gn_b
double precision, dimension(60):: h_gn

h_gn_b = 0.D0

do i = 1, eqn
    h_gn_b = h_gn_b + h_gn(i)**2
end do

h_gn_b = h_gn_b**0.5

end subroutine h_gn_betrag 



subroutine h_sd_betrag(gl,eqn,h_sd,h_sd_b)


implicit none

    type(type_gl) :: gl


integer:: i, eqn 
double precision:: h_sd_b
double precision, dimension(60):: h_sd

h_sd_b = 0.D0

do i = 1, eqn
    h_sd_b = h_sd_b + h_sd(i)**2
end do

h_sd_b = h_sd_b**0.5

end subroutine h_sd_betrag 


subroutine h_dl_betrag(gl,eqn,h_dl,h_dl_b)


implicit none

    type(type_gl) :: gl


integer:: i, eqn 
double precision:: h_dl_b
double precision, dimension(60):: h_dl

h_dl_b = 0.D0

do i = 1, eqn
    h_dl_b = h_dl_b + h_dl(i)**2
end do

h_dl_b = h_dl_b**0.5

end subroutine h_dl_betrag 


subroutine rho_rechner_DogLeg(gl,eqn,GibbsEQN_b,GibbsEQN_b_new,h_dl,g,rho)


implicit none

    type(type_gl) :: gl


integer::i,eqn
double precision:: a,b,c,rho,GibbsEQN_b,GibbsEQN_b_new
double precision, dimension(60)::h_lm,g,h_dl



a = 0.D0
b = 0.D0
c = 0.D0

!a:
a = GibbsEQN_b**2

!b:
b = GibbsEQN_b_new**2

!c:
do i = 1, eqn
    c = c + h_dl(i)* (h_dl(i)-g(i))
end do

!rho:
rho = (0.5*(a - b)) / (0.5*c)


end subroutine


subroutine q_rechner(gl,delta, h_dl_b, q)


implicit none

    type(type_gl) :: gl


double precision:: delta, h_dl_b
double precision, dimension(2):: q

q = 0.D0

q(1) = delta
q(2) = 3*h_dl_b

end subroutine 


subroutine beta_rechner(gl,eqn,h_sd,h_gn,delta,beta_lc)


implicit none

    type(type_gl) :: gl


integer:: i,j,k,eqn
double precision:: delta, beta_lc, c, ba, a
double precision, dimension(60)::h_sd,h_gn

c = 0.D0
do j = 1, eqn
    c = c + h_sd(j) * (h_gn(j) - h_sd(j))
end do

ba = 0.D0
do i = 1, eqn
    ba = ba + (h_gn(i) - h_sd(i))**2
end do

a = 0.D0
do k = 1, eqn
    a = a + h_sd(k)**2
end do


if (c <= 0) then
    beta_lc = (-c + (c**2 + ba * (delta**2 - a))**0.5) / ba

else
    beta_lc = (delta**2 - a) / (c + (c**2 + ba * (delta**2 - a))**0.5)

end if

end subroutine



subroutine ohnenamen(gl,eqn,delta,g,g_b,h_dl)


implicit none

    type(type_gl) :: gl


integer:: eqn, i
double precision:: g_b, delta
double precision, dimension(60):: g, h_dl

call g_betrag(gl,eqn,g,g_b)

do i = 1, eqn
h_dl(i) = (delta / g_b) * (-g(i))
end do

end subroutine

subroutine g_betrag(gl,eqn,g,g_b)


implicit none

    type(type_gl) :: gl


integer:: i, eqn 
double precision:: g_b
double precision, dimension(60):: g

g_b = 0.D0

do i = 1, eqn
    g_b = g_b + g(i)**2
end do

g_b = g_b**0.5

end subroutine g_betrag 

subroutine h_lm_betrag(gl,eqn,h_lm,h_lm_b)


implicit none

    type(type_gl) :: gl


integer:: i, eqn 
double precision:: h_lm_b
double precision, dimension(60):: h_lm

h_lm_b = 0.D0

do i = 1, eqn
    h_lm_b = h_lm_b + h_lm(i)**2
end do

h_lm_b = h_lm_b**0.5

end subroutine h_lm_betrag 



    end module phasenv_vbased_module
