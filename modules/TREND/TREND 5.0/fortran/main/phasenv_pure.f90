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

    ! module for file phasenv_pure.f90
    module phasenv_pure_module
    !global use inclusion
    use module_all_types
    use flash_pure_module
    use phasenv_pbased_module
    use crit_tpd_module

    contains

	


!--------------------------------------------------------------------------------
subroutine SATPLOT(gl,p_spec, T_spec, p_points, T_points, rhovap_points, rholiq_points, points, fileout, errval)
!-------------------------------------------------------------------------------------
! THIS SUBROUTINE CALCULATES THE PHASE BOUNDARIES OF A PURE FLUID AT GIVEN TEMPERATURE
! OR PRESSURE. THE RESULTS ARE RETURNED IN VECTORS AND CAN BE PRINTED INTO A FILE
!-------------------------------------------------------------------------------------
!  major changes to the code: cubic extrapolation routines added, step size control etc. - J.G. 2012
! T. Wiens, September 2012





implicit none

    type(type_gl) :: gl

!--------------------------------------------------------------------------------
integer:: iFlash                            ! indicates whether temperature or pressure is used to calculate the VLE, 1: temperature, 2: pressure
integer, parameter:: maxi = 400             ! maximum number of calculated points, lengt of the return vectors!
double precision :: p_points(maxi), T_points(maxi), rhovap_points(maxi), rholiq_points(maxi)    ! return vectors of the calculated points
character(255) :: fileout        ! optional path for output file
integer:: points, errval                    ! number of calculated points, error code
double precision:: T_spec, p_spec           ! user specific points on the saturation line. If not used, these variables must be set to zero in calling routine!
!--------------------------------------------------------------------------------

double precision:: T,p,inc,rhovap_est,rholiq_est, dpdT
integer:: error, i, iter, a(4)
character(255):: errMSG
double precision:: p_coeffs(4), T_coeffs(4), x_coeffs(30,4), rholiq_coeffs(4), rhovap_coeffs(4)
double precision:: rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, pcmax, tcmax
integer :: iterc, errvalc

p_points = 0.d0
T_points = 0.d0
rhovap_points = 0.d0
rholiq_points = 0.d0
points = 0
errval = 0
errMSG = ''

rhovap_est = 0.d0 
rholiq_est = 0.d0

T = 0.d0

if (T_spec == 0.d0) T_spec = 1.d+6
if (p_spec == 0.d0) p_spec = 1.d+6


!special case POE
if (trim(gl%components(1)) == 'pec5' .and. gl%mix_type /= 6) then
    rhoc_est_inp = 550.d0
    call find_crit_tpd(gl,gl%tc(1), rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errvalc, 1)
    pcmax = pc_eos_
    tcmax = tc_eos_

elseif (trim(gl%components(1)) == 'pec7' .and. gl%mix_type /= 6) then
    rhoc_est_inp = 410.d0
    call find_crit_tpd(gl,gl%tc(1), rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iterc, errvalc, 1)
    pcmax = pc_eos_
    tcmax = tc_eos_
else
    pcmax = gl%pc(1)
    tcmax = gl%tc(1)
end if


!--------------------------------------------------------------------------------
! T CONTROLLED STEPS
!--------------------------------------------------------------------------------
i = 1
iFlash = 1
!Changed initial temperature from tminfluid to ttp
!Andreas July 2014
T = gl%ttp(1)*1.0001d0   ! initial temperature
!exception for oil
if (gl%components(1) == "oil") T = 290.D0
inc = 2.d0*(tcmax - T)/maxi  ! initial temperature step
! calculate the first VLE point
call Flash_Pure_PhaseBoundary(gl,p, T, rhovap_est, rholiq_est, iFlash, errval, iter, 1)
if (errval == 0) then
    p_points(i) = p
    T_points(i) = T
    rhovap_points(i) = rhovap_est
    rholiq_points(i) = rholiq_est
    gl%p_pts(i) = p
    gl%T_pts(i) = T
    gl%rhovap_pts(i) = rhovap_est
    gl%rholiq_pts(i) = rholiq_est
    points = i
    gl%phasenv_pts = i
    i = i + 1
else
    if (gl%mix_type == 6) then
        do while (errval /= 0 .and. T < 0.7d0*Tcmax) 
            T = T + 5.d0
            call Flash_Pure_PhaseBoundary(gl,p, T, rhovap_est, rholiq_est, iFlash, errval, iter, 1)
        end do
        if (errval == 0) then
            p_points(i) = p
            T_points(i) = T
            rhovap_points(i) = rhovap_est
            rholiq_points(i) = rholiq_est
            gl%p_pts(i) = p
            gl%T_pts(i) = T
            gl%rhovap_pts(i) = rhovap_est
            gl%rholiq_pts(i) = rholiq_est
            points = i
            gl%phasenv_pts = i
            i = i + 1
        else
            errval = -5555    
            return
        end if
    else
        errval = -5555    
        return
    end if    
end if

T = T + inc ! increase the temperature
dpdT = 0.d0 !slope of the vapor pressure curve
do while (dpdT < 0.05)
    call Flash_Pure_PhaseBoundary(gl,p, T, rhovap_est, rholiq_est, iFlash, errval, iter, 1)
    if (errval == 0) then
        p_points(i) = p
        T_points(i) = T
        rhovap_points(i) = rhovap_est
        rholiq_points(i) = rholiq_est
        gl%p_pts(i) = p
        gl%T_pts(i) = T
        gl%rhovap_pts(i) = rhovap_est
        gl%rholiq_pts(i) = rholiq_est
        points = i
        gl%phasenv_pts = i
    else 
        errval = -5555    
        return
    end if
    inc = 2.d0*(tcmax - T)/(maxi - points)
    T = T + inc
    i = i + 1
    dpdT = (p_points(i-1) - p_points(i-2))/(T_points(i-1) - T_points(i-2))
    ! check if the specified temperature was between the last two calculated points
    if ((T_spec - T_points(i-2))*(T_spec - T_points(i-1)) < 0.d0) then
        call Flash_Pure_PhaseBoundary(gl,p, T_spec, rhovap_est, rholiq_est, iFlash, errval, iter, 1)
        if (errval == 0) then
            ! move the specified point between the two last points
            p_points(i) = p_points(i-1)
            T_points(i) = T_points(i-1)
            rhovap_points(i) = rhovap_points(i-1)
            rholiq_points(i) = rholiq_points(i-1)
            gl%p_pts(i) = p_points(i-1)
            gl%T_pts(i) = T_points(i-1)
            gl%rhovap_pts(i) = rhovap_points(i-1)
            gl%rholiq_pts(i) = rholiq_points(i-1)
            
            p_points(i-1) = p
            T_points(i-1) = T_spec
            rhovap_points(i-1) = rhovap_est
            rholiq_points(i-1) = rholiq_est
            gl%p_pts(i-1) = p
            gl%T_pts(i-1) = T_spec
            gl%rhovap_pts(i-1) = rhovap_est
            gl%rholiq_pts(i-1) = rholiq_est
            gl%pointID(i-1) = 4
            i = i + 1
        else
            gl%pointID(i-1) = -4
        end if
    end if
    ! check if the specified pressure was between the last two calculated points
    if ((p_spec - p_points(i-2))*(p_spec - p_points(i-1)) < 0.d0) then
        T = T - inc
        call Flash_Pure_PhaseBoundary(gl,p_spec, T, rhovap_est, rholiq_est, 2, errval, iter, 1)
        if (errval == 0) then
            ! move the specified point between the two last points
            p_points(i) = p_points(i-1)
            T_points(i) = T_points(i-1)
            rhovap_points(i) = rhovap_points(i-1)
            rholiq_points(i) = rholiq_points(i-1)
            gl%p_pts(i) = p_points(i-1)
            gl%T_pts(i) = T_points(i-1)
            gl%rhovap_pts(i) = rhovap_points(i-1)
            gl%rholiq_pts(i) = rholiq_points(i-1)
            
            p_points(i-1) = p_spec
            T_points(i-1) = T
            rhovap_points(i-1) = rhovap_est
            rholiq_points(i-1) = rholiq_est
            gl%p_pts(i-1) = p_spec
            gl%T_pts(i-1) = T
            gl%rhovap_pts(i-1) = rhovap_est
            gl%rholiq_pts(i-1) = rholiq_est
            gl%pointID(i-1) = 5
            i = i + 1
        else
            gl%pointID(i-1) = -5
        end if
        T = T_points(i-1) + inc        
    end if
    if (points > maxi*0.6) exit ! second exit critetion
end do

!--------------------------------------------------------------------------------
! P CONTROLLED STEPS
!--------------------------------------------------------------------------------
iFlash = 2
inc = 2.d0*(pcmax - p)/(maxi - points)
p = p + inc
do while (p < pcmax*0.9)
    call Flash_Pure_PhaseBoundary(gl,p, T, rhovap_est, rholiq_est, iFlash, errval, iter, 1)
    if (errval == 0) then
        p_points(i) = p
        T_points(i) = T
        rhovap_points(i) = rhovap_est
        rholiq_points(i) = rholiq_est
        gl%p_pts(i) = p
        gl%T_pts(i) = T
        gl%rhovap_pts(i) = rhovap_est
        gl%rholiq_pts(i) = rholiq_est
        points = i
        gl%phasenv_pts = i
    else 
        errval = -5555 
    !    exit ! in this case the points that were already calculated can be printed in a file
        return
    end if
    inc = 2.d0*(pcmax - p)/(maxi - points)
    p = p + inc
    i = i + 1
    ! check if the specified temperature was between the last two calculated points
    if ((T_spec - T_points(i-2))*(T_spec - T_points(i-1)) < 0.d0) then
        p = p - inc
        call Flash_Pure_PhaseBoundary(gl,p, T_spec, rhovap_est, rholiq_est, 1, errval, iter, 1)
        if (errval == 0) then
            ! move the specified point between the two last points
            p_points(i) = p_points(i-1)
            T_points(i) = T_points(i-1)
            rhovap_points(i) = rhovap_points(i-1)
            rholiq_points(i) = rholiq_points(i-1)
            gl%p_pts(i) = p_points(i-1)
            gl%T_pts(i) = T_points(i-1)
            gl%rhovap_pts(i) = rhovap_points(i-1)
            gl%rholiq_pts(i) = rholiq_points(i-1)
            
            p_points(i-1) = p
            T_points(i-1) = T_spec
            rhovap_points(i-1) = rhovap_est
            rholiq_points(i-1) = rholiq_est
            gl%p_pts(i-1) = p
            gl%T_pts(i-1) = T_spec
            gl%rhovap_pts(i-1) = rhovap_est
            gl%rholiq_pts(i-1) = rholiq_est
            gl%pointID(i-1) = 4
            i = i + 1
        else
            gl%pointID(i-1) = -4
        end if
        p = p_points(i-1) + inc
    end if
    ! check if the specified pressure was between the last two calculated points
    if ((p_spec - p_points(i-2))*(p_spec - p_points(i-1)) < 0.d0) then
        call Flash_Pure_PhaseBoundary(gl,p_spec, T, rhovap_est, rholiq_est, iFlash, errval, iter, 1)
        if (errval == 0) then
            ! move the specified point between the two last points
            p_points(i) = p_points(i-1)
            T_points(i) = T_points(i-1)
            rhovap_points(i) = rhovap_points(i-1)
            rholiq_points(i) = rholiq_points(i-1)
            gl%p_pts(i) = p_points(i-1)
            gl%T_pts(i) = T_points(i-1)
            gl%rhovap_pts(i) = rhovap_points(i-1)
            gl%rholiq_pts(i) = rholiq_points(i-1)
            
            p_points(i-1) = p_spec
            T_points(i-1) = T
            rhovap_points(i-1) = rhovap_est
            rholiq_points(i-1) = rholiq_est
            gl%p_pts(i-1) = p_spec
            gl%T_pts(i-1) = T
            gl%rhovap_pts(i-1) = rhovap_est
            gl%rholiq_pts(i-1) = rholiq_est
            gl%pointID(i-1) = 5
            i = i + 1
        else
            gl%pointID(i-1) = -5
        end if
    end if
    if (i > maxi*0.9) exit
end do

!--------------------------------------------------------------------------------
! CALCULATE CRITICAL REGION - D CONTROLLED STEPS
!--------------------------------------------------------------------------------
iFlash = 2
inc = (gl%rhoc(1) - rhovap_points(i-1))/(maxi - points)
do while (i < maxi)    
    a = (/i-1, i-2, i-3, i-4/)
    call Cubic_Coeffs(gl,34, a, p_coeffs, T_coeffs, x_coeffs, rholiq_coeffs, rhovap_coeffs, errval)
    ! set next density step
    rhovap_est = rhovap_est + inc
    ! calculate estimates for all variables
    p = rhovap_est**3*p_coeffs(1) + rhovap_est**2*p_coeffs(2) + rhovap_est*p_coeffs(3) + p_coeffs(4)
    T = rhovap_est**3*t_coeffs(1) + rhovap_est**2*t_coeffs(2) + rhovap_est*t_coeffs(3) + t_coeffs(4)
    rholiq_est = rhovap_est**3*rholiq_coeffs(1) + rhovap_est**2*rholiq_coeffs(2) + rhovap_est*rholiq_coeffs(3)&
                & + rholiq_coeffs(4)
    
    ! calculate next point only if estimates are good
    if ((p < pcmax) .AND. (T < tcmax) .AND. (rholiq_est > gl%rhoc(1)) .AND. (rholiq_est /= 0.d0) .AND. (errval == 0)) then
        call Flash_Pure_PhaseBoundary(gl,p, T, rhovap_est, rholiq_est, iFlash, errval, iter, 1)
        if (errval == 0) then
            p_points(i) = p
            T_points(i) = T
            rhovap_points(i) = rhovap_est
            rholiq_points(i) = rholiq_est
            gl%p_pts(i) = p
            gl%T_pts(i) = T
            gl%rhovap_pts(i) = rhovap_est
            gl%rholiq_pts(i) = rholiq_est
            points = i
            gl%phasenv_pts = i
            i = i + 1
            inc = (gl%rhoc(1) - rhovap_points(i-1))/(maxi - points)
        end if
    else
        errval = -1
    end if
    
    ! if something went wrong with te calculation, try again with smaller density step
    if (errval < 0) then
        rhovap_est = rhovap_est - inc
        inc = inc/2.d0
    end if
    
    ! exit if step size becomes too small
    if ((inc < 1.d-4) .and. (gl%components(1) /= "oil")) then
        errval = -3333
   !     exit   ! in this case the points that were already calculated can be printed in a file
    else if (inc < 1.d-9) then
        errval = -3333
        exit
    end if
    
    ! check if the specified temperature was between the last two calculated points
    if ((T_spec - T_points(i-2))*(T_spec - T_points(i-1)) < 0.d0) then
        p = p - inc
        call Flash_Pure_PhaseBoundary(gl,p, T_spec, rhovap_est, rholiq_est, 1, errval, iter, 1)
        if (errval == 0) then
            ! move the specified point between the two last points
            p_points(i) = p_points(i-1)
            T_points(i) = T_points(i-1)
            rhovap_points(i) = rhovap_points(i-1)
            rholiq_points(i) = rholiq_points(i-1)
            gl%p_pts(i) = p_points(i-1)
            gl%T_pts(i) = T_points(i-1)
            gl%rhovap_pts(i) = rhovap_points(i-1)
            gl%rholiq_pts(i) = rholiq_points(i-1)
            
            p_points(i-1) = p
            T_points(i-1) = T_spec
            rhovap_points(i-1) = rhovap_est
            rholiq_points(i-1) = rholiq_est
            gl%p_pts(i-1) = p
            gl%T_pts(i-1) = T_spec
            gl%rhovap_pts(i-1) = rhovap_est
            gl%rholiq_pts(i-1) = rholiq_est
            gl%pointID(i-1) = 4
            i = i + 1
        else
            gl%pointID(i-1) = -4
        end if
        p = p_points(i-1) + inc
    end if
    ! check if the specified pressure was between the last two calculated points
    if ((p_spec - p_points(i-2))*(p_spec - p_points(i-1)) < 0.d0) then
        call Flash_Pure_PhaseBoundary(gl,p_spec, T, rhovap_est, rholiq_est, iFlash, errval, iter, 1)
        if (errval == 0) then
            ! move the specified point between the two last points
            p_points(i) = p_points(i-1)
            T_points(i) = T_points(i-1)
            rhovap_points(i) = rhovap_points(i-1)
            rholiq_points(i) = rholiq_points(i-1)
            gl%p_pts(i) = p_points(i-1)
            gl%T_pts(i) = T_points(i-1)
            gl%rhovap_pts(i) = rhovap_points(i-1)
            gl%rholiq_pts(i) = rholiq_points(i-1)
            
            p_points(i-1) = p_spec
            T_points(i-1) = T
            rhovap_points(i-1) = rhovap_est
            rholiq_points(i-1) = rholiq_est
            gl%p_pts(i-1) = p_spec
            gl%T_pts(i-1) = T
            gl%rhovap_pts(i-1) = rhovap_est
            gl%rholiq_pts(i-1) = rholiq_est
            gl%pointID(i-1) = 5
            i = i + 1
        else
            gl%pointID(i-1) = -5
        end if
    end if
end do

!--------------------------------------------------------------------------------
! LAST POINT: CRITICAL PARAMETERS
!--------------------------------------------------------------------------------
p_points(i) = pcmax
T_points(i) = tcmax
rhovap_points(i) = gl%rhoc(1)
rholiq_points(i) = gl%rhoc(1)
gl%p_pts(i) = pcmax
gl%T_pts(i) = tcmax
gl%rhovap_pts(i) = gl%rhoc(1)
gl%rholiq_pts(i) = gl%rhoc(1)
points = i
gl%phasenv_pts = i
gl%pointID(i) = 1


!--------------------------------------------------------------------------------
! CHECK IF WRITING THE DATA TO A FILE IS WISHED
!--------------------------------------------------------------------------------
if (fileout /= '') then
    open(unit = 13, file=fileout, status='unknown', action='write', iostat=error)
    if (error /= 0) return
    write(13,*)'!---------------------------------------------------------------------------------'
    write(13,*)'! SAT DIAGRAM DATA FILE FOR ', trim(gl%components(1))
    write(13,*)'!---------------------------------------------------------------------------------'
    write(13,*)'!'
    write(13,*)'!  p/MPa           ','T/K       ','         dl/mol/l           ', 'dv/mol/l'
    if (gl%components(1) == "oil" .or. gl%components(1) == "pec5" .or. gl%components(1) == "pec7") then
        do i = 1, points
            write(13,1001) p_points(i), T_points(i), rholiq_points(i)/1000.d0, rhovap_points(i)/1000.d0
        end do
    else
        do i = 1, points
            write(13,1000) p_points(i), T_points(i), rholiq_points(i)/1000.d0, rhovap_points(i)/1000.d0
        end do
    end if
    if (errMSG /= '') write(13,*)'! ', errMSG
    close(13)
end if

1000 format(4(f15.9,3x))
1001 format(D15.9,3x,2(f15.9,3x),D15.9)

end subroutine SATPLOT


    end module phasenv_pure_module
