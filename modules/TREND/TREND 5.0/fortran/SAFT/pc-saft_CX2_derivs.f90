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

    ! module for file pc-saft_CX2_derivs.f90
    module pc_saft_CX2_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module


    contains




subroutine CX2DERIVS(gl,GETDERC)
    
! Henning Markgraf, June 2016

    ! C_1: abbreviation for the compressibility expression
    ! defined by eq. A.11 and A.31 in Gross, Sadowski 2001:
    ! C_1 = (1 + mmean*(8*eta-2*eta**2)/(1-eta)**4 + (1-mmean) &
    !       & * (20*eta-27*eta**2+12*eta**3-2*eta**4)/((1-eta)*(2-eta))**2)**(-1)
    ! dependent on T and D






implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERC
    !output: cx2_PCSAFT
    !working variables
    double precision :: recurring_factor
    double precision :: part1, part2, part3, part4, part5, part6
    double precision :: part7, part8, part9, part10, part11, part12, part13, part14
    integer :: i, xi, xj
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! III. recurring_factor occurs in all derivatives
    recurring_factor = 2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 16.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 51.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 70.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 27.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 12.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - gl%z3_PCSAFT(1)** &
            & 6 + 8.d0*gl%z3_PCSAFT(1)**5 - 27.d0*gl%z3_PCSAFT(1)**4 + 42.d0*gl%z3_PCSAFT(1)**3 - 26.d0*gl%z3_PCSAFT(1)**2 + 4.d0
    
    ! IV. X2 DERIVATIVES of C
    !calculate the derivatives C_1, C_2, C_3, ...
    ! 1: C_1
    ! requires z3_PCSAFT(1)
    if (GETDERC(1) .eq. 1) then
        part1 = 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**2/recurring_factor**3
        part2 = (4.d0*gl%z3_PCSAFT(1)**16 - 88.d0*gl%z3_PCSAFT(1)**15 + 896.d0* &
     & gl%z3_PCSAFT(1)**14 - 5552.d0*gl%z3_PCSAFT(1)**13 + 23185.d0*gl%z3_PCSAFT(1)**12 - 68282.d0*gl%z3_PCSAFT(1)** &
     & 11 + 144127.d0*gl%z3_PCSAFT(1)**10 - 216640.d0*gl%z3_PCSAFT(1)**9 + 224163.d0*gl%z3_PCSAFT(1)**8 - &
     & 145938.d0*gl%z3_PCSAFT(1)**7 + 43645.d0*gl%z3_PCSAFT(1)**6 + 8748.d0*gl%z3_PCSAFT(1)**5 - 9708.d0*gl%z3_PCSAFT(1) &
     & **4 + 864.d0*gl%z3_PCSAFT(1)**3 + 576.d0*gl%z3_PCSAFT(1)**2)
        part3 = (2.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13 - 74.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 905.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 5835.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 23098.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 59594.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 &
     & + 101491.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 110777.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 69804.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5 - 16584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 5244.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 2520.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 5.d0*gl%z3_PCSAFT(1)**13 - 65.d0*gl%z3_PCSAFT(1)**12 &
     & + 344.d0*gl%z3_PCSAFT(1)**11 - 818.d0*gl%z3_PCSAFT(1)**10 - 27.d0*gl%z3_PCSAFT(1)**9 + 5879.d0*gl%z3_PCSAFT(1)**8 &
     & - 18410.d0*gl%z3_PCSAFT(1)**7 + 30864.d0*gl%z3_PCSAFT(1)**6 - 32152.d0*gl%z3_PCSAFT(1)**5 + 20876.d0* &
     & gl%z3_PCSAFT(1)**4 - 7600.d0*gl%z3_PCSAFT(1)**3 + 864.d0*gl%z3_PCSAFT(1)**2 + 336.d0*gl%z3_PCSAFT(1) - 96.d0)
	part4 = (2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 74.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 905.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 5835.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 23098.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 &
     & - 59594.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 101491.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 110777.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**6 + 69804.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 16584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 5244.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 2520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 5.d0* &
     & gl%z3_PCSAFT(1)**13 - 65.d0*gl%z3_PCSAFT(1)**12 + 344.d0*gl%z3_PCSAFT(1)**11 - 818.d0*gl%z3_PCSAFT(1)**10 - 27.d0 &
     & *gl%z3_PCSAFT(1)**9 + 5879.d0*gl%z3_PCSAFT(1)**8 - 18410.d0*gl%z3_PCSAFT(1)**7 + 30864.d0*gl%z3_PCSAFT(1)**6 - &
     & 32152.d0*gl%z3_PCSAFT(1)**5 + 20876.d0*gl%z3_PCSAFT(1)**4 - 7600.d0*gl%z3_PCSAFT(1)**3 + 864.d0*gl%z3_PCSAFT(1)** &
     & 2 + 336.d0*gl%z3_PCSAFT(1) - 96.d0)
	part5 = (6.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 12 - 252.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 2971.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 17812.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 65442.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 160048.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**7 + 269557.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 311268.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 &
     & + 233112.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 95280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 8172.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 5040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 576.d0*gl%mmean_PCSAFT**2 + 3.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 114.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 1993.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
     & 13636.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 53901.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 139874.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**7 - 251789.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 319716.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - &
     & 282112.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 163912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 55140.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**2 + 7152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 528.d0*gl%mmean_PCSAFT - 3.d0*gl%z3_PCSAFT(1)**12 + 6.d0* &
     & gl%z3_PCSAFT(1)**11 + 256.d0*gl%z3_PCSAFT(1)**10 - 2392.d0*gl%z3_PCSAFT(1)**9 + 10737.d0*gl%z3_PCSAFT(1)**8 - &
     & 30390.d0*gl%z3_PCSAFT(1)**7 + 59442.d0*gl%z3_PCSAFT(1)**6 - 83304.d0*gl%z3_PCSAFT(1)**5 + 84096.d0* &
     & gl%z3_PCSAFT(1)**4 - 59928.d0*gl%z3_PCSAFT(1)**3 + 28632.d0*gl%z3_PCSAFT(1)**2 - 8224.d0*gl%z3_PCSAFT(1) + &
     & 1072.d0)
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
           !     if (xi .GE. xj) then
            	    gl%cx2_PCSAFT(1,xi,xj) = part1*(gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*part2 + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(1,xi)*part3 + &
         & gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(1,xj)*part4 + gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part5)
           !     end if
            end do
        end do
    end if
    
    ! 2: 1ST DERIVATIVE OF C_1 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires z3_PCSAFT(1)
    if (GETDERC(2) .eq. 1) then
        part1 = -2.d0*(gl%z3_PCSAFT(1) - 1.d0)/recurring_factor**4
        part2 = (8.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 368.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 20 + 6552.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 66320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 442018.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 2085220.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 7236292.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**15 - 18829460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 36923072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 13 - 54098356.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 57660540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - &
     & 42021172.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 17650566.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 1056408.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 2766360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 974520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 6 + 106128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 69120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 6912.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**3 + 32.d0*gl%z3_PCSAFT(1)**21 - 752.d0*gl%z3_PCSAFT(1)**20 + 8220.d0*gl%z3_PCSAFT(1)**19 - &
     & 55072.d0*gl%z3_PCSAFT(1)**18 + 250734.d0*gl%z3_PCSAFT(1)**17 - 810928.d0*gl%z3_PCSAFT(1)**16 + &
     & 1878768.d0*gl%z3_PCSAFT(1)**15 - 3004128.d0*gl%z3_PCSAFT(1)**14 + 2771940.d0*gl%z3_PCSAFT(1)**13 + &
     & 338536.d0*gl%z3_PCSAFT(1)**12 - 5894116.d0*gl%z3_PCSAFT(1)**11 + 10693032.d0*gl%z3_PCSAFT(1)**10 - &
     & 11311682.d0*gl%z3_PCSAFT(1)**9 + 7678296.d0*gl%z3_PCSAFT(1)**8 - 3045992.d0*gl%z3_PCSAFT(1)**7 + &
     & 334008.d0*gl%z3_PCSAFT(1)**6 + 276048.d0*gl%z3_PCSAFT(1)**5 - 112704.d0*gl%z3_PCSAFT(1)**4 + 1152.d0* &
     & gl%z3_PCSAFT(1)**3 + 4608.d0*gl%z3_PCSAFT(1)**2)
        part3 = (8.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 20 - 428.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 7604.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 73552.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 460622.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 2032871.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**15 + 6633593.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 16470835.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**13 + 31613041.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 47101293.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**11 + 53945491.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 46032517.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**9 + 27207753.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 9136920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **7 + 236472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 952704.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - &
     & 170856.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 34560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 3456.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**2 + 16.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 136.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 1328.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 29352.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 239044.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **16 + 1193350.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 4155786.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 10718790.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 21211698.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 32987090.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 41005718.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 41133938.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9 - 33173898.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 20899504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - &
     & 9514864.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 2521344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 13872.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**4 - 206208.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 34560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 4608.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 10.d0*gl%z3_PCSAFT(1)**20 + 175.d0*gl%z3_PCSAFT(1)**19 - 1201.d0*gl%z3_PCSAFT(1)**18 &
     & + 2941.d0*gl%z3_PCSAFT(1)**17 + 12153.d0*gl%z3_PCSAFT(1)**16 - 137337.d0*gl%z3_PCSAFT(1)**15 + 631143.d0 &
     & *gl%z3_PCSAFT(1)**14 - 1904329.d0*gl%z3_PCSAFT(1)**13 + 4228309.d0*gl%z3_PCSAFT(1)**12 - 7280186.d0* &
     & gl%z3_PCSAFT(1)**11 + 10017974.d0*gl%z3_PCSAFT(1)**10 - 11199800.d0*gl%z3_PCSAFT(1)**9 + 10201352.d0* &
     & gl%z3_PCSAFT(1)**8 - 7473672.d0*gl%z3_PCSAFT(1)**7 + 4268568.d0*gl%z3_PCSAFT(1)**6 - 1796208.d0* &
     & gl%z3_PCSAFT(1)**5 + 498352.d0*gl%z3_PCSAFT(1)**4 - 64192.d0*gl%z3_PCSAFT(1)**3 - 7488.d0*gl%z3_PCSAFT(1)**2 + &
     & 3840.d0*gl%z3_PCSAFT(1) - 384.d0)
	part4 = (8.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 428.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 7604.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 73552.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17 + 460622.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 2032871.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **15 + 6633593.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 16470835.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 &
     & + 31613041.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 47101293.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + &
     & 53945491.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 46032517.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + &
     & 27207753.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 9136920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 236472.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 952704.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 170856.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**4 - 34560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 3456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + &
     & 16.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 136.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 1328.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 18 + 29352.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 239044.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 1193350.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 4155786.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 10718790.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13 - 21211698.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 32987090.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 11 - 41005718.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 41133938.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - &
     & 33173898.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 20899504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 9514864.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 2521344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 13872.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 &
     & - 206208.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 34560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 4608.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1) - 10.d0*gl%z3_PCSAFT(1)**20 + 175.d0*gl%z3_PCSAFT(1)**19 - 1201.d0*gl%z3_PCSAFT(1)**18 + 2941.d0* &
     & gl%z3_PCSAFT(1)**17 + 12153.d0*gl%z3_PCSAFT(1)**16 - 137337.d0*gl%z3_PCSAFT(1)**15 + 631143.d0*gl%z3_PCSAFT(1) &
     & **14 - 1904329.d0*gl%z3_PCSAFT(1)**13 + 4228309.d0*gl%z3_PCSAFT(1)**12 - 7280186.d0*gl%z3_PCSAFT(1)** &
     & 11 + 10017974.d0*gl%z3_PCSAFT(1)**10 - 11199800.d0*gl%z3_PCSAFT(1)**9 + 10201352.d0*gl%z3_PCSAFT(1)** &
     & 8 - 7473672.d0*gl%z3_PCSAFT(1)**7 + 4268568.d0*gl%z3_PCSAFT(1)**6 - 1796208.d0*gl%z3_PCSAFT(1)**5 + &
     & 498352.d0*gl%z3_PCSAFT(1)**4 - 64192.d0*gl%z3_PCSAFT(1)**3 - 7488.d0*gl%z3_PCSAFT(1)**2 + 3840.d0*gl%z3_PCSAFT(1) &
     & - 384.d0)
	part5 = (24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 1464.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 25160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 227264.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**16 + 1321810.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 5488510.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**14 + 17317066.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 43098806.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**12 + 86028830.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 137489810.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**10 + 173020246.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 - 166440834.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**8 + 116897544.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - 55436568.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**6 + 15028488.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 1264680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **4 - 117792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 69120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 6912.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) + 1440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 30192.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17 + 292332.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 1754682.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **15 + 7405254.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 23620002.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 &
     & + 59526630.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 121349070.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + &
     & 201260226.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 269470806.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + &
     & 286443222.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 236174904.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + &
     & 146415528.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 65189880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + &
     & 19202184.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 3033312.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 26496.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 44928.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 4608.d0*gl%mmean_PCSAFT**2 - 18.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 342.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 11394.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 &
     & - 123306.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 776586.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 3365838.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 10940970.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 28101198.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12 + 58762224.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 101235060.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 10 + 143413836.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 165277728.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + &
     & 152487312.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 110322480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 60872400.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 24575520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 6767040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **3 - 1099008.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 64512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 4224.d0*gl%mmean_PCSAFT + &
     & 6.d0*gl%z3_PCSAFT(1)**19 - 6.d0*gl%z3_PCSAFT(1)**18 - 1294.d0*gl%z3_PCSAFT(1)**17 + 16978.d0*gl%z3_PCSAFT(1)**16 &
     & - 114842.d0*gl%z3_PCSAFT(1)**15 + 517550.d0*gl%z3_PCSAFT(1)**14 - 1729130.d0*gl%z3_PCSAFT(1)**13 + &
     & 4551286.d0*gl%z3_PCSAFT(1)**12 - 9780484.d0*gl%z3_PCSAFT(1)**11 + 17449312.d0*gl%z3_PCSAFT(1)**10 - &
     & 25906160.d0*gl%z3_PCSAFT(1)**9 + 31784848.d0*gl%z3_PCSAFT(1)**8 - 31843552.d0*gl%z3_PCSAFT(1)**7 + &
     & 25658336.d0*gl%z3_PCSAFT(1)**6 - 16317536.d0*gl%z3_PCSAFT(1)**5 + 7985984.d0*gl%z3_PCSAFT(1)**4 - &
     & 2898112.d0*gl%z3_PCSAFT(1)**3 + 734080.d0*gl%z3_PCSAFT(1)**2 - 115840.d0*gl%z3_PCSAFT(1) + 8576.d0)
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
            	    gl%cx2_PCSAFT(2,xi,xj) = part1*(gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*part2 + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(1,xi)*part3 + & 
			    & gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(1,xj)*part4 + gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part5)
                end if
            end do
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF C_1 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires z3_PCSAFT(1)
    if (GETDERC(3) .eq. 1) then
        part1 = 4.d0/recurring_factor**5
        part2 = (24.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**28 - 1440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27 + &
     & 33580.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 - 451720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 4088254.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 26946048.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 135507961.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 535343930.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 1692677365.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 4331702180.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 9019781350.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 15278925396.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + &
     & 20919521480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 22808767440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & + 19238718649.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 11839051746.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13 + 4574037537.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 436284060.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 11 - 579547620.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 282604824.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 &
     & - 11364372.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 21380976.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + &
     & 1779840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 933120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 82944.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 84.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**28 - 2160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 &
     & + 24020.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 137000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 224441.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**24 + 2882704.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 29514956.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 22 + 158761380.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 600525385.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + &
     & 1735046508.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 3974317128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + &
     & 7366109968.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 11182879285.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + &
     & 14008821640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 14525046988.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 12424411844.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 8628406707.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
     & 4656691420.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 1737867780.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
     & 271070256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 121027908.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 79902000.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 10805040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 3600000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **5 - 767232.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 110592.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 48.d0*gl%z3_PCSAFT(1) &
     & **28 + 1440.d0*gl%z3_PCSAFT(1)**27 - 20270.d0*gl%z3_PCSAFT(1)**26 + 177380.d0*gl%z3_PCSAFT(1)**25 - &
     & 1076171.d0*gl%z3_PCSAFT(1)**24 + 4766248.d0*gl%z3_PCSAFT(1)**23 - 15705399.d0*gl%z3_PCSAFT(1)**22 + &
     & 37902170.d0*gl%z3_PCSAFT(1)**21 - 60650800.d0*gl%z3_PCSAFT(1)**20 + 30583336.d0*gl%z3_PCSAFT(1)**19 + &
     & 168159120.d0*gl%z3_PCSAFT(1)**18 - 677271824.d0*gl%z3_PCSAFT(1)**17 + 1543978825.d0*gl%z3_PCSAFT(1)** &
     & 16 - 2584387840.d0*gl%z3_PCSAFT(1)**15 + 3385529885.d0*gl%z3_PCSAFT(1)**14 - 3535630038.d0* &
     & gl%z3_PCSAFT(1)**13 + 2939046826.d0*gl%z3_PCSAFT(1)**12 - 1913683280.d0*gl%z3_PCSAFT(1)**11 + &
     & 943912040.d0*gl%z3_PCSAFT(1)**10 - 329557752.d0*gl%z3_PCSAFT(1)**9 + 67410232.d0*gl%z3_PCSAFT(1)**8 + &
     & 170496.d0*gl%z3_PCSAFT(1)**7 - 5162160.d0*gl%z3_PCSAFT(1)**6 + 1776000.d0*gl%z3_PCSAFT(1)**5 - 263808.d0 &
     & *gl%z3_PCSAFT(1)**4 - 13824.d0*gl%z3_PCSAFT(1)**3 + 9216.d0*gl%z3_PCSAFT(1)**2)
        part3 = (24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27 - 1680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + 39140.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**25 - 504200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 4299782.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**23 - 26591508.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 125974859.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**21 - 473989210.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 1451755955.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**19 - 3677984500.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 7772982770.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**17 - 13721514924.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 20120913700.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**15 - 24205886460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 23406178307.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 17611040946.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + &
     & 9748346703.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 3517359420.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + &
     & 516452940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 163441872.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - &
     & 82448532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 5537376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 889920.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 466560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 41472.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**3 + 36.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27 + 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 - &
     & 27540.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 464700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 4442793.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 29011056.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 140810061.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 533937690.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 1635368130.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 4141738416.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 8811917442.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 15892338288.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - &
     & 24360661905.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 31646200440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 &
     & - 34580504781.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 31413017826.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 12 - 23312412552.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 13739660160.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**10 - 6092100960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 1779032472.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8 - 178955808.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 91465920.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6 + 34912080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 567360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **4 - 981504.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 54.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 + 900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 1545.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - &
     & 94200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 + 1321821.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 9964776.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**22 + 51850752.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 204096690.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **20 + 639504840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 1650721836.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + &
     & 3597346575.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 6727619628.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + &
     & 10884531765.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 15235004760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 18324077130.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 18713161554.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
     & 15964586892.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 11141421960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
     & 6184276320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 2619809928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + &
     & 788292096.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 140928288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 2531040.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 5851200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 1449792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **3 + 85248.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 18432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 15.d0*gl%z3_PCSAFT(1)**27 &
     & - 330.d0*gl%z3_PCSAFT(1)**26 + 2765.d0*gl%z3_PCSAFT(1)**25 - 6050.d0*gl%z3_PCSAFT(1)**24 - 87604.d0* &
     & gl%z3_PCSAFT(1)**23 + 1051470.d0*gl%z3_PCSAFT(1)**22 - 6443878.d0*gl%z3_PCSAFT(1)**21 + 27610760.d0* &
     & gl%z3_PCSAFT(1)**20 - 91268065.d0*gl%z3_PCSAFT(1)**19 + 245112266.d0*gl%z3_PCSAFT(1)**18 - &
     & 552973423.d0*gl%z3_PCSAFT(1)**17 + 1070948374.d0*gl%z3_PCSAFT(1)**16 - 1801278330.d0*gl%z3_PCSAFT(1) &
     & **15 + 2637503410.d0*gl%z3_PCSAFT(1)**14 - 3346585912.d0*gl%z3_PCSAFT(1)**13 + 3646896388.d0 &
     & *gl%z3_PCSAFT(1)**12 - 3375154496.d0*gl%z3_PCSAFT(1)**11 + 2619334480.d0*gl%z3_PCSAFT(1)**10 - &
     & 1679372000.d0*gl%z3_PCSAFT(1)**9 + 872648272.d0*gl%z3_PCSAFT(1)**8 - 357473952.d0*gl%z3_PCSAFT(1)**7 &
     & + 110280032.d0*gl%z3_PCSAFT(1)**6 - 23366400.d0*gl%z3_PCSAFT(1)**5 + 2558080.d0*gl%z3_PCSAFT(1)**4 + &
     & 140928.d0*gl%z3_PCSAFT(1)**3 - 85248.d0*gl%z3_PCSAFT(1)**2 + 8448.d0*gl%z3_PCSAFT(1))
	part4 = (24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27 - 1680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + &
     & 39140.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 - 504200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 4299782.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - 26591508.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 125974859.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - 473989210.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 1451755955.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 3677984500.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 7772982770.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 13721514924.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + &
     & 20120913700.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 24205886460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 &
     & + 23406178307.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 17611040946.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 12 + 9748346703.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 3517359420.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **10 + 516452940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 163441872.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 8 - 82448532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 5537376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + &
     & 889920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 466560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 41472.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 36.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27 + 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **26 - 27540.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 464700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - &
     & 4442793.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 29011056.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - &
     & 140810061.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 533937690.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - &
     & 1635368130.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 4141738416.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - &
     & 8811917442.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 15892338288.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 &
     & - 24360661905.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 31646200440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 14 - 34580504781.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 31413017826.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**12 - 23312412552.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 13739660160.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10 - 6092100960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 1779032472.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**8 - 178955808.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 91465920.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6 + 34912080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 567360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **4 - 981504.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 54.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 + 900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 1545.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - &
     & 94200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 + 1321821.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 9964776.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**22 + 51850752.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 204096690.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **20 + 639504840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 1650721836.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + &
     & 3597346575.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 6727619628.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + &
     & 10884531765.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 15235004760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 18324077130.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 18713161554.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
     & 15964586892.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 11141421960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
     & 6184276320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 2619809928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + &
     & 788292096.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 140928288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 2531040.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 5851200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 1449792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **3 + 85248.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 18432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 15.d0*gl%z3_PCSAFT(1)**27 &
     & - 330.d0*gl%z3_PCSAFT(1)**26 + 2765.d0*gl%z3_PCSAFT(1)**25 - 6050.d0*gl%z3_PCSAFT(1)**24 - 87604.d0* &
     & gl%z3_PCSAFT(1)**23 + 1051470.d0*gl%z3_PCSAFT(1)**22 - 6443878.d0*gl%z3_PCSAFT(1)**21 + 27610760.d0* &
     & gl%z3_PCSAFT(1)**20 - 91268065.d0*gl%z3_PCSAFT(1)**19 + 245112266.d0*gl%z3_PCSAFT(1)**18 - &
     & 552973423.d0*gl%z3_PCSAFT(1)**17 + 1070948374.d0*gl%z3_PCSAFT(1)**16 - 1801278330.d0*gl%z3_PCSAFT(1) &
     & **15 + 2637503410.d0*gl%z3_PCSAFT(1)**14 - 3346585912.d0*gl%z3_PCSAFT(1)**13 + 3646896388.d0 &
     & *gl%z3_PCSAFT(1)**12 - 3375154496.d0*gl%z3_PCSAFT(1)**11 + 2619334480.d0*gl%z3_PCSAFT(1)**10 - &
     & 1679372000.d0*gl%z3_PCSAFT(1)**9 + 872648272.d0*gl%z3_PCSAFT(1)**8 - 357473952.d0*gl%z3_PCSAFT(1)**7 &
     & + 110280032.d0*gl%z3_PCSAFT(1)**6 - 23366400.d0*gl%z3_PCSAFT(1)**5 + 2558080.d0*gl%z3_PCSAFT(1)**4 + &
     & 140928.d0*gl%z3_PCSAFT(1)**3 - 85248.d0*gl%z3_PCSAFT(1)**2 + 8448.d0*gl%z3_PCSAFT(1))
	part5 = (72.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26 - 5760.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25 + &
     & 130060.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 - 1568200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 + &
     & 12424714.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 - 72025224.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 + &
     & 326859181.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 - 1213612730.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 + &
     & 3781454365.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 - 9985308500.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 + &
     & 22326192910.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 - 42000095556.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 &
     & + 65905831100.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 - 85399406040.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 13 + 90233394229.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 - 76404304578.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**11 + 50540662497.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 - 25122947100.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**9 + 8825640660.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 - 1982688696.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**7 + 247965444.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 - 21380976.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**5 + 1779840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 + 933120.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 3 + 82944.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 - 36.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + 8640.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 - 223100.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 2820200.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**23 - 22775273.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 132844844.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**21 - 602923312.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 2236337760.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**19 - 6983565215.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 18605547996.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**17 - 42351329448.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 81982127096.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 134064069635.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + &
     & 183919955540.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 210099189056.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 12 + 198022078144.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 152066752389.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**10 + 93414878900.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 - 44626428540.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**8 + 15813916992.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - 3791660652.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**6 + 482580720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 467280.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**4 - 2309760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 1292544.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **2 - 110592.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) - 54.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 - 4320.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 140055.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 1891050.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**23 + 15690804.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 92558676.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**21 + 421749441.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 1567515930.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**19 + 4917827100.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 13245585384.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**17 + 30741479685.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 61265647614.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 104238162960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - &
     & 150534303540.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 183484158195.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 12 - 187558015014.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 159418233006.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**10 - 111294207840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 62674507200.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**8 - 27712355112.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 9229156296.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**6 - 2161272384.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 311645520.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**4 - 18673920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 813888.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**2 + 85248.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 9216.d0*gl%mmean_PCSAFT**2 + 45.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**26 + 720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - 37760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 + &
     & 560000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 4821497.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 28932480.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 132925388.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 496612360.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19 - 1568481365.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 4274726368.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**17 - 10113305576.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 20715853664.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**15 - 36534742455.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 55166107280.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13 - 70961955116.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 77365368056.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11 - 71040780376.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 54474739520.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9 - 34467972880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 17696182496.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**7 - 7198397664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 2240276032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **5 - 505055040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 74935040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - &
     & 5914176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 86784.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 8448.d0*gl%mmean_PCSAFT - 9.d0* &
     & gl%z3_PCSAFT(1)**26 + 3625.d0*gl%z3_PCSAFT(1)**24 - 61750.d0*gl%z3_PCSAFT(1)**23 + 558412.d0*gl%z3_PCSAFT(1)** &
     & 22 - 3430600.d0*gl%z3_PCSAFT(1)**21 + 15962054.d0*gl%z3_PCSAFT(1)**20 - 60118740.d0*gl%z3_PCSAFT(1)** &
     & 19 + 191389975.d0*gl%z3_PCSAFT(1)**18 - 527587752.d0*gl%z3_PCSAFT(1)**17 + 1269857253.d0* &
     & gl%z3_PCSAFT(1)**16 - 2664439958.d0*gl%z3_PCSAFT(1)**15 + 4847183030.d0*gl%z3_PCSAFT(1)**14 - &
     & 7603003520.d0*gl%z3_PCSAFT(1)**13 + 10234554284.d0*gl%z3_PCSAFT(1)**12 - 11773819024.d0* &
     & gl%z3_PCSAFT(1)**11 + 11521274256.d0*gl%z3_PCSAFT(1)**10 - 9531629600.d0*gl%z3_PCSAFT(1)**9 + &
     & 6610913760.d0*gl%z3_PCSAFT(1)**8 - 3799462080.d0*gl%z3_PCSAFT(1)**7 + 1780567744.d0*gl%z3_PCSAFT(1)** &
     & 6 - 665317504.d0*gl%z3_PCSAFT(1)**5 + 191949760.d0*gl%z3_PCSAFT(1)**4 - 40743040.d0*gl%z3_PCSAFT(1)** &
     & 3 + 5879936.d0*gl%z3_PCSAFT(1)**2 - 497664.d0*gl%z3_PCSAFT(1) + 17152.d0)
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
            	    gl%cx2_PCSAFT(3,xi,xj) = part1*(gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*part2 + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(1,xi)*part3 + gl%mPCSAFT(xi)* &
         & gl%z3x1_PCSAFT(1,xj)*part4 + gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part5)
                end if
            end do
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF C_1 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires z3_PCSAFT(1), z3_PCSAFT(4)
    if (GETDERC(4) .eq. 1) then
        part1 = -2.d0*(gl%z3_PCSAFT(1) - 1.d0)/recurring_factor**4
        part2 = (8.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 368.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 6552.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 66320.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 442018.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - &
     & 2085220.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 7236292.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4) - 18829460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 36923072.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 54098356.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + &
     & 57660540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 42021172.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) + 17650566.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 1056408.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 2766360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 974520.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 106128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - &
     & 69120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 6912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + &
     & 32.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 752.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 8220.d0*gl%z3_PCSAFT(1)** &
     & 18*gl%z3_PCSAFT(4) - 55072.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 250734.d0*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) - 810928.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 1878768.d0*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4) - 3004128.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 2771940.d0*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4) + 338536.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 5894116.d0*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) + 10693032.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 11311682.d0*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4) + 7678296.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 3045992.d0*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4) + 334008.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 276048.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) &
     & - 112704.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 1152.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 4608.d0* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4))
        part3 = (-4.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 184.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 3276.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 33160.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17 - 221009.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 1042610.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **15 - 3618146.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 9414730.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 &
     & - 18461536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 27049178.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - &
     & 28830270.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 21010586.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - &
     & 8825283.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 528204.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 1383180.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 487260.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 53064.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**4 + 34560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 3456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - &
     & 8.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 552.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 &
     & - 3920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 66462.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 448664.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**15 + 1923528.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 5830296.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 &
     & + 13008588.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 21647344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + &
     & 26719192.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 23815056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 14331758.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 4754232.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 66712.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 &
     & + 733176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 213744.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 13824.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**3 + 11520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 5.d0*gl%z3_PCSAFT(1)**20 - 110.d0*gl%z3_PCSAFT(1)**19 &
     & + 1104.d0*gl%z3_PCSAFT(1)**18 - 6534.d0*gl%z3_PCSAFT(1)**17 + 24200.d0*gl%z3_PCSAFT(1)**16 - 50794.d0* &
     & gl%z3_PCSAFT(1)**15 + 9238.d0*gl%z3_PCSAFT(1)**14 + 339894.d0*gl%z3_PCSAFT(1)**13 - 1342133.d0*gl%z3_PCSAFT(1) &
     & **12 + 3065896.d0*gl%z3_PCSAFT(1)**11 - 4858318.d0*gl%z3_PCSAFT(1)**10 + 5572200.d0*gl%z3_PCSAFT(1)** &
     & 9 - 4605528.d0*gl%z3_PCSAFT(1)**8 + 2609944.d0*gl%z3_PCSAFT(1)**7 - 854504.d0*gl%z3_PCSAFT(1)**6 + &
     & 17696.d0*gl%z3_PCSAFT(1)**5 + 120912.d0*gl%z3_PCSAFT(1)**4 - 49120.d0*gl%z3_PCSAFT(1)**3 + 4608.d0* &
     & gl%z3_PCSAFT(1)**2 + 1728.d0*gl%z3_PCSAFT(1) - 384.d0)
	part4 = (12.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 612.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 10880.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 106712.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 681631.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 3075481.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 10251739.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 25885565.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 50074577.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) - 74150471.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 82775761.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 67043103.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4) + 36033036.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 9665124.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 1146708.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 1439964.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 117792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) - 69120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 6912.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 24.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **18*gl%z3_PCSAFT(4) - 776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 33272.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 305506.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 1642014.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 6079314.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 16549086.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 34220286.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4) + 54634434.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 67724910.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 64948994.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 47505656.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 25653736.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - &
     & 9448152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 1788168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4) + 199872.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 192384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **2*gl%z3_PCSAFT(4) + 23040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4608.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) &
     & - 15.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 285.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 2305.d0*gl%z3_PCSAFT(1) &
     & **17*gl%z3_PCSAFT(4) + 9475.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 12047.d0*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) - 86543.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 621905.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) &
     & - 2244223.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 5570442.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - &
     & 10346082.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 14876292.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - &
     & 16772000.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 14806880.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - &
     & 10083616.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 5123072.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 1813904.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 377440.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 15072.d0* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 12096.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2112.d0*gl%z3_PCSAFT(4))
        part5 = (-4.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 184.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - &
     & 3276.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 33160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 221009.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 1042610.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 3618146.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**14 + 9414730.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 18461536.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**12 + 27049178.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 28830270.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**10 + 21010586.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 8825283.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8 + 528204.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 1383180.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 6 - 487260.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 53064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 34560.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 3456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 8.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 &
     & + 128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 552.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 3920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **17 + 66462.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 448664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 1923528.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 5830296.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 13008588.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12 - 21647344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 26719192.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 10 - 23815056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 14331758.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - &
     & 4754232.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 66712.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 733176.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5 - 213744.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 13824.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + &
     & 11520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 5.d0*gl%z3_PCSAFT(1)**20 - 110.d0*gl%z3_PCSAFT(1)**19 + 1104.d0* &
     & gl%z3_PCSAFT(1)**18 - 6534.d0*gl%z3_PCSAFT(1)**17 + 24200.d0*gl%z3_PCSAFT(1)**16 - 50794.d0*gl%z3_PCSAFT(1)** &
     & 15 + 9238.d0*gl%z3_PCSAFT(1)**14 + 339894.d0*gl%z3_PCSAFT(1)**13 - 1342133.d0*gl%z3_PCSAFT(1)**12 + &
     & 3065896.d0*gl%z3_PCSAFT(1)**11 - 4858318.d0*gl%z3_PCSAFT(1)**10 + 5572200.d0*gl%z3_PCSAFT(1)**9 - &
     & 4605528.d0*gl%z3_PCSAFT(1)**8 + 2609944.d0*gl%z3_PCSAFT(1)**7 - 854504.d0*gl%z3_PCSAFT(1)**6 + 17696.d0* &
     & gl%z3_PCSAFT(1)**5 + 120912.d0*gl%z3_PCSAFT(1)**4 - 49120.d0*gl%z3_PCSAFT(1)**3 + 4608.d0*gl%z3_PCSAFT(1)**2 + &
     & 1728.d0*gl%z3_PCSAFT(1) - 384.d0)
        part6 = (12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 612.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 10880.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 106712.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + &
     & 681631.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 3075481.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 &
     & *gl%z3_PCSAFT(4) + 10251739.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 25885565.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 50074577.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4) - 74150471.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 82775761.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 67043103.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + &
     & 36033036.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 9665124.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 &
     & *gl%z3_PCSAFT(4) - 1146708.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 1439964.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 117792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 69120.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 6912.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 24.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 776.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 33272.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - &
     & 305506.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 1642014.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4) - 6079314.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 16549086.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 34220286.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + &
     & 54634434.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 67724910.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) + 64948994.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 47505656.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 25653736.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 9448152.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 1788168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 199872.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 192384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & + 23040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4608.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) - 15.d0*gl%z3_PCSAFT(1)** &
     & 19*gl%z3_PCSAFT(4) + 285.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 2305.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) &
     & + 9475.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 12047.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 86543.d0* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 621905.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 2244223.d0*gl%z3_PCSAFT(1) &
     & **12*gl%z3_PCSAFT(4) + 5570442.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 10346082.d0*gl%z3_PCSAFT(1)**10 &
     & *gl%z3_PCSAFT(4) + 14876292.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 16772000.d0*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4) + 14806880.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 10083616.d0*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4) + 5123072.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 1813904.d0*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4) + 377440.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 15072.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - &
     & 12096.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2112.d0*gl%z3_PCSAFT(4))
	part7 = ( &
     & -12.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + 612.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 10880.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**17 + 106712.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 681631.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**15 + 3075481.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 10251739.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**13 + 25885565.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 50074577.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11 + 74150471.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 82775761.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**9 + 67043103.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - 36033036.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7 + 9665124.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 1146708.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **5 - 1439964.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 117792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + &
     & 69120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 + 6912.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) - 480.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**18 + 11286.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 124206.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 16 + 847815.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 4010481.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + &
     & 13915743.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 36558321.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + &
     & 73935357.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 115694667.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + &
     & 139182879.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 125996061.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + &
     & 81836268.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 34045476.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + &
     & 5779788.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 1914060.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 1111296.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 91296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 24192.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1) + 2304.d0*gl%mmean_PCSAFT**2 + 9.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 21.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 &
     & - 3117.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 43983.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 334707.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**15 + 1692093.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 6176067.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 &
     & + 16987005.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 36036618.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + &
     & 59604378.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 76895964.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 76548216.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 57246624.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 30354720.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**6 - 9812928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 737232.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + &
     & 842400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 345504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 39360.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1) + 2112.d0*gl%mmean_PCSAFT - 3.d0*gl%z3_PCSAFT(1)**19 + 33.d0*gl%z3_PCSAFT(1)**18 + 97.d0*gl%z3_PCSAFT(1)** &
     & 17 - 4279.d0*gl%z3_PCSAFT(1)**16 + 40607.d0*gl%z3_PCSAFT(1)**15 - 227921.d0*gl%z3_PCSAFT(1)**14 + &
     & 891059.d0*gl%z3_PCSAFT(1)**13 - 2592133.d0*gl%z3_PCSAFT(1)**12 + 5802520.d0*gl%z3_PCSAFT(1)**11 - &
     & 10169020.d0*gl%z3_PCSAFT(1)**10 + 14039000.d0*gl%z3_PCSAFT(1)**9 - 15211672.d0*gl%z3_PCSAFT(1)**8 + &
     & 12737344.d0*gl%z3_PCSAFT(1)**7 - 7963712.d0*gl%z3_PCSAFT(1)**6 + 3440480.d0*gl%z3_PCSAFT(1)**5 - &
     & 801536.d0*gl%z3_PCSAFT(1)**4 - 67520.d0*gl%z3_PCSAFT(1)**3 + 119552.d0*gl%z3_PCSAFT(1)**2 - 37184.d0* &
     & gl%z3_PCSAFT(1) + 4288.d0)
	part8 = (-12.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 &
     & + 612.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 10880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + 106712.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 681631.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + 3075481.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**14 - 10251739.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + 25885565.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**12 - 50074577.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 + 74150471.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**10 - 82775761.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 67043103.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**8 - 36033036.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 9665124.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **6 + 1146708.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 1439964.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + &
     & 117792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 69120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 + 6912.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1) - 480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 11286.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 &
     & - 124206.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 847815.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - &
     & 4010481.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 13915743.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - &
     & 36558321.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 73935357.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - &
     & 115694667.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 139182879.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - &
     & 125996061.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 81836268.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - &
     & 34045476.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 5779788.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + &
     & 1914060.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 1111296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 91296.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 24192.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 2304.d0*gl%mmean_PCSAFT**2 + 9.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 21.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 3117.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + &
     & 43983.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 334707.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 1692093.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14 - 6176067.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 16987005.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 &
     & - 36036618.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 59604378.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - &
     & 76895964.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 76548216.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 57246624.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 30354720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 9812928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **5 + 737232.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 842400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 345504.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 39360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 2112.d0*gl%mmean_PCSAFT - 3.d0*gl%z3_PCSAFT(1)**19 &
     & + 33.d0*gl%z3_PCSAFT(1)**18 + 97.d0*gl%z3_PCSAFT(1)**17 - 4279.d0*gl%z3_PCSAFT(1)**16 + 40607.d0*gl%z3_PCSAFT(1) &
     & **15 - 227921.d0*gl%z3_PCSAFT(1)**14 + 891059.d0*gl%z3_PCSAFT(1)**13 - 2592133.d0*gl%z3_PCSAFT(1)**12 &
     & + 5802520.d0*gl%z3_PCSAFT(1)**11 - 10169020.d0*gl%z3_PCSAFT(1)**10 + 14039000.d0*gl%z3_PCSAFT(1)**9 - &
     & 15211672.d0*gl%z3_PCSAFT(1)**8 + 12737344.d0*gl%z3_PCSAFT(1)**7 - 7963712.d0*gl%z3_PCSAFT(1)**6 + &
     & 3440480.d0*gl%z3_PCSAFT(1)**5 - 801536.d0*gl%z3_PCSAFT(1)**4 - 67520.d0*gl%z3_PCSAFT(1)**3 + 119552.d0* &
     & gl%z3_PCSAFT(1)**2 - 37184.d0*gl%z3_PCSAFT(1) + 4288.d0)
	part9 = (48.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 2688.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) &
     & + 46920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 440688.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4) + 2685072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 11639472.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 37820544.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - &
     & 94869936.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 186177984.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **10*gl%z3_PCSAFT(4) - 285790752.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 338571768.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 300527040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) + 188963616.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 74766816.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 12735072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 1615248.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 353376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) - 207360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 20736.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(4) + 2400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 52764.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 540744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - &
     & 3450312.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 15426216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(4) - 51451488.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 132643272.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 269219784.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) + 432649560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 547836564.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 538435344.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) - 399847440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 214506480.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 76749456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4) + 15374064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 810720.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 156096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 3456.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) - 36.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 384.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 17628.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 211272.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 1446000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - &
     & 6750024.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 23293104.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4) - 62075208.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 130835460.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 220443816.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 297205764.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 318374160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) + 266980560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 171031920.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 80498256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 26049984.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 5082240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - &
     & 408000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 14208.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) + 12.d0*gl%z3_PCSAFT(1)** &
     & 18*gl%z3_PCSAFT(4) - 72.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 1488.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + &
     & 25536.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 196056.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 973392.d0* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 3511248.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 9735552.d0* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 21385524.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 37787352.d0* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 53984160.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 62208192.d0* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 57318240.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 41585760.d0* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 23198496.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 9589056.d0*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4) - 2763072.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 494976.d0*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) - 41472.d0*gl%z3_PCSAFT(4))
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
            	    gl%cx2_PCSAFT(4,xi,xj) = part1*(gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*part2 + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(4,xi)*part3 + &
			    & gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(1,xi)*part4 + gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(4,xj)*part5 + & 
			    & gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(1,xj)*part6 + gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*part7 + &
			    & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*part8 + gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part9)
                end if
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF C_1 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2 
    ! requires z3_PCSAFT(1), z3_PCSAFT(4), z3_PCSAFT(5)
    if (GETDERC(5) .eq. 1) then
        part1 = -2.d0/recurring_factor**5
        part2 = (16.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(5) - 48.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **26*gl%z3_PCSAFT(4)**2 - 880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 2880.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + 20264.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) &
     & - 67160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 276200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **24*gl%z3_PCSAFT(5) + 903440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 2562084.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 8176508.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4)**2 - 17398812.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) + 53892096.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 90287438.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(5) - 271015922.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 367632998.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 1070687860.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4)**2 + 1193711530.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - &
     & 3385354730.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 3118077930.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 8663404360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 &
     & + 6569196380.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 18039562700.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - 11124192988.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(5) + 30557850792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + &
     & 14977983864.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 41839042960.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 15687599960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(5) + 45617534880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + &
     & 12242103214.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 38477437298.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 6446535814.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) &
     & + 23678103492.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 1575646146.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 9148075074.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)**2 + 534997854.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 872568120.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 573873600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(5) + 1159095240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 131352624.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 565209648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4)**2 + 36773352.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 22728744.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 19154448.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(5) + 42761952.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 907200.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 3559680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4)**2 + 933120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 1866240.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 82944.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - &
     & 165888.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 56.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27* &
     & gl%z3_PCSAFT(5) - 168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)**2 - 1640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **26*gl%z3_PCSAFT(5) + 4320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + 21976.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 48040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - &
     & 173640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 274000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(4)**2 + 843446.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 448882.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 - 2031034.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - &
     & 5765408.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 3879688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(5) + 59029912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 62163308.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 317522760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) &
     & **2 - 326625910.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 1201050770.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 1138151678.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - &
     & 3470093016.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 2957603032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **17*gl%z3_PCSAFT(5) + 7948634256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + &
     & 5965444904.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 14732219936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **15*gl%z3_PCSAFT(4)**2 - 9469367070.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + &
     & 22365758570.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + 11819551330.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 28017643280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - &
     & 11434547288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 29050093976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **12*gl%z3_PCSAFT(4)**2 + 8281748236.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - &
     & 24848823688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 4138672946.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + 17256813414.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + &
     & 1079963698.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 9313382840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 9*gl%z3_PCSAFT(4)**2 + 167460800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 3475735560.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 258356184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) &
     & - 542140512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + 78703512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 7*gl%z3_PCSAFT(5) - 242055816.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 3034320.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 159804000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 &
     & - 6868512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 21610080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4)**2 + 901440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 7200000.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 165888.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + &
     & 1534464.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 27648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(5) + 221184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 32.d0*gl%z3_PCSAFT(1)**27* &
     & gl%z3_PCSAFT(5) + 96.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)**2 + 1040.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) &
     & - 2880.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - 16108.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + &
     & 40540.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + 157580.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - &
     & 354760.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 1088146.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + &
     & 2152342.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 5614202.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - &
     & 9532496.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 22331162.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + &
     & 31410798.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 69551646.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - &
     & 75804340.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 169965900.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) &
     & + 121301600.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 320797416.d0*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(5) - 61166672.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 442978920.d0*gl%z3_PCSAFT(1)** &
     & 17*gl%z3_PCSAFT(5) - 336318240.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + 364557048.d0* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 1354543648.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + &
     & 81361550.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 3087957650.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 &
     & - 860915150.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 5168775680.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) &
     & **2 + 1638584830.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 6771059770.d0*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)**2 - 1967101162.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + 7071260076.d0*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4)**2 + 1661424768.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 5878093652.d0* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 969145860.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + &
     & 3827366560.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 331113680.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) &
     & - 1887824080.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 2950712.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) &
     & + 659115504.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 60759632.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) &
     & - 134820464.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 28279328.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) &
     & - 340992.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 2922528.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + &
     & 10324320.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 1838400.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - &
     & 3552000.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 575232.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + &
     & 527616.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 13824.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 27648.d0* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 18432.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 18432.d0*gl%z3_PCSAFT(4)**2)
        part3 = (-8.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27 + 440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **26 - 10132.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + 138100.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - &
     & 1281042.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + 8699406.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - &
     & 45143719.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + 183816499.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - &
     & 596855765.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + 1559038965.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - &
     & 3284598190.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + 5562096494.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - &
     & 7488991932.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + 7843799980.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - &
     & 6121051607.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + 3223267907.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - &
     & 787823073.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 267498927.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + &
     & 286936800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 - 65676312.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - &
     & 18386676.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 9577224.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 453600.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 466560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 41472.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**3 - 12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27 + 180.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + &
     & 1128.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - 57720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 797313.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - 6707907.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 40047441.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - 180818301.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 638285700.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 1794164886.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 4052814534.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 7372263618.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + &
     & 10747951845.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 12392896155.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 &
     & + 11001878961.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 7111849077.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 12 + 2890522962.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 283371936.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 10 - 441024480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 255569388.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 &
     & - 31862544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 19299840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + &
     & 6290064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 433440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 262656.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 13824.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 18.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 27 - 510.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 + 6507.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - 46995.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**24 + 180729.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 17727.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - &
     & 5063862.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 37216986.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 167991780.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 554465796.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 1415905185.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**17 + 2863796673.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 4625576175.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**15 + 5952877095.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 6024586500.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13 + 4653593928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 2564330712.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11 + 819289380.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 31531440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 9 - 176881392.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 78927648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - &
     & 7848432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 5432448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 1833600.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**4 - 17088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 61056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 4608.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 5.d0*gl%z3_PCSAFT(1)**27 + 155.d0*gl%z3_PCSAFT(1)**26 - 2269.d0*gl%z3_PCSAFT(1)**25 + &
     & 20665.d0*gl%z3_PCSAFT(1)**24 - 129576.d0*gl%z3_PCSAFT(1)**23 + 581070.d0*gl%z3_PCSAFT(1)**22 - &
     & 1842142.d0*gl%z3_PCSAFT(1)**21 + 3663394.d0*gl%z3_PCSAFT(1)**20 - 1237565.d0*gl%z3_PCSAFT(1)**19 - &
     & 22351321.d0*gl%z3_PCSAFT(1)**18 + 101052831.d0*gl%z3_PCSAFT(1)**17 - 272383179.d0*gl%z3_PCSAFT(1)**16 &
     & + 536685698.d0*gl%z3_PCSAFT(1)**15 - 816364400.d0*gl%z3_PCSAFT(1)**14 + 973370092.d0*gl%z3_PCSAFT(1) &
     & **13 - 904431456.d0*gl%z3_PCSAFT(1)**12 + 636388376.d0*gl%z3_PCSAFT(1)**11 - 313940656.d0* &
     & gl%z3_PCSAFT(1)**10 + 81796320.d0*gl%z3_PCSAFT(1)**9 + 15543328.d0*gl%z3_PCSAFT(1)**8 - 25154336.d0* &
     & gl%z3_PCSAFT(1)**7 + 10226688.d0*gl%z3_PCSAFT(1)**6 - 1044928.d0*gl%z3_PCSAFT(1)**5 - 704320.d0* &
     & gl%z3_PCSAFT(1)**4 + 285952.d0*gl%z3_PCSAFT(1)**3 - 21504.d0*gl%z3_PCSAFT(1)**2 - 8448.d0*gl%z3_PCSAFT(1) + &
     & 1536.d0)
	part4 = (48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 2880.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 67160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) &
     & - 903440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 8176508.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 22*gl%z3_PCSAFT(4) - 53892096.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 271015922.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 1070687860.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 3385354730.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 8663404360.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 18039562700.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) - 30557850792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 41839042960.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 45617534880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) + 38477437298.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 23678103492.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 9148075074.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) - 872568120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 1159095240.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 565209648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) - 22728744.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 42761952.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 3559680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 1866240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 165888.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) + 72.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **25*gl%z3_PCSAFT(4) - 24000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 573600.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 6587022.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + &
     & 49636032.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 272376558.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **20*gl%z3_PCSAFT(4) + 1146223140.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - &
     & 3811177740.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 10195299264.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 22174994724.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + &
     & 39402054744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 57169646190.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 67317163440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 63452706654.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 46655879412.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 25435687296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + &
     & 9089919600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 1177648560.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 692890512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 411354624.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 68178240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(4) - 10787040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 4849920.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 193536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & - 55296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 108.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) &
     & + 2640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 26370.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4) + 103980.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 454494.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **22*gl%z3_PCSAFT(4) - 9211392.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 69024396.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 342340680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + &
     & 1266535320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 3671345664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 17*gl%z3_PCSAFT(4) + 8552565030.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 16235991924.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 25290390750.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) &
     & - 32368871760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 33903584880.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 28747012368.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + &
     & 19304816376.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 9817239840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 9*gl%z3_PCSAFT(4) + 3379628160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 464263872.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 236454528.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + &
     & 161801856.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 36672000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4) - 1432320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 2189952.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 198144.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 36864.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(4) + 30.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 840.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + &
     & 10790.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 82460.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 397364.d0* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 1064640.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 322564.d0*gl%z3_PCSAFT(1) &
     & **20*gl%z3_PCSAFT(4) + 18672200.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 104902810.d0*gl%z3_PCSAFT(1)** &
     & 18*gl%z3_PCSAFT(4) + 375116264.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 1004867266.d0*gl%z3_PCSAFT(1)** &
     & 16*gl%z3_PCSAFT(4) + 2131221412.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 3669221640.d0*gl%z3_PCSAFT(1) &
     & **14*gl%z3_PCSAFT(4) + 5187460480.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 6041423488.d0* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 5772272656.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 4467659888.d0 &
     & *gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 2728924480.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 1245441440.d0 &
     & *gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 366886528.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 25196544.d0* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 35636224.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 18683520.d0* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 4056320.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 133632.d0*gl%z3_PCSAFT(1)** &
     & 2*gl%z3_PCSAFT(4) + 113664.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 16896.d0*gl%z3_PCSAFT(4))
	part5 = (24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 96.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **25*gl%z3_PCSAFT(4)**2 - 1440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 6240.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + 33580.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(5) - 145440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 451720.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 1911840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)** &
     & 2 + 4088254.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 16776072.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 26946048.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + &
     & 107075112.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 135507961.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 522965640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 &
     & - 535343930.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 2018666280.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 1692677365.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) &
     & - 6288866640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 4331702180.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 16019373360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
     & **2 + 9019781350.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 33585528240.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 15278925396.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(5) + 58000880640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + &
     & 20919521480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 82080870360.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 22808767440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(5) + 94029307800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + &
     & 19238718649.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 85289793912.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 11839051746.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) + 58900185384.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + &
     & 4574037537.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 28644768480.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 436284060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + &
     & 7907286960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 579547620.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 126189360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + &
     & 282604824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 892093392.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 11364372.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + &
     & 187625808.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 21380976.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 31687200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + &
     & 1779840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 5339520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 &
     & *gl%z3_PCSAFT(4)**2 + 933120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 2799360.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 82944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(5) - 248832.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 36.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - 240.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) &
     & **2 - 12000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 79080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4)**2 + 286800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 1503000.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 - 3293511.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(5) + 15472608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 24818016.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - 107658144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4)**2 - 136188279.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 553996680.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + 573111570.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(5) - 2214098520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - &
     & 1905588870.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 7081914000.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 5097649632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) &
     & - 18478776096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - 11087497362.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 39798829608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4)**2 + 19701027372.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - &
     & 71186731320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 28584823095.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 105890970000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4)**2 + 33658581720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - &
     & 130609564320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 31726353327.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + 132613716216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4)**2 + 23327939706.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - &
     & 109481915064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 12717843648.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 72060512400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4)**2 + 4544959800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - &
     & 36569239920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 588824280.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 13361850480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 &
     & - 346445256.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 2865174432.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 205677312.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - &
     & 53443008.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 34089120.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 251110080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - &
     & 5393520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 59037120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(4)**2 + 2424960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 3715200.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(5) + 2156544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 27648.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 248832.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**2 - 54.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **26*gl%z3_PCSAFT(5) + 216.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + 1320.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 4440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 13185.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 29460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + &
     & 51990.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 84420.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4)**2 + 227247.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 3098136.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 4605696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + &
     & 29140944.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 34512198.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 &
     & *gl%z3_PCSAFT(5) - 172725900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 171170340.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 750534060.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) &
     & **2 + 633267660.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - 2545545000.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 1835672832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + &
     & 6972789336.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + 4276282515.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **16*gl%z3_PCSAFT(5) - 15747258180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - &
     & 8117995962.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 29691231180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **14*gl%z3_PCSAFT(4)**2 + 12645195375.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - &
     & 47059454280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 16184435880.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 62838881280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + &
     & 16951792440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 70551739140.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4)**2 - 14373506184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + &
     & 66173335476.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 9652408188.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 51233990160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - &
     & 4908619920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 32100083760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 8*gl%z3_PCSAFT(4)**2 + 1689814080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - &
     & 15748180800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 232131936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **7*gl%z3_PCSAFT(5) + 5703883728.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - &
     & 118227264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 1340129664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4)**2 + 80900928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 120054720.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 18336000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) &
     & + 31609920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 716160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(5) - 10270080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 1094976.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 709632.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 99072.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 27648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**2 - 18432.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(5) + 15.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 60.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - &
     & 420.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 1500.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + 5395.d0* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 16320.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 41230.d0* &
     & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 94560.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 198682.d0* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 222156.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 532320.d0* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - 1038300.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 161282.d0* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 13210320.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + 9336100.d0* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 73893720.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 52451405.d0 &
     & *gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 287438940.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + &
     & 187558132.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 865340796.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 &
     & - 502433633.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 2110814112.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) &
     & **2 + 1065610706.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 4273118160.d0*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)**2 - 1834610820.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 7271778300.d0*gl%z3_PCSAFT(1) &
     & **13*gl%z3_PCSAFT(4)**2 + 2593730240.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 10462467300.d0* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 3020711744.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + &
     & 12734595312.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 2886136328.d0*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) - 13066065432.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 2233829944.d0* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 11217968880.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + &
     & 1364462240.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 7967593440.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 &
     & - 622720720.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 4604185440.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)** &
     & 2 + 183443264.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 2112183072.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) &
     & **2 - 12598272.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 740144448.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) &
     & **2 - 17818112.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 184923840.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) &
     & **2 + 9341760.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 28049280.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)** &
     & 2 - 2028160.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 1059840.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + &
     & 66816.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 415488.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 56832.d0* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 56832.d0*gl%z3_PCSAFT(4)**2 - 8448.d0*gl%z3_PCSAFT(5))
	part6 = (-8.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27 + 440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - &
     & 10132.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + 138100.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 1281042.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + 8699406.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 45143719.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + 183816499.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 596855765.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + 1559038965.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 3284598190.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + 5562096494.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - &
     & 7488991932.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + 7843799980.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - &
     & 6121051607.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + 3223267907.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - &
     & 787823073.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 267498927.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + &
     & 286936800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 - 65676312.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - &
     & 18386676.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 9577224.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 453600.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 466560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 41472.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**3 - 12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27 + 180.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + &
     & 1128.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - 57720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 797313.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - 6707907.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 40047441.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - 180818301.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 638285700.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 1794164886.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 4052814534.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 7372263618.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + &
     & 10747951845.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 12392896155.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 &
     & + 11001878961.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 7111849077.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 12 + 2890522962.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 283371936.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 10 - 441024480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 255569388.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 &
     & - 31862544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 19299840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + &
     & 6290064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 433440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 262656.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 13824.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 18.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 27 - 510.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 + 6507.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - 46995.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**24 + 180729.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 17727.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - &
     & 5063862.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 37216986.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 167991780.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 554465796.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 1415905185.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**17 + 2863796673.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 4625576175.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**15 + 5952877095.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 6024586500.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13 + 4653593928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 2564330712.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11 + 819289380.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 31531440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 9 - 176881392.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 78927648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - &
     & 7848432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 5432448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 1833600.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**4 - 17088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 61056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 4608.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 5.d0*gl%z3_PCSAFT(1)**27 + 155.d0*gl%z3_PCSAFT(1)**26 - 2269.d0*gl%z3_PCSAFT(1)**25 + &
     & 20665.d0*gl%z3_PCSAFT(1)**24 - 129576.d0*gl%z3_PCSAFT(1)**23 + 581070.d0*gl%z3_PCSAFT(1)**22 - &
     & 1842142.d0*gl%z3_PCSAFT(1)**21 + 3663394.d0*gl%z3_PCSAFT(1)**20 - 1237565.d0*gl%z3_PCSAFT(1)**19 - &
     & 22351321.d0*gl%z3_PCSAFT(1)**18 + 101052831.d0*gl%z3_PCSAFT(1)**17 - 272383179.d0*gl%z3_PCSAFT(1)**16 &
     & + 536685698.d0*gl%z3_PCSAFT(1)**15 - 816364400.d0*gl%z3_PCSAFT(1)**14 + 973370092.d0*gl%z3_PCSAFT(1) &
     & **13 - 904431456.d0*gl%z3_PCSAFT(1)**12 + 636388376.d0*gl%z3_PCSAFT(1)**11 - 313940656.d0* &
     & gl%z3_PCSAFT(1)**10 + 81796320.d0*gl%z3_PCSAFT(1)**9 + 15543328.d0*gl%z3_PCSAFT(1)**8 - 25154336.d0* &
     & gl%z3_PCSAFT(1)**7 + 10226688.d0*gl%z3_PCSAFT(1)**6 - 1044928.d0*gl%z3_PCSAFT(1)**5 - 704320.d0* &
     & gl%z3_PCSAFT(1)**4 + 285952.d0*gl%z3_PCSAFT(1)**3 - 21504.d0*gl%z3_PCSAFT(1)**2 - 8448.d0*gl%z3_PCSAFT(1) + &
     & 1536.d0)
	part7 = (48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 2880.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 67160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) &
     & - 903440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 8176508.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 22*gl%z3_PCSAFT(4) - 53892096.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 271015922.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 1070687860.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 3385354730.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 8663404360.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 18039562700.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) - 30557850792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 41839042960.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 45617534880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) + 38477437298.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 23678103492.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 9148075074.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) - 872568120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 1159095240.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 565209648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) - 22728744.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 42761952.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 3559680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 1866240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 165888.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) + 72.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **25*gl%z3_PCSAFT(4) - 24000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 573600.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 6587022.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + &
     & 49636032.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 272376558.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **20*gl%z3_PCSAFT(4) + 1146223140.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - &
     & 3811177740.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 10195299264.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 22174994724.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + &
     & 39402054744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 57169646190.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 67317163440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 63452706654.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 46655879412.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 25435687296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + &
     & 9089919600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 1177648560.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 692890512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 411354624.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 68178240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(4) - 10787040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 4849920.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 193536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & - 55296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 108.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) &
     & + 2640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 26370.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4) + 103980.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 454494.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **22*gl%z3_PCSAFT(4) - 9211392.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 69024396.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 342340680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + &
     & 1266535320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 3671345664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 17*gl%z3_PCSAFT(4) + 8552565030.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 16235991924.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 25290390750.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) &
     & - 32368871760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 33903584880.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 28747012368.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + &
     & 19304816376.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 9817239840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 9*gl%z3_PCSAFT(4) + 3379628160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 464263872.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 236454528.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + &
     & 161801856.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 36672000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4) - 1432320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 2189952.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 198144.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 36864.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(4) + 30.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 840.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + &
     & 10790.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 82460.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 397364.d0* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 1064640.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 322564.d0*gl%z3_PCSAFT(1) &
     & **20*gl%z3_PCSAFT(4) + 18672200.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 104902810.d0*gl%z3_PCSAFT(1)** &
     & 18*gl%z3_PCSAFT(4) + 375116264.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 1004867266.d0*gl%z3_PCSAFT(1)** &
     & 16*gl%z3_PCSAFT(4) + 2131221412.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 3669221640.d0*gl%z3_PCSAFT(1) &
     & **14*gl%z3_PCSAFT(4) + 5187460480.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 6041423488.d0* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 5772272656.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 4467659888.d0 &
     & *gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 2728924480.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 1245441440.d0 &
     & *gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 366886528.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 25196544.d0* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 35636224.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 18683520.d0* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 4056320.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 133632.d0*gl%z3_PCSAFT(1)** &
     & 2*gl%z3_PCSAFT(4) + 113664.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 16896.d0*gl%z3_PCSAFT(4))
	part8 = (24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 96.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **25*gl%z3_PCSAFT(4)**2 - 1440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 6240.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + 33580.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(5) - 145440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 451720.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 1911840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)** &
     & 2 + 4088254.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 16776072.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 26946048.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + &
     & 107075112.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 135507961.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 522965640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 &
     & - 535343930.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 2018666280.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 1692677365.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) &
     & - 6288866640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 4331702180.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 16019373360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
     & **2 + 9019781350.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 33585528240.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 15278925396.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(5) + 58000880640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + &
     & 20919521480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 82080870360.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 22808767440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(5) + 94029307800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + &
     & 19238718649.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 85289793912.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 11839051746.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) + 58900185384.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + &
     & 4574037537.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 28644768480.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 436284060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + &
     & 7907286960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 579547620.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 126189360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + &
     & 282604824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 892093392.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 11364372.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + &
     & 187625808.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 21380976.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 31687200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + &
     & 1779840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 5339520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 &
     & *gl%z3_PCSAFT(4)**2 + 933120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 2799360.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 82944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(5) - 248832.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 36.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - 240.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) &
     & **2 - 12000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 79080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4)**2 + 286800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 1503000.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 - 3293511.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(5) + 15472608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 24818016.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - 107658144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4)**2 - 136188279.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 553996680.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + 573111570.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(5) - 2214098520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - &
     & 1905588870.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 7081914000.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 5097649632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) &
     & - 18478776096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - 11087497362.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 39798829608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4)**2 + 19701027372.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - &
     & 71186731320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 28584823095.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 105890970000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4)**2 + 33658581720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - &
     & 130609564320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 31726353327.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + 132613716216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4)**2 + 23327939706.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - &
     & 109481915064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 12717843648.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 72060512400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4)**2 + 4544959800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - &
     & 36569239920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 588824280.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 13361850480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 &
     & - 346445256.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 2865174432.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 205677312.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - &
     & 53443008.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 34089120.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 251110080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - &
     & 5393520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 59037120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(4)**2 + 2424960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 3715200.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(5) + 2156544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 27648.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 248832.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**2 - 54.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **26*gl%z3_PCSAFT(5) + 216.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + 1320.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 4440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 13185.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 29460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + &
     & 51990.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 84420.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4)**2 + 227247.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 3098136.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 4605696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + &
     & 29140944.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 34512198.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 &
     & *gl%z3_PCSAFT(5) - 172725900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 171170340.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 750534060.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) &
     & **2 + 633267660.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - 2545545000.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 1835672832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + &
     & 6972789336.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + 4276282515.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **16*gl%z3_PCSAFT(5) - 15747258180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - &
     & 8117995962.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 29691231180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **14*gl%z3_PCSAFT(4)**2 + 12645195375.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - &
     & 47059454280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 16184435880.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 62838881280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + &
     & 16951792440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 70551739140.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4)**2 - 14373506184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + &
     & 66173335476.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 9652408188.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 51233990160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - &
     & 4908619920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 32100083760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 8*gl%z3_PCSAFT(4)**2 + 1689814080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - &
     & 15748180800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 232131936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **7*gl%z3_PCSAFT(5) + 5703883728.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - &
     & 118227264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 1340129664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4)**2 + 80900928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 120054720.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 18336000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) &
     & + 31609920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 716160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(5) - 10270080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 1094976.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 709632.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 99072.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 27648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**2 - 18432.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(5) + 15.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 60.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - &
     & 420.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 1500.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + 5395.d0* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 16320.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 41230.d0* &
     & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 94560.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 198682.d0* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 222156.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 532320.d0* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - 1038300.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 161282.d0* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 13210320.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + 9336100.d0* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 73893720.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 52451405.d0 &
     & *gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 287438940.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + &
     & 187558132.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 865340796.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 &
     & - 502433633.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 2110814112.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) &
     & **2 + 1065610706.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 4273118160.d0*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)**2 - 1834610820.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 7271778300.d0*gl%z3_PCSAFT(1) &
     & **13*gl%z3_PCSAFT(4)**2 + 2593730240.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 10462467300.d0* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 3020711744.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + &
     & 12734595312.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 2886136328.d0*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) - 13066065432.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 2233829944.d0* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 11217968880.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + &
     & 1364462240.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 7967593440.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 &
     & - 622720720.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 4604185440.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)** &
     & 2 + 183443264.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 2112183072.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) &
     & **2 - 12598272.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 740144448.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) &
     & **2 - 17818112.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 184923840.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) &
     & **2 + 9341760.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 28049280.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)** &
     & 2 - 2028160.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 1059840.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + &
     & 66816.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 415488.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 56832.d0* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 56832.d0*gl%z3_PCSAFT(4)**2 - 8448.d0*gl%z3_PCSAFT(5))
	part9 = (-24.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26 + 1440.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**25 - 33580.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 + 451720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 23 - 4088254.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 + 26946048.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 - &
     & 135507961.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 + 535343930.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 - &
     & 1692677365.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 + 4331702180.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 - &
     & 9019781350.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 + 15278925396.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 &
     & - 20919521480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 + 22808767440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 13 - 19238718649.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 + 11839051746.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**11 - 4574037537.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 + 436284060.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**9 + 579547620.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 - 282604824.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**7 + 11364372.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 + 21380976.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**5 - 1779840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - 933120.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 3 - 82944.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 + 12.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - 1680.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + 48020.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 710600.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**23 + 6811463.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 46753328.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**21 + 242861602.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 987461760.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**19 + 3210652355.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 8460252756.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**17 + 18200677596.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 32035944776.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**15 + 45986766905.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 53308341800.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + 48927659666.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - &
     & 34231467568.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 + 16807280589.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 &
     & - 4433228180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 - 560219220.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + &
     & 963960768.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - 290326716.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - &
     & 11723760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 21592080.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - &
     & 1249920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 573696.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 55296.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) + 18.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 25 - 21615.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 386970.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - &
     & 4065414.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 29534232.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - &
     & 159893631.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 672729450.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - &
     & 2257036140.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 6136981464.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - &
     & 13654346625.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 24969060174.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & - 37509314010.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 46014333720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13 - 45485112585.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 35318443878.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**11 - 20474713086.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 7809273120.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**9 - 1021797480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 794402568.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**7 + 522922728.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 108593856.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**5 - 12577680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 9043200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **3 - 712512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 112896.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 9216.d0* &
     & gl%mmean_PCSAFT**2 - 15.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 + 180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 2390.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 79880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 1000937.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **22 - 7942080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 45499088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - &
     & 199785880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 695272235.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - &
     & 1957144588.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 4512074102.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - &
     & 8577500984.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 13472303835.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - &
     & 17433408680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 18423939236.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - &
     & 15625903592.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 10284386896.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - &
     & 4883271200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 1330128400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + &
     & 98773984.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 276830592.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 117253568.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 16201920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 4332800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **3 + 1976256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 200448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 8448.d0*gl%mmean_PCSAFT &
     & + 3.d0*gl%z3_PCSAFT(1)**26 - 60.d0*gl%z3_PCSAFT(1)**25 + 305.d0*gl%z3_PCSAFT(1)**24 + 3790.d0*gl%z3_PCSAFT(1)**23 &
     & - 80032.d0*gl%z3_PCSAFT(1)**22 + 747520.d0*gl%z3_PCSAFT(1)**21 - 4664594.d0*gl%z3_PCSAFT(1)**20 + &
     & 21684420.d0*gl%z3_PCSAFT(1)**19 - 78917605.d0*gl%z3_PCSAFT(1)**18 + 231136332.d0*gl%z3_PCSAFT(1)**17 &
     & - 553843971.d0*gl%z3_PCSAFT(1)**16 + 1096210718.d0*gl%z3_PCSAFT(1)**15 - 1800160670.d0* &
     & gl%z3_PCSAFT(1)**14 + 2452125440.d0*gl%z3_PCSAFT(1)**13 - 2755994204.d0*gl%z3_PCSAFT(1)**12 + &
     & 2524894048.d0*gl%z3_PCSAFT(1)**11 - 1841870016.d0*gl%z3_PCSAFT(1)**10 + 1020420320.d0*gl%z3_PCSAFT(1) &
     & **9 - 381042720.d0*gl%z3_PCSAFT(1)**8 + 51970560.d0*gl%z3_PCSAFT(1)**7 + 38957312.d0*gl%z3_PCSAFT(1) &
     & **6 - 29568896.d0*gl%z3_PCSAFT(1)**5 + 8868800.d0*gl%z3_PCSAFT(1)**4 - 510080.d0*gl%z3_PCSAFT(1)**3 - &
     & 515456.d0*gl%z3_PCSAFT(1)**2 + 165888.d0*gl%z3_PCSAFT(1) - 17152.d0)
	part10 = (-24.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26 + 1440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25 - &
     & 33580.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 + 451720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 - 4088254.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 + 26946048.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 - 135507961.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 + 535343930.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 - 1692677365.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 + 4331702180.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 - 9019781350.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 + 15278925396.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 - &
     & 20919521480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 + 22808767440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 &
     & - 19238718649.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 + 11839051746.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 11 - 4574037537.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 + 436284060.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 9 + 579547620.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 - 282604824.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 + &
     & 11364372.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 + 21380976.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 - &
     & 1779840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - 933120.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 - 82944.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 + 12.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - 1680.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**25 + 48020.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 710600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 23 + 6811463.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 46753328.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + &
     & 242861602.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 987461760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + &
     & 3210652355.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 8460252756.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + &
     & 18200677596.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 32035944776.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & + 45986766905.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 53308341800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 13 + 48927659666.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 34231467568.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11 + 16807280589.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 4433228180.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**9 - 560219220.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 963960768.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7 - 290326716.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 11723760.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**5 + 21592080.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 1249920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **3 - 573696.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 55296.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) + 18.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - 21615.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**24 + 386970.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - 4065414.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **22 + 29534232.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - 159893631.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 20 + 672729450.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 2257036140.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 18 + 6136981464.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 13654346625.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **16 + 24969060174.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 37509314010.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**14 + 46014333720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 45485112585.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**12 + 35318443878.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 20474713086.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 7809273120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 1021797480.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 794402568.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 522922728.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 108593856.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 12577680.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 9043200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 712512.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**2 - 112896.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 9216.d0*gl%mmean_PCSAFT**2 - 15.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**26 + 180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 2390.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 79880.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 1000937.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 7942080.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21 + 45499088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 199785880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 19 + 695272235.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 1957144588.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + &
     & 4512074102.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 8577500984.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + &
     & 13472303835.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 17433408680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + &
     & 18423939236.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 15625903592.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + &
     & 10284386896.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 4883271200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + &
     & 1330128400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 98773984.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 276830592.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 117253568.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 16201920.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**4 - 4332800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 1976256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - &
     & 200448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 8448.d0*gl%mmean_PCSAFT + 3.d0*gl%z3_PCSAFT(1)**26 - 60.d0*gl%z3_PCSAFT(1)**25 + &
     & 305.d0*gl%z3_PCSAFT(1)**24 + 3790.d0*gl%z3_PCSAFT(1)**23 - 80032.d0*gl%z3_PCSAFT(1)**22 + 747520.d0* &
     & gl%z3_PCSAFT(1)**21 - 4664594.d0*gl%z3_PCSAFT(1)**20 + 21684420.d0*gl%z3_PCSAFT(1)**19 - 78917605.d0* &
     & gl%z3_PCSAFT(1)**18 + 231136332.d0*gl%z3_PCSAFT(1)**17 - 553843971.d0*gl%z3_PCSAFT(1)**16 + &
     & 1096210718.d0*gl%z3_PCSAFT(1)**15 - 1800160670.d0*gl%z3_PCSAFT(1)**14 + 2452125440.d0*gl%z3_PCSAFT(1) &
     & **13 - 2755994204.d0*gl%z3_PCSAFT(1)**12 + 2524894048.d0*gl%z3_PCSAFT(1)**11 - 1841870016.d0 &
     & *gl%z3_PCSAFT(1)**10 + 1020420320.d0*gl%z3_PCSAFT(1)**9 - 381042720.d0*gl%z3_PCSAFT(1)**8 + &
     & 51970560.d0*gl%z3_PCSAFT(1)**7 + 38957312.d0*gl%z3_PCSAFT(1)**6 - 29568896.d0*gl%z3_PCSAFT(1)**5 + &
     & 8868800.d0*gl%z3_PCSAFT(1)**4 - 510080.d0*gl%z3_PCSAFT(1)**3 - 515456.d0*gl%z3_PCSAFT(1)**2 + 165888.d0* &
     & gl%z3_PCSAFT(1) - 17152.d0)
	part11 = (-48.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 26 + 2880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25 - 67160.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 + &
     & 903440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 - 8176508.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 + &
     & 53892096.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 - 271015922.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 + &
     & 1070687860.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 - 3385354730.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 + &
     & 8663404360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 - 18039562700.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 &
     & + 30557850792.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 - 41839042960.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 14 + 45617534880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 - 38477437298.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**12 + 23678103492.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 - 9148075074.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)**10 + 872568120.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9 + 1159095240.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**8 - 565209648.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 + 22728744.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**6 + 42761952.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 - 3559680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **4 - 1866240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 - 165888.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 + 24.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - 3360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + 96040.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**24 - 1421200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + 13622926.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**22 - 93506656.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + 485723204.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**20 - 1974923520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + 6421304710.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**18 - 16920505512.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + 36401355192.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**16 - 64071889552.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + 91973533810.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 106616683600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + &
     & 97855319332.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 68462935136.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 &
     & + 33614561178.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 8866456360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 &
     & - 1120438440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 1927921536.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - &
     & 580653432.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 23447520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + &
     & 43184160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 2499840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - &
     & 1147392.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 110592.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) + 36.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**26 + 720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - 43230.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 &
     & + 773940.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - 8130828.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + &
     & 59068464.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - 319787262.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + &
     & 1345458900.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 4514072280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + &
     & 12273962928.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 27308693250.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 &
     & + 49938120348.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 75018628020.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 14 + 92028667440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 90970225170.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**12 + 70636887756.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 40949426172.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10 + 15618546240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 2043594960.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 1588805136.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 1045845456.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 217187712.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 25155360.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 18086400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 1425024.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**2 - 225792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 18432.d0*gl%mmean_PCSAFT**2 - 30.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 + 360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 4780.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - &
     & 159760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 2001874.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 15884160.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 90998176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 399571760.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19 + 1390544470.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 3914289176.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**17 + 9024148204.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 17155001968.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**15 + 26944607670.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 34866817360.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13 + 36847878472.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 31251807184.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11 + 20568773792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 9766542400.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9 + 2660256800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 197547968.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 7 - 553661184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 234507136.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - &
     & 32403840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 8665600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 3952512.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 400896.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 16896.d0*gl%mmean_PCSAFT + 6.d0*gl%z3_PCSAFT(1)** &
     & 26 - 120.d0*gl%z3_PCSAFT(1)**25 + 610.d0*gl%z3_PCSAFT(1)**24 + 7580.d0*gl%z3_PCSAFT(1)**23 - 160064.d0* &
     & gl%z3_PCSAFT(1)**22 + 1495040.d0*gl%z3_PCSAFT(1)**21 - 9329188.d0*gl%z3_PCSAFT(1)**20 + 43368840.d0* &
     & gl%z3_PCSAFT(1)**19 - 157835210.d0*gl%z3_PCSAFT(1)**18 + 462272664.d0*gl%z3_PCSAFT(1)**17 - &
     & 1107687942.d0*gl%z3_PCSAFT(1)**16 + 2192421436.d0*gl%z3_PCSAFT(1)**15 - 3600321340.d0*gl%z3_PCSAFT(1) &
     & **14 + 4904250880.d0*gl%z3_PCSAFT(1)**13 - 5511988408.d0*gl%z3_PCSAFT(1)**12 + 5049788096.d0 &
     & *gl%z3_PCSAFT(1)**11 - 3683740032.d0*gl%z3_PCSAFT(1)**10 + 2040840640.d0*gl%z3_PCSAFT(1)**9 - &
     & 762085440.d0*gl%z3_PCSAFT(1)**8 + 103941120.d0*gl%z3_PCSAFT(1)**7 + 77914624.d0*gl%z3_PCSAFT(1)**6 - &
     & 59137792.d0*gl%z3_PCSAFT(1)**5 + 17737600.d0*gl%z3_PCSAFT(1)**4 - 1020160.d0*gl%z3_PCSAFT(1)**3 - &
     & 1030912.d0*gl%z3_PCSAFT(1)**2 + 331776.d0*gl%z3_PCSAFT(1) - 34304.d0)
	part12 = (192.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 12480.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 290880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - &
     & 3823680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 33552144.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 21*gl%z3_PCSAFT(4) - 214150224.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 1045931280.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 4037332560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) + 12577733280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 32038746720.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 67171056480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) - 116001761280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + &
     & 164161740720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 188058615600.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 170579587824.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) &
     & - 117800370768.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 57289536960.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 15814573920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - &
     & 252378720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 1784186784.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 375251616.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 63374400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 10679040.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(4) + 5598720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 497664.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 96.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 15840.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 443040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + &
     & 6304800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 57883608.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 21*gl%z3_PCSAFT(4) + 380911704.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 1905101880.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 7509508920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) - 23885385840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 62232793296.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 133973148528.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) + 239055840240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - &
     & 352851779640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 427604437080.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 419438600376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) &
     & + 325077706104.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 190838142720.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 77869317600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) &
     & - 17011774080.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 1389619008.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 2065700448.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 456096960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 241920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 &
     & *gl%z3_PCSAFT(4) + 3939840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 1575936.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 165888.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4) - 144.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 5040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 228240.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 3721680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4) + 36496584.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 250123896.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 1288517400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 5206451400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 16952630400.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 45263898144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 100166770080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - &
     & 184693618080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 283909717800.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 362367916920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) &
     & + 380528924760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 323626210344.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 217133845920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) - 109661160960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 37601402400.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 5915968992.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4) - 1509309504.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 1112601600.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 249615360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) + 15275520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 1562112.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) + 120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 25*gl%z3_PCSAFT(4) - 600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 41400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4) + 904920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 9860112.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 71519280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - &
     & 382704960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 1593350400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 &
     & *gl%z3_PCSAFT(4) - 5329398600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 14616424008.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 33296489976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) &
     & + 63462451800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 101475307200.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 135917679840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - &
     & 151621824336.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 139357997136.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 103658938080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 60537853440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 26197786560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 7*gl%z3_PCSAFT(4) + 7248917376.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 481137024.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 575765760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 273192960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 55296000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) + 4145664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 113664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) &
     & - 24.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 360.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 840.d0*gl%z3_PCSAFT(1)** &
     & 23*gl%z3_PCSAFT(4) - 71160.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 944352.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) - 7458240.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 41990160.d0*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 181417680.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 626106600.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) - 1768915032.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 4155006888.d0*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4) - 8189712600.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 13606732560.d0*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(4) - 19064440080.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 22449767424.d0* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 22051778784.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + &
     & 17836793280.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 11630024640.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + &
     & 5886940800.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 2135900160.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + &
     & 427344384.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 45580800.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - &
     & 67438080.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 23907840.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - &
     & 4291584.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 331776.d0*gl%z3_PCSAFT(4))
	part13 = (192.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 12480.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 290880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - &
     & 3823680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 33552144.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 21*gl%z3_PCSAFT(4) - 214150224.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 1045931280.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 4037332560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) + 12577733280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 32038746720.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 67171056480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) - 116001761280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + &
     & 164161740720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 188058615600.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 170579587824.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) &
     & - 117800370768.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 57289536960.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 15814573920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - &
     & 252378720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 1784186784.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 375251616.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 63374400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 10679040.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(4) + 5598720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 497664.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 96.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 15840.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 443040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + &
     & 6304800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 57883608.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 21*gl%z3_PCSAFT(4) + 380911704.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 1905101880.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 7509508920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) - 23885385840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 62232793296.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 133973148528.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) + 239055840240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - &
     & 352851779640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 427604437080.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 419438600376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) &
     & + 325077706104.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 190838142720.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 77869317600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) &
     & - 17011774080.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 1389619008.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 2065700448.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 456096960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 241920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 &
     & *gl%z3_PCSAFT(4) + 3939840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 1575936.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 165888.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4) - 144.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 5040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 228240.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 3721680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4) + 36496584.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 250123896.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 1288517400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 5206451400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 16952630400.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 45263898144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 100166770080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - &
     & 184693618080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 283909717800.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 362367916920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) &
     & + 380528924760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 323626210344.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 217133845920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) - 109661160960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 37601402400.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 5915968992.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4) - 1509309504.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 1112601600.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 249615360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) + 15275520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 1562112.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) + 120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 25*gl%z3_PCSAFT(4) - 600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 41400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4) + 904920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 9860112.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 71519280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - &
     & 382704960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 1593350400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 &
     & *gl%z3_PCSAFT(4) - 5329398600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 14616424008.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 33296489976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) &
     & + 63462451800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 101475307200.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 135917679840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - &
     & 151621824336.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 139357997136.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 103658938080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 60537853440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 26197786560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 7*gl%z3_PCSAFT(4) + 7248917376.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 481137024.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 575765760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 273192960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 55296000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) + 4145664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 113664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) &
     & - 24.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 360.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 840.d0*gl%z3_PCSAFT(1)** &
     & 23*gl%z3_PCSAFT(4) - 71160.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 944352.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) - 7458240.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 41990160.d0*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 181417680.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 626106600.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) - 1768915032.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 4155006888.d0*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4) - 8189712600.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 13606732560.d0*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(4) - 19064440080.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 22449767424.d0* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 22051778784.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + &
     & 17836793280.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 11630024640.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + &
     & 5886940800.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 2135900160.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + &
     & 427344384.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 45580800.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - &
     & 67438080.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 23907840.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - &
     & 4291584.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 331776.d0*gl%z3_PCSAFT(4))
	part14 = (96.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 480.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 6240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + &
     & 33600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 145440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 23*gl%z3_PCSAFT(5) - 774720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 - 1911840.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) + 9880320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4)**2 + 16776072.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - 83777208.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 107075112.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(5) + 518458800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + 522965640.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 2474565000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4)**2 - 2018666280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + &
     & 9431202720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 6288866640.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 29333020560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)** &
     & 2 - 16019373360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 75384706080.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + 33585528240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(5) - 160954936080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - &
     & 58000880640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 285445862880.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + 82080870360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(5) - 418296100680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - &
     & 94029307800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + 501298508400.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 85289793912.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) - 483148526808.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - &
     & 58900185384.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 364731247200.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 28644768480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) &
     & - 206512323840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 7907286960.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 81002473920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) &
     & **2 - 126189360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 18305619120.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 892093392.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(5) + 962213472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 187625808.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 231843600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4)**2 - 31687200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 126748800.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 5339520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(5) - 21358080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 2799360.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 11197440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
     & **2 + 248832.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 995328.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(4) &
     & **2 - 48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4)**2 + 7920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 45600.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 221520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + &
     & 1236240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 3152400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **22*gl%z3_PCSAFT(5) - 16828800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - &
     & 28941804.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 147694836.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 190455852.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) &
     & - 934006440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 952550940.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 4530327180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 &
     & + 3754754460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - 17516769840.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 11942692920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(5) + 55316597400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + &
     & 31116396648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 144756177072.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 66986574264.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(5) + 316247600760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + &
     & 119527920120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 578004045120.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 176425889820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(5) + 881858164740.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + &
     & 213802218540.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 1116432101640.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 209719300188.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) + 1161220259532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + &
     & 162538853052.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 977736633360.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 95419071360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(5) + 652195229040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + &
     & 38934658800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 333701936640.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 8505887040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + &
     & 124396843680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 694809504.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 30776517504.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 &
     & + 1032850224.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 4032573840.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 228048480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - &
     & 29520000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 120960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(5) - 41765760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 1969920.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 760320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)** &
     & 2 + 787968.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 580608.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4)** &
     & 2 + 82944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(5) - 72.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + &
     & 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 2520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(5) + 18000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 114120.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 693360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 &
     & - 1860840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) + 10451520.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 18248292.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - &
     & 96243948.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 125061948.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 626296680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 &
     & + 644258700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 3100746420.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 2603225700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) &
     & + 12202475760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 8476315200.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 39226842720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4)**2 - 22631949072.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + &
     & 104745004128.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + 50083385040.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 234507806280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)**2 - 92346809040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + &
     & 441980411040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + 141954858900.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 701277133500.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)**2 - 181183958460.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + &
     & 933775773480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 190264462380.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 1037055940740.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)**2 - 161813105172.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + &
     & 951731562960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 108566922960.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 712154731680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4)**2 - 54830580480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + &
     & 426292191360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + 18800701200.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 198508224240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) &
     & **2 - 2957984496.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 68845453344.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 754654752.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(5) - 16485539040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + &
     & 556300800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 2314529280.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 124807680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - &
     & 98904960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 7637760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **2*gl%z3_PCSAFT(5) - 11289600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 781056.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 71424.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**2 + 13824.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(5) + 60.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 300.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 300.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 600.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 20700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + &
     & 153540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 452460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(5) - 2770080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 4930056.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 27361344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + &
     & 35759640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 185019360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4)**2 - 191352480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 940262520.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 796675200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(5) - 3780353760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 2664699300.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 12405215460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4)**2 + 7308212004.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 33868011576.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 16648244988.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(5) + 77795442900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + &
     & 31731225900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 151201608960.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 50737653600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + &
     & 249075491640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + 67958839920.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 347300756880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 &
     & - 75810912168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + 408319680432.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 69678998568.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - &
     & 402194923200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 51829469040.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 328830663120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + &
     & 30268926720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 220258643520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **7*gl%z3_PCSAFT(4)**2 - 13098893280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + &
     & 118671262080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 3624458688.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 50087747712.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - &
     & 240568512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 15912730560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 &
     & *gl%z3_PCSAFT(4)**2 - 287882880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 3563527680.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 136596480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) &
     & + 496128000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 27648000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 2*gl%z3_PCSAFT(5) - 30612480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 2072832.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 415488.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**2 + 56832.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(5) - 12.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 60.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + &
     & 180.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 600.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 420.d0* &
     & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 9540.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 - 35580.d0*gl%z3_PCSAFT(1) &
     & **22*gl%z3_PCSAFT(5) + 258240.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 472176.d0*gl%z3_PCSAFT(1)** &
     & 21*gl%z3_PCSAFT(5) - 2845464.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 3729120.d0*gl%z3_PCSAFT(1)** &
     & 20*gl%z3_PCSAFT(5) + 20282640.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + 20995080.d0*gl%z3_PCSAFT(1) &
     & **19*gl%z3_PCSAFT(5) - 106575240.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 90708840.d0* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 439704000.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + &
     & 313053300.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 1477157940.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)** &
     & 2 - 884457516.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 4130732904.d0*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4)**2 + 2077503444.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 9742040340.d0*gl%z3_PCSAFT(1) &
     & **14*gl%z3_PCSAFT(4)**2 - 4094856300.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 19515883680.d0* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + 6803366280.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - &
     & 33307509840.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 9532220040.d0*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(5) + 48430636320.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 11224883712.d0* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 59856655008.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - &
     & 11025889392.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 62601407520.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
     & **2 + 8918396640.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 55032395040.d0*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4)**2 - 5815012320.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 40282467840.d0*gl%z3_PCSAFT(1) &
     & **7*gl%z3_PCSAFT(4)**2 + 2943470400.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 24233623680.d0* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 1067950080.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + &
     & 11766783360.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + 213672192.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) &
     & - 4493738880.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 22790400.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) &
     & + 1298611200.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 33719040.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) &
     & - 266760960.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 11953920.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) &
     & + 34690560.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 2145792.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - &
     & 2145792.d0*gl%z3_PCSAFT(4)**2 + 165888.d0*gl%z3_PCSAFT(5))
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
            	    gl%cx2_PCSAFT(5,xi,xj) = part1*(gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*part2 &
         & + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(5,xi)*part3 + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(4,xi)*part4 + gl%mPCSAFT(xj)* &
         & gl%z3x1_PCSAFT(1,xi)*part5 + gl%mPCSAFT(xi)* &
         & gl%z3x1_PCSAFT(5,xj)*part6 + gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(4,xj)*part7 + gl%mPCSAFT(xi)* &
         & gl%z3x1_PCSAFT(1,xj)*part8 + &
         & gl%z3x1_PCSAFT(5,xj)*gl%z3x1_PCSAFT(1,xi)*part9 + gl%z3x1_PCSAFT(5,xi)* &
         & gl%z3x1_PCSAFT(1,xj)*part10 + gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(4,xi)*part11 + gl%z3x1_PCSAFT(4,xj)* &
         & gl%z3x1_PCSAFT(1,xi)*part12 + gl%z3x1_PCSAFT(4,xi)* &
         & gl%z3x1_PCSAFT(1,xj)*part13 + gl%z3x1_PCSAFT(1,xj)* &
         & gl%z3x1_PCSAFT(1,xi)*part14)
                end if
            end do
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF C_1 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires z3_PCSAFT(1), z3_PCSAFT(4)
    if (GETDERC(6) .eq. 1) then
        part1 = 2.d0/recurring_factor**5
        part2 = (32.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(4) - 2000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **26*gl%z3_PCSAFT(4) + 46896.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 627240.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 5614424.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - &
     & 36493284.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 180728484.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **21*gl%z3_PCSAFT(4) - 703054862.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + &
     & 2191643200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 5545326430.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 11470366320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - &
     & 19433657804.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 26861059096.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 29929934920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + &
     & 26235334084.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 17231567678.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 7572428928.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - &
     & 1407565974.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 585221640.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 433857024.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - &
     & 59502096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 23607504.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 6*gl%z3_PCSAFT(4) + 4466880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 933120.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 82944.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 112.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(4) - 2680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + &
     & 26064.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 100360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4) - 394564.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 7796442.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 55150224.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + &
     & 255359452.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 874424860.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 2331941338.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 4991031224.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 8766775032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
     & - 12896391500.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 16198091950.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 17615546688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 16567075452.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 13118140468.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) + 8233419142.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 3643196360.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 800496696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + &
     & 163352304.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 162838320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4) + 28478592.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 6298560.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 1700352.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 193536.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 64.d0*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(4) + 1840.d0*gl%z3_PCSAFT(1)** &
     & 26*gl%z3_PCSAFT(4) - 24432.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 197180.d0*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4) - 1064196.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 3918294.d0*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4) - 9079636.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 6252694.d0*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) + 48664300.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 259630744.d0*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) + 779297160.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 1719100696.d0*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 3006596100.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 4307860530.d0*gl%z3_PCSAFT(1)**14 &
     & *gl%z3_PCSAFT(4) + 5132474940.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 5104158914.d0*gl%z3_PCSAFT(1)** &
     & 12*gl%z3_PCSAFT(4) + 4216668884.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 2858220700.d0*gl%z3_PCSAFT(1) &
     & **10*gl%z3_PCSAFT(4) + 1556710400.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 656164792.d0*gl%z3_PCSAFT(1) &
     & **8*gl%z3_PCSAFT(4) + 195580096.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 27938336.d0*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4) - 7401792.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 5390400.d0*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4) - 1102848.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 41472.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & + 36864.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4))
        part3 = (-16.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 27 + 1000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - 23448.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + &
     & 313620.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 2807212.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + &
     & 18246642.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 90364242.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + &
     & 351527431.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 1095821600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + &
     & 2772663215.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 5735183160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + &
     & 9716828902.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 13430529548.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & + 14964967460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 13117667042.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 13 + 8615783839.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 3786214464.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **11 + 703782987.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + 292610820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **9 - 216928512.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 29751048.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 &
     & + 11803752.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 2233440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - &
     & 466560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 41472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 24.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**27 + 60.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 10872.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 &
     & - 229080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 2496198.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - &
     & 18110109.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 96140838.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - &
     & 392293269.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 1267303170.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - &
     & 3303484746.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 7034682828.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - &
     & 12328763754.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 17836871250.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & - 21265685565.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 20724474366.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13 - 16216090629.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 9827320686.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **11 - 4261587864.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 1029848760.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**9 + 90875868.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 173814768.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**7 + 53388960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 896544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **5 - 2858400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 359424.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + &
     & 41472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 36.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 - 810.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **26 + 6678.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - 4995.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 407976.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 4587969.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 29448336.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21 + 133953354.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 465275880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **19 + 1281207036.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 2860377330.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 &
     & + 5254199289.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 8019619200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + &
     & 10231558785.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 10927205940.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + &
     & 9719912256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 7088077476.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + &
     & 4089330540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 1721345520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + &
     & 409013328.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 39299616.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 73052496.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 23768448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 1117440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **4 - 1077888.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 160128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 13824.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 10.d0*gl%z3_PCSAFT(1)**27 + 265.d0*gl%z3_PCSAFT(1)**26 - 3126.d0*gl%z3_PCSAFT(1)**25 + &
     & 20565.d0*gl%z3_PCSAFT(1)**24 - 69106.d0*gl%z3_PCSAFT(1)**23 - 48750.d0*gl%z3_PCSAFT(1)**22 + 2003424.d0* &
     & gl%z3_PCSAFT(1)**21 - 12999494.d0*gl%z3_PCSAFT(1)**20 + 53688970.d0*gl%z3_PCSAFT(1)**19 - &
     & 165206811.d0*gl%z3_PCSAFT(1)**18 + 401380802.d0*gl%z3_PCSAFT(1)**17 - 793227527.d0*gl%z3_PCSAFT(1)** &
     & 16 + 1297925122.d0*gl%z3_PCSAFT(1)**15 - 1777365840.d0*gl%z3_PCSAFT(1)**14 + 2047341652.d0* &
     & gl%z3_PCSAFT(1)**13 - 1981704872.d0*gl%z3_PCSAFT(1)**12 + 1597441568.d0*gl%z3_PCSAFT(1)**11 - &
     & 1050521584.d0*gl%z3_PCSAFT(1)**10 + 540924400.d0*gl%z3_PCSAFT(1)**9 - 198986592.d0*gl%z3_PCSAFT(1)**8 &
     & + 37752608.d0*gl%z3_PCSAFT(1)**7 + 7591424.d0*gl%z3_PCSAFT(1)**6 - 8296832.d0*gl%z3_PCSAFT(1)**5 + &
     & 2732480.d0*gl%z3_PCSAFT(1)**4 - 352768.d0*gl%z3_PCSAFT(1)**3 - 35328.d0*gl%z3_PCSAFT(1)**2 + 16896.d0* &
     & gl%z3_PCSAFT(1) - 1536.d0)
	part4 = (48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) &
     & - 3360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 78280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4) - 1008400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 8599564.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 53183016.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + &
     & 251949718.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 947978420.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 2903511910.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - &
     & 7355969000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 15545965540.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 27443029848.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + &
     & 40241827400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 48411772920.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 46812356614.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - &
     & 35222081892.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 19496693406.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 7034718840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 1032905880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 326883744.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 164897064.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + &
     & 11074752.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 1779840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 &
     & *gl%z3_PCSAFT(4) + 933120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 82944.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 72.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 720.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 55080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + &
     & 929400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 8885586.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 &
     & *gl%z3_PCSAFT(4) + 58022112.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 281620122.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 1067875380.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 3270736260.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 8283476832.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 17623834884.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 31784676576.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 48721323810.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 63292400880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) - 69161009562.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 62826035652.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 46624825104.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) + 27479320320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 12184201920.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 3558064944.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) - 357911616.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 182931840.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 69824160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4) - 1134720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 1963008.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 193536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 108.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 1800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 3090.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 188400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + &
     & 2643642.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 19929552.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) + 103701504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 408193380.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 1279009680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - &
     & 3301443672.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 7194693150.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 16*gl%z3_PCSAFT(4) - 13455239256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 21769063530.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 30470009520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) + 36648154260.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 37426323108.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 31929173784.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) &
     & - 22282843920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 12368552640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **8*gl%z3_PCSAFT(4) - 5239619856.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 1576584192.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 281856576.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 5062080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 11702400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) - 2899584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 170496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & *gl%z3_PCSAFT(4) + 36864.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) + 30.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 660.d0* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 5530.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 12100.d0*gl%z3_PCSAFT(1)**23 &
     & *gl%z3_PCSAFT(4) - 175208.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 2102940.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) - 12887756.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 55221520.d0*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 182536130.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 490224532.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) - 1105946846.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 2141896748.d0*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4) - 3602556660.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 5275006820.d0*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(4) - 6693171824.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 7293792776.d0*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) - 6750308992.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 5238668960.d0* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 3358744000.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 1745296544.d0* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 714947904.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 220560064.d0* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 46732800.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 5116160.d0*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4) + 281856.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 170496.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
     & + 16896.d0*gl%z3_PCSAFT(4))
	part5 = (-16.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27 + &
     & 1000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - 23448.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + 313620.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 2807212.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + 18246642.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 90364242.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + 351527431.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 1095821600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + 2772663215.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 5735183160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + &
     & 9716828902.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 13430529548.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & + 14964967460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 13117667042.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 13 + 8615783839.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 3786214464.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **11 + 703782987.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + 292610820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **9 - 216928512.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 29751048.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 &
     & + 11803752.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 2233440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - &
     & 466560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 41472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 24.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**27 + 60.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 10872.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 &
     & - 229080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 2496198.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - &
     & 18110109.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 96140838.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - &
     & 392293269.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 1267303170.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - &
     & 3303484746.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 7034682828.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - &
     & 12328763754.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 17836871250.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & - 21265685565.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 20724474366.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13 - 16216090629.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 9827320686.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **11 - 4261587864.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 1029848760.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**9 + 90875868.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 173814768.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**7 + 53388960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 896544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **5 - 2858400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 359424.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + &
     & 41472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 36.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 - 810.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **26 + 6678.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - 4995.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 407976.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 4587969.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 29448336.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21 + 133953354.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 465275880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **19 + 1281207036.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 2860377330.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 &
     & + 5254199289.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 8019619200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + &
     & 10231558785.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 10927205940.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + &
     & 9719912256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 7088077476.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + &
     & 4089330540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 1721345520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + &
     & 409013328.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 39299616.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 73052496.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 23768448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 1117440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **4 - 1077888.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 160128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 13824.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 10.d0*gl%z3_PCSAFT(1)**27 + 265.d0*gl%z3_PCSAFT(1)**26 - 3126.d0*gl%z3_PCSAFT(1)**25 + &
     & 20565.d0*gl%z3_PCSAFT(1)**24 - 69106.d0*gl%z3_PCSAFT(1)**23 - 48750.d0*gl%z3_PCSAFT(1)**22 + 2003424.d0* &
     & gl%z3_PCSAFT(1)**21 - 12999494.d0*gl%z3_PCSAFT(1)**20 + 53688970.d0*gl%z3_PCSAFT(1)**19 - &
     & 165206811.d0*gl%z3_PCSAFT(1)**18 + 401380802.d0*gl%z3_PCSAFT(1)**17 - 793227527.d0*gl%z3_PCSAFT(1)** &
     & 16 + 1297925122.d0*gl%z3_PCSAFT(1)**15 - 1777365840.d0*gl%z3_PCSAFT(1)**14 + 2047341652.d0* &
     & gl%z3_PCSAFT(1)**13 - 1981704872.d0*gl%z3_PCSAFT(1)**12 + 1597441568.d0*gl%z3_PCSAFT(1)**11 - &
     & 1050521584.d0*gl%z3_PCSAFT(1)**10 + 540924400.d0*gl%z3_PCSAFT(1)**9 - 198986592.d0*gl%z3_PCSAFT(1)**8 &
     & + 37752608.d0*gl%z3_PCSAFT(1)**7 + 7591424.d0*gl%z3_PCSAFT(1)**6 - 8296832.d0*gl%z3_PCSAFT(1)**5 + &
     & 2732480.d0*gl%z3_PCSAFT(1)**4 - 352768.d0*gl%z3_PCSAFT(1)**3 - 35328.d0*gl%z3_PCSAFT(1)**2 + 16896.d0* &
     & gl%z3_PCSAFT(1) - 1536.d0)
	part6 = (48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) &
     & - 3360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 78280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4) - 1008400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 8599564.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 53183016.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + &
     & 251949718.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 947978420.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 2903511910.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - &
     & 7355969000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 15545965540.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 27443029848.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + &
     & 40241827400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 48411772920.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 46812356614.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - &
     & 35222081892.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 19496693406.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 7034718840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 1032905880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 326883744.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 164897064.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + &
     & 11074752.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 1779840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 &
     & *gl%z3_PCSAFT(4) + 933120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 82944.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 72.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 720.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 55080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + &
     & 929400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 8885586.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 &
     & *gl%z3_PCSAFT(4) + 58022112.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 281620122.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 1067875380.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 3270736260.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 8283476832.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 17623834884.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 31784676576.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 48721323810.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 63292400880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) - 69161009562.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 62826035652.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 46624825104.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) + 27479320320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 12184201920.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 3558064944.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) - 357911616.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 182931840.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 69824160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4) - 1134720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 1963008.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 193536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 108.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 1800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 3090.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 188400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + &
     & 2643642.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 19929552.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) + 103701504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 408193380.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 1279009680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - &
     & 3301443672.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 7194693150.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 16*gl%z3_PCSAFT(4) - 13455239256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 21769063530.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 30470009520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) + 36648154260.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 37426323108.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 31929173784.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) &
     & - 22282843920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 12368552640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **8*gl%z3_PCSAFT(4) - 5239619856.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 1576584192.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 281856576.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 5062080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 11702400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) - 2899584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 170496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & *gl%z3_PCSAFT(4) + 36864.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) + 30.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 660.d0* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 5530.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 12100.d0*gl%z3_PCSAFT(1)**23 &
     & *gl%z3_PCSAFT(4) - 175208.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 2102940.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) - 12887756.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 55221520.d0*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 182536130.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 490224532.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) - 1105946846.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 2141896748.d0*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4) - 3602556660.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 5275006820.d0*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(4) - 6693171824.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 7293792776.d0*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) - 6750308992.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 5238668960.d0* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 3358744000.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 1745296544.d0* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 714947904.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 220560064.d0* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 46732800.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 5116160.d0*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4) + 281856.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 170496.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
     & + 16896.d0*gl%z3_PCSAFT(4))
	part7 = (-48.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 26 + 3360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25 - 78280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 + &
     & 1008400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 - 8599564.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 + &
     & 53183016.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 - 251949718.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 + &
     & 947978420.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 - 2903511910.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 + &
     & 7355969000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 - 15545965540.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 &
     & + 27443029848.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 - 40241827400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 14 + 48411772920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 - 46812356614.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**12 + 35222081892.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 - 19496693406.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**10 + 7034718840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9 - 1032905880.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**8 - 326883744.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 + 164897064.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**6 - 11074752.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 - 1779840.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**4 - 933120.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 - 82944.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 &
     & + 24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - 4560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + 125480.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 1731200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + 15318878.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 96949196.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + 466827736.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 1779830940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + 5521388210.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 14195891136.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + &
     & 30585219072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 55456030568.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & + 84452356010.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 107185534940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **13 + 111863980856.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 94075917916.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11 + 61804510182.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 30068202440.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**9 + 9626325480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - 1233112032.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**7 - 452196792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 251496000.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**5 - 43063200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 529920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **3 + 359424.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 + 27648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) + 36.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 1800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - 70890.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**24 + 1086900.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - 10117464.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**22 + 65993484.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - 324471438.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**20 + 1257766800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 3962242920.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**18 + 10357986144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 22774691790.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**16 + 42408688692.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 66936230880.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 89155291020.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - &
     & 99294237210.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 91176217416.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 &
     & - 67617496788.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 39212034240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 9 - 16757106240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 4546789632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 7 - 291190704.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 339113088.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + &
     & 149963040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 25724160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + &
     & 643968.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 211968.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 18432.d0*gl%mmean_PCSAFT &
     & **2 - 30.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 60.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 15920.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**24 - 292700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 2928182.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - &
     & 19875480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 100354304.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - &
     & 397103440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 1274154830.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - &
     & 3393922828.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 7624096784.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - &
     & 14576223932.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 23793045930.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - &
     & 33092022560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 38963033696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - &
     & 38427191384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 31260695248.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - &
     & 20502384320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 10438636480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - &
     & 3822006656.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 794229696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 53375744.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 104192640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 36313600.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**3 - 6025344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 344064.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 16896.d0 &
     & *gl%mmean_PCSAFT + 6.d0*gl%z3_PCSAFT(1)**26 - 60.d0*gl%z3_PCSAFT(1)**25 - 1030.d0*gl%z3_PCSAFT(1)**24 + 28000.d0* &
     & gl%z3_PCSAFT(1)**23 - 312112.d0*gl%z3_PCSAFT(1)**22 + 2234080.d0*gl%z3_PCSAFT(1)**21 - 11665892.d0* &
     & gl%z3_PCSAFT(1)**20 + 47340000.d0*gl%z3_PCSAFT(1)**19 - 155218090.d0*gl%z3_PCSAFT(1)**18 + &
     & 422184852.d0*gl%z3_PCSAFT(1)**17 - 969815502.d0*gl%z3_PCSAFT(1)**16 + 1902434864.d0*gl%z3_PCSAFT(1)** &
     & 15 - 3203044940.d0*gl%z3_PCSAFT(1)**14 + 4627969160.d0*gl%z3_PCSAFT(1)**13 - 5712895304.d0* &
     & gl%z3_PCSAFT(1)**12 + 5976101296.d0*gl%z3_PCSAFT(1)**11 - 5234656608.d0*gl%z3_PCSAFT(1)**10 + &
     & 3774171680.d0*gl%z3_PCSAFT(1)**9 - 2181384960.d0*gl%z3_PCSAFT(1)**8 + 964008960.d0*gl%z3_PCSAFT(1)**7 &
     & - 291586816.d0*gl%z3_PCSAFT(1)**6 + 36347392.d0*gl%z3_PCSAFT(1)**5 + 15981440.d0*gl%z3_PCSAFT(1)**4 - &
     & 10933760.d0*gl%z3_PCSAFT(1)**3 + 3176704.d0*gl%z3_PCSAFT(1)**2 - 497664.d0*gl%z3_PCSAFT(1) + 34304.d0)
	part8 = (-48.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26 + 3360.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**25 - 78280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 + 1008400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **23 - 8599564.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 + 53183016.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 &
     & - 251949718.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 + 947978420.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 - &
     & 2903511910.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 + 7355969000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 - &
     & 15545965540.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 + 27443029848.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 &
     & - 40241827400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 + 48411772920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 13 - 46812356614.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 + 35222081892.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**11 - 19496693406.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 + 7034718840.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)**9 - 1032905880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 - 326883744.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**7 + 164897064.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 - 11074752.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**5 - 1779840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - 933120.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 3 - 82944.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 + 24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - 4560.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + 125480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 1731200.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**23 + 15318878.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 96949196.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**21 + 466827736.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 1779830940.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**19 + 5521388210.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 14195891136.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**17 + 30585219072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 55456030568.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + 84452356010.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - &
     & 107185534940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + 111863980856.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 12 - 94075917916.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 + 61804510182.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**10 - 30068202440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 9626325480.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**8 - 1233112032.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - 452196792.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**6 + 251496000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 43063200.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**4 + 529920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 359424.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 &
     & + 27648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) + 36.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 1800.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**25 - 70890.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 1086900.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **23 - 10117464.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 65993484.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 21 - 324471438.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 1257766800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 19 - 3962242920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 10357986144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **17 - 22774691790.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 42408688692.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**15 - 66936230880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 89155291020.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**13 - 99294237210.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 91176217416.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 67617496788.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + &
     & 39212034240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 16757106240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + &
     & 4546789632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 291190704.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - &
     & 339113088.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 149963040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - &
     & 25724160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 643968.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 211968.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 18432.d0*gl%mmean_PCSAFT**2 - 30.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 60.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**25 + 15920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 292700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + &
     & 2928182.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 19875480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 100354304.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 397103440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 1274154830.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**18 - 3393922828.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 7624096784.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**16 - 14576223932.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 23793045930.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14 - 33092022560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 38963033696.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12 - 38427191384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 31260695248.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10 - 20502384320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 10438636480.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**8 - 3822006656.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 794229696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 6 + 53375744.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 104192640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + &
     & 36313600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 6025344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 344064.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1) + 16896.d0*gl%mmean_PCSAFT + 6.d0*gl%z3_PCSAFT(1)**26 - 60.d0*gl%z3_PCSAFT(1)**25 - 1030.d0* &
     & gl%z3_PCSAFT(1)**24 + 28000.d0*gl%z3_PCSAFT(1)**23 - 312112.d0*gl%z3_PCSAFT(1)**22 + 2234080.d0* &
     & gl%z3_PCSAFT(1)**21 - 11665892.d0*gl%z3_PCSAFT(1)**20 + 47340000.d0*gl%z3_PCSAFT(1)**19 - &
     & 155218090.d0*gl%z3_PCSAFT(1)**18 + 422184852.d0*gl%z3_PCSAFT(1)**17 - 969815502.d0*gl%z3_PCSAFT(1)** &
     & 16 + 1902434864.d0*gl%z3_PCSAFT(1)**15 - 3203044940.d0*gl%z3_PCSAFT(1)**14 + 4627969160.d0* &
     & gl%z3_PCSAFT(1)**13 - 5712895304.d0*gl%z3_PCSAFT(1)**12 + 5976101296.d0*gl%z3_PCSAFT(1)**11 - &
     & 5234656608.d0*gl%z3_PCSAFT(1)**10 + 3774171680.d0*gl%z3_PCSAFT(1)**9 - 2181384960.d0*gl%z3_PCSAFT(1) &
     & **8 + 964008960.d0*gl%z3_PCSAFT(1)**7 - 291586816.d0*gl%z3_PCSAFT(1)**6 + 36347392.d0*gl%z3_PCSAFT(1) &
     & **5 + 15981440.d0*gl%z3_PCSAFT(1)**4 - 10933760.d0*gl%z3_PCSAFT(1)**3 + 3176704.d0*gl%z3_PCSAFT(1)**2 &
     & - 497664.d0*gl%z3_PCSAFT(1) + 34304.d0)
	part9 = (192.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 14880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 338400.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 4144800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4) + 33448992.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 197233464.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 905668080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 3375203880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 10466420640.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 27326586000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 60198351360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - &
     & 111443220960.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 172053489600.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 219210585000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) &
     & + 227279145072.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 188030691048.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 120578018400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) - 57280613040.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 18684187200.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 3638493648.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4) + 331033824.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 31687200.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 5339520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + &
     & 2799360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 248832.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) - 96.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 21840.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 571680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + &
     & 7371600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 60869424.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 21*gl%z3_PCSAFT(4) + 362638884.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 1672674360.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 6252506460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) - 19488518640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 51406987128.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 115287877968.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) + 219420284760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - &
     & 352580495280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 475025446020.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 532062358968.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) &
     & + 490120074204.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 365938014960.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 216897960240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) &
     & - 98879182560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 32860946016.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 7131124512.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 713665440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 42128640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4) - 5149440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 2944512.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 248832.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4) - 144.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 10440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + &
     & 351000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 4869000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 &
     & *gl%z3_PCSAFT(4) + 41499072.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 251110836.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 1167970320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 4392798660.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 13797897120.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 36849156912.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 84257651160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - &
     & 164939983920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 275412556800.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 390223898100.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) &
     & + 466262553600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 466292247444.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 386453962800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) - 261800449920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + &
     & 142106120640.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 59971499856.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 18749503296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 3983431680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 473328000.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 11623680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - &
     & 2271744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 41472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) + 120.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 1500.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - &
     & 91440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 1412700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4) - 12571176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 77740440.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 366205080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + &
     & 1390328160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 4411117560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 17*gl%z3_PCSAFT(4) + 11943375564.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 27850707936.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 56007931260.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4) - 96862530840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 143424237120.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 180886943928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4) + 193157927496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 173342256000.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 129451863360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) &
     & - 79374582240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 39214371648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **6*gl%z3_PCSAFT(4) - 15191025024.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 4427176320.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 905917440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + &
     & 113556480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 5803008.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) - 170496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) - 24.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 60.d0* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 8280.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 151500.d0*gl%z3_PCSAFT(1)** &
     & 22*gl%z3_PCSAFT(4) + 1428936.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 9095280.d0*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) + 43590000.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 167577480.d0*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) + 537998040.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 1477360356.d0*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 3509530008.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 7231314780.d0*gl%z3_PCSAFT(1)**14 &
     & *gl%z3_PCSAFT(4) + 12897411000.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 19833976200.d0*gl%z3_PCSAFT(1) &
     & **12*gl%z3_PCSAFT(4) + 26182003872.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 29523739344.d0* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 28277205120.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - &
     & 22837430880.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 15403212480.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - &
     & 8562933120.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 3852722304.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 1366982400.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 367918080.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - &
     & 70552320.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 8583168.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 497664.d0* &
     & gl%z3_PCSAFT(4))
        do xj = 1, gl%ncomp
		    do xi = 1, gl%ncomp
                if (xi .GE. xj) then
            	    gl%cx2_PCSAFT(6,xi,xj) = part1*(gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*part2 + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(4,xi)*part3 + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(1,xi)*part4 + &
         & gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(4,xj)*part5 + gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(1,xj)*part6 + gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*part7 + &
         & gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*part8 + gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part9)
                end if
            end do
        end do
    end if

    !DEC$ END IF
end subroutine CX2DERIVS



    end module pc_saft_CX2_derivs_module
