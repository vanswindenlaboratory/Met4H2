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

    ! module for file pc-saft_CX3_derivs.f90
    module pc_saft_CX3_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module


    contains



subroutine CX3DERIVS(gl,GETDERC)
    
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
    !output: cx3_PCSAFT (module variable)
    !working variables
    double precision :: recurring_factor
    double precision :: part1, part2, part3, part4, part5, part6
    double precision :: part7, part8, part9, part10, part11, part12, part13, part14, part15, part16
    double precision :: part17, part18, part19, part20, part21
    integer :: i, xi, xj, xk
    integer:: errorfld
    
!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE

    ! III. recurring_factor occurs in all derivatives
    recurring_factor = 2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 16.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 51.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 70.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 27.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 12.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - gl%z3_PCSAFT(1)** &
            & 6 + 8.d0*gl%z3_PCSAFT(1)**5 - 27.d0*gl%z3_PCSAFT(1)**4 + 42.d0*gl%z3_PCSAFT(1)**3 - 26.d0*gl%z3_PCSAFT(1)**2 + 4.d0
    
    ! IV. X3 DERIVATIVES of C
    ! 1: d^3(C_1)/(dxi dxj dxk)
    ! requires z3_PCSAFT(1)
    if (GETDERC(1) .eq. 1) then
        part1 = -2.d0*(gl%z3_PCSAFT(1) - 1.d0)/recurring_factor**4
        part2 = (24.d0*gl%z3_PCSAFT(1)**23 - 744.d0*gl%z3_PCSAFT(1)**22 + &
     & 10932.d0*gl%z3_PCSAFT(1)**21 - 100836.d0*gl%z3_PCSAFT(1)**20 + 652122.d0*gl%z3_PCSAFT(1)**19 - &
     & 3128670.d0*gl%z3_PCSAFT(1)**18 + 11492103.d0*gl%z3_PCSAFT(1)**17 - 32876319.d0*gl%z3_PCSAFT(1)**16 + &
     & 73763850.d0*gl%z3_PCSAFT(1)**15 - 129556686.d0*gl%z3_PCSAFT(1)**14 + 176018592.d0*gl%z3_PCSAFT(1)**13 &
     & - 180078852.d0*gl%z3_PCSAFT(1)**12 + 131009160.d0*gl%z3_PCSAFT(1)**11 - 58311624.d0*gl%z3_PCSAFT(1)** &
     & 10 + 6028005.d0*gl%z3_PCSAFT(1)**9 + 9501975.d0*gl%z3_PCSAFT(1)**8 - 4987764.d0*gl%z3_PCSAFT(1)**7 + &
     & 164268.d0*gl%z3_PCSAFT(1)**6 + 478224.d0*gl%z3_PCSAFT(1)**5 - 57024.d0*gl%z3_PCSAFT(1)**4 - 20736.d0* &
     & gl%z3_PCSAFT(1)**3)
        part3 = (8.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 368.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19 + 6552.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 66320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + &
     & 442018.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 2085220.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 7236292.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 18829460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 36923072.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12 - 54098356.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 57660540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 10 - 42021172.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 17650566.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - &
     & 1056408.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 2766360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 974520.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5 + 106128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 69120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 6912.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 32.d0*gl%z3_PCSAFT(1)**20 - 752.d0*gl%z3_PCSAFT(1)**19 + 8220.d0*gl%z3_PCSAFT(1)** &
     & 18 - 55072.d0*gl%z3_PCSAFT(1)**17 + 250734.d0*gl%z3_PCSAFT(1)**16 - 810928.d0*gl%z3_PCSAFT(1)**15 + &
     & 1878768.d0*gl%z3_PCSAFT(1)**14 - 3004128.d0*gl%z3_PCSAFT(1)**13 + 2771940.d0*gl%z3_PCSAFT(1)**12 + &
     & 338536.d0*gl%z3_PCSAFT(1)**11 - 5894116.d0*gl%z3_PCSAFT(1)**10 + 10693032.d0*gl%z3_PCSAFT(1)**9 - &
     & 11311682.d0*gl%z3_PCSAFT(1)**8 + 7678296.d0*gl%z3_PCSAFT(1)**7 - 3045992.d0*gl%z3_PCSAFT(1)**6 + &
     & 334008.d0*gl%z3_PCSAFT(1)**5 + 276048.d0*gl%z3_PCSAFT(1)**4 - 112704.d0*gl%z3_PCSAFT(1)**3 + 1152.d0* &
     & gl%z3_PCSAFT(1)**2 + 4608.d0*gl%z3_PCSAFT(1))
	part4 = (8.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 &
     & - 368.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 6552.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 66320.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**17 + 442018.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 2085220.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + &
     & 7236292.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 18829460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 36923072.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 54098356.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 57660540.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10 - 42021172.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 17650566.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 &
     & - 1056408.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 2766360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 974520.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 106128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 69120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 &
     & - 6912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 32.d0*gl%z3_PCSAFT(1)**20 - 752.d0*gl%z3_PCSAFT(1)**19 + 8220.d0* &
     & gl%z3_PCSAFT(1)**18 - 55072.d0*gl%z3_PCSAFT(1)**17 + 250734.d0*gl%z3_PCSAFT(1)**16 - 810928.d0*gl%z3_PCSAFT(1) &
     & **15 + 1878768.d0*gl%z3_PCSAFT(1)**14 - 3004128.d0*gl%z3_PCSAFT(1)**13 + 2771940.d0*gl%z3_PCSAFT(1)** &
     & 12 + 338536.d0*gl%z3_PCSAFT(1)**11 - 5894116.d0*gl%z3_PCSAFT(1)**10 + 10693032.d0*gl%z3_PCSAFT(1)**9 &
     & - 11311682.d0*gl%z3_PCSAFT(1)**8 + 7678296.d0*gl%z3_PCSAFT(1)**7 - 3045992.d0*gl%z3_PCSAFT(1)**6 + &
     & 334008.d0*gl%z3_PCSAFT(1)**5 + 276048.d0*gl%z3_PCSAFT(1)**4 - 112704.d0*gl%z3_PCSAFT(1)**3 + 1152.d0* &
     & gl%z3_PCSAFT(1)**2 + 4608.d0*gl%z3_PCSAFT(1))
	part5 = (12.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**19 - 612.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 10880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 &
     & - 106712.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 681631.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - &
     & 3075481.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 10251739.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - &
     & 25885565.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 50074577.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - &
     & 74150471.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 82775761.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - &
     & 67043103.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 36033036.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - &
     & 9665124.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 1146708.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 1439964.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 117792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 69120.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**2 - 6912.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 24.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 264.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 33272.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 &
     & - 305506.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 1642014.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 6079314.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 16549086.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 34220286.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11 + 54634434.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 67724910.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 &
     & + 64948994.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 47505656.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 25653736.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 9448152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 1788168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **4 + 199872.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 192384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 23040.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 4608.d0*gl%mmean_PCSAFT - 15.d0*gl%z3_PCSAFT(1)**19 + 285.d0*gl%z3_PCSAFT(1)**18 - 2305.d0 &
     & *gl%z3_PCSAFT(1)**17 + 9475.d0*gl%z3_PCSAFT(1)**16 - 12047.d0*gl%z3_PCSAFT(1)**15 - 86543.d0*gl%z3_PCSAFT(1)** &
     & 14 + 621905.d0*gl%z3_PCSAFT(1)**13 - 2244223.d0*gl%z3_PCSAFT(1)**12 + 5570442.d0*gl%z3_PCSAFT(1)**11 &
     & - 10346082.d0*gl%z3_PCSAFT(1)**10 + 14876292.d0*gl%z3_PCSAFT(1)**9 - 16772000.d0*gl%z3_PCSAFT(1)**8 + &
     & 14806880.d0*gl%z3_PCSAFT(1)**7 - 10083616.d0*gl%z3_PCSAFT(1)**6 + 5123072.d0*gl%z3_PCSAFT(1)**5 - &
     & 1813904.d0*gl%z3_PCSAFT(1)**4 + 377440.d0*gl%z3_PCSAFT(1)**3 - 15072.d0*gl%z3_PCSAFT(1)**2 - 12096.d0* &
     & gl%z3_PCSAFT(1) + 2112.d0)
	part6 = (8.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 368.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 6552.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 66320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 &
     & + 442018.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 2085220.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 7236292.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 18829460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 36923072.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12 - 54098356.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 57660540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 10 - 42021172.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 17650566.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - &
     & 1056408.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 2766360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 974520.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5 + 106128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 69120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 6912.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 32.d0*gl%z3_PCSAFT(1)**20 - 752.d0*gl%z3_PCSAFT(1)**19 + 8220.d0*gl%z3_PCSAFT(1)** &
     & 18 - 55072.d0*gl%z3_PCSAFT(1)**17 + 250734.d0*gl%z3_PCSAFT(1)**16 - 810928.d0*gl%z3_PCSAFT(1)**15 + &
     & 1878768.d0*gl%z3_PCSAFT(1)**14 - 3004128.d0*gl%z3_PCSAFT(1)**13 + 2771940.d0*gl%z3_PCSAFT(1)**12 + &
     & 338536.d0*gl%z3_PCSAFT(1)**11 - 5894116.d0*gl%z3_PCSAFT(1)**10 + 10693032.d0*gl%z3_PCSAFT(1)**9 - &
     & 11311682.d0*gl%z3_PCSAFT(1)**8 + 7678296.d0*gl%z3_PCSAFT(1)**7 - 3045992.d0*gl%z3_PCSAFT(1)**6 + &
     & 334008.d0*gl%z3_PCSAFT(1)**5 + 276048.d0*gl%z3_PCSAFT(1)**4 - 112704.d0*gl%z3_PCSAFT(1)**3 + 1152.d0* &
     & gl%z3_PCSAFT(1)**2 + 4608.d0*gl%z3_PCSAFT(1))
	part7 = (12.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**19 - 612.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 10880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 &
     & - 106712.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 681631.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - &
     & 3075481.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 10251739.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - &
     & 25885565.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 50074577.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - &
     & 74150471.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 82775761.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - &
     & 67043103.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 36033036.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - &
     & 9665124.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 1146708.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 1439964.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 117792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 69120.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**2 - 6912.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 24.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 264.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 33272.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 &
     & - 305506.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 1642014.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 6079314.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 16549086.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 34220286.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11 + 54634434.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 67724910.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 &
     & + 64948994.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 47505656.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 25653736.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 9448152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 1788168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **4 + 199872.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 192384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 23040.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 4608.d0*gl%mmean_PCSAFT - 15.d0*gl%z3_PCSAFT(1)**19 + 285.d0*gl%z3_PCSAFT(1)**18 - 2305.d0 &
     & *gl%z3_PCSAFT(1)**17 + 9475.d0*gl%z3_PCSAFT(1)**16 - 12047.d0*gl%z3_PCSAFT(1)**15 - 86543.d0*gl%z3_PCSAFT(1)** &
     & 14 + 621905.d0*gl%z3_PCSAFT(1)**13 - 2244223.d0*gl%z3_PCSAFT(1)**12 + 5570442.d0*gl%z3_PCSAFT(1)**11 &
     & - 10346082.d0*gl%z3_PCSAFT(1)**10 + 14876292.d0*gl%z3_PCSAFT(1)**9 - 16772000.d0*gl%z3_PCSAFT(1)**8 + &
     & 14806880.d0*gl%z3_PCSAFT(1)**7 - 10083616.d0*gl%z3_PCSAFT(1)**6 + 5123072.d0*gl%z3_PCSAFT(1)**5 - &
     & 1813904.d0*gl%z3_PCSAFT(1)**4 + 377440.d0*gl%z3_PCSAFT(1)**3 - 15072.d0*gl%z3_PCSAFT(1)**2 - 12096.d0* &
     & gl%z3_PCSAFT(1) + 2112.d0)
	part8 = (12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 19 - 612.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 10880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 106712.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 681631.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 3075481.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**14 + 10251739.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 25885565.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**12 + 50074577.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 74150471.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**10 + 82775761.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 67043103.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8 + 36033036.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 9665124.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **6 - 1146708.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 1439964.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - &
     & 117792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 69120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 6912.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1) + 24.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 776.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 33272.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 305506.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 15 + 1642014.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 6079314.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + &
     & 16549086.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 34220286.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 54634434.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 67724910.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 64948994.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**8 - 47505656.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 25653736.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - &
     & 9448152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 1788168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 199872.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**3 - 192384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 23040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 4608.d0* &
     & gl%mmean_PCSAFT - 15.d0*gl%z3_PCSAFT(1)**19 + 285.d0*gl%z3_PCSAFT(1)**18 - 2305.d0*gl%z3_PCSAFT(1)**17 + 9475.d0* &
     & gl%z3_PCSAFT(1)**16 - 12047.d0*gl%z3_PCSAFT(1)**15 - 86543.d0*gl%z3_PCSAFT(1)**14 + 621905.d0*gl%z3_PCSAFT(1) &
     & **13 - 2244223.d0*gl%z3_PCSAFT(1)**12 + 5570442.d0*gl%z3_PCSAFT(1)**11 - 10346082.d0*gl%z3_PCSAFT(1) &
     & **10 + 14876292.d0*gl%z3_PCSAFT(1)**9 - 16772000.d0*gl%z3_PCSAFT(1)**8 + 14806880.d0*gl%z3_PCSAFT(1) &
     & **7 - 10083616.d0*gl%z3_PCSAFT(1)**6 + 5123072.d0*gl%z3_PCSAFT(1)**5 - 1813904.d0*gl%z3_PCSAFT(1)**4 &
     & + 377440.d0*gl%z3_PCSAFT(1)**3 - 15072.d0*gl%z3_PCSAFT(1)**2 - 12096.d0*gl%z3_PCSAFT(1) + 2112.d0)
	part9 = (48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 2688.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + 46920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 440688.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**15 + 2685072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 11639472.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**13 + 37820544.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 94869936.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11 + 186177984.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 285790752.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**9 + 338571768.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - 300527040.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7 + 188963616.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 74766816.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**5 + 12735072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 1615248.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **3 - 353376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 207360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) - 20736.d0* &
     & gl%mmean_PCSAFT**3 + 2400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 52764.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + &
     & 540744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 3450312.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + &
     & 15426216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 51451488.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + &
     & 132643272.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 269219784.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + &
     & 432649560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 547836564.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + &
     & 538435344.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 399847440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + &
     & 214506480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 76749456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + &
     & 15374064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 810720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 156096.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 3456.d0*gl%mmean_PCSAFT**2 - 36.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 384.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**17 + 17628.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 211272.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + &
     & 1446000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 6750024.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 23293104.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 62075208.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 130835460.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10 - 220443816.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 297205764.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 8 - 318374160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 266980560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - &
     & 171031920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 80498256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 26049984.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 5082240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 408000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - &
     & 14208.d0*gl%mmean_PCSAFT + 12.d0*gl%z3_PCSAFT(1)**18 - 72.d0*gl%z3_PCSAFT(1)**17 - 1488.d0*gl%z3_PCSAFT(1)**16 + &
     & 25536.d0*gl%z3_PCSAFT(1)**15 - 196056.d0*gl%z3_PCSAFT(1)**14 + 973392.d0*gl%z3_PCSAFT(1)**13 - &
     & 3511248.d0*gl%z3_PCSAFT(1)**12 + 9735552.d0*gl%z3_PCSAFT(1)**11 - 21385524.d0*gl%z3_PCSAFT(1)**10 + &
     & 37787352.d0*gl%z3_PCSAFT(1)**9 - 53984160.d0*gl%z3_PCSAFT(1)**8 + 62208192.d0*gl%z3_PCSAFT(1)**7 - &
     & 57318240.d0*gl%z3_PCSAFT(1)**6 + 41585760.d0*gl%z3_PCSAFT(1)**5 - 23198496.d0*gl%z3_PCSAFT(1)**4 + &
     & 9589056.d0*gl%z3_PCSAFT(1)**3 - 2763072.d0*gl%z3_PCSAFT(1)**2 + 494976.d0*gl%z3_PCSAFT(1) - 41472.d0)
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
            	        gl%cx3_PCSAFT(1,xi,xj,xk) = part1*(gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*gl%mPCSAFT(xk)*part2 + gl%mPCSAFT(xj)*gl%mPCSAFT(xi)* &
	    & gl%z3x1_PCSAFT(1,xk)*part3 + gl%mPCSAFT(xj)*gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(1,xi)*part4 + &
	    & gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part5 + gl%mPCSAFT(xi)*gl%mPCSAFT(xk)*&
	    & gl%z3x1_PCSAFT(1,xj)*part6 + gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*part7 + &
	    & gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part8 + &
	    & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part9)
                    end if
        	    end do
            end do 
        end do
    end if
    
    ! 2: 1ST DERIVATIVE OF C_1 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires z3_PCSAFT(1)
    if (GETDERC(2) .eq. 1) then
        part1 = 4.d0/recurring_factor**5
        part2 = (24.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**28 - 1320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 + &
     & 30396.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 414300.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 3843126.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**24 - 26098218.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 135431157.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 22 - 551449497.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 1790567295.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - &
     & 4677116895.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 9853794570.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - &
     & 16686289482.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 22466975796.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - &
     & 23531399940.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 18363154821.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - &
     & 9669803721.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 2363469219.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
     & 802496781.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 860810400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
     & 197028936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 55160028.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 28731672.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 1360800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 1399680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **5 + 124416.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 132.d0*gl%z3_PCSAFT(1)**28 - 4380.d0*gl%z3_PCSAFT(1)**27 + &
     & 69312.d0*gl%z3_PCSAFT(1)**26 - 694080.d0*gl%z3_PCSAFT(1)**25 + 4922277.d0*gl%z3_PCSAFT(1)**24 - &
     & 26216823.d0*gl%z3_PCSAFT(1)**23 + 108503259.d0*gl%z3_PCSAFT(1)**22 - 355964979.d0*gl%z3_PCSAFT(1)**21 &
     & + 934979370.d0*gl%z3_PCSAFT(1)**20 - 1968039624.d0*gl%z3_PCSAFT(1)**19 + 3285634506.d0* &
     & gl%z3_PCSAFT(1)**18 - 4220456142.d0*gl%z3_PCSAFT(1)**17 + 3835754325.d0*gl%z3_PCSAFT(1)**16 - &
     & 1720034475.d0*gl%z3_PCSAFT(1)**15 - 1298004981.d0*gl%z3_PCSAFT(1)**14 + 3509697477.d0*gl%z3_PCSAFT(1) &
     & **13 - 3744449952.d0*gl%z3_PCSAFT(1)**12 + 2389775286.d0*gl%z3_PCSAFT(1)**11 - 820691040.d0* &
     & gl%z3_PCSAFT(1)**10 - 8360388.d0*gl%z3_PCSAFT(1)**9 + 140522904.d0*gl%z3_PCSAFT(1)**8 - 48796560.d0* &
     & gl%z3_PCSAFT(1)**7 - 1735344.d0*gl%z3_PCSAFT(1)**6 + 4004640.d0*gl%z3_PCSAFT(1)**5 - 290304.d0*gl%z3_PCSAFT(1) &
     & **4 - 124416.d0*gl%z3_PCSAFT(1)**3)
        part3 = (16.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **27 - 1000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 23448.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - &
     & 313620.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 2807212.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - &
     & 18246642.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 90364242.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - &
     & 351527431.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 1095821600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - &
     & 2772663215.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 5735183160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - &
     & 9716828902.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 13430529548.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & - 14964967460.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 13117667042.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13 - 8615783839.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 3786214464.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **11 - 703782987.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 292610820.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **9 + 216928512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 29751048.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 &
     & - 11803752.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 2233440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + &
     & 466560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 41472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 56.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**27 - 1340.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 + 13032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - &
     & 50180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 197282.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 3898221.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22 - 27575112.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 127679726.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 20 - 437212430.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 1165970669.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - &
     & 2495515612.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 4383387516.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - &
     & 6448195750.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 8099045975.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - &
     & 8807773344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 8283537726.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - &
     & 6559070234.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 4116709571.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - &
     & 1821598180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 400248348.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 81676152.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 81419160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 14239296.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5 + 3149280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 850176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - &
     & 96768.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 32.d0*gl%z3_PCSAFT(1)**27 + 920.d0*gl%z3_PCSAFT(1)**26 - 12216.d0* &
     & gl%z3_PCSAFT(1)**25 + 98590.d0*gl%z3_PCSAFT(1)**24 - 532098.d0*gl%z3_PCSAFT(1)**23 + 1959147.d0* &
     & gl%z3_PCSAFT(1)**22 - 4539818.d0*gl%z3_PCSAFT(1)**21 + 3126347.d0*gl%z3_PCSAFT(1)**20 + 24332150.d0* &
     & gl%z3_PCSAFT(1)**19 - 129815372.d0*gl%z3_PCSAFT(1)**18 + 389648580.d0*gl%z3_PCSAFT(1)**17 - &
     & 859550348.d0*gl%z3_PCSAFT(1)**16 + 1503298050.d0*gl%z3_PCSAFT(1)**15 - 2153930265.d0*gl%z3_PCSAFT(1) &
     & **14 + 2566237470.d0*gl%z3_PCSAFT(1)**13 - 2552079457.d0*gl%z3_PCSAFT(1)**12 + 2108334442.d0 &
     & *gl%z3_PCSAFT(1)**11 - 1429110350.d0*gl%z3_PCSAFT(1)**10 + 778355200.d0*gl%z3_PCSAFT(1)**9 - &
     & 328082396.d0*gl%z3_PCSAFT(1)**8 + 97790048.d0*gl%z3_PCSAFT(1)**7 - 13969168.d0*gl%z3_PCSAFT(1)**6 - &
     & 3700896.d0*gl%z3_PCSAFT(1)**5 + 2695200.d0*gl%z3_PCSAFT(1)**4 - 551424.d0*gl%z3_PCSAFT(1)**3 - 20736.d0* &
     & gl%z3_PCSAFT(1)**2 + 18432.d0*gl%z3_PCSAFT(1))
	part4 = (16.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**27 - 1000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 23448.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 &
     & - 313620.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 2807212.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - &
     & 18246642.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 90364242.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - &
     & 351527431.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 1095821600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - &
     & 2772663215.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 5735183160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - &
     & 9716828902.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 13430529548.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & - 14964967460.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 13117667042.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13 - 8615783839.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 3786214464.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **11 - 703782987.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 292610820.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **9 + 216928512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 29751048.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 &
     & - 11803752.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 2233440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + &
     & 466560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 41472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 56.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**27 - 1340.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 + 13032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - &
     & 50180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 197282.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 3898221.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22 - 27575112.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 127679726.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 20 - 437212430.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 1165970669.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - &
     & 2495515612.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 4383387516.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - &
     & 6448195750.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 8099045975.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - &
     & 8807773344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 8283537726.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - &
     & 6559070234.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 4116709571.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - &
     & 1821598180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 400248348.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 81676152.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 81419160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 14239296.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5 + 3149280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 850176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - &
     & 96768.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 32.d0*gl%z3_PCSAFT(1)**27 + 920.d0*gl%z3_PCSAFT(1)**26 - 12216.d0* &
     & gl%z3_PCSAFT(1)**25 + 98590.d0*gl%z3_PCSAFT(1)**24 - 532098.d0*gl%z3_PCSAFT(1)**23 + 1959147.d0* &
     & gl%z3_PCSAFT(1)**22 - 4539818.d0*gl%z3_PCSAFT(1)**21 + 3126347.d0*gl%z3_PCSAFT(1)**20 + 24332150.d0* &
     & gl%z3_PCSAFT(1)**19 - 129815372.d0*gl%z3_PCSAFT(1)**18 + 389648580.d0*gl%z3_PCSAFT(1)**17 - &
     & 859550348.d0*gl%z3_PCSAFT(1)**16 + 1503298050.d0*gl%z3_PCSAFT(1)**15 - 2153930265.d0*gl%z3_PCSAFT(1) &
     & **14 + 2566237470.d0*gl%z3_PCSAFT(1)**13 - 2552079457.d0*gl%z3_PCSAFT(1)**12 + 2108334442.d0 &
     & *gl%z3_PCSAFT(1)**11 - 1429110350.d0*gl%z3_PCSAFT(1)**10 + 778355200.d0*gl%z3_PCSAFT(1)**9 - &
     & 328082396.d0*gl%z3_PCSAFT(1)**8 + 97790048.d0*gl%z3_PCSAFT(1)**7 - 13969168.d0*gl%z3_PCSAFT(1)**6 - &
     & 3700896.d0*gl%z3_PCSAFT(1)**5 + 2695200.d0*gl%z3_PCSAFT(1)**4 - 551424.d0*gl%z3_PCSAFT(1)**3 - 20736.d0* &
     & gl%z3_PCSAFT(1)**2 + 18432.d0*gl%z3_PCSAFT(1))
	part5 = (24.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**26 - 1680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + 39140.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 24 - 504200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + 4299782.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - &
     & 26591508.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + 125974859.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - &
     & 473989210.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + 1451755955.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - &
     & 3677984500.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + 7772982770.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - &
     & 13721514924.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + 20120913700.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 &
     & - 24205886460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + 23406178307.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 12 - 17611040946.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 + 9748346703.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **10 - 3517359420.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 516452940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **8 + 163441872.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - 82448532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 &
     & + 5537376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 889920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 466560.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 41472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 + 36.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**26 + 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - 27540.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 &
     & + 464700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - 4442793.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + &
     & 29011056.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - 140810061.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + &
     & 533937690.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 1635368130.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + &
     & 4141738416.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 8811917442.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + &
     & 15892338288.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 24360661905.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 &
     & + 31646200440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 34580504781.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 12 + 31413017826.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 23312412552.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**10 + 13739660160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 6092100960.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**8 + 1779032472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 178955808.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6 - 91465920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 34912080.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**4 - 567360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 981504.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 &
     & - 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 54.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 + 900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **25 - 1545.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 94200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 1321821.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 9964776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 51850752.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**20 - 204096690.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 639504840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **18 - 1650721836.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 3597346575.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 &
     & - 6727619628.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 10884531765.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - &
     & 15235004760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 18324077130.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - &
     & 18713161554.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 15964586892.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - &
     & 11141421960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 6184276320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - &
     & 2619809928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 788292096.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - &
     & 140928288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 2531040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 5851200.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 1449792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 85248.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + &
     & 18432.d0*gl%mmean_PCSAFT + 15.d0*gl%z3_PCSAFT(1)**26 - 330.d0*gl%z3_PCSAFT(1)**25 + 2765.d0*gl%z3_PCSAFT(1)**24 - &
     & 6050.d0*gl%z3_PCSAFT(1)**23 - 87604.d0*gl%z3_PCSAFT(1)**22 + 1051470.d0*gl%z3_PCSAFT(1)**21 - 6443878.d0 &
     & *gl%z3_PCSAFT(1)**20 + 27610760.d0*gl%z3_PCSAFT(1)**19 - 91268065.d0*gl%z3_PCSAFT(1)**18 + &
     & 245112266.d0*gl%z3_PCSAFT(1)**17 - 552973423.d0*gl%z3_PCSAFT(1)**16 + 1070948374.d0*gl%z3_PCSAFT(1)** &
     & 15 - 1801278330.d0*gl%z3_PCSAFT(1)**14 + 2637503410.d0*gl%z3_PCSAFT(1)**13 - 3346585912.d0* &
     & gl%z3_PCSAFT(1)**12 + 3646896388.d0*gl%z3_PCSAFT(1)**11 - 3375154496.d0*gl%z3_PCSAFT(1)**10 + &
     & 2619334480.d0*gl%z3_PCSAFT(1)**9 - 1679372000.d0*gl%z3_PCSAFT(1)**8 + 872648272.d0*gl%z3_PCSAFT(1)**7 &
     & - 357473952.d0*gl%z3_PCSAFT(1)**6 + 110280032.d0*gl%z3_PCSAFT(1)**5 - 23366400.d0*gl%z3_PCSAFT(1)**4 &
     & + 2558080.d0*gl%z3_PCSAFT(1)**3 + 140928.d0*gl%z3_PCSAFT(1)**2 - 85248.d0*gl%z3_PCSAFT(1) + 8448.d0)
	part6 = (16.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27 - 1000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **26 + 23448.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - 313620.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + &
     & 2807212.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - 18246642.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + &
     & 90364242.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - 351527431.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + &
     & 1095821600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 2772663215.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + &
     & 5735183160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 9716828902.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + &
     & 13430529548.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 14964967460.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 &
     & + 13117667042.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 8615783839.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 12 + 3786214464.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 703782987.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 10 - 292610820.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 216928512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 &
     & - 29751048.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 11803752.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + &
     & 2233440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 466560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 41472.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 56.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 - 1340.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 &
     & + 13032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - 50180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 197282.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**23 + 3898221.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 27575112.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 &
     & + 127679726.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 437212430.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + &
     & 1165970669.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 2495515612.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + &
     & 4383387516.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 6448195750.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + &
     & 8099045975.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 8807773344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + &
     & 8283537726.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 6559070234.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + &
     & 4116709571.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 1821598180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + &
     & 400248348.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 81676152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 81419160.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 14239296.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 3149280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **4 - 850176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 96768.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 32.d0*gl%z3_PCSAFT(1) &
     & **27 + 920.d0*gl%z3_PCSAFT(1)**26 - 12216.d0*gl%z3_PCSAFT(1)**25 + 98590.d0*gl%z3_PCSAFT(1)**24 - &
     & 532098.d0*gl%z3_PCSAFT(1)**23 + 1959147.d0*gl%z3_PCSAFT(1)**22 - 4539818.d0*gl%z3_PCSAFT(1)**21 + &
     & 3126347.d0*gl%z3_PCSAFT(1)**20 + 24332150.d0*gl%z3_PCSAFT(1)**19 - 129815372.d0*gl%z3_PCSAFT(1)**18 + &
     & 389648580.d0*gl%z3_PCSAFT(1)**17 - 859550348.d0*gl%z3_PCSAFT(1)**16 + 1503298050.d0*gl%z3_PCSAFT(1)** &
     & 15 - 2153930265.d0*gl%z3_PCSAFT(1)**14 + 2566237470.d0*gl%z3_PCSAFT(1)**13 - 2552079457.d0* &
     & gl%z3_PCSAFT(1)**12 + 2108334442.d0*gl%z3_PCSAFT(1)**11 - 1429110350.d0*gl%z3_PCSAFT(1)**10 + &
     & 778355200.d0*gl%z3_PCSAFT(1)**9 - 328082396.d0*gl%z3_PCSAFT(1)**8 + 97790048.d0*gl%z3_PCSAFT(1)**7 - &
     & 13969168.d0*gl%z3_PCSAFT(1)**6 - 3700896.d0*gl%z3_PCSAFT(1)**5 + 2695200.d0*gl%z3_PCSAFT(1)**4 - &
     & 551424.d0*gl%z3_PCSAFT(1)**3 - 20736.d0*gl%z3_PCSAFT(1)**2 + 18432.d0*gl%z3_PCSAFT(1))
	part7 = (24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - 1680.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**25 + 39140.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 504200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 23 + 4299782.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 26591508.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + &
     & 125974859.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 473989210.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + &
     & 1451755955.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 3677984500.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + &
     & 7772982770.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 13721514924.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & + 20120913700.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 24205886460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 13 + 23406178307.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 17611040946.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11 + 9748346703.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 3517359420.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**9 + 516452940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 163441872.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7 - 82448532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 5537376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **5 + 889920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 466560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + &
     & 41472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 + 36.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 360.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**25 - 27540.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 464700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 23 - 4442793.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 29011056.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - &
     & 140810061.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 533937690.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - &
     & 1635368130.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 4141738416.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - &
     & 8811917442.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 15892338288.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & - 24360661905.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 31646200440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13 - 34580504781.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 31413017826.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**11 - 23312412552.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 13739660160.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**9 - 6092100960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 1779032472.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**7 - 178955808.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 91465920.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**5 + 34912080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 567360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **3 - 981504.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 54.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 + 900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - 1545.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - &
     & 94200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 1321821.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 9964776.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**21 + 51850752.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 204096690.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **19 + 639504840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 1650721836.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + &
     & 3597346575.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 6727619628.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + &
     & 10884531765.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 15235004760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + &
     & 18324077130.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 18713161554.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + &
     & 15964586892.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 11141421960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + &
     & 6184276320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 2619809928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + &
     & 788292096.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 140928288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 2531040.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 5851200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 1449792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **2 + 85248.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 18432.d0*gl%mmean_PCSAFT + 15.d0*gl%z3_PCSAFT(1)**26 - 330.d0* &
     & gl%z3_PCSAFT(1)**25 + 2765.d0*gl%z3_PCSAFT(1)**24 - 6050.d0*gl%z3_PCSAFT(1)**23 - 87604.d0*gl%z3_PCSAFT(1)**22 &
     & + 1051470.d0*gl%z3_PCSAFT(1)**21 - 6443878.d0*gl%z3_PCSAFT(1)**20 + 27610760.d0*gl%z3_PCSAFT(1)**19 - &
     & 91268065.d0*gl%z3_PCSAFT(1)**18 + 245112266.d0*gl%z3_PCSAFT(1)**17 - 552973423.d0*gl%z3_PCSAFT(1)**16 &
     & + 1070948374.d0*gl%z3_PCSAFT(1)**15 - 1801278330.d0*gl%z3_PCSAFT(1)**14 + 2637503410.d0* &
     & gl%z3_PCSAFT(1)**13 - 3346585912.d0*gl%z3_PCSAFT(1)**12 + 3646896388.d0*gl%z3_PCSAFT(1)**11 - &
     & 3375154496.d0*gl%z3_PCSAFT(1)**10 + 2619334480.d0*gl%z3_PCSAFT(1)**9 - 1679372000.d0*gl%z3_PCSAFT(1) &
     & **8 + 872648272.d0*gl%z3_PCSAFT(1)**7 - 357473952.d0*gl%z3_PCSAFT(1)**6 + 110280032.d0* &
     & gl%z3_PCSAFT(1)**5 - 23366400.d0*gl%z3_PCSAFT(1)**4 + 2558080.d0*gl%z3_PCSAFT(1)**3 + 140928.d0* &
     & gl%z3_PCSAFT(1)**2 - 85248.d0*gl%z3_PCSAFT(1) + 8448.d0)
	part8 = (24.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - 1680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + 39140.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**24 - 504200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + 4299782.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **22 - 26591508.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + 125974859.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 20 - 473989210.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + 1451755955.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 18 - 3677984500.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + 7772982770.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **16 - 13721514924.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + 20120913700.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**14 - 24205886460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + 23406178307.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**12 - 17611040946.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 + 9748346703.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 3517359420.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 516452940.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 163441872.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - 82448532.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 5537376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 889920.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**4 + 466560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 41472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 &
     & + 36.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - 27540.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**24 + 464700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - 4442793.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**22 + 29011056.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - 140810061.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**20 + 533937690.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 1635368130.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**18 + 4141738416.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 8811917442.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**16 + 15892338288.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 24360661905.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**14 + 31646200440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 34580504781.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 31413017826.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - &
     & 23312412552.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 13739660160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 &
     & - 6092100960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 1779032472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - &
     & 178955808.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 91465920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + &
     & 34912080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 567360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 981504.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 54.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 &
     & + 900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - 1545.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 94200.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**23 + 1321821.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 9964776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 &
     & + 51850752.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 204096690.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + &
     & 639504840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 1650721836.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + &
     & 3597346575.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 6727619628.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + &
     & 10884531765.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 15235004760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + &
     & 18324077130.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 18713161554.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + &
     & 15964586892.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 11141421960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + &
     & 6184276320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 2619809928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + &
     & 788292096.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 140928288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 2531040.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 5851200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 1449792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **2 + 85248.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 18432.d0*gl%mmean_PCSAFT + 15.d0*gl%z3_PCSAFT(1)**26 - 330.d0* &
     & gl%z3_PCSAFT(1)**25 + 2765.d0*gl%z3_PCSAFT(1)**24 - 6050.d0*gl%z3_PCSAFT(1)**23 - 87604.d0*gl%z3_PCSAFT(1)**22 &
     & + 1051470.d0*gl%z3_PCSAFT(1)**21 - 6443878.d0*gl%z3_PCSAFT(1)**20 + 27610760.d0*gl%z3_PCSAFT(1)**19 - &
     & 91268065.d0*gl%z3_PCSAFT(1)**18 + 245112266.d0*gl%z3_PCSAFT(1)**17 - 552973423.d0*gl%z3_PCSAFT(1)**16 &
     & + 1070948374.d0*gl%z3_PCSAFT(1)**15 - 1801278330.d0*gl%z3_PCSAFT(1)**14 + 2637503410.d0* &
     & gl%z3_PCSAFT(1)**13 - 3346585912.d0*gl%z3_PCSAFT(1)**12 + 3646896388.d0*gl%z3_PCSAFT(1)**11 - &
     & 3375154496.d0*gl%z3_PCSAFT(1)**10 + 2619334480.d0*gl%z3_PCSAFT(1)**9 - 1679372000.d0*gl%z3_PCSAFT(1) &
     & **8 + 872648272.d0*gl%z3_PCSAFT(1)**7 - 357473952.d0*gl%z3_PCSAFT(1)**6 + 110280032.d0* &
     & gl%z3_PCSAFT(1)**5 - 23366400.d0*gl%z3_PCSAFT(1)**4 + 2558080.d0*gl%z3_PCSAFT(1)**3 + 140928.d0* &
     & gl%z3_PCSAFT(1)**2 - 85248.d0*gl%z3_PCSAFT(1) + 8448.d0)
	part9 = (96.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25 - 7440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 + &
     & 169200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 - 2072400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 + &
     & 16724496.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 - 98616732.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 + &
     & 452834040.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 - 1687601940.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 + &
     & 5233210320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 - 13663293000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 &
     & + 30099175680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 - 55721610480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 14 + 86026744800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 - 109605292500.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**12 + 113639572536.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 - 94015345524.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**10 + 60289009200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9 - 28640306520.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 + 9342093600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 - 1819246824.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 + 165516912.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 - 15843600.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 + 2669760.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 + 1399680.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)**2 + 124416.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) - 48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + &
     & 10920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 285840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + 3685800.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 30434712.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + 181319442.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 836337180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + 3126253230.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 9744259320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + &
     & 25703493564.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 57643938984.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & + 109710142380.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 176290247640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **13 + 237512723010.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 266031179484.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11 + 245060037102.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 182969007480.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 108448980120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - &
     & 49439591280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 16430473008.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - &
     & 3565562256.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 356832720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + &
     & 21064320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 2574720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - &
     & 1472256.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) - 124416.d0*gl%mmean_PCSAFT**3 - 72.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 25 - 5220.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 175500.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - &
     & 2434500.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 20749536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - &
     & 125555418.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 583985160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - &
     & 2196399330.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 6898948560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - &
     & 18424578456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 42128825580.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & - 82469991960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 137706278400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **13 - 195111949050.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 233131276800.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**11 - 233146123722.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 193226981400.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 130900224960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + &
     & 71053060320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 29985749928.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + &
     & 9374751648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 1991715840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + &
     & 236664000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 5811840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - &
     & 1135872.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 20736.d0*gl%mmean_PCSAFT**2 + 60.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + &
     & 750.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 45720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 706350.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22 - 6285588.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 38870220.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 &
     & - 183102540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 695164080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - &
     & 2205558780.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 5971687782.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - &
     & 13925353968.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 28003965630.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - &
     & 48431265420.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 71712118560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - &
     & 90443471964.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 96578963748.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - &
     & 86671128000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 64725931680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - &
     & 39687291120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 19607185824.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - &
     & 7595512512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 2213588160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - &
     & 452958720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 56778240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 2901504.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 85248.d0*gl%mmean_PCSAFT - 12.d0*gl%z3_PCSAFT(1)**25 + 30.d0*gl%z3_PCSAFT(1)**24 + 4140.d0 &
     & *gl%z3_PCSAFT(1)**23 - 75750.d0*gl%z3_PCSAFT(1)**22 + 714468.d0*gl%z3_PCSAFT(1)**21 - 4547640.d0* &
     & gl%z3_PCSAFT(1)**20 + 21795000.d0*gl%z3_PCSAFT(1)**19 - 83788740.d0*gl%z3_PCSAFT(1)**18 + &
     & 268999020.d0*gl%z3_PCSAFT(1)**17 - 738680178.d0*gl%z3_PCSAFT(1)**16 + 1754765004.d0*gl%z3_PCSAFT(1)** &
     & 15 - 3615657390.d0*gl%z3_PCSAFT(1)**14 + 6448705500.d0*gl%z3_PCSAFT(1)**13 - 9916988100.d0* &
     & gl%z3_PCSAFT(1)**12 + 13091001936.d0*gl%z3_PCSAFT(1)**11 - 14761869672.d0*gl%z3_PCSAFT(1)**10 + &
     & 14138602560.d0*gl%z3_PCSAFT(1)**9 - 11418715440.d0*gl%z3_PCSAFT(1)**8 + 7701606240.d0*gl%z3_PCSAFT(1) &
     & **7 - 4281466560.d0*gl%z3_PCSAFT(1)**6 + 1926361152.d0*gl%z3_PCSAFT(1)**5 - 683491200.d0* &
     & gl%z3_PCSAFT(1)**4 + 183959040.d0*gl%z3_PCSAFT(1)**3 - 35276160.d0*gl%z3_PCSAFT(1)**2 + 4291584.d0* &
     & gl%z3_PCSAFT(1) - 248832.d0)
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then 
            	        gl%cx3_PCSAFT(2,xi,xj,xk) = part1*(gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*gl%mPCSAFT(xk)*part2 + gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(1,xk)*part3 + &
         & gl%mPCSAFT(xj)*gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(1,xi)*part4 + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part5 + &
         & gl%mPCSAFT(xi)*gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(1,xj)*part6 + gl%mPCSAFT(xi)* &
         & gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*part7 + gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)*part8 + gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
         & gl%z3x1_PCSAFT(1,xk)*part9)
                    end if
        	    end do
            end do 
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF C_1 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires z3_PCSAFT(1), z3_PCSAFT(4)
    if (GETDERC(4) .eq. 1) then
        part1 = 2.d0/recurring_factor**5
        part2 = (48.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(4) - 2640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 26*gl%z3_PCSAFT(4) + 60792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 828600.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 7686252.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - &
     & 52196436.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 270862314.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) - 1102898994.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 3581134590.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 9354233790.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) &
     & + 19707589140.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 33372578964.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 44933951592.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - &
     & 47062799880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 36726309642.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **13*gl%z3_PCSAFT(4) - 19339607442.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + &
     & 4726938438.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 1604993562.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 10*gl%z3_PCSAFT(4) - 1721620800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 394057872.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 110320056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - &
     & 57463344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 2721600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4) + 2799360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 248832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4) + 264.d0*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(4) - 8760.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) &
     & + 138624.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 1388160.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + &
     & 9844554.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 52433646.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + &
     & 217006518.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 711929958.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + &
     & 1869958740.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 3936079248.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + &
     & 6571269012.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 8440912284.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + &
     & 7671508650.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 3440068950.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - &
     & 2596009962.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 7019394954.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - &
     & 7488899904.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 4779550572.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - &
     & 1641382080.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 16720776.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + &
     & 281045808.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 97593120.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - &
     & 3470688.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 8009280.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 580608.d0* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 248832.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4))
        part3 = (-16.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27 + 880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 - &
     & 20264.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 276200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 2562084.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 17398812.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 90287438.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 367632998.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 1193711530.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 3118077930.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 6569196380.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 11124192988.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - &
     & 14977983864.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 15687599960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 &
     & - 12242103214.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 6446535814.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 12 - 1575646146.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 534997854.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 10 + 573873600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 131352624.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 &
     & - 36773352.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 19154448.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + &
     & 907200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 933120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 82944.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 56.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 + 1640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 &
     & - 21976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 173640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 843446.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**23 + 2031034.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 3879688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 &
     & - 62163308.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 326625910.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - &
     & 1138151678.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 2957603032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - &
     & 5965444904.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 9469367070.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - &
     & 11819551330.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 11434547288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - &
     & 8281748236.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 4138672946.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - &
     & 1079963698.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 167460800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + &
     & 258356184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 78703512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 3034320.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 6868512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 901440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 4 - 165888.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 27648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 32.d0*gl%z3_PCSAFT(1)** &
     & 27 - 1040.d0*gl%z3_PCSAFT(1)**26 + 16108.d0*gl%z3_PCSAFT(1)**25 - 157580.d0*gl%z3_PCSAFT(1)**24 + &
     & 1088146.d0*gl%z3_PCSAFT(1)**23 - 5614202.d0*gl%z3_PCSAFT(1)**22 + 22331162.d0*gl%z3_PCSAFT(1)**21 - &
     & 69551646.d0*gl%z3_PCSAFT(1)**20 + 169965900.d0*gl%z3_PCSAFT(1)**19 - 320797416.d0*gl%z3_PCSAFT(1)**18 &
     & + 442978920.d0*gl%z3_PCSAFT(1)**17 - 364557048.d0*gl%z3_PCSAFT(1)**16 - 81361550.d0*gl%z3_PCSAFT(1)** &
     & 15 + 860915150.d0*gl%z3_PCSAFT(1)**14 - 1638584830.d0*gl%z3_PCSAFT(1)**13 + 1967101162.d0* &
     & gl%z3_PCSAFT(1)**12 - 1661424768.d0*gl%z3_PCSAFT(1)**11 + 969145860.d0*gl%z3_PCSAFT(1)**10 - &
     & 331113680.d0*gl%z3_PCSAFT(1)**9 + 2950712.d0*gl%z3_PCSAFT(1)**8 + 60759632.d0*gl%z3_PCSAFT(1)**7 - &
     & 28279328.d0*gl%z3_PCSAFT(1)**6 + 2922528.d0*gl%z3_PCSAFT(1)**5 + 1838400.d0*gl%z3_PCSAFT(1)**4 - &
     & 575232.d0*gl%z3_PCSAFT(1)**3 - 13824.d0*gl%z3_PCSAFT(1)**2 + 18432.d0*gl%z3_PCSAFT(1))
	part4 = (48.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 2880.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 67160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 903440.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 8176508.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4) - 53892096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 271015922.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 1070687860.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 3385354730.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 8663404360.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 18039562700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) - 30557850792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 41839042960.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 45617534880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) + 38477437298.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 23678103492.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 9148075074.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) - 872568120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 1159095240.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 565209648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) - 22728744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 42761952.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 3559680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 1866240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 165888.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) + 168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 4320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(4) + 48040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 274000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4) + 448882.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 5765408.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 59029912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + &
     & 317522760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 1201050770.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 &
     & *gl%z3_PCSAFT(4) + 3470093016.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 7948634256.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 14732219936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) &
     & - 22365758570.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 28017643280.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 29050093976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + &
     & 24848823688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 17256813414.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **10*gl%z3_PCSAFT(4) + 9313382840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 3475735560.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 542140512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 242055816.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 159804000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4) + 21610080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 7200000.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 1534464.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 221184.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 96.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 2880.d0*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(4) - 40540.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 354760.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) &
     & - 2152342.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 9532496.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - &
     & 31410798.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 75804340.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - &
     & 121301600.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 61166672.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + &
     & 336318240.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 1354543648.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + &
     & 3087957650.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 5168775680.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 6771059770.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 7071260076.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + &
     & 5878093652.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 3827366560.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 1887824080.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 659115504.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 134820464.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 340992.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 10324320.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 3552000.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 527616.d0 &
     & *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 27648.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 18432.d0*gl%z3_PCSAFT(4))
	part5 = (-16.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27 + 880.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**26 - 20264.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 276200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 24 - 2562084.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 17398812.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - &
     & 90287438.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 367632998.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - &
     & 1193711530.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 3118077930.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - &
     & 6569196380.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 11124192988.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 &
     & - 14977983864.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 15687599960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 14 - 12242103214.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 6446535814.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **12 - 1575646146.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 534997854.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **10 + 573873600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 131352624.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 8 - 36773352.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 19154448.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + &
     & 907200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 933120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 82944.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 56.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 + 1640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 &
     & - 21976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 173640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 843446.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**23 + 2031034.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 3879688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 &
     & - 62163308.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 326625910.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - &
     & 1138151678.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 2957603032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - &
     & 5965444904.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 9469367070.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - &
     & 11819551330.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 11434547288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - &
     & 8281748236.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 4138672946.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - &
     & 1079963698.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 167460800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + &
     & 258356184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 78703512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 3034320.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 6868512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 901440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 4 - 165888.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 27648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 32.d0*gl%z3_PCSAFT(1)** &
     & 27 - 1040.d0*gl%z3_PCSAFT(1)**26 + 16108.d0*gl%z3_PCSAFT(1)**25 - 157580.d0*gl%z3_PCSAFT(1)**24 + &
     & 1088146.d0*gl%z3_PCSAFT(1)**23 - 5614202.d0*gl%z3_PCSAFT(1)**22 + 22331162.d0*gl%z3_PCSAFT(1)**21 - &
     & 69551646.d0*gl%z3_PCSAFT(1)**20 + 169965900.d0*gl%z3_PCSAFT(1)**19 - 320797416.d0*gl%z3_PCSAFT(1)**18 &
     & + 442978920.d0*gl%z3_PCSAFT(1)**17 - 364557048.d0*gl%z3_PCSAFT(1)**16 - 81361550.d0*gl%z3_PCSAFT(1)** &
     & 15 + 860915150.d0*gl%z3_PCSAFT(1)**14 - 1638584830.d0*gl%z3_PCSAFT(1)**13 + 1967101162.d0* &
     & gl%z3_PCSAFT(1)**12 - 1661424768.d0*gl%z3_PCSAFT(1)**11 + 969145860.d0*gl%z3_PCSAFT(1)**10 - &
     & 331113680.d0*gl%z3_PCSAFT(1)**9 + 2950712.d0*gl%z3_PCSAFT(1)**8 + 60759632.d0*gl%z3_PCSAFT(1)**7 - &
     & 28279328.d0*gl%z3_PCSAFT(1)**6 + 2922528.d0*gl%z3_PCSAFT(1)**5 + 1838400.d0*gl%z3_PCSAFT(1)**4 - &
     & 575232.d0*gl%z3_PCSAFT(1)**3 - 13824.d0*gl%z3_PCSAFT(1)**2 + 18432.d0*gl%z3_PCSAFT(1))
	part6 = (48.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 2880.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 67160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 903440.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 8176508.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4) - 53892096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 271015922.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 1070687860.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 3385354730.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 8663404360.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 18039562700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) - 30557850792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 41839042960.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 45617534880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) + 38477437298.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 23678103492.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 9148075074.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) - 872568120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 1159095240.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 565209648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) - 22728744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 42761952.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 3559680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 1866240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 165888.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) + 168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 4320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(4) + 48040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 274000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4) + 448882.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 5765408.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 59029912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + &
     & 317522760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 1201050770.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 &
     & *gl%z3_PCSAFT(4) + 3470093016.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 7948634256.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 14732219936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) &
     & - 22365758570.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 28017643280.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 29050093976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + &
     & 24848823688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 17256813414.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **10*gl%z3_PCSAFT(4) + 9313382840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 3475735560.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 542140512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 242055816.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 159804000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4) + 21610080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 7200000.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 1534464.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 221184.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 96.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 2880.d0*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(4) - 40540.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 354760.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) &
     & - 2152342.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 9532496.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - &
     & 31410798.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 75804340.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - &
     & 121301600.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 61166672.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + &
     & 336318240.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 1354543648.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + &
     & 3087957650.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 5168775680.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 6771059770.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 7071260076.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + &
     & 5878093652.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 3827366560.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 1887824080.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 659115504.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 134820464.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 340992.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 10324320.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 3552000.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 527616.d0 &
     & *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 27648.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 18432.d0*gl%z3_PCSAFT(4))
	part7 = (-24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + 1440.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**25 - 33580.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 451720.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**23 - 4088254.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 26946048.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**21 - 135507961.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 535343930.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**19 - 1692677365.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 4331702180.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**17 - 9019781350.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 15278925396.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**15 - 20919521480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 22808767440.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 19238718649.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + &
     & 11839051746.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 4574037537.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 &
     & + 436284060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 579547620.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - &
     & 282604824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 11364372.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + &
     & 21380976.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 1779840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 933120.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 82944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 36.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**26 + 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 12000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 &
     & - 286800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 3293511.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - &
     & 24818016.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 136188279.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - &
     & 573111570.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 1905588870.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - &
     & 5097649632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 11087497362.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 &
     & - 19701027372.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 28584823095.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 14 - 33658581720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 31726353327.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**12 - 23327939706.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 12717843648.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10 - 4544959800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 588824280.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**8 + 346445256.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 205677312.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**6 + 34089120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 5393520.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**4 - 2424960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 &
     & + 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 54.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 1320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **25 + 13185.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 51990.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 227247.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 4605696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 34512198.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**20 + 171170340.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 633267660.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **18 + 1835672832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 4276282515.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 &
     & + 8117995962.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 12645195375.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 16184435880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 16951792440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
     & 14373506184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 9652408188.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
     & 4908619920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 1689814080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + &
     & 232131936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 118227264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 80900928.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 18336000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 716160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **3 - 1094976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 99072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 18432.d0*gl%mmean_PCSAFT &
     & - 15.d0*gl%z3_PCSAFT(1)**26 + 420.d0*gl%z3_PCSAFT(1)**25 - 5395.d0*gl%z3_PCSAFT(1)**24 + 41230.d0*gl%z3_PCSAFT(1) &
     & **23 - 198682.d0*gl%z3_PCSAFT(1)**22 + 532320.d0*gl%z3_PCSAFT(1)**21 + 161282.d0*gl%z3_PCSAFT(1)**20 &
     & - 9336100.d0*gl%z3_PCSAFT(1)**19 + 52451405.d0*gl%z3_PCSAFT(1)**18 - 187558132.d0*gl%z3_PCSAFT(1)**17 &
     & + 502433633.d0*gl%z3_PCSAFT(1)**16 - 1065610706.d0*gl%z3_PCSAFT(1)**15 + 1834610820.d0* &
     & gl%z3_PCSAFT(1)**14 - 2593730240.d0*gl%z3_PCSAFT(1)**13 + 3020711744.d0*gl%z3_PCSAFT(1)**12 - &
     & 2886136328.d0*gl%z3_PCSAFT(1)**11 + 2233829944.d0*gl%z3_PCSAFT(1)**10 - 1364462240.d0*gl%z3_PCSAFT(1) &
     & **9 + 622720720.d0*gl%z3_PCSAFT(1)**8 - 183443264.d0*gl%z3_PCSAFT(1)**7 + 12598272.d0*gl%z3_PCSAFT(1) &
     & **6 + 17818112.d0*gl%z3_PCSAFT(1)**5 - 9341760.d0*gl%z3_PCSAFT(1)**4 + 2028160.d0*gl%z3_PCSAFT(1)**3 &
     & - 66816.d0*gl%z3_PCSAFT(1)**2 - 56832.d0*gl%z3_PCSAFT(1) + 8448.d0)
	part8 = (-24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + 1440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 - &
     & 33580.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 451720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - 4088254.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 26946048.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - 135507961.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 535343930.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 1692677365.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 4331702180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 9019781350.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 15278925396.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - &
     & 20919521480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 22808767440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 &
     & - 19238718649.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + 11839051746.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 11 - 4574037537.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + 436284060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 9 + 579547620.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - 282604824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + &
     & 11364372.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 21380976.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - &
     & 1779840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 933120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 82944.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 36.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **25 + 12000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 286800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + &
     & 3293511.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 24818016.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + &
     & 136188279.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 573111570.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + &
     & 1905588870.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 5097649632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + &
     & 11087497362.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 19701027372.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & + 28584823095.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 33658581720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13 + 31726353327.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 23327939706.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**11 + 12717843648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 4544959800.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**9 + 588824280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 346445256.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**7 - 205677312.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 34089120.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**5 + 5393520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 2424960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **3 + 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 54.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**26 - 1320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 13185.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - &
     & 51990.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 227247.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 4605696.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21 - 34512198.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 171170340.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 19 - 633267660.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 1835672832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - &
     & 4276282515.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 8117995962.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - &
     & 12645195375.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 16184435880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - &
     & 16951792440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 14373506184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - &
     & 9652408188.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 4908619920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - &
     & 1689814080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 232131936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + &
     & 118227264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 80900928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 18336000.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 716160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 1094976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 2 + 99072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 18432.d0*gl%mmean_PCSAFT - 15.d0*gl%z3_PCSAFT(1)**26 + 420.d0*gl%z3_PCSAFT(1) &
     & **25 - 5395.d0*gl%z3_PCSAFT(1)**24 + 41230.d0*gl%z3_PCSAFT(1)**23 - 198682.d0*gl%z3_PCSAFT(1)**22 + &
     & 532320.d0*gl%z3_PCSAFT(1)**21 + 161282.d0*gl%z3_PCSAFT(1)**20 - 9336100.d0*gl%z3_PCSAFT(1)**19 + &
     & 52451405.d0*gl%z3_PCSAFT(1)**18 - 187558132.d0*gl%z3_PCSAFT(1)**17 + 502433633.d0*gl%z3_PCSAFT(1)**16 &
     & - 1065610706.d0*gl%z3_PCSAFT(1)**15 + 1834610820.d0*gl%z3_PCSAFT(1)**14 - 2593730240.d0* &
     & gl%z3_PCSAFT(1)**13 + 3020711744.d0*gl%z3_PCSAFT(1)**12 - 2886136328.d0*gl%z3_PCSAFT(1)**11 + &
     & 2233829944.d0*gl%z3_PCSAFT(1)**10 - 1364462240.d0*gl%z3_PCSAFT(1)**9 + 622720720.d0*gl%z3_PCSAFT(1)** &
     & 8 - 183443264.d0*gl%z3_PCSAFT(1)**7 + 12598272.d0*gl%z3_PCSAFT(1)**6 + 17818112.d0*gl%z3_PCSAFT(1)**5 &
     & - 9341760.d0*gl%z3_PCSAFT(1)**4 + 2028160.d0*gl%z3_PCSAFT(1)**3 - 66816.d0*gl%z3_PCSAFT(1)**2 - 56832.d0 &
     & *gl%z3_PCSAFT(1) + 8448.d0)
	part9 = (96.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 25*gl%z3_PCSAFT(4) - 6240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 145440.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 1911840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + &
     & 16776072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 107075112.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **20*gl%z3_PCSAFT(4) + 522965640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - &
     & 2018666280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 6288866640.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 16019373360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + &
     & 33585528240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 58000880640.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 82080870360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 94029307800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 85289793912.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 58900185384.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + &
     & 28644768480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 7907286960.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 126189360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 892093392.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 187625808.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(4) - 31687200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 5339520.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 2799360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & + 248832.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(4) + 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 79080.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 1503000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - &
     & 15472608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 107658144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **20*gl%z3_PCSAFT(4) - 553996680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + &
     & 2214098520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 7081914000.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 18478776096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - &
     & 39798829608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 71186731320.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 105890970000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 130609564320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 132613716216.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 109481915064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) &
     & - 72060512400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 36569239920.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 13361850480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 2865174432.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 53443008.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(4) - 251110080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 59037120.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 3715200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & - 2156544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 248832.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) - &
     & 216.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 4440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - &
     & 29460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 84420.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) &
     & + 3098136.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 29140944.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) + 172725900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 750534060.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 2545545000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - &
     & 6972789336.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 15747258180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 15*gl%z3_PCSAFT(4) - 29691231180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 47059454280.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 62838881280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4) + 70551739140.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 66173335476.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 51233990160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
     & - 32100083760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 15748180800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **7*gl%z3_PCSAFT(4) - 5703883728.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 1340129664.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 120054720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - &
     & 31609920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 10270080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) - 709632.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 27648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) + &
     & 60.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 1500.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 16320.d0*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4) - 94560.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 222156.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) + 1038300.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 13210320.d0*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 73893720.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 287438940.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) + 865340796.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 2110814112.d0*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) + 4273118160.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 7271778300.d0*gl%z3_PCSAFT(1)**13 &
     & *gl%z3_PCSAFT(4) + 10462467300.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 12734595312.d0*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) + 13066065432.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 11217968880.d0* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 7967593440.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 4604185440.d0* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 2112183072.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 740144448.d0* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 184923840.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 28049280.d0* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 1059840.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 415488.d0*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) - 56832.d0*gl%z3_PCSAFT(4))
	part10 = (-16.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**27 + 880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 - 20264.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 &
     & + 276200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 2562084.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + &
     & 17398812.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 90287438.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + &
     & 367632998.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 1193711530.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + &
     & 3118077930.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 6569196380.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + &
     & 11124192988.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 14977983864.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & + 15687599960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 12242103214.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13 + 6446535814.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 1575646146.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **11 - 534997854.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 573873600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **9 - 131352624.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 36773352.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 &
     & + 19154448.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 907200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - &
     & 933120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 82944.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 56.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**27 + 1640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 21976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + &
     & 173640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 843446.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 2031034.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**22 + 3879688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 62163308.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 20 + 326625910.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 1138151678.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + &
     & 2957603032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 5965444904.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + &
     & 9469367070.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 11819551330.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 11434547288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 8281748236.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
     & 4138672946.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 1079963698.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - &
     & 167460800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 258356184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 78703512.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 3034320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 6868512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **5 - 901440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 165888.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 27648.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 32.d0*gl%z3_PCSAFT(1)**27 - 1040.d0*gl%z3_PCSAFT(1)**26 + 16108.d0*gl%z3_PCSAFT(1) &
     & **25 - 157580.d0*gl%z3_PCSAFT(1)**24 + 1088146.d0*gl%z3_PCSAFT(1)**23 - 5614202.d0*gl%z3_PCSAFT(1)** &
     & 22 + 22331162.d0*gl%z3_PCSAFT(1)**21 - 69551646.d0*gl%z3_PCSAFT(1)**20 + 169965900.d0*gl%z3_PCSAFT(1) &
     & **19 - 320797416.d0*gl%z3_PCSAFT(1)**18 + 442978920.d0*gl%z3_PCSAFT(1)**17 - 364557048.d0* &
     & gl%z3_PCSAFT(1)**16 - 81361550.d0*gl%z3_PCSAFT(1)**15 + 860915150.d0*gl%z3_PCSAFT(1)**14 - &
     & 1638584830.d0*gl%z3_PCSAFT(1)**13 + 1967101162.d0*gl%z3_PCSAFT(1)**12 - 1661424768.d0*gl%z3_PCSAFT(1) &
     & **11 + 969145860.d0*gl%z3_PCSAFT(1)**10 - 331113680.d0*gl%z3_PCSAFT(1)**9 + 2950712.d0* &
     & gl%z3_PCSAFT(1)**8 + 60759632.d0*gl%z3_PCSAFT(1)**7 - 28279328.d0*gl%z3_PCSAFT(1)**6 + 2922528.d0* &
     & gl%z3_PCSAFT(1)**5 + 1838400.d0*gl%z3_PCSAFT(1)**4 - 575232.d0*gl%z3_PCSAFT(1)**3 - 13824.d0*gl%z3_PCSAFT(1)** &
     & 2 + 18432.d0*gl%z3_PCSAFT(1))
	part11 = (48.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26* &
     & gl%z3_PCSAFT(4) - 2880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 67160.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 903440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + &
     & 8176508.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 53892096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 21*gl%z3_PCSAFT(4) + 271015922.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 1070687860.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 3385354730.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) - 8663404360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 18039562700.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 30557850792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) + 41839042960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 45617534880.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 38477437298.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4) - 23678103492.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 9148075074.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 872568120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) - 1159095240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 565209648.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 22728744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4) - 42761952.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 3559680.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 1866240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + &
     & 165888.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26* &
     & gl%z3_PCSAFT(4) - 4320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 48040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 24*gl%z3_PCSAFT(4) - 274000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 448882.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 5765408.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - &
     & 59029912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 317522760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 1201050770.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 3470093016.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 7948634256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
     & + 14732219936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 22365758570.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 28017643280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 29050093976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 24848823688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) - 17256813414.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + &
     & 9313382840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 3475735560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4) + 542140512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 242055816.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 159804000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 21610080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 7200000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) - 1534464.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 221184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & *gl%z3_PCSAFT(4) - 96.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 2880.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - &
     & 40540.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 354760.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 2152342.d0* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 9532496.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 31410798.d0* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 75804340.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 121301600.d0* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 61166672.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 336318240.d0* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 1354543648.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 3087957650.d0 &
     & *gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 5168775680.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 6771059770.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 7071260076.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + &
     & 5878093652.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 3827366560.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 1887824080.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 659115504.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 134820464.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 340992.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 10324320.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 3552000.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 527616.d0 &
     & *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 27648.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 18432.d0*gl%z3_PCSAFT(4))
	part12 = (-24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + 1440.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**25 - 33580.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 451720.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**23 - 4088254.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 26946048.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**21 - 135507961.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 535343930.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**19 - 1692677365.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 4331702180.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**17 - 9019781350.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 15278925396.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**15 - 20919521480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 22808767440.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 19238718649.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + &
     & 11839051746.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 4574037537.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 &
     & + 436284060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 579547620.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - &
     & 282604824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 11364372.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + &
     & 21380976.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 1779840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 933120.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 82944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 36.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**26 + 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 12000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 &
     & - 286800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 3293511.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - &
     & 24818016.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 136188279.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - &
     & 573111570.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 1905588870.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - &
     & 5097649632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 11087497362.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 &
     & - 19701027372.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 28584823095.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 14 - 33658581720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 31726353327.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**12 - 23327939706.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 12717843648.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10 - 4544959800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 588824280.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**8 + 346445256.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 205677312.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**6 + 34089120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 5393520.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**4 - 2424960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 &
     & + 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 54.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 1320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **25 + 13185.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 51990.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 227247.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 4605696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 34512198.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**20 + 171170340.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 633267660.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **18 + 1835672832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 4276282515.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 &
     & + 8117995962.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 12645195375.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 16184435880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 16951792440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
     & 14373506184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 9652408188.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
     & 4908619920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 1689814080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + &
     & 232131936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 118227264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 80900928.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 18336000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 716160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **3 - 1094976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 99072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 18432.d0*gl%mmean_PCSAFT &
     & - 15.d0*gl%z3_PCSAFT(1)**26 + 420.d0*gl%z3_PCSAFT(1)**25 - 5395.d0*gl%z3_PCSAFT(1)**24 + 41230.d0*gl%z3_PCSAFT(1) &
     & **23 - 198682.d0*gl%z3_PCSAFT(1)**22 + 532320.d0*gl%z3_PCSAFT(1)**21 + 161282.d0*gl%z3_PCSAFT(1)**20 &
     & - 9336100.d0*gl%z3_PCSAFT(1)**19 + 52451405.d0*gl%z3_PCSAFT(1)**18 - 187558132.d0*gl%z3_PCSAFT(1)**17 &
     & + 502433633.d0*gl%z3_PCSAFT(1)**16 - 1065610706.d0*gl%z3_PCSAFT(1)**15 + 1834610820.d0* &
     & gl%z3_PCSAFT(1)**14 - 2593730240.d0*gl%z3_PCSAFT(1)**13 + 3020711744.d0*gl%z3_PCSAFT(1)**12 - &
     & 2886136328.d0*gl%z3_PCSAFT(1)**11 + 2233829944.d0*gl%z3_PCSAFT(1)**10 - 1364462240.d0*gl%z3_PCSAFT(1) &
     & **9 + 622720720.d0*gl%z3_PCSAFT(1)**8 - 183443264.d0*gl%z3_PCSAFT(1)**7 + 12598272.d0*gl%z3_PCSAFT(1) &
     & **6 + 17818112.d0*gl%z3_PCSAFT(1)**5 - 9341760.d0*gl%z3_PCSAFT(1)**4 + 2028160.d0*gl%z3_PCSAFT(1)**3 &
     & - 66816.d0*gl%z3_PCSAFT(1)**2 - 56832.d0*gl%z3_PCSAFT(1) + 8448.d0)
	part13 = (-24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + 1440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 - &
     & 33580.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 451720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - 4088254.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 26946048.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - 135507961.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 535343930.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 1692677365.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 4331702180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 9019781350.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 15278925396.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - &
     & 20919521480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 22808767440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 &
     & - 19238718649.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + 11839051746.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 11 - 4574037537.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + 436284060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 9 + 579547620.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - 282604824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + &
     & 11364372.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 21380976.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - &
     & 1779840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 933120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 82944.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 36.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **25 + 12000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 286800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + &
     & 3293511.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 24818016.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + &
     & 136188279.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 573111570.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + &
     & 1905588870.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 5097649632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + &
     & 11087497362.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 19701027372.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & + 28584823095.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 33658581720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13 + 31726353327.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 23327939706.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**11 + 12717843648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 4544959800.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**9 + 588824280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 346445256.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**7 - 205677312.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 34089120.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**5 + 5393520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 2424960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **3 + 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 54.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**26 - 1320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 13185.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - &
     & 51990.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 227247.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 4605696.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21 - 34512198.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 171170340.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 19 - 633267660.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 1835672832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - &
     & 4276282515.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 8117995962.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - &
     & 12645195375.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 16184435880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - &
     & 16951792440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 14373506184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - &
     & 9652408188.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 4908619920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - &
     & 1689814080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 232131936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + &
     & 118227264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 80900928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 18336000.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 716160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 1094976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 2 + 99072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 18432.d0*gl%mmean_PCSAFT - 15.d0*gl%z3_PCSAFT(1)**26 + 420.d0*gl%z3_PCSAFT(1) &
     & **25 - 5395.d0*gl%z3_PCSAFT(1)**24 + 41230.d0*gl%z3_PCSAFT(1)**23 - 198682.d0*gl%z3_PCSAFT(1)**22 + &
     & 532320.d0*gl%z3_PCSAFT(1)**21 + 161282.d0*gl%z3_PCSAFT(1)**20 - 9336100.d0*gl%z3_PCSAFT(1)**19 + &
     & 52451405.d0*gl%z3_PCSAFT(1)**18 - 187558132.d0*gl%z3_PCSAFT(1)**17 + 502433633.d0*gl%z3_PCSAFT(1)**16 &
     & - 1065610706.d0*gl%z3_PCSAFT(1)**15 + 1834610820.d0*gl%z3_PCSAFT(1)**14 - 2593730240.d0* &
     & gl%z3_PCSAFT(1)**13 + 3020711744.d0*gl%z3_PCSAFT(1)**12 - 2886136328.d0*gl%z3_PCSAFT(1)**11 + &
     & 2233829944.d0*gl%z3_PCSAFT(1)**10 - 1364462240.d0*gl%z3_PCSAFT(1)**9 + 622720720.d0*gl%z3_PCSAFT(1)** &
     & 8 - 183443264.d0*gl%z3_PCSAFT(1)**7 + 12598272.d0*gl%z3_PCSAFT(1)**6 + 17818112.d0*gl%z3_PCSAFT(1)**5 &
     & - 9341760.d0*gl%z3_PCSAFT(1)**4 + 2028160.d0*gl%z3_PCSAFT(1)**3 - 66816.d0*gl%z3_PCSAFT(1)**2 - 56832.d0 &
     & *gl%z3_PCSAFT(1) + 8448.d0)
	part14 = (96.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 25*gl%z3_PCSAFT(4) - 6240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 145440.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 1911840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + &
     & 16776072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 107075112.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **20*gl%z3_PCSAFT(4) + 522965640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - &
     & 2018666280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 6288866640.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 16019373360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + &
     & 33585528240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 58000880640.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 82080870360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 94029307800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 85289793912.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 58900185384.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + &
     & 28644768480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 7907286960.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 126189360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 892093392.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 187625808.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(4) - 31687200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 5339520.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 2799360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & + 248832.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(4) + 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 79080.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 1503000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - &
     & 15472608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 107658144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **20*gl%z3_PCSAFT(4) - 553996680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + &
     & 2214098520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 7081914000.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 18478776096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - &
     & 39798829608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 71186731320.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 105890970000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 130609564320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 132613716216.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 109481915064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) &
     & - 72060512400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 36569239920.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 13361850480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 2865174432.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 53443008.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(4) - 251110080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 59037120.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 3715200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & - 2156544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 248832.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) - &
     & 216.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 4440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - &
     & 29460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 84420.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) &
     & + 3098136.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 29140944.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) + 172725900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 750534060.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 2545545000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - &
     & 6972789336.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 15747258180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 15*gl%z3_PCSAFT(4) - 29691231180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 47059454280.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 62838881280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4) + 70551739140.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 66173335476.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 51233990160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
     & - 32100083760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 15748180800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **7*gl%z3_PCSAFT(4) - 5703883728.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 1340129664.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 120054720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - &
     & 31609920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 10270080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) - 709632.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 27648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) + &
     & 60.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 1500.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 16320.d0*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4) - 94560.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 222156.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) + 1038300.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 13210320.d0*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 73893720.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 287438940.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) + 865340796.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 2110814112.d0*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) + 4273118160.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 7271778300.d0*gl%z3_PCSAFT(1)**13 &
     & *gl%z3_PCSAFT(4) + 10462467300.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 12734595312.d0*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) + 13066065432.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 11217968880.d0* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 7967593440.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 4604185440.d0* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 2112183072.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 740144448.d0* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 184923840.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 28049280.d0* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 1059840.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 415488.d0*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) - 56832.d0*gl%z3_PCSAFT(4))
	part15 = (-24.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + 1440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 - 33580.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**24 + 451720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - 4088254.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **22 + 26946048.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - 135507961.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 20 + 535343930.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 1692677365.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 18 + 4331702180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 9019781350.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **16 + 15278925396.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 20919521480.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**14 + 22808767440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 19238718649.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**12 + 11839051746.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 4574037537.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + 436284060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 579547620.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - 282604824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 11364372.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 21380976.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 1779840.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**4 - 933120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 82944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **2 - 36.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 12000.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 286800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 3293511.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**22 - 24818016.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 136188279.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**20 - 573111570.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 1905588870.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**18 - 5097649632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 11087497362.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**16 - 19701027372.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 28584823095.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 33658581720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + &
     & 31726353327.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 23327939706.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 &
     & + 12717843648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 4544959800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 &
     & + 588824280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 346445256.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - &
     & 205677312.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 34089120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + &
     & 5393520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 2424960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 96768.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 54.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 &
     & - 1320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 13185.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 51990.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**23 - 227247.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 4605696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - &
     & 34512198.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 171170340.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - &
     & 633267660.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 1835672832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - &
     & 4276282515.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 8117995962.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - &
     & 12645195375.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 16184435880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - &
     & 16951792440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 14373506184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - &
     & 9652408188.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 4908619920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - &
     & 1689814080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 232131936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + &
     & 118227264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 80900928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 18336000.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 716160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 1094976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 2 + 99072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 18432.d0*gl%mmean_PCSAFT - 15.d0*gl%z3_PCSAFT(1)**26 + 420.d0*gl%z3_PCSAFT(1) &
     & **25 - 5395.d0*gl%z3_PCSAFT(1)**24 + 41230.d0*gl%z3_PCSAFT(1)**23 - 198682.d0*gl%z3_PCSAFT(1)**22 + &
     & 532320.d0*gl%z3_PCSAFT(1)**21 + 161282.d0*gl%z3_PCSAFT(1)**20 - 9336100.d0*gl%z3_PCSAFT(1)**19 + &
     & 52451405.d0*gl%z3_PCSAFT(1)**18 - 187558132.d0*gl%z3_PCSAFT(1)**17 + 502433633.d0*gl%z3_PCSAFT(1)**16 &
     & - 1065610706.d0*gl%z3_PCSAFT(1)**15 + 1834610820.d0*gl%z3_PCSAFT(1)**14 - 2593730240.d0* &
     & gl%z3_PCSAFT(1)**13 + 3020711744.d0*gl%z3_PCSAFT(1)**12 - 2886136328.d0*gl%z3_PCSAFT(1)**11 + &
     & 2233829944.d0*gl%z3_PCSAFT(1)**10 - 1364462240.d0*gl%z3_PCSAFT(1)**9 + 622720720.d0*gl%z3_PCSAFT(1)** &
     & 8 - 183443264.d0*gl%z3_PCSAFT(1)**7 + 12598272.d0*gl%z3_PCSAFT(1)**6 + 17818112.d0*gl%z3_PCSAFT(1)**5 &
     & - 9341760.d0*gl%z3_PCSAFT(1)**4 + 2028160.d0*gl%z3_PCSAFT(1)**3 - 66816.d0*gl%z3_PCSAFT(1)**2 - 56832.d0 &
     & *gl%z3_PCSAFT(1) + 8448.d0)
	part16 = (-24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **26 + 1440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 - 33580.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + &
     & 451720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - 4088254.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + &
     & 26946048.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - 135507961.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + &
     & 535343930.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 1692677365.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + &
     & 4331702180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 9019781350.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + &
     & 15278925396.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 20919521480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 &
     & + 22808767440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 19238718649.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 12 + 11839051746.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 4574037537.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **10 + 436284060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 579547620.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 8 - 282604824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 11364372.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + &
     & 21380976.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 1779840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 933120.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 82944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 36.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**26 + 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 12000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 &
     & - 286800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 3293511.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - &
     & 24818016.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 136188279.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - &
     & 573111570.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 1905588870.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - &
     & 5097649632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 11087497362.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 &
     & - 19701027372.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 28584823095.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 14 - 33658581720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 31726353327.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**12 - 23327939706.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 12717843648.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10 - 4544959800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 588824280.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**8 + 346445256.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 205677312.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**6 + 34089120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 5393520.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**4 - 2424960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 96768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 &
     & + 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 54.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 1320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **25 + 13185.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 51990.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 227247.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 4605696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 34512198.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**20 + 171170340.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 633267660.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **18 + 1835672832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 4276282515.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 &
     & + 8117995962.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 12645195375.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 16184435880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 16951792440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
     & 14373506184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 9652408188.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
     & 4908619920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 1689814080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + &
     & 232131936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 118227264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 80900928.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 18336000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 716160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **3 - 1094976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 99072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 18432.d0*gl%mmean_PCSAFT &
     & - 15.d0*gl%z3_PCSAFT(1)**26 + 420.d0*gl%z3_PCSAFT(1)**25 - 5395.d0*gl%z3_PCSAFT(1)**24 + 41230.d0*gl%z3_PCSAFT(1) &
     & **23 - 198682.d0*gl%z3_PCSAFT(1)**22 + 532320.d0*gl%z3_PCSAFT(1)**21 + 161282.d0*gl%z3_PCSAFT(1)**20 &
     & - 9336100.d0*gl%z3_PCSAFT(1)**19 + 52451405.d0*gl%z3_PCSAFT(1)**18 - 187558132.d0*gl%z3_PCSAFT(1)**17 &
     & + 502433633.d0*gl%z3_PCSAFT(1)**16 - 1065610706.d0*gl%z3_PCSAFT(1)**15 + 1834610820.d0* &
     & gl%z3_PCSAFT(1)**14 - 2593730240.d0*gl%z3_PCSAFT(1)**13 + 3020711744.d0*gl%z3_PCSAFT(1)**12 - &
     & 2886136328.d0*gl%z3_PCSAFT(1)**11 + 2233829944.d0*gl%z3_PCSAFT(1)**10 - 1364462240.d0*gl%z3_PCSAFT(1) &
     & **9 + 622720720.d0*gl%z3_PCSAFT(1)**8 - 183443264.d0*gl%z3_PCSAFT(1)**7 + 12598272.d0*gl%z3_PCSAFT(1) &
     & **6 + 17818112.d0*gl%z3_PCSAFT(1)**5 - 9341760.d0*gl%z3_PCSAFT(1)**4 + 2028160.d0*gl%z3_PCSAFT(1)**3 &
     & - 66816.d0*gl%z3_PCSAFT(1)**2 - 56832.d0*gl%z3_PCSAFT(1) + 8448.d0)
	part17 = (96.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 6240.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 145440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - &
     & 1911840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 16776072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 21*gl%z3_PCSAFT(4) - 107075112.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 522965640.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 2018666280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) + 6288866640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 16019373360.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 33585528240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) - 58000880640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 82080870360.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 94029307800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4) + 85289793912.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 58900185384.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 28644768480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) - 7907286960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 126189360.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 892093392.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4) - 187625808.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 31687200.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 5339520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + &
     & 2799360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 248832.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) + 144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **24*gl%z3_PCSAFT(4) - 79080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 1503000.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 15472608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) + 107658144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 553996680.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 2214098520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) - 7081914000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 18478776096.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 39798829608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) + 71186731320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - &
     & 105890970000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 130609564320.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 132613716216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) &
     & + 109481915064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 72060512400.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 36569239920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - &
     & 13361850480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 2865174432.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 53443008.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 251110080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 59037120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4) + 3715200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 2156544.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 248832.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) - 216.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 4440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 29460.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 84420.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + &
     & 3098136.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 29140944.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) + 172725900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 750534060.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 2545545000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - &
     & 6972789336.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 15747258180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 15*gl%z3_PCSAFT(4) - 29691231180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 47059454280.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 62838881280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4) + 70551739140.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 66173335476.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 51233990160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
     & - 32100083760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 15748180800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **7*gl%z3_PCSAFT(4) - 5703883728.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 1340129664.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 120054720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - &
     & 31609920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 10270080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) - 709632.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 27648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) + &
     & 60.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 1500.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 16320.d0*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4) - 94560.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 222156.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) + 1038300.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 13210320.d0*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 73893720.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 287438940.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) + 865340796.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 2110814112.d0*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) + 4273118160.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 7271778300.d0*gl%z3_PCSAFT(1)**13 &
     & *gl%z3_PCSAFT(4) + 10462467300.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 12734595312.d0*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) + 13066065432.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 11217968880.d0* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 7967593440.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 4604185440.d0* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 2112183072.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 740144448.d0* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 184923840.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 28049280.d0* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 1059840.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 415488.d0*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) - 56832.d0*gl%z3_PCSAFT(4))
	part18 = ( &
     & -96.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25 + 6240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 - 145440.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**23 + 1911840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 - 16776072.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**21 + 107075112.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 - 522965640.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**19 + 2018666280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 - 6288866640.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**17 + 16019373360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 - 33585528240.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**15 + 58000880640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 - 82080870360.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 + 94029307800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 - &
     & 85289793912.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 + 58900185384.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 &
     & - 28644768480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9 + 7907286960.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 &
     & + 126189360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 - 892093392.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 + &
     & 187625808.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 + 31687200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - &
     & 5339520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 - 2799360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 - 248832.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) + 48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 - 7920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 24 + 221520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - 3152400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + &
     & 28941804.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - 190455852.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + &
     & 952550940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 3754754460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + &
     & 11942692920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 31116396648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 &
     & + 66986574264.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 119527920120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **14 + 176425889820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 213802218540.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**12 + 209719300188.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 162538853052.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + 95419071360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 - &
     & 38934658800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 8505887040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + &
     & 694809504.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 1032850224.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + &
     & 228048480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 120960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - &
     & 1969920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 787968.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) - 82944.d0*gl%mmean_PCSAFT &
     & **3 + 72.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 2520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 114120.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 1860840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 18248292.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 125061948.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 644258700.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 2603225700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 8476315200.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 22631949072.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - &
     & 50083385040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 92346809040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 &
     & - 141954858900.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 181183958460.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **12 - 190264462380.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 161813105172.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**10 - 108566922960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 54830580480.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**8 - 18800701200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 2957984496.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**6 + 754654752.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 556300800.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**4 + 124807680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 7637760.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**2 - 781056.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 13824.d0*gl%mmean_PCSAFT**2 - 60.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**25 + 300.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 + 20700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - &
     & 452460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 4930056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 35759640.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 191352480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 796675200.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**18 + 2664699300.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 7308212004.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**16 + 16648244988.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 31731225900.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14 + 50737653600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 67958839920.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12 + 75810912168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 69678998568.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10 + 51829469040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 30268926720.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**8 + 13098893280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 3624458688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **6 + 240568512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 287882880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - &
     & 136596480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 27648000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 2072832.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 56832.d0*gl%mmean_PCSAFT + 12.d0*gl%z3_PCSAFT(1)**25 - 180.d0*gl%z3_PCSAFT(1)**24 - 420.d0 &
     & *gl%z3_PCSAFT(1)**23 + 35580.d0*gl%z3_PCSAFT(1)**22 - 472176.d0*gl%z3_PCSAFT(1)**21 + 3729120.d0* &
     & gl%z3_PCSAFT(1)**20 - 20995080.d0*gl%z3_PCSAFT(1)**19 + 90708840.d0*gl%z3_PCSAFT(1)**18 - &
     & 313053300.d0*gl%z3_PCSAFT(1)**17 + 884457516.d0*gl%z3_PCSAFT(1)**16 - 2077503444.d0*gl%z3_PCSAFT(1)** &
     & 15 + 4094856300.d0*gl%z3_PCSAFT(1)**14 - 6803366280.d0*gl%z3_PCSAFT(1)**13 + 9532220040.d0* &
     & gl%z3_PCSAFT(1)**12 - 11224883712.d0*gl%z3_PCSAFT(1)**11 + 11025889392.d0*gl%z3_PCSAFT(1)**10 - &
     & 8918396640.d0*gl%z3_PCSAFT(1)**9 + 5815012320.d0*gl%z3_PCSAFT(1)**8 - 2943470400.d0*gl%z3_PCSAFT(1)** &
     & 7 + 1067950080.d0*gl%z3_PCSAFT(1)**6 - 213672192.d0*gl%z3_PCSAFT(1)**5 - 22790400.d0*gl%z3_PCSAFT(1) &
     & **4 + 33719040.d0*gl%z3_PCSAFT(1)**3 - 11953920.d0*gl%z3_PCSAFT(1)**2 + 2145792.d0*gl%z3_PCSAFT(1) - &
     & 165888.d0)
	part19 = (-96.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **25 + 6240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 - 145440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 + &
     & 1911840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 - 16776072.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 + &
     & 107075112.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 - 522965640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 + &
     & 2018666280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 - 6288866640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 + &
     & 16019373360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 - 33585528240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 &
     & + 58000880640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 - 82080870360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 13 + 94029307800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 - 85289793912.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**11 + 58900185384.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 - 28644768480.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**9 + 7907286960.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 + 126189360.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)**7 - 892093392.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 + 187625808.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**5 + 31687200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - 5339520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **3 - 2799360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 - 248832.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) + 48.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 - 7920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 221520.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**23 - 3152400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 28941804.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**21 - 190455852.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 952550940.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**19 - 3754754460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 11942692920.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**17 - 31116396648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 66986574264.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 119527920120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + &
     & 176425889820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 213802218540.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 12 + 209719300188.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 162538853052.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**10 + 95419071360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 - 38934658800.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**8 + 8505887040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 694809504.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**6 - 1032850224.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 228048480.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**4 + 120960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 1969920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 2 - 787968.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) - 82944.d0*gl%mmean_PCSAFT**3 + 72.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **25 + 2520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 114120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + &
     & 1860840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 18248292.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + &
     & 125061948.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 644258700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + &
     & 2603225700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 8476315200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + &
     & 22631949072.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 50083385040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & + 92346809040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 141954858900.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **13 + 181183958460.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 190264462380.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**11 + 161813105172.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 108566922960.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 54830580480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 18800701200.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 2957984496.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 754654752.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 556300800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 124807680.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 7637760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 781056.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1) - 13824.d0*gl%mmean_PCSAFT**2 - 60.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 300.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **24 + 20700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 452460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 4930056.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 35759640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 191352480.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19 - 796675200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 2664699300.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **17 - 7308212004.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 16648244988.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 &
     & - 31731225900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 50737653600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - &
     & 67958839920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 75810912168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - &
     & 69678998568.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 51829469040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - &
     & 30268926720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 13098893280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - &
     & 3624458688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 240568512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + &
     & 287882880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 136596480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 27648000.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 2072832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 56832.d0*gl%mmean_PCSAFT + 12.d0*gl%z3_PCSAFT(1) &
     & **25 - 180.d0*gl%z3_PCSAFT(1)**24 - 420.d0*gl%z3_PCSAFT(1)**23 + 35580.d0*gl%z3_PCSAFT(1)**22 - &
     & 472176.d0*gl%z3_PCSAFT(1)**21 + 3729120.d0*gl%z3_PCSAFT(1)**20 - 20995080.d0*gl%z3_PCSAFT(1)**19 + &
     & 90708840.d0*gl%z3_PCSAFT(1)**18 - 313053300.d0*gl%z3_PCSAFT(1)**17 + 884457516.d0*gl%z3_PCSAFT(1)**16 &
     & - 2077503444.d0*gl%z3_PCSAFT(1)**15 + 4094856300.d0*gl%z3_PCSAFT(1)**14 - 6803366280.d0* &
     & gl%z3_PCSAFT(1)**13 + 9532220040.d0*gl%z3_PCSAFT(1)**12 - 11224883712.d0*gl%z3_PCSAFT(1)**11 + &
     & 11025889392.d0*gl%z3_PCSAFT(1)**10 - 8918396640.d0*gl%z3_PCSAFT(1)**9 + 5815012320.d0*gl%z3_PCSAFT(1) &
     & **8 - 2943470400.d0*gl%z3_PCSAFT(1)**7 + 1067950080.d0*gl%z3_PCSAFT(1)**6 - 213672192.d0* &
     & gl%z3_PCSAFT(1)**5 - 22790400.d0*gl%z3_PCSAFT(1)**4 + 33719040.d0*gl%z3_PCSAFT(1)**3 - 11953920.d0* &
     & gl%z3_PCSAFT(1)**2 + 2145792.d0*gl%z3_PCSAFT(1) - 165888.d0)
	part20 = (-96.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25 + 6240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 - &
     & 145440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 + 1911840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 - &
     & 16776072.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 + 107075112.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 - &
     & 522965640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 + 2018666280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 - &
     & 6288866640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 + 16019373360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 &
     & - 33585528240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 + 58000880640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 14 - 82080870360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 + 94029307800.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**12 - 85289793912.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 + 58900185384.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**10 - 28644768480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9 + 7907286960.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 + 126189360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 - 892093392.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 + 187625808.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 + 31687200.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - 5339520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 - 2799360.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)**2 - 248832.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) + 48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 - &
     & 7920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 221520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - 3152400.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 28941804.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - 190455852.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 952550940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 3754754460.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 11942692920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - &
     & 31116396648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 66986574264.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & - 119527920120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 176425889820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **13 - 213802218540.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + 209719300188.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11 - 162538853052.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + 95419071360.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**9 - 38934658800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 8505887040.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**7 + 694809504.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 1032850224.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**5 + 228048480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 120960.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**3 - 1969920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 787968.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) - &
     & 82944.d0*gl%mmean_PCSAFT**3 + 72.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 2520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 24 - 114120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 1860840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - &
     & 18248292.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 125061948.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - &
     & 644258700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 2603225700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - &
     & 8476315200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 22631949072.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 &
     & - 50083385040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 92346809040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 14 - 141954858900.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 181183958460.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**12 - 190264462380.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 161813105172.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 108566922960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + &
     & 54830580480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 18800701200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + &
     & 2957984496.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 754654752.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - &
     & 556300800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 124807680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - &
     & 7637760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 781056.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 13824.d0*gl%mmean_PCSAFT &
     & **2 - 60.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 300.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 + 20700.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**23 - 452460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 4930056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - &
     & 35759640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 191352480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - &
     & 796675200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 2664699300.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - &
     & 7308212004.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 16648244988.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - &
     & 31731225900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 50737653600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - &
     & 67958839920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 75810912168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - &
     & 69678998568.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 51829469040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - &
     & 30268926720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 13098893280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - &
     & 3624458688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 240568512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + &
     & 287882880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 136596480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 27648000.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 2072832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 56832.d0*gl%mmean_PCSAFT + 12.d0*gl%z3_PCSAFT(1) &
     & **25 - 180.d0*gl%z3_PCSAFT(1)**24 - 420.d0*gl%z3_PCSAFT(1)**23 + 35580.d0*gl%z3_PCSAFT(1)**22 - &
     & 472176.d0*gl%z3_PCSAFT(1)**21 + 3729120.d0*gl%z3_PCSAFT(1)**20 - 20995080.d0*gl%z3_PCSAFT(1)**19 + &
     & 90708840.d0*gl%z3_PCSAFT(1)**18 - 313053300.d0*gl%z3_PCSAFT(1)**17 + 884457516.d0*gl%z3_PCSAFT(1)**16 &
     & - 2077503444.d0*gl%z3_PCSAFT(1)**15 + 4094856300.d0*gl%z3_PCSAFT(1)**14 - 6803366280.d0* &
     & gl%z3_PCSAFT(1)**13 + 9532220040.d0*gl%z3_PCSAFT(1)**12 - 11224883712.d0*gl%z3_PCSAFT(1)**11 + &
     & 11025889392.d0*gl%z3_PCSAFT(1)**10 - 8918396640.d0*gl%z3_PCSAFT(1)**9 + 5815012320.d0*gl%z3_PCSAFT(1) &
     & **8 - 2943470400.d0*gl%z3_PCSAFT(1)**7 + 1067950080.d0*gl%z3_PCSAFT(1)**6 - 213672192.d0* &
     & gl%z3_PCSAFT(1)**5 - 22790400.d0*gl%z3_PCSAFT(1)**4 + 33719040.d0*gl%z3_PCSAFT(1)**3 - 11953920.d0* &
     & gl%z3_PCSAFT(1)**2 + 2145792.d0*gl%z3_PCSAFT(1) - 165888.d0)
	part21 = (480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 33600.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 774720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - &
     & 9880320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 83777208.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 20*gl%z3_PCSAFT(4) - 518458800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 2474565000.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 9431202720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) + 29333020560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 75384706080.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 160954936080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4) - 285445862880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 418296100680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 501298508400.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 483148526808.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) &
     & - 364731247200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 206512323840.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 81002473920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 18305619120.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 962213472.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 231843600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - &
     & 126748800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 21358080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **2*gl%z3_PCSAFT(4) + 11197440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 995328.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(4) - 240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 45600.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 1236240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + &
     & 16828800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 147694836.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **20*gl%z3_PCSAFT(4) + 934006440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - &
     & 4530327180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 17516769840.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 55316597400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + &
     & 144756177072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 316247600760.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 578004045120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) &
     & - 881858164740.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 1116432101640.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 1161220259532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) + 977736633360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - &
     & 652195229040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 333701936640.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 124396843680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + &
     & 30776517504.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 4032573840.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 29520000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + &
     & 41765760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 760320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) - 580608.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4) - 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4) - 18000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 693360.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 10451520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + &
     & 96243948.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 626296680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **19*gl%z3_PCSAFT(4) + 3100746420.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - &
     & 12202475760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 39226842720.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 104745004128.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + &
     & 234507806280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 441980411040.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 701277133500.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) &
     & - 933775773480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 1037055940740.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 951731562960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) + 712154731680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - &
     & 426292191360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 198508224240.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 68845453344.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 16485539040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 2314529280.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 98904960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + &
     & 11289600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 71424.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) + 300.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - &
     & 153540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 2770080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) - 27361344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 185019360.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 940262520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + &
     & 3780353760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 12405215460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 16*gl%z3_PCSAFT(4) + 33868011576.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 77795442900.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 151201608960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) - 249075491640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 347300756880.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 408319680432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) + 402194923200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 328830663120.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 220258643520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) &
     & - 118671262080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 50087747712.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 15912730560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 3563527680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 496128000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) + 30612480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 415488.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) &
     & - 60.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 600.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 9540.d0*gl%z3_PCSAFT(1) &
     & **22*gl%z3_PCSAFT(4) - 258240.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 2845464.d0*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) - 20282640.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 106575240.d0*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) - 439704000.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 1477157940.d0*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) - 4130732904.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 9742040340.d0*gl%z3_PCSAFT(1)**14 &
     & *gl%z3_PCSAFT(4) - 19515883680.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 33307509840.d0*gl%z3_PCSAFT(1) &
     & **12*gl%z3_PCSAFT(4) - 48430636320.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 59856655008.d0* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 62601407520.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 55032395040.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 40282467840.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 24233623680.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 11766783360.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 4493738880.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 1298611200.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + &
     & 266760960.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 34690560.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2145792.d0 &
     & *gl%z3_PCSAFT(4))
        do xk = 1 , gl%ncomp
            do xj = 1 , gl%ncomp
                do xi = 1 , gl%ncomp
                    if (xi .GE. xj .AND. xj .GE. xk) then
            	        gl%cx3_PCSAFT(3,xi,xj,xk) = part1*(gl%mPCSAFT(xj)*gl%mPCSAFT(xi)*gl%mPCSAFT(xk)*part2 + gl%mPCSAFT(xj)*gl%mPCSAFT(xi)* &
         & gl%z3x1_PCSAFT(4,xk)*part3 + gl%mPCSAFT(xj)*gl%mPCSAFT(xi)* &
         & gl%z3x1_PCSAFT(1,xk)*part4 + &
         & gl%mPCSAFT(xj)*gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(4,xi)*part5 + gl%mPCSAFT(xj)*gl%mPCSAFT(xk)* &
         & gl%z3x1_PCSAFT(1,xi)*part6 + &
         & gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xk)*part7 + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(4,xk)* &
         & gl%z3x1_PCSAFT(1,xi)*part8 + gl%mPCSAFT(xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part9 + gl%mPCSAFT(xi)*gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(4,xj)*part10 + gl%mPCSAFT(xi)*gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(1,xj)*part11 + &
         & gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xk)*part12 + gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(4,xk)* &
         & gl%z3x1_PCSAFT(1,xj)*part13 + gl%mPCSAFT(xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*part14 + gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*part15 + gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*part16 + gl%mPCSAFT(xk)*gl%z3x1_PCSAFT(1,xj)* &
         & gl%z3x1_PCSAFT(1,xi)*part17 + gl%z3x1_PCSAFT(4,xj)*gl%z3x1_PCSAFT(1,xi)*gl%z3x1_PCSAFT(1,xk)*part18 + gl%z3x1_PCSAFT(4,xi)*gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xk)*part19 + gl%z3x1_PCSAFT(4,xk)*gl%z3x1_PCSAFT(1,xj)* &
         & gl%z3x1_PCSAFT(1,xi)*part20 + gl%z3x1_PCSAFT(1,xj)*gl%z3x1_PCSAFT(1,xi)* &
         & gl%z3x1_PCSAFT(1,xk)*part21)
                    end if
        	    end do
            end do 
        end do
    end if

    !DEC$ END IF
end subroutine CX3DERIVS




    end module pc_saft_CX3_derivs_module
