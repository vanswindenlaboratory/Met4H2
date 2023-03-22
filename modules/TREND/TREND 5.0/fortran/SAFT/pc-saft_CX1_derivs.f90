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

    ! module for file pc-saft_CX1_derivs.f90
    module pc_saft_CX1_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module


    contains




subroutine CX1DERIVS(gl,GETDERC)

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
    !output: cx1_PCSAFT
    !working variables
    double precision :: recurring_factor
    double precision :: part1, part2, part3, part4, part5, part6
    integer :: i, xi
    integer:: errorfld
    
!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! III. recurring_factor occurs in all derivatives
    recurring_factor = 2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 16.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 51.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 70.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 27.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 12.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - gl%z3_PCSAFT(1)** &
            & 6 + 8.d0*gl%z3_PCSAFT(1)**5 - 27.d0*gl%z3_PCSAFT(1)**4 + 42.d0*gl%z3_PCSAFT(1)**3 - 26.d0*gl%z3_PCSAFT(1)**2 + 4.d0
    
    ! IV. X1 DERIVATIVES of C
    !calculate the derivatives C_1, C_2, C_3, ...
    ! 1: C_1
    ! requires z3_PCSAFT(1)
    if (GETDERC(1) .eq. 1) then
        part1 = -(gl%z3_PCSAFT(1) - 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**3/recurring_factor**2
        part2 = (2.d0*gl%z3_PCSAFT(1)**8 - 22.d0*gl%z3_PCSAFT(1)**7 + 103.d0*gl%z3_PCSAFT(1)**6 - 255.d0*gl%z3_PCSAFT(1)**5 + 339.d0*gl%z3_PCSAFT(1)**4 - &
                    & 209.d0*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z3_PCSAFT(1)**2 + 24.d0*gl%z3_PCSAFT(1))
        part3 = (2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 52.d0* gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 230.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 372.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + &
                    & 192.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 24.d0*gl%mmean_PCSAFT + 2.d0*gl%z3_PCSAFT(1)**5 + 8.d0*gl%z3_PCSAFT(1)**4 - 70.d0*gl%z3_PCSAFT(1)**3 + &
                    & 148.d0*gl%z3_PCSAFT(1)**2 - 128.d0*gl%z3_PCSAFT(1) + 40.d0)
        do xi = 1, gl%ncomp
            gl%cx1_PCSAFT(1,xi) = part1*(gl%mPCSAFT(xi)*part2 + gl%z3x1_PCSAFT(1,xi)*part3)
        end do
    end if
    
    ! 2: 1ST DERIVATIVE OF C_1 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires z3_PCSAFT(1)
    if (GETDERC(2) .eq. 1) then
        part1 = 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**2/recurring_factor**3
        part2 = (2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 74.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + &
                    & 905.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 5835.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 23098.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
                    & **10 - 59594.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 101491.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 110777.d0* &
                    & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 69804.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 16584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - &
                    & 5244.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 2520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 &
                    & + 5.d0*gl%z3_PCSAFT(1)**14 - 65.d0*gl%z3_PCSAFT(1)**13 + 344.d0*gl%z3_PCSAFT(1)**12 - 818.d0*gl%z3_PCSAFT(1)**11 &
                    & - 27.d0*gl%z3_PCSAFT(1)**10 + 5879.d0*gl%z3_PCSAFT(1)**9 - 18410.d0*gl%z3_PCSAFT(1)**8 + 30864.d0*gl%z3_PCSAFT(1) &
                    & **7 - 32152.d0*gl%z3_PCSAFT(1)**6 + 20876.d0*gl%z3_PCSAFT(1)**5 - 7600.d0*gl%z3_PCSAFT(1)**4 + 864.d0* &
                    & gl%z3_PCSAFT(1)**3 + 336.d0*gl%z3_PCSAFT(1)**2 - 96.d0*gl%z3_PCSAFT(1))
        part3 = (4.d0*gl%mmean_PCSAFT**2* &
                    & gl%z3_PCSAFT(1)**13 - 178.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 2066.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - &
                    & 11977.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 42344.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 100454.d0* &
                    & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 168066.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 200491.d0*gl%mmean_PCSAFT**2* &
                    & gl%z3_PCSAFT(1)**6 + 163308.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 78696.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 &
                    & + 13416.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 2520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 288.d0*gl%mmean_PCSAFT &
                    & **2*gl%z3_PCSAFT(1) + 2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 91.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 1484.d0* &
                    & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 9593.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 35802.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 &
                    & + 88327.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 154444.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 197949.d0*gl%mmean_PCSAFT* &
                    & gl%z3_PCSAFT(1)**6 - 185876.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 122844.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - &
                    & 52072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 11544.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 576.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
                    & - 96.d0*gl%mmean_PCSAFT - 2.d0*gl%z3_PCSAFT(1)**13 - gl%z3_PCSAFT(1)**12 + 230.d0*gl%z3_PCSAFT(1)**11 - 1860.d0* &
                    & gl%z3_PCSAFT(1)**10 + 7654.d0*gl%z3_PCSAFT(1)**9 - 20217.d0*gl%z3_PCSAFT(1)**8 + 37734.d0*gl%z3_PCSAFT(1)**7 - &
                    & 52230.d0*gl%z3_PCSAFT(1)**6 + 54504.d0*gl%z3_PCSAFT(1)**5 - 42340.d0*gl%z3_PCSAFT(1)**4 + 23488.d0* &
                    & gl%z3_PCSAFT(1)**3 - 8624.d0*gl%z3_PCSAFT(1)**2 + 1824.d0*gl%z3_PCSAFT(1) - 160.d0)
        do xi = 1, gl%ncomp
            gl%cx1_PCSAFT(2,xi) = part1*(gl%mPCSAFT(xi)*part2 + gl%z3x1_PCSAFT(1,xi)*part3)
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF C_1 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires z3_PCSAFT(1)
    if (GETDERC(3) .eq. 1) then
        part1 = -2.d0*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)/recurring_factor**4
        part2 = (12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 612.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**19 + 10880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 106712.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **17 + 681631.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 3075481.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + &
     & 10251739.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 25885565.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + &
     & 50074577.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 74150471.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + &
     & 82775761.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 67043103.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + &
     & 36033036.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 9665124.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - &
     & 1146708.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 1439964.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 117792.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 69120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 6912.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**2 + 24.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 776.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 33272.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 305506.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 16 + 1642014.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 6079314.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 16549086.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 34220286.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 54634434.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 67724910.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 64948994.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9 - 47505656.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 25653736.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - &
     & 9448152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 1788168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 199872.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**4 - 192384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 23040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 4608.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 15.d0*gl%z3_PCSAFT(1)**20 + 285.d0*gl%z3_PCSAFT(1)**19 - 2305.d0*gl%z3_PCSAFT(1)**18 &
     & + 9475.d0*gl%z3_PCSAFT(1)**17 - 12047.d0*gl%z3_PCSAFT(1)**16 - 86543.d0*gl%z3_PCSAFT(1)**15 + 621905.d0* &
     & gl%z3_PCSAFT(1)**14 - 2244223.d0*gl%z3_PCSAFT(1)**13 + 5570442.d0*gl%z3_PCSAFT(1)**12 - 10346082.d0* &
     & gl%z3_PCSAFT(1)**11 + 14876292.d0*gl%z3_PCSAFT(1)**10 - 16772000.d0*gl%z3_PCSAFT(1)**9 + 14806880.d0* &
     & gl%z3_PCSAFT(1)**8 - 10083616.d0*gl%z3_PCSAFT(1)**7 + 5123072.d0*gl%z3_PCSAFT(1)**6 - 1813904.d0* &
     & gl%z3_PCSAFT(1)**5 + 377440.d0*gl%z3_PCSAFT(1)**4 - 15072.d0*gl%z3_PCSAFT(1)**3 - 12096.d0*gl%z3_PCSAFT(1)**2 &
     & + 2112.d0*gl%z3_PCSAFT(1))
        part3 = (24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 1464.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**18 + 25160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 227264.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**16 + 1321810.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 5488510.d0*gl%mmean_PCSAFT**3* &
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
        do xi = 1, gl%ncomp
            gl%cx1_PCSAFT(3,xi) = part1*(gl%mPCSAFT(xi)*part2 + gl%z3x1_PCSAFT(1,xi)*part3)
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF C_1 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires z3_PCSAFT(1), z3_PCSAFT(4)
    if (GETDERC(4) .eq. 1) then
        part1 = 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**2/recurring_factor**3
        part2 = (2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 74.d0*gl%mmean_PCSAFT* &
                    & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 905.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 5835.d0*gl%mmean_PCSAFT* &
                    & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 23098.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 59594.d0* &
                    & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 101491.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - &
                    & 110777.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 69804.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) &
                    & - 16584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 5244.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) &
                    & + 2520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 5.d0* &
                    & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 65.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 344.d0*gl%z3_PCSAFT(1)**11* &
                    & gl%z3_PCSAFT(4) - 818.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 27.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 5879.d0 &
                    & *gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 18410.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 30864.d0*gl%z3_PCSAFT(1)**6* &
                    & gl%z3_PCSAFT(4) - 32152.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 20876.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - &
                    & 7600.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 864.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 336.d0*gl%z3_PCSAFT(1)* &
                    & gl%z3_PCSAFT(4) - 96.d0*gl%z3_PCSAFT(4))
        part3 = (-2.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + &
                    & 74.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 905.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 5835.d0*gl%mmean_PCSAFT**2* &
                    & gl%z3_PCSAFT(1)**10 - 23098.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 59594.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 &
                    & - 101491.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 110777.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 69804.d0* &
                    & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 16584.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 5244.d0*gl%mmean_PCSAFT**2* &
                    & gl%z3_PCSAFT(1)**3 - 2520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 288.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - gl%mmean_PCSAFT &
                    & *gl%z3_PCSAFT(1)**13 - 23.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 509.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 4043.d0* &
                    & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 18099.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 51547.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 &
                    & + 97345.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 121767.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 96236.d0*gl%mmean_PCSAFT* &
                    & gl%z3_PCSAFT(1)**5 - 41068.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 3068.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 4392.d0* &
                    & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 1104.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 96.d0*gl%mmean_PCSAFT + gl%z3_PCSAFT(1)**13 - 7.d0* &
                    & gl%z3_PCSAFT(1)**12 - 26.d0*gl%z3_PCSAFT(1)**11 + 532.d0*gl%z3_PCSAFT(1)**10 - 3083.d0*gl%z3_PCSAFT(1)**9 + &
                    & 10173.d0*gl%z3_PCSAFT(1)**8 - 21708.d0*gl%z3_PCSAFT(1)**7 + 31074.d0*gl%z3_PCSAFT(1)**6 - 29592.d0* &
                    & gl%z3_PCSAFT(1)**5 + 17588.d0*gl%z3_PCSAFT(1)**4 - 5144.d0*gl%z3_PCSAFT(1)**3 - 400.d0*gl%z3_PCSAFT(1)**2 + &
                    & 752.d0*gl%z3_PCSAFT(1) - 160.d0)
        part4 = (6.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - &
                    & 252.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 2971.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
                    & gl%z3_PCSAFT(4) - 17812.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 65442.d0*gl%mmean_PCSAFT**2* &
                    & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 160048.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 269557.d0* &
                    & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 311268.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) &
                    & + 233112.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 95280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3* &
                    & gl%z3_PCSAFT(4) + 8172.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 5040.d0*gl%mmean_PCSAFT**2* &
                    & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 576.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) + 3.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12* &
                    & gl%z3_PCSAFT(4) + 114.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 1993.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10* &
                    & gl%z3_PCSAFT(4) + 13636.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 53901.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 &
                    & *gl%z3_PCSAFT(4) + 139874.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 251789.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
                    & **6*gl%z3_PCSAFT(4) + 319716.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 282112.d0*gl%mmean_PCSAFT* &
                    & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 163912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 55140.d0* &
                    & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 7152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 528.d0*gl%mmean_PCSAFT &
                    & *gl%z3_PCSAFT(4) - 3.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 6.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 256.d0* &
                    & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 2392.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 10737.d0*gl%z3_PCSAFT(1)**8* &
                    & gl%z3_PCSAFT(4) - 30390.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 59442.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - &
                    & 83304.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 84096.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 59928.d0* &
                    & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 28632.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 8224.d0*gl%z3_PCSAFT(1)* &
                    & gl%z3_PCSAFT(4) + 1072.d0*gl%z3_PCSAFT(4))
        do xi = 1, gl%ncomp
            gl%cx1_PCSAFT(4,xi) = part1*(gl%mPCSAFT(xi)*part2 + gl%z3x1_PCSAFT(4,xi)*part3 + gl%z3x1_PCSAFT(1,xi)*part4)
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF C_1 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2 
    ! requires z3_PCSAFT(1), z3_PCSAFT(4), z3_PCSAFT(5)
    if (GETDERC(5) .eq. 1) then
        part1 = 2.d0*(gl%z3_PCSAFT(1) - 1.d0)/recurring_factor**4
        part2 = (4.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 12.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 184.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 612.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 3276.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(5) - 10880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 33160.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 106712.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + &
     & 221009.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 681631.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4)**2 - 1042610.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 3075481.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + 3618146.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(5) - 10251739.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 9414730.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 25885565.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)**2 + 18461536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 50074577.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 27049178.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) + 74150471.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 28830270.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 82775761.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4)**2 - 21010586.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 67043103.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 8825283.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(5) - 36033036.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 528204.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 9665124.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4)**2 - 1383180.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 1146708.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + 487260.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(5) - 1439964.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 53064.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 117792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - &
     & 34560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 69120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)**2 - 3456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 6912.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 8.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 24.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 264.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 552.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + &
     & 776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 3920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(5) - 33272.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - 66462.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 305506.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + &
     & 448664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 1642014.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)**2 - 1923528.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 6079314.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + 5830296.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - &
     & 16549086.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 13008588.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 &
     & *gl%z3_PCSAFT(5) + 34220286.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 21647344.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 54634434.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)** &
     & 2 - 26719192.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 67724910.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 9*gl%z3_PCSAFT(4)**2 + 23815056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 64948994.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 14331758.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) &
     & + 47505656.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + 4754232.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(5) - 25653736.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 66712.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 9448152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - &
     & 733176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 1788168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4)**2 + 213744.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 199872.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 13824.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 192384.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 11520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - &
     & 23040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 4608.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**2 - 5.d0* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 15.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + 110.d0*gl%z3_PCSAFT(1)**19 &
     & *gl%z3_PCSAFT(5) - 285.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 1104.d0*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(5) + 2305.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 6534.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(5) - 9475.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - 24200.d0*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(5) + 12047.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + 50794.d0*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(5) + 86543.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 9238.d0*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(5) - 621905.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 339894.d0*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(5) + 2244223.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + 1342133.d0*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(5) - 5570442.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 3065896.d0*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) + 10346082.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 4858318.d0*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(5) - 14876292.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 5572200.d0*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(5) + 16772000.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 4605528.d0*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(5) - 14806880.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 2609944.d0*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(5) + 10083616.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 854504.d0*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(5) - 5123072.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 17696.d0*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(5) + 1813904.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 120912.d0*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(5) - 377440.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 49120.d0*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(5) + 15072.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 4608.d0*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(5) + 12096.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 1728.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - &
     & 2112.d0*gl%z3_PCSAFT(4)**2 + 384.d0*gl%z3_PCSAFT(5))
        part3 = (-4.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**20 + 184.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 3276.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + &
     & 33160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 221009.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 1042610.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 3618146.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 9414730.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**13 - 18461536.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + 27049178.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**11 - 28830270.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + 21010586.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**9 - 8825283.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 528204.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 7 + 1383180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 487260.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - &
     & 53064.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 34560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 3456.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**2 - 120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 3006.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 18 - 35376.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 258291.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - &
     & 1302792.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 4786440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - &
     & 13162656.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 27403146.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - &
     & 43125420.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 50491326.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - &
     & 42283596.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 23007675.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - &
     & 5669316.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 1656420.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 1633356.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 289464.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 84000.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**3 + 23616.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 2304.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + 3.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 18.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + &
     & 10302.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 90360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + 506538.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**15 - 2008806.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 5898522.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 &
     & - 13100283.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 22146504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - &
     & 28304742.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 26701968.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 17554992.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 6825240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 369576.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 6 - 1115328.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 519216.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 38880.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 31104.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 5952.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 384.d0 &
     & *gl%mmean_PCSAFT - gl%z3_PCSAFT(1)**20 + 16.d0*gl%z3_PCSAFT(1)**19 - 72.d0*gl%z3_PCSAFT(1)**18 - 452.d0*gl%z3_PCSAFT(1) &
     & **17 + 8230.d0*gl%z3_PCSAFT(1)**16 - 57832.d0*gl%z3_PCSAFT(1)**15 + 259468.d0*gl%z3_PCSAFT(1)**14 - &
     & 832112.d0*gl%z3_PCSAFT(1)**13 + 1994375.d0*gl%z3_PCSAFT(1)**12 - 3639048.d0*gl%z3_PCSAFT(1)**11 + &
     & 5069444.d0*gl%z3_PCSAFT(1)**10 - 5327740.d0*gl%z3_PCSAFT(1)**9 + 4081916.d0*gl%z3_PCSAFT(1)**8 - &
     & 2088256.d0*gl%z3_PCSAFT(1)**7 + 509056.d0*gl%z3_PCSAFT(1)**6 + 139664.d0*gl%z3_PCSAFT(1)**5 - 163504.d0* &
     & gl%z3_PCSAFT(1)**4 + 49408.d0*gl%z3_PCSAFT(1)**3 + 448.d0*gl%z3_PCSAFT(1)**2 - 3648.d0*gl%z3_PCSAFT(1) + 640.d0)
        part4 = (24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 1224.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 21760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 213424.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 1363262.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) - 6150962.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 20503478.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 51771130.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + &
     & 100149154.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 148300942.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 165551522.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - &
     & 134086206.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 72066072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **7*gl%z3_PCSAFT(4) - 19330248.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 2293416.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 2879928.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) &
     & - 235584.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 138240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) - 13824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 18*gl%z3_PCSAFT(4) - 22572.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 248412.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 1695630.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + &
     & 8020962.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 27831486.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(4) + 73116642.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 147870714.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 231389334.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) - 278365758.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 251992122.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 163672536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) + 68090952.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 11559576.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 3828120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 2222592.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 182592.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) - 48384.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 4608.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(4) - 18.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 42.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) + 6234.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 87966.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 16*gl%z3_PCSAFT(4) + 669414.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 3384186.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 12352134.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 33974010.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 72073236.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4) - 119208756.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 153791928.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 153096432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + &
     & 114493248.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 60709440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4) + 19625856.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 1474464.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 1684800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 691008.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 78720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 4224.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(4) + 6.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 66.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - &
     & 194.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 8558.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 81214.d0*gl%z3_PCSAFT(1) &
     & **15*gl%z3_PCSAFT(4) + 455842.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 1782118.d0*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) + 5184266.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 11605040.d0*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4) + 20338040.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 28078000.d0*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) + 30423344.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 25474688.d0*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) + 15927424.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 6880960.d0*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4) + 1603072.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 135040.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) &
     & - 239104.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 74368.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 8576.d0* &
     & gl%z3_PCSAFT(4))
        part5 = (12.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 48.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 612.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(5) + 2688.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 10880.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 46920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - &
     & 106712.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 440688.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4)**2 + 681631.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 2685072.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 3075481.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(5) + 11639472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + 10251739.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 37820544.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)**2 - 25885565.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + 94869936.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 50074577.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) - 186177984.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 74150471.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 285790752.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4)**2 + 82775761.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 338571768.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 67043103.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(5) + 300527040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + 36033036.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 188963616.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4)**2 - 9665124.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 74766816.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 1146708.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(5) - 12735072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 1439964.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 1615248.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4)**2 - 117792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 353376.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 69120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + &
     & 207360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 6912.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(5) + 20736.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4)**2 + 480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(5) - 2400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 11286.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 52764.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + &
     & 124206.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 540744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4)**2 - 847815.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 3450312.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + 4010481.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(5) - 15426216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 13915743.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 51451488.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)**2 + 36558321.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 132643272.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 73935357.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) + 269219784.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 115694667.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 432649560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4)**2 - 139182879.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 547836564.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 125996061.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(5) - 538435344.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 81836268.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 399847440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4)**2 + 34045476.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 214506480.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 5779788.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(5) + 76749456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 1914060.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 15374064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4)**2 + 1111296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 810720.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 91296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + &
     & 156096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 24192.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(5) + 3456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**2 - 2304.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(5) - &
     & 9.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 36.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - &
     & 21.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 &
     & + 3117.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 17628.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4)**2 - 43983.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 211272.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + 334707.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - &
     & 1446000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 1692093.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(5) + 6750024.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + 6176067.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 23293104.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - &
     & 16987005.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + 62075208.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4)**2 + 36036618.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 130835460.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 59604378.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(5) + 220443816.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 76895964.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 297205764.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 &
     & - 76548216.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 318374160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4)**2 + 57246624.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 266980560.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 30354720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) &
     & + 171031920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + 9812928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 &
     & *gl%z3_PCSAFT(5) - 80498256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 737232.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 26049984.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - &
     & 842400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 5082240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)**2 + 345504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 408000.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 39360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 14208.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(4)**2 - 2112.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(5) + 3.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 12.d0* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 33.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 72.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4)**2 - 97.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 1488.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
     & **2 + 4279.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 25536.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - &
     & 40607.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 196056.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + &
     & 227921.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 973392.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - &
     & 891059.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 3511248.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + &
     & 2592133.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 9735552.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - &
     & 5802520.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + 21385524.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + &
     & 10169020.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 37787352.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - &
     & 14039000.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 53984160.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + &
     & 15211672.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 62208192.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - &
     & 12737344.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 57318240.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + &
     & 7963712.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 41585760.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - &
     & 3440480.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 23198496.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + &
     & 801536.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 9589056.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + &
     & 67520.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 2763072.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - &
     & 119552.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 494976.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 37184.d0* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 41472.d0*gl%z3_PCSAFT(4)**2 - 4288.d0*gl%z3_PCSAFT(5))
        do xi = 1, gl%ncomp 
            gl%cx1_PCSAFT(5,xi) = part1*(gl%mPCSAFT(xi)*part2 + gl%z3x1_PCSAFT(5,xi)*part3 + gl%z3x1_PCSAFT(4,xi)*part4 + gl%z3x1_PCSAFT(1,xi)*part5)
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF C_1 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires z3_PCSAFT(1), z3_PCSAFT(4)
    if (GETDERC(6) .eq. 1) then
        part1 = -2.d0*(gl%z3_PCSAFT(1) - 1.d0)/recurring_factor**4
        part2 = (8.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 428.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 7604.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 73552.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 460622.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) - 2032871.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 6633593.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 16470835.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 31613041.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 47101293.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) + 53945491.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 46032517.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 27207753.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4) - 9136920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 236472.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 952704.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 170856.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 34560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - &
     & 3456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 16.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - &
     & 136.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 1328.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + &
     & 29352.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 239044.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 1193350.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 4155786.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 10718790.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 21211698.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 32987090.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4) - 41005718.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 41133938.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 33173898.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 20899504.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 9514864.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + &
     & 2521344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 13872.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) &
     & - 206208.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 34560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) + 4608.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 10.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + &
     & 175.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 1201.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 2941.d0*gl%z3_PCSAFT(1) &
     & **17*gl%z3_PCSAFT(4) + 12153.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 137337.d0*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4) + 631143.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 1904329.d0*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) + 4228309.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 7280186.d0*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4) + 10017974.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 11199800.d0*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4) + 10201352.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 7473672.d0*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) + 4268568.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 1796208.d0*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4) + 498352.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 64192.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - &
     & 7488.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 3840.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 384.d0*gl%z3_PCSAFT(4))
        part3 = (-8.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 428.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - &
     & 7604.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 73552.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 460622.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 2032871.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 6633593.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**14 + 16470835.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 31613041.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**12 + 47101293.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 53945491.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**10 + 46032517.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 - 27207753.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**8 + 9136920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - 236472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 6 - 952704.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 170856.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 34560.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 3456.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 360.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**19 + 8280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 88830.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 &
     & + 589524.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 2707689.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + &
     & 9129303.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 23395665.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + &
     & 46532211.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 72569247.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + &
     & 88691553.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 83712465.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + &
     & 58828593.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 28376160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + &
     & 7436208.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 280704.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 821832.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 175296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 576.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**2 + 6.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 39.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 2577.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**18 + 33681.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 244347.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + &
     & 1185555.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 4167261.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + 11088483.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 22936335.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + 37457874.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11 - 48591222.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + 49846248.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 &
     & - 39691632.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + 23529480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 9443352.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 1852560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 323184.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 4 - 306624.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 70464.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 3840.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1) - 384.d0*gl%mmean_PCSAFT - 2.d0*gl%z3_PCSAFT(1)**20 + 17.d0*gl%z3_PCSAFT(1)**19 + 169.d0*gl%z3_PCSAFT(1)** &
     & 18 - 3827.d0*gl%z3_PCSAFT(1)**17 + 32377.d0*gl%z3_PCSAFT(1)**16 - 170089.d0*gl%z3_PCSAFT(1)**15 + &
     & 631591.d0*gl%z3_PCSAFT(1)**14 - 1760021.d0*gl%z3_PCSAFT(1)**13 + 3808145.d0*gl%z3_PCSAFT(1)**12 - &
     & 6529972.d0*gl%z3_PCSAFT(1)**11 + 8969556.d0*gl%z3_PCSAFT(1)**10 - 9883932.d0*gl%z3_PCSAFT(1)**9 + &
     & 8655428.d0*gl%z3_PCSAFT(1)**8 - 5875456.d0*gl%z3_PCSAFT(1)**7 + 2931424.d0*gl%z3_PCSAFT(1)**6 - &
     & 941200.d0*gl%z3_PCSAFT(1)**5 + 95984.d0*gl%z3_PCSAFT(1)**4 + 70144.d0*gl%z3_PCSAFT(1)**3 - 37632.d0* &
     & gl%z3_PCSAFT(1)**2 + 7936.d0*gl%z3_PCSAFT(1) - 640.d0)
        part4 = (24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 19*gl%z3_PCSAFT(4) - 1464.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 25160.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 227264.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + &
     & 1321810.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 5488510.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 14*gl%z3_PCSAFT(4) + 17317066.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 43098806.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 86028830.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4) - 137489810.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 173020246.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 166440834.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4) + 116897544.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 55436568.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 15028488.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 1264680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 117792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) - 69120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 6912.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 1440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 30192.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 292332.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - &
     & 1754682.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 7405254.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 14*gl%z3_PCSAFT(4) - 23620002.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 59526630.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 121349070.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4) + 201260226.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 269470806.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 286443222.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4) - 236174904.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 146415528.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 65189880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4) + 19202184.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 3033312.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 26496.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 44928.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) - 18.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 342.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 11394.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 123306.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 776586.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 3365838.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + &
     & 10940970.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 28101198.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4) + 58762224.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 101235060.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 143413836.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - &
     & 165277728.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 152487312.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) - 110322480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 60872400.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 24575520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 6767040.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 1099008.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + &
     & 64512.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4224.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) + 6.d0*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 6.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 1294.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + &
     & 16978.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 114842.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 517550.d0* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 1729130.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 4551286.d0* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 9780484.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 17449312.d0* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 25906160.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 31784848.d0* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 31843552.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 25658336.d0* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 16317536.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 7985984.d0*gl%z3_PCSAFT(1) &
     & **4*gl%z3_PCSAFT(4) - 2898112.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 734080.d0*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) - 115840.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 8576.d0*gl%z3_PCSAFT(4))
        do xi = 1, gl%ncomp
            gl%cx1_PCSAFT(6,xi) = part1*(gl%mPCSAFT(xi)*part2 + gl%z3x1_PCSAFT(4,xi)*part3 + gl%z3x1_PCSAFT(1,xi)*part4)
        end do
    end if

    
    ! 7: 3RD MIXED DERIVATIVE OF C_1 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires z3_PCSAFT(1), z3_PCSAFT(4), z3_PCSAFT(5)
    if (GETDERC(7) .eq. 1) then
        part1 = -2.d0/recurring_factor**5
        part2 = (16.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(5) - 48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 &
     & *gl%z3_PCSAFT(4)**2 - 1000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 3360.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + 23448.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - &
     & 78280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 313620.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 24*gl%z3_PCSAFT(5) + 1008400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 2807212.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 8599564.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4)**2 - 18246642.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) + 53183016.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 90364242.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(5) - 251949718.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 351527431.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 947978420.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4)**2 + 1095821600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - &
     & 2903511910.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 2772663215.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 7355969000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 &
     & + 5735183160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 15545965540.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - 9716828902.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(5) + 27443029848.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + &
     & 13430529548.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 40241827400.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 14964967460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(5) + 48411772920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + &
     & 13117667042.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 46812356614.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 8615783839.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) &
     & + 35222081892.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 3786214464.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 19496693406.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)**2 - 703782987.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + &
     & 7034718840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 292610820.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 1032905880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + &
     & 216928512.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 326883744.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 29751048.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + &
     & 164897064.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 11803752.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 11074752.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + &
     & 2233440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 1779840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 &
     & *gl%z3_PCSAFT(4)**2 + 466560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 933120.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 41472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - &
     & 82944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 24.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27* &
     & gl%z3_PCSAFT(5) - 72.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)**2 - 60.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - &
     & 10872.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 55080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4)**2 + 229080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 929400.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 2496198.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(5) + 8885586.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 18110109.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 58022112.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4)**2 - 96140838.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 281620122.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 392293269.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(5) - 1067875380.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - &
     & 1267303170.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 3270736260.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 3303484746.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) &
     & - 8283476832.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 7034682828.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 17623834884.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
     & **2 + 12328763754.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 31784676576.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 17836871250.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(5) + 48721323810.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + &
     & 21265685565.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 63292400880.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 20724474366.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(5) + 69161009562.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + &
     & 16216090629.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 62826035652.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 9827320686.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) &
     & + 46624825104.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 4261587864.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 27479320320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4)**2 - 1029848760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + &
     & 12184201920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 90875868.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 3558064944.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + &
     & 173814768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 357911616.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 53388960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + &
     & 182931840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + 896544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(5) - 69824160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 2858400.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 1134720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4)**2 - 359424.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 1963008.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 41472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + &
     & 193536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 36.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27* &
     & gl%z3_PCSAFT(5) + 108.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)**2 + 810.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 26*gl%z3_PCSAFT(5) - 1800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - 6678.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 3090.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + 4995.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 188400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 &
     & + 407976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 2643642.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4)**2 - 4587969.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) + 19929552.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 29448336.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - &
     & 103701504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 133953354.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 20*gl%z3_PCSAFT(5) + 408193380.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + 465275880.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 1279009680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4)**2 - 1281207036.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 3301443672.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 2860377330.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(5) - 7194693150.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - 5254199289.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 13455239256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4)**2 + 8019619200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 21769063530.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 10231558785.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(5) + 30470009520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + &
     & 10927205940.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 36648154260.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **12*gl%z3_PCSAFT(4)**2 - 9719912256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + &
     & 37426323108.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 7088077476.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 31929173784.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - &
     & 4089330540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 22282843920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **9*gl%z3_PCSAFT(4)**2 + 1721345520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - &
     & 12368552640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 409013328.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **8*gl%z3_PCSAFT(5) + 5239619856.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 39299616.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 1576584192.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) &
     & **2 + 73052496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 281856576.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(4)**2 - 23768448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 5062080.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 1117440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - &
     & 11702400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 1077888.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(5) + 2899584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 160128.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 170496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 13824.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 36864.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**2 + 10.d0*gl%z3_PCSAFT(1)**27* &
     & gl%z3_PCSAFT(5) - 30.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)**2 - 265.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) &
     & + 660.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + 3126.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 5530.d0* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 20565.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 12100.d0* &
     & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 69106.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 175208.d0* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 48750.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 2102940.d0* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 2003424.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 12887756.d0* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 12999494.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 55221520.d0 &
     & *gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 53688970.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + &
     & 182536130.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 165206811.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) &
     & - 490224532.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 401380802.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(5) + 1105946846.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + 793227527.d0*gl%z3_PCSAFT(1) &
     & **16*gl%z3_PCSAFT(5) - 2141896748.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 1297925122.d0* &
     & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 3602556660.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + &
     & 1777365840.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 5275006820.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) &
     & **2 - 2047341652.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 6693171824.d0*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)**2 + 1981704872.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 7293792776.d0*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4)**2 - 1597441568.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + 6750308992.d0* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 1050521584.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - &
     & 5238668960.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 540924400.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) &
     & + 3358744000.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 198986592.d0*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(5) - 1745296544.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 37752608.d0*gl%z3_PCSAFT(1)**7 &
     & *gl%z3_PCSAFT(5) + 714947904.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 7591424.d0*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(5) - 220560064.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + 8296832.d0*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(5) + 46732800.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 2732480.d0*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(5) - 5116160.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 352768.d0*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(5) - 281856.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 35328.d0*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(5) + 170496.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 16896.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - &
     & 16896.d0*gl%z3_PCSAFT(4)**2 + 1536.d0*gl%z3_PCSAFT(5))
        part3 = (-16.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**27 + 1000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26 - 23448.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25 &
     & + 313620.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 - 2807212.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 + &
     & 18246642.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 - 90364242.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 + &
     & 351527431.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 - 1095821600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 + &
     & 2772663215.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 - 5735183160.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 + &
     & 9716828902.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 - 13430529548.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 &
     & + 14964967460.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 - 13117667042.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 13 + 8615783839.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 - 3786214464.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **11 + 703782987.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 + 292610820.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **9 - 216928512.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 + 29751048.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 &
     & + 11803752.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 - 2233440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 - &
     & 466560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - 41472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 + 8.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**27 - 1220.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + 34776.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 25 - 508340.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 4795114.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - &
     & 32321997.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 164706564.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - &
     & 656906812.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 2097393910.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - &
     & 5440998823.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 11573850044.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 &
     & - 20274139992.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 29225546750.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 15 - 34432325155.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 32641175388.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**13 - 24148643532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + 13095571138.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**11 - 4406466157.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + 238099340.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**9 + 582000084.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - 265953384.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**7 + 25358760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 12446208.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**5 - 2567520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 131328.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 3 - 13824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 + 12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27 + 330.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**26 - 16974.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 291585.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**24 - 2975382.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 21059748.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**22 - 111205302.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 456751023.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**20 - 1498111950.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 3992387952.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**18 - 8742162690.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 15830834583.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**16 - 23754627810.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 29447075160.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 29869230330.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + &
     & 24308313837.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 15273858102.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 &
     & + 6805123770.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 1613062320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 &
     & - 259332324.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 377279232.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - &
     & 127467552.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 8103456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + &
     & 6300000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 1524096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 49536.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 4608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 10.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 + &
     & 85.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 + 2616.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - 67140.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **24 + 782666.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 5941965.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + &
     & 32862936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 139881596.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + &
     & 473260930.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 1299237369.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + &
     & 2933689208.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 5495114348.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + &
     & 8574673798.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 11145121875.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 12001604728.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 10570972388.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
     & 7425908312.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 3951946396.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
     & 1397041120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 155930688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - &
     & 152260288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 101185136.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 25990208.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 117440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 1705088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 3 - 392832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 21504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 1536.d0*gl%mmean_PCSAFT + 2.d0* &
     & gl%z3_PCSAFT(1)**27 - 35.d0*gl%z3_PCSAFT(1)**26 + 54.d0*gl%z3_PCSAFT(1)**25 + 4615.d0*gl%z3_PCSAFT(1)**24 - &
     & 71426.d0*gl%z3_PCSAFT(1)**23 + 605880.d0*gl%z3_PCSAFT(1)**22 - 3570708.d0*gl%z3_PCSAFT(1)**21 + &
     & 15896174.d0*gl%z3_PCSAFT(1)**20 - 55790690.d0*gl%z3_PCSAFT(1)**19 + 158375661.d0*gl%z3_PCSAFT(1)**18 &
     & - 369691514.d0*gl%z3_PCSAFT(1)**17 + 717213891.d0*gl%z3_PCSAFT(1)**16 - 1163621030.d0*gl%z3_PCSAFT(1) &
     & **15 + 1582094710.d0*gl%z3_PCSAFT(1)**14 - 1798178232.d0*gl%z3_PCSAFT(1)**13 + 1694152872.d0 &
     & *gl%z3_PCSAFT(1)**12 - 1300126312.d0*gl%z3_PCSAFT(1)**11 + 785162152.d0*gl%z3_PCSAFT(1)**10 - &
     & 345574560.d0*gl%z3_PCSAFT(1)**9 + 85916320.d0*gl%z3_PCSAFT(1)**8 + 10458528.d0*gl%z3_PCSAFT(1)**7 - &
     & 20644832.d0*gl%z3_PCSAFT(1)**6 + 9101440.d0*gl%z3_PCSAFT(1)**5 - 1665600.d0*gl%z3_PCSAFT(1)**4 - &
     & 181248.d0*gl%z3_PCSAFT(1)**3 + 165632.d0*gl%z3_PCSAFT(1)**2 - 34304.d0*gl%z3_PCSAFT(1) + 2560.d0)
        part4 = (96.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 6720.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 156560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - &
     & 2016800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 17199128.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 22*gl%z3_PCSAFT(4) - 106366032.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 503899436.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 1895956840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 5807023820.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 14711938000.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 31091931080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) - 54886059696.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 80483654800.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 96823545840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) + 93624713228.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 70444163784.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 38993386812.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4) - 14069437680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 2065811760.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 653767488.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4) - 329794128.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 22149504.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 3559680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 1866240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 165888.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) - 48.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 9120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **25*gl%z3_PCSAFT(4) - 250960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 3462400.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 30637756.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4) + 193898392.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 933655472.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 3559661880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 11042776420.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 28391782272.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 61170438144.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) + 110912061136.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - &
     & 168904712020.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 214371069880.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 223727961712.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) &
     & + 188151835832.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 123609020364.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 60136404880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
     & - 19252650960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 2466224064.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 904393584.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - &
     & 502992000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 86126400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **4*gl%z3_PCSAFT(4) - 1059840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 718848.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 55296.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 72.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 3600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) &
     & + 141780.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 2173800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 23*gl%z3_PCSAFT(4) + 20234928.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 131986968.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 648942876.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) - 2515533600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 7924485840.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 20715972288.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) + 45549383580.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 84817377384.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 133872461760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4) - 178310582040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 198588474420.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 182352434832.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 135234993576.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) &
     & - 78424068480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 33514212480.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 9093579264.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 582381408.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 678226176.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(4) - 299926080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 51448320.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 1287936.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & - 423936.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 36864.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) + 60.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 31840.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 585400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - &
     & 5856364.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 39750960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) - 200708608.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 794206880.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 2548309660.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + &
     & 6787845656.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 15248193568.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 16*gl%z3_PCSAFT(4) + 29152447864.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 47586091860.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 66184045120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) - 77926067392.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 76854382768.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 62521390496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) &
     & + 41004768640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 20877272960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **8*gl%z3_PCSAFT(4) + 7644013312.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 1588459392.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 106751488.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 208385280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 72627200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) + 12050688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 688128.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 33792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) - 12.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + &
     & 120.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 2060.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 56000.d0*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4) + 624224.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 4468160.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) + 23331784.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 94680000.d0*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 310436180.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 844369704.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) + 1939631004.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 3804869728.d0*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4) + 6406089880.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 9255938320.d0*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(4) + 11425790608.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 11952202592.d0* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 10469313216.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - &
     & 7548343360.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 4362769920.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - &
     & 1928017920.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 583173632.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - &
     & 72694784.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 31962880.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + &
     & 21867520.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 6353408.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 995328.d0 &
     & *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 68608.d0*gl%z3_PCSAFT(4))
        part5 = (48.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 192.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - 3360.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 14880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4)**2 + 78280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 338400.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 1008400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(5) + 4144800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 8599564.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 33448992.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4)**2 - 53183016.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 197233464.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 251949718.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(5) - 905668080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 947978420.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 3375203880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4)**2 + 2903511910.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - &
     & 10466420640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 7355969000.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 27326586000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
     & **2 + 15545965540.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 60198351360.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 27443029848.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(5) + 111443220960.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + &
     & 40241827400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 172053489600.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 48411772920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(5) + 219210585000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + &
     & 46812356614.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 227279145072.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 35222081892.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) + 188030691048.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + &
     & 19496693406.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 120578018400.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 7034718840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) &
     & + 57280613040.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 1032905880.d0*gl%mmean_PCSAFT** &
     & 4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 18684187200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)** &
     & 2 + 326883744.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 3638493648.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 164897064.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - &
     & 331033824.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + 11074752.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 31687200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + &
     & 1779840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 5339520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 &
     & *gl%z3_PCSAFT(4)**2 + 933120.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 2799360.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 82944.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(5) - 248832.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 24.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 96.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + 4560.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 21840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4)**2 - 125480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 571680.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 1731200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(5) - 7371600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 - 15318878.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) + 60869424.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4)**2 + 96949196.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - 362638884.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 466827736.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(5) + 1672674360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + &
     & 1779830940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 6252506460.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 5521388210.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) &
     & + 19488518640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 14195891136.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 51406987128.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4)**2 - 30585219072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + &
     & 115287877968.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + 55456030568.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 219420284760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)**2 - 84452356010.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + &
     & 352580495280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + 107185534940.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 475025446020.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)**2 - 111863980856.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + &
     & 532062358968.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 94075917916.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 490120074204.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)**2 - 61804510182.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + &
     & 365938014960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 30068202440.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 216897960240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) &
     & **2 - 9626325480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 98879182560.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + 1233112032.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(5) - 32860946016.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + &
     & 452196792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 7131124512.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 251496000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - &
     & 713665440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 43063200.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 42128640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - &
     & 529920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 5149440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)**2 - 359424.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 2944512.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 27648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 248832.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4)**2 - 36.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 144.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - 1800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(5) + 10440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + 70890.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 351000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - &
     & 1086900.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 4869000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 22*gl%z3_PCSAFT(4)**2 + 10117464.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - &
     & 41499072.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 65993484.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 251110836.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 &
     & + 324471438.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 1167970320.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 1257766800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) &
     & + 4392798660.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 3962242920.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - 13797897120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) &
     & **2 - 10357986144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 36849156912.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + 22774691790.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 &
     & *gl%z3_PCSAFT(5) - 84257651160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - &
     & 42408688692.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 164939983920.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + 66936230880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(5) - 275412556800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - &
     & 89155291020.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 390223898100.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + 99294237210.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(5) - 466262553600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - &
     & 91176217416.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + 466292247444.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 67617496788.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(5) - 386453962800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - &
     & 39212034240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 261800449920.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 16757106240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) &
     & - 142106120640.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 4546789632.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 59971499856.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) &
     & **2 + 291190704.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 18749503296.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + 339113088.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(5) + 3983431680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 149963040.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 473328000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4)**2 + 25724160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 11623680.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 643968.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(5) + 2271744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 211968.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 41472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**2 - 18432.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(5) + 30.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(4)**2 + 60.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 1500.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 24*gl%z3_PCSAFT(4)**2 - 15920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 91440.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 292700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - &
     & 1412700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 - 2928182.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(5) + 12571176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 19875480.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - 77740440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)** &
     & 2 - 100354304.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 366205080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **19*gl%z3_PCSAFT(4)**2 + 397103440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - &
     & 1390328160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 1274154830.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **18*gl%z3_PCSAFT(5) + 4411117560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + &
     & 3393922828.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 11943375564.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **16*gl%z3_PCSAFT(4)**2 - 7624096784.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + &
     & 27850707936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + 14576223932.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 56007931260.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - &
     & 23793045930.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 96862530840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **13*gl%z3_PCSAFT(4)**2 + 33092022560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - &
     & 143424237120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 38963033696.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + 180886943928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 &
     & + 38427191384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 193157927496.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 31260695248.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + &
     & 173342256000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 20502384320.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 129451863360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - &
     & 10438636480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 79374582240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **7*gl%z3_PCSAFT(4)**2 + 3822006656.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - &
     & 39214371648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 794229696.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **6*gl%z3_PCSAFT(5) + 15191025024.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - &
     & 53375744.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 4427176320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4)**2 + 104192640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 905917440.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 36313600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) &
     & - 113556480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 6025344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 &
     & *gl%z3_PCSAFT(5) + 5803008.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 344064.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 170496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**2 - 16896.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(5) - 6.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 24.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + &
     & 60.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 60.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + 1030.d0*gl%z3_PCSAFT(1) &
     & **24*gl%z3_PCSAFT(5) - 8280.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 28000.d0*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(5) + 151500.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 312112.d0*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(5) - 1428936.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 2234080.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(5) + 9095280.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 11665892.d0*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(5) - 43590000.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 47340000.d0*gl%z3_PCSAFT(1)**19 &
     & *gl%z3_PCSAFT(5) + 167577480.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 155218090.d0*gl%z3_PCSAFT(1) &
     & **18*gl%z3_PCSAFT(5) - 537998040.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 422184852.d0* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 1477360356.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + &
     & 969815502.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 3509530008.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)** &
     & 2 - 1902434864.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 7231314780.d0*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)**2 + 3203044940.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 12897411000.d0* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 4627969160.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + &
     & 19833976200.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + 5712895304.d0*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(5) - 26182003872.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 5976101296.d0* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + 29523739344.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + &
     & 5234656608.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 28277205120.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
     & **2 - 3774171680.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 22837430880.d0*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4)**2 + 2181384960.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 15403212480.d0*gl%z3_PCSAFT(1) &
     & **7*gl%z3_PCSAFT(4)**2 - 964008960.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 8562933120.d0* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 291586816.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - &
     & 3852722304.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 36347392.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + &
     & 1366982400.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 15981440.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - &
     & 367918080.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 10933760.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + &
     & 70552320.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 3176704.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - &
     & 8583168.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 497664.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 497664.d0* &
     & gl%z3_PCSAFT(4)**2 - 34304.d0*gl%z3_PCSAFT(5))
        do xi = 1, gl%ncomp 
            gl%cx1_PCSAFT(7,xi) = part1*(gl%mPCSAFT(xi)*part2  + gl%z3x1_PCSAFT(5,xi)*part3 + gl%z3x1_PCSAFT(4,xi)*part4 + gl%z3x1_PCSAFT(1,xi)*part5)
        end do
    end if

    ! 8: 3RD DERIVATIVE OF C_1 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires z3_PCSAFT(1)
    if (GETDERC(8) .eq. 1) then
        part1 = 24.d0*gl%z3_PCSAFT(1)**2/recurring_factor**5
        part2 = (8.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 - 520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 &
     & + 12120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 159320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + &
     & 1398006.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 8922926.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + &
     & 43580470.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 168222190.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + &
     & 524072220.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 1334947780.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + &
     & 2798794020.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 4833406720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + &
     & 6840072530.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 7835775650.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + &
     & 7107482826.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 4908348782.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 + &
     & 2387064040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 658940580.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 - &
     & 10515780.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 74341116.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - &
     & 15635484.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 2640600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 444960.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 233280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 20736.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**2 + 12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 + 20.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - &
     & 6590.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 125250.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - 1289384.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 8971512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - 46166390.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 184508210.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 590159500.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 1539898008.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 3316569134.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 5932227610.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - &
     & 8824247500.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 10884130360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 &
     & - 11051143018.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 9123492922.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 11 - 6005042700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 3047436660.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **9 - 1113487540.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 238764536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 7 + 4453584.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 20925840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + &
     & 4919760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 309600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 179712.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 20736.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 18.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 &
     & + 370.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 - 2455.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 7035.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**23 + 258178.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 2428412.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + &
     & 14393825.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 62544505.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 212128750.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 581065778.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 1312271515.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**16 - 2474269265.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 3921621190.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14 - 5236573440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 5879311595.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12 - 5514444623.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 4269499180.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**10 - 2675006980.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 1312348400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **8 - 475323644.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 111677472.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - &
     & 10004560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 2634160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 855840.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**3 - 59136.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 2304.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 5.d0*gl%z3_PCSAFT(1) &
     & **26 - 125.d0*gl%z3_PCSAFT(1)**25 + 1360.d0*gl%z3_PCSAFT(1)**24 - 7880.d0*gl%z3_PCSAFT(1)**23 + 18513.d0 &
     & *gl%z3_PCSAFT(1)**22 + 86525.d0*gl%z3_PCSAFT(1)**21 - 1100860.d0*gl%z3_PCSAFT(1)**20 + 6157810.d0* &
     & gl%z3_PCSAFT(1)**19 - 23953245.d0*gl%z3_PCSAFT(1)**18 + 72111733.d0*gl%z3_PCSAFT(1)**17 - &
     & 175901176.d0*gl%z3_PCSAFT(1)**16 + 356093180.d0*gl%z3_PCSAFT(1)**15 - 605981525.d0*gl%z3_PCSAFT(1)** &
     & 14 + 871872275.d0*gl%z3_PCSAFT(1)**13 - 1061216276.d0*gl%z3_PCSAFT(1)**12 + 1088838786.d0* &
     & gl%z3_PCSAFT(1)**11 - 934830740.d0*gl%z3_PCSAFT(1)**10 + 663966120.d0*gl%z3_PCSAFT(1)**9 - &
     & 383682120.d0*gl%z3_PCSAFT(1)**8 + 176015256.d0*gl%z3_PCSAFT(1)**7 - 61678704.d0*gl%z3_PCSAFT(1)**6 + &
     & 15410320.d0*gl%z3_PCSAFT(1)**5 - 2337440.d0*gl%z3_PCSAFT(1)**4 + 88320.d0*gl%z3_PCSAFT(1)**3 + 34624.d0* &
     & gl%z3_PCSAFT(1)**2 - 4736.d0*gl%z3_PCSAFT(1))
        part3 = (16.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25 - &
     & 1240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 + 28200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 - 345400.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 + 2787416.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 - 16436122.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 + 75472340.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 - 281266990.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 + 872201720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 - 2277215500.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 + 5016529280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 - 9286935080.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 + 14337790800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 - &
     & 18267548750.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 + 18939928756.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 &
     & - 15669224254.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 + 10048168200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 9 - 4773384420.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 + 1557015600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 &
     & - 303207804.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 + 27586152.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 - &
     & 2640600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 + 444960.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 + 233280.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 + 20736.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) - 8.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 25 + 1820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 - 47640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + &
     & 614300.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - 5072452.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + &
     & 30219907.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 139389530.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + &
     & 521042205.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 1624043220.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + &
     & 4283915594.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 9607323164.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + &
     & 18285023730.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 29381707940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 &
     & + 39585453835.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 44338529914.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 11 + 40843339517.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 30494834580.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**9 + 18074830020.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - 8239931880.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**7 + 2738412168.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 594260376.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**5 + 59472120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 3510720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **3 - 429120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 - 245376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) - 20736.d0* &
     & gl%mmean_PCSAFT**3 - 12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 - 870.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + &
     & 29250.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - 405750.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 3458256.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - 20925903.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 97330860.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - 366066555.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 1149824760.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 3070763076.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 7021470930.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 13744998660.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + &
     & 22951046400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 32518658175.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 &
     & + 38855212800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 38857687287.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 10 + 32204496900.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 21816704160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **8 + 11842176720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 4997624988.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **6 + 1562458608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 331952640.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 4 + 39444000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 968640.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - &
     & 189312.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 3456.d0*gl%mmean_PCSAFT**2 + 10.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + &
     & 125.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 7620.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 117725.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22 - 1047598.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 + 6478370.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 &
     & - 30517090.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 115860680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - &
     & 367593130.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + 995281297.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - &
     & 2320892328.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + 4667327605.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - &
     & 8071877570.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 11952019760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - &
     & 15073911994.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 16096493958.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - &
     & 14445188000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 10787655280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - &
     & 6614548520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 3267864304.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - &
     & 1265918752.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 368931360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 75493120.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 9463040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 483584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & - 14208.d0*gl%mmean_PCSAFT - 2.d0*gl%z3_PCSAFT(1)**25 + 5.d0*gl%z3_PCSAFT(1)**24 + 690.d0*gl%z3_PCSAFT(1)**23 - &
     & 12625.d0*gl%z3_PCSAFT(1)**22 + 119078.d0*gl%z3_PCSAFT(1)**21 - 757940.d0*gl%z3_PCSAFT(1)**20 + &
     & 3632500.d0*gl%z3_PCSAFT(1)**19 - 13964790.d0*gl%z3_PCSAFT(1)**18 + 44833170.d0*gl%z3_PCSAFT(1)**17 - &
     & 123113363.d0*gl%z3_PCSAFT(1)**16 + 292460834.d0*gl%z3_PCSAFT(1)**15 - 602609565.d0*gl%z3_PCSAFT(1)** &
     & 14 + 1074784250.d0*gl%z3_PCSAFT(1)**13 - 1652831350.d0*gl%z3_PCSAFT(1)**12 + 2181833656.d0* &
     & gl%z3_PCSAFT(1)**11 - 2460311612.d0*gl%z3_PCSAFT(1)**10 + 2356433760.d0*gl%z3_PCSAFT(1)**9 - &
     & 1903119240.d0*gl%z3_PCSAFT(1)**8 + 1283601040.d0*gl%z3_PCSAFT(1)**7 - 713577760.d0*gl%z3_PCSAFT(1)**6 &
     & + 321060192.d0*gl%z3_PCSAFT(1)**5 - 113915200.d0*gl%z3_PCSAFT(1)**4 + 30659840.d0*gl%z3_PCSAFT(1)**3 &
     & - 5879360.d0*gl%z3_PCSAFT(1)**2 + 715264.d0*gl%z3_PCSAFT(1) - 41472.d0)
        do xi = 1, gl%ncomp 
            gl%cx1_PCSAFT(8,xi) = part1*(gl%mPCSAFT(xi)*part2 + gl%z3x1_PCSAFT(1,xi)*part3)
        end do
    end if

    ! 9: 3RD DERIVATIVE OF C_1 WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires z3_PCSAFT(1), z3_PCSAFT(4), z3_PCSAFT(5), z3_PCSAFT(9)
    if (GETDERC(9) .eq. 1) then
        part1 = 2.d0/recurring_factor**5
        part2 = (8.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(9) - 72.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) + 96.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**3 + 4320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 10132.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(9) - 6240.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 - 100740.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 138100.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(9) + &
     & 145440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 + 1355160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1281042.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(9) &
     & - 1911840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 - 12264762.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 8699406.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(9) + 16776072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 + 80838144.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 45143719.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **21*gl%z3_PCSAFT(9) - 107075112.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 - &
     & 406523883.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 183816499.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) + 522965640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4)**3 + 1606031790.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 596855765.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) - 2018666280.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 - 5078032095.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) - 1559038965.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(9) + &
     & 6288866640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 + 12995106540.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 3284598190.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(9) - 16019373360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**3 - &
     & 27059344050.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 5562096494.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(9) + 33585528240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4)**3 + 45836776188.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
     & + 7488991932.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) - 58000880640.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**3 - 62758564440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 7843799980.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(9) + &
     & 82080870360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 + 68426302320.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 6121051607.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(9) - 94029307800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**3 - &
     & 57716155947.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 3223267907.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) + 85289793912.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4)**3 + 35517155238.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
     & + 787823073.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) - 58900185384.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 - 13722112611.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 267498927.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) + &
     & 28644768480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**3 + 1308852180.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 286936800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(9) - 7907286960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 + &
     & 1738642860.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 65676312.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) - 126189360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) &
     & **3 - 847814472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 18386676.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(9) + 892093392.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4)**3 + 34093116.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 9577224.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) - 187625808.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(4)**3 + 64142928.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
     & - 453600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) - 31687200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **4*gl%z3_PCSAFT(4)**3 - 5339520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 466560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(9) + 5339520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 &
     & *gl%z3_PCSAFT(4)**3 - 2799360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 41472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + 2799360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4)**3 - 248832.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 248832.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 + 12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27* &
     & gl%z3_PCSAFT(9) - 108.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 180.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) + 144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) &
     & **3 + 720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1128.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(9) + 240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 + &
     & 36000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 57720.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(9) - 79080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 - &
     & 860400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 797313.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(9) + 1503000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 + &
     & 9880533.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 6707907.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(9) - 15472608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 &
     & - 74454048.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 40047441.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(9) + 107658144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4)**3 + 408564837.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 180818301.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) - 553996680.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 - 1719334710.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) - 638285700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) + 2214098520.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 + 5716766610.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 &
     & *gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1794164886.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(9) - &
     & 7081914000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 - 15292948896.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 4052814534.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(9) + 18478776096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**3 + &
     & 33262492086.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 7372263618.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(9) - 39798829608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4)**3 - 59103082116.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
     & - 10747951845.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) + 71186731320.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**3 + 85754469285.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 12392896155.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(9) - &
     & 105890970000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 - 100975745160.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 11001878961.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(9) + 130609564320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**3 + &
     & 95179059981.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 7111849077.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) - 132613716216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 &
     & *gl%z3_PCSAFT(4)**3 - 69983819118.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
     & - 2890522962.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) + 109481915064.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 + 38153530944.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 283371936.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) - &
     & 72060512400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**3 - 13634879400.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 441024480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(9) + 36569239920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 + &
     & 1766472840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 255569388.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) - 13361850480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(4)**3 + 1039335768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 31862544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(9) + 2865174432.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 - 617031936.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 19299840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) + 53443008.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3 + 102267360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 6290064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) - &
     & 251110080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 + 16180560.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 433440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(9) + 59037120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 - 7274880.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 262656.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 &
     & *gl%z3_PCSAFT(9) + 3715200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 + 290304.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 13824.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(9) - 2156544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 + 82944.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 248832.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**3 - 18.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(9) + 162.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 510.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) - 216.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 25*gl%z3_PCSAFT(4)**3 - 3960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 6507.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(9) + 4440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 + &
     & 39555.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 46995.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 24*gl%z3_PCSAFT(9) - 29460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 - 155970.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 180729.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(9) - 84420.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 - 681741.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 17727.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(9) &
     & + 3098136.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 + 13817088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 5063862.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(9) - &
     & 29140944.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 - 103536594.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 37216986.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) + &
     & 172725900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 + 513511020.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 19*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 167991780.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) - &
     & 750534060.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 - 1899802980.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 554465796.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(9) + &
     & 2545545000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 + 5507018496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1415905185.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(9) &
     & - 6972789336.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**3 - 12828847545.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 2863796673.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(9) + 15747258180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**3 + &
     & 24353987886.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 4625576175.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) - 29691231180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)**3 - 37935586125.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 5952877095.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(9) + 47059454280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **13*gl%z3_PCSAFT(4)**3 + 48553307640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 6024586500.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(9) - 62838881280.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**3 - 50855377320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 4653593928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) + &
     & 70551739140.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**3 + 43120518552.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 2564330712.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(9) - 66173335476.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 - &
     & 28957224564.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 819289380.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) + 51233990160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**3 &
     & + 14725859760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 31531440.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) - 32100083760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 - &
     & 5069442240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 176881392.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) + 15748180800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 + &
     & 696395808.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 78927648.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(9) - 5703883728.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 + &
     & 354681792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 7848432.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) + 1340129664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3 - &
     & 242702784.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 5432448.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) - 120054720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 + &
     & 55008000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1833600.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(9) - 31609920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 + &
     & 2148480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 17088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(9) + 10270080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 - 3284928.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 61056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(9) - 709632.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 + 297216.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 4608.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) - 27648.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**3 + 55296.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 5.d0*gl%z3_PCSAFT(1)**27 &
     & *gl%z3_PCSAFT(9) - 45.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 155.d0*gl%z3_PCSAFT(1)**26* &
     & gl%z3_PCSAFT(9) + 60.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**3 + 1260.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) &
     & *gl%z3_PCSAFT(5) + 2269.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(9) - 1500.d0*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4)**3 - 16185.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 20665.d0*gl%z3_PCSAFT(1) &
     & **24*gl%z3_PCSAFT(9) + 16320.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 + 123690.d0*gl%z3_PCSAFT(1)** &
     & 23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 129576.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(9) - 94560.d0* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 - 596046.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 581070.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(9) + 222156.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 + &
     & 1596960.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1842142.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(9) + 1038300.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 + 483846.d0*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 3663394.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) - 13210320.d0* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 - 28008300.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 1237565.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) + 73893720.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 + &
     & 157354215.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 22351321.d0*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(9) - 287438940.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 - 562674396.d0*gl%z3_PCSAFT(1) &
     & **17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 101052831.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(9) + &
     & 865340796.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**3 + 1507300899.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
     & *gl%z3_PCSAFT(5) + 272383179.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(9) - 2110814112.d0*gl%z3_PCSAFT(1) &
     & **15*gl%z3_PCSAFT(4)**3 - 3196832118.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 536685698.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) + 4273118160.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) &
     & **3 + 5503832460.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 816364400.d0*gl%z3_PCSAFT(1) &
     & **14*gl%z3_PCSAFT(9) - 7271778300.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 - 7781190720.d0* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 973370092.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(9) + &
     & 10462467300.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**3 + 9062135232.d0*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 904431456.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) - &
     & 12734595312.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**3 - 8658408984.d0*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 636388376.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) + &
     & 13066065432.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 + 6701489832.d0*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 313940656.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) - &
     & 11217968880.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**3 - 4093386720.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
     & *gl%z3_PCSAFT(5) - 81796320.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) + 7967593440.d0*gl%z3_PCSAFT(1)**8 &
     & *gl%z3_PCSAFT(4)**3 + 1868162160.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 15543328.d0 &
     & *gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) - 4604185440.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 - &
     & 550329792.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 25154336.d0*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(9) + 2112183072.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 + 37794816.d0*gl%z3_PCSAFT(1)** &
     & 6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 10226688.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) - 740144448.d0* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3 + 53454336.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 1044928.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) + 184923840.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 - &
     & 28025280.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 704320.d0*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(9) - 28049280.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 + 6084480.d0*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 285952.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + 1059840.d0*gl%z3_PCSAFT(1) &
     & **2*gl%z3_PCSAFT(4)**3 - 200448.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 21504.d0* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 415488.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 - 170496.d0*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 8448.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) - 56832.d0*gl%z3_PCSAFT(4)**3 + &
     & 25344.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1536.d0*gl%z3_PCSAFT(9))
        part3 = (-8.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**27 + 440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26 - 10132.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**25 + 138100.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 - 1281042.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **23 + 8699406.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 - 45143719.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 &
     & + 183816499.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 - 596855765.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 + &
     & 1559038965.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 - 3284598190.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 + &
     & 5562096494.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 - 7488991932.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 + &
     & 7843799980.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 - 6121051607.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 + &
     & 3223267907.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 - 787823073.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 - &
     & 267498927.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 + 286936800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9 - &
     & 65676312.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 - 18386676.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 + &
     & 9577224.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 + 453600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 - 466560.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - 41472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 + 4.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **27 - 460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + 13244.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 - &
     & 202260.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 2016349.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - &
     & 14431331.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 78155038.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - &
     & 330554948.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 1113258445.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - &
     & 3019253933.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 6626827552.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - &
     & 11761804784.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 16761220155.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & - 18876016645.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 16286484278.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 13 - 10082824036.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + 3711709451.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **11 - 26762023.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 798318560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 9 + 381960684.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - 24373332.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - &
     & 37082520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 9145872.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + &
     & 1317600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 - 442368.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 41472.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 + 6.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27 + 30.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 26 - 4641.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 + 95385.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - &
     & 1090032.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 8474484.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - &
     & 48688329.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 215978427.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - &
     & 758924190.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 2144593512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - &
     & 4912183935.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 9138225591.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - &
     & 13754686200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 16567258560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 &
     & - 15615882255.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 11010130041.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 12 - 5200854984.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 1004149350.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **10 + 591264840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 535070244.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 8 + 145643496.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 18873696.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - &
     & 20681136.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 + 2743200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 811584.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 - 162432.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 13824.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1) - 5.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 + 95.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 - 226.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**25 - 12740.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 + 218271.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - &
     & 2000115.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 12636152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 59904284.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 222011305.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 657907219.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**18 + 1578384894.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 3082386636.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**16 + 4897630037.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 6288286805.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14 + 6422334508.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 5054931204.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12 + 2858478584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 931324804.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **10 - 66912720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 254704672.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - &
     & 124570304.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 16068432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 9788288.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 4450240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + 271168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 3 + 192384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 29952.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 1536.d0*gl%mmean_PCSAFT + &
     & gl%z3_PCSAFT(1)**27 - 25.d0*gl%z3_PCSAFT(1)**26 + 251.d0*gl%z3_PCSAFT(1)**25 - 825.d0*gl%z3_PCSAFT(1)**24 - &
     & 8606.d0*gl%z3_PCSAFT(1)**23 + 141640.d0*gl%z3_PCSAFT(1)**22 - 1093886.d0*gl%z3_PCSAFT(1)**21 + &
     & 5788246.d0*gl%z3_PCSAFT(1)**20 - 23126915.d0*gl%z3_PCSAFT(1)**19 + 72760671.d0*gl%z3_PCSAFT(1)**18 - &
     & 184152457.d0*gl%z3_PCSAFT(1)**17 + 378996827.d0*gl%z3_PCSAFT(1)**16 - 636539640.d0*gl%z3_PCSAFT(1)** &
     & 15 + 870030730.d0*gl%z3_PCSAFT(1)**14 - 957815972.d0*gl%z3_PCSAFT(1)**13 + 830741176.d0* &
     & gl%z3_PCSAFT(1)**12 - 541743704.d0*gl%z3_PCSAFT(1)**11 + 235258168.d0*gl%z3_PCSAFT(1)**10 - &
     & 35468160.d0*gl%z3_PCSAFT(1)**9 - 33945760.d0*gl%z3_PCSAFT(1)**8 + 28498784.d0*gl%z3_PCSAFT(1)**7 - &
     & 8924064.d0*gl%z3_PCSAFT(1)**6 - 232640.d0*gl%z3_PCSAFT(1)**5 + 1155520.d0*gl%z3_PCSAFT(1)**4 - 334208.d0 &
     & *gl%z3_PCSAFT(1)**3 + 256.d0*gl%z3_PCSAFT(1)**2 + 17152.d0*gl%z3_PCSAFT(1) - 2560.d0)
        part4 = (72.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 4320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(4) + 100740.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 1355160.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 12264762.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - &
     & 80838144.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 406523883.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **20*gl%z3_PCSAFT(4) - 1606031790.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + &
     & 5078032095.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 12995106540.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 27059344050.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - &
     & 45836776188.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 62758564440.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 68426302320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 57716155947.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 35517155238.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 13722112611.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - &
     & 1308852180.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 1738642860.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 847814472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - &
     & 34093116.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 64142928.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 5*gl%z3_PCSAFT(4) + 5339520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 2799360.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 248832.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 36.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 5040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) &
     & - 144060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 2131800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 23*gl%z3_PCSAFT(4) - 20434389.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 140259984.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 728584806.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) + 2962385280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 9631957065.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 25380758268.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) - 54602032788.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 96107834328.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 137960300715.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4) + 159925025400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 146782978998.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 102694402704.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 50421841767.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + &
     & 13299684540.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 1680657660.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 2891882304.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 870980148.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 35171280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **5*gl%z3_PCSAFT(4) - 64776240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 3749760.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 1721088.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
     & + 165888.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 54.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26* &
     & gl%z3_PCSAFT(4) - 1080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 64845.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 1160910.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + &
     & 12196242.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 88602696.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **21*gl%z3_PCSAFT(4) + 479680893.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - &
     & 2018188350.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 6771108420.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 18410944392.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + &
     & 40963039875.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 74907180522.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 112527942030.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - &
     & 138043001160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 136455337755.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 105955331634.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) &
     & + 61424139258.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 23427819360.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 3065392440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + &
     & 2383207704.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 1568768184.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 325781568.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 37733040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 27129600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(4) + 2137536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 338688.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) + 45.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 &
     & *gl%z3_PCSAFT(4) - 540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 7170.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 &
     & *gl%z3_PCSAFT(4) + 239640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 3002811.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 23826240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - &
     & 136497264.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 599357640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 2085816705.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 5871433764.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 13536222306.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
     & + 25732502952.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 40416911505.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 52300226040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 55271817708.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 46877710776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) - 30853160688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + &
     & 14649813600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 3990385200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 &
     & *gl%z3_PCSAFT(4) - 296321952.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 830491776.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 351760704.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 48605760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 12998400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) - 5928768.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 601344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & *gl%z3_PCSAFT(4) + 25344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) - 9.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 180.d0* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 915.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 11370.d0*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(4) + 240096.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 2242560.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) + 13993782.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 65053260.d0*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) + 236752815.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 693408996.d0*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) + 1661531913.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 3288632154.d0*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4) + 5400482010.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 7356376320.d0*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(4) + 8267982612.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 7574682144.d0*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) + 5525610048.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 3061260960.d0* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 1143128160.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 155911680.d0* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 116871936.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 88706688.d0* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 26606400.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 1530240.d0*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4) + 1546368.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 497664.d0*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) + 51456.d0*gl%z3_PCSAFT(4))
        part5 = (72.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26* &
     & gl%z3_PCSAFT(5) - 288.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - 4320.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 18720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + &
     & 100740.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 436320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 &
     & *gl%z3_PCSAFT(4)**2 - 1355160.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 5735520.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 12264762.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(5) - 50328216.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 80838144.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 321225336.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4)**2 + 406523883.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - &
     & 1568896920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 1606031790.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 6055998840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 &
     & + 5078032095.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - 18866599920.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 12995106540.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(5) + 48058120080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + &
     & 27059344050.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 100756584720.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 45836776188.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(5) + 174002641920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + &
     & 62758564440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 246242611080.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 68426302320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(5) + 282087923400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + &
     & 57716155947.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 255869381736.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 35517155238.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(5) + 176700556152.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + &
     & 13722112611.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 85934305440.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 1308852180.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + &
     & 23721860880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 1738642860.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 378568080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + &
     & 847814472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 2676280176.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 34093116.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + &
     & 562877424.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 64142928.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 95061600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + &
     & 5339520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 16018560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 3*gl%z3_PCSAFT(4)**2 + 2799360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 8398080.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 248832.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(5) - 746496.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 36.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 144.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + 5040.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 23760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4)**2 - 144060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 664560.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 2131800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(5) - 9457200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 - 20434389.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) + 86825412.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4)**2 + 140259984.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - 571367556.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 728584806.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(5) + 2857652820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + &
     & 2962385280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 11264263380.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 9631957065.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) &
     & + 35828078760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 25380758268.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 93349189944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4)**2 - 54602032788.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + &
     & 200959722792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + 96107834328.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 358583760360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)**2 - 137960300715.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + &
     & 529277669460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + 159925025400.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 641406655620.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)**2 - 146782978998.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + &
     & 629157900564.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 102694402704.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 487616559156.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)**2 - 50421841767.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + &
     & 286257214080.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 13299684540.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 116803976400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) &
     & **2 + 1680657660.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 25517661120.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 2891882304.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(5) + 2084428512.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 870980148.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 3098550672.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4)**2 + 35171280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 684145440.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 64776240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(5) + 362880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 3749760.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 5909760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 &
     & + 1721088.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 2363904.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)**2 + 165888.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 248832.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(4)**2 - 54.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 216.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - 1080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + &
     & 7560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + 64845.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(5) - 342360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 1160910.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 5582520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4)**2 + 12196242.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 54744876.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 - 88602696.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(5) + 375185844.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 479680893.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 1932776100.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4)**2 - 2018188350.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + &
     & 7809677100.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 6771108420.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - 25428945600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)** &
     & 2 - 18410944392.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 67895847216.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + 40963039875.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(5) - 150250155120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - &
     & 74907180522.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 277040427120.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + 112527942030.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(5) - 425864576700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - &
     & 138043001160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 543551875380.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + 136455337755.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(5) - 570793387140.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - &
     & 105955331634.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + 485439315516.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 61424139258.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(5) - 325700768880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - &
     & 23427819360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 164491741440.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 3065392440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - &
     & 56402103600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + 2383207704.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 8873953488.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - &
     & 1568768184.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 2263964256.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + 325781568.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - &
     & 1668902400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 37733040.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 374423040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - &
     & 27129600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 22913280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **2*gl%z3_PCSAFT(4)**2 + 2137536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 2343168.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 338688.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - &
     & 41472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**2 + 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(5) + 45.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - 540.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - &
     & 7170.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) + 62100.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) &
     & **2 + 239640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 1357380.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 22*gl%z3_PCSAFT(4)**2 - 3002811.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) + 14790168.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 23826240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(5) - 107278920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 136497264.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 574057440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) &
     & **2 + 599357640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 2390025600.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 2085816705.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + &
     & 7994097900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 5871433764.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **17*gl%z3_PCSAFT(5) - 21924636012.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - &
     & 13536222306.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 49944734964.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **15*gl%z3_PCSAFT(4)**2 + 25732502952.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - &
     & 95193677700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 40416911505.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 152212960800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 &
     & + 52300226040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 203876519760.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 55271817708.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) + &
     & 227432736504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 46877710776.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 209036995704.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 &
     & - 30853160688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 155488407120.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 14649813600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - &
     & 90806780160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 3990385200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **8*gl%z3_PCSAFT(5) + 39296679840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - &
     & 296321952.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 10873376064.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 &
     & *gl%z3_PCSAFT(4)**2 + 830491776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 721705536.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 351760704.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) &
     & + 863648640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 48605760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 4*gl%z3_PCSAFT(5) - 409789440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 12998400.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 82944000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 &
     & - 5928768.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 6218496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)**2 + 601344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 170496.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(4)**2 + 25344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(5) - 9.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 36.d0 &
     & *gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + 180.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 540.d0*gl%z3_PCSAFT(1)** &
     & 24*gl%z3_PCSAFT(4)**2 - 915.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 1260.d0*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(4)**2 - 11370.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 106740.d0*gl%z3_PCSAFT(1)**22* &
     & gl%z3_PCSAFT(4)**2 + 240096.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 1416528.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4)**2 - 2242560.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 11187360.d0*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4)**2 + 13993782.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 62985240.d0*gl%z3_PCSAFT(1)**19 &
     & *gl%z3_PCSAFT(4)**2 - 65053260.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 272126520.d0*gl%z3_PCSAFT(1)** &
     & 18*gl%z3_PCSAFT(4)**2 + 236752815.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - 939159900.d0* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 693408996.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + &
     & 2653372548.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + 1661531913.d0*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(5) - 6232510332.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 3288632154.d0*gl%z3_PCSAFT(1) &
     & **15*gl%z3_PCSAFT(5) + 12284568900.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + 5400482010.d0* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - 20410098840.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - &
     & 7356376320.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 28596660120.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) &
     & **2 + 8267982612.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 33674651136.d0*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4)**2 - 7574682144.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + 33077668176.d0* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 5525610048.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - &
     & 26755189920.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 3061260960.d0*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(5) + 17445036960.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 1143128160.d0*gl%z3_PCSAFT(1) &
     & **8*gl%z3_PCSAFT(5) - 8830411200.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 155911680.d0* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 3203850240.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - &
     & 116871936.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 641016576.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + &
     & 88706688.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 68371200.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - &
     & 26606400.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 101157120.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + &
     & 1530240.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 35861760.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + &
     & 1546368.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 6437376.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 497664.d0 &
     & *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 497664.d0*gl%z3_PCSAFT(4)**2 + 51456.d0*gl%z3_PCSAFT(5))
        part6 = (24.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) - 288.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(9) + 480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 + 18720.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 33580.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(9) - 33600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 - 436320.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 451720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(9) + 774720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 + 5735520.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 4088254.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **22*gl%z3_PCSAFT(9) - 9880320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 - &
     & 50328216.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 26946048.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(9) + 83777208.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) &
     & **3 + 321225336.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 135507961.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) - 518458800.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 - 1568896920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) - 535343930.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) + 2474565000.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 + 6055998840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 &
     & *gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1692677365.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(9) - &
     & 9431202720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 - 18866599920.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 4331702180.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(9) + 29333020560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**3 + &
     & 48058120080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 9019781350.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(9) - 75384706080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15* &
     & gl%z3_PCSAFT(4)**3 - 100756584720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
     & - 15278925396.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) + 160954936080.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**3 + 174002641920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 20919521480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(9) - &
     & 285445862880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 - 246242611080.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 22808767440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(9) + 418296100680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**3 + &
     & 282087923400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 19238718649.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) - 501298508400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 11*gl%z3_PCSAFT(4)**3 - 255869381736.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) - 11839051746.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) + &
     & 483148526808.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 + 176700556152.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 4574037537.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 10*gl%z3_PCSAFT(9) - 364731247200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**3 - &
     & 85934305440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 436284060.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) + 206512323840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4)**3 + 23721860880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 579547620.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) - 81002473920.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 + 378568080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 282604824.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(9) + 18305619120.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 - 2676280176.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 11364372.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) - &
     & 962213472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3 + 562877424.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 21380976.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(9) - 231843600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 + 95061600.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1779840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 4*gl%z3_PCSAFT(9) - 126748800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 - &
     & 16018560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 933120.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + 21358080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 - &
     & 8398080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 82944.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 11197440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 - &
     & 746496.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 995328.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(4)**3 - 12.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) + 144.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(9) - 240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 - 23760.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 48020.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(9) + 45600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 + 664560.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 710600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23* &
     & gl%z3_PCSAFT(9) - 1236240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 - 9457200.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 6811463.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **22*gl%z3_PCSAFT(9) + 16828800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 + &
     & 86825412.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 46753328.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(9) - 147694836.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4)**3 - 571367556.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 242861602.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) + 934006440.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 + 2857652820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 987461760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) - 4530327180.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 - 11264263380.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 3210652355.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(9) &
     & + 17516769840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 + 35828078760.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 8460252756.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 17*gl%z3_PCSAFT(9) - 55316597400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**3 - &
     & 93349189944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 18200677596.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(9) + 144756177072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
     & *gl%z3_PCSAFT(4)**3 + 200959722792.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 32035944776.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) - &
     & 316247600760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**3 - 358583760360.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 45986766905.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 14*gl%z3_PCSAFT(9) + 578004045120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 + &
     & 529277669460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 53308341800.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(9) - 881858164740.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 12*gl%z3_PCSAFT(4)**3 - 641406655620.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) - 48927659666.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) + &
     & 1116432101640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**3 + 629157900564.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 34231467568.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) - 1161220259532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)**3 - 487616559156.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
     & - 16807280589.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) + 977736633360.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**3 + 286257214080.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 4433228180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) - &
     & 652195229040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 - 116803976400.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 560219220.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(9) + 333701936640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 + &
     & 25517661120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 963960768.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(9) - 124396843680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6* &
     & gl%z3_PCSAFT(4)**3 + 2084428512.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 290326716.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) + 30776517504.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3 - 3098550672.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 11723760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) - 4032573840.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 + 684145440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 21592080.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(9) + &
     & 29520000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 + 362880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1249920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + &
     & 41765760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 - 5909760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 573696.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) + &
     & 760320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 - 2363904.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 55296.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) - 580608.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4)**3 - 248832.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 18.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) + 216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) &
     & *gl%z3_PCSAFT(5) - 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(9) - 360.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 + 7560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 21615.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(9) - 18000.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 - 342360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) - 386970.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(9) + 693360.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 + 5582520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) &
     & *gl%z3_PCSAFT(5) + 4065414.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(9) - 10451520.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 - 54744876.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 29534232.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(9) + &
     & 96243948.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 + 375185844.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 159893631.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(9) - 626296680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 - &
     & 1932776100.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 672729450.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) + 3100746420.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4)**3 + 7809677100.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 2257036140.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(9) - 12202475760.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 - 25428945600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) &
     & *gl%z3_PCSAFT(5) - 6136981464.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(9) + &
     & 39226842720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**3 + 67895847216.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 13654346625.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 &
     & *gl%z3_PCSAFT(9) - 104745004128.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**3 - &
     & 150250155120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 24969060174.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) + 234507806280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 14*gl%z3_PCSAFT(4)**3 + 277040427120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 37509314010.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(9) - &
     & 441980411040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 - 425864576700.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 46014333720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 13*gl%z3_PCSAFT(9) + 701277133500.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**3 + &
     & 543551875380.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 45485112585.d0 &
     & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) - 933775773480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 11*gl%z3_PCSAFT(4)**3 - 570793387140.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) - 35318443878.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) + &
     & 1037055940740.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 + 485439315516.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 20474713086.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) - 951731562960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
     & **3 - 325700768880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 7809273120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) + 712154731680.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 + 164491741440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 1021797480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) - &
     & 426292191360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 - 56402103600.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 794402568.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(9) + 198508224240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 + &
     & 8873953488.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 522922728.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) - 68845453344.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4)**3 + 2263964256.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 108593856.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) + 16485539040.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 - 1668902400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 12577680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(9) - 2314529280.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 + 374423040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 9043200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + &
     & 98904960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 - 22913280.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 712512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(9) + 11289600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 - 2343168.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 112896.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) &
     & + 71424.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**3 - 41472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
     & + 9216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(9) + 15.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) - 180.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(9) + 300.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 + 900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 2390.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(9) + 600.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 + 62100.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 79880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(9) - 153540.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 - 1357380.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) - 1000937.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(9) + 2770080.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 + 14790168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 7942080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(9) - 27361344.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 - 107278920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) - 45499088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) + 185019360.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 + 574057440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 199785880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) - &
     & 940262520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 - 2390025600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 695272235.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(9) + &
     & 3780353760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 + 7994097900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1957144588.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(9) &
     & - 12405215460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**3 - 21924636012.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 4512074102.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(9) + 33868011576.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**3 + &
     & 49944734964.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 8577500984.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) - 77795442900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4)**3 - 95193677700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 13472303835.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(9) + 151201608960.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 + 152212960800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 17433408680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(9) - &
     & 249075491640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**3 - 203876519760.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 18423939236.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(9) + 347300756880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**3 + &
     & 227432736504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 15625903592.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) - 408319680432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)**3 - 209036995704.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 10284386896.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) + 402194923200.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**3 + 155488407120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 4883271200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) - 328830663120.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 - 90806780160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1330128400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) + &
     & 220258643520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 + 39296679840.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 98773984.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7* &
     & gl%z3_PCSAFT(9) - 118671262080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 - &
     & 10873376064.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 276830592.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) + 50087747712.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3 + &
     & 721705536.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 117253568.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) - 15912730560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 + &
     & 863648640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 16201920.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(9) + 3563527680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 - &
     & 409789440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 4332800.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) - 496128000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 + &
     & 82944000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1976256.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 30612480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 - &
     & 6218496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 200448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(9) + 415488.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**3 - 170496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) + 8448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(9) - 3.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) + 36.d0* &
     & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 60.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(9) - 60.d0* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 - 540.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 305.d0* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(9) + 600.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 - 1260.d0*gl%z3_PCSAFT(1) &
     & **23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 3790.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(9) + 9540.d0* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 + 106740.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 80032.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(9) - 258240.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 - &
     & 1416528.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 747520.d0*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(9) + 2845464.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 + 11187360.d0*gl%z3_PCSAFT(1)**20 &
     & *gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 4664594.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) - 20282640.d0* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 - 62985240.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 21684420.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) + 106575240.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 &
     & + 272126520.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 78917605.d0*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(9) - 439704000.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 - 939159900.d0*gl%z3_PCSAFT(1) &
     & **17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 231136332.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(9) + &
     & 1477157940.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**3 + 2653372548.d0*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 553843971.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(9) - 4130732904.d0 &
     & *gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**3 - 6232510332.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)* &
     & gl%z3_PCSAFT(5) - 1096210718.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) + 9742040340.d0*gl%z3_PCSAFT(1) &
     & **14*gl%z3_PCSAFT(4)**3 + 12284568900.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 1800160670.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(9) - 19515883680.d0*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4)**3 - 20410098840.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 2452125440.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(9) + 33307509840.d0*gl%z3_PCSAFT(1)**12* &
     & gl%z3_PCSAFT(4)**3 + 28596660120.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 2755994204.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) - 48430636320.d0*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4)**3 - 33674651136.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
     & 2524894048.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) + 59856655008.d0*gl%z3_PCSAFT(1)**10* &
     & gl%z3_PCSAFT(4)**3 + 33077668176.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 1841870016.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) - 62601407520.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
     & **3 - 26755189920.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1020420320.d0* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) + 55032395040.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 + &
     & 17445036960.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 381042720.d0*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(9) - 40282467840.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 - 8830411200.d0* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 51970560.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(9) + &
     & 24233623680.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 + 3203850240.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) &
     & *gl%z3_PCSAFT(5) - 38957312.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) - 11766783360.d0*gl%z3_PCSAFT(1)** &
     & 5*gl%z3_PCSAFT(4)**3 - 641016576.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 29568896.d0 &
     & *gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) + 4493738880.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 - &
     & 68371200.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 8868800.d0*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(9) - 1298611200.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 + 101157120.d0*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 510080.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + 266760960.d0* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 - 35861760.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
     & 515456.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) - 34690560.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 + &
     & 6437376.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 165888.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) + &
     & 2145792.d0*gl%z3_PCSAFT(4)**3 - 497664.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 17152.d0* &
     & gl%z3_PCSAFT(9))
        do xi = 1, gl%ncomp 
            gl%cx1_PCSAFT(9,xi) = part1*(gl%mPCSAFT(xi)*part2 + gl%z3x1_PCSAFT(9,xi)*part3 + gl%z3x1_PCSAFT(5,xi)* &
			& part4 + gl%z3x1_PCSAFT(4,xi)*part5 + gl%z3x1_PCSAFT(1,xi)*part6)
        end do
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF C_1 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires z3_PCSAFT(1), z3_PCSAFT(4)
    if (GETDERC(10) .eq. 1) then
        part1 = 4.d0/recurring_factor**5
        part2 = (24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(4) - 1680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 &
     & *gl%z3_PCSAFT(4) + 39140.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 504200.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 4299782.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - &
     & 26591508.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 125974859.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
     & **21*gl%z3_PCSAFT(4) - 473989210.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + &
     & 1451755955.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 3677984500.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 7772982770.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - &
     & 13721514924.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 20120913700.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 24205886460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + &
     & 23406178307.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - 17611040946.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 9748346703.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - &
     & 3517359420.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + 516452940.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 163441872.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - &
     & 82448532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 5537376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 &
     & *gl%z3_PCSAFT(4) + 889920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 466560.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 41472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 36.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(4) + 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - &
     & 27540.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 464700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24* &
     & gl%z3_PCSAFT(4) - 4442793.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 29011056.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 140810061.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) &
     & + 533937690.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 1635368130.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 4141738416.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - &
     & 8811917442.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 15892338288.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 24360661905.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + &
     & 31646200440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 34580504781.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 31413017826.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - &
     & 23312412552.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 13739660160.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 6092100960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 1779032472.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 178955808.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 91465920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + &
     & 34912080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 567360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4* &
     & gl%z3_PCSAFT(4) - 981504.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 96768.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 54.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(4) + 900.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 1545.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 94200.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 1321821.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - &
     & 9964776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 51850752.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
     & gl%z3_PCSAFT(4) - 204096690.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 639504840.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 1650721836.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + &
     & 3597346575.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 6727619628.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 16*gl%z3_PCSAFT(4) + 10884531765.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 15235004760.d0 &
     & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 18324077130.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13* &
     & gl%z3_PCSAFT(4) - 18713161554.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 15964586892.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 11141421960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) &
     & + 6184276320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 2619809928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 8*gl%z3_PCSAFT(4) + 788292096.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 140928288.d0*gl%mmean_PCSAFT &
     & *gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 2531040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 5851200.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 1449792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + &
     & 85248.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 18432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 15.d0 &
     & *gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(4) - 330.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 2765.d0*gl%z3_PCSAFT(1)**25* &
     & gl%z3_PCSAFT(4) - 6050.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 87604.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + &
     & 1051470.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 6443878.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + &
     & 27610760.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - 91268065.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + &
     & 245112266.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 552973423.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + &
     & 1070948374.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 1801278330.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + &
     & 2637503410.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - 3346585912.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + &
     & 3646896388.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 3375154496.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + &
     & 2619334480.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 1679372000.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + &
     & 872648272.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 357473952.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + &
     & 110280032.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 23366400.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 2558080.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 140928.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 85248.d0* &
     & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 8448.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4))
        part3 = (-24.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**27 + 1680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26 - 39140.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**25 + 504200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 - 4299782.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **23 + 26591508.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 - 125974859.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 21 + 473989210.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 - 1451755955.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 19 + 3677984500.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 - 7772982770.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
     & **17 + 13721514924.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 - 20120913700.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**15 + 24205886460.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 - 23406178307.d0*gl%mmean_PCSAFT &
     & **4*gl%z3_PCSAFT(1)**13 + 17611040946.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 - 9748346703.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 + 3517359420.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 - 516452940.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9 - 163441872.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 + 82448532.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 - 5537376.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 - 889920.d0*gl%mmean_PCSAFT**4 &
     & *gl%z3_PCSAFT(1)**5 - 466560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - 41472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 &
     & + 12.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27 - 2280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26 + 62740.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**25 - 865600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 7659439.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**23 - 48474598.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 233413868.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**21 - 889915470.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 2760694105.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**19 - 7097945568.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 15292609536.d0*gl%mmean_PCSAFT** &
     & 3*gl%z3_PCSAFT(1)**17 - 27728015284.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 42226178005.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 53592767470.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + &
     & 55931990428.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 47037958958.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 &
     & + 30902255091.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 - 15034101220.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
     & 10 + 4813162740.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 - 616556016.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 &
     & - 226098396.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 125748000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - &
     & 21531600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 264960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 179712.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 13824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 + 18.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**27 + 900.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26 - 35445.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25 &
     & + 543450.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 5058732.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + &
     & 32996742.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 162235719.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + &
     & 628883400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 1981121460.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + &
     & 5178993072.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 11387345895.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 &
     & + 21204344346.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 33468115440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 15 + 44577645510.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 49647118605.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**13 + 45588108708.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 33808748394.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**11 + 19606017120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 8378553120.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 + 2273394816.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 145595352.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 169556544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 74981520.d0* &
     & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 12862080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 321984.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)**3 + 105984.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 9216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - &
     & 15.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27 - 30.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26 + 7960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 &
     & - 146350.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 + 1464091.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 9937740.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 50177152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 198551720.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**20 + 637077415.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 1696961414.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **18 + 3812048392.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 7288111966.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 &
     & + 11896522965.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 16546011280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
     & 19481516848.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 19213595692.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
     & 15630347624.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 10251192160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
     & 5219318240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 1911003328.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + &
     & 397114848.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 26687872.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 52096320.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 18156800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 3012672.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **3 + 172032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 8448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 3.d0*gl%z3_PCSAFT(1)**27 - &
     & 30.d0*gl%z3_PCSAFT(1)**26 - 515.d0*gl%z3_PCSAFT(1)**25 + 14000.d0*gl%z3_PCSAFT(1)**24 - 156056.d0*gl%z3_PCSAFT(1) &
     & **23 + 1117040.d0*gl%z3_PCSAFT(1)**22 - 5832946.d0*gl%z3_PCSAFT(1)**21 + 23670000.d0*gl%z3_PCSAFT(1) &
     & **20 - 77609045.d0*gl%z3_PCSAFT(1)**19 + 211092426.d0*gl%z3_PCSAFT(1)**18 - 484907751.d0* &
     & gl%z3_PCSAFT(1)**17 + 951217432.d0*gl%z3_PCSAFT(1)**16 - 1601522470.d0*gl%z3_PCSAFT(1)**15 + &
     & 2313984580.d0*gl%z3_PCSAFT(1)**14 - 2856447652.d0*gl%z3_PCSAFT(1)**13 + 2988050648.d0*gl%z3_PCSAFT(1) &
     & **12 - 2617328304.d0*gl%z3_PCSAFT(1)**11 + 1887085840.d0*gl%z3_PCSAFT(1)**10 - 1090692480.d0 &
     & *gl%z3_PCSAFT(1)**9 + 482004480.d0*gl%z3_PCSAFT(1)**8 - 145793408.d0*gl%z3_PCSAFT(1)**7 + 18173696.d0 &
     & *gl%z3_PCSAFT(1)**6 + 7990720.d0*gl%z3_PCSAFT(1)**5 - 5466880.d0*gl%z3_PCSAFT(1)**4 + 1588352.d0* &
     & gl%z3_PCSAFT(1)**3 - 248832.d0*gl%z3_PCSAFT(1)**2 + 17152.d0*gl%z3_PCSAFT(1))
        part4 = (72.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) - 5760.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) &
     & + 130060.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 1568200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
     & 23*gl%z3_PCSAFT(4) + 12424714.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 72025224.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 326859181.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) - 1213612730.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 3781454365.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 9985308500.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) + 22326192910.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 42000095556.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 65905831100.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4) - 85399406040.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 90233394229.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 76404304578.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11* &
     & gl%z3_PCSAFT(4) + 50540662497.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 25122947100.d0 &
     & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 8825640660.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8* &
     & gl%z3_PCSAFT(4) - 1982688696.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 247965444.d0* &
     & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 21380976.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5* &
     & gl%z3_PCSAFT(4) + 1779840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 933120.d0*gl%mmean_PCSAFT**4* &
     & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 82944.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 36.d0*gl%mmean_PCSAFT &
     & **3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 8640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - &
     & 223100.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) + 2820200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 &
     & *gl%z3_PCSAFT(4) - 22775273.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 132844844.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - 602923312.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) + 2236337760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) - 6983565215.d0* &
     & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 18605547996.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17* &
     & gl%z3_PCSAFT(4) - 42351329448.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) + 81982127096.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 134064069635.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14* &
     & gl%z3_PCSAFT(4) + 183919955540.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 210099189056.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 198022078144.d0*gl%mmean_PCSAFT**3 &
     & *gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) - 152066752389.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) &
     & + 93414878900.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 44626428540.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) + 15813916992.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - &
     & 3791660652.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 482580720.d0*gl%mmean_PCSAFT**3* &
     & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - 467280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 2309760.d0 &
     & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 1292544.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2* &
     & gl%z3_PCSAFT(4) - 110592.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 54.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
     & 26*gl%z3_PCSAFT(4) - 4320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) + 140055.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 1891050.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + &
     & 15690804.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) - 92558676.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **21*gl%z3_PCSAFT(4) + 421749441.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) - &
     & 1567515930.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 4917827100.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) - 13245585384.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + &
     & 30741479685.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) - 61265647614.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 104238162960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) - &
     & 150534303540.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 183484158195.d0*gl%mmean_PCSAFT**2 &
     & *gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) - 187558015014.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) &
     & + 159418233006.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 111294207840.d0*gl%mmean_PCSAFT &
     & **2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 62674507200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) &
     & - 27712355112.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 9229156296.d0*gl%mmean_PCSAFT**2* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 2161272384.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + &
     & 311645520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 18673920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
     & **3*gl%z3_PCSAFT(4) - 813888.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 85248.d0*gl%mmean_PCSAFT** &
     & 2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 9216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4) + 45.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26* &
     & gl%z3_PCSAFT(4) + 720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4) - 37760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 &
     & *gl%z3_PCSAFT(4) + 560000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) - 4821497.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4) + 28932480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) - &
     & 132925388.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) + 496612360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19* &
     & gl%z3_PCSAFT(4) - 1568481365.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) + 4274726368.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) - 10113305576.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
     & + 20715853664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) - 36534742455.d0*gl%mmean_PCSAFT* &
     & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4) + 55166107280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) - &
     & 70961955116.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4) + 77365368056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
     & **11*gl%z3_PCSAFT(4) - 71040780376.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) + &
     & 54474739520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) - 34467972880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
     & 8*gl%z3_PCSAFT(4) + 17696182496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) - 7198397664.d0* &
     & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) + 2240276032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) - &
     & 505055040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) + 74935040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(4) - 5914176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 86784.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)* &
     & gl%z3_PCSAFT(4) + 8448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4) - 9.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) + 3625.d0* &
     & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) - 61750.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4) + 558412.d0*gl%z3_PCSAFT(1)** &
     & 22*gl%z3_PCSAFT(4) - 3430600.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) + 15962054.d0*gl%z3_PCSAFT(1)**20* &
     & gl%z3_PCSAFT(4) - 60118740.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) + 191389975.d0*gl%z3_PCSAFT(1)**18* &
     & gl%z3_PCSAFT(4) - 527587752.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) + 1269857253.d0*gl%z3_PCSAFT(1)**16* &
     & gl%z3_PCSAFT(4) - 2664439958.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) + 4847183030.d0*gl%z3_PCSAFT(1)**14 &
     & *gl%z3_PCSAFT(4) - 7603003520.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4) + 10234554284.d0*gl%z3_PCSAFT(1)** &
     & 12*gl%z3_PCSAFT(4) - 11773819024.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4) + 11521274256.d0* &
     & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) - 9531629600.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) + 6610913760.d0* &
     & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) - 3799462080.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) + 1780567744.d0* &
     & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4) - 665317504.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) + 191949760.d0* &
     & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) - 40743040.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 5879936.d0*gl%z3_PCSAFT(1) &
     & **2*gl%z3_PCSAFT(4) - 497664.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 17152.d0*gl%z3_PCSAFT(4))
        do xi = 1, gl%ncomp 
            gl%cx1_PCSAFT(10,xi) = part1*(gl%mPCSAFT(xi)*part2 + gl%z3x1_PCSAFT(4,xi)*part3 + gl%z3x1_PCSAFT(1,xi)*part4)
        end do
    end if

    !DEC$ END IF
end subroutine CX1DERIVS



    end module pc_saft_CX1_derivs_module
