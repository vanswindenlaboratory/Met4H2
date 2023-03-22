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


    module sub_file_input_helper
    use module_all_types
    implicit none
    contains
    subroutine bwr_recalc(gl,nrsubst)


    implicit none

    type(type_gl) :: gl



    integer :: i, l, count, nrsubst, pol, pol_add,loop


    count = 19           !for MBWR exp-coeffs 20-40
    pol = 0              !number of polynomial terms
    pol_add = 0          !number of additional polynomial terms


    !only pol di
    if (gl%eos_coeff%nreg(nrsubst) == 32) then            !MBWR
        do l=1, 5
            gl%eos_coeff%di(l,nrsubst) = 1.d0
        enddo
        do l= 6, 9
            gl%eos_coeff%di(l,nrsubst) = 2.d0
        enddo
        do l=10, 12
            gl%eos_coeff%di(l,nrsubst) = 3.d0
        enddo
        gl%eos_coeff%di(13,nrsubst) = 4.d0
        do l=14, 15
            gl%eos_coeff%di(l,nrsubst) = 5.d0
        enddo
        gl%eos_coeff%di(16,nrsubst) = 6.d0
        do l=17, 18
            gl%eos_coeff%di(l,nrsubst) = 7.d0
        enddo
        gl%eos_coeff%di(19,nrsubst) = 8.d0

        do l=1, 6
            gl%eos_coeff%di(l+count,nrsubst) = 2.d0*l
            gl%eos_coeff%di(l+count+1,nrsubst) = 2.d0*l
            count=count+1
        enddo
        gl%eos_coeff%di(32,nrsubst) = 12.d0


        !only pol ti
        gl%eos_coeff%ti(1,nrsubst) = 0.d0
        gl%eos_coeff%ti(2,nrsubst) = 0.5d0
        gl%eos_coeff%ti(3,nrsubst) = 1.d0
        gl%eos_coeff%ti(4,nrsubst) = 2.d0
        gl%eos_coeff%ti(5,nrsubst) = 3.d0
        gl%eos_coeff%ti(6,nrsubst) = 0.d0
        gl%eos_coeff%ti(7,nrsubst) = 1.d0
        gl%eos_coeff%ti(8,nrsubst) = 2.d0
        gl%eos_coeff%ti(9,nrsubst) = 3.d0
        gl%eos_coeff%ti(10,nrsubst) = 0.d0
        gl%eos_coeff%ti(11,nrsubst) = 1.d0
        gl%eos_coeff%ti(12,nrsubst) = 2.d0
        gl%eos_coeff%ti(13,nrsubst) = 1.d0
        gl%eos_coeff%ti(14,nrsubst) = 2.d0
        gl%eos_coeff%ti(15,nrsubst) = 3.d0
        gl%eos_coeff%ti(16,nrsubst) = 2.d0
        gl%eos_coeff%ti(17,nrsubst) = 2.d0
        gl%eos_coeff%ti(18,nrsubst) = 3.d0
        gl%eos_coeff%ti(19,nrsubst) = 3.d0
        gl%eos_coeff%ti(20,nrsubst) = 3.d0
        gl%eos_coeff%ti(21,nrsubst) = 4.d0
        gl%eos_coeff%ti(22,nrsubst) = 3.d0
        gl%eos_coeff%ti(23,nrsubst) = 5.d0
        gl%eos_coeff%ti(24,nrsubst) = 3.d0
        gl%eos_coeff%ti(25,nrsubst) = 4.d0
        gl%eos_coeff%ti(26,nrsubst) = 3.d0
        gl%eos_coeff%ti(27,nrsubst) = 5.d0
        gl%eos_coeff%ti(28,nrsubst) = 3.d0
        gl%eos_coeff%ti(29,nrsubst) = 4.d0
        gl%eos_coeff%ti(30,nrsubst) = 3.d0
        gl%eos_coeff%ti(31,nrsubst) = 4.d0
        gl%eos_coeff%ti(32,nrsubst) = 5.d0

    elseif (gl%eos_coeff%nreg(nrsubst) == 19) then        !Bender
        do l=1, 5
            gl%eos_coeff%di(l,nrsubst) = 1.d0
        enddo
        do l=6, 8
            gl%eos_coeff%di(l,nrsubst) = 2.d0
        enddo
        do l=9, 10
            gl%eos_coeff%di(l,nrsubst) = 3.d0
        enddo
        do l=11, 12
            gl%eos_coeff%di(l,nrsubst) = 4.d0
        enddo
        gl%eos_coeff%di(13,nrsubst) = 5.d0
        do l=14, 16
            gl%eos_coeff%di(l,nrsubst) = 2.d0
        enddo
        do l=17,19
            gl%eos_coeff%di(l,nrsubst) = 4.d0
        enddo


        gl%eos_coeff%ti(1,nrsubst) = 0.d0
        gl%eos_coeff%ti(2,nrsubst) = 1.d0
        gl%eos_coeff%ti(3,nrsubst) = 2.d0
        gl%eos_coeff%ti(4,nrsubst) = 3.d0
        gl%eos_coeff%ti(5,nrsubst) = 4.d0
        gl%eos_coeff%ti(6,nrsubst) = 0.d0
        gl%eos_coeff%ti(7,nrsubst) = 1.d0
        gl%eos_coeff%ti(8,nrsubst) = 2.d0
        gl%eos_coeff%ti(9,nrsubst) = 0.d0
        gl%eos_coeff%ti(10,nrsubst) = 1.d0
        gl%eos_coeff%ti(11,nrsubst) = 0.d0
        gl%eos_coeff%ti(12,nrsubst) = 1.d0
        gl%eos_coeff%ti(13,nrsubst) = 1.d0
        gl%eos_coeff%ti(14,nrsubst) = 3.d0
        gl%eos_coeff%ti(15,nrsubst) = 4.d0
        gl%eos_coeff%ti(16,nrsubst) = 5.d0
        gl%eos_coeff%ti(17,nrsubst) = 3.d0
        gl%eos_coeff%ti(18,nrsubst) = 4.d0
        gl%eos_coeff%ti(19,nrsubst) = 5.d0

    elseif (gl%eos_coeff%nreg(nrsubst) == 12) then        !Starling
        do l=1, 5
            gl%eos_coeff%di(l,nrsubst) = 1.d0
        enddo
        do l=6, 8
            gl%eos_coeff%di(l,nrsubst) = 2.d0
        enddo
        do l=9, 10
            gl%eos_coeff%di(l,nrsubst) = 5.d0
        enddo
        gl%eos_coeff%di(11,nrsubst) = 2.d0
        gl%eos_coeff%di(12,nrsubst) = 4.d0


        gl%eos_coeff%ti(1,nrsubst) = 0.d0
        gl%eos_coeff%ti(2,nrsubst) = 1.d0
        gl%eos_coeff%ti(3,nrsubst) = 3.d0
        gl%eos_coeff%ti(4,nrsubst) = 4.d0
        gl%eos_coeff%ti(5,nrsubst) = 5.d0
        gl%eos_coeff%ti(6,nrsubst) = 0.d0
        gl%eos_coeff%ti(7,nrsubst) = 1.d0
        gl%eos_coeff%ti(8,nrsubst) = 2.d0
        gl%eos_coeff%ti(9,nrsubst) = 1.d0
        gl%eos_coeff%ti(10,nrsubst) = 2.d0
        gl%eos_coeff%ti(11,nrsubst) = 3.d0
        gl%eos_coeff%ti(12,nrsubst) = 3.d0

    endif


    !n2i to ni  (Multiparameter Equations of State (R.Span, 2000) page 26 Eq[3.28] to Eq[3.26])
    if (((gl%Req(nrsubst) /= 1.d0) .or. ((gl%Req(nrsubst) == 1.d0) .and. gl%vir)) .and. (.not.gl%sta_read) .and. (.not.gl%ben_read)) then !(gl%Req(nrsubst) == 1.d0) .and. gl%vir) Exception for mbwr equations in reduced units
        do i=1, gl%eos_coeff%nreg(nrsubst)                                      !
            gl%ncoeff(nrsubst,i) = gl%ncoeff(nrsubst,i) * ((gl%rhored(nrsubst)/gl%Factor)**gl%eos_coeff%di(i,nrsubst)) / &
                (gl%tred(nrsubst)**(gl%eos_coeff%ti(i,nrsubst))) /(gl%Req(nrsubst)/gl%factorrbwr)
        enddo
    elseif ((gl%sta_read) .or. (gl%ben_read)) then
        do i=1, gl%eos_coeff%nreg(nrsubst)
            gl%ncoeff(nrsubst,i) = gl%ncoeff(nrsubst,i) * ((gl%rhored(nrsubst))**gl%eos_coeff%di(i,nrsubst)) / &
                (gl%tred(nrsubst)**(gl%eos_coeff%ti(i,nrsubst))) /gl%Req(nrsubst)
        enddo
    endif

    !term number after integration
    if (gl%bwr_read) then
        gl%eos_coeff%nreg(nrsubst) = 40
        pol = 22
        pol_add = 3
    elseif (gl%ben_read) then
        gl%eos_coeff%nreg(nrsubst) = 22
        pol = 16
        pol_add = 3
    elseif (gl%sta_read) then
        gl%eos_coeff%nreg(nrsubst) = 13
        pol = 11
        pol_add = 1
    endif


    do l=1, pol
        gl%eos_coeff%p_i(l,nrsubst) = 0
    enddo
    do l=(pol+1), gl%eos_coeff%nreg(nrsubst)
        gl%eos_coeff%p_i(l,nrsubst) = 2
    enddo


    !poly-term coefficients transformation
    Do i=1, (pol-pol_add)
        gl%eos_coeff%ni(i,nrsubst) = gl%ncoeff(nrsubst,i)/gl%eos_coeff%di(i,nrsubst)
    enddo


    !gamma calculation
    do i=1, (pol-pol_add)
        gl%eos_coeff%gama(i,nrsubst)=  0.d0
    enddo
    do i=(pol-pol_add+1), gl%eos_coeff%nreg(nrsubst)
        gl%eos_coeff%gama(i,nrsubst)= ((gl%rhoc(nrsubst)/gl%Factor)/gl%gama_bwr(nrsubst))**2.d0
    enddo



    !integrated ni-terms      from Multiparameter Equations of State(R.Span) page 28 / Table 3.4 and 3.5
    !MBWR terms
    if (gl%bwr_read) then
        gl%eos_coeff%ni(20,nrsubst) = gl%ncoeff(nrsubst,20)/(2.d0*gl%eos_coeff%gama(20,nrsubst)) + gl%ncoeff(nrsubst,22)/(2.d0*gl%eos_coeff%gama(20,nrsubst)**2) &
            + gl%ncoeff(nrsubst,24)/gl%eos_coeff%gama(20,nrsubst)**3 + 3.d0*gl%ncoeff(nrsubst,26)/gl%eos_coeff%gama(20,nrsubst)**4 &
            + 12.d0*gl%ncoeff(nrsubst,28)/gl%eos_coeff%gama(20,nrsubst)**5 + 60.d0*gl%ncoeff(nrsubst,30)/gl%eos_coeff%gama(20,nrsubst)**6
        gl%eos_coeff%ni(21,nrsubst) = gl%ncoeff(nrsubst,21)/(2.d0*gl%eos_coeff%gama(21,nrsubst)) + gl%ncoeff(nrsubst,25)/(gl%eos_coeff%gama(21,nrsubst)**3) &
            + 12.d0*gl%ncoeff(nrsubst,29)/(gl%eos_coeff%gama(21,nrsubst)**5) + 60.d0*gl%ncoeff(nrsubst,31)/gl%eos_coeff%gama(21,nrsubst)**6
        gl%eos_coeff%ni(22,nrsubst) = (gl%ncoeff(nrsubst,23)/(2.d0*gl%eos_coeff%gama(22,nrsubst)**2)) + (3.d0*gl%ncoeff(nrsubst,27)/(gl%eos_coeff%gama(22,nrsubst)**4)) &
            + (60.d0*gl%ncoeff(nrsubst,32)/(gl%eos_coeff%gama(22,nrsubst)**6))
        gl%eos_coeff%ni(23,nrsubst) = -gl%eos_coeff%ni(20,nrsubst)
        gl%eos_coeff%ni(24,nrsubst) = -gl%eos_coeff%ni(21,nrsubst)
        gl%eos_coeff%ni(25,nrsubst) = -gl%eos_coeff%ni(22,nrsubst)
        gl%eos_coeff%ni(26,nrsubst) = -gl%ncoeff(nrsubst,22)/(2.d0*gl%eos_coeff%gama(26,nrsubst)) - gl%ncoeff(nrsubst,24)/(gl%eos_coeff%gama(26,nrsubst)**2) &
            - 3.d0*gl%ncoeff(nrsubst,26)/(gl%eos_coeff%gama(26,nrsubst)**3) - 12.d0*gl%ncoeff(nrsubst,28)/(gl%eos_coeff%gama(26,nrsubst)**4) &
            - 60.d0*gl%ncoeff(nrsubst,30)/(gl%eos_coeff%gama(26,nrsubst)**5)
        gl%eos_coeff%ni(27,nrsubst) = -gl%ncoeff(nrsubst,25)/(gl%eos_coeff%gama(27,nrsubst)**2) - 12.d0*gl%ncoeff(nrsubst,29)/(gl%eos_coeff%gama(27,nrsubst)**4) &
            - 60.d0*gl%ncoeff(nrsubst,31)/(gl%eos_coeff%gama(27,nrsubst)**5)
        gl%eos_coeff%ni(28,nrsubst) = -gl%ncoeff(nrsubst,23)/(2.d0*gl%eos_coeff%gama(28,nrsubst)) - 3.d0*gl%ncoeff(nrsubst,27)/(gl%eos_coeff%gama(28,nrsubst)**3) &
            - 60.d0*gl%ncoeff(nrsubst,32)/(gl%eos_coeff%gama(28,nrsubst)**5)
        gl%eos_coeff%ni(29,nrsubst) = -gl%ncoeff(nrsubst,24)/(2.d0*gl%eos_coeff%gama(29,nrsubst)) - 3.d0*gl%ncoeff(nrsubst,26)/(2.d0*gl%eos_coeff%gama(29,nrsubst)**2) &
            - 6.d0*gl%ncoeff(nrsubst,28)/(gl%eos_coeff%gama(29,nrsubst)**3) - 30.d0*gl%ncoeff(nrsubst,30)/(gl%eos_coeff%gama(29,nrsubst)**4)
        gl%eos_coeff%ni(30,nrsubst) = -gl%ncoeff(nrsubst,25)/(2.d0*gl%eos_coeff%gama(30,nrsubst)) - 6.d0*gl%ncoeff(nrsubst,29)/(gl%eos_coeff%gama(30,nrsubst)**3) &
            - 30.d0*gl%ncoeff(nrsubst,31)/(gl%eos_coeff%gama(30,nrsubst)**4)
        gl%eos_coeff%ni(31,nrsubst) = -3.d0*gl%ncoeff(nrsubst,27)/(2.d0*gl%eos_coeff%gama(31,nrsubst)**2) - 30.d0*gl%ncoeff(nrsubst,32)/(gl%eos_coeff%gama(31,nrsubst)**4)
        gl%eos_coeff%ni(32,nrsubst) = -gl%ncoeff(nrsubst,26)/(2.d0*gl%eos_coeff%gama(32,nrsubst)) - 2.d0*gl%ncoeff(nrsubst,28)/(gl%eos_coeff%gama(32,nrsubst)**2) &
            - 10.d0*gl%ncoeff(nrsubst,30)/(gl%eos_coeff%gama(32,nrsubst)**3)
        gl%eos_coeff%ni(33,nrsubst) = -2.d0*gl%ncoeff(nrsubst,29)/(gl%eos_coeff%gama(33,nrsubst)**2) - 10.d0*gl%ncoeff(nrsubst,31)/(gl%eos_coeff%gama(33,nrsubst)**3)
        gl%eos_coeff%ni(34,nrsubst) = -gl%ncoeff(nrsubst,27)/(2.d0*gl%eos_coeff%gama(34,nrsubst)) - 10.d0*gl%ncoeff(nrsubst,32)/(gl%eos_coeff%gama(34,nrsubst)**3)
        gl%eos_coeff%ni(35,nrsubst) = -gl%ncoeff(nrsubst,28)/(2.d0*gl%eos_coeff%gama(35,nrsubst)) - 5.d0*gl%ncoeff(nrsubst,30)/(2.d0*(gl%eos_coeff%gama(35,nrsubst)**2))
        gl%eos_coeff%ni(36,nrsubst) = -gl%ncoeff(nrsubst,29)/(2.d0*gl%eos_coeff%gama(36,nrsubst)) - 5.d0*gl%ncoeff(nrsubst,31)/(2.d0*(gl%eos_coeff%gama(36,nrsubst)**2))
        gl%eos_coeff%ni(37,nrsubst) = -5.d0*gl%ncoeff(nrsubst,32)/(2.d0*(gl%eos_coeff%gama(37,nrsubst)**2))
        gl%eos_coeff%ni(38,nrsubst) = -gl%ncoeff(nrsubst,30)/(2.d0*gl%eos_coeff%gama(38,nrsubst))
        gl%eos_coeff%ni(39,nrsubst) = -gl%ncoeff(nrsubst,31)/(2.d0*gl%eos_coeff%gama(39,nrsubst))
        gl%eos_coeff%ni(40,nrsubst) = -gl%ncoeff(nrsubst,32)/(2.d0*gl%eos_coeff%gama(40,nrsubst))

        !Bender terms
    elseif (gl%ben_read) then
        gl%eos_coeff%ni(14,nrsubst) = gl%ncoeff(nrsubst,14)/(2.d0*gl%eos_coeff%gama(14,nrsubst)) + gl%ncoeff(nrsubst,17)/(2.d0*(gl%eos_coeff%gama(14,nrsubst)**2))
        gl%eos_coeff%ni(15,nrsubst) = gl%ncoeff(nrsubst,15)/(2.d0*gl%eos_coeff%gama(15,nrsubst)) + gl%ncoeff(nrsubst,18)/(2.d0*(gl%eos_coeff%gama(15,nrsubst)**2))
        gl%eos_coeff%ni(16,nrsubst) = gl%ncoeff(nrsubst,16)/(2.d0*gl%eos_coeff%gama(16,nrsubst)) + gl%ncoeff(nrsubst,19)/(2.d0*(gl%eos_coeff%gama(16,nrsubst)**2))
        gl%eos_coeff%ni(17,nrsubst) = -gl%eos_coeff%ni(14,nrsubst)
        gl%eos_coeff%ni(18,nrsubst) = -gl%eos_coeff%ni(15,nrsubst)
        gl%eos_coeff%ni(19,nrsubst) = -gl%eos_coeff%ni(16,nrsubst)
        gl%eos_coeff%ni(20,nrsubst) = -gl%ncoeff(nrsubst,17)/(2.d0*gl%eos_coeff%gama(20,nrsubst))
        gl%eos_coeff%ni(21,nrsubst) = -gl%ncoeff(nrsubst,18)/(2.d0*gl%eos_coeff%gama(21,nrsubst))
        gl%eos_coeff%ni(22,nrsubst) = -gl%ncoeff(nrsubst,19)/(2.d0*gl%eos_coeff%gama(22,nrsubst))

        !Starling terms
    elseif (gl%sta_read) then
        gl%eos_coeff%ni(11,nrsubst) = gl%ncoeff(nrsubst,11)/(2.d0*gl%eos_coeff%gama(11,nrsubst)) + gl%ncoeff(nrsubst,12)/(2.d0*(gl%eos_coeff%gama(11,nrsubst)**2))
        gl%eos_coeff%ni(12,nrsubst) = -gl%eos_coeff%ni(11,nrsubst)
        gl%eos_coeff%ni(13,nrsubst) = -gl%ncoeff(nrsubst,12)/(2.d0*gl%eos_coeff%gama(13,nrsubst))
    endif


    !integrated di coeff
    if (.not.gl%sta_read) then
        do l=(pol-2), (pol+3)
            gl%eos_coeff%di(l,nrsubst) = 0.d0
        enddo
        do l=(pol+4), (pol+6)
            gl%eos_coeff%di(l,nrsubst) = 2.d0
        enddo
        do l=29, 31
            gl%eos_coeff%di(l,nrsubst) = 4.d0
        enddo

        do l=32, 34
            gl%eos_coeff%di(l,nrsubst) = 6.d0
        enddo

        do l=35, 37
            gl%eos_coeff%di(l,nrsubst) = 8.d0
        enddo
        do l=38, 40
            gl%eos_coeff%di(l,nrsubst) = 10.d0
        enddo
    else
        gl%eos_coeff%di(11,nrsubst) = 0.d0
        gl%eos_coeff%di(12,nrsubst) = 0.d0
        gl%eos_coeff%di(13,nrsubst) = 2.d0
    endif


    !integrated ti coeff
    if (.not.gl%sta_read) then
        if (gl%bwr_read) then
            loop = 7
            count = pol-pol_add
        elseif (gl%ben_read) then
            loop = 3
            count = pol
        endif
        do l=1, loop
            gl%eos_coeff%ti((l+count),nrsubst) = 3.d0
            gl%eos_coeff%ti((l+count+1),nrsubst) = 4.d0
            gl%eos_coeff%ti((l+count+2),nrsubst) = 5.d0
            count = count+2
        enddo
    else
        do l=11, 13
            gl%eos_coeff%ti(l,nrsubst) = 3.d0
        enddo
    endif



    end subroutine

    end module sub_file_input_helper