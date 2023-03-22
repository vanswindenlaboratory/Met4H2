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
    !          TREND. Thermodynamic Reference and Engineering Data 5.0.s
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

    module unit_convertion

    use module_all_types
    use setup_module
    use seawater_module
    use electrolytes

    implicit none

    contains



    subroutine convert_flash_spec(gl, prop_overall,wm_mix)
    !converting flash properties for output handle
    !Benedikt May 2018
    !
    implicit none

    type(type_gl) :: gl

    integer :: i, j, nrofphases, IS_ERROR
    double precision :: sum_weight
    double precision :: wm_mix, wm_water
    double precision, dimension (30) :: prop_overall
    double precision, dimension (30,5) :: prop_phase, x_phase
    character(255) :: errorcodes

    if(gl%seawater) then
        wm_water = gl%wm(1)
        gl%wm(1) = gl%sea%wm_sea
        call wm_mix_calc(gl, wm_mix)
    elseif(gl%el_present) then
        wm_water = gl%wm(1)
        gl%wm(1) = wm_brine(gl)
        call wm_mix_calc(gl, wm_mix)
    end if

    do i=3, 10
        !if (IS_ERROR(int(prop_overall(i))) == 1).and.((prop_overall(i) < -1000.d0) .and. (dabs(prop_overall(i) - int(prop_overall(i))) < 1.d-14))) then
        !    prop_overall(i) = prop_overall(i)
        if(i == 3) then
            prop_overall(i) = prop_overall(i) * wm_mix
        elseif (( i .gt. 3 ) .and. (i .lt. 11)) then
            prop_overall(i) = prop_overall(i) / wm_mix
        end if
    end do

    if(gl%seawater .or. gl%el_present) then
        gl%wm(1) = wm_water
    end if


    end subroutine
    !***************************************************************************************************
    !***************************************************************************************************


    !
    !***************************************************************************************************
    subroutine convert_fractions(gl, converttype, wm_phase, x_spec, n_zero)
    !converting phase composition form molar to specific or from specific to molar
    !if converttype==1, then convert specific to molar; if converttype ==2: convert molar to specific
    implicit none

    type(type_gl) :: gl

    integer :: i,j, converttype, n_zero_comp
    double precision :: wm_mix, sum_weight, wm_phase      !Mixed molar weight for (phase) composition
    ! double precision, dimension (5) :: wm_phase
    double precision, dimension (30) :: x_molar, wm_comp
    double precision, dimension (30) :: x_spec
    integer, dimension(gl%ncomp) :: loc_comp
    integer, optional :: n_zero
    !integer, allocatable :: zero_pos(:)
    ! double precision, dimension (30,5) :: x_phase

    !call wm_mix_calc(gl, wm_mix)
    sum_weight = 0.d0
    x_molar =0.d0
    wm_comp =0.d0
    loc_comp = 0.d0
    !if(present(n_zero))then
    !    if(.not. allocated(zero_pos)) allocate(zero_pos(n_zero))
    !    n_zero_comp = n_zero
    !    zero_pos = findloc(x_spec(:), 0.d0)
    !end if
    if ( (gl%zero_comp) .and. (present(n_zero)) ) then
        j = 1
        do i=1,gl%ncomp+n_zero
            if(x_spec(i) .ge. 1.d-14) then
                loc_comp(j) = i
                j = j + 1
            end if
        end do
    end if

    if(converttype == 1) then    !convert specific input to moles
        if( (gl%zero_comp) .and. (present(n_zero)) ) then
            do i=1, gl%ncomp
                wm_comp(loc_comp(i)) = gl%wm(i)
            end do
            do i=1,gl%ncomp+n_zero                                                                     !this loop will calculate the the moles from specific input
                if(wm_comp(i) .le. 1.d-10) then
                    sum_weight = sum_weight
                else
                    sum_weight = sum_weight + x_spec(i) * (wm_comp(i)**(-1.d0))
                end if
            end do

            do i=1,gl%ncomp+n_zero
                if(wm_comp(i) .le. 1.d-10) then
                    x_spec(i) = 0.d0
                else
                    x_molar(i) = (x_spec(i) * (wm_comp(i)**(-1.d0)))  / (sum_weight)
                    !gl%molfractions(i) = x_molar(i)
                    x_spec(i) = x_molar(i)
                end if
            end do
        else
            do i=1,gl%ncomp                                                                     !this loop will calculate the the moles from specific input
                sum_weight = sum_weight + gl%molfractions(i) * (gl%wm(i)**(-1.d0))
            end do

            do i=1,gl%ncomp
                x_molar(i) = (gl%molfractions(i) * (gl%wm(i)**(-1.d0)))  / (sum_weight)
                gl%molfractions(i) = x_molar(i)
                x_spec(i) = x_molar(i)
            end do
        end if

    elseif (converttype == 2) then  !convert calculated molar compositopns from TREND to specific output
        !call wm_mix_calc(gl, wm_mix)
        if ( (gl%zero_comp) .and. (present(n_zero)) ) then
            do i=1, gl%ncomp
                wm_comp(loc_comp(i)) = gl%wm(i)
            end do
            do i = 1,gl%ncomp+n_zero
                if(wm_comp(i) .le. 1.d-10) then
                    x_molar(i) = 0.d0
                else
                    x_molar(i) = x_spec(i) * (wm_comp(i) / wm_phase)
                end if
            end do
        else
            do i=1, gl%ncomp
                x_molar(i) = x_spec(i) * (gl%wm(i) / wm_phase)
            end do
        end if
        x_spec = x_molar
    end if

    end subroutine
    !************************************************************************************************************
    !***************************************************************************************************
    !
    !
    !***************************************************************************************************
    subroutine convert_inputprop(gl, inputflag, prop)
    !converting specific input to molar props for TREND calculations
    !Benedikt May 2018
    !
    !

    implicit none

    type(type_gl) :: gl

    integer :: inputflag
    double precision :: wm_mix
    double precision :: prop

    if(gl%already_converted == 1) then
        prop = prop
    else

        call wm_mix_calc(gl, wm_mix)
        gl%already_converted = 1
        if (inputflag == 1) then        ! density input kg/m^3 --> mol / m^3
            prop = prop / wm_mix
        elseif(inputflag == 2) then     ! entropy, enthalpy,....
            prop = prop * wm_mix
        end if
    end if
    !futher routines for estimated densities in flash expert?!

    end subroutine
    !**************************************************************************************************
    !**************************************************************************************************


    !***************************************************************************************************
    subroutine convert_single_prop(gl, inputflag, prop, nrsubst)
    !converting single prop molar<--> specific for single component
    !Benedikt June 2018


    implicit none

    type(type_gl) :: gl

    integer :: inputflag, nrsubst
    double precision :: wm, wm_sea
    double precision :: prop


    if(nrsubst .eq. 9999)then
        !wm_sea = molar_sea(gl)
        wm_sea = gl%wm(1)
        !wm_sea = 0.0314038218d0
        ! wm_Sea = 1.d0
        if (inputflag == 1) then        ! density input kg/m^3 --> mol / m^3
            prop = prop / wm_sea
        elseif(inputflag == 2) then     ! entropy, enthalpy,....
            prop = prop * wm_sea
        end if
    else
        if (inputflag == 1) then        ! density input kg/m^3 --> mol / m^3
            prop = prop / gl%wm(nrsubst)
        elseif(inputflag == 2) then     ! entropy, enthalpy,....
            prop = prop * gl%wm(nrsubst)
        end if
    end if
    end subroutine
    !***************************************************************************************************
    !***************************************************************************************************



    !***************************************************************************************************
    subroutine convert_units_out(gl, calctype_internal, unitdefinition, errorflag, moles, TREND_EOS_RESULT)
    implicit none
    type(type_gl) :: gl

    integer :: calctype_internal, errorflag
    double precision :: TREND_EOS_RESULT, prop2, wm_mix
    double precision, dimension(30) :: moles
    character(12) :: unitdefinition


    if((calctype_internal .ge. 12) .and. (gl%seawater .or. gl%el_present)) then
        errorflag = -12800.d0
        trend_eos_result = -12800.d0
        return
    end if


    if (errorflag /= 0) then
        TREND_EOS_RESULT = errorflag
    else


        call moles_incheck(gl,moles, errorflag)

        if(errorflag ==0) then
            gl%molfractions = moles



            if(gl%unitin == 1)   then

                TREND_EOS_RESULT =  TREND_EOS_RESULT

            elseif((gl%unitin == 2) .and. (gl%converttype == 1)) then
                if(gl%seawater) then
                    gl%wm(1) = gl%sea%wm_sea
                elseif(gl%el_present) then
                    gl%wm(1) = wm_brine(gl)
                end if
                call wm_mix_calc(gl, wm_mix)

                TREND_EOS_RESULT =  TREND_EOS_RESULT / wm_mix !*1000

            elseif((gl%unitin == 2) .and. (gl%converttype ==2))then
                if(gl%seawater) then
                    gl%wm(1) = gl%sea%wm_sea
                elseif(gl%el_present) then
                    gl%wm(1) = wm_brine(gl)
                end if
                call wm_mix_calc(gl, wm_mix)

                TREND_EOS_RESULT =  TREND_EOS_RESULT * wm_mix !* 1000

            elseif((gl%unitin == 2) .and. (gl%converttype ==3))then

                call wm_mix_calc(gl, wm_mix)

                TREND_EOS_RESULT =  TREND_EOS_RESULT / wm_mix !* 1000

            elseif((gl%unitin == 2) .and. (gl%converttype ==4))then

                call wm_mix_calc(gl, wm_mix)

                TREND_EOS_RESULT =  TREND_EOS_RESULT * (wm_mix)**(-2.d0) !* 1e6

            elseif((gl%unitin == 2) .and. (gl%converttype ==5))then

                call wm_mix_calc(gl, wm_mix)

                TREND_EOS_RESULT =  TREND_EOS_RESULT * (wm_mix)**(-4.d0) !* 1e12

            end if

        else
            TREND_EOS_RESULT = errorflag
        end if



        if(trim(unitdefinition)=="specific") then
            gl%unitin=2
            prop2 = gl%inputprops(2)
            gl%molfractions = gl%moles_read
            moles = gl%molfractions
            !if(gl%seawater) gl%wm(1) = molar_sea(gl)
        end if

    end if
    end subroutine
    !***************************************************************************************************
    !***************************************************************************************************



    !***************************************************************************************************
    subroutine moles_incheck(gl,moles,  errorflag)
    !Subroutine to check validity of input compositions (for example of FLASH_EXPERT or FLASH_NC_3P_EXPERT)
    !convert if neccessary

    implicit none

    type(type_gl) :: gl


    double precision, dimension(30) :: moles, molesave
    integer:: l, e, j, t, i, c, b, k
    integer:: errorflag, converttype
    double precision:: sum
    character(12) :: dd
    double precision :: sum_weight, wm_mix
    double precision, dimension (30) :: x_molar

    errorflag = 0
    i=0

    ! check if all entries are within 0 < x_i < 1 and fill the vector to the size of 30 with zeros
    sum = 0.d0

    do k = 1, gl%ncomp
        if (((0.d0 >= moles(k)) .OR. (moles(k) > 1.d0)) .AND. (k == 1)) then    !check if 0 < xi <= 1
            errorflag = -9951   ! there is a valid molefraction given in the first entry
            return
        end if
        if (((0.d0 > moles(k)) .OR. (moles(k) >= 1.d0)) .AND. (k > 1)) then    !check if 0 < xi < 1
            errorflag = -9951
            return
        end if

        if (moles(k) == 0.d0) then  ! check whether there are any further molefractions given -- >  error
            do c = k, 30
                if (moles(c) /= 0.d0) then
                    errorflag = -9951
                    return
                end if
            end do
            exit
        end if
        sum = sum + moles(k)
    end do

    if(gl%seawater) then
        sum = sum - gl%sea%salinity + moles(gl%ncomp+1)
    end if

    if (abs(sum - 1.d0) > 1.d-10) then  !check if sum /= 1
        errorflag = -9952
        return
    end if

    !if(gl%unitin == 2) then
    !
    !    molesave = gl%molfractions
    !    gl%molfractions = moles
    !    converttype = 1
    !    call wm_mix_calc(gl, wm_mix)
    !    call convert_fractions(gl, converttype, wm_mix, moles)
    !
    !    gl%molfractions = molesave
    !
    !end if



    End subroutine moles_incheck
    !********************************************************************************************
    !********************************************************************************************

    !*******************************************************************************************
    integer function comp_convert_flag(unit_char)
    character(*), intent(in) :: unit_char

    if(trim(unit_char) == "specific") then
        comp_convert_flag = 2
    elseif(trim(unit_char) == "molar") then
        comp_convert_flag = 1
    else
        comp_convert_flag = -123456
        return
    end if

    end function comp_convert_flag
    !*********************************************************************************************




    end module
