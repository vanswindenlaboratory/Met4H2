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
    
    module trend_calc_m

    use module_all_types
    use calc_functions
    use calc_func_ptr
    use prop_unit_m
    use calctype_setter_m
    use unit_convertion
    use controlling

    implicit none

    contains
    !-------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------
    subroutine TREND_CALC_INTERFACE(T,D,nrsubst,prop_list,results,prop_name,errorflag,gl_handle_arr)

    type(type_gl), pointer:: gl
    type(c_ptr), intent(inout) :: gl_handle_arr
    integer, intent(in):: nrsubst
    integer, dimension(100), intent(in) :: prop_list
    double precision, dimension(100):: results
    character(30), dimension(100), intent(out) :: prop_name
    double precision, intent(in) :: T,D
    integer:: loop_props,errorflag
    integer, dimension(20):: calc_three_args !(gl, T,  nrsubst)
    integer, dimension(20):: calc_int_args !(gl, T, integer, nrsubst)
    integer, dimension(10):: calc_fld_nr !(gl, fluidnr)
    INTEGER :: fluidnr !enter fluid in nrsubst
    integer, dimension(100):: calc_not_exist
    integer:: loop_size,max_props
    double precision :: wm_mix
    double precision , dimension(30) :: x_spec
    integer :: converttype
    ! to know the convert types for each calctypes
    type(calctype_t) :: c


    if(.not.c_associated(gl_handle_arr))then
        errorflag = -5122  !no gl_handle: return
        return
    end if

    !get info about gl
    CALL C_F_POINTER(gl_handle_arr,gl)

    !spinprop put  it into D - it is converted
    max_props = 81 !no dspin support
    calc_three_args = 0
    calc_int_args = 0
    calc_not_exist = 0
    calc_fld_nr = 0
    results = 0.d0
    calc_three_args(1:15) = (/12,13,14,41,57,58,73,74,75,76,77,78,79,80,81/) !calc functions have only 3 arguments
    calc_not_exist(1:16) = (/1,2,15,17,18,36,37,38,39,41,45,51,52,53,54,64/) !this calctypes are not available
    calc_fld_nr(1:3) = (/68,69,70/) !tboyle etc.
    !calc_int_args(1) = 64 !DSPIN (spinprop is integer) !no dspin support (only single available)
    loop_size = count(prop_list/=0)

    !unitin = get_unitin_m(unitdefinition)
    !if(unitin < 0) then
    !    errorflag = unitin
    !    return
    !end if

    ! If specific input -> convert the density

    x_spec = 0d0

    if(gl%unitin .eq. 2) then
        
        if(gl%trend_calc_called) then
            converttype = 1
            call convert_fractions(gl, converttype, 1234.d0,  x_spec)
        end if

        call convert_inputprop(gl, 1, D)

        call wm_mix_calc(gl, wm_mix)
        
        if(.not.gl%trend_calc_called) then
            call convert_fractions(gl, gl%unitin, wm_mix, gl%molfractions)
        end if
        
        gl%already_converted = 0
    else
        wm_mix = 1d0
    end if

    do loop_props = 1,loop_size



        !check if the prop exits at the moment 64 possible positions (not all are set)
        if( (any(calc_not_exist==prop_list(loop_props)).and.prop_list(loop_props).ne.0) &
            .or.prop_list(loop_props)>max_props                                   &
            .or.prop_list(loop_props)<0               &                                       
            .or.prop_list(loop_props)==0) then

            if(prop_list(loop_props) == 1) then
                results(loop_props) = T
            elseif(prop_list(loop_props) == 2) then
                results(loop_props) = D * wm_mix
            else
                results(loop_props) = -5000
            end if

            if(prop_list(loop_props) > 0) then
                prop_name(loop_props) = act_func_ptrs%function_ptrs(prop_list(loop_props))%propname
            else
                prop_name(loop_props) = 'Invalid proplist'
                results(loop_props) = -5111
            end if
        else

            if(any(calc_int_args==prop_list(loop_props))) then
                results(loop_props) = act_func_ptrs%function_ptrs(prop_list(loop_props))%ptr3 (gl,T,INT(D),nrsubst)

            elseif(any(calc_three_args==prop_list(loop_props))) then
                results(loop_props) = act_func_ptrs%function_ptrs(prop_list(loop_props))%ptr2 (gl,T,nrsubst)
            elseif(any(calc_fld_nr==prop_list(loop_props))) then
                fluidnr = nrsubst
                results(loop_props) = act_func_ptrs%function_ptrs(prop_list(loop_props))%ptr4 (gl, fluidnr)
            else
                results(loop_props) = act_func_ptrs%function_ptrs(prop_list(loop_props))%ptr (gl,T,D,nrsubst)
            end if

            prop_name(loop_props) = act_func_ptrs%function_ptrs(prop_list(loop_props))%propname

            if(gl%unitin .ne. 1) then
                gl%converttype = c%converttype(prop_list(loop_props))
                call convert_units_out(gl, prop_list(loop_props), gl%unitdefinition, errorflag, gl%molfractions, results(loop_props))
            end if

        end if



    end do

    gl%trend_calc_called = .true.


    end subroutine


    !-------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------
    subroutine CONTROL_FLUID_handle_INTERFACE (gl,input,fluids, moles, EOS_indicator, MIX_indicator, path, unitdefinition,errorflag, gl_handle_control, gl_handle_arr)


    logical :: fluid_present
    integer:: i,nr_fld,nr_fld_prev
    integer :: unitin
    double precision :: prop1, prop2, wm_mix
    character (12) :: input,unit
    character (30), dimension(30) :: fluids
    double precision, dimension (30) :: moles
    character (255) :: path
    integer, dimension (30) :: EOS_indicator
    integer :: MIX_indicator
    integer :: errorflag
    character (20) :: unitdefinition
    integer, intent(inout):: gl_handle_control      ! 1 or 0 (1 -> create fluid, 0 -> destroy fluid)
    type(type_gl), pointer:: gl
    type(c_ptr), intent(inout):: gl_handle_arr

    errorflag = 0
    ! set input to td& to irgnore bounds
    prop1 = 1d0
    prop2 = 1d0

    ! No fluid is associated with the handle and a new one should be created
    if(gl_handle_control == 1 .and. .not. c_associated(gl_handle_arr))  then

        ! Set all input chars to lower chars
        !call uppertolower_char_list((/unitdefinition,path,input/))
        call uppertolower_char(unitdefinition,len(unitdefinition))
        call uppertolower_char(path,len(path))
        call uppertolower_char(input,len(input))


        ! Get the unitin form unitdefinition
        unitin = get_unitin_m(unitdefinition)
        if(unitin < 0) then
            errorflag = unitin
            return
        end if

        call control_fluids(gl,input,fluids,moles,EOS_indicator,unitin,gl_handle_arr,errorflag)
        call set_func_pointers()

        gl%unitin = unitin
        gl%unitdefinition = unitdefinition

        !dummy setup call to guarantte memory allocation
        call setup (gl,input, prop1, prop2, fluids, moles, path, EOS_indicator, MIX_indicator, errorflag)

        gl%trend_calc_called = .false.

        ! A fluid is associated with the handle and the fluid should be destroyed
    elseif(gl_handle_control == 0 .and. c_associated(gl_handle_arr)) then
        call DESTROY_FLUID (gl_handle_arr)
    end if


    end subroutine



    end module
