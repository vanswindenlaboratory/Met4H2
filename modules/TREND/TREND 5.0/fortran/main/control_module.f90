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
    !!          Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (2020):
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
    !
    !
    module controlling
    use utility_module
    use phase_properties
    contains

    subroutine control_fluids(gl,input,fluids,moles,EOS_indicator,unitin,handle,errorflag)
!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "control_fluids" :: control_fluids

    use module_all_types
    use, intrinsic :: iso_c_binding


    implicit none

    !Handle for gl
    !-------------------------------------------------------------------------------------------------------------------------------------------------
    type(type_gl), pointer, intent(inout) :: gl
    type(c_ptr),intent(inout):: handle
    type(type_phase_properties), dimension(:), pointer :: ph_prop
    type(type_gl), pointer :: gl_check

    character (12) :: input                      ! TD, DT, TP, PT, PH, HP, PS, SP
    character (30), dimension(30) :: fluids, fluids_orig
    double precision, dimension(30) :: moles, moles_orig
    integer, dimension (30) :: EOS_indicator, comp_map, eos_indicator_orig !if molfractions of a mixture are zero this array contains the positions of the original mixture of the fluids that are present
    logical:: fluid_present, solids, lj, holdlimits, zero_comp !zero_comp: false: all specified molfractions are > 1d-14, true: one or more specified molfraction is < 1.d-14
    integer:: errorflag,i,nr_fld,nr_fld_prev, unitin

    !if the handle is no adress, its negative
    if(.not. c_associated(handle)) then !create a new fluid and save the handle
        call CREATE_FLUID(handle)
        fluid_present = .TRUE.
        CALL C_F_POINTER(handle,gl)
        call check_zero_comp(gl, fluids, moles, eos_indicator, comp_map, zero_comp, nr_fld, fluids_orig, moles_orig, eos_indicator_orig, errorflag)
        if(errorflag /=0) return
        gl%unitin = unitin
    elseif (c_associated(handle)) then
        CALL C_F_POINTER(handle,gl)
        call check_zero_comp(gl, fluids, moles, eos_indicator, comp_map, zero_comp, nr_fld, fluids_orig, moles_orig, eos_indicator_orig, errorflag)

        if(errorflag /=0) return
        fluid_present = .TRUE.
        if (index(input,'+') /= 0) then
            solids = .true.
        else
            solids = .false.
        endif
        !
        !if (zero_comp) then
        !    fluids_orig(:) = gl%fluids_zero_orig(:)
        !    moles_orig(:) = gl%moles_zero_orig(:)
        !    eos_indicator_orig(:) = gl%eos_indicator_zero_orig(:)
        !end if

        if (gl%check_solid .neqv. solids) fluid_present = .FALSE.
        if (fluid_present) then
            if ((index(input,'&') /= 0).or.(solids)) then
                holdlimits = .false.
            else
                holdlimits = .true.
            endif
            if (gl%hold_limits .neqv. holdlimits) fluid_present = .FALSE.
            if (fluid_present) then
                if (gl%unitin /= unitin) fluid_present = .FALSE.
                if (fluid_present) then
                    !get number of fluids
                    nr_fld = count(fluids(:).ne.' ')
                    nr_fld_prev = count (gl%components(:).ne.'')
                    if (nr_fld /= nr_fld_prev) fluid_present = .FALSE.
                    if (fluid_present) then
                        !loop over fluids
                        if((.not. gl%seawatercalled) .or. (.not. gl%el_present)) then
                            do i = 1, nr_fld
                                call uppertolower_char(fluids(i),len(fluids(i)))
                                if (sum(gl%mapping) == 0) then
                                    if((gl%components(i).eq.fluids(i)).and.(EOS_indicator(i).eq.gl%Eq_type(i)).and.(nr_fld.eq.nr_fld_prev).and..not.solids) then !same fluid and EQ_Type
                                        gl%same_components = .TRUE.
                                    else
                                        !the actual input differs from the previous connected with the handle
                                        fluid_present = .FALSE.
                                        exit
                                    end if
                                else
                                    if((gl%components(i).eq.fluids(i)).and.(EOS_indicator(i).eq.gl%Eq_type(i)).and.(nr_fld.eq.nr_fld_prev).and..not.solids) then !same fluid and EQ_Type
                                        gl%same_components = .TRUE.
                                    elseif((gl%components(gl%mapping(i)).eq.fluids(i)).and.(EOS_indicator(gl%mapping(i)).eq.gl%Eq_type(i)).and.(nr_fld.eq.nr_fld_prev).and.solids) then !same fluid and EQ_Type
                                        gl%same_components = .TRUE.
                                    else
                                        !the actual input differs from the previous connected with the handle
                                        fluid_present = .FALSE.
                                        exit
                                    end if
                                end if
                            end do
                        endif
                    endif
                endif
            endif
        endif


        if(.not.fluid_present) then
            if (c_associated(gl%s_p%ph_prop_handle)) then
                call DESTROY_PH_PROP(gl%s_p%ph_prop_handle, gl%s_p%nr_checks)
            endif
            if (c_associated(gl%s_p%gl_check_handle)) then
                call DESTROY_FLUID(gl%s_p%gl_check_handle)
            endif
            !destroy the actual handle and fluid
            call DESTROY_FLUID(handle)
            !create a new fluid
            call CREATE_FLUID(handle)
            CALL C_F_POINTER(handle,gl)
            gl%unitin = unitin
        end if

    end if
    if (zero_comp) then
        gl%zero_comp = zero_comp
        if (.not. allocated(gl%comp_map)) allocate(gl%comp_map(nr_fld))
        gl%comp_map(1:nr_fld) = comp_map(1:nr_fld)
        if (.not.allocated(gl%fluids_zero_orig))allocate(gl%fluids_zero_orig(30))
        gl%fluids_zero_orig = ''
        gl%fluids_zero_orig(:) = fluids_orig(:)
        if (.not.allocated(gl%moles_zero_orig))allocate(gl%moles_zero_orig(30))
        gl%moles_zero_orig = 0.d0
        gl%moles_zero_orig(:) = moles_orig(:)
        if (.not.allocated(gl%eos_indicator_zero_orig))allocate(gl%eos_indicator_zero_orig(30))
        gl%eos_indicator_zero_orig = 0
        gl%eos_indicator_zero_orig(:) = eos_indicator_orig(:)
    endif
    !save information if fluid has changed
    gl%s_p%fluid_present = fluid_present
    if(gl%seawatercalled .or. gl%el_present) then
        call destroy_fluid(handle)
        call CREATE_FLUID(handle)
        call c_f_pointer(handle,gl)
    end if
    !end of checking fluid is present or not



    end subroutine


    !Function to create a new fluid
    !The return value of the function is the memory adress associated with the global derived type gl
    !TYPE(C_PTR) function CREATE_FLUID()
    subroutine create_fluid(handle)
    use, intrinsic :: iso_c_binding
    use module_all_types
    implicit none
    type(c_ptr) :: handle
    type(type_gl) , pointer:: gl

    allocate(gl)

    handle = C_LOC(gl)

    end subroutine CREATE_FLUID


    !Subroutine to destroy the created fluid
    subroutine DESTROY_FLUID (handle)
    use, intrinsic :: iso_c_binding
    use module_all_types
    implicit none

    TYPE(C_PTR), intent(in)  :: handle
    type(type_gl),pointer :: gl_ptr

    CALL C_F_POINTER(handle,gl_ptr)

    deallocate(gl_ptr)

    end subroutine DESTROY_FLUID

    !Subroutine to destroy the created fluid
    subroutine DESTROY_PH_PROP (handle, dim)
    use, intrinsic :: iso_c_binding
    use module_all_types
    implicit none

    TYPE(C_PTR), intent(in) :: handle
    type(type_phase_properties), dimension(:),pointer :: ph_prop
    integer :: dim, test

    !vorher nur:
    !CALL C_F_POINTER(handle,ph_prop, [dim])
    !test = C_LOC(ph_prop)
    if (dim .gt. 1) then
        CALL C_F_POINTER(handle,ph_prop, [dim])
    else
        CALL C_F_POINTER(handle,ph_prop)
    end if
    allocate(ph_prop(dim))
    deallocate(ph_prop)

    end subroutine DESTROY_PH_PROP

    subroutine check_zero_comp(gl, fluids, moles, eos_indicator, comp_map, zero_comp, ncomps, fluids_orig, moles_orig, eos_indicator_orig, errorflag)
    use module_all_types

    type(type_gl) :: gl
    character (30), dimension(30) :: fluids, fluids_orig
    double precision, dimension(30) :: moles, moles_orig
    integer, dimension (30) :: EOS_indicator, comp_map, eos_indicator_orig !if molfractions of a mixture are zero this array contains the positions of the original mixture of the fluids that are present
    integer :: ncomps, i, j, errorflag
    logical:: zero_comp !zero_comp: false: all specified molfractions are > 1d-14, true: one or more specified molfraction is < 1.d-14

    errorflag = 0
    ncomps = count(fluids /= '')
    if (count(eos_indicator /= 0) /= ncomps)then
        errorflag = -9953
        return
    endif

    comp_map = 0
    fluids_orig = ''
    moles_orig = 0d0
    eos_indicator_orig = 0
    zero_comp = .false.
    j = 1
    do i = 1,ncomps
        if (moles(i) > 1.d-14) then
            comp_map(j) = i
            j = j + 1
        endif
    end do
    if ((j-1) < ncomps) then

        !Moved to control_fluid
        !store original mixture with molfraction == 0
        !if (.not.allocated(gl%fluids_zero_orig))allocate(gl%fluids_zero_orig(30))
        !gl%fluids_zero_orig = ''
        !gl%fluids_zero_orig(1:j) = fluids(1:j)
        !if (.not.allocated(gl%moles_zero_orig))allocate(gl%moles_zero_orig(30))
        !gl%moles_zero_orig = 0.d0
        !gl%moles_zero_orig(1:j) = moles(1:j)
        !if (.not.allocated(gl%eos_indicator_zero_orig))allocate(gl%eos_indicator_zero_orig(30))
        !gl%eos_indicator_zero_orig = 0
        !gl%eos_indicator_zero_orig(1:j) = eos_indicator(1:j)

        fluids_orig = fluids(1:j)
        moles_orig = moles(1:j)
        eos_indicator_orig = eos_indicator(1:j)


        zero_comp = .true.
        ncomps = j - 1
        fluids(1:ncomps) = fluids(comp_map(1:ncomps))
        fluids(j:30) = ''
        moles(1:ncomps) = moles(comp_map(1:ncomps))
        moles(j:30) = 0.d0
        eos_indicator(1:ncomps) = eos_indicator(comp_map(1:ncomps))
        eos_indicator(j:30) = 0
    endif

    end subroutine check_zero_comp

    subroutine revert_zero_inpt(gl, fluids, moles, eos_indicator)
    use module_all_types

    type(type_gl) :: gl
    character(30), dimension(30) :: fluids
    double precision, dimension(30) :: moles
    integer, dimension(30) :: eos_indicator

    if (allocated(gl%fluids_zero_orig)) fluids = gl%fluids_zero_orig
    if (allocated(gl%moles_zero_orig)) then
        if(gl%zero_spec) then
            moles = gl%moles_zero_orig_spec
        else
            moles = gl%moles_zero_orig
        end if
    end if
    if (allocated(gl%EOS_indicator_zero_orig)) EOS_indicator = gl%EOS_indicator_zero_orig

    end subroutine revert_zero_inpt

    subroutine revert_zero_vector(gl, argument)
    use module_all_types

    type(type_gl) :: gl
    double precision, dimension(30) :: argument, temp

    temp = argument
    argument = 0.d0
    argument(gl%comp_map) = temp(1:gl%ncomp)

    end subroutine revert_zero_vector

    end module
