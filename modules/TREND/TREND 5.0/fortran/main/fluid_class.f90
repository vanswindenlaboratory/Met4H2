    ! FLUID CLASS

    module fluid
    use iso_c_binding
    use IFPORT
    use utility_module
    use module_all_types
    use interface_helper
    use setup_module
    use controlling
    use calc_func_ptr
    implicit none

    type fld
        double precision, dimension(30) :: moles
        character (12) :: input,calctype
        character(20) :: unit
        character (30), dimension(30) :: fluids
        integer, dimension(30) :: eos_ind
        integer:: mix_ind,unitin
        character(255) :: path
        type(c_ptr) :: handle
        integer:: handle_control
    end type

    !************************************************************
    !interface for standard constructor
    interface fld
    module procedure init
    end interface

    !************************************************************
    contains


    !###################################################################
    !constructor function for the fluid
    type(fld) function init(input_in,calctype_in,unit_in,fld_in,moles_in,eos_ind_in,mix_ind_in,path_in)
    !DEC$ ATTRIBUTES DLLEXPORT:: init
    !##################################################################
    type(type_gl), pointer :: gl
    character(12) :: input_in,calctype_in
    character(20) :: unit_in
    character(30), dimension(30) :: fld_in
    double precision, dimension(30) :: moles_in
    integer, dimension(30) :: eos_ind_in
    integer :: mix_ind_in
    character(255) :: path_in
    integer , dimension(3) :: sizes_in
    integer:: errorflag,i
    double precision :: prop1,prop2
    !###################################################################
    ! check if the sizes are correct
    sizes_in(1) = count(fld_in(:).ne.' ')
    sizes_in(2) = count(moles_in(:).ne.0d0)
    sizes_in(3) = count(eos_ind_in(:).ne.0)

    !###################################################################
    !all sizes fit
    if(all(sizes_in(:) .eq. sizes_in(1))) then

        !Set all variables to the class variables
        init%moles = moles_in
        init%input = input_in
        init%calctype = calctype_in
        init%unit = unit_in
        init%mix_ind = mix_ind_in
        init%eos_ind =eos_ind_in
        init%fluids = fld_in
        init%handle = c_null_ptr
        init%path =  path_in
        init%handle_control = 0

    else ! an error occured
        print *,'! ERROR !Moles, Fluid and EOS_ind must have the same size!'
        pause
    end if

    !###################################################################
    !init the fluid
    call uppertolower_char(init%unit,len(init%unit))
    call uppertolower_char(init%calctype,len(init%calctype))
    call uppertolower_char(init%input,len(init%input))
    call uppertolower_char(init%path,len(init%path))

    do i=1,sizes_in(1)
        call uppertolower_char(init%fluids(i),len(init%fluids(i)))
    end do

    init%unitin = get_unitin(init%unit)
    
    call control_fluids(gl,init%input,init%fluids, init%moles &
        & , init%eos_ind ,init%unitin,init%handle,errorflag )

    prop1 = 1d0
    prop2 = 1d0
    call setup(gl,init%input, prop1, prop2, init%fluids, init%moles, init%path, init%eos_ind, init%mix_ind, errorflag)

    call set_func_pointers()


    end function init
    !###################################################################

    end module
    ! END FLUID CLASS