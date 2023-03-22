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

    ! module for file setup.f90
    module setup_module
    !global use inclusion
    use module_all_types

    
! interface for routines with circular dependencies (implemented in submodule in file <this-file-name>_impl)    
    interface
    




    module subroutine setup (gl,input, prop1, prop2, fluids, moles, path, eqtype, mixtype, errorflag)
!DEC$ ATTRIBUTES DLLEXPORT, decorate, alias : "setup" :: setup
    implicit none
    type(type_gl) :: gl
    character (12), intent(inout) :: input                        ! TD, DT, TP, PT, PH, HP, PS, SP
    double precision, intent(inout) :: prop1, prop2               ! input values
    double precision, dimension (30), intent(inout) :: moles      ! molefraction vector
    double precision, dimension (30) :: x_spec
    character (255), intent(inout) :: path                                     ! path where to find fluid and mix files
    integer :: errorflag                        ! information about warnings or errors
    character (30), dimension (30), intent(inout) :: fluids                    ! fluid vector
    integer, dimension(30), intent(inout):: eqtype
    integer, intent(inout) :: mixtype
    end subroutine setup



    module subroutine setpropinput (gl,input, prop1, prop2, errorflag)
    implicit none
    type(type_gl) :: gl
    character (12), intent(inout) :: input                      ! TD, DT, TP, PT, PH, HP, PS, SP
    double precision, intent(inout) :: prop1, prop2             ! input values
    integer, intent(out) :: errorflag                           ! information about warnings or errors
    end subroutine setpropinput




    module subroutine setfluid (gl, fluids, nf, errorflag)
    implicit none
    type(type_gl) :: gl
    integer, intent(out) :: errorflag                       ! information about warnings or errors
    character (30), dimension (30) :: fluids                ! fluid vector
    integer :: nf                                           ! amount of fluid in vector
    end subroutine setfluid


    module subroutine setmoles (gl, moles, nm, errorflag)
    implicit none
    type(type_gl) :: gl
    double precision, dimension (30), intent(inout) :: moles      ! molefraction vector
    integer, intent(out) :: errorflag                        ! information about warnings or errors
    integer:: nm
    end subroutine setmoles


    module subroutine setseawater(gl)
    implicit none
    type(type_gl) :: gl
    end subroutine setseawater


    module subroutine setsolid(gl,fluids, nf, moles, resorted)
    implicit none
    type(type_gl) :: gl
    integer :: nf                                               ! amount of fluid in vector
    character (30), dimension (30) :: fluids                    ! fluid vector
    double precision, dimension (30), intent(out) :: moles      ! molefraction vector
    logical :: resorted                                         ! flag whether fluid list has been resorted  Theresa
    end subroutine setsolid

    module subroutine seteqtype (gl, eqtype, mixtype, ne, errorflag)
    implicit none
    type(type_gl) :: gl
    integer, intent(out) :: errorflag                        ! information about warnings or errors
    integer :: ne                                              ! amount of fluid in vector
    integer, dimension(30), intent(inout) :: eqtype
    integer, intent(inout):: mixtype
    end subroutine seteqtype


    module subroutine setpath (gl,path)
    implicit none
    type(type_gl) :: gl
    character (255), intent(inout) :: path                                     ! path where to find fluid and mix files
    end subroutine setpath



    module subroutine wm_mix_calc(gl,wm_mix)
    implicit none
    type(type_gl) :: gl
    double precision :: wm_mix !Output of Subroutine
    end subroutine wm_mix_calc

    module subroutine wm_mix_reac_calc(gl, wm_comps, molfrac, wm_reac_mix)
    implicit none
    type(type_gl) :: gl
    double precision :: wm_reac_mix
    double precision, dimension(30), intent(in) :: wm_comps, molfrac
    end subroutine wm_mix_reac_calc

    module subroutine R_mix_calc(gl,R_mix)
    implicit none
    type(type_gl) :: gl
    double precision :: R_mix !Output of Subroutine
    end subroutine R_mix_calc


    module subroutine ref_state(gl,errorflag)
    implicit none
    type(type_gl) :: gl
    integer :: errorflag                        ! information about warnings or errors
    end subroutine ref_state

    

    module subroutine check_limits (gl,input, t, p, d, dliq, dvap, nrofphases, x_phase, phasetype, ILIMITS)
    implicit none
    type(type_gl) :: gl
    character(12) :: input
    double precision :: t, d, dliq, dvap, p
    double precision :: x_Phase(30,5)               ! vector containing the compositions of all phases
    integer :: nrofphases
    integer :: phasetype(5)                          !phasetype contains the phase indicator number
    integer :: ILIMITS                           ! states whether the limits of the EOS have been hold, ONE VALUE BACK

    end subroutine check_limits


    module subroutine  set_limits_trans(gl)       ! will not be called with new version of interface
    implicit none
    type(type_gl) :: gl
    end subroutine set_limits_trans

    module subroutine check_limits_trans (gl,input, t, p, d, dliq, dvap, nrofphases, x_phase, phasetype, prop, ILIMITS)
    implicit none
    type(type_gl) :: gl
    character(12) :: input
    double precision :: t, p, d, dliq, dvap
    double precision :: x_Phase(30,5)               ! vector containing the compositions of all phases
    integer :: nrofphases
    integer :: phasetype(5)                          !phasetype contains the phase indicator number
    integer :: ILIMITS                           ! states whether the limits of the EOS have been hold, ONE VALUE BACK
    integer :: prop
    end subroutine check_limits_trans



    module subroutine setup_one_fluid_model(gl)
    implicit none
    type(type_gl) :: gl
    end subroutine

    module subroutine setup_HelmgE(gl, errval)
    implicit none
    type(type_gl) :: gl
    integer :: errval
    end subroutine


    !module subroutine convert_flash_spec(gl, prop_overall,wm_mix)
    !implicit none
    !type(type_gl) :: gl
    !double precision :: wm_mix
    !double precision, dimension (30), intent(inout) :: prop_overall
    !end subroutine
    !
    !module subroutine convert_fractions(gl, converttype, wm_phase, x_spec)!, wm_mix)
    !implicit none
    !type(type_gl) :: gl
    !integer :: converttype
    !double precision :: wm_phase      !Mixed molar weight for (phase) composition
    !double precision, dimension (30), intent(inout) :: x_spec
    !end subroutine
    !
    !module subroutine convert_inputprop(gl, inputflag, prop)
    !implicit none
    !type(type_gl) :: gl
    !integer :: inputflag
    !double precision, intent(inout) :: prop
    !end subroutine
    !
    !module subroutine convert_single_prop(gl, inputflag, prop, nrsubst)
    !implicit none
    !type(type_gl) :: gl
    !
    !integer :: inputflag, nrsubst
    !double precision, intent(inout) :: prop
    !end subroutine



    module subroutine allocate_subtypes(gl, errorflag)
    implicit none
    type(type_gl), intent(inout):: gl
    integer:: errorflag
    end subroutine allocate_subtypes


    end interface
    
    contains

    end module setup_module
