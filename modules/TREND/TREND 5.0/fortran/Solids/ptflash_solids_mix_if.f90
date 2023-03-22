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
    !Cite as: Span, R.; Beckmuller, R.; Eckermann, T.; Herrig, S.; Hielscher, S.; 
	!          Jager, A.; Mickoleit, E.; Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (2020): 	
    !          TREND. Thermodynamic Reference and Engineering Data 5.0. 
    !          Lehrstuhl fur Thermodynamik, Ruhr-Universitat Bochum.

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

    ! module for file ptflash_solids_mix.f90
    module ptflash_solids_mix_module
    !global use inclusion
    use module_all_types


    interface

    module subroutine ptflash_solid_NC_4P(gl,press, Temp, x_known, rho, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
        & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iFlash, iphase, Phasefrac_0, errval, iter)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp
    double precision, dimension(30):: x_known
    double precision, dimension(4):: rho
    double precision, dimension(30):: x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_hyd, x_sol
    double precision:: rhofluid1_est, rhofluid2_est, rhofluid3_est, rhofluid4_est
    double precision, dimension(4) :: Phasefrac
    integer:: errval
    integer:: iter
    integer:: iFlash
    integer:: iphase
    integer:: Phasefrac_0
    end subroutine ptflash_solid_NC_4P



    !**************************************************************************
    module subroutine SysOfEqs_solid_NC_4P (gl,press, Temp, x_known, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, rhofluid1_est, &
        & rhofluid2_est, rhofluid3_est, rhofluid4_est, Phasefrac, iphase_given, nr_of_fluidphases, &
        & fug_guest, rho, GibbsEQN, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp
    double precision, dimension(30):: x_known, x_fluid1, x_fluid2, x_fluid3, x_fluid4
    double precision, dimension(30):: x_sol, x_hyd
    double precision, dimension(60):: GibbsEQN
    double precision, dimension(30):: fug_guest
    double precision:: rhofluid1_est, rhofluid2_est, rhofluid3_est, rhofluid4_est
    double precision, dimension(4):: Phasefrac, rho
    integer:: errval, iphase_given, nr_of_fluidphases


    end subroutine SysOfEqs_solid_NC_4P


    
    module subroutine Jacobi_solid_NC_4P (gl,press, Temp, x_known, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd, &
        & Phasefrac, iFlash, nr_of_fluidphases, fug_guest, mapp_index, JacMatrix, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: Temp, press
    double precision, dimension(30):: x_known, x_fluid1, x_fluid2, x_fluid3, x_fluid4, x_sol, x_hyd
    double precision, dimension(60, 60):: JacMatrix
    double precision, dimension(4):: Phasefrac
    integer:: iFlash, nr_of_fluidphases
    integer, dimension(3):: mapp_index
    integer:: errval
    double precision, dimension(30):: fug_guest
    end subroutine Jacobi_solid_NC_4P

    module subroutine ptflash_solid_NC_3P(gl,press, Temp, x_known, rho, x_fluid1, x_fluid2, x_sol, x_hyd, rhofluid1_est, &
        & rhofluid2_est, Phasefrac, iFlash, iphase, iter, errval)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30)::  x_known
    double precision:: rhofluid1_est, rhofluid2_est
    integer:: iFlash
    integer:: iPhase
    double precision:: press, temp
    double precision, dimension(30):: x_fluid1, x_fluid2, x_sol, x_hyd
    double precision, dimension(3):: Phasefrac
    double precision, dimension(3) :: rho
    integer:: errval, iter
    end subroutine



    module subroutine SysOfEqs_solid_NC_3P(gl,P, T, x_known, x_fluid1, x_fluid2, x_sol, x_hyd, rhofluid1_est, &
        & rhofluid2_est, fug_save, twofluidphases, iFlash, iphaseIn, Phasefrac, GibbsEQN, errval)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30):: x_hyd
    double precision:: P, T
    double precision, dimension(30):: x_known, x_fluid1, x_fluid2, x_sol
    double precision:: rhofluid1_est, rhofluid2_est
    double precision, dimension(3):: Phasefrac
    integer:: iFlash
    logical:: twofluidphases
    integer:: iPhaseIn
    double precision, dimension(30):: fug_save

    double precision, dimension(60):: GibbsEQN
    integer:: errval
    end subroutine

    
    module subroutine Jacobi_solid_NC_3P (gl,P, T, x_fluid1, x_fluid2, x_sol, x_hyd, fug_save, twofluidphases, iFlash, Phasefrac, JacMatrix, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: T, p
    double precision, dimension(30):: x_fluid1, x_fluid2, x_sol, x_hyd
    double precision, dimension(3):: Phasefrac
    integer :: iFlash
    logical :: twofluidphases
    double precision, dimension(60, 60), intent(out):: JacMatrix
    integer :: errval
    double precision, dimension(30):: fug_save
    end subroutine

    
    module subroutine ptflash_solid_NC_2P (gl,press, Temp, rho, x_known, x_solid, x_hyd, x_fluid, rhofluid_est, beta_loc, iFlash, &
        & iPhase, errval, iter)

    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, beta_loc
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_fluid, x_solid, x_hyd! oder doch kein intent?!
    double precision, dimension(30):: z, x_fluid_new!, occup_ls, occup_ld, occup_sms, occup_smd
    double precision:: rhofluid_est
    integer:: iFlash
    integer:: errval, iter
    integer:: iPhase
    double precision, dimension(3) :: rho
    end subroutine

    
    module subroutine SysOfEqs_solid_NC_2P (gl,P, T, x_fluid, x_solid, x_hyd, rhofluid_est, beta_loc, fug_save, &
        & iFlash, iPhaseIn, GibbsEQN, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: P, T, beta_loc
    double precision, dimension(30):: x_fluid, x_solid
    double precision, dimension(30):: x_hyd
    double precision, dimension(60):: GibbsEQN
    double precision:: rhofluid_est
    integer:: iFlash
    integer:: errval
    integer:: iPhaseIn
    double precision, dimension(30):: fug_save
    end subroutine

    
    module subroutine Jacobi_solid_NC_2P (gl,P, T, x_fluid, x_solid , x_hyd, beta_loc, fug_save, iFlash, JacMatrix, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: T, p, beta_loc
    double precision, dimension(30):: x_fluid, x_solid, x_hyd
    double precision, dimension(60, 60):: JacMatrix
    integer:: errval
    integer:: iFlash
    integer:: i, j, k, h
    double precision:: rho_fluid
    double precision, dimension(30, 30):: dChempoti_dxj
    double precision, dimension(30):: dChempoti_dT, d2nadnidT, dChempoti_dp, z, dnadni
    double precision:: rhoredmix_orig, tredmix_orig, Rmix, dummy
    double precision, dimension(30):: fug_save

    end subroutine

    end interface
    
    contains

    end module ptflash_solids_mix_module
