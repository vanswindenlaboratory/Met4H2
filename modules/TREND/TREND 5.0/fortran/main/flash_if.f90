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

    ! module for file flash.f90
    module flash_module
    !global use inclusion
    use module_all_types

	use module_regula_falsi_support


    interface

    module subroutine Flash_PhaseBoundary(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
        &  Nr_x_given, errval, iter)
    !DEC$ ATTRIBUTES DLLEXPORT :: Flash_PhaseBoundary
        
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision:: rholiq_est, rhovap_est
    integer:: errval, iter, iFlash, Nr_x_given
    End subroutine Flash_PhaseBoundary

    module subroutine Flash_PhaseBoundary_calc(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash,&
        & iPhase_try, Nr_x_given, errval, iter)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter
    end subroutine Flash_PhaseBoundary_calc

    module subroutine SysOfEqs_PhaseBoundary(gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, iPhase_try, GibbsEQN, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: P, T
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, iPhase_try
    integer:: errval
    double precision, dimension(30):: z, lnfi_liq, lnfi_vap
    double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq
    integer:: i, errorflag

    end subroutine SysOfEqs_PhaseBoundary

    module subroutine Jacobi_PhaseBoundary (gl,P, T, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, JacMatrix, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: T, p, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JacMatrix
    integer:: errval
    integer:: iFlash, Nr_x_given
    end subroutine Jacobi_PhaseBoundary

    module subroutine PTX_startvals_PhaseBoundary (gl,P, T, x_known, x_vap, x_liq, vapfrac, iFlash, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: P, T, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    integer:: iFlash
    integer:: errval

    end subroutine PTX_startvals_PhaseBoundary

    module subroutine Flash_pT(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, errval, iter)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision:: rholiq_est, rhovap_est
    integer:: errval, iter
    End subroutine Flash_pT

    module subroutine Flash_pT_calc(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, errval, iter)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(30):: z, x_vap_new, x_liq_new, x_save
    double precision:: rholiq_est, rhovap_est, Delta_X_min
    integer:: iPhase_try
    integer:: errval, iter
    end subroutine Flash_pT_calc

    module subroutine SysOfEqs_pT(gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iPhase_try, GibbsEQN, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: P, T
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est
    integer:: iPhase_try
    integer:: errval
    end subroutine SysOfEqs_pT

    module subroutine Jacobi_pT (gl,P, T, x_vap, x_liq, vapfrac, JacMatrix, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: T, p, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JacMatrix
    integer:: errval
    end subroutine Jacobi_pT

    module subroutine PTX_startvals_pT (gl,P, T, x_known, x_vap, x_liq, vapfrac, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: P, T, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    integer:: errval
    end subroutine PTX_startvals_pT

    module subroutine Flash_ph(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, h_spec, iPhase_try, errval, iter)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(30):: x_vap_new, x_liq_new
    double precision:: rholiq_est, rhovap_est, h_spec
    integer:: iPhase_try
    integer:: errval, iter

    end subroutine Flash_ph

    module subroutine SysOfEqs_ph(gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, h_spec, iPhase_try, GibbsEQN, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: P, T
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est, h_spec
    integer::  iPhase_try
    integer:: errval
    end subroutine SysOfEqs_ph

    module subroutine Jacobi_ph (gl,P, T, x_vap, x_liq, vapfrac, JacMatrix, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: T, p, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JacMatrix
    integer:: errval
    end subroutine Jacobi_ph

    module subroutine Flash_ps(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, s_spec, iPhase_try, errval, iter)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(30):: x_vap_new, x_liq_new
    double precision:: rholiq_est, rhovap_est, s_spec
    integer:: iPhase_try
    integer:: errval, iter
    end subroutine Flash_ps

    module subroutine SysOfEqs_ps(gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, s_spec, iPhase_try, GibbsEQN, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: P, T
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est, s_spec
    integer::  iPhase_try
    integer:: errval
    end subroutine SysOfEqs_ps

    module subroutine Jacobi_ps (gl,P, T, x_vap, x_liq, vapfrac, JacMatrix, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: T, p, vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JacMatrix
    integer:: errval
    end subroutine Jacobi_ps

    module subroutine RachRice (gl,K_Val, x_known, x_vap, x_liq, vapfrac, errorflag)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30) :: x_known, K_val
    double precision, dimension(30) :: x_vap, x_liq
    double precision :: vapfrac
    integer:: errorflag
    End Subroutine RachRice

    double precision module function rac_func(gl,frac, parameters)
    implicit none
    type(type_gl) :: gl
    double precision:: frac
    type(type_additional_parameters) :: Parameters            ! Inputvariables
    End function rac_func

    double precision module function Tsat_iter (gl,P, x, IPhase)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30):: x
    double precision:: P
    integer:: IPhase
    end function Tsat_iter

    Double Precision module Function RacRice_div(gl,temp, Parameters)
    implicit none
    type(type_gl) :: gl
    Double Precision :: temp              ! Inputvariablen
    type(type_additional_parameters) :: Parameters            ! Inputvariables
    End Function RacRice_div

    module subroutine newvars (gl,P, T, x_vap, x_liq, rhovap_est, rholiq_est, iPhase_try, errorflag)
    implicit none
    type(type_gl) :: gl
    double precision:: T, P, rhovap_est, rholiq_est
    integer:: iPhase_try
    double precision, dimension(30):: x_vap, x_liq, x_liq_reac, x_vap_reac
    integer:: errorflag, iPhaseVap, iPhaseLiq
    end subroutine newvars

    module subroutine LUdecomp (gl,n,aMatrix,cMatrix,ierr,herr)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(60, 60):: amatrix
    double precision, dimension(60):: cmatrix,  ctemp, sdecomp
    character(255):: herr
    integer:: ierr, i, j, k, n
    END subroutine LUdecomp

    module subroutine Pivot (gl,n,j,iord,aMatrix,sdecomp)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(60,60)::amatrix
    integer, dimension(60):: iord
    double precision, dimension(60):: sdecomp
    integer:: j, n
    END SUBROUTINE Pivot

    module subroutine Succ_Sub(gl,press, Temp, x_known, x_vap, x_liq, rhovap_est, rholiq_est, vapfrac, iFlash, errval, Nr_of_iter, converged)
    implicit none
    type(type_gl) :: gl
    double precision:: Temp, press, vapfrac
    double precision:: rhovap_est, rholiq_est
    double precision, dimension(30):: x_known
    double precision, dimension(30):: x_vap, x_liq
    integer:: iFlash, Nr_of_iter
    integer:: errval
    logical:: converged
    end subroutine Succ_Sub

    module subroutine Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_A, errorflag)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(60,60):: MatrixA
    integer:: rankA
    double precision, dimension(60):: vectorb
    double precision, dimension(60):: vectorx
    double precision:: det_A
    integer:: errorflag

    end subroutine Gauss_algorithm

    module subroutine Adjugate(gl,MatrixA, rankA, adj_A, errorflag)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(60,60):: MatrixA
    integer:: rankA
    double precision, dimension(60,60):: adj_A
    integer:: errorflag

    end subroutine Adjugate

    module subroutine Trace(gl,MatrixA, rankA, trace_A, errorflag)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(60,60):: MatrixA
    integer:: rankA
    double precision:: trace_A
    integer:: errorflag
    end subroutine Trace

    module subroutine Mat_mult(gl,MatrixA, rowA, colA, MatrixB, rowB, colB, MatrixC, errorflag)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(60,60):: MatrixA
    double precision, dimension(60,60):: MatrixB
    double precision, dimension(60,60):: MatrixC
    integer:: rowA, rowB, colA, colB
    integer:: errorflag
    end subroutine Mat_mult

    module subroutine crit_pt(gl,z_given, T_crit, p_crit, d_crit, nr_crit_pts, errorflag)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30):: z_given
    double precision, dimension(10):: T_crit
    double precision, dimension(10):: p_crit
    double precision, dimension(10):: d_crit
    integer:: nr_crit_pts
    integer:: errorflag
    end subroutine crit_pt

    module subroutine dAdjL_dX(gl,MatrixL, dL_dtau, dL_ddel, rankL, dadjL_dtau, dadjL_ddel, errorflag)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(60,60):: MatrixL
    double precision, dimension(60,60):: dL_dtau
    double precision, dimension(60,60):: dL_ddel
    integer:: rankL
    double precision, dimension(60,60):: dadjL_dtau
    double precision, dimension(60,60):: dadjL_ddel
    integer:: errorflag

    end subroutine dAdjL_dX

    
    module subroutine Flash_PhaseBoundary_Vbased_calc(gl,press, Temp, rhovap, rholiq, x_known, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, errval, iter)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac, rhovap, rholiq, rhovap_new, rholiq_new, GibbsEQN_b!, p_calc
    double precision, dimension(30):: x_vap, x_liq, x_known
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter
    end subroutine Flash_PhaseBoundary_Vbased_calc

    module subroutine SysOfEqs_PhaseBoundary_Vbased(gl,P, T, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: P, T,  p_vap, p_liq,rhovap, rholiq !P_CALC
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, iPhase_try
    integer:: errval
    end subroutine SysOfEqs_PhaseBoundary_Vbased

    module subroutine Jacobi_PhaseBoundary_Vbased (gl,P, T, rhovap, rholiq, x_vap, x_liq,vapfrac, iFlash, Nr_x_given, JT, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: T, p, vapfrac, dPdT_liq, dPdT_vap, dPdrho_liq, dPdrho_vap, rhovap, rholiq
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JT
    integer:: errval
    integer:: iFlash, Nr_x_given
    end subroutine Jacobi_PhaseBoundary_Vbased

    module subroutine Flash_PhaseBoundary_Vbased_calc_LevMar(gl,press, Temp, rhovap, rholiq, x_known, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, errval, iter)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac, rhovap, rholiq, rhovap_new, rholiq_new !, P_CALC
    double precision, dimension(30):: x_vap, x_liq, x_known
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter
    end subroutine Flash_PhaseBoundary_Vbased_calc_LevMar

    module subroutine SysOfEqs_PhaseBoundary_Vbased_LevMar(gl,P, T, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: P, T, p_vap, p_liq,rhovap, rholiq !P_CALC
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, iPhase_try
    integer:: errval
    double precision, dimension(30):: z, lnfi_liq, lnfi_vap
    double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq
    integer:: i, errorflag
    end subroutine SysOfEqs_PhaseBoundary_Vbased_LevMar

    module subroutine Jacobi_PhaseBoundary_Vbased_LevMar(gl,P, T, rhovap, rholiq, x_vap, x_liq,vapfrac, iFlash, Nr_x_given, JT, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: T, p, vapfrac, dPdT_liq, dPdT_vap, dPdrho_liq, dPdrho_vap, rhovap, rholiq
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JT
    integer:: errval
    integer:: iFlash, Nr_x_given
    end subroutine Jacobi_PhaseBoundary_Vbased_LevMar

    
    module subroutine Flash_PhaseBoundary_Vbased_calc_LevMar_mix(gl,press, Temp, rhovap, rholiq, x_known, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, errval, iter, GibbsEQN_b)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac, rhovap, rholiq, rhovap_new, rholiq_new!, P_CALC
    double precision, dimension(30):: x_vap, x_liq, x_known
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter
    double precision:: GibbsEQN_b
    end subroutine Flash_PhaseBoundary_Vbased_calc_LevMar_mix


    module subroutine Flash_PhaseBoundary_Vbased_calc_DogLeg(gl,press, Temp, rhovap, rholiq, x_known, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, errval, iter)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac, rhovap, rholiq, rhovap_new, rholiq_new !, p_calc
    double precision, dimension(30):: x_vap, x_liq, x_known
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter
    end subroutine Flash_PhaseBoundary_Vbased_calc_DogLeg

    module subroutine SysOfEqs_PhaseBoundary_Vbased_DogLeg(gl,P, T, x_vap, x_liq, rhovap, rholiq, vapfrac, iFlash, GibbsEQN, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: P, T, p_vap, p_liq,rhovap, rholiq !P_CALC
    double precision:: vapfrac
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60):: GibbsEQN
    double precision:: rholiq_est, rhovap_est
    integer:: iFlash, iPhase_try
    integer:: errval
    double precision, dimension(30):: z, lnfi_liq, lnfi_vap
    double precision:: rhoredmix_orig, tredmix_orig, d_vap, d_liq
    integer:: i, errorflag
    end subroutine SysOfEqs_PhaseBoundary_Vbased_DogLeg

    module subroutine Jacobi_PhaseBoundary_Vbased_DogLeg(gl,P, T, rhovap, rholiq, x_vap, x_liq,vapfrac, iFlash, Nr_x_given, JT, errval)
    implicit none
    type(type_gl) :: gl
    double precision:: T, p, vapfrac, dPdT_liq, dPdT_vap, dPdrho_liq, dPdrho_vap, rhovap, rholiq
    double precision, dimension(30):: x_vap, x_liq
    double precision, dimension(60, 60):: JT
    integer:: errval
    integer:: iFlash, Nr_x_given
    end subroutine Jacobi_PhaseBoundary_Vbased_DogLeg

    module subroutine Flash_PhaseBoundary_Vbased_calc_DogLeg_mix(gl,press, Temp, rhovap, rholiq, x_known, x_vap, x_liq, vapfrac, iFlash, Nr_x_given, errval, iter)
    implicit none
    type(type_gl) :: gl
    double precision:: press, temp, vapfrac, rhovap, rholiq, rhovap_new, rholiq_new !, p_calc
    double precision, dimension(30):: x_vap, x_liq, x_known
    double precision, dimension(30):: z, x_vap_new, x_liq_new
    integer:: iFlash, Nr_x_given, iPhase_try
    integer:: errval, iter
    end subroutine Flash_PhaseBoundary_Vbased_calc_DogLeg_mix


    end interface

    contains


    end module flash_module
