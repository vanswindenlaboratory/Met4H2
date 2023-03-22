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

    ! module for file cubic_eos.f90
    module cubic_eos_module
    !global use inclusion
    use module_all_types

    
! interface for routines with circular dependencies (implemented in submodule in file <this-file-name>_impl)    
    interface
    
    
    module subroutine init_SRK(gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: temp
    integer :: nrsubst
    end subroutine init_SRK
    
    module subroutine da_dxi_dtau_cubic_all (gl,Temp, dens, get_a_deriv, dnadtaun, dnadxidtaum, dnadxi2dtaum)
    implicit none
    type(type_gl) :: gl
    double precision:: TEMP, DENS
    integer get_a_deriv                             ! array specifier to indicate, which derivatives are needed 
    double precision, dimension(15)::dnadtaun       ! array with the computed values for the derivatives
    double precision, dimension(15,30)::dnadxidtaum       ! array with the computed values for the derivatives
    double precision, dimension(15,30,30) :: dnadxi2dtaum       ! array with the computed values for the derivatives
    end subroutine da_dxi_dtau_cubic_all
    


    double precision module Function SRK(gl,T,P0_lc,iphase, nrsubst) 
    implicit none
    type(type_gl) :: gl
    real *8 :: t,p,P0_lc ! input T and P
    integer :: iphase !input iphase: 1-liquid; 2-vapor
    integer :: nrsubst
    End Function


    double precision module function rho_SRK (gl,Temp, press_in, Iphase, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: Temp, press, press_in
    integer:: Iphase, nrsubst
    end function rho_SRK

    
    module subroutine da_SRK_dxi(gl,daSRKdxi)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30):: daSRKdxi
    end subroutine da_SRK_dxi

    
    double precision module function da_SRK_dtau(gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: temp
    integer:: nrsubst
    end function da_SRK_dtau

    
    double precision module function d2a_SRK_dtau2(gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision :: Temp
    integer:: nrsubst
    end function d2a_SRK_dtau2


    double precision module function d3a_SRK_dtau3(gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision :: Temp
    integer:: nrsubst
    end function d3a_SRK_dtau3


    module subroutine d2a_SRK_dxidtau(gl,Temp,d2aSRKdxidtau)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30):: d2aSRKdxidtau
    double precision:: temp
    end subroutine d2a_SRK_dxidtau

    
    module subroutine d2a_SRK_dxidxj(gl,d2aSRKdxidxj)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30, 30):: d2aSRKdxidxj
    end subroutine d2a_SRK_dxidxj


    module subroutine db_SRK_dxi(gl,dbSRKdxi)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30):: dbSRKdxi
    end subroutine db_SRK_dxi

    
    
    module subroutine d2b_SRK_dxidxj(gl,d2bSRKdxidxj)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30, 30):: d2bSRKdxidxj
    end subroutine d2b_SRK_dxidxj
    

    Double Precision module Function PSRK(gl,T,P0_lc,iphase, nrsubst) 
    implicit none
    type(type_gl) :: gl
    integer :: iphase                                ! input iphase: 1-liquid; 2-vapor
    integer:: nrsubst                             ! if 0, then mixture, otherwise pure fluid
    double precision :: T,P0_lc,p                         ! input: T,p0  transfered to p
    End Function


    module subroutine CUBIC_NT(gl,para,volume)
    implicit none
    type(type_gl) :: gl
    double precision :: para(3)                ! a,b,c input parameters
    double precision :: volume(3)              ! output [m³/mol]
    end  

    
    module subroutine init_PR(gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: temp
    integer :: nrsubst
    end subroutine init_PR
    
    
    double precision module function rho_PR (gl,Temp, press_in, Iphase, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: Temp, press_in
    integer:: Iphase, nrsubst
    end function rho_PR

    
    double precision module function da_PR_dtau (gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: temp
    integer:: nrsubst
    end function da_PR_dtau 

    
    double precision module function d2a_PR_dtau2(gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision :: Temp
    integer:: nrsubst
    end function d2a_PR_dtau2

    
    double precision module function d3a_PR_dtau3(gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision :: Temp
    integer:: nrsubst
    end function d3a_PR_dtau3



    module subroutine da_PR_dxi(gl,daPRdxi)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30):: daPRdxi
    end subroutine da_PR_dxi

    
    module subroutine d2a_PR_dxidxj(gl,d2aPRdxidxj)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30, 30):: d2aPRdxidxj
    end subroutine d2a_PR_dxidxj

    
    module subroutine d2a_PR_dxidtau(gl,Temp,d2aPRdxidtau)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30):: d2aPRdxidtau
    double precision:: temp
    end subroutine d2a_PR_dxidtau


    module subroutine db_PR_dxi(gl,dbPRdxi)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30):: dbPRdxi
    end subroutine db_PR_dxi


    module subroutine d2b_PR_dxidxj(gl,d2bPRdxidxj)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(30, 30):: d2bPRdxidxj
    end subroutine d2b_PR_dxidxj

    
    module subroutine dar_dxi_PR(gl,Temp, dens, dardxi_PR)
    implicit none
    type(type_gl) :: gl
    double precision:: TEMP, DENS 
    double precision, dimension(30):: dardxi_PR
    end subroutine dar_dxi_PR

    
    module subroutine d2ar_dxidxj_PR(gl,Temp, dens, d2ardxidxj_PR)
    implicit none
    type(type_gl) :: gl
    double precision:: TEMP, DENS 
    double precision, dimension(30,30) :: d2ardxidxj_PR
    end subroutine d2ar_dxidxj_PR

    
    
    module subroutine d2ar_dxidtau_PR(gl,Temp, dens, d2ardxidtau_PR)
    implicit none
    type(type_gl) :: gl
    double precision:: TEMP, DENS 
    double precision, dimension(30)::d2ardxidtau_PR
    end subroutine d2ar_dxidtau_PR

    
    module subroutine d2ar_dxiddel_PR(gl,Temp, dens, d2ardxiddel_PR)
    implicit none
    type(type_gl) :: gl
    double precision:: TEMP, DENS 
    double precision, dimension(30)::d2ardxiddel_PR
    end subroutine d2ar_dxiddel_PR

    
    double precision module function d4a_SRK_dtau4(gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    Double precision:: temp
    integer :: nrsubst
    end function d4a_SRK_dtau4


    double precision module function d4a_PR_dtau4(gl,Temp, nrsubst)
    implicit none
    type(type_gl) :: gl
    Double precision:: temp
    integer :: nrsubst
    end function d4a_PR_dtau4


    module subroutine MIXDERIVSFNR_HIGHER_CUBIC (gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)
    implicit none
    type(type_gl) :: gl
    double precision:: TEMPERATURE, DENSITY  
    integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed 
    double precision, dimension(15)::MIXDERFNR  ! array with the computed values for the derivatives
    End subroutine MIXDERIVSFNR_HIGHER_CUBIC


    module subroutine MIXDERIVSFNR_dxi_CUBIC (gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR_dxi)
    implicit none
    type(type_gl) :: gl
    double precision:: TEMPERATURE, DENSITY  
    integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed 
    double precision, dimension(15,30)::MIXDERFNR_dxi  ! array with the computed values for the derivatives
    End subroutine MIXDERIVSFNR_dxi_CUBIC


    module subroutine MIXDERIVSFNR_dxidxj_CUBIC (gl,TEMPERATURE, DENSITY, GETDER)
    implicit none
    type(type_gl) :: gl
    double precision:: TEMPERATURE, DENSITY  
    integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed 
    End subroutine MIXDERIVSFNR_dxidxj_CUBIC


    module subroutine MIXDERIVSFNR_dxidxjdxk_CUBIC (gl,TEMPERATURE, DENSITY, GETDER)
    implicit none
    type(type_gl) :: gl
    double precision:: TEMPERATURE, DENSITY  
    integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed 
    End subroutine MIXDERIVSFNR_dxidxjdxk_CUBIC

    
    end interface    
    
    
    
    contains

    end module cubic_eos_module
