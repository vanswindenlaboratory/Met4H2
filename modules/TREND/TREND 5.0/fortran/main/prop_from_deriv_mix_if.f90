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



    MODULE calc_functions
    use module_all_types
    use module_regula_falsi
    use module_regula_falsi_support

    ! interface for routines with circular dependencies (implemented in submodule in file <this-file-name>_impl)
    interface

    DOUBLE PRECISION module FUNCTION OMEGA_CALC(gl,nrsubst)
    implicit none
    type(type_gl) :: gl
    INTEGER::nrsubst
    end FUNCTION OMEGA_CALC

    DOUBLE PRECISION module FUNCTION HSAT_CALC(gl,TSAT, IPHASE, nrsubst)
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION:: TSAT
    integer:: IPHASE, nrsubst
    end FUNCTION HSAT_CALC

    Double precision module function DSPIN_CALC (gl,T,spinprop,nrsubst) !mol/m³
    implicit none
    type(type_gl) :: gl
    double precision:: T
    integer:: nrsubst, spinprop
    end FUNCTION DSPIN_CALC


    DOUBLE PRECISION module FUNCTION P_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION P_CALC


    DOUBLE PRECISION module FUNCTION P_CALC_FROM_DER(gl,T,D, nrsubst, FNRDER)
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION P_CALC_FROM_DER

    DOUBLE PRECISION module FUNCTION Z_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION Z_CALC

    DOUBLE PRECISION module FUNCTION DPDD_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION DPDD_CALC

    DOUBLE PRECISION module FUNCTION DPDD_CALC_FROM_DER(gl,T,D, nrsubst, FNRDER)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(nderivs):: FNRDER
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION DPDD_CALC_FROM_DER

    double precision module function D2PDD2_CALC (gl,T, D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION D2PDD2_CALC

    double precision module function D2PDD2_CALC_FROM_DER (gl,T, D, nrsubst, der_res)
    implicit none
    type(type_gl) :: gl
    double precision, dimension(nderivs):: der_res
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION D2PDD2_CALC_FROM_DER

    DOUBLE PRECISION module FUNCTION DPDT_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION DPDT_CALC

    DOUBLE PRECISION module FUNCTION DDDT_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION DDDT_CALC

    DOUBLE PRECISION module FUNCTION D2PDTT_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION D2PDTT_CALC

    DOUBLE PRECISION module FUNCTION D2PDDT_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION D2PDDT_CALC

    double precision module function FUGCOPURE_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION FUGCOPURE_CALC

    double precision module function CPOTR_CALC(gl,T,D,nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION CPOTR_CALC

    double precision module function CPOTI_CALC(gl,T,D,nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION CPOTI_CALC

    double precision module function DFugcoefpureDD_CALC (gl,T, D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION DFugcoefpureDD_CALC

    DOUBLE PRECISION module FUNCTION S_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION S_CALC

    DOUBLE PRECISION module FUNCTION U_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION U_CALC

    DOUBLE PRECISION module FUNCTION CV_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION CV_CALC

    DOUBLE PRECISION module FUNCTION H_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION H_CALC

    DOUBLE PRECISION module FUNCTION CP_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION CP_CALC

    DOUBLE PRECISION module FUNCTION PIP_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION PIP_CALC

    DOUBLE PRECISION module FUNCTION GRUEN_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION GRUEN_CALC

    DOUBLE PRECISION module FUNCTION CP0_CALC(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION CP0_CALC

    DOUBLE PRECISION module FUNCTION G_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION G_CALC

    DOUBLE PRECISION module FUNCTION GR_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION GR_CALC

    DOUBLE PRECISION module FUNCTION WS_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION WS_CALC

    DOUBLE PRECISION module FUNCTION JTCO_CALC(gl,T,D,nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION JTCO_CALC

    DOUBLE PRECISION module FUNCTION B_CALC(gl,T, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION B_CALC

    DOUBLE PRECISION module FUNCTION DBDT_CALC(gl,T, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION DBDT_CALC

    DOUBLE PRECISION module FUNCTION C_CALC(gl,T, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION C_CALC

    DOUBLE PRECISION module FUNCTION DCDT_CALC(gl,T, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION DCDT_CALC

    DOUBLE PRECISION module FUNCTION D2CD2T_CALC(gl,T, nrsubst)

    implicit none

    type(type_gl) :: gl
    integer:: nrsubst                                 ! Flag = 0: mixture property; Flag > 0:  Fluid of the same number will be calculated
    DOUBLE PRECISION :: T
    end function D2CD2T_CALC

    DOUBLE PRECISION module FUNCTION D_CALC(gl,T, nrsubst)
    implicit none

    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION D_CALC

    DOUBLE PRECISION module FUNCTION VOLEXP_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION VOLEXP_CALC

    DOUBLE PRECISION module FUNCTION COMPT_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION COMPT_CALC

    DOUBLE PRECISION module FUNCTION COMPS_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION COMPS_CALC

    DOUBLE PRECISION module FUNCTION THROT_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION THROT_CALC

    DOUBLE PRECISION module FUNCTION EXPANS_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION EXPANS_CALC

    DOUBLE PRECISION module FUNCTION EXPANT_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION EXPANT_CALC

    DOUBLE PRECISION module FUNCTION DUDV_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION DUDV_CALC

    DOUBLE PRECISION module FUNCTION A_CALC(gl,T, D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A_CALC

    DOUBLE PRECISION module FUNCTION AI_CALC(gl,T, D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION AI_CALC

    DOUBLE PRECISION module FUNCTION BETA_CALC(gl,T, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION BETA_CALC

    DOUBLE PRECISION module FUNCTION UR_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION UR_CALC

    DOUBLE PRECISION module FUNCTION SR_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION SR_CALC

    DOUBLE PRECISION module FUNCTION CVR_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION CVR_CALC

    DOUBLE PRECISION module FUNCTION HR_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION HR_CALC

    DOUBLE PRECISION module FUNCTION AR_CALC(gl,T, D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION AR_CALC

    DOUBLE PRECISION module FUNCTION PR_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION PR_CALC

    DOUBLE PRECISION module FUNCTION A00_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A00_CALC

    DOUBLE PRECISION module FUNCTION A10_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A10_CALC

    DOUBLE PRECISION module FUNCTION A01_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A01_CALC

    DOUBLE PRECISION module FUNCTION A11_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A11_CALC

    DOUBLE PRECISION module FUNCTION A20_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A20_CALC

    DOUBLE PRECISION module FUNCTION A02_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A02_CALC

    DOUBLE PRECISION module FUNCTION A12_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A12_CALC

    DOUBLE PRECISION module FUNCTION A03_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A03_CALC

    DOUBLE PRECISION module FUNCTION A30_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A30_CALC

    DOUBLE PRECISION module FUNCTION A40_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A40_CALC

    DOUBLE PRECISION module FUNCTION A21_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION A21_CALC


    double precision module function RIEM_CALC (gl,t,d,nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION RIEM_CALC

    DOUBLE PRECISION module FUNCTION GAMMAGD_CALC(gl,T, D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION GAMMAGD_CALC

    DOUBLE PRECISION module FUNCTION DCVDT_CALC(gl,T, D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION DCVDT_CALC

    DOUBLE PRECISION module FUNCTION DCVRDT_CALC(gl,T, D, nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION DCVRDT_CALC

    DOUBLE PRECISION module FUNCTION D4PD4_CALC(gl,T,D,nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION D4PD4_CALC



    Double Precision module Function TBoyle_calc(gl,fluidnr)
    implicit none
    type(type_gl) :: gl
    INTEGER::fluidnr
    end FUNCTION TBoyle_calc

    Double Precision module Function TJT_calc(gl,fluidnr)
    implicit none
    type(type_gl) :: gl
    INTEGER::fluidnr
    end FUNCTION TJT_calc

    Double Precision module Function TJTINV_calc(gl,fluidnr)
    implicit none
    type(type_gl) :: gl
    INTEGER::fluidnr
    end FUNCTION TJTINV_calc

    DOUBLE PRECISION module FUNCTION H_REF(gl,T,D,nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION H_REF

    DOUBLE PRECISION module FUNCTION S_REF(gl,T,D,nrsubst)
    implicit none
    type(type_gl) :: gl
    double precision:: T,D
    INTEGER::nrsubst
    end FUNCTION S_REF

    module Subroutine FUGCO_CALC_MIX(gl,T,D, fugcoef_mix)
    !DEC$ ATTRIBUTES DLLEXPORT :: FUGCO_CALC_MIX
    implicit none
    type(type_gl) :: gl
    double precision :: T, D
    double precision, dimension(30), intent(out):: fugcoef_mix
    end Subroutine FUGCO_CALC_MIX

    DOUBLE PRECISION module FUNCTION DlnFugcoef_pureDD (gl,T, D, i)
    implicit none
    type(type_gl) :: gl
    double precision :: T, D
    integer :: i
    end FUNCTION DlnFugcoef_pureDD


    DOUBLE PRECISION module FUNCTION DFugcoef_pureDT (gl,T, D, i)
    implicit none
    type(type_gl) :: gl
    double precision :: T, D
    integer :: i
    end function DFugcoef_pureDT

    module Subroutine FUGCO_CALC_PURE_TP(gl,T,P, fugcoef_pure)
    implicit none
    type(type_gl) :: gl
    double precision :: T, P
    double precision, dimension(30), intent(out)::fugcoef_pure
    END Subroutine FUGCO_CALC_PURE_TP

    module Subroutine Chempot_CALC(gl,T,D, Chempot, oir)
    implicit none
    type(type_gl) :: gl
    double precision:: T, D
    integer :: oir
    double precision, dimension(30), intent(out):: chempot
    end subroutine Chempot_CALC

    DOUBLE PRECISION module FUNCTION S_CALC_TP(gl,T,P)
    implicit none
    type(type_gl) :: gl
    double precision :: T,P
    END FUNCTION S_CALC_TP
    DOUBLE PRECISION module FUNCTION H_CALC_TP(gl,T,P)
    implicit none
    type(type_gl) :: gl
    double precision :: T,P
    END FUNCTION H_CALC_TP

    DOUBLE PRECISION module FUNCTION Res_TBoyle(gl,T_akt,parameters)
    implicit none
    type(type_gl) :: gl
    type(type_additional_parameters) :: parameters
    DOUBLE PRECISION :: T_akt
    END function Res_TBoyle

    DOUBLE PRECISION module FUNCTION Res_TJT(gl,T_akt,parameters)
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION :: T_akt
    type(type_additional_parameters) :: parameters
    END function Res_TJT

    DOUBLE PRECISION module FUNCTION Res_TJTINV(gl,T_akt,parameters)
    implicit none
    type(type_gl) :: gl
    type(type_additional_parameters) :: parameters
    DOUBLE PRECISION :: T_akt
    END function Res_TJTINV

    DOUBLE PRECISION module FUNCTION TDENSMAX_calc(gl,press,nrsubst)
    implicit none
    type(type_gl) :: gl
    INTEGER :: nrsubst
    DOUBLE PRECISION :: press
    End Function TDENSMAX_calc

    DOUBLE PRECISION module FUNCTION Res_TDENMAX(gl,T_akt,parameters)
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION :: T_akt
    type(type_additional_parameters) :: parameters
    END function Res_TDENMAX

    DOUBLE PRECISION module FUNCTION Res_Dspin(gl,D_akt,parameters)
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION :: D_akt
    type(type_additional_parameters) :: parameters
    END function Res_Dspin

    module Subroutine GAMMATF_CALC(gl,T, D, GAMMATF)
    implicit none
    type(type_gl) :: gl
    double precision:: T, D
    double precision, dimension(30, 30), intent(out):: GAMMATF
    END SUBROUTINE GAMMATF_CALC




    DOUBLE PRECISION module FUNCTION D_PIP_DD_CALC(gl,T,D, nrsubst)

    use module_all_types
    implicit none

    type(type_gl) :: gl
    DOUBLE PRECISION :: T, D
    INTEGER :: nrsubst
    END FUNCTION D_PIP_DD_CALC




    DOUBLE PRECISION module FUNCTION D2_PIP_D2D_CALC(gl,T,D, nrsubst)

    use module_all_types
    implicit none

    type(type_gl) :: gl
    DOUBLE PRECISION :: T, D
    INTEGER :: nrsubst
    END FUNCTION D2_PIP_D2D_CALC

    DOUBLE PRECISION module FUNCTION D_PIP_DT_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION :: T, D
    INTEGER :: nrsubst
    END FUNCTION D_PIP_DT_CALC

    DOUBLE PRECISION module FUNCTION d2DdT2_CALC(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    END FUNCTION d2DdT2_CALC

    DOUBLE PRECISION module FUNCTION d2BdT2_CALC(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    END FUNCTION d2BdT2_CALC

    DOUBLE PRECISION module FUNCTION d2CdT2_CALC(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    END FUNCTION d2CdT2_CALC


    DOUBLE PRECISION module FUNCTION d3BdT3_CALC(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    END FUNCTION d3BdT3_CALC


    DOUBLE PRECISION module FUNCTION d3CdT3_CALC(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    END FUNCTION d3CdT3_CALC


    DOUBLE PRECISION module FUNCTION d3DdT3_CALC(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    END FUNCTION d3DdT3_CALC


    DOUBLE PRECISION module FUNCTION d4BdT4_CALC(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    END FUNCTION d4BdT4_CALC


    DOUBLE PRECISION module FUNCTION d4CdT4_CALC(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    END FUNCTION d4CdT4_CALC


    DOUBLE PRECISION module FUNCTION d4DdT4_CALC(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    END FUNCTION d4DdT4_CALC

    double precision module function b_calc_num(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    end function b_calc_num
    double precision module function c_calc_num(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    end function c_calc_num
    double precision module function d_calc_num(gl,T,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    double precision:: T
    end function d_calc_num
    double precision module function n_eff_calc(gl,T,D,nrsubst)
    implicit none
    type(type_gl):: gl
    integer:: nrsubst
    DOUBLE PRECISION :: T,D
    end function
    double precision module function zeta_calc(gl,T,D,nrsubst)
    implicit none
    type(type_gl):: gl
    integer:: nrsubst
    DOUBLE PRECISION :: T,D
    end function
    double precision module function d2zetadd2_calc(gl,T,D,nrsubst)
    implicit none
    type(type_gl):: gl
    integer:: nrsubst
    DOUBLE PRECISION :: T,D
    end function

    double precision module function GRUEN_WO_CV0_CALC(gl,T,D, nrsubst)
    implicit none
    type(type_gl) :: gl
    DOUBLE PRECISION :: T, D
    INTEGER :: nrsubst
    end function

    DOUBLE PRECISION module FUNCTION GRUEN_WO_CV0_CALC_DT(gl,T,D, nrsubst)
    type(type_gl) :: gl
    DOUBLE PRECISION :: T, D
    INTEGER :: nrsubst
    end function

    end interface


    contains
    !*************************************************************************************
    !				TREND Version 3.0
    !		   Thermodynamic Reference & Engineering Data
    !
    !- software for the calculation of thermodynamic and other properties -
    !
    !Copyright (C) 2016,  Prof. Dr.-Ing. R.Span
    !                     Lehrstuhl fuer Thermodynamik
    !                     Ruhr-Universitaet Bochum
    !                     Universitaetsstr. 150
    !                     D-44892 Bochum
    !
    !Cite as: Span, R.; Eckermann, T.; Herrig, S.; Hielscher, S.; Jäger, A.; Thol, M. (2016):
    !	     TREND. Thermodynamic Reference and Engineering Data 3.0. Lehrstuhl fuer
    !	     Thermodynamik, Ruhr-Universitaet Bochum.
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





    end module calc_functions



    module calc_func_ptr

    use calc_functions
    use module_all_types
    !abstract interface for all the functions
    abstract interface

    double precision function custom_name(gl,T,D,nrsubst)

    use module_all_types

    implicit none

    !--------------------------------------------------------------------
    double precision :: T
    double precision  :: D   !input temperature and density
    integer :: nrsubst
    type(type_gl) :: gl
    !--------------------------------------------------------------------

    end function custom_name

    double precision function custom_name2(gl,T,nrsubst)

    use module_all_types

    implicit none

    !--------------------------------------------------------------------
    double precision  :: T
    integer :: nrsubst
    type(type_gl) :: gl
    !--------------------------------------------------------------------

    end function custom_name2

    double precision function custom_name3(gl,T,int_prop,nrsubst)

    use module_all_types

    implicit none

    !--------------------------------------------------------------------
    double precision  :: T
    integer :: nrsubst,int_prop
    type(type_gl) :: gl
    !--------------------------------------------------------------------

    end function custom_name3

    double precision function custom_name4(gl, fluidnr)
    use module_all_types
    implicit none
    type(type_gl) :: gl
    INTEGER*4 :: errTBoyle,fluidnr
    end function custom_name4


    end interface

    !container type to create the pointers of the procedures(functions)
    type container
        procedure(custom_name), pointer, nopass :: ptr !(gl, T, D, nrsubst)
        procedure(custom_name2), pointer, nopass :: ptr2 !(gl, T, nrsubst)
        procedure(custom_name3), pointer, nopass :: ptr3 !(gl, T, integer, nrsubst)
        procedure(custom_name4), pointer, nopass :: ptr4 !(gl, fluidnr)
        character(20):: propname
    end type
    !
    type func_ptrs
        type(container), dimension(100):: function_ptrs
    end type


    type(func_ptrs), pointer:: act_func_ptrs

    !procedure( P_CALC ) , pointer:: ptr_P_CALC
    !procedure( Z_CALC ) , pointer:: ptr_Z_CALC
    !procedure( DPDD_CALC ) , pointer:: ptr_DPDD_CALC
    !procedure( D2PDD2_CALC ) , pointer:: ptr_D2PDD2_CALC
    !procedure( DPDT_CALC ) , pointer:: ptr_DPDT_CALC
    !procedure( DDDT_CALC ) , pointer:: ptr_DDDT_CALC
    !procedure( D2PDTT_CALC ) , pointer:: ptr_D2PDTT_CALC
    !procedure( D2PDDT_CALC ) , pointer:: ptr_D2PDDT_CALC
    !!procedure( FUGCO_CALC_MIX ) , pointer:: ptr_FUGCO_CALC_MIX
    !procedure( FUGCOPURE_CALC ) , pointer:: ptr_FUGCOPURE_CALC
    !procedure( CPOTR_CALC ) , pointer:: ptr_CPOTR_CALC
    !procedure( CPOTI_CALC ) , pointer:: ptr_CPOTI_CALC
    !procedure( DlnFugcoef_pureDD ) , pointer:: ptr_DlnFugcoef_pureDD
    !procedure( DFugcoefpureDD_CALC ) , pointer:: ptr_DFugcoefpureDD_CALC
    !procedure( DFugcoef_pureDT ) , pointer:: ptr_DFugcoef_pureDT
    !procedure( FUGCO_CALC_PURE_TP ) , pointer:: ptr_FUGCO_CALC_PURE_TP
    !procedure( S_CALC ) , pointer:: ptr_S_CALC
    !procedure( U_CALC ) , pointer:: ptr_U_CALC
    !procedure( CV_CALC ) , pointer:: ptr_CV_CALC
    !procedure( H_CALC ) , pointer:: ptr_H_CALC
    !procedure( CP_CALC ) , pointer:: ptr_CP_CALC
    !procedure( PIP_CALC ) , pointer:: ptr_PIP_CALC
    !procedure( GRUEN_CALC ) , pointer:: ptr_GRUEN_CALC
    !procedure( CP0_CALC ) , pointer:: ptr_CP0_CALC
    !procedure( OMEGA_CALC ) , pointer:: ptr_OMEGA_CALC
    !procedure( G_CALC ) , pointer:: ptr_G_CALC
    !procedure( GR_CALC ) , pointer:: ptr_GR_CALC
    !procedure( Chempot_CALC ) , pointer:: ptr_Chempot_CALC
    !procedure( WS_CALC ) , pointer:: ptr_WS_CALC
    !procedure( JTCO_CALC ) , pointer:: ptr_JTCO_CALC
    !procedure( B_CALC ) , pointer:: ptr_B_CALC
    !procedure( DBDT_CALC ) , pointer:: ptr_DBDT_CALC
    !procedure( C_CALC ) , pointer:: ptr_C_CALC
    !procedure( DCDT_CALC ) , pointer:: ptr_DCDT_CALC
    !procedure( D_CALC ) , pointer:: ptr_D_CALC
    !procedure( VOLEXP_CALC ) , pointer:: ptr_VOLEXP_CALC
    !procedure( COMPT_CALC ) , pointer:: ptr_COMPT_CALC
    !procedure( COMPS_CALC ) , pointer:: ptr_COMPS_CALC
    !procedure( THROT_CALC ) , pointer:: ptr_THROT_CALC
    !procedure( EXPANS_CALC ) , pointer:: ptr_EXPANS_CALC
    !procedure( EXPANT_CALC ) , pointer:: ptr_EXPANT_CALC
    !procedure( DUDV_CALC ) , pointer:: ptr_DUDV_CALC
    !procedure( A_CALC ) , pointer:: ptr_A_CALC
    !procedure( AI_CALC ) , pointer:: ptr_AI_CALC
    !procedure( BETA_CALC ) , pointer:: ptr_BETA_CALC
    !procedure( HSAT_CALC ) , pointer:: ptr_HSAT_CALC
    !procedure( H_REF ) , pointer:: ptr_H_REF
    !procedure( S_REF ) , pointer:: ptr_S_REF
    !procedure( UR_CALC ) , pointer:: ptr_UR_CALC
    !procedure( SR_CALC ) , pointer:: ptr_SR_CALC
    !procedure( CVR_CALC ) , pointer:: ptr_CVR_CALC
    !procedure( HR_CALC ) , pointer:: ptr_HR_CALC
    !procedure( AR_CALC ) , pointer:: ptr_AR_CALC
    !procedure( PR_CALC ) , pointer:: ptr_PR_CALC
    !procedure( A00_CALC ) , pointer:: ptr_A00_CALC
    !procedure( A10_CALC ) , pointer:: ptr_A10_CALC
    !procedure( A01_CALC ) , pointer:: ptr_A01_CALC
    !procedure( A11_CALC ) , pointer:: ptr_A11_CALC
    !procedure( A20_CALC ) , pointer:: ptr_A20_CALC
    !procedure( A02_CALC ) , pointer:: ptr_A02_CALC
    !procedure( A12_CALC ) , pointer:: ptr_A12_CALC
    !procedure( A03_CALC ) , pointer:: ptr_A03_CALC
    !procedure( A30_CALC ) , pointer:: ptr_A30_CALC
    !procedure( A21_CALC ) , pointer:: ptr_A21_CALC
    !procedure( S_CALC_TP ) , pointer:: ptr_S_CALC_TP
    !procedure( H_CALC_TP ) , pointer:: ptr_H_CALC_TP
    !procedure( RIEM_CALC ) , pointer:: ptr_RIEM_CALC
    !procedure( TBoyle_calc ) , pointer:: ptr_TBoyle_calc
    !procedure( Res_TBoyle ) , pointer:: ptr_Res_TBoyle
    !procedure( TJT_calc ) , pointer:: ptr_TJT_calc
    !procedure( Res_TJT ) , pointer:: ptr_Res_TJT
    !procedure( TJTINV_calc ) , pointer:: ptr_TJTINV_calc
    !procedure( Res_TJTINV ) , pointer:: ptr_Res_TJTINV
    !procedure( TDENMAX_calc ) , pointer:: ptr_TDENMAX_calc
    !procedure( Res_TDENMAX ) , pointer:: ptr_Res_TDENMAX
    !procedure( DSPIN_CALC ) , pointer:: ptr_DSPIN_CALC
    !procedure( Res_Dspin ) , pointer:: ptr_Res_Dspin
    !procedure( GAMMAGD_CALC ) , pointer:: ptr_GAMMAGD_CALC
    !procedure( GAMMATF_CALC ) , pointer:: ptr_GAMMATF_CALC
    !procedure( DCVDT_CALC ) , pointer:: ptr_DCVDT_CALC
    !procedure( DCVRDT_CALC ) , pointer:: ptr_DCVRDT_CALC

    contains

    !in this subroutine the pointers for the necessary functions are set
    subroutine set_func_pointers()
    use, intrinsic :: iso_c_binding

    implicit none
    !   prop    property
    !    1       T
    !type(c_ptr), intent(inout):: handle_procedures
    !Pointer to all functions

    allocate(act_func_ptrs)


    act_func_ptrs%function_ptrs(1)%propname = "T"
    act_func_ptrs%function_ptrs(2)%propname = "D"
    act_func_ptrs%function_ptrs(3)%ptr => P_CALC                     !    3       P
    act_func_ptrs%function_ptrs(3)%propname = "P"
    act_func_ptrs%function_ptrs(4)%ptr => U_CALC                     !    4       U
    act_func_ptrs%function_ptrs(4)%propname = "U"
    act_func_ptrs%function_ptrs(5)%ptr => H_CALC                     !    5       H
    act_func_ptrs%function_ptrs(5)%propname = "H"
    act_func_ptrs%function_ptrs(6)%ptr => S_CALC                     !    6       S
    act_func_ptrs%function_ptrs(6)%propname = "S"
    act_func_ptrs%function_ptrs(7)%ptr => G_CALC                     !    7       G
    act_func_ptrs%function_ptrs(7)%propname = "G"
    act_func_ptrs%function_ptrs(8)%ptr => A_CALC                     !    8       A
    act_func_ptrs%function_ptrs(8)%propname = "A"
    act_func_ptrs%function_ptrs(9)%ptr => CP_CALC                    !    9       CP
    act_func_ptrs%function_ptrs(9)%propname = "CP"
    act_func_ptrs%function_ptrs(10)%ptr => CV_CALC                   !   10       CV
    act_func_ptrs%function_ptrs(10)%propname = "CV"
    act_func_ptrs%function_ptrs(11)%ptr => WS_CALC                   !   11       WS
    act_func_ptrs%function_ptrs(11)%propname = "WS"
    act_func_ptrs%function_ptrs(12)%ptr2 => B_CALC                    !   12       B
    act_func_ptrs%function_ptrs(12)%propname = "BVIR"
    act_func_ptrs%function_ptrs(13)%ptr2 => C_CALC                    !   13       C
    act_func_ptrs%function_ptrs(13)%propname = "CVIR"
    act_func_ptrs%function_ptrs(14)%ptr2 => CP0_CALC                  !   14       CP0
    act_func_ptrs%function_ptrs(14)%propname = "CP0"
    act_func_ptrs%function_ptrs(15)%propname = "Q"                   !   15       Q
    act_func_ptrs%function_ptrs(16)%ptr => Z_CALC                    !   16       Z
    act_func_ptrs%function_ptrs(16)%propname = "Z"
    act_func_ptrs%function_ptrs(17)%propname = "PNUMER"             !   17       PNUMER
    act_func_ptrs%function_ptrs(18)%propname = "WSNUMER"            !   18       WSNUMER
    act_func_ptrs%function_ptrs(19)%ptr => A00_CALC                   !   19       A00
    act_func_ptrs%function_ptrs(19)%propname = "A00"
    act_func_ptrs%function_ptrs(20)%ptr => A01_CALC                   !   20       A01
    act_func_ptrs%function_ptrs(20)%propname = "A01"
    act_func_ptrs%function_ptrs(21)%ptr => A02_CALC                   !   21       A02
    act_func_ptrs%function_ptrs(21)%propname = "A02"
    act_func_ptrs%function_ptrs(22)%ptr => A03_CALC                   !   22       A03
    act_func_ptrs%function_ptrs(22)%propname = "A03"
    act_func_ptrs%function_ptrs(23)%ptr => A10_CALC                   !   23       A10
    act_func_ptrs%function_ptrs(23)%propname = "A10"
    act_func_ptrs%function_ptrs(24)%ptr => A11_CALC                   !   24       A11
    act_func_ptrs%function_ptrs(24)%propname = "A11"
    act_func_ptrs%function_ptrs(25)%ptr => A12_CALC                   !   25       A12
    act_func_ptrs%function_ptrs(25)%propname = "A12"
    act_func_ptrs%function_ptrs(26)%ptr => A20_CALC                   !   26       A20
    act_func_ptrs%function_ptrs(26)%propname = "A20"
    act_func_ptrs%function_ptrs(27)%ptr => A21_CALC                   !   27       A21
    act_func_ptrs%function_ptrs(27)%propname = "A21"
    act_func_ptrs%function_ptrs(28)%ptr => A30_CALC                   !   28       A30
    act_func_ptrs%function_ptrs(28)%propname = "A30"
    act_func_ptrs%function_ptrs(29)%ptr => UR_CALC                    !   29       UR
    act_func_ptrs%function_ptrs(29)%propname = "UR"
    act_func_ptrs%function_ptrs(30)%ptr => HR_CALC                    !   30       HR
    act_func_ptrs%function_ptrs(30)%propname = "HR"
    act_func_ptrs%function_ptrs(31)%ptr => SR_CALC                    !   31       SR
    act_func_ptrs%function_ptrs(31)%propname = "SR"
    act_func_ptrs%function_ptrs(32)%ptr => CVR_CALC                   !   32       CVR
    act_func_ptrs%function_ptrs(32)%propname = "CVR"
    act_func_ptrs%function_ptrs(33)%ptr => CPOTR_CALC                 !   33       CPOTR
    act_func_ptrs%function_ptrs(33)%propname = "CPOTR"
    act_func_ptrs%function_ptrs(34)%ptr => GRUEN_CALC                 !   34       GRUEN
    act_func_ptrs%function_ptrs(34)%propname = "GRUEN"
    act_func_ptrs%function_ptrs(35)%ptr => PIP_CALC                   !   35       PIP
    act_func_ptrs%function_ptrs(35)%propname = "PIP"
    act_func_ptrs%function_ptrs(36)%propname = "DE"                     !   36       DE
    act_func_ptrs%function_ptrs(37)%propname = "ST"                     !   37       ST
    act_func_ptrs%function_ptrs(38)%propname = "VISDYN"                 !   38       VISDYN
    act_func_ptrs%function_ptrs(39)%propname = "TCX"                   !   39       TCX
    act_func_ptrs%function_ptrs(40)%ptr => VOLEXP_CALC                 !   40       VOLEXP
    act_func_ptrs%function_ptrs(40)%propname = "VOLEXP"
    act_func_ptrs%function_ptrs(41)%ptr2 => D_CALC                     !   41       DVIR
    act_func_ptrs%function_ptrs(41)%propname = "DVIR"
    act_func_ptrs%function_ptrs(42)%ptr => RIEM_CALC                   !   42       RIEM
    act_func_ptrs%function_ptrs(42)%propname = "RIEM"
    act_func_ptrs%function_ptrs(43)%ptr => COMPT_CALC                  !   43       COMPT
    act_func_ptrs%function_ptrs(43)%propname = "COMPT"
    act_func_ptrs%function_ptrs(44)%ptr => COMPS_CALC                  !   44       COMPS
    act_func_ptrs%function_ptrs(44)%propname = "COMPS"
    act_func_ptrs%function_ptrs(45)%propname = "throt"                  !  45
    act_func_ptrs%function_ptrs(46)%ptr => EXPANS_CALC                 !   46       EXPANS
    act_func_ptrs%function_ptrs(46)%propname = "EXPANS"
    act_func_ptrs%function_ptrs(47)%ptr => EXPANT_CALC                 !   47       EXPANT
    act_func_ptrs%function_ptrs(47)%propname = "EXPANT"
    act_func_ptrs%function_ptrs(48)%ptr => DUDV_CALC                   !   48       DUDV
    act_func_ptrs%function_ptrs(48)%propname = "DUDV"
    act_func_ptrs%function_ptrs(49)%ptr => JTCO_CALC                   !   49       JTCO
    act_func_ptrs%function_ptrs(49)%propname = "JTCO"
    act_func_ptrs%function_ptrs(50)%ptr => DPDT_CALC                   !   50       DPDT
    act_func_ptrs%function_ptrs(50)%propname = "DPDT"
    act_func_ptrs%function_ptrs(51)%propname = "VE"                     !   51       VE
    act_func_ptrs%function_ptrs(52)%propname = "HE"                      !   52       HE
    act_func_ptrs%function_ptrs(53)%propname = "GE"                        !   53       GE
    act_func_ptrs%function_ptrs(54)%propname = "B12"                        !  54       B12
    act_func_ptrs%function_ptrs(55)%ptr => GAMMAGD_CALC                !   55       GAMMAGD
    act_func_ptrs%function_ptrs(55)%propname = "GAMMAGD"
    act_func_ptrs%function_ptrs(56)%ptr => PR_CALC                     !   56       Residual pressure
    act_func_ptrs%function_ptrs(56)%propname = "Residual pressure"
    act_func_ptrs%function_ptrs(57)%ptr2 => dBdT_CALC                   !   57       dBdT (derivative of the second virial coeffcient)
    act_func_ptrs%function_ptrs(57)%propname = "dBdT"
    act_func_ptrs%function_ptrs(58)%ptr2 => dCdT_CALC                   !   58       dCdT (derivative of the third virial coeffcient)
    act_func_ptrs%function_ptrs(58)%propname = "dCdT"
    act_func_ptrs%function_ptrs(59)%ptr => dPdD_CALC                   !   59       dPdD (1st derivative of the pressure with respect to density at constant temperature)
    act_func_ptrs%function_ptrs(59)%propname = "dPdD"
    act_func_ptrs%function_ptrs(60)%ptr => d2PdD2_CALC                   !   60       d2PdD2 (2nd derivative of the pressure with respect to density at constant temperature)
    act_func_ptrs%function_ptrs(60)%propname = "d2PdD2"
    act_func_ptrs%function_ptrs(61)%ptr => D2PDDT_CALC                   !   61      d2PdDdT (1st mixed derivatives of the pressure with respect to density and temperature)
    act_func_ptrs%function_ptrs(61)%propname = "d2PdDdT"
    act_func_ptrs%function_ptrs(62)%ptr => d2PdTT_CALC                   !   62     d2PdT2 (2nd derivative of the pressure with respect to temperature at constant density)
    act_func_ptrs%function_ptrs(62)%propname = "d2PdDdT"
    act_func_ptrs%function_ptrs(63)%ptr => dDdT_CALC                   !   63       dDdT (1st derivative of the density with respect to temperature at constant pressure)
    act_func_ptrs%function_ptrs(63)%propname = "dDdT"
    act_func_ptrs%function_ptrs(64)%ptr3 => DSPIN_CALC                  !   64       DSPIN (Calculation of Spinodals)
    act_func_ptrs%function_ptrs(64)%propname = "DSPIN"
    act_func_ptrs%function_ptrs(68)%ptr4 => TBoyle_calc                  !   65    TBoyle_calc
    act_func_ptrs%function_ptrs(68)%propname = "TBoyle"
    act_func_ptrs%function_ptrs(69)%ptr4 => TJT_calc                  !   66    TJT_calc
    act_func_ptrs%function_ptrs(69)%propname = "TJT"
    act_func_ptrs%function_ptrs(70)%ptr4 => TJTINV_calc                  !   67    TJTINV_calc
    act_func_ptrs%function_ptrs(70)%propname = "TJTINV"

    act_func_ptrs%function_ptrs(71)%ptr =>  D_PIP_DD_CALC    !  71 First PIP derivative with respect to density (at constant temperature) (numerical derivative)
    act_func_ptrs%function_ptrs(71)%propname = "D_PIP_DD_CALC"

    act_func_ptrs%function_ptrs(72)%ptr => D2_PIP_D2D_CALC    !   72 Second PIP derivative with respect to density^2 (at constant temperature) (numerical derivative)
    act_func_ptrs%function_ptrs(72)%propname = "D2_PIP_D2D_CALC"

    act_func_ptrs%function_ptrs(73)%ptr2 =>  d2BdT2_CALC    ! 73  Second derivative of second virial coeff. B with respect to T^2
    act_func_ptrs%function_ptrs(73)%propname = "d2BdT2_CALC"

    act_func_ptrs%function_ptrs(74)%ptr2 => d2CdT2_CALC      ! 74 Second derivative of third virial coeff. C with respect to T^2
    act_func_ptrs%function_ptrs(74)%propname = "d2CdT2_CALC"

    act_func_ptrs%function_ptrs(75)%ptr2 =>  d2DdT2_CALC    !  75 Second derivative of fourth virial coeff. D with respect to T^2
    act_func_ptrs%function_ptrs(75)%propname = "d2DdT2_CALC"

    act_func_ptrs%function_ptrs(76)%ptr2 =>   d3BdT3_CALC   !  76 Third derivative of second virial coeff. B with respect to T^3
    act_func_ptrs%function_ptrs(76)%propname = "d3BdT3_CALC"

    act_func_ptrs%function_ptrs(77)%ptr2 => d3CdT3_CALC    !  77 Third derivative of third virial coeff. C with respect to T^3
    act_func_ptrs%function_ptrs(77)%propname = "d3CdT3_CALC"

    act_func_ptrs%function_ptrs(78)%ptr2 =>    d3DdT3_CALC  !  78 Third derivative of fourth virial coeff. D with respect to T^3
    act_func_ptrs%function_ptrs(78)%propname = "d3DdT3_CALC"

    act_func_ptrs%function_ptrs(79)%ptr2 => d4BdT4_CALC     !  79 Fourth derivative of second virial coeff. B with respect to T^4
    act_func_ptrs%function_ptrs(79)%propname = "d4BdT4_CALC"

    act_func_ptrs%function_ptrs(80)%ptr2 =>  d4CdT4_CALC    !  80 Fourth derivative of third virial coeff. C with respect to T^4
    act_func_ptrs%function_ptrs(80)%propname = "d4CdT4_CALC"

    act_func_ptrs%function_ptrs(81)%ptr2 =>  d4DdT4_CALC    !  81 Fourth derivative of fourth virial coeff. D with respect to T^4
    act_func_ptrs%function_ptrs(81)%propname = "d4DdT4_CALC"



    end subroutine set_func_pointers

    end module calc_func_ptr