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
#include "access_macros.fi"

    ! module for file vle_derivs.f90
    module vle_derivs_module
    !global use inclusion
    use module_all_types
    use calc_functions
    use rhomix_pt_module
    use flash_module
    use reduced_parameters_calc_module
    use cubic_eos_module
    use seawater_module
    use pc_saft_ARX_derivs_module
    use fnrderivs_module
    use setup_module
    use fniderivs_module

    contains




    !**************************************************************************
    !           --------------------------------------------------
    !           Routines for the calculation of the mixture
    !           derivatives of the pressure with respect to
    !           T, V, RHO and ni
    !           These routines are needed for phase equilibrium
    !           calculations.
    !
    !           A. Jäger, J. Gernert Aug. 2010 Boulder
    !           --------------------------------------------------
    !**************************************************************************



    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    subroutine dP_drho (gl,TEMPERATURE, DENSITY, dPdD)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE PRESSURE WITH RESPECT TO DENSITY
    ! AT CONSTANT TEMPERATURE AND MOLE NUMBERS / CONSTANT TEMPERATURE AND MOLAR COMPOSITION AS PUBLISHED BY
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 108 table 7.2)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dPdD        - Variable the derivative is stored in
    !--------------------------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision:: dPdD

    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivs):: getder_res
    double precision:: Fr_d_delta, Fr_d_delta_d_delta, Rmix


    getder_res = (/0,1,1,0,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.d0


    dPdD = 0.D0                         !Initialize return variable
    Fr_d_delta = 0.D0                   !Initialize
    Fr_d_delta_d_delta = 0.D0             !Initialize


    Call MIXDERIVSFNR (gl,temperature, density, getder_res, der_res)
    Fr_d_delta = der_res(2)             !Derivative of the residual Helmholtz with respect to delta
    Fr_d_delta_d_delta = der_res(3)     !Second derivative of the residual Helmholtz with respect to delta

    Call R_mix_calc(gl,Rmix)

    dPdD = Temperature *  Rmix * (1.d0 +  2.D0*Fr_d_delta + Fr_d_delta_d_delta)

    end subroutine


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    double precision function d2P_drho2 (gl,TEMPERATURE, DENSITY)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE PRESSURE WITH RESPECT TO DENSITY
    ! AT CONSTANT TEMPERATURE AND MOLE NUMBERS / CONSTANT TEMPERATURE AND MOLAR COMPOSITION AS PUBLISHED BY
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 108 table 7.2)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! ddPddD      - Variable the derivative is stored in
    !--------------------------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY


    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivs):: getder_res
    double precision:: dar_dd, d2ar_dd2, d3ar_dd3, Rmix


    getder_res = (/0,1,1,0,0,0,0,1,0,0,0,0,0,0,0/)
    der_res = 0.d0


    d2P_drho2 = 0.D0                             !Initialize return variable
    dar_dd = 0.D0                       !Initialize
    d2ar_dd2 = 0.D0               !Initialize
    d3ar_dd3= 0.D0        !Initialize

    Call MIXDERIVSFNR (gl,temperature, density, getder_res, der_res)
    dar_dd = der_res(2)                 !Derivative of the residual Helmholtz with respect to delta
    d2ar_dd2 = der_res(3)         !Second derivative of the residual Helmholtz with respect to delta
    d3ar_dd3 = der_res(8) !Third derivative of the residual Helmholtz with respect to delta


    Call R_mix_calc(gl,Rmix)

    d2P_drho2 = Temperature *  Rmix / density * (2.d0*dar_dd + 4.d0 * d2ar_dd2 + d3ar_dd3)
    end function d2P_drho2



    subroutine DEPFUNCFNR_caller(gl,temperature, density, getder, depfuncder, fld1_arg, fld2_arg, der_index, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
    type(type_gl) :: gl
    double precision:: temperature, density
    integer, dimension(nderivs):: getder            ! array specifier to indicate, which derivative is needed
    integer :: fld1, fld2, fld1_arg, fld2_arg                      ! line index of the array DFCOEFF which contains the coefficients of the binary departure functions
    integer :: der_index    ! which derivative !!! ATTENTION ONLY ONE SINGLE DERIVATIVE IS SUPPORTED
    double precision, dimension(nderivs)::depfuncder ! array with the computed values for the derivatives
    double precision, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij
    logical, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij_calculated

    ! ensure fld1 < fld2
    if (fld1_arg < fld2_arg) then
        fld1 = fld1_arg
        fld2 = fld2_arg
    else
        fld1 = fld2_arg
        fld2 = fld1_arg
    endif

    if (DEPFUNCFNR_ij_calculated(fld1, fld2)) then
        depfuncder(der_index) = DEPFUNCFNR_ij(fld1, fld2)
    else
        call DEPFUNCFNR(gl,temperature, density, getder, depfuncder, fld1_arg, fld2_arg)
        DEPFUNCFNR_ij(fld1, fld2) = depfuncder(der_index)
        DEPFUNCFNR_ij_calculated(fld1, fld2) = .true.
    endif
    end subroutine DEPFUNCFNR_caller
    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    subroutine dar_dxi (gl,TEMPERATURE, DENSITY, DERIVFX)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE MOLEFRACTION Xi FOR EACH COMPONENT AT CONSTANT  del, tau
    ! ACCORDING TO THE ALGORITHM BY KEN. HALL, THE MOLEFRACTION OF THE LAST
    ! COMPONENT Xn IS REPLACED BY 1 - SUM OF ALL OTHER MOLEFRACTIONS.
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A MATRIX (30)IN WHICH THE DERIVATIVE WITH RESPECT TO COMPOSITION FOR
    !               EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------










    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVFX

    double precision, dimension(nderivs):: der_res, der_dep
    integer, dimension(nderivs):: getder_res, getder_dep
    integer:: i, j, k
    double precision:: alphaij_r, alphaiN_r, alphajN_r, F0i_r, F0N_r, delr

    !Variable for SRK AND PR
    double precision, dimension(30):: dardxi
    double precision, dimension(15,30) :: MIXDERFNR_dxi

    !Variables for LKP
    double precision :: B_0, C_0, D_0, B_ref, C_ref, D_ref
    double precision, dimension(30) :: daccenLKP_dxi, dzcLKP_dxi
    double precision :: zcLKP2, zcLKP3, zcLKP4, zcLKP5,zcLKP6, zcLKP8, tau, tau2, tau3, del, del2, del3, del4, del5
    double precision, dimension(30) :: part_id, part_ref,part_dar_dxi_0, part_dar_dxi_ref

    !New variables for quadratic mixing rules for the residual reduced Helmholtz energy
    !Andreas, February 2016
    double precision:: alpha_r_ij            !Combination term alpha^r_(ij) = 0.5 * (alpha^r_i + alpha^r_j) (1-k_(ij))
    double precision:: alpha_r_jN, alpha_r_iN, alpha_r_NN
    double precision, dimension(30):: alpha_r_i         !reduced residual Helmholtz energies of the pure fluids

    !New variables for numerical derivative and analytical derivative of the excess based departure function
    !Andreas Jäger, July 2016
    double precision:: delta_x                          !Step in x for numerical derivative
    double precision:: alpha_dep_p, alpha_dep_m         !Value of departure function at increased and decreased xi (and decreased or increased xN, respectively)
    double precision, dimension(30):: molfractions_orig !Original value of molfractions
    double precision:: Temp_p, Temp_m, dens_p, dens_m   !Decrease and increased temperature and density. When xi is changed and tau and delta stay constant, T and rho change
    double precision:: tredmix_orig, rhoredmix_orig
    double precision:: rhomix_ref_copy
    double precision, dimension(30):: rho_i_ref_copy
    !double precision, dimension(30):: DERIVFX_num
    !Andreas Jäger, June 2017
    double precision:: const_A1
    double precision:: help_h, help_f
    double precision, dimension(30):: dhelph_dxa, dhelpf_dxa
    double precision, dimension(30):: alpha_oi_r_mix, alpha_oi_r_i
    integer, dimension(nderivs):: GETDER_i                         ! array specifier to indicate, which derivative is needed for the pure fluid residual Helmholtz energies
    double precision, dimension(nderivs)::FNRDER_i                 ! array with the computed values for the derivatives of the residual Helmholtz energy of the pure fluids
    double precision:: aE
    double precision, dimension(30):: dalpha_oi_r_mix_ddelref, dalpha_oi_r_i_dtaui
    double precision, dimension(30):: ddelref_dxa, dtaui_dxa
    double precision, dimension(30):: dTred_dxi, drhored_dxi
    double precision, dimension(30):: dbHelmgE_dxi
    double precision, dimension(30):: daE_dxi
    !Variables for UNIFAC
    !double precision, dimension(nderivs) :: gl%ge%gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gl%ge%gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau
    integer:: C_or_R
    integer:: errval
    !double precision, dimension(nderivs, 30) :: gl%ge%gE_C_dxa        !Combinatorial part of gE and derivatives with respect to tau and delta and xa
    !double precision, dimension(nderivs, 30) :: gl%ge%gE_R_dxa        !Residual part of gE and derivatives with respect to tau and delta and xa
    !double precision, allocatable :: gl%ge%ln_gamma_C_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to delta and tau and xa
    !double precision, allocatable :: gl%ge%ln_gamma_R_dxa(:,:,:)       !Residual activity coefficients and derivatives with respect to delta and tau and xa



    !New variables for PC-SAFT equation type
    double precision, dimension (10,gl%ncomp) :: x1

    !Erik, April 2018
    !New variables for COSMO-SAC
    double precision, dimension(nderivs) :: gE_C_p, gE_C_m, gE_R_p, gE_R_m
    double precision, dimension(:,:), allocatable :: ln_gamma_C_p, ln_gamma_C_m, ln_gamma_R_p, ln_gamma_R_m

    double precision, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij
    logical, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij_calculated

    if (.not. allocated(ln_gamma_C_p)) allocate(ln_gamma_C_p(nderivs, 30))
    allocate(ln_gamma_C_m, ln_gamma_R_p, ln_gamma_R_m, mold=ln_gamma_C_p)

    getder_res = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    getder_dep = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.d0
    der_dep = 0.d0

    DERIVFX = 0.D0               !Initialize return matrix DERIVFX

    alphaij_r = 0.D0                    !Initialize departure function ij
    alphaiN_r = 0.D0                    !Initialize departure function iN
    alphajN_r = 0.D0                    !Initialize departure function jN
    F0i_r = 0.D0                        !Initialize residual helmholtz pure of component i
    F0N_r = 0.D0                        !Initialize residual helmholtz pure of component N
    aE = 0.D0

    if (gl%mix_Type  ==  1) then !Lorentz Berthelot or modified mixing rules used Andreas March 2012

        DEPFUNCFNR_ij_calculated = .false.

        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, gl%NCOMP)
        F0N_r = der_res(1)                  !Get residual helmh. pure of comp. N

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components
            Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
            F0i_r = der_res(1)                      !residual part of component i: F0i_r

            Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP, 1, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
            alphaiN_r = der_dep(1)                  !Get departure Function for Comp. i and N
            DERIVFX(i) = F0i_r - F0N_r + (1.D0 - 2.D0 * gl%molfractions(i)) * gl%Fij(i,gl%NCOMP)* alphaiN_r
            Do j = 1, (gl%NCOMP-1)                         !Sum over all components with departure function
                if(i  /=  j) then
                    if((gl%Fij(i,j)  /=  0.d0).OR.(gl%Fij(i,gl%NCOMP)  /=  0.d0).OR.(gl%Fij(j,gl%NCOMP)  /=  0.d0)) then
                        Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, j, 1, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                        alphaij_r = der_dep(1)
                        Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, j, gl%NCOMP, 1, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                        alphajN_r = der_dep(1)          !Get departure Function for Comp. i and N

                        DERIVFX(i) = DERIVFX(i)+ gl%molfractions(j)*(gl%Fij(i,j) * alphaij_r  - &
                            & gl%Fij(i,gl%NCOMP) * alphaiN_r - gl%Fij(j,gl%NCOMP) * alphajN_r)
                    end if
                end if
            end do
        end do

        !SRK with SRK mixing rules used
    else if ((gl%mix_type  ==  2) .or. (gl%mix_type == 21) .or. (gl%mix_type == 22)) then
        !Replaced with new cubics routine
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_dxi_CUBIC (gl,Temperature, density, getder_res, MIXDERFNR_dxi)
        DERIVFX = MIXDERFNR_dxi(1,:)
        !Old version
        !call dar_dxi_SRK(gl,Temperature, density, dardxi)
        !DERIVFX = dardxi

        !PR with PR mixing rules used, not yet included!!!
    else if ((gl%mix_type  ==  3) .or. (gl%mix_type == 31)) then
        !Replaced with new cubics routine
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_dxi_CUBIC (gl,Temperature, density, getder_res, MIXDERFNR_dxi)
        DERIVFX = MIXDERFNR_dxi(1,:)
        !Old version
        !call dar_dxi_PR(Temperature, density, dardxi)
        !DERIVFX = dardxi

        !Stefan & Andreas July 2014, edited by Monika 10/2014
        !LKP with LKP mixing rules
    else if (gl%mix_type  ==  4) then

        tau = gl%tredmix/temperature
        del = density/gl%rhoredmix

        del2=del*del
        del3=del2*del
        del4=del*del3
        del5=del4*del
        zcLKP2=gl%zcLKP*gl%zcLKP
        zcLKP3=zcLKP2*gl%zcLKP
        zcLKP4=zcLKP2*zcLKP2
        zcLKP5=gl%zcLKP*zcLKP4
        zcLKP6=gl%zcLKP*zcLKP5
        zcLKP8=zcLKP6*zcLKP2
        tau2=tau*tau
        tau3=tau*tau2


        do i = 1, gl%ncomp - 1
            daccenLKP_dxi(i) = gl%accen(i) - gl%accen(gl%ncomp)
            dzcLKP_dxi(i)=-0.085d0 * daccenLKP_dxi(i)
        end do

        B_0 = gl%lkp_b1_0 - gl%lkp_b2_0 * tau - gl%lkp_b3_0 * tau2 - gl%lkp_b4_0 * tau3
        C_0 = gl%lkp_c1_0 - gl%lkp_c2_0 * tau + gl%lkp_c3_0 * tau3
        D_0 = gl%lkp_d1_0 + gl%lkp_d2_0 * tau
        B_ref = gl%lkp_b1_ref - gl%lkp_b2_ref * tau - gl%lkp_b3_ref * tau2 - gl%lkp_b4_ref * tau3
        C_ref = gl%lkp_c1_ref - gl%lkp_c2_ref * tau + gl%lkp_c3_ref * tau3
        D_ref = gl%lkp_d1_ref + gl%lkp_d2_ref * tau

        do i = 1, gl%ncomp - 1

            !Helmholtz energy
            !Ideal fluid contribution
            part_id(i) = B_0 / gl%zcLKP * del + 0.5d0 * C_0 / zcLKP2 * del2 + 0.2d0 * D_0 / zcLKP5 * del5 -  &
                & gl%lkp_c4_0 / (2.d0 * gl%lkp_gamma_0) * tau3 * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0 + 1.d0) * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2)  &
                & + gl%lkp_c4_0 / (2.d0 * gl%lkp_gamma_0) * tau3 * (gl%lkp_beta_0 + 1.d0)

            !Reference fluid contribution
            part_ref(i) = B_ref / gl%zcLKP * del + 0.5d0 * C_ref / zcLKP2 * del2 + 0.2d0 * D_ref / zcLKP5 * del5 -  &
                & gl%lkp_c4_ref / (2.d0 * gl%lkp_gamma_ref) * tau3 * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref + 1.d0) * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2)  &
                & + gl%lkp_c4_ref / (2.d0 * gl%lkp_gamma_ref) * tau3 * (gl%lkp_beta_ref + 1.d0)


            !partial derivative
            !Ideal fluid contribution
            part_dar_dxi_0(i) = (- B_0 / zcLKP2 * del - C_0 / zcLKP3 * del2 - D_0 / zcLKP6 * del5) * dzcLKP_dxi(i)  &
                - gl%lkp_c4_0 * tau3 / (2.d0 * gl%lkp_gamma_0) * dzcLKP_dxi(i) * 2.d0 * del2 * gl%lkp_gamma_0 / zcLKP3 * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) &
                * (gl%lkp_gamma_0 / zcLKP2 *del2 + gl%lkp_beta_0)

            !Reference fluid contribution
            part_dar_dxi_ref(i) = (- B_ref / zcLKP2 * del - C_ref / zcLKP3 * del2 - D_ref / zcLKP6 * del5) * dzcLKP_dxi(i) &
                - gl%lkp_c4_ref * tau3 / (2.d0 * gl%lkp_gamma_ref) * dzcLKP_dxi(i) * 2.d0 * del2 * gl%lkp_gamma_ref / zcLKP3 * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) &
                * (gl%lkp_gamma_ref / zcLKP2 *del2 + gl%lkp_beta_ref)

            !calculate Dalpha_r_Dx_i(i)
            DERIVFX(i) = - daccenLKP_dxi(i) / gl%lkp_w_ref * part_id(i) + (1.d0 - gl%accenLKPMix / gl%lkp_w_ref) * part_dar_dxi_0(i) &
                + daccenLKP_dxi(i) / gl%lkp_w_ref * part_ref(i) + gl%accenLKPMix / gl%lkp_w_ref * part_dar_dxi_ref(i)
        end do

        !PC-SAFT
    else if (gl%mix_type == 6) then

        call ARX1DERIVS(gl,temperature, density, getder_res, x1)
        DERIVFX(1:gl%ncomp) = x1(1,:)

        !transformation to usual writing
        do i = 1, gl%ncomp
            DERIVFX(i) = DERIVFX(i)- DERIVFX(gl%ncomp)
        end do


        !Quadratic mixing rules for the residual reduced Helmholtz energy, Andreas February 2015
    else if (gl%mix_type  ==  11) then

        DEPFUNCFNR_ij_calculated = .false.

        !First, calculate all reduced residual Helmholtz energies of the pure fluids
        der_res = 0.D0
        do i = 1,gl%ncomp
            call FNRDERIVS(gl,temperature, density, getder_res, der_res, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_r_i(i) = der_res(1)    !Residual reduced Helmholtz energy of pure component i in the mixture
        end do

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xncomp-1
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components

            !Contribution of the quadratic mixing term for the pure reduced residual Helmholtz energies
            alpha_r_NN = alpha_r_i(gl%ncomp)
            alpha_r_iN = 0.5D0 * (alpha_r_i(i) + alpha_r_i(gl%ncomp)) * (1.D0 - ACCESS_KIJ_HELM(i,gl%ncomp))
            DERIVFX(i) = DERIVFX(i) + 2.D0 * gl%molfractions(gl%ncomp) * (alpha_r_iN - alpha_r_NN)

            !Contribution of the departure function
            Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP, 1, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
            alphaiN_r = der_dep(1)                  !Get departure Function for Comp. i and N
            DERIVFX(i) = DERIVFX(i) + (1.D0 - 2.D0 * gl%molfractions(i)) * gl%Fij(i,gl%NCOMP)* alphaiN_r

            Do j = 1, (gl%NCOMP-1)                         !Sum over all components with departure function
                !Contribution of the quadratic mixing term for the pure reduced residual Helmholtz energies
                if (i == j) then
                    alpha_r_ij = alpha_r_i(i)
                end if
                if (i /= j) then
                    alpha_r_ij = 0.5D0 * (alpha_r_i(i) + alpha_r_i(j)) * (1.D0 - ACCESS_KIJ_HELM(i,j))
                end if
                alpha_r_jN = 0.5D0 * (alpha_r_i(j) + alpha_r_i(gl%ncomp)) * (1.D0 - ACCESS_KIJ_HELM(j,gl%ncomp))
                DERIVFX(i) = DERIVFX(i) + 2.D0 * gl%molfractions(j) * (alpha_r_ij - alpha_r_jN)

                !Contributions of the departure function
                if(i  /=  j) then
                    !if((Fij(i,j)  /=  0.d0).OR.(Fij(i,NCOMP)  /=  0.d0).OR.(Fij(j,NCOMP)  /=  0.d0)) then
                    Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, j, 1, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                    alphaij_r = der_dep(1)
                    Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, j, gl%NCOMP, 1, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                    alphajN_r = der_dep(1)          !Get departure Function for Comp. i and N

                    DERIVFX(i) = DERIVFX(i)+ gl%molfractions(j)*(gl%Fij(i,j) * alphaij_r  - &
                        & gl%Fij(i,gl%NCOMP) * alphaiN_r - gl%Fij(j,gl%NCOMP) * alphajN_r)
                    !end if
                end if
            end do
        end do

        !Excess based departure function, Andreas Jäger July 2016, modified April 2018, Erik
    else if ((gl%mix_type  ==  12) .or. (gl%mix_type  ==  13)) then
        if(.not. allocated(gl%ge%ln_gamma_C_dxa)) allocate(gl%ge%ln_gamma_C_dxa(nderivs, 30, 30))
        if(.not. allocated(gl%ge%ln_gamma_R_dxa)) allocate(gl%ge%ln_gamma_R_dxa(nderivs, 30, 30))

        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, gl%NCOMP)
        F0N_r = der_res(1)                  !Get residual helmh. pure of comp. N

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components
            Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
            F0i_r = der_res(1)                      !residual part of component i: F0i_r
            DERIVFX(i) = F0i_r - F0N_r
        end do
        !DERIVFX_num = DERIVFX

        tredmix_orig = gl%tredmix
        rhoredmix_orig = gl%rhoredmix
        tau = gl%tredmix / temperature
        del = density / gl%rhoredmix
        !Calculate analytical derivative of the excess based departure function with respect to xi
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !Compute auxiliary function f(del,x)
        !----
        const_A1 = -dlog(1.D0 / gl%u_pack + 1.D0)
        help_f = -dlog(gl%b_HelmgE * density + 1.D0) / const_A1
        !----

        !Compute derivative of auxiliary function f(del,x) with respect to xa
        !----
        call dYr_dxi(gl,dTred_dxi, drhored_dxi)
        call db_HelmhE_dxi(gl,dbHelmgE_dxi)
        do j = 1, gl%ncomp-1
            dhelpf_dxa(j) = -del * (gl%rhoredmix * dbHelmgE_dxi(j) + drhored_dxi(j) * gl%b_HelmgE) &
                & / const_A1 / (1.D0 + del * gl%b_HelmgE * gl%rhoredmix)
        end do
        !----

        !Compute auxiliary function h(tau,x)
        !----
        GETDER_i = 0
        GETDER_i(1) = 1 !The Helmholtz energy itself is needed
        GETDER_i(2) = 1 !Get the first derivative of the Helmholtz energy with respect to deltaref
        do i = 1, gl%ncomp
            !Get the reduced Helmholtz energies of the pure components at mixture conditions
            rhomix_ref_copy = gl%rho_mix_ref
            call FNRDERIVS(gl,temperature, rhomix_ref_copy, GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_oi_r_mix(i) = FNRDER_i(1)
            dalpha_oi_r_mix_ddelref(i) = FNRDER_i(2) / (gl%rho_mix_ref/gl%rhoredmix)
        end do
        GETDER_i(2) = 0
        GETDER_i(4) = 1 !Get the first derivative of the Helmholtz energy with respect to taui of the pure fluid part
        do i = 1, gl%ncomp
            !Get the reduced Helmholtz energies of the pure components at pure fluid conditions
            gl%tredmix = gl%tc(i)         !Needed to evaluate the pure fluid at correct tau (tau_i)
            gl%rhoredmix = gl%rhoc(i)     !Needed to evaluate the pure fluid at correct delta (delta_i_ref)
            rho_i_ref_copy(i) = gl%rho_i_ref(i)
            call FNRDERIVS(gl,temperature, rho_i_ref_copy(i), GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_oi_r_i(i) = FNRDER_i(1)
            dalpha_oi_r_i_dtaui(i) = FNRDER_i(4) / (gl%tc(i)/temperature)
        end do
        gl%tredmix = tredmix_orig
        gl%rhoredmix = rhoredmix_orig

        C_or_R = 2  !Get residual part of UNIFAC/COSMO-SAC, modified April 2018, Erik
        !call gE_UNIFAC_MIXDERIVS(gl,temperature, getder_dep, gl%ge%gE_C, gl%ge%gE_R, gl%ge%ln_gamma_C, gl%ge%ln_gamma_R, C_or_R, errval)
        if (gl%mix_type  ==  12) then
            call gE_UNIFAC_MIXDERIVS(gl,temperature, getder_dep, C_or_R, errval)
        elseif (gl%mix_type  ==  13) then
            call gE_COSMO_SAC_MIXDERIVS(gl,temperature, getder_dep, C_or_R, errval)
        end if
        aE = gl%ge%gE_C(1) + gl%ge%gE_R(1)    !Specify which parts of gE are considered (combinatorial, residual, or both)

        help_h = aE / R_HelmgE / temperature
        do i = 1, gl%ncomp
            !Andreas Jäger, December 2017
            !term "sum(xi ln(bi/b))" deleted from the excess based departure function, because it is assumed that it cancels out with the combinatorial part of gE
            !help_h = help_h - molfractions(i) * dlog(bi_HelmgE(i) / b_HelmgE) - molfractions(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
            help_h = help_h - gl%molfractions(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
        end do
        !----

        !Compute derivative of auxiliary function h(tau,x) with respect to xa
        !----
        !call gE_UNIFAC_MIXDERIVS_dxa(gl,temperature, getder_dep, gl%ge%gE_C_dxa, gl%ge%gE_R_dxa, gl%ge%ln_gamma_C_dxa, gl%ge%ln_gamma_R_dxa, C_or_R, errval)
        if (gl%mix_type  ==  12) then
            call gE_UNIFAC_MIXDERIVS_dxa(gl,temperature, getder_dep, C_or_R, errval)
            !Erik, April 2018
        elseif (gl%mix_type  ==  13) then
            !!Analytical derivative of the excess based departure function with respect to xi at constant tau and del
            if ((gl.cosmo.COSMO_ver == 1) .and. (gl.cosmo.analytical)) then
                call gE_COSMO_SAC_MIXDERIVS_dxa(gl,temperature, getder_dep, C_or_R, errval)
            else
                !!Numerical derivative of the excess based departure function with respect to xi at constant tau and del
                delta_x = 1.D-4
                molfractions_orig = gl%molfractions
                tredmix_orig = gl%tredmix
                rhoredmix_orig = gl%rhoredmix
                do i =1, gl%ncomp-1
                    !Increase xi
                    gl%molfractions(i) = molfractions_orig(i) + delta_x
                    gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) - delta_x
                    call reduced_parameters_calc(gl, Temperature)
                    Temp_p = gl%tredmix/tredmix_orig * Temperature
                    !Dens_p = rhoredmix/rhoredmix_orig * Density
                    call gE_COSMO_SAC_MIXDERIVS(gl, Temp_p, getder_dep, C_or_R, errval)
                    gE_C_p = gl%ge%gE_C
                    gE_R_p = gl%ge%gE_R
                    ln_gamma_C_p = gl%ge%ln_gamma_C
                    ln_gamma_R_p = gl%ge%ln_gamma_R
                    !test
                    !seggamma_p = gl.cosmo.seggamma_gl
                    !    alpha_dep_p = der_dep(1)
                    !Decrease xi
                    gl%molfractions(i) = molfractions_orig(i) - delta_x
                    gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) + delta_x
                    call reduced_parameters_calc(gl, Temperature)
                    Temp_m = gl%tredmix/tredmix_orig * Temperature
                    !    Dens_m = rhoredmix/rhoredmix_orig * Density
                    call gE_COSMO_SAC_MIXDERIVS(gl, Temp_m, getder_dep, C_or_R, errval)
                    gE_C_m = gl%ge%gE_C
                    gE_R_m = gl%ge%gE_R
                    ln_gamma_C_m = gl%ge%ln_gamma_C
                    ln_gamma_R_m = gl%ge%ln_gamma_R
                    !test
                    !seggamma_m = gl.cosmo.seggamma_gl
                    !    alpha_dep_m = der_dep(1)
                    !Calculate the numerical derivative
                    gl%ge%gE_C_dxa(1,i) = (gE_C_p(1) - gE_C_m(1)) / (2.D0 * delta_x)
                    gl%ge%gE_R_dxa(1,i) = (gE_R_p(1) - gE_R_m(1)) / (2.D0 * delta_x)
                    !test
                    !seggamma_dxa = (seggamma_p - seggamma_m) / (2.D0 * delta_x)
                    gl.molfractions = molfractions_orig
                end do

                call reduced_parameters_calc(gl, Temperature)

                !test
                !!Derivatives of tredmix with respect to xa are required for the following derivatives
                !call dYr_dxi(gl,dTred_dxa, drhored_dxa)
                !GETDER = 1
                !call gE_COSMO_SAC_MIXDERIVS(gl, Temperature, GETDER, C_or_R, errval)
                !gE_R_Trho(4) = - tau**2 / gl.tredmix * gl.gE_R(4)
                !Do i = 1,gl.ncomp-1
                !    gE_R_Trho_dxa(1,i) =  gl.gE_R_dxa(1,i) - dTred_dxa(i) / tau * gE_R_Trho(4)
                !end do
            end if
        end if
        daE_dxi = gl%ge%gE_C_dxa(1,:) + gl%ge%gE_R_dxa(1,:)
        do j = 1, gl%ncomp-1
            !Andreas Jäger, December 2017
            !term "sum(xi ln(bi/b))" deleted from the excess based departure function, because it is assumed that it cancels out with the combinatorial part of gE
            !dhelph_dxa(j) = tau / R_HelmgE / tredmix**2 * (tredmix * daE_dxi(j) - aE * dTred_dxi(j)) &
            !            & + dlog(bi_HelmgE(ncomp)/bi_HelmgE(j)) + dbHelmgE_dxi(j) / b_HelmgE &
            !            & - (alpha_oi_r_mix(j) - alpha_oi_r_i(j)) + (alpha_oi_r_mix(ncomp) - alpha_oi_r_i(ncomp))
            dhelph_dxa(j) = tau / R_HelmgE / gl%tredmix**2 * (gl%tredmix * daE_dxi(j) - aE * dTred_dxi(j)) &
                & - (alpha_oi_r_mix(j) - alpha_oi_r_i(j)) + (alpha_oi_r_mix(gl%ncomp) - alpha_oi_r_i(gl%ncomp))
            ddelref_dxa(j) = - 1.D0 / gl%u_pack * (drhored_dxi(j) / gl%b_HelmgE / gl%rhoredmix**2 + dbHelmgE_dxi(j) / gl%b_HelmgE**2 / gl%rhoredmix)
            do i = 1, gl%ncomp
                dtaui_dxa(j) =  - tau * gl%tc(i) / gl%tredmix**2 * dTred_dxi(j)
                dhelph_dxa(j) = dhelph_dxa(j) - gl%molfractions(i) * dalpha_oi_r_mix_ddelref(i) * ddelref_dxa(j) &
                    & + gl%molfractions(i) * dalpha_oi_r_i_dtaui(i) * dtaui_dxa(j)
            end do
        end do
        !----

        DERIVFX = DERIVFX + dhelpf_dxa * help_h + help_f * dhelph_dxa
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !!Numerical derivative of the excess based departure function with respect to xi at constant tau and del
        !delta_x = 1.D-4
        !molfractions_orig = molfractions
        !tredmix_orig = tredmix
        !rhoredmix_orig = rhoredmix
        !do i =1, ncomp-1
        !    !Increase xi
        !    molfractions(i) = molfractions_orig(i) + delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !    call reduced_parameters_calc(gl,Temperature)
        !    Temp_p = tredmix/tredmix_orig * Temperature
        !    Dens_p = rhoredmix/rhoredmix_orig * Density
        !    call DEPFUNC_GE_BASED (Temp_p, Dens_p, getder_dep, der_dep)
        !    alpha_dep_p = der_dep(1)
        !    !Decrease xi
        !    molfractions(i) = molfractions_orig(i) - delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !    call reduced_parameters_calc(gl,Temperature)
        !    Temp_m = tredmix/tredmix_orig * Temperature
        !    Dens_m = rhoredmix/rhoredmix_orig * Density
        !    call DEPFUNC_GE_BASED (Temp_m, Dens_m, getder_dep, der_dep)
        !    alpha_dep_m = der_dep(1)
        !    !Calculate the numerical derivative
        !    DERIVFX_num(i) = DERIVFX_num(i) + (alpha_dep_p - alpha_dep_m) / (2.D0 * delta_x)
        !    molfractions = molfractions_orig
        !end do
        !call reduced_parameters_calc(gl,Temperature)
        !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        !One-fluid model for multiparameter EOS, Andreas Jäger July 2016
    else if (gl%mix_type == 19) then

        !!Numerical derivative of the one-fluid model at constant tau and del (same as constant T and rho)
        !delta_x = 1.D-4
        !molfractions_orig = molfractions
        !
        !do i =1, ncomp-1
        !    !Increase xi
        !    molfractions(i) = molfractions_orig(i) + delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !    call MIXDERIVSFNR (Temperature, density, getder_dep, der_dep)
        !    alpha_dep_p = der_dep(1)
        !    !Decrease xi
        !    molfractions(i) = molfractions_orig(i) - delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !    call MIXDERIVSFNR (Temperature, density, getder_dep, der_dep)
        !    alpha_dep_m = der_dep(1)
        !    !Calculate the numerical derivative
        !    DERIVFX(i) = DERIVFX(i) + (alpha_dep_p - alpha_dep_m) / (2.D0 * delta_x)
        !    molfractions = molfractions_orig
        !end do


        !Analytical derivative of the one-fluid model at constant tau and del (same as constant T and rho)
        gl%tredmix = gl%tc(gl%ncomp)
        gl%rhoredmix = gl%rhoc(gl%ncomp)
        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, gl%NCOMP)
        gl%tredmix = 1.D0
        gl%rhoredmix = 1.D0
        F0N_r = der_res(1)                  !Get residual helmh. pure of comp. N

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components
            gl%tredmix = gl%tc(i)
            gl%rhoredmix = gl%rhoc(i)
            Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
            gl%tredmix = 1.D0
            gl%rhoredmix = 1.D0
            F0i_r = der_res(1)                      !residual part of component i: F0i_r
            Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP)
            alphaiN_r = der_dep(1)                  !Get departure Function for Comp. i and N
            DERIVFX(i) = 2.D0 * gl%molfractions(i) * F0i_r - 2.D0 * gl%molfractions(gl%ncomp) * F0N_r + 2.D0 * gl%molfractions(gl%ncomp) * alphaiN_r
            Do j = 1, (gl%NCOMP-1)
                Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, j, gl%NCOMP)
                alphajN_r = der_dep(1)          !Get departure Function for Comp. i and N
                DERIVFX(i) = DERIVFX(i) - 2.D0 * gl%molfractions(j) * alphajN_r
                if(i  /=  j) then
                    Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, i, j)
                    alphaij_r = der_dep(1)
                    DERIVFX(i) = DERIVFX(i) + 2.D0 * gl%molfractions(j) * alphaij_r
                end if
            end do
        end do


    end if

    end subroutine dar_dxi



    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    subroutine dar_dxi_TV (gl,TEMPERATURE, DENSITY, DERIVFXTV)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE MOLEFRACTION Xi FOR EACH COMPONENT.
    ! ACCORDING TO THE ALGORITHM BY KEN. HALL, THE MOLEFRACTION OF THE LAST
    ! COMPONENT Xn IS REPLACED BY 1 - SUM OF ALL OTHER MOLEFRACTIONS.
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A MATRIX (30)IN WHICH THE DERIVATIVE WITH RESPECT TO COMPOSITION FOR
    !               EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVFXTV

    double precision, dimension(nderivs):: der_mix
    integer, dimension(nderivs):: getder_mix

    double precision, dimension(30):: d_rhored_d_xi, d_tred_d_xi, dardxi
    double precision:: dar_ddel, dar_dtau, del, tau

    integer:: i

    getder_mix = (/0,1,0,1,0,0,0,0,0,0,0,0,0,0,0/)
    der_mix = 0.d0

    DERIVFXTV = 0.D0                 !Initialize return matrix DERIVFX

    dar_ddel = 0.D0                          !Initialize residual helmholtz pure of component i
    dar_dtau = 0.D0                      !Initialize residual helmholtz pure of component N

    Call MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, getder_mix, der_mix)
    dar_ddel = der_mix(2)                     !Get residial helmh. pure of comp. N
    dar_dtau = der_mix(4)

    Call dYr_dxi(gl,d_tred_d_xi, d_rhored_d_xi)
    Call dar_dxi (gl,TEMPERATURE, DENSITY, dardxi)
    del = Density / gl%rhoredmix
    tau = gl%tredmix / temperature
    !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
    Do i = 1, (gl%NCOMP-1)                         !Loop over all components
        DERIVFXTV(i) = dar_ddel  / del * Density* d_rhored_d_xi(i) / (-gl%rhoredmix**2) + &
            & dar_dtau / tau * d_tred_d_xi(i) / Temperature + dardxi(i)
    end do

    end subroutine dar_dxi_TV

    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    subroutine d2ar_dxi2 (gl,TEMPERATURE, DENSITY, DERIVFX2)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE MOLEFRACTION Xi FOR EACH COMPONENT.
    ! ACCORDING TO THE ALGORITHM BY KEN. HALL, THE MOLEFRACTION OF THE LAST
    ! COMPONENT Xn IS REPLACED BY 1 - SUM OF ALL OTHER MOLEFRACTIONS.
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX2    - A MATRIX (30)IN WHICH THE SECOND DERIVATIVE WITH RESPECT TO COMPOSITION FOR
    !               EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------









    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVFX2

    double precision, dimension(nderivs):: der_dep, der_res
    integer, dimension(nderivs):: getder_dep, getder_res

    double precision:: alphaiN_r, delr

    integer::i, j

    !Variable for SRK and PR
    double precision, dimension(30,30):: d2ardxidxj
    !SH variable is part of gl
    !double precision, allocatable ::MIXDERFNR_dxidxj(:,:,:)

    !Variables for LKP
    double precision :: B_0, C_0, D_0, B_ref, C_ref, D_ref
    double precision, dimension(30) :: daccenLKP_dxi, dzcLKP_dxi
    double precision :: zcLKP2, zcLKP3, zcLKP4, zcLKP5,zcLKP6, zcLKP7, zcLKP8, tau, tau2, tau3, del, del2, del4, del5
    !double precision :: help_11, help_12, help_21, help_22
    !double precision, dimension(30,30) :: help_1, help_2, help_3
    !double precision, dimension(30,30) :: help_1_0, help_2_0, help_3_0

    double precision, dimension(30) :: part_dar_dxi_ref, part_dar_dxi_0, help_0, help_ref
    double precision, dimension(30,30) :: part_2dar_dxidxj_0, part_2dar_dxidxj_ref, help_ref_d2x, help_0_d2x

    !New variables for quadratic mixing rules for the residual reduced Helmholtz energy
    !Andreas, February 2016
    double precision:: alpha_r_ij            !Combination term alpha^r_(ij) = 0.5 * (alpha^r_i + alpha^r_j) (1-k_(ij))
    double precision:: alpha_r_jN, alpha_r_iN, alpha_r_NN
    double precision, dimension(30):: alpha_r_i         !reduced residual Helmholtz energies of the pure fluids

    !New variables for numerical derivative of the excess based departure function
    !Andreas Jäger, July 2016
    double precision:: delta_x                          !Step in x for numerical derivative
    double precision:: alpha_dep_p, alpha_dep_m         !Value of departure function at increased and decreased xi (and decreased or increased xN, respectively)
    double precision:: alpha_dep                        !Value of the departure function at original T, rho and x
    double precision, dimension(30):: molfractions_orig !Original value of molfractions
    double precision:: Temp_p, Temp_m, dens_p, dens_m   !Decrease and increased temperature and density. When xi is changed and tau and delta stay constant, T and rho change
    double precision:: tredmix_orig, rhoredmix_orig

    !New variables for PC-SAFT equation type
    double precision, allocatable  :: x2(:,:,:)
    double precision, dimension(gl%ncomp,gl%ncomp):: DERIVFX2_TRANS

    !double precision, dimension(30):: DERIVFX2_NUM

    !Erik, April 2018
    !New variables for COSMO-SAC
    integer:: errval, C_or_R
    double precision, dimension(nderivs) :: gE_C_p, gE_C_m, gE_R_p, gE_R_m
    double precision, dimension(nderivs, 30) :: ln_gamma_C_p, ln_gamma_C_m, ln_gamma_R_p, ln_gamma_R_m
    double precision, dimension(nderivs, 30) :: gE_C_dxa        !Combinatorial part of gE and derivatives with respect to tau and delta and xa
    double precision, dimension(nderivs, 30) :: gE_R_dxa        !Residual part of gE and derivatives with respect to tau and delta and xa
    !double precision, dimension(nderivs) :: gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau

    !test
    !double precision :: gE_R_dxa_p, gE_R_dxa_m
    !double precision, dimension(30) :: gE_R_dxadxb

    if (.not. allocated(gl%MIXDERIVFNR_dxidxj)) allocate(gl%MIXDERIVFNR_dxidxj(15,30,30))

    getder_res = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    getder_dep = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.D0
    der_dep = 0.D0
    DERIVFX2 = 0.d0

    if (gl%mix_Type  ==  1) then !Lorentz Berthelot or modified mixing rules used Andreas March 2012
        Do i = 1, gl%NCOMP-1
            Call DEPFUNCFNR (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP)
            alphaiN_r = der_dep(1)                          !Get departure Function for Comp. i and N
            DERIVFX2(i) = -2.D0*gl%Fij(i,gl%NCOMP)* alphaiN_r
        End do

    else if ((gl%mix_type  ==  2) .or. (gl%mix_type == 21) .or. (gl%mix_type == 22)) then !SRK with SRK mixing rules used
        !Replaced with new cubics routine
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_dxidxj_CUBIC (gl,temperature, density, getder_res)
        Do i = 1, gl%ncomp - 1
            DERIVFX2(i) = gl%MIXDERIVFNR_dxidxj(1,i,i)
        End do
        !Old version
        !call d2ar_dxidxj_SRK(gl,temperature, density, d2ardxidxj)
        !Do i = 1, ncomp - 1
        !    DERIVFX2(i) = d2ardxidxj(i,i)
        !End do

    else if ((gl%mix_type  ==  3) .or. (gl%mix_type == 31)) then !PR with PR mixing rules
        !Replaced with new cubics routine
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_dxidxj_CUBIC (gl,temperature, density, getder_res)
        Do i = 1, gl%ncomp - 1
            DERIVFX2(i) = gl%MIXDERIVFNR_dxidxj(1,i,i)
        End do
        !Old version
        !call d2ar_dxidxj_PR(temperature, density, d2ardxidxj)
        !Do i = 1, ncomp - 1
        !    DERIVFX2(i) = d2ardxidxj(i,i)
        !End do

    else if (gl%mix_type  ==  4) then !LKP with LKP mixing rules


        tau = gl%tredmix/temperature
        del = density/gl%rhoredmix

        del2=del*del
        del4=del2*del2
        del5=del4*del
        zcLKP2=gl%zcLKP*gl%zcLKP
        zcLKP3=zcLKP2*gl%zcLKP
        zcLKP4=zcLKP2*zcLKP2
        zcLKP5=gl%zcLKP*zcLKP4
        zcLKP6=gl%zcLKP*zcLKP5
        zcLKP7=zcLKP6*gl%zcLKP
        zcLKP8=zcLKP7*gl%zcLKP
        tau2=tau*tau
        tau3=tau*tau2


        do i = 1, gl%ncomp - 1
            daccenLKP_dxi(i) = gl%accen(i) - gl%accen(gl%ncomp)
            dzcLKP_dxi(i)=-0.085d0 * daccenLKP_dxi(i)
        end do

        B_0 = gl%lkp_b1_0 - gl%lkp_b2_0 * tau - gl%lkp_b3_0 * tau2 - gl%lkp_b4_0 * tau3
        C_0 = gl%lkp_c1_0 - gl%lkp_c2_0 * tau + gl%lkp_c3_0 * tau3
        D_0 = gl%lkp_d1_0 + gl%lkp_d2_0 * tau
        B_ref = gl%lkp_b1_ref - gl%lkp_b2_ref * tau - gl%lkp_b3_ref * tau2 - gl%lkp_b4_ref * tau3
        C_ref = gl%lkp_c1_ref - gl%lkp_c2_ref * tau + gl%lkp_c3_ref * tau3
        D_ref = gl%lkp_d1_ref + gl%lkp_d2_ref * tau

        do i = 1, gl%ncomp - 1
            j=i !Preliminary solution
            !1st partial derivative WRT xi
            !Ideal fluid contribution
            part_dar_dxi_0(i) = (- B_0 / zcLKP2 * del - C_0 / zcLKP3 * del2 - D_0 / zcLKP6 * del5) * dzcLKP_dxi(i)  &
                - gl%lkp_c4_0 * tau3 / (2.d0 * gl%lkp_gamma_0) * dzcLKP_dxi(i) * 2.d0 * del2 * gl%lkp_gamma_0 / zcLKP3 * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) &
                * (gl%lkp_gamma_0 / zcLKP2 *del2 + gl%lkp_beta_0)

            !Reference fluid contribution
            part_dar_dxi_ref(i) = (- B_ref / zcLKP2 * del - C_ref / zcLKP3 * del2 - D_ref / zcLKP6 * del5) * dzcLKP_dxi(i) &
                - gl%lkp_c4_ref * tau3 / (2.d0 * gl%lkp_gamma_ref) * dzcLKP_dxi(i) * 2.d0 * del2 * gl%lkp_gamma_ref / zcLKP3 * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) &
                * (gl%lkp_gamma_ref / zcLKP2 *del2 + gl%lkp_beta_ref)

            !2nd partial derivative WRT xi,xj
            !ideal fluid
            part_2dar_dxidxj_0(i,j) = (2.d0 * B_0 * del / zcLKP3 + 3.d0 * C_0 * del2 / zcLKP4 + 6.d0 * D_0 * del5 / zcLKP7) * dzcLKP_dxi(i) * dzcLKP_dxi(j) &
                - gl%lkp_c4_0 * tau3 / gl%lkp_gamma_0 * dzcLKP_dxi(i) * dzcLKP_dxi(j) * del2 * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) &
                * (2.d0 * (gl%lkp_gamma_0**3 / zcLKP8 * del4 + gl%lkp_gamma_0**2 / zcLKP6 * del2 * gl%LKP_beta_0) &
                -5.d0 * gl%lkp_gamma_0**2 / zcLKP6 * del2 - 3.d0 * gl%lkp_gamma_0 / zcLKP4 * gl%LKP_beta_0)

            !reference fluid
            part_2dar_dxidxj_ref(i,j) = (2.d0 * B_ref * del / zcLKP3 + 3.d0 * C_ref * del2 / zcLKP4 + 6.d0 * D_ref * del5 / zcLKP7) * dzcLKP_dxi(i) * dzcLKP_dxi(j) &
                - gl%lkp_c4_ref * tau3 / gl%lkp_gamma_ref * dzcLKP_dxi(i) * dzcLKP_dxi(j) * del2 * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) &
                * (2.d0 * (gl%lkp_gamma_ref**3 / zcLKP8 * del4 + gl%lkp_gamma_ref**2 / zcLKP6 * del2 * gl%LKP_beta_ref) &
                -5.d0 * gl%lkp_gamma_ref**2 / zcLKP6 * del2 - 3.d0 * gl%lkp_gamma_ref / zcLKP4 * gl%LKP_beta_ref)

            !calculate DERIVFXiXj(i,j)
            DERIVFX2(i) = - daccenLKP_dxi(i) / gl%lkp_w_ref * part_dar_dxi_0(j) - daccenLKP_dxi(j) / gl%lkp_w_ref * part_dar_dxi_0(i) + (1.d0 - gl%accenLKPMix / gl%lkp_w_ref) * part_2dar_dxidxj_0(i,j) &
                + daccenLKP_dxi(i) / gl%lkp_w_ref * part_dar_dxi_ref(j) + daccenLKP_dxi(j) / gl%lkp_w_ref * part_dar_dxi_ref(i) + gl%accenLKPMix / gl%lkp_w_ref * part_2dar_dxidxj_ref(i,j)

        end do

        !PC-SAFT
    else if (gl%mix_type == 6) then

        if(.not. allocated(x2)) allocate(x2(6,gl%ncomp,gl%ncomp))
        call ARX2DERIVS(gl,temperature, density, getder_res, x2)
        DERIVFX2_TRANS(1:gl%ncomp,1:gl%ncomp) = x2(1,:,:)
        !transformation to usual writing
        do j = 1, gl%ncomp
            do i = 1, gl%ncomp
                DERIVFX2_TRANS(i,j) = DERIVFX2_TRANS(i,j)- DERIVFX2_TRANS(gl%ncomp,j)- DERIVFX2_TRANS(gl%ncomp,i)+ DERIVFX2_TRANS(gl%ncomp,gl%ncomp)
            end do
        end do

        do i = 1, gl%ncomp
            DERIVFX2(i) = DERIVFX2_TRANS(i,i)
        end do

        !Quadratic mixing rules for the residual reduced Helmholtz energy, Andreas February 2015
    elseif (gl%mix_Type  ==  11) then !Lorentz Berthelot or modified mixing rules used Andreas March 2012

        !First, calculate all reduced residual Helmholtz energies of the pure fluids
        der_res = 0.D0
        do i = 1,gl%ncomp
            call FNRDERIVS(gl,temperature, density, getder_res, der_res, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_r_i(i) = der_res(1)    !Residual reduced Helmholtz energy of pure component i in the mixture
        end do

        alpha_r_NN = alpha_r_i(gl%ncomp)
        Do i = 1, gl%NCOMP-1
            alpha_r_ij = alpha_r_i(i)
            alpha_r_iN = 0.5D0 * (alpha_r_i(i) + alpha_r_i(gl%ncomp)) * (1.D0 - ACCESS_KIJ_HELM(i,gl%ncomp))
            alpha_r_jN = alpha_r_iN
            !Contribution of the quadratic mixing term for the pure reduced residual Helmholtz energies
            DERIVFX2(i) = DERIVFX2(i) + 2.D0 * (alpha_r_ij - alpha_r_jN - alpha_r_iN + alpha_r_NN)

            !Contribution of the departure function
            Call DEPFUNCFNR (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP)
            alphaiN_r = der_dep(1)                          !Get departure Function for Comp. i and N
            DERIVFX2(i) = DERIVFX2(i) + (-2.D0)*gl%Fij(i,gl%NCOMP)* alphaiN_r
        End do


        !Excess based departure function, Andreas Jäger July 2016
    else if (gl%mix_type  ==  12) then

        !Preliminary solution, because this routine is in principle not needed
        !For mix_type = 12 the diagonal (so d2ar/dxidxj for i=j) is already filled in the routine d2ar_dxidxj
        !which is simply called here
        !Andreas Jäger, July 2017
        call d2ar_dxidxj (gl,Temperature, Density, d2ardxidxj)
        do i = 1, gl%ncomp-1
            DERIVFX2(i) = d2ardxidxj(i,i)
        end do

        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !!Second numerical derivative of the excess based departure function with respect to xi^2 at constant tau and del
        !delta_x = 1.D-4
        !molfractions_orig = molfractions
        !tredmix_orig = tredmix
        !rhoredmix_orig = rhoredmix
        !
        !!Value at unchanged temperature, density and composition
        !call DEPFUNC_GE_BASED (Temperature, Density, getder_dep, der_dep)
        !alpha_dep = der_dep(1)
        !
        !do i =1, ncomp-1
        !
        !    !Increase xi
        !    molfractions(i) = molfractions_orig(i) + delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !    call reduced_parameters_calc(gl,Temperature)
        !    Temp_p = tredmix/tredmix_orig * Temperature
        !    Dens_p = rhoredmix/rhoredmix_orig * Density
        !    call DEPFUNC_GE_BASED (Temp_p, Dens_p, getder_dep, der_dep)
        !    alpha_dep_p = der_dep(1)
        !    !Decrease xi
        !    molfractions(i) = molfractions_orig(i) - delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !    call reduced_parameters_calc(gl,Temperature)
        !    Temp_m = tredmix/tredmix_orig * Temperature
        !    Dens_m = rhoredmix/rhoredmix_orig * Density
        !    call DEPFUNC_GE_BASED (Temp_m, Dens_m, getder_dep, der_dep)
        !    alpha_dep_m = der_dep(1)
        !    !Calculate the numerical derivative
        !    DERIVFX2_NUM(i) = (alpha_dep_p - 2.D0 * alpha_dep + alpha_dep_m) / delta_x**2
        !    molfractions = molfractions_orig
        !
        !end do
        !call reduced_parameters_calc(gl,Temperature)
        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !Erik, April 2018
    else if (gl%mix_type == 13) then
        !!!Second numerical derivative of the excess based departure function with respect to xi^2 at constant tau and del
        !delta_x = 1.D-4
        !molfractions_orig = molfractions
        !tredmix_orig = tredmix
        !rhoredmix_orig = rhoredmix
        !!
        !!!Value at unchanged temperature, density and composition
        !call gE_COSMO_SAC_MIXDERIVS(Temperature, getder_dep, gE_C, gE_R, ln_gamma_C, ln_gamma_R, C_or_R, errval)
        !!alpha_dep = der_dep(1)
        !!
        !do i =1, ncomp-1
        !!
        !    !Increase xi
        !    molfractions(i) = molfractions_orig(i) + delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !    call reduced_parameters_calc(Temperature)
        !    Temp_p = tredmix/tredmix_orig * Temperature
        !    Dens_p = rhoredmix/rhoredmix_orig * Density
        !!    call DEPFUNC_GE_BASED (Temp_p, Dens_p, getder_dep, der_dep)
        !    call gE_COSMO_SAC_MIXDERIVS(Temp_p, getder_dep, gE_C_p, gE_R_p, ln_gamma_C_p, ln_gamma_R_p, C_or_R, errval)
        !!    alpha_dep_p = der_dep(1)
        !    !Decrease xi
        !    molfractions(i) = molfractions_orig(i) - delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !    call reduced_parameters_calc(Temperature)
        !    Temp_m = tredmix/tredmix_orig * Temperature
        !    Dens_m = rhoredmix/rhoredmix_orig * Density
        !!    call DEPFUNC_GE_BASED (Temp_m, Dens_m, getder_dep, der_dep)
        !    call gE_COSMO_SAC_MIXDERIVS(Temp_m, getder_dep, gE_C_m, gE_R_m, ln_gamma_C_m, ln_gamma_R_m, C_or_R, errval)
        !!    alpha_dep_m = der_dep(1)
        !    !Calculate the numerical derivative
        !!    DERIVFX2_NUM(i) = (alpha_dep_p - 2.D0 * alpha_dep + alpha_dep_m) / delta_x**2
        !    gE_C_dxa(1,i) = (gE_C_p(1) - 2.D0 * gE_C(1) + gE_C_m(1)) / delta_x**2
        !    gE_R_dxa(1,i) = (gE_R_p(1) - 2.D0 * gE_R(1) + gE_R_m(1)) / delta_x**2
        !    molfractions = molfractions_orig
        !!
        !end do
        !call reduced_parameters_calc(Temperature)

        if ((gl.cosmo.COSMO_ver == 1) .and. (gl.cosmo.analytical)) then
            !analytical derivatives
            call d2ar_dxidxj (gl,Temperature, Density, d2ardxidxj)
            do i = 1, gl.ncomp-1
                DERIVFX2(i) = d2ardxidxj(i,i)
            end do
        else
            !Second numerical derivative of the excess based departure function with respect to xi^2 at constant tau and del
            delta_x = 1.D-4
            molfractions_orig = gl%molfractions
            tredmix_orig = gl%tredmix
            rhoredmix_orig = gl%rhoredmix

            !Value at unchanged temperature, density and composition
            call DEPFUNC_GE_BASED (gl, Temperature, Density, getder_dep, der_dep)
            alpha_dep = der_dep(1)

            do i =1, gl%ncomp-1

                !Increase xi
                gl%molfractions(i) = molfractions_orig(i) + delta_x
                gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) - delta_x
                call reduced_parameters_calc(gl, Temperature)
                Temp_p = gl%tredmix/tredmix_orig * Temperature
                Dens_p = gl%rhoredmix/rhoredmix_orig * Density
                call DEPFUNC_GE_BASED (gl, Temp_p, Dens_p, getder_dep, der_dep)
                alpha_dep_p = der_dep(1)
                !Decrease xi
                gl%molfractions(i) = molfractions_orig(i) - delta_x
                gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) + delta_x
                call reduced_parameters_calc(gl, Temperature)
                Temp_m = gl%tredmix/tredmix_orig * Temperature
                Dens_m = gl%rhoredmix/rhoredmix_orig * Density
                call DEPFUNC_GE_BASED (gl, Temp_m, Dens_m, getder_dep, der_dep)
                alpha_dep_m = der_dep(1)
                !Calculate the numerical derivative
                DERIVFX2(i) = (alpha_dep_p - 2.D0 * alpha_dep + alpha_dep_m) / delta_x**2
                gl%molfractions = molfractions_orig

            end do
            call reduced_parameters_calc(gl, Temperature)
        end if

        !-------------------------------------------------------------------------------
        !!numerical test of second derivative of gE_R
        !C_or_R = 2
        !do i =1, gl.ncomp-1
        !
        !    !Increase xi
        !    gl.molfractions(i) = molfractions_orig(i) + delta_x
        !    gl.molfractions(gl.ncomp) = molfractions_orig(gl.ncomp) - delta_x
        !    call reduced_parameters_calc(gl, Temperature)
        !    Temp_p = gl.tredmix/tredmix_orig * Temperature
        !    Dens_p = gl.rhoredmix/rhoredmix_orig * Density
        !    call gE_COSMO_SAC_MIXDERIVS_dxa(gl,Temp_p, getder_dep, C_or_R, errval)
        !    !alpha_dep_p = der_dep(1)
        !    gE_R_dxa_p = gl.gE_R_dxa(1,i)
        !    !Decrease xi
        !    gl.molfractions(i) = molfractions_orig(i) - delta_x
        !    gl.molfractions(gl.ncomp) = molfractions_orig(gl.ncomp) + delta_x
        !    call reduced_parameters_calc(gl, Temperature)
        !    Temp_m = gl.tredmix/tredmix_orig * Temperature
        !    Dens_m = gl.rhoredmix/rhoredmix_orig * Density
        !    call gE_COSMO_SAC_MIXDERIVS_dxa(gl,Temp_m, getder_dep, C_or_R, errval)
        !    !alpha_dep_m = der_dep(1)
        !    gE_R_dxa_m = gl.gE_R_dxa(1,i)
        !    !Calculate the numerical derivative
        !    !DERIVFX2(i) = (alpha_dep_p - 2.D0 * alpha_dep + alpha_dep_m) / delta_x**2
        !    gE_R_dxadxb(i) = (gE_R_dxa_p - gE_R_dxa_m) / (2.D0 *delta_x)
        !    gl.cosmo.gE_R_dxadxb_num(i) = gE_R_dxadxb(i)
        !    gl.molfractions = molfractions_orig
        !
        !end do
        !call reduced_parameters_calc(gl, Temperature)
        !------------------------------------------------------------------------


        !One-fluid model for multiparameter EOS, Andreas Jäger July 2016
    else if (gl%mix_type == 19) then

        !!Second numerical derivative of the non-corresponding states based model with respect to xi^2 at constant tau and del
        !delta_x = 1.D-4
        !molfractions_orig = molfractions
        !
        !!Value at unchanged temperature, density and composition
        !call MIXDERIVSFNR  (Temperature, Density, getder_dep, der_dep)
        !alpha_dep = der_dep(1)
        !
        !do i =1, ncomp-1
        !
        !    !Increase xi
        !    molfractions(i) = molfractions_orig(i) + delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !    call MIXDERIVSFNR  (Temperature, Density, getder_dep, der_dep)
        !    alpha_dep_p = der_dep(1)
        !    !Decrease xi
        !    molfractions(i) = molfractions_orig(i) - delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !    call MIXDERIVSFNR  (Temperature, Density, getder_dep, der_dep)
        !    alpha_dep_m = der_dep(1)
        !    !Calculate the numerical derivative
        !    DERIVFX2(i) = (alpha_dep_p - 2.D0 * alpha_dep + alpha_dep_m) / delta_x**2
        !    molfractions = molfractions_orig
        !
        !end do

        !Analytical derivative of the non-corresponding states based model with respect to xi and xj at constant tau and delta
        !First, calculate all reduced residual Helmholtz energies of the pure fluids
        der_res = 0.D0
        do i = 1,gl%ncomp
            gl%tredmix = gl%tc(i)
            gl%rhoredmix = gl%rhoc(i)
            call FNRDERIVS(gl,temperature, density, getder_res, der_res, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_r_i(i) = der_res(1)    !Residual reduced Helmholtz energy of pure component i in the mixture
            gl%tredmix = 1.D0
            gl%rhoredmix = 1.D0
        end do

        Do i = 1, (gl%NCOMP - 1)
            Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, i,gl%NCOMP)
            alphaiN_r = der_dep(1)
            DERIVFX2(i) = 2.D0 * alpha_r_i(i) + 2.D0 * alpha_r_i(gl%ncomp) - 4.D0 * alphaiN_r
        end do

    end if

    end subroutine d2ar_dxi2
    !**************************************************************************

    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    subroutine ndP_dni_TV (gl,TEMPERATURE, DENSITY, dPdni)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE PRESSURE WITH RESPECT TO AMOUNT OF SUBSTANCE
    ! AT CONSTANT TEMPERATURE AND VOLUME.
    ! THE DERIVATIVE IS GIVEN TIMES THE TOTAL NUMBER OF MOLES n. AS PUBLISHED BY
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 108 table 7.2)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOLE/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dPdni       - Vector(30) the derivatives for each component are stored in
    !--------------------------------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: dPdni

    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivs):: getder_res
    double precision:: Fr_d_delta, Rmix
    double precision, dimension(30):: d_rhoredmix_d_ni, d_tredmix_d_ni, Fr_d_delta_d_ni

    integer:: i


    getder_res = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.D0
    dPdni = 0.D0                                        !Initialize return variable
    Fr_d_delta = 0.D0                                       !Initialize

    Call R_mix_calc(gl,Rmix)

    Call MIXDERIVSFNR (gl,temperature, density, getder_res, der_res)
    Fr_d_delta = der_res(2)                                 !Derivative of the residual Helmholtz with respect to delta
    call ndYr_dni(gl,d_tredmix_d_ni, d_rhoredmix_d_ni)       !Get derivatives of reduced parameters with respect to ni
    Call nd2ar_dniddel (gl,TEMPERATURE, DENSITY, Fr_d_delta_d_ni) !Second Derivative of the residual Helmholtz with respect to delta and ni

    Do i = 1, gl%NCOMP
        dPdni(i) = Density * Temperature *  Rmix * (1.D0+Fr_d_delta*(2.D0-d_rhoredmix_d_ni(i)/gl%rhoredmix)+Fr_d_delta_d_ni(i))
    end do

    end subroutine ndP_dni_TV
    !**************************************************************************

    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine nd2ar_dniddel (gl,TEMPERATURE, DENSITY, Danidel)
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE RESIDUAL HELMHOLTZ ENERGY
    ! WITH RESPECT TO THE REDUCED DENSITY DELTA AND THE MOLAR NUMBER OF COMPONENT i IN THE
    ! MIXTURE. THE DERIVATIVE IS GIVEN TIMES THE TOTAL NUMBER OF MOLES n ADN DELTA. THE
    ! CALCULATION IS BASED ON THE FORMULATION OF THE HELMHOLTZ FREE ENERGY AS PUBLISHED BY:
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 123 EQ. (7.64))
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! Temperature
    ! Density
    !
    ! OUTPUT PARAMETERS:
    ! Danidel - n x 1 vector where n is the number of components in the mixture. The derivative
    ! n * delta * d(Fr)_d(ni)_d(del) is stored in this vector.
    !-------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: Danidel

    integer, dimension(nderivs):: getder_res
    double precision, dimension(nderivs):: der_res
    double precision, dimension(30):: dda_dxdd

    !derivatives of the reducingfunctions with respect to ni multiplied by n
    double precision, dimension(30):: drhoredmix_dni, dtredmix_dni
    !second derivatives of Fr with respect to delta, tau, x and combinations of those are needed
    double precision:: dda_dddd, dda_dtdd, dda_dtdxi, dda_dtdxj
    integer:: i, j

    getder_res = (/0,0,1,0,0,1,0,0,0,0,0,0,0,0,0/)                                !Derivatives with respect to del^2 and to del and tau are needed
    der_res = 0.D0                                                  !Initialize                                                    !Initialize
    dda_dxdd = 0.D0
    Danidel = 0.D0

    Call ndYr_dni(gl,dtredmix_dni, drhoredmix_dni)                         !Get derivatives of reduced parameters with respect to ni
    Call MIXDERIVSFNR (gl,Temperature, Density, getder_res, der_res)       !Second (mixed) derivatives of Fr with respect to delta and tau
    Call d2ar_dxiddel (gl,Temperature, Density, dda_dxdd)                  !Second Derivative of F_res with respect to xi and delta
    dda_dddd = der_res(3)
    dda_dtdd = der_res(6)

    Do i = 1, gl%NCOMP
        Danidel(i) = dda_dddd*(1.d0-1.d0/gl%rhoredmix*drhoredmix_dni(i))+ &
            & dda_dtdd/gl%tredmix*dtredmix_dni(i)
        if (i  <  gl%NCOMP) then
            dda_dtdxi = dda_dxdd(i)
            Danidel(i) = Danidel(i) + dda_dtdxi
        end if
        Do j = 1, (gl%NCOMP-1)
            dda_dtdxj = dda_dxdd(j)
            Danidel(i) = Danidel(i) - gl%molfractions(j) * dda_dtdxj
        end do
    end do

    end subroutine nd2ar_dniddel
    !**************************************************************************

    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    subroutine d2ar_dxidxj (gl,TEMPERATURE, DENSITY, DERIVFXiXj)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE MIXED SECOND DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE MOLEFRACTION Xi and Xj FOR EACH COMPONENT.
    ! ACCORDING TO THE ALGORITHM BY KEN. HALL, THE MOLEFRACTION OF THE LAST
    ! COMPONENT Xn IS REPLACED BY 1 - SUM OF ALL OTHER MOLEFRACTIONS.
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A MATRIX (30,30)IN WHICH THE MIXED DERIVATIVE WITH RESPECT TO MOLFRACTION i and j FOR
    !              EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------









    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30,30):: DERIVFXiXj, DERIVFXiXj_transform

    double precision, dimension(nderivs):: der_dep
    double precision, dimension(nderivs) :: der_res
    integer, dimension(nderivs):: getder_dep
    integer, dimension(nderivs) :: getder_res
    integer:: i, j, k
    double precision:: Fij_r, FiN_r, FjN_r, delr

    !Variable for SRK and PR
    double precision, dimension(30,30):: d2ardxidxj
    !SH variable is part of gl
    !double precision,allocatable :: MIXDERFNR_dxidxj(:,:,:)

    !Variables for LKP
    double precision :: B_0, C_0, D_0, B_ref, C_ref, D_ref
    double precision, dimension(30) :: daccenLKP_dxi, dzcLKP_dxi
    double precision :: zcLKP2, zcLKP3, zcLKP4, zcLKP5,zcLKP6, zcLKP7, zcLKP8, tau, tau2, tau3, del, del2, del4, del5

    double precision, dimension(30) :: part_dar_dxi_ref, part_dar_dxi_0, help_0, help_ref
    double precision, dimension(:,:), allocatable :: part_2dar_dxidxj_0, part_2dar_dxidxj_ref

    !New variables for quadratic mixing rules for the residual reduced Helmholtz energy
    !Andreas, February 2016
    double precision:: alpha_r_ij            !Combination term alpha^r_(ij) = 0.5 * (alpha^r_i + alpha^r_j) (1-k_(ij))
    double precision:: alpha_r_jN, alpha_r_iN, alpha_r_NN
    double precision, dimension(30):: alpha_r_i         !reduced residual Helmholtz energies of the pure fluids

    !New variables for numerical derivative of the excess based departure function
    !Andreas Jäger, July 2016
    double precision:: delta_x                          !Step in x for numerical derivative
    double precision:: alpha_dep_ipjp                   !Value of departure function at increased xi and increased xj (and decreased xN)
    double precision:: alpha_dep_ipjm                   !Value of departure function at increased xi and decreased xj
    double precision:: alpha_dep_imjp                   !Value of departure function at decreased xi and increased xj
    double precision:: alpha_dep_imjm                   !Value of departure function at decreased xi and decreased xj (and increased xN)
    double precision, dimension(30):: molfractions_orig !Original value of molfractions
    double precision:: Temp_pm, dens_pm                 !Decrease and increased temperature and density. When xi and xj are changed and tau and delta stay constant, T and rho change
    double precision:: tredmix_orig, rhoredmix_orig
    double precision:: rhomix_ref_copy
    double precision, dimension(30):: rho_i_ref_copy
    !double precision, dimension(30,30):: DERIVFXiXj_num
    !Andreas Jäger, June 2017
    double precision:: const_A1
    double precision:: help_h, help_f
    double precision, dimension(30):: dhelph_dxa, dhelpf_dxa
    double precision, dimension(:,:), allocatable :: d2helpf_dxadxb, d2helph_dxadxb
    double precision, dimension(30):: alpha_oi_r_mix, alpha_oi_r_i
    integer, dimension(nderivs):: GETDER_i                         ! array specifier to indicate, which derivative is needed for the pure fluid residual Helmholtz energies
    double precision, dimension(nderivs)::FNRDER_i                 ! array with the computed values for the derivatives of the residual Helmholtz energy of the pure fluids
    double precision:: aE
    double precision, dimension(30):: dalpha_oi_r_mix_ddelref, dalpha_oi_r_i_dtaui
    double precision, dimension(30):: d2alpha_oi_r_mix_ddelref2, d2alpha_oi_r_i_dtaui2
    double precision, dimension(30):: ddelref_dxa
    double precision, dimension(:,:), allocatable :: dtaui_dxa
    double precision, dimension(:,:), allocatable :: d2delref_dxadxb
    double precision, allocatable :: d2taui_dxadxb(:,:,:)
    double precision, dimension(30):: dTred_dxi, drhored_dxi, d2Tred_dxi2, d2rhored_dxi2
    double precision, dimension(:,:), allocatable :: d2Tred_dxidxj, d2rhored_dxidxj
    double precision, dimension(30):: dbHelmgE_dxi
    double precision, dimension(:,:), allocatable :: d2bHelmgE_dxidxj
    double precision, dimension(30):: daE_dxi
    double precision, dimension(:,:), allocatable :: d2aE_dxidxj
    !Variables for UNIFAC
    !double precision, dimension(nderivs) :: gl%ge%gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gl%ge%gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau
    integer:: C_or_R
    integer:: errval
    !double precision, dimension(nderivs, 30) :: gl%ge%gE_C_dxa        !Combinatorial part of gE and derivatives with respect to tau and delta and xa
    !double precision, dimension(nderivs, 30) :: gl%ge%gE_R_dxa        !Residual part of gE and derivatives with respect to tau and delta and xa
    !double precision, allocatable :: gl%ge%ln_gamma_C_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to delta and tau and xa
    !double precision, allocatable :: gl%ge%ln_gamma_R_dxa(:,:,:)       !Residual activity coefficients and derivatives with respect to delta and tau and xa
    !double precision, allocatable :: gl%ge%gE_C_dxadxb(:,:,:)                 !Combinatorial part of gE and derivatives with respect to tau and delta and xa and xb
    !double precision, allocatable :: gl%ge%gE_R_dxadxb(:,:,:)                 !Residual part of gE and derivatives with respect to tau and delta and xa and xb
    !double precision, allocatable :: gl%ge%ln_gamma_C_dxadxb(:,:,:,:)   !Combinatorial activity coefficients and derivatives with respect to delta and tau and xa and xb
    !double precision, allocatable :: gl%ge%ln_gamma_R_dxadxb(:,:,:,:)           !Residual activity coefficients and derivatives with respect to delta and tau and xa and xb


    !New variables for PC-SAFT equation type
    double precision, allocatable :: x2(:,:,:)

    !New variables for COSMO-SAC
    double precision :: Temp_p, Temp_m
    double precision, dimension(nderivs) :: gE_C_p, gE_C_m, gE_R_p, gE_R_m, gE_C_ipjp, gE_R_ipjp, gE_C_imjm, gE_R_imjm
    double precision, dimension(:, :), allocatable :: ln_gamma_C_p, ln_gamma_C_m, ln_gamma_R_p, ln_gamma_R_m, ln_gamma_C_ipjp, ln_gamma_R_ipjp, ln_gamma_C_imjm, ln_gamma_R_imjm
    double precision, dimension(nderivs) :: gE_C_ipjm, gE_R_ipjm, gE_C_imjp, gE_R_imjp
    double precision, dimension(:, :), allocatable :: ln_gamma_C_ipjm, ln_gamma_R_ipjm, ln_gamma_C_imjp, ln_gamma_R_imjp


    if(.not. allocated(gl%MIXDERIVFNR_dxidxj)) allocate(gl%MIXDERIVFNR_dxidxj(15,30,30))
    if (.not. allocated(part_2dar_dxidxj_0)) then
        allocate(part_2dar_dxidxj_0(30,30))
        allocate(part_2dar_dxidxj_ref, d2helpf_dxadxb, d2helph_dxadxb, dtaui_dxa, d2delref_dxadxb, d2Tred_dxidxj, d2rhored_dxidxj, &
        & d2bHelmgE_dxidxj, d2aE_dxidxj, mold=part_2dar_dxidxj_0)
    endif

    getder_res = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    getder_dep = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.D0
    der_dep = 0.D0
    DERIVFXiXj = 0.D0              !Initialize return matrix DERIVFXiXj
    Fij_r = 0.D0                        !Initialize departure function
    DERIVFXiXj_transform = 0.d0

    if (gl%mix_Type  ==  1) then !Lorentz Berthelot or modified mixing rules used Andreas March 2012
        Do i = 1, (gl%NCOMP - 1)
            Call DEPFUNCFNR (gl,temperature, density, getder_dep, der_dep, i,gl%NCOMP)
            FiN_r = der_dep(1)
            Do j = (i+1), (gl%NCOMP-1)
                if((gl%Fij(i,j)  /=  0.d0).OR.(gl%Fij(i,gl%NCOMP)  /=  0.d0).OR.(gl%Fij(j,gl%NCOMP)  /=  0.d0)) then
                    Call DEPFUNCFNR (gl,temperature, density, getder_dep, der_dep, i,j)
                    Fij_r = der_dep(1)
                    Call DEPFUNCFNR (gl,temperature, density, getder_dep, der_dep, j,gl%NCOMP)
                    FjN_r = der_dep(1)
                    DERIVFXiXj(i,j) = gl%Fij(i,j) * Fij_r -  gl%Fij(i,gl%NCOMP) * FiN_r - gl%Fij(j,gl%NCOMP) * FjN_r
                    DERIVFXiXj(j,i) = DERIVFXiXj(i,j)
                end if
            end do
        end do
    else if ((gl%mix_Type  ==  2) .or. (gl%mix_type == 21) .or. (gl%mix_type == 22)) then !SRK with SRK mixing rules used
        !Replaced with new cubics routine
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_dxidxj_CUBIC (gl,temperature, density, getder_res)
        DERIVFXiXj(:,:) = gl%MIXDERIVFNR_dxidxj(1,:,:)
        !Old version
        !call d2ar_dxidxj_SRK(gl,temperature, density, d2ardxidxj)
        !DERIVFXiXj = d2ardxidxj

    else if ((gl%mix_Type  ==  3) .or. (gl%mix_type == 31)) then !PR with PR mixing rules used
        !Replaced with new cubics routine
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_dxidxj_CUBIC (gl,temperature, density, getder_res)
        DERIVFXiXj(:,:) = gl%MIXDERIVFNR_dxidxj(1,:,:)
        !Old version
        !call d2ar_dxidxj_PR(temperature, density, d2ardxidxj)
        !DERIVFXiXj = d2ardxidxj

    else if (gl%mix_Type  ==  4) then !LKP with LKP mixing rules used
        !edited by Monika 10/2014

        tau = gl%tredmix/temperature
        del = density/gl%rhoredmix

        del2=del*del
        del4=del2*del2
        del5=del4*del
        zcLKP2=gl%zcLKP*gl%zcLKP
        zcLKP3=zcLKP2*gl%zcLKP
        zcLKP4=zcLKP2*zcLKP2
        zcLKP5=gl%zcLKP*zcLKP4
        zcLKP6=gl%zcLKP*zcLKP5
        zcLKP7=zcLKP6*gl%zcLKP
        zcLKP8=zcLKP7*gl%zcLKP
        tau2=tau*tau
        tau3=tau*tau2


        do i = 1, gl%ncomp - 1
            daccenLKP_dxi(i) = gl%accen(i) - gl%accen(gl%ncomp)
            dzcLKP_dxi(i)=-0.085d0 * daccenLKP_dxi(i)
        end do

        B_0 = gl%lkp_b1_0 - gl%lkp_b2_0 * tau - gl%lkp_b3_0 * tau2 - gl%lkp_b4_0 * tau3
        C_0 = gl%lkp_c1_0 - gl%lkp_c2_0 * tau + gl%lkp_c3_0 * tau3
        D_0 = gl%lkp_d1_0 + gl%lkp_d2_0 * tau
        B_ref = gl%lkp_b1_ref - gl%lkp_b2_ref * tau - gl%lkp_b3_ref * tau2 - gl%lkp_b4_ref * tau3
        C_ref = gl%lkp_c1_ref - gl%lkp_c2_ref * tau + gl%lkp_c3_ref * tau3
        D_ref = gl%lkp_d1_ref + gl%lkp_d2_ref * tau

        do i = 1, gl%ncomp - 1
            do j = 1, gl%ncomp - 1
                !1st partial derivative WRT xi
                !Ideal fluid contribution
                part_dar_dxi_0(i) = (- B_0 / zcLKP2 * del - C_0 / zcLKP3 * del2 - D_0 / zcLKP6 * del5) * dzcLKP_dxi(i)  &
                    - gl%lkp_c4_0 * tau3 / (2.d0 * gl%lkp_gamma_0) * dzcLKP_dxi(i) * 2.d0 * del2 * gl%lkp_gamma_0 / zcLKP3 * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) &
                    * (gl%lkp_gamma_0 / zcLKP2 *del2 + gl%lkp_beta_0)

                !Reference fluid contribution
                part_dar_dxi_ref(i) = (- B_ref / zcLKP2 * del - C_ref / zcLKP3 * del2 - D_ref / zcLKP6 * del5) * dzcLKP_dxi(i) &
                    - gl%lkp_c4_ref * tau3 / (2.d0 * gl%lkp_gamma_ref) * dzcLKP_dxi(i) * 2.d0 * del2 * gl%lkp_gamma_ref / zcLKP3 * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) &
                    * (gl%lkp_gamma_ref / zcLKP2 *del2 + gl%lkp_beta_ref)

                !2nd partial derivative WRT xi,xj
                !ideal fluid
                part_2dar_dxidxj_0(i,j) = (2.d0 * B_0 * del / zcLKP3 + 3.d0 * C_0 * del2 / zcLKP4 + 6.d0 * D_0 * del5 / zcLKP7) * dzcLKP_dxi(i) * dzcLKP_dxi(j) &
                    - gl%lkp_c4_0 * tau3 / gl%lkp_gamma_0 * dzcLKP_dxi(i) * dzcLKP_dxi(j) * del2 * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) &
                    * (2.d0 * (gl%lkp_gamma_0**3 / zcLKP8 * del4 + gl%lkp_gamma_0**2 / zcLKP6 * del2 * gl%LKP_beta_0) &
                    -5.d0 * gl%lkp_gamma_0**2 / zcLKP6 * del2 - 3.d0 * gl%lkp_gamma_0 / zcLKP4 * gl%LKP_beta_0)

                !reference fluid
                part_2dar_dxidxj_ref(i,j) = (2.d0 * B_ref * del / zcLKP3 + 3.d0 * C_ref * del2 / zcLKP4 + 6.d0 * D_ref * del5 / zcLKP7) * dzcLKP_dxi(i) * dzcLKP_dxi(j) &
                    - gl%lkp_c4_ref * tau3 / gl%lkp_gamma_ref * dzcLKP_dxi(i) * dzcLKP_dxi(j) * del2 * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) &
                    * (2.d0 * (gl%lkp_gamma_ref**3 / zcLKP8 * del4 + gl%lkp_gamma_ref**2 / zcLKP6 * del2 * gl%LKP_beta_ref) &
                    -5.d0 * gl%lkp_gamma_ref**2 / zcLKP6 * del2 - 3.d0 * gl%lkp_gamma_ref / zcLKP4 * gl%LKP_beta_ref)

                !calculate DERIVFXiXj(i,j)
                DERIVFXiXj(i,j) = - daccenLKP_dxi(i) / gl%lkp_w_ref * part_dar_dxi_0(j) - daccenLKP_dxi(j) / gl%lkp_w_ref * part_dar_dxi_0(i) + (1.d0 - gl%accenLKPMix / gl%lkp_w_ref) * part_2dar_dxidxj_0(i,j) &
                    + daccenLKP_dxi(i) / gl%lkp_w_ref * part_dar_dxi_ref(j) + daccenLKP_dxi(j) / gl%lkp_w_ref * part_dar_dxi_ref(i) + gl%accenLKPMix / gl%lkp_w_ref * part_2dar_dxidxj_ref(i,j)

            end do
        end do

        !PC-SAFT
    else if (gl%mix_type == 6) then

        if(.not. allocated(x2)) allocate(x2(6,gl%ncomp,gl%ncomp))
        call ARX2DERIVS(gl,temperature, density, getder_res, x2)
        DERIVFXiXj(1:gl%ncomp,1:gl%ncomp) = x2(1,:,:)
        !transformation to usual writing
        do j = 1, gl%ncomp
            do i = 1, gl%ncomp
                DERIVFXiXj_transform(i,j) = DERIVFXiXj(i,j)- DERIVFXiXj(gl%ncomp,j)- DERIVFXiXj(gl%ncomp,i)+ DERIVFXiXj(gl%ncomp,gl%ncomp)
            end do
        end do

        DERIVFXiXj = DERIVFXiXj_transform

        do i = 1, gl%ncomp
            DERIVFXiXj(i,i) = 0.d0
        end do


        !Quadratic mixing rules for the residual reduced Helmholtz energy, Andreas February 2015
    elseif (gl%mix_Type  ==  11) then !Lorentz Berthelot or modified mixing rules used Andreas March 2012

        !First, calculate all reduced residual Helmholtz energies of the pure fluids
        der_res = 0.D0
        do i = 1,gl%ncomp
            call FNRDERIVS(gl,temperature, density, getder_res, der_res, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_r_i(i) = der_res(1)    !Residual reduced Helmholtz energy of pure component i in the mixture
        end do


        !Contribution of the quadratic mixing term for the pure reduced residual Helmholtz energies
        alpha_r_NN = alpha_r_i(gl%ncomp)
        Do i = 1, gl%NCOMP-1
            alpha_r_iN = 0.5D0 * (alpha_r_i(i) + alpha_r_i(gl%ncomp)) * (1.D0 - ACCESS_KIJ_HELM(i,gl%ncomp))
            Do j = 1,gl%NCOMP-1
                if (i == j) then
                    alpha_r_ij = alpha_r_i(i)
                end if
                if (i /= j) then
                    alpha_r_ij = 0.5D0 * (alpha_r_i(i) + alpha_r_i(j)) * (1.D0 - ACCESS_KIJ_HELM(i,j))
                end if
                alpha_r_jN = 0.5D0 * (alpha_r_i(j) + alpha_r_i(gl%ncomp)) * (1.D0 - ACCESS_KIJ_HELM(j,gl%ncomp))
                DERIVFXiXj(i,j) = DERIVFXiXj(i,j) + 2.D0 * (alpha_r_ij - alpha_r_jN - alpha_r_iN + alpha_r_NN)
            End do
        End do

        !Contribution of the departure function
        Do i = 1, (gl%NCOMP - 1)
            Call DEPFUNCFNR (gl,temperature, density, getder_dep, der_dep, i,gl%NCOMP)
            FiN_r = der_dep(1)
            Do j = (i+1), (gl%NCOMP-1)
                !if((Fij(i,j)  /=  0.d0).OR.(Fij(i,NCOMP)  /=  0.d0).OR.(Fij(j,NCOMP)  /=  0.d0)) then
                Call DEPFUNCFNR (gl,temperature, density, getder_dep, der_dep, i,j)
                Fij_r = der_dep(1)
                Call DEPFUNCFNR (gl,temperature, density, getder_dep, der_dep, j,gl%NCOMP)
                FjN_r = der_dep(1)
                DERIVFXiXj(i,j) = DERIVFXiXj(i,j) + gl%Fij(i,j) * Fij_r -  gl%Fij(i,gl%NCOMP) * FiN_r - gl%Fij(j,gl%NCOMP) * FjN_r
                DERIVFXiXj(j,i) = DERIVFXiXj(i,j)
                !end if
            end do
        end do

        !Excess based departure function, Andreas Jäger July 2016
    else if ((gl%mix_type  ==  12) .or. (gl%mix_type  ==  13)) then
        if(.not. allocated(gl%ge%ln_gamma_C_dxa)) allocate(gl%ge%ln_gamma_C_dxa(nderivs, 30, 30))
        if(.not. allocated(gl%ge%ln_gamma_R_dxa)) allocate(gl%ge%ln_gamma_R_dxa(nderivs, 30, 30))
        if(.not. allocated(gl%ge%gE_C_dxadxb)) allocate(gl%ge%gE_C_dxadxb(nderivs, 30, 30))
        if(.not. allocated(gl%ge%gE_R_dxadxb)) allocate(gl%ge%gE_R_dxadxb(nderivs, 30, 30))
        if(.not. allocated(gl%ge%ln_gamma_C_dxadxb)) allocate(gl%ge%ln_gamma_C_dxadxb(nderivs, 30, 30, 30))
        if(.not. allocated(gl%ge%ln_gamma_R_dxadxb)) allocate(gl%ge%ln_gamma_R_dxadxb(nderivs, 30, 30, 30))
        if(.not. allocated(d2taui_dxadxb)) allocate(d2taui_dxadxb(30,30,30))

        tredmix_orig = gl%tredmix
        rhoredmix_orig = gl%rhoredmix
        tau = gl%tredmix / temperature
        del = density / gl%rhoredmix
        !Calculate analytical derivative of the excess based departure function with respect to xi and xj
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !Compute auxiliary function f(del,x)
        !----
        const_A1 = -dlog(1.D0 / gl%u_pack + 1.D0)
        help_f = -dlog(gl%b_HelmgE * density + 1.D0) / const_A1
        !----

        !Compute derivative of auxiliary function f(del,x) with respect to xa
        !----
        call dYr_dxi(gl,dTred_dxi, drhored_dxi)
        call db_HelmhE_dxi(gl,dbHelmgE_dxi)
        do j = 1, gl%ncomp-1
            dhelpf_dxa(j) = -del * (gl%rhoredmix * dbHelmgE_dxi(j) + drhored_dxi(j) * gl%b_HelmgE) &
                & / const_A1 / (1.D0 + del * gl%b_HelmgE * gl%rhoredmix)
        end do
        !----

        !Compute the second derivative of the auxiliary function f(del,x) with respect to xa and xb
        !----
        call d2Yr_dxi2(gl,d2Tred_dxi2, d2rhored_dxi2)
        call d2Yr_dxidxj(gl,d2Tred_dxidxj, d2rhored_dxidxj)
        do i = 1, gl%ncomp-1
            d2Tred_dxidxj(i,i) = d2Tred_dxi2(i)
            d2rhored_dxidxj(i,i) = d2rhored_dxi2(i)
        end do
        call d2b_HelmhE_dxidxj(gl,d2bHelmgE_dxidxj)
        do i = 1, gl%ncomp-1
            do j = 1, gl%ncomp-1
                d2helpf_dxadxb(i,j) = -1.D0 / const_A1 * ((1.D0 + del * gl%b_HelmgE * gl%rhoredmix) * del &
                    & * (drhored_dxi(j) * dbHelmgE_dxi(i) + gl%rhoredmix * d2bHelmgE_dxidxj(i,j) &
                    & + d2rhored_dxidxj(i,j) * gl%b_HelmgE + drhored_dxi(i) * dbHelmgE_dxi(j)) &
                    & - del**2 * (gl%rhoredmix * dbHelmgE_dxi(i) + drhored_dxi(i) * gl%b_HelmgE) &
                    & * (gl%rhoredmix * dbHelmgE_dxi(j) + drhored_dxi(j) * gl%b_HelmgE)) &
                    &  / (1.D0 + del * gl%b_HelmgE * gl%rhoredmix)**2
            end do
        end do
        !----




        !Compute auxiliary function h(tau,x)
        !----
        GETDER_i = 0
        GETDER_i(1) = 1 !The Helmholtz energy itself is needed
        GETDER_i(2) = 1 !Get the first derivative of the Helmholtz energy with respect to deltaref
        GETDER_i(3) = 1 !Get the second derivative of the Helmholtz energy with respect to deltaref^2
        do i = 1, gl%ncomp
            !Get the reduced Helmholtz energies of the pure components at mixture conditions
            rhomix_ref_copy = gl%rho_mix_ref
            call FNRDERIVS(gl,temperature, rhomix_ref_copy, GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_oi_r_mix(i) = FNRDER_i(1)
            dalpha_oi_r_mix_ddelref(i) = FNRDER_i(2) / (gl%rho_mix_ref/gl%rhoredmix)
            d2alpha_oi_r_mix_ddelref2(i) = FNRDER_i(3) / (gl%rho_mix_ref/gl%rhoredmix)**2
        end do
        GETDER_i(2) = 0
        GETDER_i(3) = 0
        GETDER_i(4) = 1 !Get the first derivative of the Helmholtz energy with respect to taui of the pure fluid part
        GETDER_i(5) = 1 !Get the second derivative of the Helmholtz energy with respect to taui^2 of the pure fluid part
        do i = 1, gl%ncomp
            !Get the reduced Helmholtz energies of the pure components at pure fluid conditions
            gl%tredmix = gl%tc(i)         !Needed to evaluate the pure fluid at correct tau (tau_i)
            gl%rhoredmix = gl%rhoc(i)     !Needed to evaluate the pure fluid at correct delta (delta_i_ref)
            rho_i_ref_copy(i) = gl%rho_i_ref(i)
            call FNRDERIVS(gl,temperature, rho_i_ref_copy(i), GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_oi_r_i(i) = FNRDER_i(1)
            dalpha_oi_r_i_dtaui(i) = FNRDER_i(4) / (gl%tc(i)/temperature)
            d2alpha_oi_r_i_dtaui2(i) = FNRDER_i(5) / (gl%tc(i)/temperature)**2
        end do
        gl%tredmix = tredmix_orig
        gl%rhoredmix = rhoredmix_orig

        C_or_R = 2  !Get residual part of UNIFAC
        !Erik, April 2018
        if (gl%mix_Type == 12) then
            call gE_UNIFAC_MIXDERIVS(gl,temperature, getder_dep, C_or_R, errval)
        elseif (gl%mix_Type == 13) then
            call gE_COSMO_SAC_MIXDERIVS(gl,temperature, getder_dep, C_or_R, errval)
        end if
        aE = gl%ge%gE_C(1) + gl%ge%gE_R(1)    !Specify which parts of gE are considered (combinatorial, residual, or both)

        help_h = aE / R_HelmgE / temperature
        do i = 1, gl%ncomp
            !Andreas Jäger, December 2017
            !term "sum(xi ln(bi/b))" deleted from the excess based departure function, because it is assumed that it cancels out with the combinatorial part of gE
            !help_h = help_h - molfractions(i) * dlog(bi_HelmgE(i) / b_HelmgE) - molfractions(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
            help_h = help_h - gl%molfractions(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
        end do
        !----

        !Compute derivative of auxiliary function h(tau,x) with respect to xa
        !----
        !Erik ,April 2018
        if (gl%mix_Type == 12) then
            call gE_UNIFAC_MIXDERIVS_dxa(gl,temperature, getder_dep, C_or_R, errval)
        elseif (gl%mix_Type == 13) then
            if (.not.allocated(ln_gamma_C_p)) then
                allocate(ln_gamma_C_p(nderivs, 30))
                allocate(ln_gamma_C_m, ln_gamma_R_p, ln_gamma_R_m, ln_gamma_C_ipjp, ln_gamma_R_ipjp, ln_gamma_C_imjm, ln_gamma_R_imjm, &
                    ln_gamma_C_ipjm, ln_gamma_R_ipjm, ln_gamma_C_imjp, ln_gamma_R_imjp, mold = ln_gamma_C_p)
            end if
            !!Analytical derivative of the excess based departure function with respect to xi at constant tau and del
            if ((gl.cosmo.COSMO_ver == 1) .and. (gl.cosmo.analytical)) then
                call gE_COSMO_SAC_MIXDERIVS_dxa(gl,temperature, getder_dep, C_or_R, errval)
            else
                !!Numerical derivative of the excess based departure function with respect to xi at constant tau and del
                delta_x = 1.D-4
                molfractions_orig = gl%molfractions
                tredmix_orig = gl%tredmix
                rhoredmix_orig = gl%rhoredmix
                do i =1, gl%ncomp-1
                    !Increase xi
                    gl%molfractions(i) = molfractions_orig(i) + delta_x
                    gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) - delta_x
                    call reduced_parameters_calc(gl, Temperature)
                    Temp_p = gl%tredmix/tredmix_orig * Temperature
                    !Dens_p = rhoredmix/rhoredmix_orig * Density
                    call gE_COSMO_SAC_MIXDERIVS(gl, Temp_p, getder_dep, C_or_R, errval)
                    gE_C_p = gl%ge%gE_C
                    gE_R_p = gl%ge%gE_R
                    ln_gamma_C_p = gl%ge%ln_gamma_C
                    ln_gamma_R_p = gl%ge%ln_gamma_R
                    !    alpha_dep_p = der_dep(1)
                    !Decrease xi
                    gl%molfractions(i) = molfractions_orig(i) - delta_x
                    gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) + delta_x
                    call reduced_parameters_calc(gl, Temperature)
                    Temp_m = gl%tredmix/tredmix_orig * Temperature
                    !    Dens_m = rhoredmix/rhoredmix_orig * Density
                    call gE_COSMO_SAC_MIXDERIVS(gl, Temp_m, getder_dep, C_or_R, errval)
                    gE_C_m = gl%ge%gE_C
                    gE_R_m = gl%ge%gE_R
                    ln_gamma_C_m = gl%ge%ln_gamma_C
                    ln_gamma_R_m = gl%ge%ln_gamma_R
                    !    alpha_dep_m = der_dep(1)
                    !Calculate the numerical derivative
                    gl%ge%gE_C_dxa(1,i) = (gE_C_p(1) - gE_C_m(1)) / (2.D0 * delta_x)
                    gl%ge%gE_R_dxa(1,i) = (gE_R_p(1) - gE_R_m(1)) / (2.D0 * delta_x)
                    gl%molfractions = molfractions_orig
                end do
                call reduced_parameters_calc(gl, Temperature)
            end if
        end if
        daE_dxi = gl%ge%gE_C_dxa(1,:) + gl%ge%gE_R_dxa(1,:)
        do j = 1, gl%ncomp-1
            !Andreas Jäger, December 2017
            !term "sum(xi ln(bi/b))" deleted from the excess based departure function, because it is assumed that it cancels out with the combinatorial part of gE
            !dhelph_dxa(j) = tau / R_HelmgE / tredmix**2 * (tredmix * daE_dxi(j) - aE * dTred_dxi(j)) &
            !            & + dlog(bi_HelmgE(ncomp)/bi_HelmgE(j)) + dbHelmgE_dxi(j) / b_HelmgE &
            !            & - (alpha_oi_r_mix(j) - alpha_oi_r_i(j)) + (alpha_oi_r_mix(ncomp) - alpha_oi_r_i(ncomp))
            dhelph_dxa(j) = tau / R_HelmgE / gl%tredmix**2 * (gl%tredmix * daE_dxi(j) - aE * dTred_dxi(j)) &
                & - (alpha_oi_r_mix(j) - alpha_oi_r_i(j)) + (alpha_oi_r_mix(gl%ncomp) - alpha_oi_r_i(gl%ncomp))
            ddelref_dxa(j) = - 1.D0 / gl%u_pack * (drhored_dxi(j) / gl%b_HelmgE / gl%rhoredmix**2 + dbHelmgE_dxi(j) / gl%b_HelmgE**2 / gl%rhoredmix)
            do i = 1, gl%ncomp
                dtaui_dxa(i,j) =  - tau * gl%tc(i) / gl%tredmix**2 * dTred_dxi(j)
                dhelph_dxa(j) = dhelph_dxa(j) - gl%molfractions(i) * dalpha_oi_r_mix_ddelref(i) * ddelref_dxa(j) &
                    & + gl%molfractions(i) * dalpha_oi_r_i_dtaui(i) * dtaui_dxa(i,j)
            end do
        end do
        !----

        !Compute the second derivative of the auxiliary function h(tau,x) with respect to xa and xb
        !----
        !Erik ,April 2018
        if (gl%mix_Type == 12) then
            call gE_UNIFAC_MIXDERIVS_dxadxb(gl,temperature, getder_dep, C_or_R, errval)
        elseif (gl%mix_Type == 13) then
            if(allocated(gl%ge%gE_C_dxadxb))gl%ge%gE_C_dxadxb = 0.D0
            !!Analytical derivative of the excess based departure function with respect to xi at constant tau and del
            if ((gl.cosmo.COSMO_ver == 1) .and. (gl.cosmo.analytical)) then
                call gE_COSMO_SAC_MIXDERIVS_dxadxb(gl,temperature, getder_dep, C_or_R, errval)
            else
                !!Second numerical derivative of the excess based departure function with respect to xi and xj at constant tau and del
                delta_x = 1.D-4
                molfractions_orig = gl%molfractions
                tredmix_orig = gl%tredmix
                rhoredmix_orig = gl%rhoredmix
                !
                Do i = 1, gl%ncomp-1
                    do j = i+1, gl%ncomp-1
                        !Increase xi, increase xj
                        gl%molfractions(i) = molfractions_orig(i) + delta_x
                        gl%molfractions(j) = molfractions_orig(j) + delta_x
                        gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) - 2.D0 * delta_x
                        call reduced_parameters_calc(gl, Temperature)
                        Temp_pm = gl%tredmix/tredmix_orig * Temperature
                        Dens_pm = gl%rhoredmix/rhoredmix_orig * Density
                        call gE_COSMO_SAC_MIXDERIVS(gl, Temp_pm, getder_dep, C_or_R, errval)
                        gE_C_ipjp = gl%ge%gE_C
                        gE_R_ipjp = gl%ge%gE_R
                        ln_gamma_C_ipjp = gl%ge%ln_gamma_C
                        ln_gamma_R_ipjp = gl%ge%ln_gamma_R
                        !        alpha_dep_ipjp = der_dep(1)
                        !
                        !Increase xi, decrease xj
                        gl%molfractions(i) = molfractions_orig(i) + delta_x
                        gl%molfractions(j) = molfractions_orig(j) - delta_x
                        gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp)
                        call reduced_parameters_calc(gl, Temperature)
                        Temp_pm = gl%tredmix/tredmix_orig * Temperature
                        Dens_pm = gl%rhoredmix/rhoredmix_orig * Density
                        call gE_COSMO_SAC_MIXDERIVS(gl, Temp_pm, getder_dep, C_or_R, errval)
                        gE_C_ipjm = gl%ge%gE_C
                        gE_R_ipjm = gl%ge%gE_R
                        ln_gamma_C_ipjm = gl%ge%ln_gamma_C
                        ln_gamma_R_ipjm = gl%ge%ln_gamma_R
                        !        alpha_dep_ipjm = der_dep(1)
                        !
                        !        !Decrease xi, increase xj
                        gl%molfractions(i) = molfractions_orig(i) - delta_x
                        gl%molfractions(j) = molfractions_orig(j) + delta_x
                        gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp)
                        call reduced_parameters_calc(gl, Temperature)
                        Temp_pm = gl%tredmix/tredmix_orig * Temperature
                        Dens_pm = gl%rhoredmix/rhoredmix_orig * Density
                        call gE_COSMO_SAC_MIXDERIVS(gl, Temp_pm, getder_dep, C_or_R, errval)
                        gE_C_imjp = gl%ge%gE_C
                        gE_R_imjp = gl%ge%gE_R
                        ln_gamma_C_imjp = gl%ge%ln_gamma_C
                        ln_gamma_R_imjp = gl%ge%ln_gamma_R
                        !        alpha_dep_imjp = der_dep(1)
                        !
                        !         !decrease xi, decrease xj
                        gl%molfractions(i) = molfractions_orig(i) - delta_x
                        gl%molfractions(j) = molfractions_orig(j) - delta_x
                        gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) + 2.D0 * delta_x
                        call reduced_parameters_calc(gl, Temperature)
                        Temp_pm = gl%tredmix/tredmix_orig * Temperature
                        Dens_pm = gl%rhoredmix/rhoredmix_orig * Density
                        call gE_COSMO_SAC_MIXDERIVS(gl, Temp_pm, getder_dep, C_or_R, errval)
                        gE_C_imjm = gl%ge%gE_C
                        gE_R_imjm = gl%ge%gE_R
                        ln_gamma_C_imjm = gl%ge%ln_gamma_C
                        ln_gamma_R_imjm = gl%ge%ln_gamma_R
                        !        alpha_dep_imjm = der_dep(1)
                        !
                        !        DERIVFXiXj_num(i,j) = (alpha_dep_ipjp + alpha_dep_imjm - alpha_dep_ipjm - alpha_dep_imjp) / 4.D0 / delta_x**2
                        gl%ge%gE_C_dxadxb(1,i,j) = (gE_C_ipjp(1) + gE_C_imjm(1) - gE_C_ipjm(1) - gE_C_imjp(1)) / 4.D0 / delta_x**2
                        gl%ge%gE_R_dxadxb(1,i,j) = (gE_R_ipjp(1) + gE_R_imjm(1) - gE_R_ipjm(1) - gE_R_imjp(1)) / 4.D0 / delta_x**2
                        !        DERIVFXiXj_num(j,i) = DERIVFXiXj_num(i,j)
                        gl%ge%gE_C_dxadxb(1,j,i) = gl%ge%gE_C_dxadxb(1,i,j)
                        gl%ge%gE_R_dxadxb(1,j,i) = gl%ge%gE_R_dxadxb(1,i,j)
                        gl%molfractions = molfractions_orig
                    end do
                end do
                call reduced_parameters_calc(gl, Temperature)
            end if
        end if
        d2aE_dxidxj =  gl%ge%gE_C_dxadxb(1,:,:) + gl%ge%gE_R_dxadxb(1,:,:)
        !Calculate help variables d2taui_dxadxb and d2delref_dxadxb
        do i = 1, gl%ncomp-1
            do j = 1, gl%ncomp-1
                d2delref_dxadxb(i,j) = - 1.D0 / gl%u_pack * (d2rhored_dxidxj(i,j) / gl%b_HelmgE / gl%rhoredmix**2 &
                    & + d2bHelmgE_dxidxj(i,j) / gl%b_HelmgE**2 / gl%rhoredmix &
                    & - 2.D0 * drhored_dxi(i) * drhored_dxi(j) / gl%b_HelmgE / gl%rhoredmix**3 &
                    & - 2.D0 * dbHelmgE_dxi(i) * dbHelmgE_dxi(j) / gl%b_HelmgE**3 / gl%rhoredmix &
                    & - (drhored_dxi(i) * dbHelmgE_dxi(j) + drhored_dxi(j) * dbHelmgE_dxi(i)) &
                    & / gl%b_HelmgE**2 / gl%rhoredmix**2)
                do k = 1, gl%ncomp
                    d2taui_dxadxb(k,i,j) = -tau * gl%tc(k) &
                        & * (d2Tred_dxidxj(i,j) / gl%tredmix**2 - 2.D0 * dTred_dxi(i) * dTred_dxi(j) / gl%tredmix**3)
                end do
            end do
        end do
        do i = 1, gl%ncomp-1
            do j = 1, gl%ncomp-1
                !Andreas Jäger, December 2017
                !term "sum(xi ln(bi/b))" deleted from the excess based departure function, because it is assumed that it cancels out with the combinatorial part of gE
                !d2helph_dxadxb(i,j) = tau / R_HelmgE / tredmix**2 * (tredmix * d2aE_dxidxj(i,j) - aE * d2Tred_dxidxj(i,j) &
                !                  & - daE_dxi(j) * dTred_dxi(i) - daE_dxi(i) * dTred_dxi(j)) &
                !                  & + 2.D0 * tau / R_HelmgE * aE * dTred_dxi(i) * dTred_dxi(j) / tredmix**3 &
                !                  & + d2bHelmgE_dxidxj(i,j) / b_HelmgE - dbHelmgE_dxi(i) * dbHelmgE_dxi(j) / b_HelmgE**2 &
                !                  & - (dalpha_oi_r_mix_ddelref(i) * ddelref_dxa(j) - dalpha_oi_r_i_dtaui(i) * dtaui_dxa(i,j)) &
                !                  & + (dalpha_oi_r_mix_ddelref(ncomp) * ddelref_dxa(j) - dalpha_oi_r_i_dtaui(ncomp) * dtaui_dxa(ncomp,j)) &
                !                  & - (dalpha_oi_r_mix_ddelref(j) * ddelref_dxa(i) - dalpha_oi_r_i_dtaui(j) * dtaui_dxa(j,i)) &
                !                  & + (dalpha_oi_r_mix_ddelref(ncomp) * ddelref_dxa(i) - dalpha_oi_r_i_dtaui(ncomp) * dtaui_dxa(ncomp,i))
                d2helph_dxadxb(i,j) = tau / R_HelmgE / gl%tredmix**2 * (gl%tredmix * d2aE_dxidxj(i,j) - aE * d2Tred_dxidxj(i,j) &
                    & - daE_dxi(j) * dTred_dxi(i) - daE_dxi(i) * dTred_dxi(j)) &
                    & + 2.D0 * tau / R_HelmgE * aE * dTred_dxi(i) * dTred_dxi(j) / gl%tredmix**3 &
                    & - (dalpha_oi_r_mix_ddelref(i) * ddelref_dxa(j) - dalpha_oi_r_i_dtaui(i) * dtaui_dxa(i,j)) &
                    & + (dalpha_oi_r_mix_ddelref(gl%ncomp) * ddelref_dxa(j) - dalpha_oi_r_i_dtaui(gl%ncomp) * dtaui_dxa(gl%ncomp,j)) &
                    & - (dalpha_oi_r_mix_ddelref(j) * ddelref_dxa(i) - dalpha_oi_r_i_dtaui(j) * dtaui_dxa(j,i)) &
                    & + (dalpha_oi_r_mix_ddelref(gl%ncomp) * ddelref_dxa(i) - dalpha_oi_r_i_dtaui(gl%ncomp) * dtaui_dxa(gl%ncomp,i))
                do k = 1, gl%ncomp
                    d2helph_dxadxb(i,j) = d2helph_dxadxb(i,j) &
                        & - gl%molfractions(k) * (d2alpha_oi_r_mix_ddelref2(k) * ddelref_dxa(i) * ddelref_dxa(j) &
                        & + dalpha_oi_r_mix_ddelref(k) * d2delref_dxadxb(i,j)) &
                        & + gl%molfractions(k) * (d2alpha_oi_r_i_dtaui2(k) * dtaui_dxa(k,i) * dtaui_dxa(k,j) &
                        & + dalpha_oi_r_i_dtaui(k) * d2taui_dxadxb(k,i,j))
                end do
            end do
        end do
        !----

        do i = 1, gl%ncomp-1
            do j = 1, gl%ncomp-1
                DERIVFXiXj(i,j) = d2helpf_dxadxb(i,j) * help_h + help_f * d2helph_dxadxb(i,j)&
                    & + dhelpf_dxa(i) * dhelph_dxa(j) + dhelpf_dxa(j) * dhelph_dxa(i)
            end do
        end do
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !!Second numerical derivative of the excess based departure function with respect to xi and xj at constant tau and del
        !delta_x = 1.D-4
        !molfractions_orig = molfractions
        !tredmix_orig = tredmix
        !rhoredmix_orig = rhoredmix
        !
        !Do i = 1, ncomp-1
        !    do j = i+1, ncomp-1
        !        !Increase xi, increase xj
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(j) = molfractions_orig(j) + delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) - 2.D0 * delta_x
        !        call reduced_parameters_calc(gl,Temperature)
        !        Temp_pm = tredmix/tredmix_orig * Temperature
        !        Dens_pm = rhoredmix/rhoredmix_orig * Density
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_ipjp = der_dep(1)
        !
        !        !Increase xi, decrease xj
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(j) = molfractions_orig(j) - delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp)
        !        call reduced_parameters_calc(gl,Temperature)
        !        Temp_pm = tredmix/tredmix_orig * Temperature
        !        Dens_pm = rhoredmix/rhoredmix_orig * Density
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_ipjm = der_dep(1)
        !
        !        !Decrease xi, increase xj
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(j) = molfractions_orig(j) + delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp)
        !        call reduced_parameters_calc(gl,Temperature)
        !        Temp_pm = tredmix/tredmix_orig * Temperature
        !        Dens_pm = rhoredmix/rhoredmix_orig * Density
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_imjp = der_dep(1)
        !
        !         !decrease xi, decrease xj
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(j) = molfractions_orig(j) - delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) + 2.D0 * delta_x
        !        call reduced_parameters_calc(gl,Temperature)
        !        Temp_pm = tredmix/tredmix_orig * Temperature
        !        Dens_pm = rhoredmix/rhoredmix_orig * Density
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_imjm = der_dep(1)
        !
        !        DERIVFXiXj_num(i,j) = (alpha_dep_ipjp + alpha_dep_imjm - alpha_dep_ipjm - alpha_dep_imjp) / 4.D0 / delta_x**2
        !        DERIVFXiXj_num(j,i) = DERIVFXiXj_num(i,j)
        !        molfractions = molfractions_orig
        !    end do
        !end do
        !call reduced_parameters_calc(gl,Temperature)
        !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !One-fluid model for multiparameter EOS, Andreas Jäger July 2016
    else if (gl%mix_type == 19) then

        !!Second numerical derivative of the non-corresponding states model with respect to xi and xj at constant tau and del
        !delta_x = 1.D-4
        !molfractions_orig = molfractions
        !
        !Do i = 1, ncomp-1
        !    do j = i+1, ncomp-1
        !        !Increase xi, increase xj
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(j) = molfractions_orig(j) + delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) - 2.D0 * delta_x
        !        call MIXDERIVSFNR  (Temperature, Density, getder_dep, der_dep)
        !        alpha_dep_ipjp = der_dep(1)
        !
        !        !Increase xi, decrease xj
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(j) = molfractions_orig(j) - delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp)
        !        call MIXDERIVSFNR  (Temperature, Density, getder_dep, der_dep)
        !        alpha_dep_ipjm = der_dep(1)
        !
        !        !Decrease xi, increase xj
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(j) = molfractions_orig(j) + delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp)
        !        call MIXDERIVSFNR  (Temperature, Density, getder_dep, der_dep)
        !        alpha_dep_imjp = der_dep(1)
        !
        !         !decrease xi, decrease xj
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(j) = molfractions_orig(j) - delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) + 2.D0 * delta_x
        !        call MIXDERIVSFNR  (Temperature, Density, getder_dep, der_dep)
        !        alpha_dep_imjm = der_dep(1)
        !
        !        DERIVFXiXj(i,j) = (alpha_dep_ipjp + alpha_dep_imjm - alpha_dep_ipjm - alpha_dep_imjp) / 4.D0 / delta_x**2
        !        DERIVFXiXj(j,i) = DERIVFXiXj(i,j)
        !        molfractions = molfractions_orig
        !    end do
        !end do

        !Analytical derivative of the non-corresponding states based model with respect to xi and xj at constant tau and delta
        !First, calculate the residual Helmholtz energy of the last component
        der_res = 0.D0
        gl%tredmix = gl%tc(gl%ncomp)
        gl%rhoredmix = gl%rhoc(gl%ncomp)
        call FNRDERIVS(gl,temperature, density, getder_res, der_res, gl%ncomp) ! call the calculation routine for the derivatives of the fluid i
        alpha_r_i(gl%ncomp) = der_res(1)    !Residual reduced Helmholtz energy of pure component i in the mixture
        gl%tredmix = 1.D0
        gl%rhoredmix = 1.D0

        DERIVFXiXj = 2.D0 * alpha_r_i(gl%ncomp)

        Do i = 1, (gl%NCOMP - 1)
            Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, i,gl%NCOMP)
            FiN_r = der_dep(1)
            Do j = (i+1), (gl%NCOMP-1)
                Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, i,j)
                Fij_r = der_dep(1)
                Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, j,gl%NCOMP)
                FjN_r = der_dep(1)
                DERIVFXiXj(i,j) = DERIVFXiXj(i,j) + 2.D0 * Fij_r - 2.D0 * FiN_r - 2.D0 * FjN_r
                DERIVFXiXj(j,i) = DERIVFXiXj(i,j)
            end do
        end do

    end if

    end subroutine d2ar_dxidxj


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    subroutine d2ar_dxiddel (gl,TEMPERATURE, DENSITY, DERIVFXDEL)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE COMPOSITION Xi AND DELTA FOR EACH COMPONENT
    ! AT CONST. TAU AND Xj MULTIPLIED BY DELTA
    ! ACCORDING TO THE ALGORITHM BY KEN. HALL, THE MOLEFRACTION OF THE LAST
    ! COMPONENT Xn IS REPLACED BY 1 - SUM OF ALL OTHER MOLEFRACTIONS.
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A MATRIX (30)IN WHICH THE DERIVATIVE WITH RESPECT TO COMPOSITION AND DEL FOR
    !              EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------









    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVFXDEL

    double precision, dimension(nderivs):: der_res, der_dep
    integer, dimension(nderivs):: getder_res, getder_dep
    integer:: i, j
    double precision:: d_alphaijr_d_del, d_alphaiNr_d_del, d_alphajNr_d_del, d_F0ir_d_del, d_F0Nr_d_del

    !Variable for SRK and PR
    double precision, dimension(30):: d2ardxiddel
    double precision, dimension(15,30)::MIXDERFNR_dxi

    !Variables for LKP
    double precision :: B_0, C_0, D_0, B_ref, C_ref, D_ref
    double precision, dimension(30) :: daccenLKP_dxi, dzcLKP_dxi
    double precision :: zcLKP2, zcLKP3, zcLKP4, zcLKP5,zcLKP6, zcLKP7, zcLKP8, tau, tau2, tau3, del, del2, del3, del4, del5
    double precision, dimension(30) :: help_1, help_11, help_12, help_2, help_21, help_22, help_3, help_31, help_32

    !New variables for quadratic mixing rules for the residual reduced Helmholtz energy
    !Andreas, February 2016
    double precision:: dalpha_r_ij_ddel            !First delta derivative of combination term alpha^r_(ij) = 0.5 * (alpha^r_i + alpha^r_j) (1-k_(ij))
    double precision:: dalpha_r_jN_ddel, dalpha_r_iN_ddel, dalpha_r_NN_ddel
    double precision, dimension(30):: dalpha_r_i_ddel         !First derivative of reduced residual Helmholtz energies of the pure fluids wrt delta

    !New variables for numerical derivative of the excess based departure function
    !Andreas Jäger, July 2016
    double precision:: delta_x                          !Step in x for numerical derivative
    double precision:: delta_del                          !Step in del for numerical derivative
    double precision:: alpha_dep_ipdelp                   !Value of departure function at increased xi and increased del (and decreased xN)
    double precision:: alpha_dep_ipdelm                   !Value of departure function at increased xi and decreased del (and increased xN)
    double precision:: alpha_dep_imdelp                   !Value of departure function at decreased xi and increased del (and decreased xN)
    double precision:: alpha_dep_imdelm                   !Value of departure function at decreased xi and decreased del (and increased xN)
    double precision, dimension(30):: molfractions_orig !Original value of molfractions
    Double precision:: del_p, del_m
    double precision:: Temp_pm, dens_pm                 !Decrease and increased temperature and density. When xi and xj are changed and tau and delta stay constant, T and rho change
    double precision:: tredmix_orig, rhoredmix_orig
    double precision:: rhomix_ref_copy
    double precision, dimension(30):: rho_i_ref_copy
    !double precision, dimension(30):: DERIVFXDEL_num
    !Andreas Jäger, June 2017
    double precision:: const_A1
    double precision:: help_h, df_ddel
    double precision, dimension(30):: dhelph_dxa, d2helpf_dxaddel
    double precision, dimension(30):: alpha_oi_r_mix, alpha_oi_r_i
    integer, dimension(nderivs):: GETDER_i                         ! array specifier to indicate, which derivative is needed for the pure fluid residual Helmholtz energies
    double precision, dimension(nderivs)::FNRDER_i                 ! array with the computed values for the derivatives of the residual Helmholtz energy of the pure fluids
    double precision:: aE
    double precision, dimension(30):: dalpha_oi_r_mix_ddelref, dalpha_oi_r_i_dtaui
    double precision, dimension(30):: ddelref_dxa, dtaui_dxa
    double precision, dimension(30):: dTred_dxi, drhored_dxi
    double precision, dimension(30):: dbHelmgE_dxi
    double precision, dimension(30):: daE_dxi
    !Variables for UNIFAC
    !double precision, dimension(nderivs) :: gl%ge%gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gl%ge%gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau
    integer:: C_or_R
    integer:: errval
    !double precision, dimension(nderivs, 30) :: gl%ge%gE_C_dxa        !Combinatorial part of gE and derivatives with respect to tau and delta and xa
    !double precision, dimension(nderivs, 30) :: gl%ge%gE_R_dxa        !Residual part of gE and derivatives with respect to tau and delta and xa
    !double precision, allocatable :: gl%ge%ln_gamma_C_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to delta and tau and xa
    !double precision, allocatable :: gl%ge%ln_gamma_R_dxa(:,:,:)       !Residual activity coefficients and derivatives with respect to delta and tau and xa


    !New variables for PC-SAFT equation type
    double precision, dimension (10,gl%ncomp) :: x1

    !New variables for COSMO-SAC
    double precision :: Temp_p, Temp_m
    double precision, dimension(nderivs) :: gE_C_p, gE_C_m, gE_R_p, gE_R_m
    double precision, dimension(nderivs, 30) :: ln_gamma_C_p, ln_gamma_C_m, ln_gamma_R_p, ln_gamma_R_m

    double precision, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij
    logical, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij_calculated

    getder_res = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    getder_dep = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.D0
    der_dep = 0.D0
    DERIVFXDEL = 0.D0                       !Initialize return matrix DERIVFX
    d_alphaijr_d_del = 0.D0                    !Initialize departure function ij
    d_alphaiNr_d_del = 0.D0                    !Initialize departure function iN
    d_alphajNr_d_del = 0.D0                    !Initialize departure function jN
    d_F0ir_d_del = 0.D0                        !Initialize residual helmholtz pure of component i
    d_F0Nr_d_del = 0.D0                        !Initialize residual helmholtz pure of component N


    if (gl%mix_Type  ==  1) then !Lorentz Berthelot or modified mixing rules used Andreas March 2012
        DEPFUNCFNR_ij_calculated = .false.

        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, gl%NCOMP)
        d_F0Nr_d_del = der_res(2)                  !Get residial helmh. pure of comp. N

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components
            Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
            d_F0ir_d_del = der_res(2)               !residual part of component i: F0i_r
            Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP, 2, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
            d_alphaiNr_d_del = der_dep(2)           !Get departure Function for Comp. i and N
            DERIVFXDEL(i) = d_F0ir_d_del - d_F0Nr_d_del + (1.d0 - 2.d0 * gl%molfractions(i)) * gl%Fij(i,gl%NCOMP)* d_alphaiNr_d_del
            Do j = 1, (gl%NCOMP -1)                        !Sum over all components with departure function
                if(i  /=  j) then
                    if((gl%Fij(i,j)  /=  0.d0).OR.(gl%Fij(i,gl%NCOMP)  /=  0.d0).OR.(gl%Fij(j,gl%NCOMP)  /=  0.d0)) then
                        Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, j, 2, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                        d_alphaijr_d_del = der_dep(2)           !Get departure Function for Comp. i and N
                        Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, j, gl%NCOMP, 2, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                        d_alphajNr_d_del = der_dep(2)           !Get departure Function for Comp. j and N
                        DERIVFXDEL(i) = DERIVFXDEL(i)+ gl%molfractions(j)*(gl%Fij(i,j) * d_alphaijr_d_del  - &
                            & gl%Fij(i,gl%NCOMP) * d_alphaiNr_d_del - gl%Fij(j,gl%NCOMP) * d_alphajNr_d_del)
                    end if
                end if
            end do
        end do

        !SRK with SRK mixing rules used
    else if ((gl%mix_Type  ==  2) .or. (gl%mix_type == 21) .or. (gl%mix_type == 22)) then
        !Replaced with new cubics routine
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_dxi_CUBIC (gl,Temperature, density, getder_res, MIXDERFNR_dxi)
        DERIVFXDEL = MIXDERFNR_dxi(2,:)
        !Old version
        !call d2ar_dxiddel_SRK(gl,temperature, density, d2ardxiddel)
        !DERIVFXDEL = d2ardxiddel

        !PR with PR mixing rules
    else if ((gl%mix_Type  ==  3) .or. (gl%mix_type == 31)) then
        !Replaced with new cubics routine
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_dxi_CUBIC (gl,Temperature, density, getder_res, MIXDERFNR_dxi)
        DERIVFXDEL = MIXDERFNR_dxi(2,:)
        !Old version
        !call d2ar_dxiddel_PR(temperature, density, d2ardxiddel)
        !DERIVFXDEL = d2ardxiddel

        !LKP with LKP mixing rules
    else if (gl%mix_Type  ==  4) then


        tau = gl%tredmix/temperature
        del = density/gl%rhoredmix

        del2=del*del
        del3=del2*del
        del4=del*del3
        del5=del4*del
        zcLKP2=gl%zcLKP*gl%zcLKP
        zcLKP3=zcLKP2*gl%zcLKP
        zcLKP4=zcLKP2*zcLKP2
        zcLKP5=gl%zcLKP*zcLKP4
        zcLKP6=gl%zcLKP*zcLKP5
        zcLKP7=zcLKP6*gl%zcLKP
        zcLKP8=zcLKP7*gl%zcLKP
        tau2=tau*tau
        tau3=tau*tau2


        do i = 1, gl%ncomp - 1
            daccenLKP_dxi(i) = gl%accen(i) - gl%accen(gl%ncomp)
            dzcLKP_dxi(i)=-0.085d0 * daccenLKP_dxi(i)
        end do

        B_0 = gl%lkp_b1_0 - gl%lkp_b2_0 * tau - gl%lkp_b3_0 * tau2 - gl%lkp_b4_0 * tau3
        C_0 = gl%lkp_c1_0 - gl%lkp_c2_0 * tau + gl%lkp_c3_0 * tau3
        D_0 = gl%lkp_d1_0 + gl%lkp_d2_0 * tau
        B_ref = gl%lkp_b1_ref - gl%lkp_b2_ref * tau - gl%lkp_b3_ref * tau2 - gl%lkp_b4_ref * tau3
        C_ref = gl%lkp_c1_ref - gl%lkp_c2_ref * tau + gl%lkp_c3_ref * tau3
        D_ref = gl%lkp_d1_ref + gl%lkp_d2_ref * tau

        do i = 1, gl%ncomp - 1
            !Berechne help_1(i)
            help_11(i) = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
            help_12(i) = del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0
            help_1(i) = 2.d0 * gl%lkp_c4_0 * del * tau3 * help_11(i) * dzcLKP_dxi(i) / zcLKP3 - 2.d0 * C_0 * del * dzcLKP_dxi(i) / zcLKP3 &
                & - 5.d0 * D_0 * del4 * dzcLKP_dxi(i) / zcLKP6 - B_0 * dzcLKP_dxi(i) / zcLKP2 &
                & + 2.d0 * gl%lkp_c4_0 * del * tau3 * help_11(i) * (gl%lkp_beta_0 + 1.d0) * dzcLKP_dxi(i) / zcLKP3 &
                & - 2.d0 * gl%lkp_c4_0 * del3 * gl%lkp_gamma_0 * tau3 * help_11(i) * dzcLKP_dxi(i) / zcLKP5 &
                & - 4.d0 * gl%lkp_c4_0 * gl%lkp_gamma_0 * del * tau3 * help_11(i) * help_12(i) * dzcLKP_dxi(i) / zcLKP5 &
                & + 2.d0 * gl%lkp_c4_0 * gl%lkp_gamma_0 ** 2 * del3 * tau3 * help_11(i) * help_12(i) * dzcLKP_dxi(i) / zcLKP7
            !Berechne help_2(i)
            help_21(i) = dexp(-gl%lkp_gamma_0 * del2 / zcLKP2)
            help_22(i) = del2 + (gl%lkp_beta_0 + 1.d0) * zcLKP2 / gl%lkp_gamma_0
            help_2(i) = daccenLKP_dxi(i) / gl%lkp_w_ref * (1.d0 / del + B_0 / gl%zcLKP + C_0 * del / zcLKP2 &
                & + D_0 * del4 / zcLKP5 - gl%lkp_c4_0 * tau3 * del * help_21(i) / zcLKP2 &
                & + gl%lkp_c4_0 * del * gl%lkp_gamma_0 * tau3 * help_21(i) * help_22(i) / zcLKP4) &
                & - gl%accenLKPMix / gl%lkp_w_ref * (B_0 * dzcLKP_dxi(i) / zcLKP2 + 2.d0 * C_0 * del * dzcLKP_dxi(i) / zcLKP3 &
                & + 5.d0 * D_0 * del4 * dzcLKP_dxi(i) / zcLKP6 &
                & - 2.d0 * gl%lkp_c4_0 * del * tau3 * help_21(i) * dzcLKP_dxi(i) / zcLKP3 &
                & - 2.d0 * gl%lkp_c4_0 * del * tau3 * help_21(i) * (gl%lkp_beta_0 + 1.d0) * dzcLKP_dxi(i) / zcLKP3 &
                & + 2.d0 * gl%lkp_c4_0 * del3 * gl%lkp_gamma_0 * tau3 * help_21(i) * dzcLKP_dxi(i) / zcLKP5 &
                & + 4.d0 * gl%lkp_c4_0 * gl%lkp_gamma_0 * del * tau3 * help_21(i) * help_22(i) * dzcLKP_dxi(i) / zcLKP5 &
                & - 2.d0 * gl%lkp_c4_0 * gl%lkp_gamma_0 ** 2 * del3 * tau3 * help_21(i) * help_22(i) * dzcLKP_dxi(i) / zcLKP7)
            !Berechne help_3(i)
            help_31(i) = dexp(-gl%lkp_gamma_ref * del2 / zcLKP2)
            help_32(i) = del2 + (gl%lkp_beta_ref + 1.d0) * zcLKP2 / gl%lkp_gamma_ref
            help_3(i) = daccenLKP_dxi(i) / gl%lkp_w_ref * (1.d0 / del + B_ref / gl%zcLKP + C_ref * del / zcLKP2 &
                & + D_ref * del4 / zcLKP5 - gl%lkp_c4_ref * tau3 * del * help_31(i) / zcLKP2 &
                & + gl%lkp_c4_ref * del * gl%lkp_gamma_ref * tau3 * help_31(i) * help_32(i) / zcLKP4) &
                & - gl%accenLKPMix / gl%lkp_w_ref * (B_ref * dzcLKP_dxi(i) / zcLKP2 + 2.d0 * C_ref * del * dzcLKP_dxi(i) / zcLKP3 &
                & + 5.d0 * D_ref * del4 * dzcLKP_dxi(i) / zcLKP6 &
                & - 2.d0 * gl%lkp_c4_ref * del * tau3 * help_31(i) * dzcLKP_dxi(i) / zcLKP3 &
                & - 2.d0 * gl%lkp_c4_ref * del * tau3 * help_31(i) * (gl%lkp_beta_ref + 1) * dzcLKP_dxi(i) / zcLKP3 &
                & + 2.d0 * gl%lkp_c4_ref * del3 * gl%lkp_gamma_ref * tau3 * help_31(i) * dzcLKP_dxi(i) / zcLKP5 &
                & + 4.d0 * gl%lkp_c4_ref * gl%lkp_gamma_ref * del * tau3 * help_31(i) * help_32(i) * dzcLKP_dxi(i) / zcLKP5 &
                & - 2.d0 * gl%lkp_c4_ref * gl%lkp_gamma_ref ** 2 * del3 * tau3 * help_31(i) * help_32(i) * dzcLKP_dxi(i) / zcLKP7)
            !Berechne D2alpha_r_Dx_i_Ddelta
            DERIVFXDEL(i) = (help_1(i) - help_2(i) + help_3(i))*del
        end do

        !PC-SAFT
    else if (gl%mix_Type == 6) then

        call ARX1DERIVS(gl,temperature, density, getder_res, x1)
        DERIVFXDEL(1:gl%ncomp) = x1(2,:)

        !transformation to usual writing
        do i = 1, gl%ncomp
            DERIVFXDEL(i) = DERIVFXDEL(i)- DERIVFXDEL(gl%ncomp)
        end do


        !Quadratic mixing rules for the residual reduced Helmholtz energy, Andreas February 2015
    else if (gl%mix_type  ==  11) then

        DEPFUNCFNR_ij_calculated = .false.

        !First, calculate all reduced residual Helmholtz energies of the pure fluids
        der_res = 0.D0
        do i = 1,gl%ncomp
            call FNRDERIVS(gl,temperature, density, getder_res, der_res, i) ! call the calculation routine for the derivatives of the fluid i
            dalpha_r_i_ddel(i) = der_res(2)    !Residual reduced Helmholtz energy of pure component i in the mixture
        end do

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xncomp-1
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components

            !Contribution of the quadratic mixing term for the pure reduced residual Helmholtz energies
            dalpha_r_NN_ddel = dalpha_r_i_ddel(gl%ncomp)
            dalpha_r_iN_ddel = 0.5D0 * (dalpha_r_i_ddel(i) + dalpha_r_i_ddel(gl%ncomp)) * (1.D0 - ACCESS_KIJ_HELM(i,gl%ncomp))
            DERIVFXDEL(i) = DERIVFXDEL(i) + 2.D0 * gl%molfractions(gl%ncomp) * (dalpha_r_iN_ddel - dalpha_r_NN_ddel)

            !Contribution of the departure function
            Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP, 2, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
            d_alphaiNr_d_del = der_dep(2)           !Get departure Function for Comp. i and N
            DERIVFXDEL(i) = DERIVFXDEL(i) + (1.d0 - 2.d0 * gl%molfractions(i)) * gl%Fij(i,gl%NCOMP)* d_alphaiNr_d_del

            Do j = 1, (gl%NCOMP-1)                         !Sum over all components with departure function
                !Contribution of the quadratic mixing term for the pure reduced residual Helmholtz energies
                if (i == j) then
                    dalpha_r_ij_ddel = dalpha_r_i_ddel(i)
                end if
                if (i /= j) then
                    dalpha_r_ij_ddel = 0.5D0 * (dalpha_r_i_ddel(i) + dalpha_r_i_ddel(j)) * (1.D0 - ACCESS_KIJ_HELM(i,j))
                end if
                dalpha_r_jN_ddel = 0.5D0 * (dalpha_r_i_ddel(j) + dalpha_r_i_ddel(gl%ncomp)) * (1.D0 - ACCESS_KIJ_HELM(j,gl%ncomp))
                DERIVFXDEL(i) = DERIVFXDEL(i) + 2.D0 * gl%molfractions(j) * (dalpha_r_ij_ddel - dalpha_r_jN_ddel)

                !Contributions of the departure function
                if(i  /=  j) then
                    !if((Fij(i,j)  /=  0.d0).OR.(Fij(i,NCOMP)  /=  0.d0).OR.(Fij(j,NCOMP)  /=  0.d0)) then
                    Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, j, 2, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                    d_alphaijr_d_del = der_dep(2)           !Get departure Function for Comp. i and N
                    Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, j, gl%NCOMP, 2, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                    d_alphajNr_d_del = der_dep(2)           !Get departure Function for Comp. j and N
                    DERIVFXDEL(i) = DERIVFXDEL(i)+ gl%molfractions(j)*(gl%Fij(i,j) * d_alphaijr_d_del  - &
                        & gl%Fij(i,gl%NCOMP) * d_alphaiNr_d_del - gl%Fij(j,gl%NCOMP) * d_alphajNr_d_del)
                    !end if
                end if
            end do
        end do


        !Excess based departure function, Andreas Jäger July 2016, modified April 2018 by Erik
    else if ((gl%mix_type  ==  12) .or. (gl%mix_type  ==  13))then
        if(.not. allocated(gl%ge%ln_gamma_C_dxa)) allocate(gl%ge%ln_gamma_C_dxa(nderivs, 30, 30))
        if(.not. allocated(gl%ge%ln_gamma_R_dxa)) allocate(gl%ge%ln_gamma_R_dxa(nderivs, 30, 30))

        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, gl%NCOMP)
        d_F0Nr_d_del = der_res(2)                  !Get residial helmh. pure of comp. N

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components
            Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
            d_F0ir_d_del = der_res(2)               !residual part of component i: F0i_r
            DERIVFXDEL(i) = d_F0ir_d_del - d_F0Nr_d_del
        end do

        !DERIVFXDEL_num = DERIVFXDEL

        tredmix_orig = gl%tredmix
        rhoredmix_orig = gl%rhoredmix
        tau = gl%tredmix / temperature
        del = density / gl%rhoredmix
        !Calculate analytical derivative of the excess based departure function with respect to xi and del
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !Compute auxiliary function f(del,x)
        !----
        const_A1 = -dlog(1.D0 / gl%u_pack + 1.D0)
        df_ddel = - (gl%b_HelmgE * gl%rhoredmix) / (1.D0 + gl%b_HelmgE * del * gl%rhoredmix) / const_A1
        !----

        !Compute second derivative of auxiliary function f(del,x) with respect to xa and delta
        !----
        call dYr_dxi(gl,dTred_dxi, drhored_dxi)
        call db_HelmhE_dxi(gl,dbHelmgE_dxi)
        do j = 1, gl%ncomp-1
            d2helpf_dxaddel(j) = -(gl%rhoredmix * dbHelmgE_dxi(j) + drhored_dxi(j) * gl%b_HelmgE) &
                & / const_A1 / (1.D0 + del * gl%b_HelmgE * gl%rhoredmix)**2
        end do
        !----

        !Compute auxiliary function h(tau,x)
        !----
        GETDER_i = 0
        GETDER_i(1) = 1 !The Helmholtz energy itself is needed
        GETDER_i(2) = 1 !Get the first derivative of the Helmholtz energy with respect to deltaref
        do i = 1, gl%ncomp
            !Get the reduced Helmholtz energies of the pure components at mixture conditions
            rhomix_ref_copy = gl%rho_mix_ref
            call FNRDERIVS(gl,temperature, rhomix_ref_copy, GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_oi_r_mix(i) = FNRDER_i(1)
            dalpha_oi_r_mix_ddelref(i) = FNRDER_i(2) / (gl%rho_mix_ref/gl%rhoredmix)
        end do
        GETDER_i(2) = 0
        GETDER_i(4) = 1 !Get the first derivative of the Helmholtz energy with respect to taui of the pure fluid part
        do i = 1, gl%ncomp
            !Get the reduced Helmholtz energies of the pure components at pure fluid conditions
            gl%tredmix = gl%tc(i)         !Needed to evaluate the pure fluid at correct tau (tau_i)
            gl%rhoredmix = gl%rhoc(i)     !Needed to evaluate the pure fluid at correct delta (delta_i_ref)
            rho_i_ref_copy(i) = gl%rho_i_ref(i)
            call FNRDERIVS(gl,temperature, rho_i_ref_copy(i), GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
            alpha_oi_r_i(i) = FNRDER_i(1)
            dalpha_oi_r_i_dtaui(i) = FNRDER_i(4) / (gl%tc(i)/temperature)
        end do
        gl%tredmix = tredmix_orig
        gl%rhoredmix = rhoredmix_orig

        C_or_R = 2  !Get residual part of UNIFAC
        getder_dep = 0
        getder_dep(1) = 1 !gE itself is needed, not the derivative of gE with respect to delta (which is 0). The same is true for dgE/dxa
        !Erik, April 2018
        if (gl%mix_type  ==  12) then
            !call gE_UNIFAC_MIXDERIVS(gl,temperature, getder_dep, gl%ge%gE_C, gl%ge%gE_R, gl%ge%ln_gamma_C, gl%ge%ln_gamma_R, C_or_R, errval)
            call gE_UNIFAC_MIXDERIVS(gl,temperature, getder_dep, C_or_R, errval)
        elseif (gl%mix_type  ==  13) then
            call gE_COSMO_SAC_MIXDERIVS(gl,temperature, getder_dep, C_or_R, errval)
        end if
        aE = gl%ge%gE_C(1) + gl%ge%gE_R(1)    !Specify which parts of gE are considered (combinatorial, residual, or both)

        help_h = aE / R_HelmgE / temperature
        do i = 1, gl%ncomp
            !Andreas Jäger, December 2017
            !term "sum(xi ln(bi/b))" deleted from the excess based departure function, because it is assumed that it cancels out with the combinatorial part of gE
            !help_h = help_h - molfractions(i) * dlog(bi_HelmgE(i) / b_HelmgE) - molfractions(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
            help_h = help_h - gl%molfractions(i) * (alpha_oi_r_mix(i) - alpha_oi_r_i(i))
        end do
        !----

        !Compute derivative of auxiliary function h(tau,x) with respect to xa
        !----
        !Erik, April 2018
        if (gl%mix_type  ==  12) then
            !call gE_UNIFAC_MIXDERIVS_dxa(gl,temperature, getder_dep, gl%ge%gE_C_dxa, gl%ge%gE_R_dxa, gl%ge%ln_gamma_C_dxa, gl%ge%ln_gamma_R_dxa, C_or_R, errval)
            call gE_UNIFAC_MIXDERIVS_dxa(gl,temperature, getder_dep, C_or_R, errval)
        elseif (gl%mix_type  ==  13) then
            !!Numerical derivative of the excess based departure function with respect to xi at constant tau and del
            delta_x = 1.D-4
            molfractions_orig = gl%molfractions
            tredmix_orig = gl%tredmix
            rhoredmix_orig = gl%rhoredmix
            do i =1, gl%ncomp-1
                !Increase xi
                gl%molfractions(i) = molfractions_orig(i) + delta_x
                gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) - delta_x
                call reduced_parameters_calc(gl, Temperature)
                Temp_p = gl%tredmix/tredmix_orig * Temperature
                !Dens_p = gl%rhoredmix/rhoredmix_orig * Density
                call gE_COSMO_SAC_MIXDERIVS(gl, Temp_p, getder_dep, C_or_R, errval)
                gE_C_p = gl%ge%gE_C
                gE_R_p = gl%ge%gE_R
                ln_gamma_C_p = gl%ge%ln_gamma_C
                ln_gamma_R_p = gl%ge%ln_gamma_R
                !    alpha_dep_p = der_dep(1)
                !Decrease xi
                gl%molfractions(i) = molfractions_orig(i) - delta_x
                gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) + delta_x
                call reduced_parameters_calc(gl, Temperature)
                Temp_m = gl%tredmix/tredmix_orig * Temperature
                !    Dens_m = gl%rhoredmix/rhoredmix_orig * Density
                call gE_COSMO_SAC_MIXDERIVS(gl, Temp_m, getder_dep, C_or_R, errval)
                gE_C_m = gl%ge%gE_C
                gE_R_m = gl%ge%gE_R
                ln_gamma_C_m = gl%ge%ln_gamma_C
                ln_gamma_R_m = gl%ge%ln_gamma_R
                !    alpha_dep_m = der_dep(1)
                !Calculate the numerical derivative
                gl%ge%gE_C_dxa(1,i) = (gE_C_p(1) - gE_C_m(1)) / (2.D0 * delta_x)
                gl%ge%gE_R_dxa(1,i) = (gE_R_p(1) - gE_R_m(1)) / (2.D0 * delta_x)
                gl%molfractions = molfractions_orig
            end do
            call reduced_parameters_calc(gl, Temperature)
        end if
        daE_dxi = gl%ge%gE_C_dxa(1,:) + gl%ge%gE_R_dxa(1,:)
        do j = 1, gl%ncomp-1
            !Andreas Jäger, December 2017
            !term "sum(xi ln(bi/b))" deleted from the excess based departure function, because it is assumed that it cancels out with the combinatorial part of gE
            !dhelph_dxa(j) = tau / R_HelmgE / tredmix**2 * (tredmix * daE_dxi(j) - aE * dTred_dxi(j)) &
            !            & + dlog(bi_HelmgE(ncomp)/bi_HelmgE(j)) + dbHelmgE_dxi(j) / b_HelmgE &
            !            & - (alpha_oi_r_mix(j) - alpha_oi_r_i(j)) + (alpha_oi_r_mix(ncomp) - alpha_oi_r_i(ncomp))
            dhelph_dxa(j) = tau / R_HelmgE / gl%tredmix**2 * (gl%tredmix * daE_dxi(j) - aE * dTred_dxi(j)) &
                & - (alpha_oi_r_mix(j) - alpha_oi_r_i(j)) + (alpha_oi_r_mix(gl%ncomp) - alpha_oi_r_i(gl%ncomp))
            ddelref_dxa(j) = - 1.D0 / gl%u_pack * (drhored_dxi(j) / gl%b_HelmgE / gl%rhoredmix**2 + dbHelmgE_dxi(j) / gl%b_HelmgE**2 / gl%rhoredmix)
            do i = 1, gl%ncomp
                dtaui_dxa(j) =  - tau * gl%tc(i) / gl%tredmix**2 * dTred_dxi(j)
                dhelph_dxa(j) = dhelph_dxa(j) - gl%molfractions(i) * dalpha_oi_r_mix_ddelref(i) * ddelref_dxa(j) &
                    & + gl%molfractions(i) * dalpha_oi_r_i_dtaui(i) * dtaui_dxa(j)
            end do
        end do
        !----

        DERIVFXDEL = DERIVFXDEL + (d2helpf_dxaddel * help_h + df_ddel * dhelph_dxa) * del
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !!Second numerical derivative of the excess based departure function with respect to xi and tau at constant delta
        !delta_x = 1.D-6
        !delta_del = 1.D-6
        !molfractions_orig = molfractions
        !tredmix_orig = tredmix
        !rhoredmix_orig = rhoredmix
        !del = Density / rhoredmix
        !getder_res = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        !getder_dep = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        !
        !Do i = 1, ncomp-1
        !        !Increase xi, increase del
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !        call reduced_parameters_calc(gl,Temperature)
        !        del_p = del + delta_del
        !        Dens_pm = rhoredmix * del_p
        !        Temp_pm = tredmix/tredmix_orig * Temperature
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_ipdelp = der_dep(1)
        !
        !        !Increase xi, decrease del
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !        call reduced_parameters_calc(gl,Temperature)
        !        del_p = del - delta_del
        !        Dens_pm = rhoredmix * del_p
        !        Temp_pm = tredmix/tredmix_orig * Temperature
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_ipdelm = der_dep(1)
        !
        !        !Decrease xi, increase del
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !        call reduced_parameters_calc(gl,Temperature)
        !        del_p = del + delta_del
        !        Dens_pm = rhoredmix * del_p
        !        Temp_pm = tredmix/tredmix_orig * Temperature
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_imdelp = der_dep(1)
        !
        !         !decrease xi, decrease del
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !        call reduced_parameters_calc(gl,Temperature)
        !        del_p = del - delta_del
        !        Dens_pm = rhoredmix * del_p
        !        Temp_pm = tredmix/tredmix_orig * Temperature
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_imdelm = der_dep(1)
        !
        !        DERIVFXDEL_num(i) = DERIVFXDEL_num(i) + (alpha_dep_ipdelp + alpha_dep_imdelm - alpha_dep_ipdelm - alpha_dep_imdelp) / 4.D0 / delta_x / delta_del * del
        !        molfractions = molfractions_orig
        !end do
        !call reduced_parameters_calc(gl,Temperature)
        !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        !One-fluid model for multiparameter EOS, Andreas Jäger July 2016
    else if (gl%mix_type == 19) then

        !!Second numerical derivative of the non-corresponding states based model with respect to xi and tau at constant delta
        !delta_x = 1.D-4
        !molfractions_orig = molfractions
        !getder_dep = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        !
        !Do i = 1, ncomp-1
        !        !Increase xi
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !        call MIXDERIVSFNR (Temperature, Density, getder_dep, der_dep)
        !        alpha_dep_ipdelp = der_dep(2)
        !
        !        !Decrease xi
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !        call MIXDERIVSFNR (Temperature, Density, getder_dep, der_dep)
        !        alpha_dep_imdelp = der_dep(2)
        !
        !        DERIVFXDEL(i) = DERIVFXDEL(i) + (alpha_dep_ipdelp - alpha_dep_imdelp) / (2.D0 * delta_x)
        !        molfractions = molfractions_orig
        !end do


        !Analytical derivative of the one-fluid model at constant tau and del (same as constant T and rho)
        gl%tredmix = gl%tc(gl%ncomp)
        gl%rhoredmix = gl%rhoc(gl%ncomp)
        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, gl%NCOMP)
        gl%tredmix = 1.D0
        gl%rhoredmix = 1.D0
        d_F0Nr_d_del = der_res(2)                  !Get residual helmh. pure of comp. N

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components
            gl%tredmix = gl%tc(i)
            gl%rhoredmix = gl%rhoc(i)
            Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
            gl%tredmix = 1.D0
            gl%rhoredmix = 1.D0
            d_F0ir_d_del = der_res(2)                      !residual part of component i: F0i_r
            Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP)
            d_alphaiNr_d_del = der_dep(2)                  !Get departure Function for Comp. i and N
            DERIVFXDEL(i) = 2.D0 * gl%molfractions(i) * d_F0ir_d_del - 2.D0 * gl%molfractions(gl%ncomp) * d_F0Nr_d_del + 2.D0 * gl%molfractions(gl%ncomp) * d_alphaiNr_d_del
            Do j = 1, (gl%NCOMP-1)
                Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, j, gl%NCOMP)
                d_alphajNr_d_del = der_dep(2)          !Get departure Function for Comp. i and N
                DERIVFXDEL(i) = DERIVFXDEL(i) - 2.D0 * gl%molfractions(j) * d_alphajNr_d_del
                if(i  /=  j) then
                    Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, i, j)
                    d_alphaijr_d_del = der_dep(2)
                    DERIVFXDEL(i) = DERIVFXDEL(i) + 2.D0 * gl%molfractions(j) * d_alphaijr_d_del
                end if
            end do
        end do


    end if

    end subroutine d2ar_dxiddel

    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    subroutine d2ar_dxiddel_TV (gl,TEMPERATURE, DENSITY, DERIVFXDELTV)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE COMPOSITION Xi AND DELTA FOR EACH COMPONENT
    ! AT CONST. T, V and Xj
    ! ACCORDING TO THE ALGORITHM BY KEN. HALL, THE MOLEFRACTION OF THE LAST
    ! COMPONENT Xn IS REPLACED BY 1 - SUM OF ALL OTHER MOLEFRACTIONS.
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A MATRIX (30) IN WHICH THE DERIVATIVE WITH RESPECT TO COMPOSITION AND DEL FOR
    !               EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVFXDELTV

    double precision, dimension(nderivs):: der_mix
    integer, dimension(nderivs):: getder_mix

    double precision, dimension(30):: d_rhored_d_xi, d_tred_d_xi, d2ar_ddel_dxi
    double precision:: d2ar_ddel2, d2ar_ddel_dtau, del, tau

    integer:: i

    getder_mix = (/0,0,1,0,0,1,0,0,0,0,0,0,0,0,0/)
    der_mix = 0.D0
    DERIVFXDELTV = 0.D0                        !Initialize return matrix DERIVFX
    d2ar_ddel2 = 0.D0                          !Initialize residual helmholtz pure of component i
    d2ar_ddel_dtau = 0.D0                      !Initialize residual helmholtz pure of component N

    Call MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, getder_mix, der_mix)
    d2ar_ddel2 = der_mix(3)                     !Get residial helmh. pure of comp. N
    d2ar_ddel_dtau = der_mix(6)

    Call dYr_dxi(gl,d_tred_d_xi, d_rhored_d_xi)
    Call d2ar_dxiddel (gl,TEMPERATURE, DENSITY, d2ar_ddel_dxi)
    del = Density / gl%rhoredmix
    tau = gl%tredmix / temperature

    !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
    Do i = 1, (gl%NCOMP-1)                         !Loop over all components
        DERIVFXDELTV(i) = d2ar_ddel2 / del * Density* d_rhored_d_xi(i) / (-gl%rhoredmix**2) + &
            & d2ar_ddel_dtau / tau * d_tred_d_xi(i) / Temperature + d2ar_ddel_dxi(i)
    end do

    end subroutine d2ar_dxiddel_TV

    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    subroutine d2ar_dxidtau (gl,TEMPERATURE, DENSITY, DERIVFXTAU)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE COMPOSITION Xi AND TAU FOR EACH COMPONENT
    ! AT CONSTANT DEL AND Xj
    ! ACCORDING TO THE ALGORITHM BY KEN. HALL, THE MOLEFRACTION OF THE LAST
    ! COMPONENT Xn IS REPLACED BY 1 - SUM OF ALL OTHER MOLEFRACTIONS.
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 111 table 7.5)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A MATRIX (30)IN WHICH THE DERIVATIVE WITH RESPECT TO COMPOSITION AND TAU FOR
    !              EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------









    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVFXTAU

    double precision, dimension(nderivs):: der_res, der_dep
    integer, dimension(nderivs):: getder_res, getder_dep
    integer:: i, j
    double precision:: d_alphaijr_d_tau, d_alphaiNr_d_tau, d_alphajNr_d_tau, d_F0ir_d_tau, d_F0Nr_d_tau, delr

    !Variable for SRK and PR
    double precision, dimension(30):: d2ardxidtau
    double precision, dimension(15,30)::MIXDERFNR_dxi

    !Variables for LKP
    double precision :: DB_0_Dtau, DC_0_Dtau, DD_0_Dtau, DB_ref_Dtau, DC_ref_Dtau, DD_ref_Dtau
    double precision, dimension(30) :: daccenLKP_dxi, dzcLKP_dxi
    double precision :: zcLKP2, zcLKP3, zcLKP4, zcLKP5,zcLKP6, tau, tau2, del, del2, del5

    double precision, dimension(30) :: part_d2ar_dtau_dxi_0, part_d2ar_dtau_dxi_ref
    double precision :: part_id, part_ref
    double precision :: del3, del4, tau3, zcLKP8

    !New variables for quadratic mixing rules for the residual reduced Helmholtz energy
    !Andreas, February 2016
    double precision:: dalpha_r_ij_dtau            !First tau derivative of combination term alpha^r_(ij) = 0.5 * (alpha^r_i + alpha^r_j) (1-k_(ij))
    double precision:: dalpha_r_jN_dtau, dalpha_r_iN_dtau, dalpha_r_NN_dtau
    double precision, dimension(30):: dalpha_r_i_dtau         !First derivative of reduced residual Helmholtz energies of the pure fluids wrt tau

    !New variables for numerical derivative of the excess based departure function
    !Andreas Jäger, July 2016
    double precision:: delta_x                          !Step in x for numerical derivative
    double precision:: delta_tau                          !Step in tau for numerical derivative
    double precision:: alpha_dep_iptaup                   !Value of departure function at increased xi and increased tau (and decreased xN)
    double precision:: alpha_dep_iptaum                   !Value of departure function at increased xi and decreased tau (and increased xN)
    double precision:: alpha_dep_imtaup                   !Value of departure function at decreased xi and increased tau (and decreased xN)
    double precision:: alpha_dep_imtaum                   !Value of departure function at decreased xi and decreased tau (and increased xN)
    double precision, dimension(30):: molfractions_orig !Original value of molfractions
    Double precision:: tau_p, tau_m
    double precision:: Temp_pm, dens_pm                 !Decrease and increased temperature and density. When xi and xj are changed and tau and delta stay constant, T and rho change
    double precision:: tredmix_orig, rhoredmix_orig
    double precision:: rhomix_ref_copy
    double precision, dimension(30):: rho_i_ref_copy
    !double precision, dimension(30):: DERIVFXTAU_num
    !Andreas Jäger, June 2017
    double precision:: const_A1
    double precision:: dh_dtau, help_f
    double precision, dimension(30):: d2helph_dxadtau, dhelpf_dxa
    double precision, dimension(30):: dalpha_oi_r_mix_dtau, dalpha_oi_r_i_dtaui
    integer, dimension(nderivs):: GETDER_i                         ! array specifier to indicate, which derivative is needed for the pure fluid residual Helmholtz energies
    double precision, dimension(nderivs)::FNRDER_i                 ! array with the computed values for the derivatives of the residual Helmholtz energy of the pure fluids
    double precision:: aE, daE_dtau
    double precision, dimension(30):: d2alpha_oi_r_mix_ddelrefdtau
    double precision, dimension(30):: d2alpha_oi_r_i_dtaui2
    double precision, dimension(30):: ddelref_dxa, dtaui_dxa, d2taui_dxadtau
    double precision, dimension(30):: dTred_dxi, drhored_dxi
    double precision, dimension(30):: dbHelmgE_dxi
    double precision, dimension(30):: daE_dxi, d2aE_dxidtau
    !Variables for UNIFAC
    !double precision, dimension(nderivs) :: gl%ge%gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gl%ge%gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau
    integer:: C_or_R
    integer:: errval
    !double precision, dimension(nderivs, 30) :: gl%ge%gE_C_dxa        !Combinatorial part of gE and derivatives with respect to tau and delta and xa
    !double precision, dimension(nderivs, 30) :: gl%ge%gE_R_dxa        !Residual part of gE and derivatives with respect to tau and delta and xa
    !double precision, allocatable :: gl%ge%ln_gamma_C_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to delta and tau and xa
    !double precision, allocatable :: gl%ge%ln_gamma_R_dxa(:,:,:)       !Residual activity coefficients and derivatives with respect to delta and tau and xa


    !New variables for PC-SAFT equation type
    double precision, dimension (10,gl%ncomp) :: x1

    !New variables for COSMO-SAC
    double precision :: Temp_p, Temp_m
    double precision, dimension(nderivs) :: gE_C_p, gE_C_m, gE_R_p, gE_R_m
    double precision, dimension(nderivs, 30) :: ln_gamma_C_p, ln_gamma_C_m, ln_gamma_R_p, ln_gamma_R_m
    double precision, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij
    logical, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij_calculated

    !test
    !double precision, dimension(15,30) :: gE_R_Trho_dxa, gE_R_Test
    !double precision, dimension(30) :: dTred_dxa, drhored_dxa, gE_R_Trho
    !integer, dimension(nderivs):: GETDER                !array specifier to indicate, which derivative is needed

    getder_res = (/0,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)
    getder_dep = (/0,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.D0
    der_dep = 0.D0
    DERIVFXTAU = 0.D0                       !Initialize return matrix DERIVFX
    d_alphaijr_d_tau = 0.D0                    !Initialize departure function ij
    d_alphaiNr_d_tau = 0.D0                    !Initialize departure function iN
    d_alphajNr_d_tau = 0.D0                    !Initialize departure function jN
    d_F0ir_d_tau = 0.D0                        !Initialize residual helmholtz pure of component i
    d_F0Nr_d_tau = 0.D0                        !Initialize residual helmholtz pure of component N

    if (gl%mix_Type  ==  1) then !Lorentz Berthelot or modified mixing rules used Andreas March 2012
        DEPFUNCFNR_ij_calculated = .false.

        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, gl%NCOMP)
        d_F0Nr_d_tau = der_res(4)                  !Get residial helmh. pure of comp. N

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components
            Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
            d_F0ir_d_tau = der_res(4)               !residual part of component i: F0i_r
            Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP, 4, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
            d_alphaiNr_d_tau = der_dep(4)           !Get departure Function for Comp. i and N
            DERIVFXTAU(i) = d_F0ir_d_tau - d_F0Nr_d_tau + (1.D0 - 2.D0 * gl%molfractions(i)) * gl%Fij(i,gl%NCOMP)* d_alphaiNr_d_tau
            Do j = 1, (gl%NCOMP-1)                     !Sum over all components with departure function
                if(i  /=  j) then
                    if((gl%Fij(i,j)  /=  0.d0).OR.(gl%Fij(i,gl%NCOMP)  /=  0.d0).OR.(gl%Fij(j,gl%NCOMP)  /=  0.d0)) then
                        Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, j, 4, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                        d_alphaijr_d_tau = der_dep(4)           !Get departure Function for Comp. i and N
                        Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, j, gl%NCOMP, 4, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                        d_alphajNr_d_tau = der_dep(4)           !Get departure Function for Comp. j and N
                        DERIVFXTAU(i) = DERIVFXTAU(i)+ gl%molfractions(j)*(gl%Fij(i,j) * d_alphaijr_d_tau  - &
                            & gl%Fij(i,gl%NCOMP) * d_alphaiNr_d_tau - gl%Fij(j,gl%NCOMP) * d_alphajNr_d_tau)
                    end if
                end if
            end do
        end do

    else if ((gl%mix_Type  ==  2) .or. (gl%mix_type == 21) .or. (gl%mix_type == 22)) then !SRK with SRK mixing rules used
        !Replaced with new cubics routine
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_dxi_CUBIC (gl,Temperature, density, getder_res, MIXDERFNR_dxi)
        DERIVFXTAU = MIXDERFNR_dxi(4,:)
        !Old version
        !call d2ar_dxidtau_SRK(gl,temperature, density, d2ardxidtau)
        !DERIVFXTAU = d2ardxidtau

    else if ((gl%mix_Type  ==  3) .or. (gl%mix_type == 31)) then !PR with PR mixing rules used
        !Replaced with new cubics routine
        !Andreas Jäger, March 2017
        call MIXDERIVSFNR_dxi_CUBIC (gl,Temperature, density, getder_res, MIXDERFNR_dxi)
        DERIVFXTAU = MIXDERFNR_dxi(4,:)
        !Old version
        !call d2ar_dxidtau_PR(temperature, density, d2ardxidtau)
        !DERIVFXTAU = d2ardxidtau

        !LKP with LKP mixing rules
    else if (gl%mix_Type  ==  4) then

        tau = gl%tredmix/temperature
        del = density/gl%rhoredmix

        del2=del*del
        del5=del2*del2*del
        zcLKP2=gl%zcLKP*gl%zcLKP
        zcLKP3=zcLKP2*gl%zcLKP
        zcLKP5=zcLKP2*zcLKP3
        zcLKP6=gl%zcLKP*zcLKP5
        tau2=tau*tau

        do i = 1, gl%ncomp - 1
            daccenLKP_dxi(i) = gl%accen(i) - gl%accen(gl%ncomp)
            dzcLKP_dxi(i)=-0.085d0 * daccenLKP_dxi(i)
        end do


        DB_0_Dtau = -gl%lkp_b2_0 - 2.d0 * gl%lkp_b3_0 * tau - 3.d0 * gl%lkp_b4_0 * tau2
        DC_0_Dtau = -gl%lkp_c2_0 + 3.d0 * gl%lkp_c3_0 * tau2
        DD_0_Dtau = gl%lkp_d2_0
        DB_ref_Dtau = -gl%lkp_b2_ref - 2.d0 * gl%lkp_b3_ref * tau - 3.d0 * gl%lkp_b4_ref * tau2
        DC_ref_Dtau = -gl%lkp_c2_ref + 3.d0 * gl%lkp_c3_ref * tau2
        DD_ref_Dtau = gl%lkp_d2_ref

        do i = 1, gl%ncomp - 1

            !1st derivative WRT tau
            !Ideal fluid contribution
            part_id = DB_0_Dtau / gl%zcLKP * del + 0.5d0 * DC_0_Dtau / zcLKP2 * del2 + 0.2d0 * DD_0_Dtau / zcLKP5 * del5 &
                & -1.5d0 * gl%lkp_c4_0 * tau2 / gl%lkp_gamma_0 * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0 + 1.d0) &
                & * dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) + 1.5d0 * gl%lkp_c4_0 * tau2 / gl%lkp_gamma_0 * (gl%lkp_beta_0 + 1.d0)
            !Reference fluid contribution
            part_ref = DB_ref_Dtau / gl%zcLKP * del + 0.5d0 * DC_ref_Dtau / zcLKP2 * del2 + 0.2d0 * DD_ref_Dtau / zcLKP5 * del5 &
                & -1.5d0 * gl%lkp_c4_ref * tau2 / gl%lkp_gamma_ref * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref + 1.d0) &
                & * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) + 1.5d0 * gl%lkp_c4_ref * tau2 / gl%lkp_gamma_ref * (gl%lkp_beta_ref + 1.d0)

            !2nd derivative WRT tau and xi
            !Ideal fluid contribution
            part_d2ar_dtau_dxi_0(i) = (-DB_0_Dtau * del / zcLKP2 - DC_0_Dtau * del2 / zcLKP3 - DD_0_Dtau * del5 / zcLKP6) * dzcLKP_dxi(i)  &
                - 3.d0 * gl%lkp_c4_0 * tau2  / zcLKP3 * del2 * dzcLKP_dxi(i)* dexp(-gl%lkp_gamma_0 / zcLKP2 * del2) &
                * (gl%lkp_gamma_0 / zcLKP2 * del2 + gl%lkp_beta_0)

            !Reference fluid contribution
            part_d2ar_dtau_dxi_ref(i) = (-DB_ref_Dtau * del / zcLKP2 - DC_ref_Dtau * del2 / zcLKP3 - DD_ref_Dtau * del5 / zcLKP6) * dzcLKP_dxi(i)  &
                - 3.d0 * gl%lkp_c4_ref * tau2 * dzcLKP_dxi(i) * del2 / zcLKP3 * dexp(-gl%lkp_gamma_ref / zcLKP2 * del2) &
                * (gl%lkp_gamma_ref / zcLKP2 * del2 + gl%lkp_beta_ref)

            !Berechne DERIVFXTAU(i)
            DERIVFXTAU(i) = (- daccenLKP_dxi(i) / gl%lkp_w_ref * part_id + (1.d0 - gl%accenLKPMix / gl%lkp_w_ref) * part_d2ar_dtau_dxi_0(i) &
                + daccenLKP_dxi(i) / gl%lkp_w_ref * part_ref + gl%accenLKPMix / gl%lkp_w_ref * part_d2ar_dtau_dxi_ref(i)) * tau
        end do


        !PC-SAFT
    else if (gl%mix_type == 6) then

        call ARX1DERIVS(gl,temperature, density, getder_res, x1)
        DERIVFXTAU(1:gl%ncomp) = x1(4,:)

        !transformation to usual writing
        do i = 1, gl%ncomp
            DERIVFXTAU(i) = DERIVFXTAU(i)- DERIVFXTAU(gl%ncomp)
        end do


        !Quadratic mixing rules for the residual reduced Helmholtz energy, Andreas February 2015
    else if (gl%mix_type  ==  11) then
        DEPFUNCFNR_ij_calculated = .false.

        !First, calculate all reduced residual Helmholtz energies of the pure fluids
        der_res = 0.D0
        do i = 1,gl%ncomp
            call FNRDERIVS(gl,temperature, density, getder_res, der_res, i) ! call the calculation routine for the derivatives of the fluid i
            dalpha_r_i_dtau(i) = der_res(4)    !Residual reduced Helmholtz energy of pure component i in the mixture
        end do

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xncomp-1
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components

            !Contribution of the quadratic mixing term for the pure reduced residual Helmholtz energies
            dalpha_r_NN_dtau = dalpha_r_i_dtau(gl%ncomp)
            dalpha_r_iN_dtau = 0.5D0 * (dalpha_r_i_dtau(i) + dalpha_r_i_dtau(gl%ncomp)) * (1.D0 - ACCESS_KIJ_HELM(i,gl%ncomp))
            DERIVFXTAU(i) = DERIVFXTAU(i) + 2.D0 * gl%molfractions(gl%ncomp) * (dalpha_r_iN_dtau - dalpha_r_NN_dtau)

            !Contribution of the departure function
            Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP, 4, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
            d_alphaiNr_d_tau = der_dep(4)
            DERIVFXTAU(i) = DERIVFXTAU(i) + (1.D0 - 2.D0 * gl%molfractions(i)) * gl%Fij(i,gl%NCOMP)* d_alphaiNr_d_tau

            Do j = 1, (gl%NCOMP-1)                         !Sum over all components with departure function
                !Contribution of the quadratic mixing term for the pure reduced residual Helmholtz energies
                if (i == j) then
                    dalpha_r_ij_dtau = dalpha_r_i_dtau(i)
                end if
                if (i /= j) then
                    dalpha_r_ij_dtau = 0.5D0 * (dalpha_r_i_dtau(i) + dalpha_r_i_dtau(j)) * (1.D0 - ACCESS_KIJ_HELM(i,j))
                end if
                dalpha_r_jN_dtau = 0.5D0 * (dalpha_r_i_dtau(j) + dalpha_r_i_dtau(gl%ncomp)) * (1.D0 - ACCESS_KIJ_HELM(j,gl%ncomp))
                DERIVFXTAU(i) = DERIVFXTAU(i) + 2.D0 * gl%molfractions(j) * (dalpha_r_ij_dtau - dalpha_r_jN_dtau)

                !Contributions of the departure function
                if(i  /=  j) then
                    !if((Fij(i,j)  /=  0.d0).OR.(Fij(i,NCOMP)  /=  0.d0).OR.(Fij(j,NCOMP)  /=  0.d0)) then
                    Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, j, 4, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                    d_alphaijr_d_tau = der_dep(4)           !Get departure Function for Comp. i and N
                    Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, j, gl%NCOMP, 4, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                    d_alphajNr_d_tau = der_dep(4)           !Get departure Function for Comp. j and N
                    DERIVFXTAU(i) = DERIVFXTAU(i)+ gl%molfractions(j)*(gl%Fij(i,j) * d_alphaijr_d_tau  - &
                        & gl%Fij(i,gl%NCOMP) * d_alphaiNr_d_tau - gl%Fij(j,gl%NCOMP) * d_alphajNr_d_tau)
                    !end if
                end if
            end do
        end do


        !Excess based departure function, Andreas Jäger July 2016
    else if ((gl%mix_type  ==  12) .or. (gl%mix_type  ==  13)) then
        if(.not. allocated(gl%ge%ln_gamma_C_dxa)) allocate(gl%ge%ln_gamma_C_dxa(nderivs, 30, 30))
        if(.not. allocated(gl%ge%ln_gamma_R_dxa)) allocate(gl%ge%ln_gamma_R_dxa(nderivs, 30, 30))

        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, gl%NCOMP)
        d_F0Nr_d_tau = der_res(4)                  !Get residial helmh. pure of comp. N

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components
            Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
            d_F0ir_d_tau = der_res(4)               !residual part of component i: F0i_r
            DERIVFXTAU(i) = d_F0ir_d_tau - d_F0Nr_d_tau
        end do
        !DERIVFXTAU_NUM = DERIVFXTAU

        tredmix_orig = gl%tredmix
        rhoredmix_orig = gl%rhoredmix
        tau = gl%tredmix / temperature
        del = density / gl%rhoredmix
        !Calculate analytical derivative of the excess based departure function with respect to xi and tau
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !Compute auxiliary function f(del,x)
        !----
        const_A1 = -dlog(1.D0 / gl%u_pack + 1.D0)
        help_f = -dlog(gl%b_HelmgE * density + 1.D0) / const_A1
        !----

        !Compute derivative of auxiliary function f(del,x) with respect to xa
        !----
        call dYr_dxi(gl,dTred_dxi, drhored_dxi)
        call db_HelmhE_dxi(gl,dbHelmgE_dxi)
        do j = 1, gl%ncomp-1
            dhelpf_dxa(j) = -del * (gl%rhoredmix * dbHelmgE_dxi(j) + drhored_dxi(j) * gl%b_HelmgE) &
                & / const_A1 / (1.D0 + del * gl%b_HelmgE * gl%rhoredmix)
        end do
        !----

        !Compute the derivative of the auxiliary function h(tau,x) with respect to tau
        !----
        GETDER_i = 0
        GETDER_i(4) = 1 !The derivative of the residual Helmholtz energy with respect to tau is needed
        GETDER_i(6) = 1 !The second derivative of the residual Helmholtz energy with respect to tau and delta is needed
        do i = 1, gl%ncomp
            !Get the reduced Helmholtz energies of the pure components at mixture conditions
            rhomix_ref_copy = gl%rho_mix_ref
            call FNRDERIVS(gl,temperature, rhomix_ref_copy, GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
            dalpha_oi_r_mix_dtau(i) = FNRDER_i(4) / tau
            d2alpha_oi_r_mix_ddelrefdtau(i) = FNRDER_i(6) / (gl%rho_mix_ref/gl%rhoredmix) / tau
        end do
        GETDER_i(5) = 1 !For the pure component part, the second derivative of the residual Helmholtz energy with respect to taui is needed
        GETDER_i(6) = 0 !This derivative is not needed for the pure part
        do i = 1, gl%ncomp
            !Get the reduced Helmholtz energies of the pure components at pure fluid conditions
            gl%tredmix = gl%tc(i)         !Needed to evaluate the pure fluid at correct tau (tau_i)
            gl%rhoredmix = gl%rhoc(i)     !Needed to evaluate the pure fluid at correct delta (delta_i_ref)
            rho_i_ref_copy(i) = gl%rho_i_ref(i)
            call FNRDERIVS(gl,temperature, rho_i_ref_copy(i), GETDER_i, FNRDER_i, i) ! call the calculation routine for the derivatives of the fluid i
            dalpha_oi_r_i_dtaui(i) = FNRDER_i(4) / (gl%tc(i)/temperature)
            d2alpha_oi_r_i_dtaui2(i) = FNRDER_i(5) / (gl%tc(i)/temperature)**2
        end do
        gl%tredmix = tredmix_orig
        gl%rhoredmix = rhoredmix_orig

        C_or_R = 2  !Get residual part of UNIFAC
        getder_dep = 0
        getder_dep(1) = 1 !gE itself is needed. The same is true for dgE/dxa
        getder_dep(4) = 1 !The derivative of gE with respect to tau is needed. The same is true for d2gE/dxadtau
        !Erik, April 2018
        if (gl%mix_type  ==  12) then
            !call gE_UNIFAC_MIXDERIVS(gl,temperature, getder_dep, gl%ge%gE_C, gl%ge%gE_R, gl%ge%ln_gamma_C, gl%ge%ln_gamma_R, C_or_R, errval)
            call gE_UNIFAC_MIXDERIVS(gl,temperature, getder_dep, C_or_R, errval)
        elseif (gl%mix_type  ==  13) then
            call gE_COSMO_SAC_MIXDERIVS(gl,temperature, getder_dep, C_or_R, errval)
        end if
        aE = gl%ge%gE_C(1) + gl%ge%gE_R(1)    !Specify which parts of gE are considered (combinatorial, residual, or both)
        daE_dtau = gl%ge%gE_C(4) + gl%ge%gE_R(4)

        dh_dtau = (tau * daE_dtau + aE)/ R_HelmgE / gl%tredmix
        do i = 1, gl%ncomp
            dh_dtau = dh_dtau - gl%molfractions(i) * (dalpha_oi_r_mix_dtau(i) - dalpha_oi_r_i_dtaui(i) * gl%tc(i)/gl%tredmix)
        end do
        !----

        !Compute derivative of auxiliary function h(tau,x) with respect to xa
        !----
        !Erik, April 2018
        if (gl%mix_type  ==  12) then
            !call gE_UNIFAC_MIXDERIVS_dxa(gl,temperature, getder_dep, gl%ge%gE_C_dxa, gl%ge%gE_R_dxa, gl%ge%ln_gamma_C_dxa, gl%ge%ln_gamma_R_dxa, C_or_R, errval)
            call gE_UNIFAC_MIXDERIVS_dxa(gl,temperature, getder_dep, C_or_R, errval)
        elseif (gl%mix_type  ==  13) then
            if ((gl.cosmo.COSMO_ver == 1) .and. (gl.cosmo.analytical)) then
                !Analytical derivative of the excess based departure function with respect to xi at constant tau and del
                call gE_COSMO_SAC_MIXDERIVS_dxa(gl,temperature, getder_dep, C_or_R, errval)
            else
                !!Numerical derivative of the excess based departure function with respect to xi at constant tau and del
                delta_x = 1.D-4
                molfractions_orig = gl%molfractions
                tredmix_orig = gl%tredmix
                rhoredmix_orig = gl%rhoredmix
                gl%ge%gE_C_dxa = 0.D0
                gl%ge%gE_R_dxa = 0.D0
                do i =1, gl%ncomp-1
                    !Increase xi
                    gl%molfractions(i) = molfractions_orig(i) + delta_x
                    gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) - delta_x
                    call reduced_parameters_calc(gl, Temperature)
                    Temp_p = gl%tredmix/tredmix_orig * Temperature
                    !Temp_p = temperature
                    !Dens_p = gl.rhoredmix/rhoredmix_orig * Density
                    call gE_COSMO_SAC_MIXDERIVS(gl, Temp_p, getder_dep, C_or_R, errval)
                    gE_C_p = gl%ge%gE_C
                    gE_R_p = gl%ge%gE_R
                    ln_gamma_C_p = gl%ge%ln_gamma_C
                    ln_gamma_R_p = gl%ge%ln_gamma_R
                    !    alpha_dep_p = der_dep(1)
                    !Decrease xi
                    gl%molfractions(i) = molfractions_orig(i) - delta_x
                    gl%molfractions(gl%ncomp) = molfractions_orig(gl%ncomp) + delta_x
                    call reduced_parameters_calc(gl, Temperature)
                    Temp_m = gl%tredmix/tredmix_orig * Temperature
                    !Temp_m = temperature
                    !    Dens_m = gl.rhoredmix/rhoredmix_orig * Density
                    call gE_COSMO_SAC_MIXDERIVS(gl, Temp_m, getder_dep, C_or_R, errval)
                    gE_C_m = gl%ge%gE_C
                    gE_R_m = gl%ge%gE_R
                    ln_gamma_C_m = gl%ge%ln_gamma_C
                    ln_gamma_R_m = gl%ge%ln_gamma_R
                    !    alpha_dep_m = der_dep(1)
                    !Calculate the numerical derivative
                    if ((C_or_R == 1) .or. (C_or_R == 0)) then
                        gl%ge%gE_C_dxa(1,i) = (gE_C_p(1) - gE_C_m(1)) / (2.D0 * delta_x)
                        gl%ge%gE_C_dxa(4,i) = (gE_C_p(4) - gE_C_m(4)) / (2.D0 * delta_x)

                    elseif ((C_or_R == 2) .or. (C_or_R == 0)) then
                        gl%ge%gE_R_dxa(1,i) = (gE_R_p(1) - gE_R_m(1)) / (2.D0 * delta_x)
                        gl%ge%gE_R_dxa(4,i) = (gE_R_p(4) - gE_R_m(4)) / (2.D0 * delta_x)
                    end if
                    gl%molfractions = molfractions_orig
                end do
                call reduced_parameters_calc(gl, Temperature)
            end if
        end if
        daE_dxi = gl%ge%gE_C_dxa(1,:) + gl%ge%gE_R_dxa(1,:)
        d2aE_dxidtau = gl%ge%gE_C_dxa(4,:) + gl%ge%gE_R_dxa(4,:)
        do j = 1, gl%ncomp-1
            d2helph_dxadtau(j) = (gl%tredmix * daE_dxi(j) - aE * dTred_dxi(j)) / R_HelmgE / gl%tredmix**2   &
                &  + (gl%tredmix * d2aE_dxidtau(j) - daE_dtau * dTred_dxi(j)) * tau / R_HelmgE / gl%tredmix**2   &
                &  - (dalpha_oi_r_mix_dtau(j) - dalpha_oi_r_i_dtaui(j) * gl%tc(j)/gl%tredmix) &
                &  + (dalpha_oi_r_mix_dtau(gl%ncomp) - dalpha_oi_r_i_dtaui(gl%ncomp) * gl%tc(gl%ncomp)/gl%tredmix)
            ddelref_dxa(j) = - 1.D0 / gl%u_pack * (drhored_dxi(j) / gl%b_HelmgE / gl%rhoredmix**2 + dbHelmgE_dxi(j) / gl%b_HelmgE**2 / gl%rhoredmix)
            do i = 1, gl%ncomp
                dtaui_dxa(j) = - tau * gl%tc(i) / gl%tredmix**2 * dTred_dxi(j)
                d2taui_dxadtau(j) = - gl%tc(i) / gl%tredmix**2 * dTred_dxi(j)
                d2helph_dxadtau(j) = d2helph_dxadtau(j) - gl%molfractions(i) * d2alpha_oi_r_mix_ddelrefdtau(i) * ddelref_dxa(j) &
                    & + gl%molfractions(i) * d2alpha_oi_r_i_dtaui2(i) * gl%tc(i)/gl%tredmix * dtaui_dxa(j) &
                    & + gl%molfractions(i) * dalpha_oi_r_i_dtaui(i) * d2taui_dxadtau(j)
            end do
        end do
        !----

        DERIVFXTAU = DERIVFXTAU + (dhelpf_dxa * dh_dtau + help_f * d2helph_dxadtau) * tau
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !!Second numerical derivative of the excess based departure function with respect to xi and tau at constant delta
        !delta_x = 1.D-4
        !delta_tau = 1.D-4
        !molfractions_orig = molfractions
        !tredmix_orig = tredmix
        !rhoredmix_orig = rhoredmix
        !tau = tredmix / Temperature
        !getder_res = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        !getder_dep = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        !
        !Do i = 1, ncomp-1
        !        !Increase xi, increase tau
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !        call reduced_parameters_calc(gl,Temperature)
        !        tau_p = tau + delta_tau
        !        Temp_pm = tredmix/tau_p
        !        Dens_pm = rhoredmix/rhoredmix_orig * Density
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_iptaup = der_dep(1)
        !
        !        !Increase xi, decrease tau
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !        call reduced_parameters_calc(gl,Temperature)
        !        tau_p = tau - delta_tau
        !        Temp_pm = tredmix/tau_p
        !        Dens_pm = rhoredmix/rhoredmix_orig * Density
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_iptaum = der_dep(1)
        !
        !        !Decrease xi, increase tau
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !        call reduced_parameters_calc(gl,Temperature)
        !        tau_p = tau + delta_tau
        !        Temp_pm = tredmix/tau_p
        !        Dens_pm = rhoredmix/rhoredmix_orig * Density
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_imtaup = der_dep(1)
        !
        !         !decrease xi, decrease xj
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !        call reduced_parameters_calc(gl,Temperature)
        !        tau_p = tau - delta_tau
        !        Temp_pm = tredmix/tau_p
        !        Dens_pm = rhoredmix/rhoredmix_orig * Density
        !        call DEPFUNC_GE_BASED (Temp_pm, Dens_pm, getder_dep, der_dep)
        !        alpha_dep_imtaum = der_dep(1)
        !
        !        DERIVFXTAU_NUM(i) = DERIVFXTAU_NUM(i) + (alpha_dep_iptaup + alpha_dep_imtaum - alpha_dep_iptaum - alpha_dep_imtaup) / 4.D0 / delta_x / delta_tau * tau
        !        molfractions = molfractions_orig
        !end do
        !call reduced_parameters_calc(gl,Temperature)
        !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        !One-fluid model for multiparameter EOS, Andreas Jäger July 2016
    else if (gl%mix_type == 19) then

        !Second numerical derivative of the non-corresponding states model with respect to xi and tau at constant delta
        !delta_x = 1.D-4
        !molfractions_orig = molfractions
        !getder_dep = (/0,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)
        !
        !Do i = 1, ncomp-1
        !        !Increase xi
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !        call MIXDERIVSFNR (Temperature, Density, getder_dep, der_dep)
        !        alpha_dep_iptaup = der_dep(4)
        !
        !        !Decrease xi
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !        call MIXDERIVSFNR (Temperature, Density, getder_dep, der_dep)
        !        alpha_dep_imtaup = der_dep(4)
        !
        !        DERIVFXTAU(i) = DERIVFXTAU(i) + (alpha_dep_iptaup - alpha_dep_imtaup) / (2.D0 * delta_x)
        !        molfractions = molfractions_orig
        !end do

        !Analytical derivative of the one-fluid model at constant tau and del (same as constant T and rho)
        gl%tredmix = gl%tc(gl%ncomp)
        gl%rhoredmix = gl%rhoc(gl%ncomp)
        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, gl%NCOMP)
        gl%tredmix = 1.D0
        gl%rhoredmix = 1.D0
        d_F0Nr_d_tau = der_res(4)                  !Get residual helmh. pure of comp. N

        !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
        Do i = 1, (gl%NCOMP-1)                         !Loop over all components
            gl%tredmix = gl%tc(i)
            gl%rhoredmix = gl%rhoc(i)
            Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
            gl%tredmix = 1.D0
            gl%rhoredmix = 1.D0
            d_F0ir_d_tau = der_res(4)                      !residual part of component i: F0i_r
            Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, i, gl%NCOMP)
            d_alphaiNr_d_tau = der_dep(4)                  !Get departure Function for Comp. i and N
            DERIVFXTAU(i) = 2.D0 * gl%molfractions(i) * d_F0ir_d_tau - 2.D0 * gl%molfractions(gl%ncomp) * d_F0Nr_d_tau + 2.D0 * gl%molfractions(gl%ncomp) * d_alphaiNr_d_tau
            Do j = 1, (gl%NCOMP-1)
                Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, j, gl%NCOMP)
                d_alphajNr_d_tau = der_dep(4)          !Get departure Function for Comp. i and N
                DERIVFXTAU(i) = DERIVFXTAU(i) - 2.D0 * gl%molfractions(j) * d_alphajNr_d_tau
                if(i  /=  j) then
                    Call DEPFUNC_NON_COR_STATE (gl,temperature, density, getder_dep, der_dep, i, j)
                    d_alphaijr_d_tau = der_dep(4)
                    DERIVFXTAU(i) = DERIVFXTAU(i) + 2.D0 * gl%molfractions(j) * d_alphaijr_d_tau
                end if
            end do
        end do

    end if


    end subroutine d2ar_dxidtau
    !**************************************************************************


    !**************************************************************************
    subroutine dnda_dnidxj(gl,TEMPERATURE, DENSITY, dndadndx)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF alpha_r W.R.T. xj AND
    ! ni - d(n*d(ar)/d(ni))dxj. According to the algorithm by K. Hall the last component
    ! x_n is replaced by 1 - Sum(xi)

    ! A. Jäger , J. Gernert, Aug. 2010
    !--------------------------------------------------------------------------------------------------
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 121 EQ. 7.52)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! dndadndx     - A MATRIX OF(30 x 30)IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30,30):: dndadndx

    integer, dimension(nderivs)::GETDER
    double precision, dimension(nderivs)::MIXDERFNR
    double precision, dimension(30)::DERIVFXDEL, DERIVFXTAU, ndTred_dni, ndrhored_dni
    double precision, dimension(30)::dTred_dxi, drhored_dxi
    double precision, dimension(:,:), allocatable ::dndTr_dnidxj, dndrhor_dnidxj
    double precision, dimension(30):: DERIVFX, DERIVFX2
    !warnings (Moni)
    !double precision, dimension(30):: D2fx2
    double precision, dimension(:,:), allocatable :: DERIVFXiXj
    double precision:: da_dd, da_dt, dda_dxjdd, dda_dxjdt, da_dxj, dda_dxidxj, xm, dda_dxjdxm, dda_dxi2, dda_dxm2
    integer:: i, j, m

    dndadndx = 0.D0

    if (.not. allocated(dndTr_dnidxj)) then
        allocate(dndTr_dnidxj(30,30))
        allocate(dndrhor_dnidxj, DERIVFXiXj, mold=dndTr_dnidxj)
    endif
    !get the derivatives tau*d(ar)/d(tau) and del*d(ar)/d(del)
    GETDER = (/0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)

    ! get the derivatives del*d^2(ar)/d(xi)d(del) for all xi
    call d2ar_dxiddel (gl,TEMPERATURE, DENSITY, DERIVFXDEL)

    ! get the derivatives tau*d^2(ar)/d(xi)d(tau) for all xi
    call d2ar_dxidtau (gl,TEMPERATURE, DENSITY, DERIVFXTAU)

    ! get the derivatives n*d(Tred)/d(ni) and n*d(rhored)/d(ni) for all ni
    call ndYr_dni(gl,ndTred_dni, ndrhored_dni)

    ! get the derivatives d(Tred)/d(xi) and d(rhored)/d(xi) for all xi
    call dYr_dxi(gl,dTred_dxi, drhored_dxi)

    ! get the derivatives d(n*d(Tred)/d(ni))dxj and d(n*d(rhored)/d(ni))dxj for all ni, xj
    call dndYr_dnidxj(gl,dndTr_dnidxj, dndrhor_dnidxj)

    ! get the derivatives d(ar)/d(xi) for all xi
    call dar_dxi (gl,TEMPERATURE, DENSITY, DERIVFX)

    ! get the derivatives d^2(ar)/d(xi)d(xj) for all xi, xj
    call d2ar_dxidxj (gl,TEMPERATURE, DENSITY, DERIVFXiXj)

    Call d2ar_dxi2 (gl,TEMPERATURE, DENSITY, DERIVFX2)

    da_dd = MIXDERFNR(2)
    da_dt = MIXDERFNR(4)
    do j = 1, (gl%ncomp-1)
        do i = 1, (gl%ncomp)
            dda_dxjdd = DERIVFXDEL(j)
            dda_dxjdt = DERIVFXTAU(j)
            da_dxj = DERIVFX(j)
            dda_dxi2 = DERIVFX2(j)

            dndadndx(j, i) = dda_dxjdd*(1.D0 - ndrhored_dni(i)/gl%rhoredmix) + dda_dxjdt*ndTred_dni(i)/gl%tredmix &
                & - da_dd/gl%rhoredmix*(dndrhor_dnidxj(j, i) - ndrhored_dni(i)/gl%rhoredmix*drhored_dxi(j)) &
                & + da_dt/gl%tredmix*(dndtr_dnidxj(j, i) - ndtred_dni(i)/gl%tredmix*dtred_dxi(j)) &
                & - da_dxj
            if (j  ==  i) then
                dndadndx(j, i) = dndadndx(j, i) + dda_dxi2
            else if (i  <  gl%ncomp) then
                dda_dxidxj = DERIVFXiXj(j, i)
                dndadndx(j, i) = dndadndx(j, i) + dda_dxidxj
            end if
            do m = 1, (gl%ncomp-1)
                xm = gl%molfractions(m)
                if (j  /=  m) then
                    dda_dxjdxm = DERIVFXiXj(j, m)
                    dndadndx(j, i) = dndadndx(j, i) - xm * dda_dxjdxm
                else
                    dda_dxm2 = DERIVFX2(m)
                    dndadndx(j, i) = dndadndx(j, i) - xm * dda_dxm2
                end if
            end do
        end do
    end do

    end subroutine dnda_dnidxj

    !**************************************************************************
    subroutine dnda_dnidxj_TV(gl,TEMPERATURE, DENSITY, dndadndxTV)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF alpha_r W.R.T. xj AND
    ! ni - d(n*d(ar)/d(ni))dxj. According to the algorithm by K. Hall the last component
    ! x_n is replaced by 1 - Sum(xi)

    ! A. Jäger , J. Gernert, Aug. 2010
    !--------------------------------------------------------------------------------------------------
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 121 EQ. 7.52)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! dndadndx     - A MATRIX OF(30 x 30)IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30,30):: dndadndxTV
    double precision, dimension(30):: d_rhored_d_xi, d_tred_d_xi, dar_dniddel, dar_dnidtau
    double precision, dimension(30,30):: dndar_dndx
    double precision:: del, tau
    integer:: i, j

    dndadndxTV = 0.D0

    Call dYr_dxi(gl,d_tred_d_xi, d_rhored_d_xi)
    Call dnda_dnidxj(gl,TEMPERATURE, DENSITY, dndar_dndx)
    Call dndar_dniddel (gl,TEMPERATURE, DENSITY, dar_dniddel)
    Call dndar_dnidtau (gl,TEMPERATURE, DENSITY, dar_dnidtau)

    del = Density / gl%rhoredmix
    tau = gl%tredmix / temperature

    Do j = 1, (gl%NCOMP-1)                         !Loop over all components
        Do i = 1, (gl%NCOMP)
            dndadndxTV(j, i) = dar_dniddel(i)  / del * Density* d_rhored_d_xi(j) / (-gl%rhoredmix**2) + &
                & dar_dnidtau(i) / tau * d_tred_d_xi(j) / Temperature + dndar_dndx(j, i)
        End Do
    end do

    end subroutine dnda_dnidxj_TV

    !**************************************************************************
    subroutine d2na_dnidxj_TV(gl,TEMPERATURE, DENSITY, d2nadndxTV)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF alpha_r W.R.T. xj AND
    ! ni - d(n*d(ar)/d(ni))dxj. According to the algorithm by K. Hall the last component
    ! x_n is replaced by 1 - Sum(xi)

    ! A. Jäger , J. Gernert, Aug. 2010
    !--------------------------------------------------------------------------------------------------
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 121 EQ. 7.52)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! dndadndx     - A MATRIX OF(30 x 30)IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision,  dimension(:,:), allocatable :: d2nadndxTV

    double precision,  dimension(:,:), allocatable :: dnda_dndx_TV
    double precision, dimension(30):: dar_dxi

    integer:: i, j
    
    if (.not.allocated(d2nadndxTV)) allocate(d2nadndxTV(30,30))
if (.not.allocated(dnda_dndx_TV)) allocate(dnda_dndx_TV(30,30))

    d2nadndxTV = 0.D0

    Call dnda_dnidxj_TV(gl,TEMPERATURE, DENSITY, dnda_dndx_TV)
    Call dar_dxi_TV(gl,TEMPERATURE, DENSITY,dar_dxi)

    Do j = 1, (gl%NCOMP - 1)                         !Loop over all components
        Do i = 1, (gl%NCOMP)
            d2nadndxTV(j, i) = dnda_dndx_TV(j, i) + dar_dxi(j)
        End Do
    end do

    end subroutine d2na_dnidxj_TV

    !
    !**************************************************************************
    subroutine dlnphii_dxj_TP(gl,TEMPERATURE, DENSITY, dphiidxj)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE FUGACITY COEFFICIENT
    ! OF COMPONENT i W.R.T. xj AT CONSTANT T AND P
    !                            --- d(ln(phi_i)/d(xj)) ---
    ! A. Jäger , J. Gernert, Nov. 2010
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! dphiidxj    - A MATRIX OF(30 x 30)IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(:,:), allocatable :: dphiidxj

    double precision,  dimension(:,:), allocatable ::d2nadndxTV
    double precision, dimension(30)::dPdni, dPdXi
    double precision::  dPdV, Rmix
    !warnings (Moni)
    !double precision:: xk,xj
    !integer::k
    integer::i, j
    
    if (.not.allocated(dphiidxj)) allocate(dphiidxj(30,30), d2nadndxTV(30,30))
    ! get the derivative d^2(n*ar)/d(ni)d(xj) at const. T, V
    call d2na_dnidxj_TV(gl,TEMPERATURE, DENSITY, d2nadndxTV)
    ! get the derivative n*d(p)/d(ni) at const. T, V
    call ndP_dni_TV (gl,TEMPERATURE, DENSITY, dPdni)
    ! get the derivative n*d(p)/d(V) at const. T, x
    call ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    ! get the mixture Molar Gas Constant
    Call R_mix_calc(gl,Rmix)
    ! get the derivative d(p)/d(xj) at const. T, V
    call dP_dxi_TV (gl,TEMPERATURE, DENSITY, dPdXi)

    dphiidxj = 0.d0

    do j = 1, gl%ncomp-1
        do i = 1, gl%ncomp
            dphiidxj(j, i) = d2nadndxTV(j, i) + dPdni(i)/(Rmix*temperature)*dPdXi(j)/dpdV
        end do
    end do

    end subroutine dlnphii_dxj_TP
    !**************************************************************************

    !**************************************************************************
    subroutine dndar_dniddel(gl,TEMPERATURE, DENSITY, DERIVS)
    !**************************************************************************
    !subroutine for the calculation of the derivative of n*d(ar)/d(ni) with
    !respect to del, multiplied with del - del*d(n*d(ar)/d(ni))/d(del)
    !J. Gernert, Aug. 2010
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 121, Equ. 7.50)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    ! OUTPUT PARAMETERS:
    ! DERIVS      - A VECTOR (30)IN WHICH THE DERIVATIVE WITH RESPECT TO THE AMOUNT OF SUBSTANCE OF
    !               COMPONENT i AND DELTA IS STORED FOR EACH COMPONENT OF THE MIXTURE.
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVS
    double precision, dimension(30):: ndTred_dni, ndrhored_dni, DERIVFXDEL
    double precision, dimension(nderivs):: FNRDER
    integer, dimension(nderivs):: GETDER
    double precision:: dadd, d2add2, d2adddt, d2addeldxi, d2addeldxk, xk
    integer:: i, k


    ! get derivatives of alpha_r w.r.t. del and tau
    GETDER = (/0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR(gl,TEMPERATURE, DENSITY, GETDER, FNRDER)
    ! 1st der. of alpha_r w.r.t del, multiplied with del - del*d(ar)/d(del)
    dadd = FNRDER(2)
    ! 2nd der. of alpha_r w.r.t del, multiplied with del^2 - del^2*d^2(ar)/d(del)^2
    d2add2 = FNRDER(3)
    ! 1st mixed der. of alpha_r w.r.t. del and tau, multiplied with del and tau - del*tau*d^2(ar)/d(del)d(tau)
    d2adddt = FNRDER(6)

    ! get derivatives of the reducing parameters Tred and rhored w.r.t all ni
    call ndYr_dni(gl,ndTred_dni, ndrhored_dni)

    !get derivatives of alpha_r w.r.t del and all xi
    call d2ar_dxiddel (gl,TEMPERATURE, DENSITY, DERIVFXDEL)
    DERIVS = 0.d0
    ! calculate the derivative
    do i = 1, gl%ncomp
        DERIVS(i)= 0.D0
        DERIVS(i) = ((dadd + d2add2)*(1.D0 - ndrhored_dni(i)/gl%rhoredmix) &
            & + d2adddt*ndTred_dni(i)/gl%tredmix)
        if (i  <  gl%ncomp) then
            d2addeldxi = DERIVFXDEL(i)
            DERIVS(i) = DERIVS(i) + d2addeldxi
        end if
        ! inner summation over all components ...
        do k = 1, gl%ncomp-1
            xk = gl%MOLFRACTIONS(k)
            d2addeldxk = DERIVFXDEL(k)
            DERIVS(i) = DERIVS(i) - xk*d2addeldxk
        end do
    end do

    end subroutine dndar_dniddel
    !**************************************************************************

    !**************************************************************************
    subroutine dndar_dnidtau (gl,TEMPERATURE, DENSITY, DANITAU)
    !**************************************************************************
    ! A. Jäger, Aug. 2010
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF n* d(ar)/d(ni) WITH RESPECT TO tau, MULTIPLIED BY tau
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 121 EQ. 7.51)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A VECTOR (30)IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DANITAU

    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivs):: getder_res
    integer:: i,j
    double precision, dimension(30):: drhoredmix_dni, dtredmix_dni, dda_dxdt
    !dda_dtdd = second deriv of a with respect to del and tau
    !da_dt = first deriv of a with respect to tau
    !dda_dtdt = second deriv of a with respect to tau
    !dda_dxdt = second deriv of a with respect to tau and xi
    double precision:: dda_dtdd, da_dt, dda_dtdt, dda_dxidt, dda_dxjdt

    getder_res = (/0,0,0,1,1,1,0,0,0,0,0,0,0,0,0/)
    der_res = 0.D0
    DANITAU = 0.D0               !Initialize return vector
    dda_dtdd = 0.D0                     !Initialize
    da_dt  = 0.D0                       !Initialize
    dda_dtdt= 0.D0                      !Initialize

    call ndYr_dni(gl,dtredmix_dni, drhoredmix_dni)                     !Get derivatives of reduced parameters with respect to ni
    Call d2ar_dxidtau(gl,temperature, density, dda_dxdt)               !Get derivatives of a with respect to tau and xi

    Call MIXDERIVSFNR (gl,temperature, density, getder_res, der_res)   !Get derivatives of ar
    dda_dtdd = der_res(6)
    da_dt  = der_res(4)
    dda_dtdt= der_res(5)

    Do i = 1, gl%NCOMP                   !Loop over all components
        DANITAU(i) = (dda_dtdd*(1.d0-drhoredmix_dni(i)/gl%rhoredmix)+(da_dt+dda_dtdt)*dtredmix_dni(i)/gl%tredmix)
        if (i  <  gl%ncomp) then
            dda_dxidt = dda_dxdt(i)
            DANITAU(i) = DANITAU(i) + dda_dxidt
        end if
        Do j = 1, gl%NCOMP-1
            dda_dxjdt = dda_dxdt(j)
            DANITAU(i) = DANITAU(i) - gl%molfractions(j) * dda_dxjdt
        end do
    end do

    end subroutine dndar_dnidtau

    !**************************************************************************
    subroutine dP_dxi_TV (gl,TEMPERATURE, DENSITY, dPdXi)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE PRESSURE WITH RESPECT TO THE MOL FRACTION
    ! xi AT CONSTANT TEMPERATURE AND VOLUME
    !                                           d(p)/d(xi) at const. T, V, xj
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOLE/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dPdXi       - Vector(30) the derivatives for each component are stored in
    !--------------------------------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: dPdxi
    double precision, dimension(30):: d2ar_dxiddel, dTreddxi, drhoreddxi
    double precision, dimension(nderivs):: FNRDERIVS
    integer, dimension(nderivs):: getderivs
    integer:: i
    double precision:: Rmix, dar_ddel

    dPdxi = 0.d0

    ! get the derivative del*d^2(ar)/d(xi)d(del) at const. T, V
    call d2ar_dxiddel_TV (gl,TEMPERATURE, DENSITY, d2ar_dxiddel)
    ! get the 1st derivative of alpha-r w.r.t. density
    getderivs = (/0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR(gl,Temperature, Density, getderivs, FNRDERIVS)
    ! get the 1st derivative of the reducing functions w.r.t xi
    call dYr_dxi(gl,dTreddxi, drhoreddxi)
    ! get the mixture Molar Gas Constant
    Call R_mix_calc(gl,Rmix)

    dar_ddel = FNRDERIVS(2)

    do i = 1, gl%ncomp-1
        ! calculate the derivative d(p)/d(xi) at const. T, V
        dPdxi(i) = density*Rmix*Temperature*( -drhoreddxi(i)/gl%rhoredmix*dar_ddel  + d2ar_dxiddel(i))
    end do

    end subroutine dP_dxi_TV

    !**************************************************************************
    Subroutine lnf_mix(gl,T, D, p, lnfi)
    !**************************************************************************
    ! Calculation of the logarithm of the mixture fugacities (the fugacities are expressed in MPa)
    !-------------------------------------
    !input:
    !-------------------------------------
    ! T: Temperatuee (K)
    ! D: density (mol/m^3)
    ! p: pressure (MPa)
    ! output:
    !-------------------------------------
    ! lnfi: vector log(f_i) of all components
    !-------------------------------------
    !!DEC$ ATTRIBUTES DLLEXPORT :: lnf_mix






    implicit none

    type(type_gl) :: gl


    double precision:: T, D, p
    double precision, dimension(30), intent(out) :: lnfi

    integer,DIMENSION(nderivs) :: GETDERR != (/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    DOUBLE PRECISION,DIMENSION(nderivs) :: FNRDER              !Vector gives the values of the calculated derivatives
    integer:: i, oir
    DOUBLE PRECISION :: D_ARD, Rmix
    double precision, dimension(30):: chempot_res


    GETDERR = (/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    lnfi = 0.D0

    CALL MIXDERIVSFNR(gl,T,D,GETDERR,FNRDER)   !subroutine calculates derivatives of residual part

    D_ARD = FNRDER(2)

    oir = 2                    !The residual chemical potential will be calculated (0 = overall, 1 = ideal)
    call dna_dni(gl,T, D, chempot_res, OIR)
    call R_mix_calc(gl,Rmix)

    !Errorhandling, catch negative (or zero) values in the logarithm
    if ((D_ARD  <=  -1.D0) .or. (gl%molfractions(1)  <  1.D-24).or. (gl%molfractions(2)  <  1.D-24)) then
        return
    end if

    Do i = 1, gl%NCOMP
        lnfi(i) = dlog(gl%molfractions(i)*p) + chempot_res(i) - dlog(1.D0 + D_ARD)
    End Do

    END subroutine lnf_mix
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine dna_dni (gl,TEMPERATURE, DENSITY, CHEMPOT, OIR)
    ! SUBROUTINE FOR THE CALCULATION OF THE REDUCED CHEMICAL POTENTIALS. THE CALCULATION
    ! IS BASED ON THE FORMULATION OF THE HELMHOLTZ FREE ENERGY AS PUBLISHED BY:
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 58 Eq. 5.36, page 109 Table 7.4)
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! Temperature
    ! Density
    ! OIR       - 0: Overall Chemical potential is given    d(nF)/d(ni)
    !           - 1: Ideal chemical potential is given      d(nFi)/d(ni)
    !           - 2: Residual chemical potential is given   d(nFr)/d(ni)
    !
    ! OUTPUT PARAMETERS:
    ! CHEMPOT - n x 1 vector where n is the number of components in the mixture
    !-------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30), intent(out):: CHEMPOT
    integer:: OIR

    double precision, dimension(30):: chempot_ideal, chempot_res
    integer, dimension(nderivs):: getder_res
    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivsi):: getder_ideal
    double precision, dimension(nderivsi):: der_ideal
    double precision, dimension(30):: der_FRX

    !derivatives of the reducingfunctions with respect to ni multiplied by n
    double precision, dimension(30):: d_rhoredmix_d_ni, d_tredmix_d_ni
    double precision:: F0i_0, F_res, F_res_d_delta, F_res_d_tau, F_res_d_xi, F_res_d_x
    integer:: i, j


    getder_res = (/1,1,0,1,0,0,0,0,0,0,0,0,0,0,0/)                                !F_res, d_F_res_d_delta and d_F_res_d_tau are needed
    getder_ideal = (/1,0,0,0,0,0,0,0,0,0/)                                  !F_0 for each pure component is needed

    der_res = 0.D0    !Initialize
    der_ideal = 0.D0    !Initialize
    der_FRX = 0.D0
    CHEMPOT = 0.D0
    chempot_ideal = 0.D0
    chempot_res = 0.D0

    if (gl%NCOMP  ==  1) Then                                                    !PURE COMPONENT

        if (OIR  /=  2) then      !If only the residual part is needed, calculating the ideal part is not necessary
            Call FNIDERIVS (gl,Temperature, Density, getder_ideal, der_ideal,1)    !Get F of of the pure component
            F0i_0 = der_ideal(1)
            chempot_ideal(1) = 1.D0 + F0i_0                                        !Ideal chemical potential for pure component
        End if

        if (OIR  /=  1) then     !If only the ideal part is needed, calculating the residual part is not necessary
            Call FNRDERIVS (gl,Temperature, Density, getder_res, der_res, 1)
            F_res = der_res(1)
            F_res_d_delta = der_res(2)
            chempot_res(1) = F_res + F_res_d_delta                              !Residual chemical potential for pure component
        End if

    else
        !MIXTURE
        if (OIR  /=  2) then      !If only the residual part is needed, calculating the ideal part is not necessary
            Do i = 1, gl%NCOMP
                Call FNIDERIVS (gl,Temperature, Density, getder_ideal, der_ideal,i)    !Get F of the pure component
                F0i_0 = der_ideal(1)
                if (gl%molfractions(i)  >  0.D0) then
                    chempot_ideal(i) = 1.d0 + F0i_0 + dlog(gl%molfractions(i))             !Ideal chem.pot of comp i = F0i + ln(xi)
                end if
            end do
        End if

        if (OIR  /=  1) then     !If only the ideal part is needed, calculating the residual part is not necessary

            call ndYr_dni(gl,d_tredmix_d_ni, d_rhoredmix_d_ni)                       !Get derivatives of reduced parameters with respect to ni
            Call MIXDERIVSFNR (gl,Temperature, Density, getder_res, der_res)           !F_res and derivatives with respect to delta and tau
            Call dar_dxi (gl,Temperature, Density, der_FRX)                             !Derivative of F_res with respect to xi
            F_res = der_res(1)
            F_res_d_delta = der_res(2)
            F_res_d_tau = der_res(4)

            Do i = 1, gl%NCOMP
                chempot_res(i) = F_res + F_res_d_delta*(1.d0-1.d0/gl%rhoredmix*d_rhoredmix_d_ni(i))+ &
                    & F_res_d_tau/gl%tredmix*d_tredmix_d_ni(i)
                if (i  <  gl%NCOMP)then
                    F_res_d_xi = der_FRX(i)
                    chempot_res(i) = chempot_res(i) + F_res_d_xi
                end if
                Do j = 1, (gl%NCOMP-1)
                    F_res_d_x = der_FRX(j)
                    chempot_res(i) = chempot_res(i) - gl%molfractions(j) * F_res_d_x
                end do
            end do
        End if

    end if
    chempot = chempot_ideal + chempot_res

    end subroutine dna_dni
    !**************************************************************************

    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE PRESSURE WITH RESPECT TO VOLUME AT
    ! CONSTANT TEMPERATURE AND MOLE NUMBERS / CONSTANT TEMPERATURE AND MOLAR COMPOSITION.
    ! THE DERIVATIVE IS GIVEN TIMES THE TOTAL NUMBER OF MOLES n. AS PUBLISHED BY
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 108 table 7.2)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOLE/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dPdV        - Variable the derivative is stored in
    !--------------------------------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision:: dPdV

    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivs):: getder_res
    double precision:: Fr_d_delta, Fr_d_delta_d_delta, Rmix


    getder_res = (/0,1,1,0,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.D0
    DPDV = 0.D0                         !Initialize return variable
    Fr_d_delta = 0.D0                   !Initialize
    Fr_d_delta_d_delta = 0.D0           !Initialize


    Call MIXDERIVSFNR (gl,temperature, density, getder_res, der_res)
    Fr_d_delta = der_res(2)             !Derivative of the residual Helmholtz with respect to delta
    Fr_d_delta_d_delta = der_res(3)     !Second derivative of the residual Helmholtz with respect to delta

    Call R_mix_calc(gl,Rmix)
    dPdV = - Density**2 * Temperature *  Rmix * (1.d0 + 2.D0*Fr_d_delta + Fr_d_delta_d_delta)

    end subroutine ndP_dV

    !**************************************************************************
    subroutine dlnfi_dxj_TP (gl,TEMPERATURE, DENSITY, dlnfidXj)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE LOGARITHM OF THE FUGACITY
    ! WITH RESPECT TO THE MOLE FRACTION xj AT CONSTANT TEMPERATURE AND PRESSURE
    !                                           d(ln(fi))/d(xj) at const. T, P, xk
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - p   MPa
    !
    ! OUTPUT PARAMETERS:
    ! dlnfidXi       - MAtrix(30, 30) the derivatives for each component are stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30,30) :: dlnfidXj
    double precision, dimension(:, :), allocatable:: dlnphiidxj 
    
    double precision:: xj, xN
    integer:: i, j

    if (.not.allocated(dlnphiidxj)) allocate(dlnphiidxj(30,30))
    ! get the derivative d(ln(phii))/d(xj) at const. T, P
    call dlnphii_dxj_TP(gl,TEMPERATURE, DENSITY, dlnphiidxj)

    dlnfidXj = 0.D0
    do j = 1, gl%ncomp-1
        do i = 1, gl%ncomp-1
            if (j  /=  i) then
                dlnfidXj(j, i) = dlnphiidxj(j, i)
            else
                xj = gl%molfractions(j)
                dlnfidXj(j, i) = 1.D0/xj + dlnphiidxj(j, i)
            end if
        end do
        xN = gl%molfractions(gl%NCOMP)
        dlnfidXj(j,gl%NCOMP) = -1.d0/xN + dlnphiidxj(j, gl%NCOMP)
    end do

    end subroutine dlnfi_dxj_TP
    !**************************************************************************

    !**************************************************************************
    subroutine dlnphii_dP(gl,TEMPERATURE, DENSITY, dphiidP)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE FUGACITY COEFFICIENT
    ! OF COMPONENT i W.R.T. T AT CONSTANT T and x
    !                            --- d(ln(phi_i)/d(P)) ---
    ! J. Gernert, Dec. 2010
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! dphiidP     - A VECTOR OF LENGTH 30 IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: dphiidP
    double precision, dimension(30)::dPdni
    double precision:: dPdV, Rmix, Pressure
    integer::i

    ! get the derivative n*d(p)/d(ni) at const. T, V
    call ndP_dni_TV (gl,TEMPERATURE, DENSITY, dPdni)
    ! get the derivative n*d(p)/d(V) at const. T, x
    call ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    ! get the mixture Molar Gas Constant
    Call R_mix_calc(gl,Rmix)
    ! get the pressure
    Pressure = P_calc(gl,Temperature, Density, 0)

    dphiidP = 0.D0
    do i = 1, gl%ncomp
        dphiidP(i) = - dPdni(i)/(dPdV*Rmix*Temperature) - 1.D-6/Pressure
    end do

    end subroutine dlnphii_dP
    !**************************************************************************


    !**************************************************************************
    subroutine dlnphii_dT(gl,TEMPERATURE, DENSITY, dphiidT)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE FUGACITY COEFFICIENT
    ! OF COMPONENT i W.R.T. T AT CONSTANT P and x
    !                            --- d(ln(phi_i)/d(T)) ---
    ! J. Gernert, Dec. 2010
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! dphiidT     - A VECTOR OF LENGTH 30 IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: dphiidT

    double precision, dimension(30):: d2nadndT
    double precision, dimension(30):: dPdni
    double precision:: dPdV, dPdT, Rmix
    integer:: i

    ! get the derivative d^2(n*ar)/d(ni)d(T) at const. V, n
    call d2na_dnidT_V(gl,TEMPERATURE, DENSITY, d2nadndT)
    ! get the derivative n*d(p)/d(ni) at const. T, V
    call ndP_dni_TV (gl,TEMPERATURE, DENSITY, dPdni)
    ! get the derivative n*d(p)/d(V) at const. T, x
    call ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    ! get the derivative d(p)/d(T) at const. V, n
    call dP_dT (gl,TEMPERATURE, DENSITY, dPdT)
    ! get the mixture Molar Gas Constant
    Call R_mix_calc(gl,Rmix)

    dphiidT = 0.D0
    do i = 1, gl%ncomp
        dphiidT(i) = d2nadndT(i) + 1.D0/Temperature + dPdni(i)*dPdT/(dPdV*Rmix*Temperature)
    end do

    end subroutine dlnphii_dT
    !**************************************************************************


    !*********************************************************************************
    subroutine chempot_reac_TP(gl, T, d ,p , chempot, reactpos)
    !***********************************************************************
    !subroutine for the calculation of the chemical potential for reactive mixtures and gE models with helmholtz mixtures.
    !Benedikt 10/2018
    !---------------------------
    ! Input T: Temperature [K]
    !       d: density [mol/m^3]
    !       p: pressure [MPa] as needed for gibbs
    !output dchempotidxj: derivative of chemical potential of component i with respect to composition change of comp. j
    !       reactpos: here are the reactive components for which additive terms have to be considered.


    implicit none

    type(type_gl) :: gl


    double precision:: T, D, p, gibbsphase, sumnui
    integer :: oir
    double precision, dimension(30), intent(out):: chempot
    double precision, dimension(30) :: chempotge, dnadni

    double precision::  Rmix, wmmix

    integer:: i, nreac, reactpos

    !oir = 0                    !The overall chemical potential will be calculated (1 = ideal, 2 = residual)
    chempot = 0                !Initialize chempot
    chempotge = 0.d0
    call R_mix_calc(gl,Rmix)
    call wm_mix_calc(gl,wmmix)

    call dna_dni(gl,T, D, dnadni, 0)

    if(reactpos .eq. 1) then
        call chempot_reactive(gl, t, d, gl%sea%seap, chempotge, 1)
    end if
    sumnui = 0.d0
    Do i = 1, gl%NCOMP
        if((i .eq. 1) .and. (gl%seawater) .and. (reactpos .eq.1)) then
            chempot(i) = dnadni(i)* Rmix * T /gl%wm(i) + chempotge(i)*1.d0
            chempot(i) = chempot(i) *gl%wm(1) !
        else
            Chempot(i)= dnadni(i) * Rmix * T !/ gl%wm(i)
        end if


        write(*,*) chempot(1)

        sumnui = sumnui + gl%molfractions(i) * chempot(i)
    end do
    if(d .ge. 1.d4) then
        gl%seacalc = .true.
        gibbsphase = g_calc(gl, t, d, 0)
        gl%seacalc = .false.
        write(*,*) gibbsphase - sumnui, 'diff gibbs'
        write(*,*) (gl%molfractions(2) * chempot(2)) - gibbsphase
    end if
    end subroutine chempot_reac_TP
    !*********************************************************************************



    !********************************************************************************
    subroutine chempot_reactive(gl, t, d, p, chempotge, nreac)
    !*********************************************************
    !subroutine calculating the part of the chemical potential resulting form gibbs or gE model
    !----------------------------------------------
    !


    implicit none
    type(type_gl) :: gl
    double precision :: t, p, d, v, a, g, dgds, d2gdpds
    double precision, dimension(30), intent(out) :: chempotge
    double precision, dimension(30) :: dPdXi
    integer :: nreac

    chempotge = 0.d0
    if(gl%seawater) gl%seacalc = .true.
    if( (nreac .eq. 1) .and. (gl%seacalc) ) then

        !p = p_calc(gl, t, d, 0)
        !p = p_calc(gl, t, d, 1)
        p = gl%sea%seap
        v = v_saline(gl, t, p)
        a = a_saline(gl, t, p)
        g = g_saline(gl, t, p) - p*1.d6*v
        dgds = DGDS_saline(gl, t, p)
        chempotge(nreac) = -((gl%sea%salinity*  DGDS) -  (p * 1.d6 * v) - a )-g_saline(gl,t,p)
        chempotge(nreac) = chempotge(nreac) * (-10.d0)

        !chempotge(nreac) = 6500d0

        chempotge(2) = gl%sea%salinity
        gl%sea%salinity = gl%sea%salinity * (1.d0+1.d-6)
        chempotge(5) = gl%sea%salinity
        chempotge(3) = a_saline(gl, t, p)
        gl%sea%salinity = chempotge(2) * (1.d0 - 1.d-6)
        chempotge(6) = gl%sea%salinity
        chempotge(4) = a_saline(gl,t,p)
        chempotge(7)=  (chempotge(3)-chempotge(4))/(chempotge(5)-chempotge(6))
        gl%sea%salinity = chempotge(2)
        chempotge(1) = chempotge(7)
        gl%sea%salinity = chempotge(2)
        chempotge(nreac) = a !+  chempotge(7)

        !new part here (nach gerg)


        gl%seacalc = .false.
        call dP_dxi_TV (gl,t, d, dPdXi)
        gl%seacalc = .true.

        g = g_saline(gl, t, p)! - p*1.d6*v
        dgds = DGDS_saline(gl, t, p)
        d2gdpds = d2gdsdp_saline(gl, t, p)
        a = a_saline(gl, t, p)

        chempotge(nreac) = a + chempotge(7)*gl%sea%salinity! - p*d2gdpds - dpdxi(1) * v


        !end new part



        if(gl%sea%salinity .le. 1.d-9) then
            chempotge(nreac) = 0.d0
            gl%sea%salinity = 0.d0
        end if

    end if
    gl%seacalc = .false.
    end subroutine chempot_reactive
    !*******************************************************************************


    !**************************************************************************
    subroutine d2na_dnidT_V(gl,TEMPERATURE, DENSITY, dnadnidT)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF alpha_r W.R.T. ni AND
    ! the Temperature - d(d(n*ar)/d(ni))dT at const. V
    ! J. Gernert, Dec.. 2010
    !--------------------------------------------------------------------------------------------------
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 120 EQ. 7.44)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! dnadnidT     - A VECTOR OF 30 IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30)::dnadnidT
    double precision, dimension(30):: DANITAU
    double precision, dimension(nderivs):: MIXDERFNR
    integer, dimension(nderivs) :: GETDER
    double precision:: dar_dtau
    integer:: i

    dnadnidT = 0.D0
    ! get the derivative d(n*d(ar)/d(ni))/d(tau) at const. del, x
    call dndar_dnidtau (gl,TEMPERATURE, DENSITY, DANITAU)
    ! get the derivative tau*d(ar)/d(tau) at const. del, x
    GETDER = (/0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)
    dar_dtau = MIXDERFNR(4)

    do i = 1, gl%ncomp
        dnadnidT(i) = -1.d0/Temperature*(dar_dtau + DANITAU(i))
    end do

    end subroutine d2na_dnidT_V
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************

    subroutine dP_dT (gl,TEMPERATURE, DENSITY, dPdT)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE PRESSURE WITH RESPECT TO TEMPERATURE
    ! AT CONSTANT VOLUME AND MOLE NUMBERS / CONSTANT DENSITY AND MOLAR COMPOSITION AS PUBLISHED BY
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 108 table 7.2)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dPdT        - Variable the derivative is stored in    Unit: [Pa / K]
    !--------------------------------------------------------------------------------------------------



    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision:: dPdT

    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivs):: getder_res
    double precision:: Fr_d_delta, Fr_d_delta_d_tau, Rmix


    getder_res = (/0,1,0,0,0,1,0,0,0,0,0,0,0,0,0/)
    der_res = 0.D0
    DPDT = 0.D0                         !Initialize return variable
    Fr_d_delta = 0.D0                   !Initialize
    Fr_d_delta_d_tau = 0.D0             !Initialize

    Call MIXDERIVSFNR (gl,temperature, density, getder_res, der_res)
    Fr_d_delta = der_res(2)             !Derivative of the residual Helmholtz with respect to delta
    Fr_d_delta_d_tau = der_res(6)       !Second derivative of the residual Helmholtz with respect to delta and tau

    Call R_mix_calc(gl,Rmix)

    dPdT = Density *  Rmix * (1.d0 +  Fr_d_delta - Fr_d_delta_d_tau)

    end subroutine dP_dT
    !**************************************************************************

    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine d2na_dnidT_P (gl,TEMPERATURE, DENSITY, DCHPOTDT, OIR)
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE CHEMICAL POTENTIAL W.R.T.
    ! Temperature T at constant p and composition x.
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! Temperature
    ! Density
    ! OIR       - 0: Derivative of the overall chemical potential   d(d(nF)/d(ni))/dT
    !           - 1: Derivative of the ideal chemical potential     d(d(nFi)/d(ni)/dT
    !           - 2: Derivative of the residual chemical potential  d(d(nFr)/d(ni)/dT
    !
    ! OUTPUT PARAMETERS:
    ! DCHPOTDT - n x 1 vector where n is the number of components in the mixture
    !-------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DCHPOTDT
    integer:: OIR

    double precision, dimension(30):: dchpotdt_ideal, dchpotdt_res

    integer, dimension(nderivsi):: getder_ideal
    double precision, dimension(nderivsi):: der_ideal

    double precision:: ndaidV, daidT, daiddel, daidtau
    double precision, dimension(30):: ndnardnidV, dnardnidT
    double precision:: dVdT, dpdT, dpdV

    integer:: i
    !warnings (Moni)
    !integer:: j

    getder_ideal = (/0,1,0,0,1,0,0,0,0,0/)
    der_ideal = 0.D0             !Initialize
    DCHPOTDT = 0.d0

    call dP_dT (gl,TEMPERATURE, DENSITY, dPdT)
    call ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    dVdT = - dPdT / dPdV

    Do i = 1, gl%NCOMP
        Call FNIDERIVS (gl,Temperature, Density, getder_ideal, der_ideal,i)
        daiddel = der_ideal(2)
        daidtau = der_ideal(5)
        ndaidV = -1.d0 * Density * daiddel
        daidT = -1.d0 / Temperature * daidtau
        dchpotdt_ideal(i) = ndaidV * dVdT + daidT
    end do

    call nd2na_dnidV_T(gl,TEMPERATURE, DENSITY, ndnardnidV)
    call d2na_dnidT_V(gl,TEMPERATURE, DENSITY, dnardnidT)

    Do i = 1, gl%NCOMP
        dchpotdt_res(i) = ndnardnidV(i) * dVdT + dnardnidT(i)
    end do

    select case (OIR)
    case (0)
        Do i = 1, gl%NCOMP
            DCHPOTDT(i) = dchpotdt_ideal(i) + dchpotdt_res(i)                  !Overall chemical potential
        end do
    case (1)
        Do i = 1, gl%NCOMP
            DCHPOTDT(i) = dchpotdt_ideal(i)                                   !Ideal chemical potential
        end do
    case (2)
        Do i = 1, gl%NCOMP
            DCHPOTDT(i) = dchpotdt_res(i)                                     !Residual chemical potential
        end do
    end select

    end subroutine d2na_dnidT_P


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine d2na_dnidp_T (gl,TEMPERATURE, DENSITY, DCHPOTDP, OIR)
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE CHEMICAL POTENTIAL W.R.T.
    ! pressure p at constant T and composition x.
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! Temperature
    ! Density
    ! OIR       - 0: Derivative of the overall chemical potential   d(d(nF)/d(ni))/dp
    !           - 1: Derivative of the ideal chemical potential     d(d(nFi)/d(ni)/dp
    !           - 2: Derivative of the residual chemical potential  d(d(nFr)/d(ni)/dp
    !
    ! OUTPUT PARAMETERS:
    ! DCHPOTDT - n x 1 vector where n is the number of components in the mixture
    !-------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DCHPOTDP
    integer:: OIR
    double precision, dimension(30):: dchpotdp_ideal, dchpotdp_res
    integer, dimension(nderivsi):: getder_ideal
    double precision, dimension(nderivsi):: der_ideal
    double precision:: ndaidV, daiddel
    double precision, dimension(30):: ndnardnidV
    double precision:: dVdp, dpdV
    integer:: i
    !warnings (Moni)
    !integer:: j

    getder_ideal = (/0,1,0,0,0,0,0,0,0,0/)
    der_ideal = 0.D0              !Initialize
    DCHPOTDP = 0.d0

    call ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    dVdp = 1.D0 / dPdV

    Do i = 1, gl%NCOMP
        Call FNIDERIVS (gl,Temperature, Density, getder_ideal, der_ideal,i)
        daiddel = der_ideal(2)
        ndaidV = -1.d0 * Density * daiddel
        dchpotdp_ideal(i) = ndaidV * dVdp
    end do

    call nd2na_dnidV_T(gl,TEMPERATURE, DENSITY, ndnardnidV)

    Do i = 1, gl%NCOMP
        dchpotdp_res(i) = ndnardnidV(i) * dVdp
    end do

    select case (OIR)
    case (0)
        Do i = 1, gl%NCOMP
            DCHPOTDP(i) = dchpotdp_ideal(i) + dchpotdp_res(i)                  !Overall chemical potential
        end do
    case (1)
        Do i = 1, gl%NCOMP
            DCHPOTDP(i) = dchpotdp_ideal(i)                                   !Ideal chemical potential
        end do
    case (2)
        Do i = 1, gl%NCOMP
            DCHPOTDP(i) = dchpotdp_res(i)                                     !Residual chemical potential
        end do
    end select

    end subroutine d2na_dnidp_T


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine d2na_dnidxj_PT (gl,TEMPERATURE, DENSITY, DCHPOTDX, OIR)
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE CHEMICAL POTENTIAL W.R.T.
    ! the molfractions xj at constant p and T
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! Temperature
    ! Density
    ! OIR       - 0: Derivative of the overall chemical potential   d(d(nF)/d(ni))/dxj
    !           - 1: Derivative of the ideal chemical potential     d(d(nFi)/d(ni)/dxj
    !           - 2: Derivative of the residual chemical potential  d(d(nFr)/d(ni)/dxj
    !
    ! OUTPUT PARAMETERS:
    ! DCHPOTDX - n x n Matrix where n is the number of components in the mixture
    !-------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30, 30):: DCHPOTDX
    integer:: OIR
    double precision, dimension(:,:), allocatable :: dchpotdx_ideal, dchpotdx_res
    integer, dimension(nderivsi):: getder_ideal
    double precision, dimension(nderivsi):: der_ideal
    double precision:: ndaidV, daiddel
    double precision, dimension(30):: ndnardnidV, dPdxj, dVdxj
    double precision, dimension(:,:), allocatable :: d2nardnidxj
    double precision:: dpdV
    integer:: i, j

    if (.not.allocated(dchpotdx_ideal)) allocate(dchpotdx_res(30, 30),d2nardnidxj(30, 30),dchpotdx_ideal(30,30))
    
    DCHPOTDX = 0.D0
    dchpotdx_ideal = 0.D0
    dchpotdx_res = 0.D0

    !Monika, March 2017: not needed, see below
    !getder_ideal = (/0,1,0,0,0,0,0,0,0,0/)
    !der_ideal = 0.D0

    call ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    call dP_dxi_TV (gl,TEMPERATURE, DENSITY, dPdXj)

    Do i = 1, gl%Ncomp-1
        dVdxj(i) = - dPdxj(i) / dPdV
    end do

    Do j = 1, gl%NCOMP - 1
        Do i = 1, gl%NCOMP
            !Monika, March 2017: The 1st derivative of the ideal part wrt density is always 1, so that it is not needed to call this routine
            !Call FNIDERIVS (Temperature, Density, getder_ideal, der_ideal,i)
            !daiddel = der_ideal(2)
            !ndaidV = -1.d0 * Density * daiddel
            !dchpotdx_ideal(j,i) = ndaidV * dVdxj(j)

            dchpotdx_ideal(j,i) = -Density * dVdxj(j)
            if (j  ==  i) then
                dchpotdx_ideal(j,i) = dchpotdx_ideal(j,i) + 1.D0 / gl%molfractions(i)
            end if
        end do
        dchpotdx_ideal(j,gl%NCOMP) = dchpotdx_ideal(j,gl%NCOMP) - 1.D0 / gl%molfractions(gl%NCOMP)
    End do

    call nd2na_dnidV_T(gl,TEMPERATURE, DENSITY, ndnardnidV)
    call d2na_dnidxj_TV(gl,TEMPERATURE, DENSITY, d2nardnidxj)

    Do j = 1, gl%NCOMP - 1
        Do i = 1, gl%NCOMP
            dchpotdx_res(j,i) = ndnardnidV(i) * dVdxj(j) + d2nardnidxj(j,i)
        end do
    end do

    select case (OIR)
    case (0)
        DCHPOTDX = dchpotdx_ideal + dchpotdx_res                !Overall chemical potential
    case (1)
        DCHPOTDX = dchpotdx_ideal                               !Ideal chemical potential
    case (2)
        DCHPOTDX = dchpotdx_res                                 !Residual chemical potential
    end select

    end subroutine d2na_dnidxj_PT
    !********************************************************************************


    !*****************************************************************************
    subroutine d2na_dnidxj_PT_reac(gl, t,p, d, dchempot_reac_dx, reacpos)
    !subroutine for the calcoulation fo the derivative of the chemical potential wrt. tp xj, reactive part only
    !Benedikt 10/2018
    !-------------------------

    implicit none

    type(type_gl) :: gl

    double precision :: t, p, d, d2gds2, dgds, d2gds2_num
    integer :: reacpos, i, j
    double precision, dimension(30,30) :: dchempot_reac_dx

    dchempot_reac_dx = 0.d0

    Do j = 1, gl%NCOMP - 1
        Do i = 1, gl%NCOMP
            !dchpotdx_res(j,i) = ndnardnidV(i) * dVdxj(j) + d2nardnidxj(j,i)
            if((j .ne. 1.) .and. (gl%seacalc) .and. (i .eq. 1)) then
                !dchempot_reac_dx(j,i) = (g_saline(gl, t, p) - gl%sea%salinity * dgds_saline(gl, t, p) )*gl%wm(1)   !hier war noch was falsch?    !convert from spec to mol!!!!!!!
                dchempot_reac_dx(j,i) = 0.d0!(g_saline(gl, t, p) - gl%sea%salinity * dgds_saline(gl, t, p) )*gl%wm(1)
            elseif((j .eq. 1) .and. (i .eq. 1).and. (gl%seacalc)) then !vapor water
                !dchempot_reac_dx(j,i) = gl%sea%salinity * d2gds2_saline(gl, t, p) * gl%wm(1) !hier war noch was falsch?
                d2gds2 = d2gds2_saline(gl, t, p)!D2GDS2_saline(gl, t, p)
                dgds = dgds_saline(gl, t, p)
                d2gds2_num = d2gds2_saline_num(gl, t, p)!D2GDS2_saline(gl, t, p)
                dgds = dgds_saline(gl, t, p)
                dchempot_reac_dx(j,i) = -  (- gl%sea%salinity * d2gds2 - dgds - d2gds2 )!    gl%sea%salinity * d2gds2_saline(gl, t, p) * gl%wm(1)  !old version
                dchempot_reac_dx(j,i) = - dchempot_reac_dx(j,i) !* gl%wm(1)
            end if
        end do
    end do

    end subroutine d2na_dnidxj_PT_reac
    !***************************************************************************************


    !**************************************************************************
    subroutine nd2na_dnidV_T(gl,TEMPERATURE, DENSITY, ndnadnidV)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF alpha_r W.R.T. ni AND
    ! the Volume V - nd(d(n*ar)/d(ni))dV at const. T and compositon x
    ! A. Jäger, Mar.. 2011
    !--------------------------------------------------------------------------------------------------
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 120 EQ. 7.45)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! dnadnidV     - A VECTOR OF 30 IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30)::ndnadnidV
    double precision, dimension(30):: DANIDEL
    double precision, dimension(nderivs):: MIXDERFNR
    integer, dimension(nderivs) :: GETDER
    double precision:: dar_ddel
    integer:: i

    ndnadnidV = 0.D0
    ! get the derivative d(n*d(ar)/d(ni))/d(tau) at const. del, x
    call dndar_dniddel(gl,TEMPERATURE, DENSITY, DANIDEL)
    ! get the derivative tau*d(ar)/d(tau) at const. del, x
    GETDER = (/0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)
    dar_ddel = MIXDERFNR(2)

    do i = 1, gl%ncomp
        ndnadnidV(i) = -DENSITY*(dar_ddel + DANIDEL(i))
    end do

    end subroutine nd2na_dnidV_T
    !**************************************************************************

    !**************************************************************************
    !   Johannes Gernert
    !**************************************************************************

    double precision function dh_dT_px (gl,TEMPERATURE, DENSITY)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE ENTHALPY WITH RESPECT TO TEMPERATURE
    ! AT CONSTANT PRESSURE AND MOLE FRACTIONS
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dHdT        - Variable the derivative is stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision :: TEMPERATURE, DENSITY
    double precision :: dPdV, dPdT, dVdT, dhdV, dhdT_Vx, R
    double precision ::  ai_deltau, ai_tautau, ar_del, ar_deltau, ar_tautau, ar_deldel
    !warnings (Moni)
    !double precision :: ai_tau,ai_del, ar_tau
    double precision, dimension(nderivs):: FNRDER
    double precision, dimension(nderivsi)::FNIDER
    integer,dimension(nderivs):: GETDERFNR
    integer,dimension(nderivsi):: GETDERFNI

    call ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    call dP_dT (gl,TEMPERATURE, DENSITY, dPdT)

    GETDERFNI = (/0, 1, 0, 1, 1, 1,0,0,0,0/)
    call MIXDERIVSFNI(gl,TEMPERATURE, DENSITY, GETDERFNI, FNIDER)
    !warnings (Moni)
    !ai_del = FNIDER(2)
    !ai_tau = FNIDER(5)
    ai_deltau = FNIDER(4)
    ai_tautau = FNIDER(6)

    GETDERFNR = (/0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, GETDERFNR, FNRDER)
    ar_del = FNRDER(2)
    ar_deldel = FNRDER(3)
    !warnings (Moni)
    !ar_tau = FNRDER(4)
    ar_tautau = FNRDER(5)
    ar_deltau = FNRDER(6)

    call R_mix_calc(gl,R)

    dhdV = -Density*R*Temperature*(ai_deltau + ar_deltau + ar_del + ar_deldel)
    dhdT_Vx = R*(1.d0 + ar_del - ai_tautau - ar_tautau - ar_deltau)
    dVdT = - dPdT/dPdV

    dh_dT_px = dhdV*dVdT + dhdT_Vx

    end function dh_dT_px
    !**************************************************************************

    !**************************************************************************
    !   Johannes Gernert, November 2011
    !**************************************************************************

    subroutine dh_dx_TP (gl,TEMPERATURE, DENSITY, dhdx)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE ENTHALPY WITH RESPECT TO THE COMPOSITION
    ! OF THE COMPONENT i AT CONSTANT PRESSURE AND TEMPERATURE
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dHdT        - Variable the derivative is stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: dhdx
    double precision:: dPdV, dhdV, dh_ddel, dh_dtau, h, R
    double precision, dimension(30):: dVdx, dPdXi, d2a0dtaudx
    double precision:: ai_tau, ai_deltau, ai_tautau, ar_del, ar_tau, ar_deltau, ar_tautau, ar_deldel
    !warnings (Moni)
    !double precision:: ai_del, daiN_dtau
    double precision, dimension(nderivs):: FNRDER
    double precision, dimension(nderivsi)::FNIDER
    integer,dimension(nderivs):: GETDERFNR
    integer,dimension(nderivsi):: GETDERFNI
    double precision, dimension(30):: dTred_dxi, drhored_dxi, d2ardxiddel, d2ardxidtau
    integer:: i

    !warnings (Moni)
    !integer:: k

    dPdXi = 0.D0
    dVdx = 0.D0
    d2a0dtaudx = 0.D0
    dhdx = 0.d0

    call ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    call dP_dxi_TV (gl,TEMPERATURE, DENSITY, dPdXi)

    GETDERFNI = (/0, 1, 0, 1, 1, 1,0,0,0,0/)
    call MIXDERIVSFNI(gl,TEMPERATURE, DENSITY, GETDERFNI, FNIDER)
    !warnings
    !ai_del = FNIDER(2)
    ai_deltau = FNIDER(4)
    ai_tau = FNIDER(5)
    ai_tautau = FNIDER(6)

    GETDERFNR = (/0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, GETDERFNR, FNRDER)
    ar_del = FNRDER(2)
    ar_deldel = FNRDER(3)
    ar_tau = FNRDER(4)
    ar_tautau = FNRDER(5)
    ar_deltau = FNRDER(6)

    ! get the derivatives of the reducing parameters w.r.t. all x_i
    call dYr_dxi(gl,dTred_dxi, drhored_dxi)

    !get the derivatives of the residual part of alpha w.r.t. del and all x_i
    call d2ar_dxiddel (gl,TEMPERATURE, DENSITY, d2ardxiddel)

    !get the derivatives of the residual part of alpha w.r.t. tau and all x_i
    call d2ar_dxidtau (gl,TEMPERATURE, DENSITY, d2ardxidtau)

    !get the derivatives of the ideal-gas part of alpha w.r.t. tau and all x_i
    call d2a0_dtaudx (gl,TEMPERATURE, DENSITY, d2a0dtaudx)

    call R_mix_calc(gl,R)

    dhdV = -Density*R*Temperature*(ai_deltau + ar_deltau + ar_del + ar_deldel)
    dVdx = - dPdXi / dPdV
    dh_ddel = -R*Temperature*(ai_deltau + ar_deltau + ar_del + ar_deldel)   !dh_ddel multiplied by -1 and delta (coming from ddel_dxi)
    dh_dtau = R*Temperature*(-1.D0 - ar_del + ai_tautau + ar_tautau + ar_deltau) !dh_dtau multiplied by tau (coming from dtau_dxi)
    h = (1.d0 + ai_tau + ar_tau + ar_del)

    do i = 1, gl%ncomp - 1
        dhdx(i) = dhdV*dVdx(i) + dh_ddel*drhored_dxi(i)/gl%rhoredmix + dh_dtau*dTred_dxi(i)/gl%Tredmix &
            & + h*R*Temperature*dTred_dxi(i)/gl%Tredmix &
            & + R*Temperature*(d2a0dtaudx(i) + d2ardxidtau(i) + d2ardxiddel(i))
    end do

    end subroutine dh_dx_TP
    !**************************************************************************


    !**************************************************************************
    !   Johannes Gernert, November 2011
    !**************************************************************************

    subroutine d2a0_dtaudx (gl,TEMPERATURE, DENSITY, d2a0dtaudx)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE IDEAL-GAS PART OF THE HELMHOLTZ
    ! ENERGY WITH RESPECT TO THE COMPOSITION OF THE COMPONENT i AND del AT CONSTANT del AND tau.
    ! THIS DERIVATIVE IS GIVEN MULTIPLIED BY TAU: d2a0dtaudx * tau
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! d2a0dtaudx        - Variable the derivative is stored in (Vector, size 30)
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: d2a0dtaudx
    double precision, dimension(30):: dTred_dxi, drhored_dxi, da0i_dtau, da0i_dtau_ddel, da0i_dtau2
    double precision, dimension(nderivsi):: SETDERIVS
    integer, dimension(nderivsi):: GETDERIVS
    double precision:: xi
    !warnings (Moni)
    !double precision:: del, tau
    integer:: i, k

    ! get the derivatives of T_r and rho_r with respect to all x_i
    call dYr_dxi(gl,dTred_dxi, drhored_dxi)

    !warnings (Moni)
    !del = density / rhoredmix
    !tau = tredmix / temperature

    GETDERIVS = (/0,0,0,1,1,1,0,0,0,0/)
    SETDERIVS = 0.d0
    da0i_dtau = 0.d0
    da0i_dtau2 = 0.d0
    da0i_dtau_ddel = 0.d0
    !get all derivatives da0i_dtau, da0i_dtau^2, and da0i_ddel_dtau
    do i = 1, gl%ncomp
        call FNIDERIVS(gl,TEMPERATURE, DENSITY, GETDERIVS, SETDERIVS, i)
        da0i_dtau_ddel(i) = setderivs(4)
        da0i_dtau(i) = setderivs(5)
        da0i_dtau2(i) = setderivs(6)
    end do

    !2020-11-20, Andy: Do we have to consider different gas constants of the components in the mixture?
    d2a0dtaudx = 0.d0
    do k = 1, gl%ncomp-1
        d2a0dtaudx(k) = da0i_dtau(k) - da0i_dtau(gl%ncomp)
        do i = 1, gl%ncomp-1
            xi = gl%molfractions(i)
            d2a0dtaudx(k) = d2a0dtaudx(k) + xi*(da0i_dtau_ddel(i)*drhored_dxi(k)/gl%rhoredmix &
                & - da0i_dtau(i)*dTred_dxi(k)/gl%Tredmix - da0i_dtau2(i)*dTred_dxi(k)/gl%Tredmix &
                & - (da0i_dtau_ddel(gl%ncomp)*drhored_dxi(k)/gl%rhoredmix &
                & - da0i_dtau(gl%ncomp)*dTred_dxi(k)/gl%Tredmix - da0i_dtau2(gl%ncomp)*dTred_dxi(k)/gl%Tredmix))
        end do
        d2a0dtaudx(k) = d2a0dtaudx(k) + da0i_dtau_ddel(gl%ncomp)*drhored_dxi(k)/gl%rhoredmix &
            & - da0i_dtau(gl%ncomp)*dTred_dxi(k)/gl%Tredmix - da0i_dtau2(gl%ncomp)*dTred_dxi(k)/gl%Tredmix
    end do

    end subroutine d2a0_dtaudx
    !**************************************************************************

    !**************************************************************************
    !   Johannes Gernert, November 2011
    !**************************************************************************

    subroutine da0_dx (gl,TEMPERATURE, DENSITY, da0dx)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE IDEAL-GAS PART OF THE HELMHOLTZ
    ! ENERGY WITH RESPECT TO THE COMPOSITION OF THE COMPONENT AT CONSTANT del AND tau
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! da0dx        - Variable the derivative is stored in (Vector, size 30)
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: da0dx

    double precision, dimension(30):: dTred_dxi, drhored_dxi, da0i_dtau, a0i, da0i_ddel
    double precision, dimension(nderivsi):: SETDERIVS
    integer, dimension(nderivsi):: GETDERIVS
    double precision:: xi, xk, xn
    !warnings (Moni)
    !double precision:: del,tau
    integer:: i, k

    ! get the derivatives of T_r and rho_r with respect to all x_i
    call dYr_dxi(gl,dTred_dxi, drhored_dxi)

    !warnings (Moni)
    !del = density / rhoredmix
    !tau = tredmix / temperature

    GETDERIVS = (/1,1,0,0,1,0,0,0,0,0/)
    SETDERIVS = 0.d0
    da0i_dtau = 0.d0
    da0i_ddel = 0.d0
    a0i = 0.d0

    !get all derivatives da0i_dtau, da0i_dtau^2, and da0i_ddel_dtau
    do i = 1, gl%ncomp
        call FNIDERIVS(gl,TEMPERATURE, DENSITY, GETDERIVS, SETDERIVS, i)
        da0i_dtau(i) = setderivs(5)
        da0i_ddel(i) = setderivs(2)
        a0i(i) = setderivs(1)
    end do

    
    !2020-11-20, Andy: Do we have to consider different gas constants of the components in the mixture?
    xn = gl%molfractions(gl%ncomp)
    da0dx = 0.d0
    do k = 1, gl%ncomp-1
        xk = gl%molfractions(k)
        da0dx(k) = a0i(k) - a0i(gl%ncomp) + dlog(xk/xn) - da0i_dtau(gl%ncomp)*dTred_dxi(k)/gl%Tredmix &
            & + da0i_ddel(gl%ncomp)*drhored_dxi(k)/gl%rhoredmix
        do i = 1, gl%ncomp-1
            xi = gl%molfractions(i)
            da0dx(k) = da0dx(k) + xi*(-da0i_dtau(i)*dTred_dxi(k)/gl%Tredmix + da0i_ddel(i)*drhored_dxi(k)/gl%rhoredmix &
                & + da0i_dtau(gl%ncomp)*dTred_dxi(k)/gl%Tredmix - da0i_ddel(gl%ncomp)*drhored_dxi(k)/gl%rhoredmix)
        end do
    end do

    end subroutine da0_dx
    !**************************************************************************


    !**************************************************************************
    !   Johannes Gernert, November 2011
    !**************************************************************************
    double precision function ds_dT_px (gl,TEMPERATURE, DENSITY)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE ENTROPY WITH RESPECT TO TEMPERATURE
    ! AT CONSTANT PRESSURE AND MOLE FRACTIONS
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! ds_dT_px        - Variable the derivative is stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision :: TEMPERATURE, DENSITY
    double precision :: dPdV, dPdT, dVdT, dsdV, dsdT_Vx, R
    double precision :: ai_del, ai_tautau, ar_del, ar_tautau, ar_deltau
    double precision, dimension(nderivs):: FNRDER
    double precision, dimension(nderivsi)::FNIDER
    integer,dimension(nderivs):: GETDERFNR
    integer,dimension(nderivsi):: GETDERFNI

    call ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    call dP_dT (gl,TEMPERATURE, DENSITY, dPdT)

    GETDERFNI = (/0, 1, 0, 0, 0, 1,0,0,0,0/)
    call MIXDERIVSFNI(gl,TEMPERATURE, DENSITY, GETDERFNI, FNIDER)
    ai_del = FNIDER(2)
    ai_tautau = FNIDER(6)

    GETDERFNR = (/0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, GETDERFNR, FNRDER)
    ar_del = FNRDER(2)
    ar_deltau = FNRDER(6)
    ar_tautau = FNRDER(5)

    call R_mix_calc(gl,R)

    dsdV = -R*Density*(ar_deltau - ai_del - ar_del) !multiplied by -delta*dens coming from ddel_dV
    dsdT_Vx = - R/Temperature*(ai_tautau + ar_tautau) !multiplied by -tau/T coming from dtau_dT
    dVdT = - dPdT/dPdV

    ds_dT_px = dsdV*dVdT + dsdT_Vx

    end function ds_dT_px
    !**************************************************************************

    !**************************************************************************
    !   Johannes Gernert, November 2011
    !**************************************************************************
    subroutine ds_dx_TP (gl,TEMPERATURE, DENSITY, dsdx)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE ENTROPY WITH RESPECT TO THE COMPOSITION
    ! OF THE COMPONENT i AT CONSTANT PRESSURE AND TEMPERATURE
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dsdx        - Variable the derivative is stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: dsdx

    double precision:: dPdV, dsdV, ds_ddel, ds_dtau, R
    double precision, dimension(30):: dVdx, dPdXi
    double precision:: ai_del, ai_tautau, ar_del, ar_deltau, ar_tautau
    double precision, dimension(nderivs):: FNRDER
    double precision, dimension(nderivsi)::FNIDER
    integer,dimension(nderivs):: GETDERFNR
    integer,dimension(nderivsi):: GETDERFNI
    double precision, dimension(30):: dTred_dxi, drhored_dxi, dardxi, d2ardxidtau, da0dx, d2a0dtaudx
    integer:: i


    dPdXi = 0.D0
    dVdx = 0.D0
    d2a0dtaudx = 0.D0
    dsdx = 0.d0

    call ndP_dV (gl,TEMPERATURE, DENSITY, dPdV)
    call dP_dxi_TV (gl,TEMPERATURE, DENSITY, dPdXi)

    GETDERFNI = (/0, 1, 0, 0, 1, 1,0,0,0,0/)
    call MIXDERIVSFNI(gl,TEMPERATURE, DENSITY, GETDERFNI, FNIDER)
    ai_del = FNIDER(2)
    ai_tautau = FNIDER(6)

    GETDERFNR = (/0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, GETDERFNR, FNRDER)
    ar_del = FNRDER(2)
    ar_tautau = FNRDER(5)
    ar_deltau = FNRDER(6)

    ! get the derivatives of the reducing parameters w.r.t. all x_i
    call dYr_dxi(gl,dTred_dxi, drhored_dxi)

    !get the derivatives of the residual part of alpha w.r.t. all x_i
    call dar_dxi (gl,TEMPERATURE, DENSITY, dardxi)

    !get the derivatives of the residual part of alpha w.r.t. tau and all x_i
    call d2ar_dxidtau (gl,TEMPERATURE, DENSITY, d2ardxidtau)

    !get the derivatives of the ideal-gas part of alpha w.r.t. all x_i
    call da0_dx (gl,TEMPERATURE, DENSITY, da0dx)

    !get the derivatives of the ideal-gas part of alpha w.r.t. tau and all x_i
    call d2a0_dtaudx (gl,TEMPERATURE, DENSITY, d2a0dtaudx)

    call R_mix_calc(gl,R)

    dsdV = -R*Density*(ar_deltau - ai_del - ar_del)
    dVdx = - dPdXi / dPdV
    ds_ddel = -R*(ar_deltau - ai_del - ar_del)  !Multiplied by -delta coming from drhoredmix_dxi
    ds_dtau = R*(ai_tautau + ar_tautau) !Multiplied by tau coming from dtau_dxi

    do i = 1, gl%ncomp - 1
        dsdx(i) = dsdV*dVdx(i) + ds_ddel*drhored_dxi(i)/gl%rhoredmix + ds_dtau*dTred_dxi(i)/gl%Tredmix &
            & +R*(d2a0dtaudx(i) + d2ardxidtau(i) - da0dx(i) - dardxi(i))
    end do

    end subroutine ds_dx_TP
    !**************************************************************************



    !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    !NEW DERIVATIVES NEEDED FOR CRITICAL POINT CALCULATION OF A MIXTURE
    !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

    !**************************************************************************
    subroutine ndlnfi_dnj(gl,Temp, Dens, ndlnfidnj)    !CONSTANT T, V, and nm
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE LOGARITHM OF
    ! THE FUGACITY OF COMPONENT i IN THE MIXTURE WITH RESPECT
    ! TO MOLE NUMBER nj MULTIPLIED BY TOTAL NUMBER OF MOLES n
    ! AT CONSTANT T AND V AND nm.
    ! n*d(lnfi)/dnj)_TV
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 1, Eq. (5)
    !
    !These derivatives are needed to built matrix L (page 3, Eq. 4) and the
    !rows 1...ncomp-1 (all except for the last row) of matrix M (page3, Eq. 6)
    !
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMP - Temperature in   K
    ! DENS - Density in   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! ndlnfidnj    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas, Feb 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens
    Double precision, dimension(30,30):: ndlnfidnj

    Double precision, dimension(30,30):: nd2nar_dnidnj
    integer:: i,j

    call nddna_dnidnj(gl,Temp, dens, nd2nar_dnidnj)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            if (i == j) then
                ndlnfidnj(i,j) = (nd2nar_dnidnj(i,j) + 1.D0 / gl%molfractions(i))
            else
                ndlnfidnj(i,j) = nd2nar_dnidnj(i,j)
            end if
        end do
    end do

    end subroutine
    !**************************************************************************


    !**************************************************************************
    subroutine n2d2lnfi_dnjdnk(gl,Temp, Dens) !CONSTANT T, V, and nm
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE LOGARITHM
    ! OF THE FUGACITY OF COMPONENT i WITH RESPECT TO MOLE NUMBERS nj AND nk
    ! AT CONSTANT T, V and other mole numbers
    ! THE WHOLE TERM IS MULTIPLIED BY TOTAL NUMBER OF MOLES n^2
    ! n^2 * d(lnfi)/d(nj)d(nk) at const. T, V, nm
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 1, Eq. (9)
    !
    !These derivatives are needed to fill the last row of matrix M
    ! (page3, Eq. 6 of the main article)
    !
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMP - Temperature in   K
    ! DENS - Density in   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! n2d2lnfidnjdnk    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas, Feb 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens


    integer:: i,j,k,m


    !Other derivatives needed to calculate the derivative n * d/dnk (n*d(lnfi)/dnj)_TV)
    double precision, dimension(30,30):: ndlnfidnj

    if(.not. allocated(gl%n2d2lnfidnjdnk)) allocate(gl%n2d2lnfidnjdnk(30,30,30))
    if(.not. allocated(gl%ndndlnfidnjdnk)) allocate(gl%ndndlnfidnjdnk(30,30,30))

    gl%n2d2lnfidnjdnk = 0.D0

    call ndlnfi_dnj(gl,Temp, Dens, ndlnfidnj)
    call ndndlnfi_dnjdnk(gl,Temp, Dens)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            do k = 1, gl%ncomp
                gl%n2d2lnfidnjdnk(i,j,k) = gl%ndndlnfidnjdnk(i,j,k) - ndlnfidnj(i,j)
            end do
        end do
    end do

    end subroutine n2d2lnfi_dnjdnk
    !**************************************************************************


    !**************************************************************************
    subroutine dn2d2lnfi_dnjdnkddel(gl,Temp, Dens) !CONSTANT T, V, and nm
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF
    ! n^2 * d(lnfi)/d(nj)d(nk) at const. T, V, nm
    ! WRT DELTA AT CONSTANT TAU AND X
    ! TERM IS MULTIPLIED BY DELTA!!
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 2, Eq. (13)
    !
    !
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMP - Temperature in   K
    ! DENS - Density in   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dn2d2lnfidnjdnkddel    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens

    integer:: i,j,k,m


    !Other derivatives needed to calculate the derivative n * d/dnk (n*d(lnfi)/dnj)_TV)
    double precision, dimension(30,30):: dndlnfidnjddel



    if(.not. allocated(gl%dn2d2lnfidnjdnkddel)) allocate(gl%dn2d2lnfidnjdnkddel(30,30,30))
    if(.not. allocated(gl%dndndlnfidnjdnkddel)) allocate(gl%dndndlnfidnjdnkddel(30,30,30))

    gl%dn2d2lnfidnjdnkddel = 0.D0

    call dndlnfi_dnjddel(gl,Temp, Dens, dndlnfidnjddel)
    call dndndlnfi_dnjdnkddel(gl,Temp, Dens)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            do k = 1, gl%ncomp
                gl%dn2d2lnfidnjdnkddel(i,j,k) = gl%dndndlnfidnjdnkddel(i,j,k) - dndlnfidnjddel(i,j)
            end do
        end do
    end do

    end subroutine dn2d2lnfi_dnjdnkddel
    !**************************************************************************



    !**************************************************************************
    subroutine dn2d2lnfi_dnjdnkdtau(gl,Temp, Dens) !CONSTANT T, V, and nm
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF
    ! n^2 * d(lnfi)/d(nj)d(nk) at const. T, V, nm
    ! WRT TAU AT CONSTANT DEL AND X
    ! TERM IS MULTIPLIED BY TAU!!
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 2, Eq. (14)
    !
    !
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMP - Temperature in   K
    ! DENS - Density in   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dn2d2lnfidnjdnkdtau    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens

    integer:: i,j,k,m


    !Other derivatives needed to calculate the derivative n * d/dnk (n*d(lnfi)/dnj)_TV)
    double precision, dimension(30,30):: dndlnfidnjdtau


    if(.not. allocated(gl%dn2d2lnfidnjdnkdtau)) allocate(gl%dn2d2lnfidnjdnkdtau(30,30,30))
    if(.not. allocated(gl%dndndlnfidnjdnkdtau)) allocate(gl%dndndlnfidnjdnkdtau(30,30,30))

    gl%dn2d2lnfidnjdnkdtau = 0.D0

    call dndlnfi_dnjdtau(gl,Temp, Dens, dndlnfidnjdtau)
    call dndndlnfi_dnjdnkdtau(gl,Temp, Dens)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            do k = 1, gl%ncomp
                gl%dn2d2lnfidnjdnkdtau(i,j,k) = gl%dndndlnfidnjdnkdtau(i,j,k) - dndlnfidnjdtau(i,j)
            end do
        end do
    end do

    end subroutine dn2d2lnfi_dnjdnkdtau
    !**************************************************************************



    !**************************************************************************
    subroutine ndndlnfi_dnjdnk(gl,Temp, Dens) !CONSTANT T, V, and nm
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE DERIVATIVE OF
    ! THE LOGARITHM OF THE FUGACITY OF COMPONENT i IN THE MIXTURE WITH RESPECT
    ! TO MOLE NUMBER nj MULTIPLIED BY TOTAL NUMBER OF MOLES n WITH RESPECT TO
    ! MOLE NUMBER nk AT CONSTANT T AND V AND nm.
    ! THE WHOLE TERM IS MULTIPLIED BY TOTAL NUMBER OF MOLES n
    ! n * d/dnk (n*d(lnfi)/dnj)_TV)
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 2, Eq. (10)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMP - Temperature in   K
    ! DENS - Density in   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! ndndlnfidnjdnk    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas, Feb 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens
    integer:: i,j,k,m

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative n * d/dnk (n*d(lnfi)/dnj)_TV)
    Double precision, dimension(30):: nddel_dnk, ndtau_dnk
    Double precision, dimension(30,30):: d_n_dlnfi_dnj_dtau, d_n_dlnfi_dnj_ddel
    double precision, dimension(30):: ndTred_dnk, ndrhored_dnk



    if(.not. allocated(gl%ndndlnfidnjdnk)) allocate(gl%ndndlnfidnjdnk(30,30,30))
    if(.not. allocated(gl%d_n_dlnfi_dnj_dxm)) allocate(gl%d_n_dlnfi_dnj_dxm(30,30,30))

    gl%ndndlnfidnjdnk = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    !Calculate the derivatives nddel_dnk and ndtau_dnk
    ! get derivatives of the reducing parameters Tred and rhored w.r.t all ni
    call ndYr_dni(gl,ndTred_dnk, ndrhored_dnk)
    ndtau_dnk = ndTred_dnk / gl%tredmix                   !Divided by tau, because all tau derivatives are already multiplied by tau
    nddel_dnk = 1.D0 - ndrhored_dnk / gl%rhoredmix        !Divided by del, because all del derivatives are already multiplied by del

    !Get tau-derivative of n*d(lnfi)/dnj
    call dndlnfi_dnjdtau(gl,Temp, Dens, d_n_dlnfi_dnj_dtau)
    !Get del-derivative of n*d(lnfi)/dnj
    call dndlnfi_dnjddel(gl,Temp, Dens, d_n_dlnfi_dnj_ddel)
    !Get x-derivatives of n*d(lnfi)/dnj
    call dndlnfi_dnjdxk(gl,Temp, Dens)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            do k = 1, gl%ncomp
                gl%ndndlnfidnjdnk(i,j,k) = d_n_dlnfi_dnj_dtau(i,j) * ndtau_dnk(k) + d_n_dlnfi_dnj_ddel(i,j) * nddel_dnk(k)
                do m= 1, gl%ncomp-1
                    if (k == m) then
                        gl%ndndlnfidnjdnk(i,j,k) = gl%ndndlnfidnjdnk(i,j,k) + gl%d_n_dlnfi_dnj_dxm(i,j,m) * (1.D0 - gl%molfractions(m))
                    else
                        gl%ndndlnfidnjdnk(i,j,k) = gl%ndndlnfidnjdnk(i,j,k) + gl%d_n_dlnfi_dnj_dxm(i,j,m) * (-gl%molfractions(m))
                    end if
                end do
            end do
        end do
    end do


    end subroutine ndndlnfi_dnjdnk
    !**************************************************************************


    !**************************************************************************
    subroutine dndndlnfi_dnjdnkdtau (gl,Temp, Dens) !CONSTANT T, V, and nm
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE TERM
    ! n * d/dnk (n*d(lnfi)/dnj)_TV)
    ! WITH RESPECT TO TAU AT CONSTANT DELTA AND X
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 2, Eq. (15)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMP - Temperature in   K
    ! DENS - Density in   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dndndlnfidnjdnkdtau    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens

    integer:: i,j,k,m

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative n * d/dnk (n*d(lnfi)/dnj)_TV)
    Double precision, dimension(30):: nddel_dnk, ndtau_dnk, dnddel_dnkddel, dndtau_dnkdtau
    Double precision, dimension(30,30):: d_n_dlnfi_dnj_dtau, d2_n_dlnfi_dnj_dtauddel, d2_n_dlnfi_dnj_dtau2
    double precision, dimension(30):: ndTred_dnk, ndrhored_dnk


    if(.not. allocated(gl%dndndlnfidnjdnkdtau)) allocate(gl%dndndlnfidnjdnkdtau(30,30,30))
    if(.not. allocated(gl%d2_n_dlnfi_dnj_dxmdtau)) allocate(gl%d2_n_dlnfi_dnj_dxmdtau(30,30,30))


    gl%dndndlnfidnjdnkdtau = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    !Calculate the derivatives nddel_dnk and ndtau_dnk
    ! get derivatives of the reducing parameters Tred and rhored w.r.t all ni
    call ndYr_dni(gl,ndTred_dnk, ndrhored_dnk)
    ndtau_dnk = ndTred_dnk / gl%tredmix                   !Divided by tau, because all tau derivatives are already multiplied by tau
    nddel_dnk = 1.D0 - ndrhored_dnk / gl%rhoredmix        !Divided by del, because all del derivatives are already multiplied by del
    dndtau_dnkdtau = (ndTred_dnk / gl%tredmix)! / tau                   !Divided by tau, because all tau derivatives are already multiplied by tau
    dnddel_dnkddel = (1.D0 - ndrhored_dnk / gl%rhoredmix) !/ delta      !Divided by del, because all del derivatives are already multiplied by del

    !Get tau-derivative of n*d(lnfi)/dnj
    call dndlnfi_dnjdtau(gl,Temp, Dens, d_n_dlnfi_dnj_dtau)
    !Get tau-del-derivative of n*d(lnfi)/dnj
    call d2ndlnfi_dnjdtauddel(gl,Temp, Dens, d2_n_dlnfi_dnj_dtauddel)
    !Get tau-tau-derivative of n*d(lnfi)/dnj
    call d2ndlnfi_dnjdtau2(gl,Temp, Dens, d2_n_dlnfi_dnj_dtau2)
    !Get tau-x-derivatives of n*d(lnfi)/dnj
    call d2ndlnfi_dnjdxkdtau(gl,Temp, Dens)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            do k = 1, gl%ncomp
                gl%dndndlnfidnjdnkdtau(i,j,k) = d2_n_dlnfi_dnj_dtau2(i,j) * ndtau_dnk(k) &
                    & + d_n_dlnfi_dnj_dtau(i,j) * dndtau_dnkdtau(k) &
                    & + d2_n_dlnfi_dnj_dtauddel(i,j) * nddel_dnk(k)
                do m= 1, gl%ncomp-1
                    if (k == m) then
                        gl%dndndlnfidnjdnkdtau(i,j,k) = gl%dndndlnfidnjdnkdtau(i,j,k) + gl%d2_n_dlnfi_dnj_dxmdtau(i,j,m) * (1.D0 - gl%molfractions(m))
                    else
                        gl%dndndlnfidnjdnkdtau(i,j,k) = gl%dndndlnfidnjdnkdtau(i,j,k) + gl%d2_n_dlnfi_dnj_dxmdtau(i,j,m) * (-gl%molfractions(m))
                    end if
                end do
            end do
        end do
    end do


    end subroutine dndndlnfi_dnjdnkdtau
    !**************************************************************************



    !**************************************************************************
    subroutine dndndlnfi_dnjdnkddel (gl,Temp, Dens) !CONSTANT T, V, and nm
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE TERM
    ! n * d/dnk (n*d(lnfi)/dnj)_TV)
    ! WITH RESPECT TO DELTA AT CONSTANT TAU AND X
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 2, Eq. (16)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMP - Temperature in   K
    ! DENS - Density in   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dndndlnfidnjdnkddel    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens

    integer:: i,j,k,m

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative n * d/dnk (n*d(lnfi)/dnj)_TV)
    Double precision, dimension(30):: nddel_dnk, ndtau_dnk, dnddel_dnkddel, dndtau_dnkdtau
    Double precision, dimension(30,30):: d_n_dlnfi_dnj_ddel, d2_n_dlnfi_dnj_dtauddel, d2_n_dlnfi_dnj_ddel2
    double precision, dimension(30):: ndTred_dnk, ndrhored_dnk


    if(.not. allocated(gl%dndndlnfidnjdnkddel)) allocate(gl%dndndlnfidnjdnkddel(30,30,30))
    if(.not. allocated(gl%d2_n_dlnfi_dnj_dxmddel)) allocate(gl%d2_n_dlnfi_dnj_dxmddel(30,30,30))

    gl%dndndlnfidnjdnkddel = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    !Calculate the derivatives nddel_dnk and ndtau_dnk
    ! get derivatives of the reducing parameters Tred and rhored w.r.t all ni
    call ndYr_dni(gl,ndTred_dnk, ndrhored_dnk)
    ndtau_dnk = ndTred_dnk / gl%tredmix                   !Divided by tau, because all tau derivatives are already multiplied by tau
    nddel_dnk = 1.D0 - ndrhored_dnk / gl%rhoredmix        !Divided by del, because all del derivatives are already multiplied by del
    dndtau_dnkdtau = (ndTred_dnk / gl%tredmix)! / tau                   !Divided by tau, because all tau derivatives are already multiplied by tau
    dnddel_dnkddel = (1.D0 - ndrhored_dnk / gl%rhoredmix)! / delta      !Divided by del, because all del derivatives are already multiplied by del

    !Get del-derivative of n*d(lnfi)/dnj
    call dndlnfi_dnjddel(gl,Temp, Dens, d_n_dlnfi_dnj_ddel)
    !Get tau-del-derivative of n*d(lnfi)/dnj
    call d2ndlnfi_dnjdtauddel(gl,Temp, Dens, d2_n_dlnfi_dnj_dtauddel)
    !Get del-del-derivative of n*d(lnfi)/dnj
    call d2ndlnfi_dnjddel2(gl,Temp, Dens, d2_n_dlnfi_dnj_ddel2)
    !Get tau-x-derivatives of n*d(lnfi)/dnj
    call d2ndlnfi_dnjdxkddel(gl,Temp, Dens)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            do k = 1, gl%ncomp
                gl%dndndlnfidnjdnkddel(i,j,k) = d2_n_dlnfi_dnj_dtauddel(i,j) * ndtau_dnk(k) &
                    & + d2_n_dlnfi_dnj_ddel2(i,j) * nddel_dnk(k) &
                    & + d_n_dlnfi_dnj_ddel(i,j) * dnddel_dnkddel(k)
                do m= 1, gl%ncomp-1
                    if (k == m) then
                        gl%dndndlnfidnjdnkddel(i,j,k) = gl%dndndlnfidnjdnkddel(i,j,k) + gl%d2_n_dlnfi_dnj_dxmddel(i,j,m) * (1.D0 - gl%molfractions(m))
                    else
                        gl%dndndlnfidnjdnkddel(i,j,k) = gl%dndndlnfidnjdnkddel(i,j,k) + gl%d2_n_dlnfi_dnj_dxmddel(i,j,m) * (-gl%molfractions(m))
                    end if
                end do
            end do
        end do
    end do


    end subroutine dndndlnfi_dnjdnkddel
    !**************************************************************************



    !**************************************************************************
    subroutine dndlnfi_dnjdtau(gl,Temp, Dens, dndlnfidnjdtau)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE DERIVATIVE OF
    ! THE LOGARITHM OF THE FUGACITY OF COMPONENT i IN THE MIXTURE WITH RESPECT
    ! TO MOLE NUMBER nj MULTIPLIED BY TOTAL NUMBER OF MOLES n WITH RESPECT TO
    ! TAU AT CONSTANT del AND x.
    ! THE RESULT IS MULTIPLIED WITH TAU!!!
    ! d/dtau (n*d(lnfi)/dnj)_del,x)
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 3, Eq. (20)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dndlnfidnjdtau    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas, Feb 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens
    Double precision, dimension(30,30):: dndlnfidnjdtau

    integer:: i,j

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative d/dtau (n*d(lnfi)/dnj)_del,x)
    Double precision, dimension(15,30)::dnardnjALL
    integer, dimension(15):: GETDER

    if(.not. allocated(gl%dndnardnidnjALL))allocate(gl%dndnardnidnjALL(15,30,30))

    dndlnfidnjdtau = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    GETDER = 0
    !Get the tau derivative of ndar_dnj
    GETDER(4) = 1
    call ndar_dni_ALL(gl,Temp, dens, GETDER, dnardnjALL)
    !Get the tau derivative of ndnar_dnidnj
    call ndndar_dnidnj_ALL(gl,Temp, dens, GETDER)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            dndlnfidnjdtau(i,j) = dnardnjALL(4,j) + gl%dndnardnidnjALL(4,i,j)
        end do
    end do

    end subroutine dndlnfi_dnjdtau
    !**************************************************************************



    !**************************************************************************
    subroutine dndlnfi_dnjddel(gl,Temp, Dens, dndlnfidnjddel)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE DERIVATIVE OF
    ! THE LOGARITHM OF THE FUGACITY OF COMPONENT i IN THE MIXTURE WITH RESPECT
    ! TO MOLE NUMBER nj MULTIPLIED BY TOTAL NUMBER OF MOLES n WITH RESPECT TO
    ! DEL AT CONSTANT tau AND x.
    ! THE RESULT IS MULTIPLIED WITH DEL!!!
    ! d/ddel (n*d(lnfi)/dnj)_tau,x)
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 3, Eq. (21)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dndlnfidnjddel    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens
    Double precision, dimension(30,30):: dndlnfidnjddel

    integer:: i,j

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative d/dtau (n*d(lnfi)/dnj)_del,x)
    Double precision, dimension(15,30)::dnardnjALL
    integer, dimension(15):: GETDER

    if(.not. allocated(gl%dndnardnidnjALL))allocate(gl%dndnardnidnjALL(15,30,30))

    dndlnfidnjddel = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    GETDER = 0
    !Get the del derivative of ndar_dnj
    GETDER(2) = 1
    call ndar_dni_ALL(gl,Temp, dens, GETDER, dnardnjALL)
    !Get the del derivative of ndnar_dnidnj
    call ndndar_dnidnj_ALL(gl,Temp, dens, GETDER)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            dndlnfidnjddel(i,j) = dnardnjALL(2,j) + gl%dndnardnidnjALL(2,i,j)
        end do
    end do

    end subroutine dndlnfi_dnjddel
    !**************************************************************************


    !**************************************************************************
    subroutine d2ndlnfi_dnjdtauddel(gl,Temp, Dens, d2ndlnfidnjdtauddel)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE DERIVATIVE OF
    ! THE LOGARITHM OF THE FUGACITY OF COMPONENT i IN THE MIXTURE WITH RESPECT
    ! TO MOLE NUMBER nj MULTIPLIED BY TOTAL NUMBER OF MOLES n WITH RESPECT TO
    ! DEL and TAU AT CONSTANT tau AND x.
    ! THE RESULT IS MULTIPLIED WITH DEL AND TAU!!!
    ! d2/ddel dtau (n*d(lnfi)/dnj)_tau,x)
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 4, Eq. (25)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! d2ndlnfidnjdtauddel    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens
    Double precision, dimension(30,30):: d2ndlnfidnjdtauddel

    integer:: i,j

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative d/dtau (n*d(lnfi)/dnj)_del,x)
    Double precision, dimension(15,30)::dnardnjALL
    integer, dimension(15):: GETDER


    d2ndlnfidnjdtauddel = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    GETDER = 0
    !Get the del and tau derivative of ndar_dnj
    GETDER(6) = 1
    call ndar_dni_ALL(gl,Temp, dens, GETDER, dnardnjALL)
    !Get the del derivative of ndnar_dnidnj
    call ndndar_dnidnj_ALL(gl,Temp, dens, GETDER)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            d2ndlnfidnjdtauddel(i,j) = dnardnjALL(6,j) + gl%dndnardnidnjALL(6,i,j)
        end do
    end do

    end subroutine d2ndlnfi_dnjdtauddel
    !**************************************************************************


    !**************************************************************************
    subroutine d2ndlnfi_dnjdtau2(gl,Temp, Dens, d2ndlnfidnjdtau2)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE DERIVATIVE OF
    ! THE LOGARITHM OF THE FUGACITY OF COMPONENT i IN THE MIXTURE WITH RESPECT
    ! TO MOLE NUMBER nj MULTIPLIED BY TOTAL NUMBER OF MOLES n WITH RESPECT TO
    ! TAU and TAU AT CONSTANT tau AND x.
    ! THE RESULT IS MULTIPLIED WITH TAU AND TAU!!!
    ! d2/dtau2 (n*d(lnfi)/dnj)_tau,x)
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 3, Eq. (23)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! d2ndlnfidnjdtau2    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens
    Double precision, dimension(30,30):: d2ndlnfidnjdtau2

    integer:: i,j

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative d/dtau (n*d(lnfi)/dnj)_del,x)
    Double precision, dimension(15,30)::dnardnjALL
    integer, dimension(15):: GETDER

    if(.not. allocated(gl%dndnardnidnjALL))allocate(gl%dndnardnidnjALL(15,30,30))
    d2ndlnfidnjdtau2 = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    GETDER = 0
    !Get the tau and tau derivative of ndar_dnj
    GETDER(5) = 1
    call ndar_dni_ALL(gl,Temp, dens, GETDER, dnardnjALL)
    !Get the tau and tau derivative of ndnar_dnidnj
    call ndndar_dnidnj_ALL(gl,Temp, dens, GETDER)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            d2ndlnfidnjdtau2(i,j) = dnardnjALL(5,j) + gl%dndnardnidnjALL(5,i,j)
        end do
    end do

    end subroutine d2ndlnfi_dnjdtau2
    !**************************************************************************


    !**************************************************************************
    subroutine d2ndlnfi_dnjddel2(gl,Temp, Dens, d2ndlnfidnjddel2)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE DERIVATIVE OF
    ! THE LOGARITHM OF THE FUGACITY OF COMPONENT i IN THE MIXTURE WITH RESPECT
    ! TO MOLE NUMBER nj MULTIPLIED BY TOTAL NUMBER OF MOLES n WITH RESPECT TO
    ! DEL and DEL AT CONSTANT tau AND x.
    ! THE RESULT IS MULTIPLIED WITH DEL AND DEL!!!
    ! d2/ddel2 (n*d(lnfi)/dnj)_tau,x)
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 3, Eq. (24)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! d2ndlnfidnjddel2    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens
    Double precision, dimension(30,30):: d2ndlnfidnjddel2

    integer:: i,j

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative d/dtau (n*d(lnfi)/dnj)_del,x)
    Double precision, dimension(15,30)::dnardnjALL
    integer, dimension(15):: GETDER

    if(.not. allocated(gl%dndnardnidnjALL))allocate(gl%dndnardnidnjALL(15,30,30))

    d2ndlnfidnjddel2 = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    GETDER = 0
    !Get the del and tau derivative of ndar_dnj
    GETDER(3) = 1
    call ndar_dni_ALL(gl,Temp, dens, GETDER, dnardnjALL)
    !Get the del derivative of ndnar_dnidnj
    call ndndar_dnidnj_ALL(gl,Temp, dens, GETDER)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            d2ndlnfidnjddel2(i,j) = dnardnjALL(3,j) + gl%dndnardnidnjALL(3,i,j)
        end do
    end do

    end subroutine d2ndlnfi_dnjddel2
    !**************************************************************************


    !**************************************************************************
    subroutine dndlnfi_dnjdxk(gl,Temp, Dens)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE DERIVATIVE OF
    ! THE LOGARITHM OF THE FUGACITY OF COMPONENT i IN THE MIXTURE WITH RESPECT
    ! TO MOLE NUMBER nj MULTIPLIED BY TOTAL NUMBER OF MOLES n WITH RESPECT TO
    ! MOLFRACTIONS xk AT CONSTANT tau, del, and xm.
    ! d/dxk (n*d(lnfi)/dnj)_del,tau,xm)
    !
    !The derivative is given in the supplementary material of the article by
    !Bell and Jäger (XXX), page 3, Eq. (21)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dndlnfidnjdxk    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens

    integer:: i,j,k

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative d/dtau (n*d(lnfi)/dnj)_del,x)
    integer, dimension(15):: GETDER

    if(.not. allocated(gl%dndlnfidnjdxk))allocate(gl%dndlnfidnjdxk(30,30,30))
    if(.not. allocated(gl%ndardnidxi_all))allocate(gl%ndardnidxi_all(15,30,30))
    if(.not. allocated(gl%dndndardnidnjdxk_ALL))allocate(gl%dndndardnidnjdxk_ALL(15,30,30,30))

    gl%dndlnfidnjdxk = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    GETDER = 0
    !Get the xk derivatives of ndar_dnj
    GETDER(1) = 1
    call dndar_dni_dxi_ALL(gl,Temp, dens, GETDER)
    !Get the xk derivatives of ndnar_dnidnj
    call dndndar_dnidnj_dxk_ALL(gl,Temp, dens, GETDER)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                gl%dndlnfidnjdxk(i,j,k) = gl%ndardnidxi_all(1,j,k) + gl%dndndardnidnjdxk_ALL(1,i,j,k)
                if ((i == gl%ncomp) .and. (j == gl%ncomp)) then
                    gl%dndlnfidnjdxk(i,j,k) = gl%dndlnfidnjdxk(i,j,k) + 1.D0 / gl%molfractions(i)**2
                else
                    if ((k == i) .and. (i == j)) then
                        gl%dndlnfidnjdxk(i,j,k) = gl%dndlnfidnjdxk(i,j,k) - 1.D0 / gl%molfractions(i)**2
                    end if
                end if
            end do
        end do
    end do

    end subroutine dndlnfi_dnjdxk
    !**************************************************************************


    !**************************************************************************
    subroutine d2ndlnfi_dnjdxkddel(gl,Temp, Dens)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE DERIVATIVE OF
    ! THE LOGARITHM OF THE FUGACITY OF COMPONENT i IN THE MIXTURE WITH RESPECT
    ! TO MOLE NUMBER nj MULTIPLIED BY TOTAL NUMBER OF MOLES n WITH RESPECT TO
    ! MOLFRACTIONS xk and del AT CONSTANT tau, and xm.
    ! d2 / dxk ddel (n*d(lnfi)/dnj))
    ! THE RESULT IS MULTIPLIED BY delta
    !
    ! The derivative is given in the supplementary material of the article by
    ! Bell and Jäger (XXX), page 3, Eq. (??) - NOT YET THERE
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! d2ndlnfidnjdxkddel    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens

    integer:: i,j,k

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative d/dtau (n*d(lnfi)/dnj)_del,x)
    integer, dimension(15):: GETDER

    if(.not. allocated(gl%d2ndlnfidnjdxkddel))allocate(gl%d2ndlnfidnjdxkddel(15,30,30))
    if(.not. allocated(gl%ndardnidxi_all))allocate(gl%ndardnidxi_all(15,30,30))
    if(.not. allocated(gl%dndndardnidnjdxk_ALL))allocate(gl%dndndardnidnjdxk_ALL(15,30,30,30))

    gl%d2ndlnfidnjdxkddel = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    GETDER = 0
    !Get the xk derivatives of d (ndar_dnj) / d del
    GETDER(2) = 1
    call dndar_dni_dxi_ALL(gl,Temp, dens, GETDER)
    !Get the xk derivatives of d (ndnar_dnidnj) / d del
    call dndndar_dnidnj_dxk_ALL(gl,Temp, dens, GETDER)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                gl%d2ndlnfidnjdxkddel(i,j,k) = gl%ndardnidxi_all(2,j,k) + gl%dndndardnidnjdxk_ALL(2,i,j,k)
            end do
        end do
    end do

    end subroutine d2ndlnfi_dnjdxkddel
    !**************************************************************************


    !**************************************************************************
    subroutine d2ndlnfi_dnjdxkdtau(gl,Temp, Dens)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE DERIVATIVE OF
    ! THE LOGARITHM OF THE FUGACITY OF COMPONENT i IN THE MIXTURE WITH RESPECT
    ! TO MOLE NUMBER nj MULTIPLIED BY TOTAL NUMBER OF MOLES n WITH RESPECT TO
    ! MOLFRACTIONS xk and tau AT CONSTANT del, and xm.
    ! d2 / dxk dtau (n*d(lnfi)/dnj))
    ! THE RESULT IS MULTIPLIED BY tau
    !
    ! The derivative is given in the supplementary material of the article by
    ! Bell and Jäger (XXX), page 3, Eq. (??) - NOT YET THERE
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    !
    ! OUTPUT PARAMETERS:
    ! d2ndlnfidnjdxkdtau    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, Dens

    integer:: i,j,k

    double precision:: delta, tau

    !Other derivatives needed to calculate the derivative d/dtau (n*d(lnfi)/dnj)_del,x)
    integer, dimension(15):: GETDER

    if(.not. allocated(gl%d2ndlnfidnjdxkdtau))allocate(gl%d2ndlnfidnjdxkdtau(15,30,30))
    if(.not. allocated(gl%ndardnidxi_all))allocate(gl%ndardnidxi_all(15,30,30))
    if(.not. allocated(gl%dndndardnidnjdxk_ALL))allocate(gl%dndndardnidnjdxk_ALL(15,30,30,30))

    gl%d2ndlnfidnjdxkdtau = 0.D0

    !Calculate delta and tau
    delta = Dens / gl%rhoredmix
    tau = gl%tredmix / temp

    GETDER = 0
    !Get the xk derivatives of d (ndar_dnj) / d tau
    GETDER(4) = 1
    call dndar_dni_dxi_ALL(gl,Temp, dens, GETDER)
    !Get the xk derivatives of d (ndnar_dnidnj) / d tau
    call dndndar_dnidnj_dxk_ALL(gl,Temp, dens, GETDER)

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                gl%d2ndlnfidnjdxkdtau(i,j,k) = gl%ndardnidxi_all(4,j,k) + gl%dndndardnidnjdxk_ALL(4,i,j,k)
            end do
        end do
    end do

    end subroutine d2ndlnfi_dnjdxkdtau
    !**************************************************************************


    !**************************************************************************
    subroutine ndar_dni_ALL(gl,Temp, dens, GETDER, ndardni_all)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF ALL (UP TO THIRD) DELTA AND TAU
    ! DERIVATIVES OF ONE TERM OF THE RESIDUAL PART OF THE CHEMICAL POTENTIAL:
    ! (n*d(alpha^r)/dni at constant T,V, and n_m)
    !
    ! The calculated derivatives can be written:
    ! d(n+q)/d(del^n)d(tau^q) (n*d(alpha^r)/dni_constT,V,nm)
    ! AT CONSTANT TAU, DEL, AND x
    !
    ! The derivatives are given in the supplementary material of the article by
    ! Bell and Jäger (XXX), page 4+5, Eq. (34-42)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    ! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. RESIDUAL PART OF THE CHEMICAL POTENTIAL (CPR)
    !                2. 1ST DERIVATIVE OF CPR WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
    !                3. 2ND DERIVATIVE OF CPR WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
    !                4. 1ST DERIVATIVE OF CPR WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
    !                5: 2ND DERIVATIVE OF CPR WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
    !                6: 1ST MIXED DERIVATIVE OF CPR WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
    !                7: 3RD MIXED DERIVATIVE OF CPR WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
    !                8: 3RD DERIVATIVE OF CPR WITH RESPECT TO D, MULTIPLIED BY DEL^3
    !                9: 3RD DERIVATIVE OF CPR WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !               10: 3RD MIXED DERIVATIVE OF CPR WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    !               -- NOT YET IMPLEMENTED (NOT NEEDED SO FAR)
    !               11: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^4, MULTIPLIED BY DEL^4
    !               12: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^3 and TAU, MULTIPLIED BY DEL^3*TAU
    !               13: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^2 and TAU^2, MULTIPLIED BY DEL^2*TAU^2
    !               14: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL and TAU^3, MULTIPLIED BY DEL*TAU^3
    !               15: 4TH DERIVATIVE OF CPR WITH RESPECT TO TAU^4, MULTIPLIED BY TAU^4
    !
    ! OUTPUT PARAMETERS:
    ! ndardni_all    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, dens
    Double precision, dimension(15,30)::ndardni_all
    integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed
    integer, dimension(15):: GETDER_help       ! Specifies which help derivatives are calculated

    !Other derivatives needed for the calculation of the derivatives of the residual chemical potential wrt delta and tau
    double precision, dimension(15):: MIXDERIVFNR
    double precision, dimension(15,30):: MIXDERIVFNR_dxi

    double precision:: delta, del2, del3, del4, tau, tau2, tau3, tau4
    integer:: i,k

    if(.not. allocated(gl%d2PSI_Y_dxjdxk)) allocate(gl%d2PSI_Y_dxjdxk(2,30,30,30))

    call Help_PSI_Derivs(gl)

    !FOR CUBICS ONLY AT THE MOMENT!!!
    GETDER_help = 1 !Calculate all help derivatives
    call MIXDERIVSFNR_HIGHER_CUBIC (gl,Temp, Dens, GETDER_help, MIXDERIVFNR)
    call MIXDERIVSFNR_dxi_CUBIC (gl,Temp, Dens, GETDER_help, MIXDERIVFNR_dxi)


    delta = dens / gl%rhoredmix
    del2 = delta*delta
    del3 = del2*delta
    del4 = del2*del2

    tau = gl%tredmix / temp
    tau2 = tau * tau
    tau3 = tau2 * tau
    tau4 = tau2 * tau2


    ndardni_all = 0.D0

    !alphar
    if (GETDER(1) == 1) then
        ndardni_all(1,:) = MIXDERIVFNR(2) * gl%PSI_Y(2,:) + MIXDERIVFNR(4) * gl%PSI_Y(1,:) + MIXDERIVFNR_dxi(1,:)
        do i = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                ndardni_all(1,i) = ndardni_all(1,i) - gl%molfractions(k) * MIXDERIVFNR_dxi(1,k)
            end do
        end do
    end if

    !d_alphar_d_delta * delta
    if (GETDER(2) == 1) then
        ndardni_all(2,:) = (MIXDERIVFNR(2) + MIXDERIVFNR(3))  * gl%PSI_Y(2,:) + MIXDERIVFNR(6) * gl%PSI_Y(1,:) + MIXDERIVFNR_dxi(2,:)
        do i = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                ndardni_all(2,i) = ndardni_all(2,i) - gl%molfractions(k) * MIXDERIVFNR_dxi(2,k)
            end do
        end do
    end if

    !d2_alphar_d_delta2 * delta^2
    if (GETDER(3) == 1) then
        ndardni_all(3,:) = (2.D0 * MIXDERIVFNR(3) + MIXDERIVFNR(8))  * gl%PSI_Y(2,:) + MIXDERIVFNR(10) * gl%PSI_Y(1,:) + MIXDERIVFNR_dxi(3,:)
        do i = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                ndardni_all(3,i) = ndardni_all(3,i) - gl%molfractions(k) * MIXDERIVFNR_dxi(3,k)
            end do
        end do
    end if

    !d_alphar_d_tau * tau
    if (GETDER(4) == 1) then
        ndardni_all(4,:) = MIXDERIVFNR(6) * gl%PSI_Y(2,:) + (MIXDERIVFNR(5) + MIXDERIVFNR(4)) * gl%PSI_Y(1,:) + MIXDERIVFNR_dxi(4,:)
        do i = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                ndardni_all(4,i) = ndardni_all(4,i) - gl%molfractions(k) * MIXDERIVFNR_dxi(4,k)
            end do
        end do
    end if

    !d2_alphar_d_tau2 * tau^2
    if (GETDER(5) == 1) then
        ndardni_all(5,:) = MIXDERIVFNR(7) * gl%PSI_Y(2,:) + (2.D0 * MIXDERIVFNR(5) + MIXDERIVFNR(9)) * gl%PSI_Y(1,:) + MIXDERIVFNR_dxi(5,:)
        do i = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                ndardni_all(5,i) = ndardni_all(5,i) - gl%molfractions(k) * MIXDERIVFNR_dxi(5,k)
            end do
        end do
    end if

    !d2_alphar_d_delta_d_tau * delta * tau
    if (GETDER(6) == 1) then
        ndardni_all(6,:) = (MIXDERIVFNR(6) + MIXDERIVFNR(10)) * gl%PSI_Y(2,:) + (MIXDERIVFNR(7) + MIXDERIVFNR(6)) * gl%PSI_Y(1,:) + MIXDERIVFNR_dxi(6,:)
        do i = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                ndardni_all(6,i) = ndardni_all(6,i) - gl%molfractions(k) * MIXDERIVFNR_dxi(6,k)
            end do
        end do
    end if

    !d3_alphar_d_delta_d_tau2 * delta * tau^2
    if (GETDER(7) == 1) then
        ndardni_all(7,:) = (MIXDERIVFNR(7) + MIXDERIVFNR(13)) * gl%PSI_Y(2,:) + (MIXDERIVFNR(14) + MIXDERIVFNR(7)) * gl%PSI_Y(1,:) + MIXDERIVFNR_dxi(7,:)
        do i = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                ndardni_all(7,i) = ndardni_all(7,i) - gl%molfractions(k) * MIXDERIVFNR_dxi(7,k)
            end do
        end do
    end if

    !d3_alphar_d_delta3 * delta^3
    if (GETDER(8) == 1) then
        ndardni_all(8,:) = (3.D0 * MIXDERIVFNR(8) + MIXDERIVFNR(11))  * gl%PSI_Y(2,:) + MIXDERIVFNR(12) * gl%PSI_Y(1,:) + MIXDERIVFNR_dxi(8,:)
        do i = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                ndardni_all(8,i) = ndardni_all(8,i) - gl%molfractions(k) * MIXDERIVFNR_dxi(8,k)
            end do
        end do
    end if

    !d3_alphar_d_tau3 * tau^3
    if (GETDER(9) == 1) then
        ndardni_all(9,:) = MIXDERIVFNR(14) * gl%PSI_Y(2,:) + (3.D0 * MIXDERIVFNR(9) + MIXDERIVFNR(15)) * gl%PSI_Y(1,:) + MIXDERIVFNR_dxi(9,:)
        do i = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                ndardni_all(9,i) = ndardni_all(9,i) - gl%molfractions(k) * MIXDERIVFNR_dxi(9,k)
            end do
        end do
    end if

    !d3_alphar_d_delta2_dtau * delta^2 * tau
    if (GETDER(10) == 1) then
        ndardni_all(10,:) = (2.D0 * MIXDERIVFNR(10) + MIXDERIVFNR(12)) * gl%PSI_Y(2,:) + (MIXDERIVFNR(13) + MIXDERIVFNR(10)) * gl%PSI_Y(1,:) + MIXDERIVFNR_dxi(10,:)
        do i = 1, gl%ncomp
            do k = 1, gl%ncomp-1
                ndardni_all(10,i) = ndardni_all(10,i) - gl%molfractions(k) * MIXDERIVFNR_dxi(10,k)
            end do
        end do
    end if



    !Following derivatives not yet implemented:

    !d4_alphar_d_delta4 * delta^4
    if (GETDER(11) == 1) then
        ndardni_all(11,:) = 0.D0
    end if

    !d4_alphar_d_delta3_dtau * delta^3*tau
    if (GETDER(12) == 1) then
        ndardni_all(12,:) = 0.D0
    end if

    !d4_alphar_d_delta2_dtau2 * delta^2*tau^2
    if (GETDER(13) == 1) then
        ndardni_all(13,:) = 0.D0
    end if

    !d4_alphar_d_delta_dtau3 * delta*tau^3
    if (GETDER(14) == 1) then
        ndardni_all(14,:) = 0.D0
    end if

    !d4_alphar_d_tau4 * tau^4
    if (GETDER(15) == 1) then
        ndardni_all(15,:) = 0.D0
    end if


    end subroutine ndar_dni_ALL
    !**************************************************************************


    !**************************************************************************
    subroutine dndar_dni_dxi_ALL(gl,Temp, dens, GETDER)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF ALL (UP TO SECOND) DELTA AND TAU
    ! DERIVATIVES OF ONE TERM OF THE RESIDUAL PART OF THE CHEMICAL POTENTIAL:
    ! (n*d(alpha^r)/dni at constant T,V, and n_m)
    ! ALL DERIVATIVES ALSO WITH RESPECT TO MOLFRACTIONS xi
    !
    ! The calculated derivatives can be written:
    ! d d(n+q)/d(del^n)d(tau^q) (n*d(alpha^r)/dni_constT,V,nm) / dxi
    ! AT CONSTANT TAU, DEL, AND x
    !
    ! The derivatives are given in the supplementary material of the article by
    ! Bell and Jäger (XXX), page 5+6, Eq. (43-48)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    ! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
    !               ALL DERIVATIVES ALSO WITH RESPECT TO xj
    !                1. RESIDUAL PART OF THE CHEMICAL POTENTIAL (CPR)
    !                2. 1ST DERIVATIVE OF CPR WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
    !                3. 2ND DERIVATIVE OF CPR WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
    !                4. 1ST DERIVATIVE OF CPR WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
    !                5: 2ND DERIVATIVE OF CPR WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
    !                6: 1ST MIXED DERIVATIVE OF CPR WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
    !               -- NOT YET IMPLEMENTED (NOT NEEDED SO FAR)
    !                7: 3RD MIXED DERIVATIVE OF CPR WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
    !                8: 3RD DERIVATIVE OF CPR WITH RESPECT TO D, MULTIPLIED BY DEL^3
    !                9: 3RD DERIVATIVE OF CPR WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !               10: 3RD MIXED DERIVATIVE OF CPR WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    !               11: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^4, MULTIPLIED BY DEL^4
    !               12: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^3 and TAU, MULTIPLIED BY DEL^3*TAU
    !               13: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^2 and TAU^2, MULTIPLIED BY DEL^2*TAU^2
    !               14: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL and TAU^3, MULTIPLIED BY DEL*TAU^3
    !               15: 4TH DERIVATIVE OF CPR WITH RESPECT TO TAU^4, MULTIPLIED BY TAU^4
    !
    ! OUTPUT PARAMETERS:
    ! ndardnidxi_all    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, dens
    integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed
    integer, dimension(15):: GETDER_help       ! Specifies which help derivatives are calculated

    !Other derivatives needed for the calculation of the derivatives of the residual chemical potential wrt delta and tau
    double precision, dimension(15):: MIXDERIVFNR
    double precision, dimension(15,30):: MIXDERIVFNR_dxi

    double precision:: delta, del2, del3, del4, tau, tau2, tau3, tau4
    integer:: i,j,k

    if(.not. allocated(gl%ndardnidxi_all)) allocate(gl%ndardnidxi_all(15,30,30))
    if(.not. allocated(gl%MIXDERIVFNR_dxidxj)) allocate(gl%MIXDERIVFNR_dxidxj(15,30,30))
    if(.not. allocated(gl%dPSI_Y_dxj)) allocate(gl%dPSI_Y_dxj(2,30,30))
    if(.not. allocated(gl%d2PSI_Y_dxjdxk)) allocate(gl%d2PSI_Y_dxjdxk(2,30,30,30))


    call Help_PSI_Derivs(gl)

    !FOR CUBICS ONLY AT THE MOMENT!!!
    GETDER_help = 1 !Calculate all help derivatives
    call MIXDERIVSFNR_HIGHER_CUBIC (gl,Temp, Dens, GETDER_help, MIXDERIVFNR)
    call MIXDERIVSFNR_dxi_CUBIC (gl,Temp, Dens, GETDER_help, MIXDERIVFNR_dxi)
    call MIXDERIVSFNR_dxidxj_CUBIC (gl,Temp, Dens, GETDER_help)

    delta = dens / gl%rhoredmix
    del2 = delta*delta
    del3 = del2*delta
    del4 = del2*del2

    tau = gl%tredmix / temp
    tau2 = tau * tau
    tau3 = tau2 * tau
    tau4 = tau2 * tau2


    gl%ndardnidxi_all = 0.D0

    !alphar
    if (GETDER(1) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp-1
                gl%ndardnidxi_all(1,i,j) = MIXDERIVFNR_dxi(2,j) * gl%PSI_Y(2,i) + MIXDERIVFNR(2) * gl%dPSI_Y_dxj(2,i,j) &
                    & + MIXDERIVFNR_dxi(4,j) * gl%PSI_Y(1,i) + MIXDERIVFNR(4) * gl%dPSI_Y_dxj(1,i,j) &
                    & + gl%MIXDERIVFNR_dxidxj(1,i,j) - MIXDERIVFNR_dxi(1,j)
                do k = 1, gl%ncomp-1
                    gl%ndardnidxi_all(1,i,j) = gl%ndardnidxi_all(1,i,j) - gl%molfractions(k) * gl%MIXDERIVFNR_dxidxj(1,j,k)
                end do
            end do
        end do
    end if

    !d_alphar_d_delta * delta
    if (GETDER(2) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp-1
                gl%ndardnidxi_all(2,i,j) = (MIXDERIVFNR_dxi(2,j) + MIXDERIVFNR_dxi(3,j))  * gl%PSI_Y(2,i) + &
                    &    (MIXDERIVFNR(2) + MIXDERIVFNR(3)) * gl%dPSI_Y_dxj(2,i,j) &
                    & + MIXDERIVFNR_dxi(6,j) * gl%PSI_Y(1,i) + MIXDERIVFNR(6) * gl%dPSI_Y_dxj(1,i,j) &
                    & + gl%MIXDERIVFNR_dxidxj(2,i,j) - MIXDERIVFNR_dxi(2,j)
                do k = 1, gl%ncomp-1
                    gl%ndardnidxi_all(2,i,j) = gl%ndardnidxi_all(2,i,j) - gl%molfractions(k) * gl%MIXDERIVFNR_dxidxj(2,j,k)
                end do
            end do
        end do
    end if

    !d2_alphar_d_delta2 * delta^2
    if (GETDER(3) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp-1
                gl%ndardnidxi_all(3,i,j) = (2.D0 * MIXDERIVFNR_dxi(3,j) + MIXDERIVFNR_dxi(8,j))  * gl%PSI_Y(2,i) + &
                    &    (2.D0 * MIXDERIVFNR(3) + MIXDERIVFNR(8)) * gl%dPSI_Y_dxj(2,i,j) &
                    & + MIXDERIVFNR_dxi(10,j) * gl%PSI_Y(1,i) + MIXDERIVFNR(10) * gl%dPSI_Y_dxj(1,i,j) &
                    & + gl%MIXDERIVFNR_dxidxj(3,i,j) - MIXDERIVFNR_dxi(3,j)
                do k = 1, gl%ncomp-1
                    gl%ndardnidxi_all(3,i,j) = gl%ndardnidxi_all(3,i,j) - gl%molfractions(k) * gl%MIXDERIVFNR_dxidxj(3,j,k)
                end do
            end do
        end do
    end if

    !d_alphar_d_tau * tau
    if (GETDER(4) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp-1
                gl%ndardnidxi_all(4,i,j) = MIXDERIVFNR_dxi(6,j)  * gl%PSI_Y(2,i) + MIXDERIVFNR(6) * gl%dPSI_Y_dxj(2,i,j) &
                    & + (MIXDERIVFNR_dxi(5,j) + MIXDERIVFNR_dxi(4,j)) * gl%PSI_Y(1,i) &
                    & + (MIXDERIVFNR(5) + MIXDERIVFNR(4))  * gl%dPSI_Y_dxj(1,i,j) &
                    & + gl%MIXDERIVFNR_dxidxj(4,i,j) - MIXDERIVFNR_dxi(4,j)
                do k = 1, gl%ncomp-1
                    gl%ndardnidxi_all(4,i,j) = gl%ndardnidxi_all(4,i,j) - gl%molfractions(k) * gl%MIXDERIVFNR_dxidxj(4,j,k)
                end do
            end do
        end do
    end if

    !d2_alphar_d_tau2 * tau^2
    if (GETDER(5) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp-1
                gl%ndardnidxi_all(5,i,j) = MIXDERIVFNR_dxi(7,j)  * gl%PSI_Y(2,i) + MIXDERIVFNR(7) * gl%dPSI_Y_dxj(2,i,j) &
                    & + (MIXDERIVFNR_dxi(9,j) + 2.D0 * MIXDERIVFNR_dxi(5,j)) * gl%PSI_Y(1,i) &
                    & + (MIXDERIVFNR(9) + 2.D0 * MIXDERIVFNR(5))  * gl%dPSI_Y_dxj(1,i,j) &
                    & + gl%MIXDERIVFNR_dxidxj(5,i,j) - MIXDERIVFNR_dxi(5,j)
                do k = 1, gl%ncomp-1
                    gl%ndardnidxi_all(5,i,j) = gl%ndardnidxi_all(5,i,j) - gl%molfractions(k) * gl%MIXDERIVFNR_dxidxj(5,j,k)
                end do
            end do
        end do
    end if

    !d2_alphar_d_delta_d_tau * delta * tau
    if (GETDER(6) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp-1
                gl%ndardnidxi_all(6,i,j) = (MIXDERIVFNR_dxi(6,j) + MIXDERIVFNR_dxi(10,j))  * gl%PSI_Y(2,i) + &
                    &    (MIXDERIVFNR(6) + MIXDERIVFNR(10)) * gl%dPSI_Y_dxj(2,i,j) &
                    & + (MIXDERIVFNR_dxi(7,j) + MIXDERIVFNR_dxi(6,j)) * gl%PSI_Y(1,i) &
                    & + (MIXDERIVFNR(7) + MIXDERIVFNR(6))  * gl%dPSI_Y_dxj(1,i,j) &
                    & + gl%MIXDERIVFNR_dxidxj(6,i,j) - MIXDERIVFNR_dxi(6,j)
                do k = 1, gl%ncomp-1
                    gl%ndardnidxi_all(6,i,j) = gl%ndardnidxi_all(6,i,j) - gl%molfractions(k) * gl%MIXDERIVFNR_dxidxj(6,j,k)
                end do
            end do
        end do
    end if




    !Following derivatives not yet implemented:

    !d3_alphar_d_delta_d_tau2 * delta * tau^2
    if (GETDER(7) == 1) then
        gl%ndardnidxi_all(7,:,:) = 0.D0
    end if

    !d3_alphar_d_delta3 * delta^3
    if (GETDER(8) == 1) then
        gl%ndardnidxi_all(8,:,:) = 0.D0
    end if

    !d3_alphar_d_tau3 * tau^3
    if (GETDER(9) == 1) then
        gl%ndardnidxi_all(9,:,:) = 0.D0
    end if

    !d3_alphar_d_delta2_dtau * delta^2 * tau
    if (GETDER(10) == 1) then
        gl%ndardnidxi_all(10,:,:) = 0.D0
    end if

    !d4_alphar_d_delta4 * delta^4
    if (GETDER(11) == 1) then
        gl%ndardnidxi_all(11,:,:) = 0.D0
    end if

    !d4_alphar_d_delta3_dtau * delta^3*tau
    if (GETDER(12) == 1) then
        gl%ndardnidxi_all(12,:,:) = 0.D0
    end if

    !d4_alphar_d_delta2_dtau2 * delta^2*tau^2
    if (GETDER(13) == 1) then
        gl%ndardnidxi_all(13,:,:) = 0.D0
    end if

    !d4_alphar_d_delta_dtau3 * delta*tau^3
    if (GETDER(14) == 1) then
        gl%ndardnidxi_all(14,:,:) = 0.D0
    end if

    !d4_alphar_d_tau4 * tau^4
    if (GETDER(15) == 1) then
        gl%ndardnidxi_all(15,:,:) = 0.D0
    end if


    end subroutine dndar_dni_dxi_ALL
    !**************************************************************************


    !**************************************************************************
    subroutine d2ndar_dni_dxidxj_ALL(gl,Temp, dens, GETDER)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF ALL (UP TO FIRST) DELTA AND TAU
    ! DERIVATIVES OF ONE TERM OF THE RESIDUAL PART OF THE CHEMICAL POTENTIAL:
    ! (n*d(alpha^r)/dni at constant T,V, and n_m)
    ! ALL DERIVATIVES ALSO WITH RESPECT TO MOLFRACTIONS xi AND xj
    !
    ! The calculated derivatives can be written:
    ! d2 d(n+q)/d(del^n)d(tau^q) (n*d(alpha^r)/dni_constT,V,nm) / dxidxj
    ! AT CONSTANT TAU, DEL, AND x
    !
    ! The derivatives are given in the supplementary material of the article by
    ! Bell and Jäger (XXX), page 6+7, Eq. (49-51)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    ! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
    !               ALL DERIVATIVES ALSO WITH RESPECT TO xj
    !                1. RESIDUAL PART OF THE CHEMICAL POTENTIAL (CPR)
    !                2. 1ST DERIVATIVE OF CPR WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
    !                3. 2ND DERIVATIVE OF CPR WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2  - NOT YET IMPLEMENTED
    !                4. 1ST DERIVATIVE OF CPR WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
    !               -- NOT YET IMPLEMENTED (NOT NEEDED SO FAR)
    !                5: 2ND DERIVATIVE OF CPR WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
    !                6: 1ST MIXED DERIVATIVE OF CPR WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
    !                7: 3RD MIXED DERIVATIVE OF CPR WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
    !                8: 3RD DERIVATIVE OF CPR WITH RESPECT TO D, MULTIPLIED BY DEL^3
    !                9: 3RD DERIVATIVE OF CPR WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !               10: 3RD MIXED DERIVATIVE OF CPR WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    !               11: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^4, MULTIPLIED BY DEL^4
    !               12: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^3 and TAU, MULTIPLIED BY DEL^3*TAU
    !               13: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^2 and TAU^2, MULTIPLIED BY DEL^2*TAU^2
    !               14: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL and TAU^3, MULTIPLIED BY DEL*TAU^3
    !               15: 4TH DERIVATIVE OF CPR WITH RESPECT TO TAU^4, MULTIPLIED BY TAU^4
    !
    ! OUTPUT PARAMETERS:
    ! d2ndardnidxidxj_all    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, dens
    integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed
    integer, dimension(15):: GETDER_help       ! Specifies which help derivatives are calculated

    !Other derivatives needed for the calculation of the derivatives of the residual chemical potential wrt delta and tau
    double precision, dimension(15):: MIXDERIVFNR
    double precision, dimension(15,30):: MIXDERIVFNR_dxi


    double precision:: delta, del2, del3, del4, tau, tau2, tau3, tau4
    integer:: i,j,k,m

    if(.not. allocated(gl%MIXDERIVFNR_dxidxj)) allocate(gl%MIXDERIVFNR_dxidxj(15,30,30))
    if(.not. allocated(gl%MIXDERIVFNR_dxidxjdxk)) allocate(gl%MIXDERIVFNR_dxidxjdxk(15,30,30,30))
    if(.not. allocated(gl%PSI_Y)) allocate(gl%PSI_Y(2,30))
    if(.not. allocated(gl%dPSI_Y_dxj)) allocate(gl%dPSI_Y_dxj(2,30,30))
    if(.not. allocated(gl%d2PSI_Y_dxjdxk)) allocate(gl%d2PSI_Y_dxjdxk(2,30,30,30))

    call Help_PSI_Derivs(gl)

    !FOR CUBICS ONLY AT THE MOMENT!!!
    GETDER_help = 1 !Calculate all help derivatives
    call MIXDERIVSFNR_HIGHER_CUBIC (gl,Temp, Dens, GETDER_help, MIXDERIVFNR)
    call MIXDERIVSFNR_dxi_CUBIC (gl,Temp, Dens, GETDER_help, MIXDERIVFNR_dxi)
    call MIXDERIVSFNR_dxidxj_CUBIC (gl,Temp, Dens, GETDER_help)
    call MIXDERIVSFNR_dxidxjdxk_CUBIC (gl,Temp, Dens, GETDER_help)


    delta = dens / gl%rhoredmix
    del2 = delta*delta
    del3 = del2*delta
    del4 = del2*del2

    tau = gl%tredmix / temp
    tau2 = tau * tau
    tau3 = tau2 * tau
    tau4 = tau2 * tau2


    gl%d2ndardnidxidxj_all = 0.D0

    !alphar
    if (GETDER(1) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp-1
                do k = 1, gl%ncomp-1
                    gl%d2ndardnidxidxj_all(1,i,j,k) = MIXDERIVFNR_dxi(2,j) * gl%dPSI_Y_dxj(2,i,k) + gl%MIXDERIVFNR_dxidxj(2,j,k) * gl%PSI_Y(2,i) &
                        & + MIXDERIVFNR(2) * gl%d2PSI_Y_dxjdxk(2,i,j,k) + MIXDERIVFNR_dxi(2,k) * gl%dPSI_Y_dxj(2,i,j) &
                        & + MIXDERIVFNR_dxi(4,j) * gl%dPSI_Y_dxj(1,i,k) + gl%MIXDERIVFNR_dxidxj(4,j,k) * gl%PSI_Y(1,i) &
                        & + MIXDERIVFNR(4) * gl%d2PSI_Y_dxjdxk(1,i,j,k) + MIXDERIVFNR_dxi(4,k) * gl%dPSI_Y_dxj(1,i,j) &
                        & + gl%MIXDERIVFNR_dxidxjdxk(1,i,j,k) - 2.D0 * gl%MIXDERIVFNR_dxidxj(1,j,k)
                    do m = 1, gl%ncomp-1
                        gl%d2ndardnidxidxj_all(1,i,j,k) = gl%d2ndardnidxidxj_all(1,i,j,k) &
                            & - gl%molfractions(m) * gl%MIXDERIVFNR_dxidxjdxk(1,j,m,k)
                    end do
                end do
            end do
        end do
    end if

    !d_alphar_d_delta * delta
    if (GETDER(2) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp-1
                do k = 1, gl%ncomp-1
                    gl%d2ndardnidxidxj_all(2,i,j,k) = (MIXDERIVFNR_dxi(3,j) + MIXDERIVFNR_dxi(2,j)) * gl%dPSI_Y_dxj(2,i,k) &
                        & + (gl%MIXDERIVFNR_dxidxj(3,j,k) + gl%MIXDERIVFNR_dxidxj(2,j,k)) * gl%PSI_Y(2,i) &
                        & + (MIXDERIVFNR(3) + MIXDERIVFNR(2)) * gl%d2PSI_Y_dxjdxk(2,i,j,k) &
                        & + (MIXDERIVFNR_dxi(3,k) + MIXDERIVFNR_dxi(2,k)) * gl%dPSI_Y_dxj(2,i,j) &
                        & + MIXDERIVFNR_dxi(6,j) * gl%dPSI_Y_dxj(1,i,k) + gl%MIXDERIVFNR_dxidxj(6,j,k) * gl%PSI_Y(1,i) &
                        & + MIXDERIVFNR(6) * gl%d2PSI_Y_dxjdxk(1,i,j,k) + MIXDERIVFNR_dxi(6,k) * gl%dPSI_Y_dxj(1,i,j) &
                        & + gl%MIXDERIVFNR_dxidxjdxk(2,i,j,k) - 2.D0 * gl%MIXDERIVFNR_dxidxj(2,j,k)
                    do m = 1, gl%ncomp-1
                        gl%d2ndardnidxidxj_all(2,i,j,k) = gl%d2ndardnidxidxj_all(2,i,j,k) &
                            & - gl%molfractions(m) * gl%MIXDERIVFNR_dxidxjdxk(2,j,m,k)
                    end do
                end do
            end do
        end do
    end if

    !d2_alphar_d_delta2 * delta^2
    if (GETDER(3) == 1) then
        !Derivative not yet implemented
        gl%d2ndardnidxidxj_all(3,:,:,:) = 0.D0
    end if

    !d_alphar_d_tau * tau
    if (GETDER(4) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp-1
                do k = 1, gl%ncomp-1
                    gl%d2ndardnidxidxj_all(4,i,j,k) = MIXDERIVFNR_dxi(6,j) * gl%dPSI_Y_dxj(2,i,k) + gl%MIXDERIVFNR_dxidxj(6,j,k) * gl%PSI_Y(2,i) &
                        & + MIXDERIVFNR(6) * gl%d2PSI_Y_dxjdxk(2,i,j,k) + MIXDERIVFNR_dxi(6,k) * gl%dPSI_Y_dxj(2,i,j) &
                        & + (MIXDERIVFNR_dxi(5,j) + MIXDERIVFNR_dxi(4,j)) * gl%dPSI_Y_dxj(1,i,k) &
                        & + (gl%MIXDERIVFNR_dxidxj(5,j,k) + gl%MIXDERIVFNR_dxidxj(4,j,k)) * gl%PSI_Y(1,i) &
                        & + (MIXDERIVFNR(5) + MIXDERIVFNR(4)) * gl%d2PSI_Y_dxjdxk(1,i,j,k) &
                        & + (MIXDERIVFNR_dxi(5,k) + MIXDERIVFNR_dxi(4,k)) * gl%dPSI_Y_dxj(1,i,j) &
                        & + gl%MIXDERIVFNR_dxidxjdxk(4,i,j,k) - 2.D0 * gl%MIXDERIVFNR_dxidxj(4,j,k)
                    do m = 1, gl%ncomp-1
                        gl%d2ndardnidxidxj_all(4,i,j,k) = gl%d2ndardnidxidxj_all(4,i,j,k) &
                            & - gl%molfractions(m) * gl%MIXDERIVFNR_dxidxjdxk(4,j,m,k)
                    end do
                end do
            end do
        end do
    end if



    !Following derivatives not yet implemented:

    !d2_alphar_d_tau2 * tau^2
    if (GETDER(5) == 1) then
        gl%d2ndardnidxidxj_all(5,:,:,:) = 0.D0
    end if

    !d2_alphar_d_delta_d_tau * delta * tau
    if (GETDER(6) == 1) then
        gl%d2ndardnidxidxj_all(6,:,:,:) = 0.D0
    end if

    !d3_alphar_d_delta_d_tau2 * delta * tau^2
    if (GETDER(7) == 1) then
        gl%d2ndardnidxidxj_all(7,:,:,:) = 0.D0
    end if

    !d3_alphar_d_delta3 * delta^3
    if (GETDER(8) == 1) then
        gl%d2ndardnidxidxj_all(8,:,:,:) = 0.D0
    end if

    !d3_alphar_d_tau3 * tau^3
    if (GETDER(9) == 1) then
        gl%d2ndardnidxidxj_all(9,:,:,:) = 0.D0
    end if

    !d3_alphar_d_delta2_dtau * delta^2 * tau
    if (GETDER(10) == 1) then
        gl%d2ndardnidxidxj_all(10,:,:,:) = 0.D0
    end if

    !d4_alphar_d_delta4 * delta^4
    if (GETDER(11) == 1) then
        gl%d2ndardnidxidxj_all(11,:,:,:) = 0.D0
    end if

    !d4_alphar_d_delta3_dtau * delta^3*tau
    if (GETDER(12) == 1) then
        gl%d2ndardnidxidxj_all(12,:,:,:) = 0.D0
    end if

    !d4_alphar_d_delta2_dtau2 * delta^2*tau^2
    if (GETDER(13) == 1) then
        gl%d2ndardnidxidxj_all(13,:,:,:) = 0.D0
    end if

    !d4_alphar_d_delta_dtau3 * delta*tau^3
    if (GETDER(14) == 1) then
        gl%d2ndardnidxidxj_all(14,:,:,:) = 0.D0
    end if

    !d4_alphar_d_tau4 * tau^4
    if (GETDER(15) == 1) then
        gl%d2ndardnidxidxj_all(15,:,:,:) = 0.D0
    end if


    end subroutine d2ndar_dni_dxidxj_ALL
    !**************************************************************************


    !**************************************************************************
    subroutine ndndar_dnidnj_ALL(gl,Temp, dens, GETDER)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF ALL (UP TO SECOND) DELTA AND TAU
    ! DERIVATIVES OF Eq. 7.47 of the GERG-MONOGRAPH (Kunz et al. (2007))
    ! (n*d(n*dalpha^r/dni)/dnj at constant T,V, and n_m)
    !
    ! The calculated derivatives can be written:
    ! d(n+q)/d(del^n)d(tau^q) (n*d(n*dalpha^r/dni)/dnj)
    ! AT CONSTANT TAU, DEL, AND x
    !
    ! The derivatives are given in the supplementary material of the article by
    ! Bell and Jäger (XXX), page 7-9, Eq. (52-57)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    ! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. RESIDUAL PART OF THE CHEMICAL POTENTIAL (CPR)
    !                2. 1ST DERIVATIVE OF CPR WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
    !                3. 2ND DERIVATIVE OF CPR WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
    !                4. 1ST DERIVATIVE OF CPR WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
    !                5: 2ND DERIVATIVE OF CPR WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
    !                6: 1ST MIXED DERIVATIVE OF CPR WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
    !               -- NOT YET IMPLEMENTED (NOT NEEDED SO FAR)
    !                7: 3RD MIXED DERIVATIVE OF CPR WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
    !                8: 3RD DERIVATIVE OF CPR WITH RESPECT TO D, MULTIPLIED BY DEL^3
    !                9: 3RD DERIVATIVE OF CPR WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !               10: 3RD MIXED DERIVATIVE OF CPR WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    !               11: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^4, MULTIPLIED BY DEL^4
    !               12: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^3 and TAU, MULTIPLIED BY DEL^3*TAU
    !               13: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^2 and TAU^2, MULTIPLIED BY DEL^2*TAU^2
    !               14: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL and TAU^3, MULTIPLIED BY DEL*TAU^3
    !               15: 4TH DERIVATIVE OF CPR WITH RESPECT TO TAU^4, MULTIPLIED BY TAU^4
    !
    ! OUTPUT PARAMETERS:
    ! ndndardnidnj_ALL    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, dens
    integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed
    integer, dimension(15):: GETDER_help       ! Specifies which help derivatives are calculated

    !Other derivatives needed for the calculation of the derivatives of Eq. 7.47 of GERG
    Double precision, dimension(15,30):: ndardni_ALL
    Double precision, dimension(30):: nddel_dnj, ndtau_dnj
    Double precision, dimension(30):: dnddel_dnjddel, dndtau_dnjdtau
    double precision, dimension(30):: ndTred_dni, ndrhored_dni


    double precision:: delta, del2, del3, del4, tau, tau2, tau3, tau4
    integer:: i,j,m

    if(.not. allocated(gl%ndndardnidnj_ALL)) allocate(gl%ndndardnidnj_ALL(15,30,30))
    if(.not. allocated(gl%dndardnidxi_ALL)) allocate(gl%dndardnidxi_ALL(15,30,30))


    !FOR CUBICS ONLY AT THE MOMENT!!!
    GETDER_help = 1 !Calculate all help derivatives
    call ndar_dni_ALL(gl,Temp, Dens, GETDER_help, ndardni_ALL)
    call dndar_dni_dxi_ALL(gl,Temp, Dens, GETDER_help)
    call ndYr_dni(gl,ndTred_dni, ndrhored_dni)


    delta = dens / gl%rhoredmix
    del2 = delta*delta
    del3 = del2*delta
    del4 = del2*del2

    tau = gl%tredmix / temp
    tau2 = tau * tau
    tau3 = tau2 * tau
    tau4 = tau2 * tau2



    ndtau_dnj = (tau / gl%tredmix * ndTred_dni) / tau                      !Divided by tau, because all tau derivatives are already multiplied by tau
    nddel_dnj = (delta - delta / gl%rhoredmix * ndrhored_dni) / delta      !Divided by delta, because all delta derivatives are already multiplied by delta

    dndtau_dnjdtau = 1.D0 / gl%tredmix * ndTred_dni
    dnddel_dnjddel = 1.D0 - 1.D0 / gl%rhoredmix * ndrhored_dni

    gl%ndndardnidnj_ALL = 0.D0

    !alphar
    if (GETDER(1) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                gl%ndndardnidnj_ALL(1,i,j) = ndardni_ALL(2,i) * nddel_dnj(j) + ndardni_ALL(4,i) * ndtau_dnj(j) &
                    & + gl%dndardnidxi_ALL(1,i,j)
                do m = 1,gl%ncomp-1
                    gl%ndndardnidnj_ALL(1,i,j) = gl%ndndardnidnj_ALL(1,i,j) - gl%molfractions(m) * gl%dndardnidxi_ALL(1,i,m)
                end do
            end do
        end do
    end if

    !d_alphar_d_delta * delta
    if (GETDER(2) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                gl%ndndardnidnj_ALL(2,i,j) = ndardni_ALL(3,i) * nddel_dnj(j) + ndardni_ALL(6,i) * ndtau_dnj(j) &
                    & + ndardni_ALL(2,i) * dnddel_dnjddel(j) &
                    & + gl%dndardnidxi_ALL(2,i,j)
                do m = 1,gl%ncomp-1
                    gl%ndndardnidnj_ALL(2,i,j) = gl%ndndardnidnj_ALL(2,i,j) - gl%molfractions(m) * gl%dndardnidxi_ALL(2,i,m)
                end do
            end do
        end do
    end if

    !d2_alphar_d_delta2 * delta^2
    if (GETDER(3) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                gl%ndndardnidnj_ALL(3,i,j) = ndardni_ALL(8,i) * nddel_dnj(j) + 2.D0 * ndardni_ALL(3,i) * dnddel_dnjddel(j) &
                    & + ndardni_ALL(10,i) * ndtau_dnj(j) &
                    & + gl%dndardnidxi_ALL(3,i,j)
                do m = 1,gl%ncomp-1
                    gl%ndndardnidnj_ALL(3,i,j) = gl%ndndardnidnj_ALL(3,i,j) - gl%molfractions(m) * gl%dndardnidxi_ALL(3,i,m)
                end do
            end do
        end do
    end if

    !d_alphar_d_tau * tau
    if (GETDER(4) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                gl%ndndardnidnj_ALL(4,i,j) = ndardni_ALL(6,i) * nddel_dnj(j) + ndardni_ALL(5,i) * ndtau_dnj(j) &
                    & + ndardni_ALL(4,i) * dndtau_dnjdtau(j) &
                    & + gl%dndardnidxi_ALL(4,i,j)
                do m = 1,gl%ncomp-1
                    gl%ndndardnidnj_ALL(4,i,j) = gl%ndndardnidnj_ALL(4,i,j) - gl%molfractions(m) * gl%dndardnidxi_ALL(4,i,m)
                end do
            end do
        end do
    end if

    !d2_alphar_d_tau2 * tau^2
    if (GETDER(5) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                gl%ndndardnidnj_ALL(5,i,j) = ndardni_ALL(7,i) * nddel_dnj(j) + 2.D0 * ndardni_ALL(5,i) * dndtau_dnjdtau(j) &
                    & + ndardni_ALL(9,i) * ndtau_dnj(j) &
                    & + gl%dndardnidxi_ALL(4,i,j)
                do m = 1,gl%ncomp-1
                    gl%ndndardnidnj_ALL(5,i,j) = gl%ndndardnidnj_ALL(5,i,j) - gl%molfractions(m) * gl%dndardnidxi_ALL(5,i,m)
                end do
            end do
        end do
    end if

    !d2_alphar_d_delta_d_tau * delta * tau
    if (GETDER(6) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                gl%ndndardnidnj_ALL(6,i,j) = ndardni_ALL(10,i) * nddel_dnj(j) + ndardni_ALL(6,i) * dnddel_dnjddel(j) &
                    & + ndardni_ALL(6,i) * dndtau_dnjdtau(j) + ndardni_ALL(7,i) * ndtau_dnj(j) &
                    & + gl%dndardnidxi_ALL(6,i,j)
                do m = 1,gl%ncomp-1
                    gl%ndndardnidnj_ALL(6,i,j) = gl%ndndardnidnj_ALL(6,i,j) - gl%molfractions(m) * gl%dndardnidxi_ALL(6,i,m)
                end do
            end do
        end do
    end if




    !Following derivatives not yet implemented:

    !d3_alphar_d_delta_d_tau2 * delta * tau^2
    if (GETDER(7) == 1) then
        gl%ndndardnidnj_ALL(7,:,:) = 0.D0
    end if

    !d3_alphar_d_delta3 * delta^3
    if (GETDER(8) == 1) then
        gl%ndndardnidnj_ALL(8,:,:) = 0.D0
    end if

    !d3_alphar_d_tau3 * tau^3
    if (GETDER(9) == 1) then
        gl%ndndardnidnj_ALL(9,:,:) = 0.D0
    end if

    !d3_alphar_d_delta2_dtau * delta^2 * tau
    if (GETDER(10) == 1) then
        gl%ndndardnidnj_ALL(10,:,:) = 0.D0
    end if

    !d4_alphar_d_delta4 * delta^4
    if (GETDER(11) == 1) then
        gl%ndndardnidnj_ALL(11,:,:) = 0.D0
    end if

    !d4_alphar_d_delta3_dtau * delta^3*tau
    if (GETDER(12) == 1) then
        gl%ndndardnidnj_ALL(12,:,:) = 0.D0
    end if

    !d4_alphar_d_delta2_dtau2 * delta^2*tau^2
    if (GETDER(13) == 1) then
        gl%ndndardnidnj_ALL(13,:,:) = 0.D0
    end if

    !d4_alphar_d_delta_dtau3 * delta*tau^3
    if (GETDER(14) == 1) then
        gl%ndndardnidnj_ALL(14,:,:) = 0.D0
    end if

    !d4_alphar_d_tau4 * tau^4
    if (GETDER(15) == 1) then
        gl%ndndardnidnj_ALL(15,:,:) = 0.D0
    end if


    end subroutine ndndar_dnidnj_ALL
    !**************************************************************************



    !**************************************************************************
    subroutine dndndar_dnidnj_dxk_ALL(gl,Temp, dens, GETDER)
    !--------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF ALL (UP TO FIRST) DELTA AND TAU
    ! DERIVATIVES OF Eq. 7.47 of the GERG-MONOGRAPH (Kunz et al. (2007))
    ! (n*d(n*dalpha^r/dni)/dnj at constant T,V, and n_m)
    ! ALL DERIVATIVES ALSO WITH RESPECT TO xk
    !
    ! The calculated derivatives can be written:
    ! d(n+q)/d(del^n)d(tau^q) d(n*d(n*dalpha^r/dni)/dnj)/dxk
    ! AT CONSTANT TAU, DEL, AND x
    !
    ! The derivatives are given in the supplementary material of the article by
    ! Bell and Jäger (XXX), page 9-11, Eq. (58-60)
    !--------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    ! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
    !               ALL DERIVATIVES ALSO WITH RESPECT TO xj
    !                1. RESIDUAL PART OF THE CHEMICAL POTENTIAL (CPR)
    !                2. 1ST DERIVATIVE OF CPR WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
    !                3. 2ND DERIVATIVE OF CPR WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2 --NOT YET IMPLEMENTED
    !                4. 1ST DERIVATIVE OF CPR WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
    !               -- NOT YET IMPLEMENTED (NOT NEEDED SO FAR)
    !                5: 2ND DERIVATIVE OF CPR WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
    !                6: 1ST MIXED DERIVATIVE OF CPR WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
    !                7: 3RD MIXED DERIVATIVE OF CPR WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
    !                8: 3RD DERIVATIVE OF CPR WITH RESPECT TO D, MULTIPLIED BY DEL^3
    !                9: 3RD DERIVATIVE OF CPR WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !               10: 3RD MIXED DERIVATIVE OF CPR WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    !               11: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^4, MULTIPLIED BY DEL^4
    !               12: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^3 and TAU, MULTIPLIED BY DEL^3*TAU
    !               13: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL^2 and TAU^2, MULTIPLIED BY DEL^2*TAU^2
    !               14: 4TH DERIVATIVE OF CPR WITH RESPECT TO DEL and TAU^3, MULTIPLIED BY DEL*TAU^3
    !               15: 4TH DERIVATIVE OF CPR WITH RESPECT TO TAU^4, MULTIPLIED BY TAU^4
    !
    ! OUTPUT PARAMETERS:
    ! dndndardnidnjdxk_ALL    - multidimensional matrix containing the derivatives
    !--------------------------------------------------------------------------
    !Andreas Jäger, March 2016





    implicit none

    type(type_gl) :: gl


    Double precision:: Temp, dens
    integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed
    integer, dimension(15):: GETDER_help       ! Specifies which help derivatives are calculated

    !Other derivatives needed for the calculation of the derivatives of Eq. 7.47 of GERG (also with respect to xj)
    Double precision, dimension(15,30):: ndardni_ALL
    Double precision, dimension(30):: nddel_dnj, ndtau_dnj
    Double precision, dimension(30):: dnddel_dnjddel, dndtau_dnjdtau
    double precision, dimension(30):: ndTred_dni, ndrhored_dni
    double precision, dimension(30,30):: dndTr_dnidxj, dndrhor_dnidxj
    double precision, dimension(30,30):: dndtau_dnidxj, dnddel_dnidxj
    double precision, dimension(30,30):: d2ndtau_dnidxjdtau, d2nddel_dnidxjddel
    double precision, dimension(30):: dTred_dxi, drhored_dxi

    double precision:: delta, del2, del3, del4, tau, tau2, tau3, tau4
    integer:: i,j,k,m

    if(.not. allocated(gl%dndardnidxi_ALL)) allocate(gl%dndardnidxi_ALL(15,30,30))
    if(.not. allocated(gl%d2ndardnidxidxj_ALL)) allocate(gl%d2ndardnidxidxj_ALL(15,30,30,30))


    !FOR CUBICS ONLY AT THE MOMENT!!!
    GETDER_help = 1 !Calculate all help derivatives
    call ndar_dni_ALL(gl,Temp, Dens, GETDER_help, ndardni_ALL)
    call dndar_dni_dxi_ALL(gl,Temp, Dens, GETDER_help)
    call d2ndar_dni_dxidxj_ALL(gl,Temp, Dens, GETDER_help)
    call ndYr_dni(gl,ndTred_dni, ndrhored_dni)
    call dndYr_dnidxj(gl,dndTr_dnidxj, dndrhor_dnidxj)
    call dYr_dxi(gl,dTred_dxi, drhored_dxi)


    delta = dens / gl%rhoredmix
    del2 = delta*delta
    del3 = del2*delta
    del4 = del2*del2

    tau = gl%tredmix / temp
    tau2 = tau * tau
    tau3 = tau2 * tau
    tau4 = tau2 * tau2



    ndtau_dnj = (tau / gl%tredmix * ndTred_dni) / tau                      !Divided by tau, because all tau derivatives are already multiplied by tau
    nddel_dnj = (delta - delta / gl%rhoredmix * ndrhored_dni) / delta      !Divided by delta, because all delta derivatives are already multiplied by delta

    !Derivatives of ndtau_dnj and nddel_dnj with respect to tau and delta
    dndtau_dnjdtau = 1.D0 / gl%tredmix * ndTred_dni
    dnddel_dnjddel = 1.D0 - 1.D0 / gl%rhoredmix * ndrhored_dni
    !Derivatives of ndtau_dnj and nddel_dnj with respect to xk
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp-1
            dnddel_dnidxj(i,j) = -delta / gl%rhoredmix * (dndrhor_dnidxj(i,j) - 1.D0 / gl%rhoredmix * drhored_dxi(j) * ndrhored_dni(i))
            dndtau_dnidxj(i,j) = -tau / gl%tredmix * (dndTr_dnidxj(i,j) - 1.D0 / gl%tredmix * dTred_dxi(j) * ndTred_dni(i))
        end do
    end do
    dnddel_dnidxj = dnddel_dnidxj / delta           !Divided by delta, because all delta derivatives are already multiplied by delta
    dndtau_dnidxj = dndtau_dnidxj / tau             !Divided by tau, because all tau derivatives are already multiplied by tau
    !Derivatives of ndtau_dnj and nddel_dnj with respect to xk and tau and delta
    d2ndtau_dnidxjdtau = dndtau_dnidxj / tau        !Divided by tau, because all tau derivatives are already multiplied by tau
    d2nddel_dnidxjddel = dnddel_dnidxj / delta      !Divided by delta, because all delta derivatives are already multiplied by delta

    gl%dndndardnidnjdxk_ALL = 0.D0

    !alphar
    if (GETDER(1) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                do k = 1, gl%ncomp-1
                    gl%dndndardnidnjdxk_ALL(1,i,j,k) = ndardni_ALL(2,i) * dnddel_dnidxj(j,k) + gl%dndardnidxi_ALL(2,i,k) * nddel_dnj(j) &
                        & + ndardni_ALL(4,i) * dndtau_dnidxj(j,k) + gl%dndardnidxi_ALL(4,i,k) * ndtau_dnj(j) &
                        & + gl%d2ndardnidxidxj_ALL(1,i,j,k) - gl%dndardnidxi_ALL(1,i,k)
                    do m = 1, gl%ncomp-1
                        gl%dndndardnidnjdxk_ALL(1,i,j,k) = gl%dndndardnidnjdxk_ALL(1,i,j,k) - gl%molfractions(m) * gl%d2ndardnidxidxj_ALL(1,i,k,m)
                    end do
                end do
            end do
        end do
    end if

    !d_alphar_d_delta * delta
    if (GETDER(2) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                do k = 1, gl%ncomp-1
                    gl%dndndardnidnjdxk_ALL(2,i,j,k) = ndardni_ALL(2,i) * d2nddel_dnidxjddel(j,k) + ndardni_ALL(3,i) * dnddel_dnidxj(j,k) &
                        & + gl%dndardnidxi_ALL(2,i,k) * dnddel_dnjddel(j) + gl%dndardnidxi_ALL(3,i,k) * nddel_dnj(j) &
                        & + gl%dndardnidxi_ALL(6,i,k) * dndtau_dnidxj(j,k) + gl%dndardnidxi_ALL(6,i,k) * ndtau_dnj(j) &
                        & + gl%d2ndardnidxidxj_ALL(2,i,j,k) - gl%dndardnidxi_ALL(2,i,k)
                    do m = 1, gl%ncomp-1
                        gl%dndndardnidnjdxk_ALL(2,i,j,k) = gl%dndndardnidnjdxk_ALL(2,i,j,k) - gl%molfractions(m) * gl%d2ndardnidxidxj_ALL(2,i,k,m)
                    end do
                end do
            end do
        end do
    end if

    !d2_alphar_d_delta2 * delta^2
    if (GETDER(3) == 1) then
        !NOT YET IMPLEMENTED
        gl%dndndardnidnjdxk_ALL(3,:,:,:) = 0.D0
    end if

    !d_alphar_d_tau * tau
    if (GETDER(4) == 1) then
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                do k = 1, gl%ncomp-1
                    gl%dndndardnidnjdxk_ALL(4,i,j,k) = ndardni_ALL(6,i) * dnddel_dnidxj(j,k) + gl%dndardnidxi_ALL(6,i,k) * nddel_dnj(j) &
                        & + ndardni_ALL(4,i) * d2ndtau_dnidxjdtau(j,k) + ndardni_ALL(5,i) * dndtau_dnidxj(j,k) &
                        & + gl%dndardnidxi_ALL(4,i,k) * dndtau_dnjdtau(j) + gl%dndardnidxi_ALL(5,i,k) * ndtau_dnj(j) &
                        & + gl%d2ndardnidxidxj_ALL(4,i,j,k) - gl%dndardnidxi_ALL(4,i,k)
                    do m = 1, gl%ncomp-1
                        gl%dndndardnidnjdxk_ALL(4,i,j,k) = gl%dndndardnidnjdxk_ALL(4,i,j,k) - gl%molfractions(m) * gl%d2ndardnidxidxj_ALL(4,i,k,m)
                    end do
                end do
            end do
        end do
    end if





    !Following derivatives not yet implemented:

    !d2_alphar_d_tau2 * tau^2
    if (GETDER(5) == 1) then
        gl%dndndardnidnjdxk_ALL(5,:,:,:) = 0.D0
    end if

    !d2_alphar_d_delta_d_tau * delta * tau
    if (GETDER(6) == 1) then
        gl%dndndardnidnjdxk_ALL(6,:,:,:) = 0.D0
    end if

    !d3_alphar_d_delta_d_tau2 * delta * tau^2
    if (GETDER(7) == 1) then
        gl%dndndardnidnjdxk_ALL(7,:,:,:) = 0.D0
    end if

    !d3_alphar_d_delta3 * delta^3
    if (GETDER(8) == 1) then
        gl%dndndardnidnjdxk_ALL(8,:,:,:) = 0.D0
    end if

    !d3_alphar_d_tau3 * tau^3
    if (GETDER(9) == 1) then
        gl%dndndardnidnjdxk_ALL(9,:,:,:) = 0.D0
    end if

    !d3_alphar_d_delta2_dtau * delta^2 * tau
    if (GETDER(10) == 1) then
        gl%dndndardnidnjdxk_ALL(10,:,:,:) = 0.D0
    end if

    !d4_alphar_d_delta4 * delta^4
    if (GETDER(11) == 1) then
        gl%dndndardnidnjdxk_ALL(11,:,:,:) = 0.D0
    end if

    !d4_alphar_d_delta3_dtau * delta^3*tau
    if (GETDER(12) == 1) then
        gl%dndndardnidnjdxk_ALL(12,:,:,:) = 0.D0
    end if

    !d4_alphar_d_delta2_dtau2 * delta^2*tau^2
    if (GETDER(13) == 1) then
        gl%dndndardnidnjdxk_ALL(13,:,:,:) = 0.D0
    end if

    !d4_alphar_d_delta_dtau3 * delta*tau^3
    if (GETDER(14) == 1) then
        gl%dndndardnidnjdxk_ALL(14,:,:,:) = 0.D0
    end if

    !d4_alphar_d_tau4 * tau^4
    if (GETDER(15) == 1) then
        gl%dndndardnidnjdxk_ALL(15,:,:,:) = 0.D0
    end if


    end subroutine dndndar_dnidnj_dxk_ALL
    !**************************************************************************


    !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    !END OF NEW DERIVATIVES NEEDED FOR CRITICAL POINT CALCULATION OF A MIXTURE
    !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





    !Derivatives of UNIFAC with respect to mole frations
    !Needed for phase equilibrium routines for PSRK(gl,and VTPR and new Helmholtz+gE) model
    !------------------------------------------------------------------------------------------------------------------------------
    !subroutine gE_UNIFAC_MIXDERIVS_dxa(gl,Temp, GETDER, gl%ge%gE_C_dxa, gl%ge%gE_R_dxa, gl%ge%ln_gamma_C_dxa, gl%ge%ln_gamma_R_dxa, C_or_R, errval)
    subroutine gE_UNIFAC_MIXDERIVS_dxa(gl,Temp, GETDER, C_or_R, errval)
    !------------------------------------------------------------------------------------------------------------------------------
    !Subroutine to calculate the excess Gibbs energy and the activity coefficients of the
    !components in the mixture according to the UNIFAC model as well as tau- and delta-
    !derivatives of the model.
    !All quantities are given as derivatives with respect to mole fractions xa for all components
    !in the mixture. The convention that the "last" mole fraction xN of component N is replaced with
    !1-sum(xk) is used here.
    !
    !INPUT
    !Temp       -    Temperature in K
    !Note that the model makes use of the global variable "molfractions" defined in the module
    !"module_fluid_parameters"
    ! GETDER      - AN ARRAY WITH nderivs ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED. !!!ALL DERIVATIVES ARE ADDITIONALLY DERIVATIVES WITH RESPECT TO xa!!!:
    !                1. gE from UNIFAC as a function of temperature and composition
    !                2. 1ST DERIVATIVE OF gE WITH RESPECT TO DEL AT CONSTANT TAU
    !                3. 2ND DERIVATIVE OF gE WITH RESPECT TO DEL AT CONSTANT TAU
    !                4. 1ST DERIVATIVE OF gE WITH RESPECT TO TAU AT CONSTANT DEL
    !                5: 2ND DERIVATIVE OF gE WITH RESPECT TO TAU AT CONSTANT DEL
    !                6: 2ND MIXED DERIVATIVE OF gE WITH RESPECT TO DEL AND TAU
    !                7: 3RD MIXED DERIVATIVE OF gE WITH RESPECT TO DEL, TAU, AND TAU
    !                8: 3RD DERIVATIVE OF gE WITH RESPECT TO DEL
    !                9: 3RD DERIVATIVE OF gE WITH RESPECT TO TAU
    !               10: 3RD MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, DEL, AND DEL
    !               11: 4TH DERIVATIVE OF gE WITH RESPECT TO DEL
    !               12: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, DEL, DEL AND DEL
    !               13: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, TAU, DEL, AND DEL
    !               14: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, TAU, TAU AND DEL
    !               15: 4TH DERIVATIVE OF gE WITH RESPECT TO TAU
    !C_or_R       - Integer that specifies if the combinatorial, the residual or both parts shall be calculated
    !               C_or_R = 0: Calculate both parts
    !               C_or_R = 1: Calculate the combinatorial part only
    !               C_or_R = 2: Calculate the the residual part only
    !
    !OUTPUT
    !gE_C         -    Matrix with molar excess Gibbs energy of the combinatorial part (C) in J/molK and derivatives with respect to tau and delta and xa
    !gE_R         -    Matrix with molar excess Gibbs energy of the residual part (R) in J/molK and derivatives with respect to tau and delta and xa
    !ln_gamma_C   -    Matrix containing the natural logarithm of the activity coefficients of the combinatorial part for each component and derivatives with respect to tau and delta and xa
    !ln_gamma_R   -    Matrix containing the natural logarithm of the activity coefficients of the residual part for each component and derivatives with respect to tau and delta and xa
    !errval       -    Indicates if an error occurred during calculations
    !
    !Andreas Jäger, January 2017
    !------------------------------------------------------------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl


    double precision :: Temp
    integer, dimension(nderivs):: GETDER                !array specifier to indicate, which derivative is needed
    !double precision, dimension(nderivs, 30) :: gl%ge%gE_C_dxa        !Combinatorial part of gE and derivatives with respect to tau and delta and xa
    !double precision, dimension(nderivs, 30) :: gl%ge%gE_R_dxa        !Residual part of gE and derivatives with respect to tau and delta and xa
    !double precision, allocatable:: gl%ge%ln_gamma_C_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to delta and tau and xa
    !double precision, allocatable:: gl%ge%ln_gamma_R_dxa(:,:,:)       !Residual activity coefficients and derivatives with respect to delta and tau and xa
    integer:: C_or_R
    integer:: errval

    double precision, dimension(30) :: r_i              !Molecular vdW volume of components i
    double precision, dimension(30) :: q_i              !Molecular vdW surface area of components i

    double precision :: R_const                         !universal gas constant

    !Help variables
    double precision, dimension(30) :: x
    double precision, dimension(30) :: phi_i
    double precision, dimension(30) :: Theta_i
    !double precision, dimension(100) :: lnGamma_k
    !double precision, dimension(100,30) :: lnGamma_ki   !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    double precision, dimension(100) :: Theta_m
    double precision, dimension(:,:), allocatable :: Theta_mi     !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    double precision, dimension(100) :: X_m             !Mole fraction of groups in the mixture
    double precision, dimension(:,:), allocatable :: X_mi         !Mole fraction of groups in pure component i. !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    double precision, dimension(:,:), allocatable :: Psi_nm
    double precision, dimension(100):: Q_m              !Help variable to list all group areas in one vector in order to not have to map them back to Q_ik(i, k) for derivative calculation

    !More help variables
    double precision:: sum_rjxj
    double precision:: sum_qjxj
    double precision:: sum_QnXn
    double precision, dimension(30):: sum_QndXn_dxa
    double precision:: sumsum_vnxj
    !double precision:: sum_xjLj
    double precision:: sum_ThetamPsimk
    !double precision:: sumsum_ThetamPsikm
    double precision:: sum_ThetanPsinm

    !Derivatives of the combinatorial part with respect to "natural" variables
    double precision, dimension(nderivs, 30) :: ln_gamma_C_Trho       !Combinatorial activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs, 30) :: ln_gamma_R_Trho       !Residual activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs) :: gE_C_Trho        !Combinatorial part of gE and derivatives with respect to T and rho
    double precision, dimension(nderivs) :: gE_R_Trho        !Residual part of gE and derivatives with respect to T and rho
    double precision, dimension(:,:,:), allocatable:: ln_gamma_C_Trho_dxa       !Combinatorial activity coefficients and derivatives with respect to rho and T, and xa
    double precision, dimension(:,:,:), allocatable:: ln_gamma_R_Trho_dxa       !Residual activity coefficients and derivatives with respect to rho and T, and xa
    double precision, dimension(nderivs, 30) :: gE_C_Trho_dxa        !Combinatorial part of gE and derivatives with respect to T and rho and xa
    double precision, dimension(nderivs, 30) :: gE_R_Trho_dxa        !Residual part of gE and derivatives with respect to T and rho and xa

    !Derivatives of the residual part with respect to "natural" variables
    !For combinatorial part
    double precision, dimension(:,:), allocatable :: dphi_i_dxa
    double precision, dimension(:,:), allocatable :: dTheta_i_dxa
    !For residual part
    double precision, dimension(:,:), allocatable :: dlnGamma_k_dxa
    double precision, dimension(:,:), allocatable :: d2lnGamma_k_dxadT
    double precision, dimension(:,:), allocatable :: dPsi_nm_dT
    double precision, dimension(:,:), allocatable :: dTheta_m_dxa
    double precision, dimension(:,:), allocatable :: dX_m_dxa
    double precision:: sum_dThetam_dxa_Psimk
    double precision:: sum_dThetan_dxa_Psinm
    double precision:: sumsum_ThetamPsikm_dxa          !Variable for double sum in xa derivative
    double precision:: sum_Thetam_dPsimk_dT
    double precision:: sum_Thetan_dPsinm_dT
    double precision:: sum_dThetam_dxa_dPsimk_dT
    double precision:: sum_dThetan_dxa_dPsinm_dT
    double precision:: sumsum_ThetamPsikm_dxadT        !Variable for double sum in xa and temperature derivative

    double precision, dimension(30):: dTred_dxa, drhored_dxa    !Variables for the derivatives of the reducing functions with respect to molfractions xa

    !Required derivatives of activities and excess Gibbs energy only with respect to tau and delta (not xa, however, needed for the calculation of derivatives with respect to xa)
    integer, dimension(nderivs):: GETDER_no_xa                !array specifier to indicate, which derivative is needed
    !double precision, dimension(nderivs) :: gl%ge%gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gl%ge%gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau


    !Corresponding states relevant variables
    double precision:: tau, del

    !Summation variables
    integer:: a, i, j, k, m, n, help_k

    !Variable that indicates if derivatives need to be recalculated
    logical:: recalc


    !allocation of large arrays
    if(.not. allocated(gl%ge%ln_gamma_C_dxa))   allocate(gl%ge%ln_gamma_C_dxa(15,30,30))                                                                                                                  ! allocate(gl%ge%ln_gamma_C_dxa(nderivs,gl%ncomp,gl%ncomp))
    if(.not. allocated(gl%ge%ln_gamma_R_dxa))   allocate(gl%ge%ln_gamma_R_dxa(15,30,30))                                                                                                                  ! allocate(gl%ge%ln_gamma_R_dxa(nderivs,gl%ncomp,gl%ncomp))
    if(.not. allocated(ln_gamma_C_Trho_dxa))    allocate(ln_gamma_C_Trho_dxa(15,30,30))                                                                                                                   ! allocate(ln_gamma_C_Trho_dxa(nderivs,gl%ncomp,gl%ncomp))
    if(.not. allocated(ln_gamma_R_Trho_dxa))    allocate(ln_gamma_R_Trho_dxa(15,30,30))                                                                                                                   ! allocate(ln_gamma_R_Trho_dxa(nderivs,gl%ncomp,gl%ncomp))
    !
    if(.not. allocated(Theta_mi))               allocate(Theta_mi(100,30))                                                                        ! allocate(Theta_mi((gl%ncomp*gl%nr_of_groups_mix-gl%nr_of_groups_mix)*gl%nr_of_groups_mix,gl%ncomp))
    if(.not. allocated(X_mi))                   allocate(X_mi(100,30))                                                                            ! allocate(X_mi((gl%ncomp*gl%nr_of_groups_mix-gl%nr_of_groups_mix)*gl%nr_of_groups_mix,gl%ncomp))
    if(.not. allocated(Psi_nm))                 allocate(Psi_nm(100,100))            ! allocate(Psi_nm((gl%ncomp*gl%nr_of_groups_mix-gl%nr_of_groups_mix)*gl%nr_of_groups_mix,(gl%ncomp*gl%nr_of_groups_mix-gl%nr_of_groups_mix)*gl%nr_of_groups_mix))
    !
    if(.not. allocated(dPsi_nm_dT))             allocate(dPsi_nm_dT(100,100))        ! allocate(dPsi_nm_dT((gl%ncomp*gl%nr_of_groups_mix-gl%nr_of_groups_mix)*gl%nr_of_groups_mix,(gl%ncomp*gl%nr_of_groups_mix-gl%nr_of_groups_mix)*gl%nr_of_groups_mix))
    !
    if(.not. allocated(dphi_i_dxa))             allocate(dphi_i_dxa(30,30))                                                                                                                                    ! allocate(dphi_i_dxa(gl%ncomp,gl%ncomp))
    if(.not. allocated(dTheta_i_dxa))           allocate(dTheta_i_dxa(30,30))                                                                                                                                  ! allocate(dTheta_i_dxa(gl%ncomp,gl%ncomp))
    !
    !
    !
    if(.not. allocated(dlnGamma_k_dxa))         allocate(dlnGamma_k_dxa(100,30))                                                                  ! allocate(dlnGamma_k_dxa((gl%ncomp*gl%nr_of_groups_mix-gl%nr_of_groups_mix)*gl%nr_of_groups_mix,gl%ncomp))
    if(.not. allocated(d2lnGamma_k_dxadT))      allocate(d2lnGamma_k_dxadT(100,30))                                                               ! allocate(d2lnGamma_k_dxadT((gl%ncomp*gl%nr_of_groups_mix-gl%nr_of_groups_mix)*gl%nr_of_groups_mix,gl%ncomp))
    if(.not. allocated(dTheta_m_dxa))           allocate(dTheta_m_dxa(100,30))                                                                    ! allocate(dTheta_m_dxa((gl%ncomp*gl%nr_of_groups_mix-gl%nr_of_groups_mix)*gl%nr_of_groups_mix,gl%ncomp))
    if(.not. allocated(dX_m_dxa))               allocate(dX_m_dxa(100,30))                                                                        ! allocate(dX_m_dxa((gl%ncomp*gl%nr_of_groups_mix-gl%nr_of_groups_mix)*gl%nr_of_groups_mix,gl%ncomp))





    errval = 0

    !Check for wrong input to C_or_R
    if ((C_or_R < 0) .or. (C_or_R > 2)) then
        errval = -1111
        gl%ge%gE_C_dxa = errval
        gl%ge%gE_R_dxa = errval
        gl%ge%ln_gamma_C_dxa = errval
        gl%ge%ln_gamma_R_dxa = errval
        return
    end if

    recalc = .true.

    !Check whether calculations have to be redone or not.
    !If the temperature, (density), composition, equation of state and, the required part of UNIFAC stayed the same, no need to recalculate
    !if(dabs(gl%ge%Temp_dxa_prev - Temp) > 1.D-16) then
    !   recalc = .true.
    !end if
    !do i=1,gl%ncomp
    !    if(dabs(gl%ge%molfrac_dxa_prev(i) - gl%molfractions(i)) > 1.D-16) then
    !       recalc = .true.
    !    end if
    !    if(gl%ge%Eq_type_dxa_prev(i) .ne. gl%Eq_type(i)) then
    !       recalc = .true.
    !    end if
    !end do
    !if (gl%ge%mixtype_dxa_prev .ne. gl%mix_type) then
    !    recalc = .true.
    !end if
    !if (gl%ge%C_or_R_dxa_prev .ne. C_or_R) then
    !    recalc = .true.
    !end if


    !Initialize variables
    ln_gamma_C_Trho = 0.D0
    ln_gamma_R_Trho = 0.D0
    gE_C_Trho = 0.D0
    gE_R_Trho = 0.d0
    ln_gamma_C_Trho_dxa = 0.D0
    ln_gamma_R_Trho_dxa = 0.d0
    gE_C_Trho_dxa = 0.D0
    gE_R_Trho_dxa = 0.D0
    gl%ge%ln_gamma_C_dxa = 0.D0
    gl%ge%ln_gamma_R_dxa = 0.D0
    gl%ge%gE_C_dxa = 0.D0
    gl%ge%gE_R_dxa = 0.D0

    !Value for the ideal gas constant
    R_const = 8.3144598D0

    !Get mole fractions of the mixture from module variable
    x = gl%molfractions

    !Calculate tau
    tau = gl%tredmix / Temp

    !Set dummy value for delta (Not needed because UNIFAC is not a function of density. However, already implemented here for possible later modifications)
    del = 1.D0

    !Derivatives of tredmix with respect to xa are required for the following derivatives
    call dYr_dxi(gl,dTred_dxa, drhored_dxa)


    !Get the required derivatives of activities and excess Gibbs energies with respect tau (and not with respect to xa, however, needed for the derivatives here)
    GETDER_no_xa = 0
    if (GETDER(1) .eq. 1)  then
        GETDER_no_xa(1) = 1
        GETDER_no_xa(4) = 1
    end if
    if (GETDER(4) .eq. 1)  then
        GETDER_no_xa(1) = 1
        GETDER_no_xa(4) = 1
        GETDER_no_xa(5) = 1
    end if
    !call gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER_no_xa, gl%ge%gE_C, gl%ge%gE_R, gl%ge%ln_gamma_C, gl%ge%ln_gamma_R, C_or_R, errval)
    call gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER_no_xa, C_or_R, errval)
    if (GETDER(1) .eq. 1) then
        !Transform the tau derivatives into the needed temperature derivatives
        ln_gamma_C_Trho(1,:) = gl%ge%ln_gamma_C(1,:)
        gE_C_Trho(1) = gl%ge%gE_C(1)
        ln_gamma_R_Trho(1,:) = gl%ge%ln_gamma_R(1,:)
        gE_R_Trho(1) = gl%ge%gE_R(1)
        ln_gamma_C_Trho(4,:) = - tau**2 / gl%tredmix * gl%ge%ln_gamma_C(4,:)
        gE_C_Trho(4) = - tau**2 / gl%tredmix * gl%ge%gE_C(4)
        ln_gamma_R_Trho(4,:) = - tau**2 / gl%tredmix * gl%ge%ln_gamma_R(4,:)
        gE_R_Trho(4) = - tau**2 / gl%tredmix * gl%ge%gE_R(4)
    end if
    if (GETDER(4) .eq. 1) then
        !Transform the tau derivatives into the needed temperature derivatives
        ln_gamma_C_Trho(1,:) = gl%ge%ln_gamma_C(1,:)
        gE_C_Trho(1) = gl%ge%gE_C(1)
        ln_gamma_R_Trho(1,:) = gl%ge%ln_gamma_R(1,:)
        gE_R_Trho(1) = gl%ge%gE_R(1)
        ln_gamma_C_Trho(4,:) = - tau**2 / gl%tredmix * gl%ge%ln_gamma_C(4,:)
        gE_C_Trho(4) = - tau**2 / gl%tredmix * gl%ge%gE_C(4)
        ln_gamma_R_Trho(4,:) = - tau**2 / gl%tredmix * gl%ge%ln_gamma_R(4,:)
        gE_R_Trho(4) = - tau**2 / gl%tredmix * gl%ge%gE_R(4)
        ln_gamma_C_Trho(5,:) = tau**4 / gl%tredmix**2 * gl%ge%ln_gamma_C(5,:) - 2.D0 * tau / gl%tredmix * ln_gamma_C_Trho(4,:)
        gE_C_Trho(5) = tau**4 / gl%tredmix**2 * gl%ge%gE_C(5) - 2.D0 * tau / gl%tredmix * gE_C_Trho(4)
        ln_gamma_R_Trho(5,:) = tau**4 / gl%tredmix**2 * gl%ge%ln_gamma_R(5,:) - 2.D0 * tau / gl%tredmix * ln_gamma_R_Trho(4,:)
        gE_R_Trho(5) = tau**4 / gl%tredmix**2 * gl%ge%gE_R(5) - 2.D0 * tau / gl%tredmix * gE_R_Trho(4)
    end if


    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 1)) then
        !Calculate the combinatorial part of the activity coefficients
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !ln_gamma_C is always needed for the following calculations, thus it is always calculated

        if ((gl%ge%GETDER_dxa_prev(1) == 1) .and. (recalc .eqv. .false.)) then

            gl%ge%ln_gamma_C_dxa(1,:,:) = gl%ge%ln_gamma_C_dxa_prev(1,:,:)

            !Transformation of variables from tau,delta, and x to T,rho,x
            do a = 1, gl%ncomp-1
                do i = 1, gl%ncomp
                    ln_gamma_C_Trho_dxa(1,i,a) =  gl%ge%ln_gamma_C_dxa(1,i,a) - dTred_dxa(a) / tau * ln_gamma_C_Trho(4,i)
                end do
            end do

        else

            !Calculate help variables
            r_i = 0.D0
            q_i = 0.D0
            Do i = 1, gl%ncomp
                Do k = 1, gl%nr_of_groups_i(i)
                    r_i(i) = r_i(i) + gl%v_ik(i, k) * gl%R_ik(i, k)
                    q_i(i) = q_i(i) + gl%v_ik(i, k) * gl%Q_ik(i, k)
                End do
            End do

            sum_rjxj = 0.D0
            sum_qjxj = 0.D0
            Do j = 1, gl%ncomp
                sum_rjxj = sum_rjxj + r_i(j) * x(j)
                sum_qjxj = sum_qjxj + q_i(j) * x(j)
            End do
            Do i = 1, gl%ncomp
                phi_i(i) = r_i(i) / sum_rjxj
                Theta_i(i) = q_i(i) / sum_qjxj
            End do
            !Calculation of all activity coefficients of the combinatorial part
            !(Not needed anmyore because called from gE_UNIFAC_MIXDERIVS above)
            !Do i = 1, ncomp
            !    ln_gamma_C_Trho(1,i) = 1.D0 - phi_i(i) + dlog(phi_i(i)) - &
            !                    &  5.D0 * q_i(i) * (1 - phi_i(i) / Theta_i(i) + dlog(phi_i(i) / Theta_i(i)))
            !End do

            !Calculate the derivativs of phi_i and Theta_i with respect to xa. These derivatives are also needed for all tau and delta derivatives
            Do a = 1, gl%ncomp-1
                Do i = 1, gl%ncomp
                    dphi_i_dxa(i,a) = -r_i(i) * (r_i(a) - r_i(gl%ncomp)) / sum_rjxj**2
                    dTheta_i_dxa(i,a) = -q_i(i) * (q_i(a) - q_i(gl%ncomp)) / sum_qjxj**2
                End do
            End do

            !Compute the mole fraction derivatives of all activity coefficients of the combinatorial part
            do a = 1, gl%ncomp-1
                Do i = 1, gl%ncomp
                    ln_gamma_C_Trho_dxa(1,i,a) = (1.D0 / phi_i(i) - 1.D0) * dphi_i_dxa(i,a) &
                        & - 5.D0 * q_i(i) * ((Theta_i(i)/phi_i(i) - 1.D0) &
                        & * (Theta_i(i)*dphi_i_dxa(i,a) - phi_i(i)*dTheta_i_dxa(i,a)) / Theta_i(i)**2)
                End do
            End do

            !Transformation of variables from T,rho,x, to tau,delta, and x
            do a = 1, gl%ncomp-1
                do i = 1, gl%ncomp
                    gl%ge%ln_gamma_C_dxa(1,i,a) = dTred_dxa(a) / tau * ln_gamma_C_Trho(4,i) + ln_gamma_C_Trho_dxa(1,i,a)  !The ln_gamma_C_Trho(4,i) is the first derivative of ln(gammaiC) with respect to temperature. As gammaiC is not a function of temperature, this derivative is 0
                end do
            end do

        end if

        if ((GETDER(1) .eq. 1) .or. (GETDER(4) .eq. 1) .or. (GETDER(5) .eq. 1)) then

            if ((gl%ge%GETDER_dxa_prev(1) == 1) .and. (recalc .eqv. .false.)) then

                gl%ge%gE_C_dxa(1,:) = gl%ge%gE_C_dxa_prev(1,:)
                !Transformation of variables from tau,delta,x to T,rho,x
                do a=1,gl%ncomp-1
                    gE_C_Trho_dxa(1,a) = gl%ge%gE_C_dxa(1,a) - dTred_dxa(a) / tau * gE_C_Trho(4)
                end do

            else

                Do a = 1,gl%ncomp-1
                    Do i = 1, gl%ncomp
                        gE_C_Trho_dxa(1,a) = gE_C_Trho_dxa(1,a) + R_const * Temp * x(i) * ln_gamma_C_Trho_dxa(1,i,a)
                    end do
                    gE_C_Trho_dxa(1,a) = gE_C_Trho_dxa(1,a) + R_const * Temp * (ln_gamma_C_Trho(1,a) - ln_gamma_C_Trho(1,gl%ncomp))
                End do
                !Calculate the derivative of gE_C with respect to temperature T at constant rho and x. Also needed for derivatives with respect to tau and xa
                !(Not needed anmyore because called from gE_UNIFAC_MIXDERIVS above)
                !do i=1,ncomp
                !    gE_C_Trho(4) = gE_C_Trho(4) + x(i) * ln_gamma_C_Trho(1,i) * R_const
                !end do
                !Transformation of variables from T,rho,x, to tau,delta, and x
                do a=1,gl%ncomp-1
                    gl%ge%gE_C_dxa(1,a) =  dTred_dxa(a) / tau * gE_C_Trho(4) + gE_C_Trho_dxa(1,a)
                end do

            end if
        end if

        !All other derivatives of the combinatorial part are 0 (except for the tau-only derivatives of gE), because the combinatorial part is a function of composition only
        !Derivative wrt delta
        if (GETDER(2) .eq. 1) then
            ln_gamma_C_Trho_dxa(2,:,:) = 0.D0
            gE_C_Trho_dxa(2,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(2,:,:) = 0.d0
            gl%ge%gE_C_dxa(2,:) = 0.D0
        end if
        !Second derivative wrt delta
        if (GETDER(3) .eq. 1) then
            ln_gamma_C_Trho_dxa(3,:,:) = 0.D0
            gE_C_Trho_dxa(3,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(3,:,:) = 0.d0
            gl%ge%gE_C_dxa(3,:) = 0.D0
        end if

        !Derivative wrt tau
        if ((GETDER(4) .eq. 1) .or. (GETDER(5) .eq. 1)) then

            if ((gl%ge%GETDER_dxa_prev(4) == 1) .and. (recalc .eqv. .false.)) then

                gl%ge%ln_gamma_C_dxa(4,:,:) = gl%ge%ln_gamma_C_dxa_prev(4,:,:)
                gl%ge%gE_C_dxa(4,:) = gl%ge%gE_C_dxa_prev(4,:)

                !Transformation of variables from tau,delta,x to T,rho,x
                do a = 1, gl%ncomp -1
                    do i = 1, gl%ncomp
                        ln_gamma_C_Trho_dxa(4,i,a) = -dTred_dxa(a) / gl%tredmix * ln_gamma_C_Trho(4,i) &
                            & - dTred_dxa(a) / tau * ln_gamma_C_Trho(5,i) &
                            & - tau**2 / gl%tredmix * gl%ge%ln_gamma_C_dxa(4,i,a)
                    end do
                    gE_C_Trho_dxa(4,a) = -dTred_dxa(a) / gl%tredmix * gE_C_Trho(4) &
                        & - dTred_dxa(a) / tau * dTred_dxa(a)*gE_C_Trho(5) &
                        & - tau**2 / gl%tredmix * gl%ge%gE_C_dxa(4,a)
                end do

            else

                !Derivatives with respect to natural variables (T,rho,x)
                ln_gamma_C_Trho_dxa(4,:,:) = 0.D0   !ln_gamma_C is not a function of temperature
                !Compute second derivative of the combinatorial part of the excess Gibbs energy with respect to xa and temperature
                Do a = 1,gl%ncomp-1
                    Do i = 1, gl%ncomp
                        gE_C_Trho_dxa(4,a) = gE_C_Trho_dxa(4,a) + R_const * x(i) * ln_gamma_C_Trho_dxa(1,i,a) &
                            & +  R_const * Temp * x(i) * ln_gamma_C_Trho_dxa(4,i,a) !This part of the sum is 0, as lngammaC is not a function of temperature
                    end do
                    gE_C_Trho_dxa(4,a) = gE_C_Trho_dxa(4,a) + R_const * (ln_gamma_C_Trho(1,a) - ln_gamma_C_Trho(1,gl%ncomp)) &
                        &  + R_const * Temp * (ln_gamma_C_Trho(4,a) - ln_gamma_C_Trho(4,gl%ncomp)) !This part of the sum is 0, as lngammaC is not a function of temperature
                End do

                !Calculate the derivative of gE_C with respect to temperature T at constant rho and x. Also needed for all higher tau derivatives of gE
                !gE_C_Trho(5) = 0.D0

                do a = 1, gl%ncomp -1
                    do i = 1, gl%ncomp
                        gl%ge%ln_gamma_C_dxa(4,i,a) = -1.D0/tau**2* dTred_dxa(a)*ln_gamma_C_Trho(4,i) &
                            & - gl%tredmix/tau**3*dTred_dxa(a)*ln_gamma_C_Trho(5,i) &
                            & - gl%tredmix/tau**2*ln_gamma_C_Trho_dxa(4,i,a)       !Whole term is 0, because all ln_gamma_C derivatives with respect to temperature are 0
                    end do
                    gl%ge%gE_C_dxa(4,a) = -1.D0/tau**2* dTred_dxa(a)*gE_C_Trho(4) &
                        & - gl%tredmix/tau**3*dTred_dxa(a)*gE_C_Trho(5) &
                        & - gl%tredmix/tau**2*gE_C_Trho_dxa(4,a)
                end do

            end if

        end if

        !Second derivative wrt tau
        if (GETDER(5) .eq. 1) then
            ln_gamma_C_Trho_dxa(5,:,:) = 0.D0
            gE_C_Trho_dxa(5,:) = 0.D0
            gl%ge%ln_gamma_C_dxa(5,:,:) = 0.d0
            gl%ge%gE_C_dxa(5,:) = 0.D0
        end if
        !Second derivative wrt delta and tau
        if (GETDER(6) .eq. 1) then
            ln_gamma_C_Trho_dxa(6,:,:) = 0.D0
            gE_C_Trho_dxa(6,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(6,:,:) = 0.d0
            gl%ge%gE_C_dxa(6,:) = 0.D0
        end if
        !Third derivative wrt delta, tau, and tau
        if (GETDER(7) .eq. 1) then
            ln_gamma_C_Trho_dxa(7,:,:) = 0.D0
            gE_C_Trho_dxa(7,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(7,:,:) = 0.d0
            gl%ge%gE_C_dxa(7,:) = 0.D0
        end if
        !Third derivative wrt delta
        if (GETDER(8) .eq. 1) then
            ln_gamma_C_Trho_dxa(8,:,:) = 0.D0
            gE_C_Trho_dxa(8,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(8,:,:) = 0.d0
            gl%ge%gE_C_dxa(8,:) = 0.D0
        end if
        !Third derivative wrt tau
        if (GETDER(9) .eq. 1) then
            ln_gamma_C_Trho_dxa(9,:,:) = 0.D0
            gE_C_Trho_dxa(9,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(9,:,:) = 0.d0
            gl%ge%gE_C_dxa(9,:) = 0.D0
        end if
        !Third derivative wrt tau, delta, and delta
        if (GETDER(10) .eq. 1) then
            ln_gamma_C_Trho_dxa(10,:,:) = 0.D0
            gE_C_Trho_dxa(10,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(10,:,:) = 0.d0
            gl%ge%gE_C_dxa(10,:) = 0.D0
        end if
        !Fourth derivative wrt delta
        if (GETDER(11) .eq. 1) then
            ln_gamma_C_Trho_dxa(11,:,:) = 0.D0
            gE_C_Trho_dxa(11,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(11,:,:) = 0.d0
            gl%ge%gE_C_dxa(11,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, delta and tau
        if (GETDER(12) .eq. 1) then
            ln_gamma_C_Trho_dxa(12,:,:) = 0.D0
            gE_C_Trho_dxa(12,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(12,:,:) = 0.d0
            gl%ge%gE_C_dxa(12,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, tau, tau
        if (GETDER(13) .eq. 1) then
            ln_gamma_C_Trho_dxa(13,:,:) = 0.D0
            gE_C_Trho_dxa(13,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(13,:,:) = 0.d0
            gl%ge%gE_C_dxa(13,:) = 0.D0
        end if
        !Fourth derivative wrt delta, tau, tau, tau
        if (GETDER(14) .eq. 1) then
            ln_gamma_C_Trho_dxa(14,:,:) = 0.D0
            gE_C_Trho_dxa(14,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(14,:,:) = 0.d0
            gl%ge%gE_C_dxa(14,:) = 0.D0
        end if
        !Fourth derivative wrt tau
        if (GETDER(15) .eq. 1) then
            ln_gamma_C_Trho_dxa(15,:,:) = 0.D0
            gE_C_Trho_dxa(15,:) = 0.d0
            gl%ge%ln_gamma_C_dxa(15,:,:) = 0.d0
            gl%ge%gE_C_dxa(15,:) = 0.D0
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end if


    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 2)) then
        !Calculate the residual part of the activity coefficients
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !ln_gamma_R is always needed for the following calculations, thus it is always calculated

        !----------------------------------------------------------------------
        !Calculate the mole fractions of all groups in the mixture
        !----------------------------------------------------------------------
        sumsum_vnxj = 0.D0
        Do j = 1, gl%ncomp                     !sum over all components in the mixture
            Do n = 1, gl%nr_of_groups_i(j)     !Sum over all groups of component j
                sumsum_vnxj = sumsum_vnxj + gl%v_ik(j, n) * x(j)
            End do
        End do

        m = 1
        Do j = 1, gl%ncomp                     !Sum over all components in the mixture
            Do k = 1, gl%nr_of_groups_i(j)     !Sum over all groups of component j
                X_m(m) = gl%v_ik(j, k) * x(j) / sumsum_vnxj
                m = m + 1
            End do
        End do
        !----------------------------------------------------------------------


        !----------------------------------------------------------------------
        !Calculate Theta_m for all groups in the mixture
        !----------------------------------------------------------------------
        m = 1
        sum_QnXn = 0.D0
        Do j = 1, gl%ncomp                      !Sum over all components in the mixture
            Do k = 1, gl%nr_of_groups_i(j)      !Sum over all groups of component j
                sum_QnXn = sum_QnXn + gl%Q_ik(j, k) * X_m(m)
                m = m + 1
            End do
        End do

        m = 1
        Do j = 1, gl%ncomp                      !Sum over all components in the mixture
            Do k = 1, gl%nr_of_groups_i(j)      !Sum over all groups of component j
                Theta_m(m) = gl%Q_ik(j, k) * X_m(m) / sum_QnXn
                !Set help variable Q_m to list all group areas in one vector in order to not have to map them back to Q_ik(i, k) for derivative calculation
                Q_m(m) = gl%Q_ik(j, k)
                m = m + 1
            End do
        End do
        !----------------------------------------------------------------------


        !----------------------------------------------------------------------
        !Calculate the derivative of the group mole fractions Xm with respect
        !to xa. This Derivative is needed for all following derivatives
        Do a = 1, gl%ncomp-1
            Do m = 1, gl%nr_of_groups_mix
                !Three cases have to be distinguished:
                !Case 1: Group m is part of molecule a
                !Case 2: Group m is part of molecule N (Last component)
                !Case 3: Group m is neither part of molecule a nor N

                !Case 1
                if (m >= gl%group_start_index(a) .and. m < gl%group_start_index(a+1)) then
                    !The 0.D0 in the following equation stands for the number of groups of type m in molecule ncomp, which is 0, because here m is a group of molecule a
                    !v_ik(a, m+1-group_start_index(a)) refers to the number of groups of type m in molecule of type a
                    dX_m_dxa(m,a) = ((gl%v_ik(a, m+1-gl%group_start_index(a)) - 0.D0)*sumsum_vnxj &
                        & - (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(a) * gl%v_ik(a, m+1-gl%group_start_index(a))) &
                        &  / sumsum_vnxj**2

                    !Case 2
                elseif (m >= gl%group_start_index(gl%ncomp)) then
                    !The 0.D0 in the following equation stands for the number of groups of type m in molecule a, which is 0, because here m is a group of molecule N
                    !v_ik(ncomp, m+1-group_start_index(ncomp)) refers to the number of groups of type m in molecule of type ncomp
                    dX_m_dxa(m,a) = ((0.D0 - gl%v_ik(gl%ncomp, m+1-gl%group_start_index(gl%ncomp))) * sumsum_vnxj &
                        & - (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(gl%ncomp) * gl%v_ik(gl%ncomp, m+1-gl%group_start_index(gl%ncomp))) &
                        &  / sumsum_vnxj**2

                    !Case 3
                else
                    !The first 0.D0 stands for the number of groups of type m in molecule a and the second 0.D0 stands for the number of groups of type m in molecule N.
                    !Both are 0, because the group belongs to another molecule of type k which needs to be determined
                    !Determine which molecule k the group m belongs to
                    do k=1,gl%ncomp-1
                        if (m >= gl%group_start_index(k) .and. m < gl%group_start_index(k+1)) then
                            exit
                        end if
                    end do
                    dX_m_dxa(m,a) = ((0.D0 - 0.D0) * sumsum_vnxj &
                        & - (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(k) * gl%v_ik(k, m+1-gl%group_start_index(k))) &
                        &  / sumsum_vnxj**2

                end if

            end do
        end do
        !----------------------------------------------------------------------

        !----------------------------------------------------------------------
        !Calculate the derivative of Theta_m with respect to xa
        !This Derivative is needed for all following derivatives
        sum_QndXn_dxa = 0.D0
        do a = 1,gl%ncomp-1
            Do m = 1, gl%nr_of_groups_mix
                sum_QndXn_dxa(a) = sum_QndXn_dxa(a) + Q_m(m) * dX_m_dxa(m,a)
            end do
        end do

        Do a = 1, gl%ncomp-1
            Do m = 1, gl%nr_of_groups_mix
                dTheta_m_dxa(m,a) = (Q_m(m) * dX_m_dxa(m,a) * sum_QnXn - Q_m(m) * X_m(m) * sum_QndXn_dxa(a)) / sum_QnXn**2
            end do
        end do
        !----------------------------------------------------------------------


        !----------------------------------------------------------------------
        !Calculate all group interaction parameters Psi_nm
        !----------------------------------------------------------------------
        Do n = 1, gl%nr_of_groups_mix
            Do m = 1, gl%nr_of_groups_mix
                Psi_nm(n, m) = dexp(-(gl%a_nm(n, m) + gl%b_nm(n, m) * Temp + gl%c_nm(n, m) * Temp**2) / Temp)
            End do
        End do
        !----------------------------------------------------------------------


        !Loop over all molefactions xa
        Do a = 1, gl%ncomp-1
            help_k = 0     !Help index, because for Psi_nm a different index for k is needed (k runs through the groups of each molecule, help_k runs through ALL groups)
            !Loop over all components in the mixture
            Do i = 1, gl%ncomp

                !----------------------------------------------------------------------
                !Calculate dlnGamma_k_dxa
                !Also needed for the following calculations
                !----------------------------------------------------------------------
                Do k = 1, gl%nr_of_groups_i(i)
                    sum_ThetamPsimk = 0.D0
                    sum_dThetam_dxa_Psimk = 0.D0
                    sumsum_ThetamPsikm_dxa = 0.D0
                    help_k = help_k + 1
                    Do m = 1, gl%nr_of_groups_mix
                        sum_ThetamPsimk = sum_ThetamPsimk + Theta_m(m) * Psi_nm(m, help_k)
                        sum_dThetam_dxa_Psimk = sum_dThetam_dxa_Psimk + dTheta_m_dxa(m,a) * Psi_nm(m, help_k)
                        sum_ThetanPsinm = 0.D0
                        sum_dThetan_dxa_Psinm = 0.D0
                        Do n = 1, gl%nr_of_groups_mix
                            sum_ThetanPsinm = sum_ThetanPsinm + Theta_m(n) * Psi_nm(n, m)
                            sum_dThetan_dxa_Psinm = sum_dThetan_dxa_Psinm + dTheta_m_dxa(n,a) * Psi_nm(n, m)
                        End do
                        sumsum_ThetamPsikm_dxa = sumsum_ThetamPsikm_dxa &
                            & + (dTheta_m_dxa(m,a) * Psi_nm(help_k, m) * sum_ThetanPsinm &
                            & - Theta_m(m) * Psi_nm(help_k, m) * sum_dThetan_dxa_Psinm) &
                            & / sum_ThetanPsinm**2
                    End do
                    dlnGamma_k_dxa(k,a) = -gl%Q_ik(i, k) * sum_dThetam_dxa_Psimk / sum_ThetamPsimk &
                        & -gl%Q_ik(i, k) * sumsum_ThetamPsikm_dxa
                End do
                !----------------------------------------------------------------------


                !Calculate the derivative of the logarithm of the residual activity coefficients with respect to xa
                !!!NOTE THAT LNGAMMA_ki OF PURE COMPONENT i IS NOT A FUNCTION OF COMPOSITION, HENCE dlnGamma_ki_dxa(k,i,:) = 0!!! This is indicated with the 0.D0 in the following equation
                Do k = 1, gl%nr_of_groups_i(i)
                    ln_gamma_R_Trho_dxa(1,i,a) = ln_gamma_R_Trho_dxa(1,i,a) + gl%v_ik(i, k) * (dlnGamma_k_dxa(k,a) - 0.D0)
                End do


            End do
        end do

        !Transformation of variables from T,rho,x, to tau,delta, and x
        do a=1,gl%ncomp-1
            Do i = 1, gl%ncomp
                gl%ge%ln_gamma_R_dxa(1,i,a) =  dTred_dxa(a) / tau * ln_gamma_R_Trho(4,i) + ln_gamma_R_Trho_dxa(1,i,a)
            end do
        end do

        !Note in the following that all derivatives with respect to delta at constant tau and x are 0, because the residual part is not a function of density
        if ((GETDER(1) .eq. 1) .or. (GETDER(4) .eq. 1) .or. (GETDER(5) .eq. 1)) then

            if ((gl%ge%GETDER_dxa_prev(1) == 1) .and. (recalc .eqv. .false.)) then

                gl%ge%gE_R_dxa(1,:) = gl%ge%gE_R_dxa_prev(1,:)
                !Transformation of variables from tau,delta,x to T,rho,x
                do a=1,gl%ncomp-1
                    gE_R_Trho_dxa(1,a) = gl%ge%gE_R_dxa(1,a) - dTred_dxa(a) / tau * gE_R_Trho(4)
                end do

            else

                Do a = 1,gl%ncomp-1
                    Do i = 1, gl%ncomp
                        gE_R_Trho_dxa(1,a) = gE_R_Trho_dxa(1,a) + R_const * Temp * x(i) * ln_gamma_R_Trho_dxa(1,i,a)
                    end do
                    gE_R_Trho_dxa(1,a) = gE_R_Trho_dxa(1,a) + R_const * Temp * (ln_gamma_R_Trho(1,a) - ln_gamma_R_Trho(1,gl%ncomp))
                End do
                !Calculate the derivative of gE_C with respect to temperature T at constant rho and x. Also needed for derivatives with respect to tau and xa
                !(Not needed anmyore because called from gE_UNIFAC_MIXDERIVS above)
                !do i=1,ncomp
                !    gE_R_Trho(4) = gE_R_Trho(4) + x(i) * ln_gamma_R_Trho(1,i) * R_const
                !end do
                !Transformation of variables from T,rho,x, to tau,delta, and x
                do a=1,gl%ncomp-1
                    gl%ge%gE_R_dxa(1,a) =  dTred_dxa(a) / tau * gE_R_Trho(4) + gE_R_Trho_dxa(1,a)
                end do

            end if

        end if

        !Second derivative wrt delta and tau
        if (GETDER(2) .eq. 1) then
            ln_gamma_R_Trho_dxa(2,:,:) = 0.D0
            gE_R_Trho_dxa(2,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(2,:,:) = 0.d0
            gl%ge%gE_R_dxa(2,:) = 0.D0
        end if

        !Third derivative wrt delta, tau, and tau
        if (GETDER(3) .eq. 1) then
            ln_gamma_R_Trho_dxa(3,:,:) = 0.D0
            gE_R_Trho_dxa(3,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(3,:,:) = 0.d0
            gl%ge%gE_R_dxa(3,:) = 0.D0
        end if

        !Derivative wrt tau, times tau
        if ((GETDER(4) .eq. 1) .or. (GETDER(5) .eq. 1)) then
            !---------------------------------------------------------------------------------------
            !Calculate the first temperature derivative of all group interaction parameters Psi_nm
            !Do this always because it might be needed for the following derivatives, too
            !---------------------------------------------------------------------------------------
            Do n = 1, gl%nr_of_groups_mix
                Do m = 1, gl%nr_of_groups_mix
                    dPsi_nm_dT(n, m) = (gl%a_nm(n, m) / Temp**2 - gl%c_nm(n, m)) * &
                        & dexp(-(gl%a_nm(n, m) + gl%b_nm(n, m) * Temp + gl%c_nm(n, m) * Temp**2) / Temp)
                End do
            End do

            if ((gl%ge%GETDER_dxa_prev(4) == 1) .and. (recalc .eqv. .false.)) then

                gl%ge%ln_gamma_R_dxa(4,:,:) = gl%ge%ln_gamma_R_dxa_prev(4,:,:)
                gl%ge%gE_R_dxa(4,:) = gl%ge%gE_R_dxa_prev(4,:)

                !Transformation of variables from tau,delta,x to T,rho,x
                do a = 1, gl%ncomp -1
                    do i = 1, gl%ncomp
                        ln_gamma_R_Trho_dxa(4,i,a) = -dTred_dxa(a) / gl%tredmix * ln_gamma_R_Trho(4,i) &
                            & - dTred_dxa(a) / tau * ln_gamma_R_Trho(5,i) &
                            & - tau**2 / gl%tredmix * gl%ge%ln_gamma_R_dxa(4,i,a)
                    end do
                    gE_R_Trho_dxa(4,a) = -dTred_dxa(a) / gl%tredmix * gE_R_Trho(4) &
                        & - dTred_dxa(a) / tau * dTred_dxa(a)*gE_R_Trho(5) &
                        & - tau**2 / gl%tredmix * gl%ge%gE_R_dxa(4,a)
                end do

            else

                !Loop over all molefactions xa
                Do a = 1, gl%ncomp-1
                    help_k = 0     !Help index, because for Psi_nm a different index for k is needed (k runs through the groups of each molecule, help_k runs through ALL groups)
                    !Loop over all components in the mixture
                    Do i = 1, gl%ncomp

                        !----------------------------------------------------------------------
                        !Calculate dlnGamma_k_dxa
                        !Also needed for the following calculations
                        !----------------------------------------------------------------------
                        Do k = 1, gl%nr_of_groups_i(i)
                            sum_ThetamPsimk = 0.D0
                            sum_dThetam_dxa_Psimk = 0.D0
                            sum_Thetam_dPsimk_dT = 0.D0
                            sum_dThetam_dxa_dPsimk_dT = 0.D0
                            sumsum_ThetamPsikm_dxadT = 0.D0
                            help_k = help_k + 1
                            Do m = 1, gl%nr_of_groups_mix
                                sum_ThetamPsimk = sum_ThetamPsimk + Theta_m(m) * Psi_nm(m, help_k)
                                sum_dThetam_dxa_Psimk = sum_dThetam_dxa_Psimk + dTheta_m_dxa(m,a) * Psi_nm(m, help_k)
                                sum_Thetam_dPsimk_dT = sum_Thetam_dPsimk_dT + Theta_m(m) * dPsi_nm_dT(m, help_k)
                                sum_dThetam_dxa_dPsimk_dT = sum_dThetam_dxa_dPsimk_dT + dTheta_m_dxa(m,a) * dPsi_nm_dT(m, help_k)
                                sum_ThetanPsinm = 0.D0
                                sum_dThetan_dxa_Psinm = 0.D0
                                sum_Thetan_dPsinm_dT = 0.D0
                                sum_dThetan_dxa_dPsinm_dT = 0.D0
                                Do n = 1, gl%nr_of_groups_mix
                                    sum_ThetanPsinm = sum_ThetanPsinm + Theta_m(n) * Psi_nm(n, m)
                                    sum_dThetan_dxa_Psinm = sum_dThetan_dxa_Psinm + dTheta_m_dxa(n,a) * Psi_nm(n, m)
                                    sum_Thetan_dPsinm_dT = sum_Thetan_dPsinm_dT + Theta_m(n) * dPsi_nm_dT(n, m)
                                    sum_dThetan_dxa_dPsinm_dT = sum_dThetan_dxa_dPsinm_dT + dTheta_m_dxa(n,a) * dPsi_nm_dT(n, m)
                                End do
                                sumsum_ThetamPsikm_dxadT = sumsum_ThetamPsikm_dxadT &
                                    & + (dTheta_m_dxa(m,a) * dPsi_nm_dT(help_k, m) * sum_ThetanPsinm &
                                    & + Theta_m(m) * dPsi_nm_dT(help_k, m) * sum_dThetan_dxa_Psinm) &
                                    & / sum_ThetanPsinm**2 &
                                    & - (dTheta_m_dxa(m,a) * Psi_nm(help_k, m) * sum_Thetan_dPsinm_dT &
                                    & + Theta_m(m) * Psi_nm(help_k, m) * sum_dThetan_dxa_dPsinm_dT) &
                                    & / sum_ThetanPsinm**2 &
                                    & - 2.D0 * sum_dThetan_dxa_Psinm  &
                                    & * (Theta_m(m) * dPsi_nm_dT(help_k, m) * sum_ThetanPsinm &
                                    & - Theta_m(m) * Psi_nm(help_k, m) * sum_Thetan_dPsinm_dT) &
                                    & / sum_ThetanPsinm**3
                            End do
                            d2lnGamma_k_dxadT(k,a) = gl%Q_ik(i, k) * sum_dThetam_dxa_Psimk * sum_Thetam_dPsimk_dT / sum_ThetamPsimk**2 &
                                & -gl%Q_ik(i, k) * sum_dThetam_dxa_dPsimk_dT / sum_ThetamPsimk &
                                & -gl%Q_ik(i, k) * sumsum_ThetamPsikm_dxadT
                        End do
                        !----------------------------------------------------------------------


                        !Calculate the derivative of the logarithm of the residual activity coefficients with respect to xa and temperature
                        !!!NOTE THAT LNGAMMA_ki OF PURE COMPONENT i IS NOT A FUNCTION OF COMPOSITION, HENCE d2lnGamma_ki_dxadT(k,i,:) = 0!!! This is indicated with the 0.D0 in the following equation
                        Do k = 1, gl%nr_of_groups_i(i)
                            ln_gamma_R_Trho_dxa(4,i,a) = ln_gamma_R_Trho_dxa(4,i,a) + gl%v_ik(i, k) * (d2lnGamma_k_dxadT(k,a) - 0.D0)
                        End do


                    End do
                end do


                !Compute the derivative of gE_R with respect to xa and T (natural variables)
                Do a = 1,gl%ncomp-1
                    Do i = 1, gl%ncomp
                        gE_R_Trho_dxa(4,a) = gE_R_Trho_dxa(4,a) + R_const * x(i) * ln_gamma_R_Trho_dxa(1,i,a) &
                            & +  R_const * Temp * x(i) * ln_gamma_R_Trho_dxa(4,i,a)
                    end do
                    gE_R_Trho_dxa(4,a) = gE_R_Trho_dxa(4,a) + R_const * (ln_gamma_R_Trho(1,a) - ln_gamma_R_Trho(1,gl%ncomp)) &
                        &  + R_const * Temp * (ln_gamma_R_Trho(4,a) - ln_gamma_R_Trho(4,gl%ncomp))
                End do

                !Transformation of variables from T,rho,x, to tau,delta, and x
                do a=1,gl%ncomp-1
                    Do i = 1, gl%ncomp
                        gl%ge%ln_gamma_R_dxa(4,i,a) = - dTred_dxa(a) / tau**2 * ln_gamma_R_Trho(4,i) &
                            & - gl%tredmix / tau**3 * dTred_dxa(a) * ln_gamma_R_Trho(5,i) &
                            & - gl%tredmix / tau**2 * ln_gamma_R_Trho_dxa(4,i,a)
                    end do
                    gl%ge%gE_R_dxa(4,a) = - dTred_dxa(a) / tau**2 * gE_R_Trho(4) &
                        & - gl%tredmix / tau**3 * dTred_dxa(a) * gE_R_Trho(5) &
                        & - gl%tredmix / tau**2 * gE_R_Trho_dxa(4,a)
                end do

            end if

        end if


        !Second derivative wrt tau NOT YET IMPLEMENTED
        if (GETDER(5) .eq. 1) then
            ln_gamma_R_Trho_dxa(5,:,:) = 0.D0
            gE_R_Trho_dxa(5,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(5,:,:) = 0.d0
            gl%ge%gE_R_dxa(5,:) = 0.D0
        end if
        !Second derivative wrt delta and tau
        if (GETDER(6) .eq. 1) then
            ln_gamma_R_Trho_dxa(6,:,:) = 0.D0
            gE_R_Trho_dxa(6,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(6,:,:) = 0.d0
            gl%ge%gE_R_dxa(6,:) = 0.D0
        end if
        !Third derivative wrt delta, tau, and tau
        if (GETDER(7) .eq. 1) then
            ln_gamma_R_Trho_dxa(7,:,:) = 0.D0
            gE_R_Trho_dxa(7,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(7,:,:) = 0.d0
            gl%ge%gE_R_dxa(7,:) = 0.D0
        end if
        !Third derivative wrt delta
        if (GETDER(8) .eq. 1) then
            ln_gamma_R_Trho_dxa(8,:,:) = 0.D0
            gE_R_Trho_dxa(8,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(8,:,:) = 0.d0
            gl%ge%gE_R_dxa(8,:) = 0.D0
        end if
        !Third derivative wrt tau
        if (GETDER(9) .eq. 1) then
            ln_gamma_R_Trho_dxa(9,:,:) = 0.D0
            gE_R_Trho_dxa(9,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(9,:,:) = 0.d0
            gl%ge%gE_R_dxa(9,:) = 0.D0
        end if
        !Third derivative wrt tau, delta, and delta
        if (GETDER(10) .eq. 1) then
            ln_gamma_R_Trho_dxa(10,:,:) = 0.D0
            gE_R_Trho_dxa(10,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(10,:,:) = 0.d0
            gl%ge%gE_R_dxa(10,:) = 0.D0
        end if
        !Fourth derivative wrt delta
        if (GETDER(11) .eq. 1) then
            ln_gamma_R_Trho_dxa(11,:,:) = 0.D0
            gE_R_Trho_dxa(11,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(11,:,:) = 0.d0
            gl%ge%gE_R_dxa(11,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, delta and tau
        if (GETDER(12) .eq. 1) then
            ln_gamma_R_Trho_dxa(12,:,:) = 0.D0
            gE_R_Trho_dxa(12,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(12,:,:) = 0.d0
            gl%ge%gE_R_dxa(12,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, tau, tau
        if (GETDER(13) .eq. 1) then
            ln_gamma_R_Trho_dxa(13,:,:) = 0.D0
            gE_R_Trho_dxa(13,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(13,:,:) = 0.d0
            gl%ge%gE_R_dxa(13,:) = 0.D0
        end if
        !Fourth derivative wrt delta, tau, tau, tau
        if (GETDER(14) .eq. 1) then
            ln_gamma_R_Trho_dxa(14,:,:) = 0.D0
            gE_R_Trho_dxa(14,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(14,:,:) = 0.d0
            gl%ge%gE_R_dxa(14,:) = 0.D0
        end if
        !Fourth derivative wrt tau
        if (GETDER(15) .eq. 1) then
            ln_gamma_R_Trho_dxa(15,:,:) = 0.D0
            gE_R_Trho_dxa(15,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(15,:,:) = 0.d0
            gl%ge%gE_R_dxa(15,:) = 0.D0
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    End if


    !Write new values in save variables
    !gl%ge%Temp_dxa_prev = Temp
    !gl%ge%GETDER_dxa_prev = GETDER
    !gl%ge%molfrac_dxa_prev = gl%molfractions
    !gl%ge%Eq_type_dxa_prev = gl%Eq_type
    !gl%ge%mixtype_dxa_prev = gl%mix_type
    !gl%ge%C_or_R_dxa_prev = C_or_R
    !gl%ge%gE_C_dxa_prev = gl%ge%gE_C_dxa
    !gl%ge%gE_R_dxa_prev = gl%ge%gE_R_dxa
    !gl%ge%ln_gamma_C_dxa_prev = gl%ge%ln_gamma_C_dxa
    !gl%ge%ln_gamma_R_dxa_prev = gl%ge%ln_gamma_R_dxa


    end subroutine
    !--------------------------------------------------------------------------------------




    !Derivatives of UNIFAC with respect to mole frations xa and xb
    !Needed for phase equilibrium routines for PSRK(gl,and VTPR and new Helmholtz+gE) model
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine gE_UNIFAC_MIXDERIVS_dxadxb(gl,Temp, GETDER, C_or_R, errval)
    !------------------------------------------------------------------------------------------------------------------------------
    !Subroutine to calculate the excess Gibbs energy and the activity coefficients of the
    !components in the mixture according to the UNIFAC model as well as tau- and delta-
    !derivatives of the model.
    !All quantities are given as derivatives with respect to mole fractions xa and xb for all components
    !in the mixture. The convention that the "last" mole fraction xN of component N is replaced with
    !1-sum(xk) is used here.
    !
    !INPUT
    !Temp       -    Temperature in K
    !Note that the model makes use of the global variable "molfractions" defined in the module
    !"module_fluid_parameters"
    ! GETDER      - AN ARRAY WITH nderivs ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED. !!!ALL DERIVATIVES ARE ADDITIONALLY DERIVATIVES WITH RESPECT TO xa and xb!!!:
    !                1. gE from UNIFAC as a function of temperature and composition
    !                2. 1ST DERIVATIVE OF gE WITH RESPECT TO DEL AT CONSTANT TAU
    !                3. 2ND DERIVATIVE OF gE WITH RESPECT TO DEL AT CONSTANT TAU
    !                4. 1ST DERIVATIVE OF gE WITH RESPECT TO TAU AT CONSTANT DEL
    !                5: 2ND DERIVATIVE OF gE WITH RESPECT TO TAU AT CONSTANT DEL
    !                6: 2ND MIXED DERIVATIVE OF gE WITH RESPECT TO DEL AND TAU
    !                7: 3RD MIXED DERIVATIVE OF gE WITH RESPECT TO DEL, TAU, AND TAU
    !                8: 3RD DERIVATIVE OF gE WITH RESPECT TO DEL
    !                9: 3RD DERIVATIVE OF gE WITH RESPECT TO TAU
    !               10: 3RD MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, DEL, AND DEL
    !               11: 4TH DERIVATIVE OF gE WITH RESPECT TO DEL
    !               12: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, DEL, DEL AND DEL
    !               13: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, TAU, DEL, AND DEL
    !               14: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, TAU, TAU AND DEL
    !               15: 4TH DERIVATIVE OF gE WITH RESPECT TO TAU
    !C_or_R       - Integer that specifies if the combinatorial, the residual or both parts shall be calculated
    !               C_or_R = 0: Calculate both parts
    !               C_or_R = 1: Calculate the combinatorial part only
    !               C_or_R = 2: Calculate the the residual part only
    !
    !OUTPUT
    !gE_C_dxadxb         -    Matrix with molar excess Gibbs energy of the combinatorial part (C) in J/molK and derivatives with respect to tau and delta and xa
    !gE_R_dxadxb         -    Matrix with molar excess Gibbs energy of the residual part (R) in J/molK and derivatives with respect to tau and delta and xa
    !ln_gamma_C_dxadxb   -    Matrix containing the natural logarithm of the activity coefficients of the combinatorial part for each component and derivatives with respect to tau and delta and xa
    !ln_gamma_R_dxadxb   -    Matrix containing the natural logarithm of the activity coefficients of the residual part for each component and derivatives with respect to tau and delta and xa
    !errval              -    Indicates if an error occurred during calculations
    !
    !Andreas Jäger, April 2017
    !------------------------------------------------------------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl


    double precision :: Temp
    integer, dimension(nderivs):: GETDER                                    !array specifier to indicate, which derivative is needed
    !double precision, allocatable :: gl%ge%gE_C_dxadxb(:,:,:)                 !Combinatorial part of gE and derivatives with respect to tau and delta and xa and xb
    !double precision, allocatable :: gl%ge%gE_R_dxadxb(:,:,:)                 !Residual part of gE and derivatives with respect to tau and delta and xa and xb
    !double precision, allocatable :: gl%ge%ln_gamma_C_dxadxb(:,:,:,:)     !Combinatorial activity coefficients and derivatives with respect to delta and tau and xa and xb
    !double precision, allocatable :: gl%ge%ln_gamma_R_dxadxb(:,:,:,:)     !Residual activity coefficients and derivatives with respect to delta and tau and xa and xb
    integer:: C_or_R
    integer:: errval

    double precision, dimension(30) :: r_i              !Molecular vdW volume of components i
    double precision, dimension(30) :: q_i              !Molecular vdW surface area of components i

    double precision :: R_const                         !universal gas constant

    !Help variables
    double precision, dimension(30) :: x
    double precision, dimension(30) :: phi_i
    double precision, dimension(30) :: Theta_i
    !double precision, dimension(100) :: lnGamma_k
    !double precision, dimension(100,30) :: lnGamma_ki   !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    double precision, dimension(100) :: Theta_m
    double precision, dimension(100,30) :: Theta_mi     !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    double precision, dimension(100) :: X_m             !Mole fraction of groups in the mixture
    double precision, dimension(100,30) :: X_mi         !Mole fraction of groups in pure component i. !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    double precision, dimension(:,:), allocatable :: Psi_nm
    double precision, dimension(100):: Q_m              !Help variable to list all group areas in one vector in order to not have to map them back to Q_ik(i, k) for derivative calculation

    !More help variables
    double precision:: sum_rjxj
    double precision:: sum_qjxj
    double precision:: sum_QnXn
    double precision, dimension(30):: sum_QndXn_dxa
    double precision, dimension(30,30):: sum_QndXn_dxadxb
    double precision:: sumsum_vnxj
    !double precision:: sum_xjLj
    double precision:: sum_ThetamPsimk
    !double precision:: sumsum_ThetamPsikm
    double precision:: sum_ThetanPsinm

    !Derivatives of the combinatorial part with respect to "natural" variables
    double precision, dimension(nderivs, 30) :: ln_gamma_C_Trho       !Combinatorial activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs, 30) :: ln_gamma_R_Trho       !Residual activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs) :: gE_C_Trho        !Combinatorial part of gE and derivatives with respect to T and rho
    double precision, dimension(nderivs) :: gE_R_Trho        !Residual part of gE and derivatives with respect to T and rho
    double precision, allocatable  :: ln_gamma_C_Trho_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to rho and T, and xa
    double precision, allocatable  :: ln_gamma_R_Trho_dxa(:,:,:)       !Residual activity coefficients and derivatives with respect to rho and T, and xa
    double precision, dimension(nderivs, 30) :: gE_C_Trho_dxa        !Combinatorial part of gE and derivatives with respect to T and rho and xa
    double precision, dimension(nderivs, 30) :: gE_R_Trho_dxa        !Residual part of gE and derivatives with respect to T and rho and xa
    double precision, allocatable :: gE_C_Trho_dxadxb(:,:,:)                 !Combinatorial part of gE and derivatives with respect to T and rho and xa and xb
    double precision, allocatable :: gE_R_Trho_dxadxb(:,:,:)                 !Residual part of gE and derivatives with respect to T and rho and xa and xb
    double precision, allocatable :: ln_gamma_C_Trho_dxadxb(:,:,:,:)       !Combinatorial activity coefficients and derivatives with respect to rho and T and xa and xb
    double precision, allocatable :: ln_gamma_R_Trho_dxadxb(:,:,:,:)       !Residual activity coefficients and derivatives with respect to rho and T and xa and xb

    !Derivatives of the residual part with respect to "natural" variables
    !For combinatorial part
    double precision, dimension(30,30):: dphi_i_dxa
    double precision, dimension(30,30):: dTheta_i_dxa
    double precision, allocatable:: d2phi_i_dxadxb(:,:,:)
    double precision, allocatable:: d2Theta_i_dxadxb(:,:,:)
    !For residual part
    !double precision, dimension(100,30):: dlnGamma_k_dxa
    !double precision, dimension(100,30):: d2lnGamma_k_dxadT
    double precision, allocatable:: d2lnGamma_k_dxadxb(:,:,:)
    !double precision, dimension(100,100) :: dPsi_nm_dT
    double precision, dimension(100,30):: dTheta_m_dxa
    double precision, allocatable:: d2Theta_m_dxadxb(:,:,:)
    double precision, dimension(100,30):: dX_m_dxa
    double precision, allocatable:: d2X_m_dxadxb(:,:,:)
    double precision:: sum_dThetam_dxa_Psimk
    double precision:: sum_dThetam_dxb_Psimk
    double precision:: sum_d2Thetam_dxadxb_Psimk
    double precision:: sum_dThetan_dxa_Psinm
    double precision:: sum_dThetan_dxb_Psinm
    double precision:: sum_d2Thetan_dxadxb_Psinm
    double precision:: sumsum_ThetamPsikm_dxadxb          !Variable for double sum in xa and xb derivative
    !double precision:: sum_Thetam_dPsimk_dT
    !double precision:: sum_Thetan_dPsinm_dT
    !double precision:: sum_dThetam_dxa_dPsimk_dT
    !double precision:: sum_dThetan_dxa_dPsinm_dT
    !double precision:: sumsum_ThetamPsikm_dxadT        !Variable for double sum in xa and temperature derivative

    double precision, dimension(30):: dTred_dxa, drhored_dxa            !Variables for the derivatives of the reducing functions with respect to molfractions xa
    double precision, dimension(30):: d2Tred_dxa2, drhored_dxa2         !Variables for the second derivatives of the reducing functions with respect to molfractions xa
    double precision, dimension(30,30):: d2Tred_dxadxb, drhored_dxadxb  !Variables for the second derivatives of the reducing functions with respect to molfractions xa and xb

    !Required derivatives of activities and excess Gibbs energy only with respect to tau and delta (not xa and xb, however, needed for the calculation of derivatives with respect to xa)
    integer, dimension(nderivs):: GETDER_no_xaxb        !array specifier to indicate, which derivative is needed
    !double precision, dimension(nderivs) :: gl%ge%gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gl%ge%gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau

    !Required derivatives of activities and excess Gibbs energy only with respect to tau and delta and xa (not xb, however, needed for the calculation of derivatives with respect to xa and xb)
    integer, dimension(nderivs):: GETDER_no_xb                      !array specifier to indicate, which derivative is needed
    !double precision, dimension(nderivs,30) :: gl%ge%gE_C_dxa                !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs,30) :: gl%ge%gE_R_dxa                !Residual part of gE and derivatives with respect to tau and delta
    !double precision, allocatable:: gl%ge%ln_gamma_C_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, allocatable:: gl%ge%ln_gamma_R_dxa(:,:,:)      !Residual activity coefficients and derivatives with respect to delta and tau


    !Corresponding states relevant variables
    double precision:: tau, del

    !Summation variables
    integer:: a, b, i, j, k, m, n, help_k

    !Variable that indicates if derivatives need to be recalculated
    logical:: recalc

    !allocation of large arrays
    if(.not. allocated(gl%ge%gE_C_dxadxb)) allocate(gl%ge%gE_C_dxadxb(nderivs, gl%ncomp, gl%ncomp))
    if(.not. allocated(gl%ge%gE_R_dxadxb)) allocate(gl%ge%gE_R_dxadxb(nderivs, gl%ncomp, gl%ncomp))
    if(.not. allocated(gl%ge%ln_gamma_C_dxadxb)) allocate(gl%ge%ln_gamma_C_dxadxb(nderivs,gl%ncomp, gl%ncomp, gl%ncomp))
    if(.not. allocated(gl%ge%ln_gamma_R_dxadxb)) allocate(gl%ge%ln_gamma_R_dxadxb(nderivs,gl%ncomp, gl%ncomp, gl%ncomp))
    if(.not. allocated(ln_gamma_C_Trho_dxa)) allocate(ln_gamma_C_Trho_dxa(nderivs, gl%ncomp, gl%ncomp))
    if(.not. allocated(ln_gamma_R_Trho_dxa)) allocate(ln_gamma_R_Trho_dxa(nderivs, gl%ncomp, gl%ncomp))
    if(.not. allocated(gE_C_Trho_dxadxb)) allocate(gE_C_Trho_dxadxb(nderivs, gl%ncomp, gl%ncomp))
    if(.not. allocated(gE_R_Trho_dxadxb)) allocate(gE_R_Trho_dxadxb(nderivs, gl%ncomp, gl%ncomp))
    if(.not. allocated(ln_gamma_C_Trho_dxadxb)) allocate(ln_gamma_C_Trho_dxadxb(nderivs,gl%ncomp, gl%ncomp, gl%ncomp))
    if(.not. allocated(ln_gamma_R_Trho_dxadxb)) allocate(ln_gamma_R_Trho_dxadxb(nderivs,gl%ncomp, gl%ncomp, gl%ncomp))
    if(.not. allocated(d2phi_i_dxadxb)) allocate(d2phi_i_dxadxb(gl%ncomp,gl%ncomp,gl%ncomp))
    if(.not. allocated(d2Theta_i_dxadxb)) allocate(d2Theta_i_dxadxb(gl%ncomp,gl%ncomp,gl%ncomp))

    if(.not. allocated(d2lnGamma_k_dxadxb)) allocate(d2lnGamma_k_dxadxb(gl%nr_of_groups_mix,gl%ncomp,gl%ncomp))
    if(.not. allocated(d2Theta_m_dxadxb)) allocate(d2Theta_m_dxadxb(gl%nr_of_groups_mix,gl%ncomp,gl%ncomp))
    if(.not. allocated(d2X_m_dxadxb)) allocate(d2X_m_dxadxb(gl%nr_of_groups_mix,gl%ncomp,gl%ncomp))

    if(.not. allocated(gl%ge%ln_gamma_C_dxa)) allocate(gl%ge%ln_gamma_C_dxa(nderivs, gl%ncomp, gl%ncomp))
    if(.not. allocated(gl%ge%ln_gamma_R_dxa)) allocate(gl%ge%ln_gamma_R_dxa(nderivs, gl%ncomp, gl%ncomp))

    if(.not. allocated(Psi_nm)) allocate(Psi_nm(gl%nr_of_groups_mix*gl%ncomp, gl%nr_of_groups_mix*gl%ncomp))
    !end of allocation


    errval = 0

    !Check for wrong input to C_or_R
    if ((C_or_R < 0) .or. (C_or_R > 2)) then
        errval = -1111
        gl%ge%gE_C_dxadxb = errval
        gl%ge%gE_R_dxadxb = errval
        gl%ge%ln_gamma_C_dxadxb = errval
        gl%ge%ln_gamma_R_dxadxb = errval
        return
    end if

    recalc = .true.

    !Check whether calculations have to be redone or not.
    !If the temperature, (density), composition, equation of state and, the required part of UNIFAC stayed the same, no need to recalculate
    if(dabs(gl%ge%Temp_dxadxb_prev - Temp) > 1.D-16) then
        recalc = .true.
    end if
    do i=1,gl%ncomp
        if(dabs(gl%ge%molfrac_dxadxb_prev(i) - gl%molfractions(i)) > 1.D-16) then
            recalc = .true.
        end if
        if(gl%ge%Eq_type_dxadxb_prev(i) .ne. gl%Eq_type(i)) then
            recalc = .true.
        end if
    end do
    if (gl%ge%mixtype_dxadxb_prev .ne. gl%mix_type) then
        recalc = .true.
    end if
    if (gl%ge%C_or_R_dxadxb_prev .ne. C_or_R) then
        recalc = .true.
    end if


    !Initialize variables
    ln_gamma_C_Trho = 0.D0
    ln_gamma_R_Trho = 0.D0
    gE_C_Trho = 0.D0
    gE_R_Trho = 0.d0
    ln_gamma_C_Trho_dxa = 0.D0
    ln_gamma_R_Trho_dxa = 0.d0
    gE_C_Trho_dxa = 0.D0
    gE_R_Trho_dxa = 0.D0
    gl%ge%ln_gamma_C_dxa = 0.D0
    gl%ge%ln_gamma_R_dxa = 0.D0
    gl%ge%gE_C_dxa = 0.D0
    gl%ge%gE_R_dxa = 0.D0
    gl%ge%gE_C_dxadxb = 0.D0
    gl%ge%gE_R_dxadxb = 0.D0
    gl%ge%ln_gamma_C_dxadxb = 0.D0
    gl%ge%ln_gamma_R_dxadxb = 0.D0
    gE_C_Trho_dxadxb = 0.D0
    gE_R_Trho_dxadxb = 0.D0
    ln_gamma_C_Trho_dxadxb = 0.D0
    ln_gamma_R_Trho_dxadxb = 0.D0

    !Value for the ideal gas constant
    R_const = 8.3144598D0

    !Get mole fractions of the mixture from module variable
    x = gl%molfractions

    !Calculate tau
    tau = gl%tredmix / Temp

    !Set dummy value for delta (Not needed because UNIFAC is not a function of density. However, already implemented here for possible later modifications)
    del = 1.D0

    !Derivatives of tredmix with respect to xa are required for the following derivatives
    call dYr_dxi(gl,dTred_dxa, drhored_dxa)

    !Derivatives of tredmix with respect to xa and xb are required for the following derivatives
    call d2Yr_dxi2(gl,d2Tred_dxa2, drhored_dxa2)
    call d2Yr_dxidxj(gl,d2Tred_dxadxb, drhored_dxadxb)
    !write the derivatives with respect to two times the same mole fraction also in the variable d2Tred_dxadxb
    do i = 1, gl%ncomp-1
        d2Tred_dxadxb(i,i) = d2Tred_dxa2(i)
    end do

    !Get the required derivatives of activities and excess Gibbs energies with respect tau (and not with respect to xa and xb, however, needed for the derivatives here)
    GETDER_no_xaxb = 0
    if (GETDER(1) .eq. 1)  then
        GETDER_no_xaxb(1) = 1
        GETDER_no_xaxb(4) = 1
        GETDER_no_xaxb(5) = 1
    end if
    !call gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER_no_xaxb, gl%ge%gE_C, gl%ge%gE_R, gl%ge%ln_gamma_C, gl%ge%ln_gamma_R, C_or_R, errval)
    call gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER_no_xaxb, C_or_R, errval)
    if (GETDER(1) .eq. 1) then
        !Transform the tau derivatives into the needed temperature derivatives
        ln_gamma_C_Trho(1,:) = gl%ge%ln_gamma_C(1,:)
        gE_C_Trho(1) = gl%ge%gE_C(1)
        ln_gamma_R_Trho(1,:) = gl%ge%ln_gamma_R(1,:)
        gE_R_Trho(1) = gl%ge%gE_R(1)
        ln_gamma_C_Trho(4,:) = - tau**2 / gl%tredmix * gl%ge%ln_gamma_C(4,:)
        gE_C_Trho(4) = - tau**2 / gl%tredmix * gl%ge%gE_C(4)
        ln_gamma_R_Trho(4,:) = - tau**2 / gl%tredmix * gl%ge%ln_gamma_R(4,:)
        gE_R_Trho(4) = - tau**2 / gl%tredmix * gl%ge%gE_R(4)
        ln_gamma_C_Trho(5,:) = tau**4 / gl%tredmix**2 * gl%ge%ln_gamma_C(5,:) - 2.D0 * tau / gl%tredmix * ln_gamma_C_Trho(4,:)
        gE_C_Trho(5) = tau**4 / gl%tredmix**2 * gl%ge%gE_C(5) - 2.D0 * tau / gl%tredmix * gE_C_Trho(4)
        ln_gamma_R_Trho(5,:) = tau**4 / gl%tredmix**2 * gl%ge%ln_gamma_R(5,:) - 2.D0 * tau / gl%tredmix * ln_gamma_R_Trho(4,:)
        gE_R_Trho(5) = tau**4 / gl%tredmix**2 * gl%ge%gE_R(5) - 2.D0 * tau / gl%tredmix * gE_R_Trho(4)
    end if



    !Get the required derivatives of activities and excess Gibbs energies with respect tau (and not with respect to xa and xb, however, needed for the derivatives here)
    GETDER_no_xb = 0
    if (GETDER(1) .eq. 1)  then
        GETDER_no_xb(1) = 0
        GETDER_no_xb(4) = 1
        !GETDER_no_xb(5) = 1
    end if
    call gE_UNIFAC_MIXDERIVS_dxa(gl,Temp, GETDER_no_xb, C_or_R, errval)
    if (GETDER(1) .eq. 1) then
        !Transform the tau derivatives into the needed temperature derivatives
        do a = 1, gl%ncomp-1
            do i = 1, gl%ncomp
                ln_gamma_C_Trho_dxa(1,i,a) =  gl%ge%ln_gamma_C_dxa(1,i,a) - dTred_dxa(a) / tau * ln_gamma_C_Trho(4,i)
                ln_gamma_R_Trho_dxa(1,i,a) =  gl%ge%ln_gamma_R_dxa(1,i,a) - dTred_dxa(a) / tau * ln_gamma_R_Trho(4,i)
                ln_gamma_C_Trho_dxa(4,i,a) = - tau**2 / gl%tredmix * gl%ge%ln_gamma_C_dxa(4,i,a) &
                    & - dTred_dxa(a) / gl%tredmix * ln_gamma_C_Trho(4,i) &
                    & - dTred_dxa(a) / tau * ln_gamma_C_Trho(5,i)
                ln_gamma_R_Trho_dxa(4,i,a) = - tau**2 / gl%tredmix * gl%ge%ln_gamma_R_dxa(4,i,a) &
                    & - dTred_dxa(a) / gl%tredmix * ln_gamma_R_Trho(4,i) &
                    & - dTred_dxa(a) / tau * ln_gamma_R_Trho(5,i)
            end do
            gE_C_Trho_dxa(1,a) = gl%ge%gE_C_dxa(1,a) - dTred_dxa(a) / tau * gE_C_Trho(4)
            gE_R_Trho_dxa(1,a) = gl%ge%gE_R_dxa(1,a) - dTred_dxa(a) / tau * gE_R_Trho(4)
            gE_C_Trho_dxa(4,a) = - tau**2 / gl%tredmix * gl%ge%gE_C_dxa(4,a) &
                & - dTred_dxa(a) / gl%tredmix * gE_C_Trho(4) &
                & - dTred_dxa(a) / tau * gE_C_Trho(5)
            gE_R_Trho_dxa(4,a) = - tau**2 / gl%tredmix * gl%ge%gE_R_dxa(4,a) &
                & - dTred_dxa(a) / gl%tredmix * gE_R_Trho(4) &
                & - dTred_dxa(a) / tau * gE_R_Trho(5)
        end do
    end if



    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 1)) then
        !Calculate the combinatorial part of the activity coefficients
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !ln_gamma_C is always needed for the following calculations, thus it is always calculated

        if ((gl%ge%GETDER_dxadxb_prev(1) == 1) .and. (recalc .eqv. .false.)) then

            gl%ge%ln_gamma_C_dxadxb(1,:,:,:) = gl%ge%ln_gamma_C_dxadxb_prev(1,:,:,:)

            !Transformation of variables from tau,delta, and x to T,rho,x
            do b = 1, gl%ncomp-1
                do a = 1, gl%ncomp-1
                    do i = 1, gl%ncomp
                        ln_gamma_C_Trho_dxadxb(1,i,a,b) = gl%ge%ln_gamma_C_dxadxb(1,i,a,b) &
                            & - dTred_dxa(a) * dTred_dxa(b) / tau**2 * ln_gamma_C_Trho(5,i) &
                            & - dTred_dxa(a) / tau * ln_gamma_C_Trho_dxa(4,i,b) &
                            & - dTred_dxa(b) / tau * ln_gamma_C_Trho_dxa(4,i,a) &
                            & - d2Tred_dxadxb(a,b) / tau * ln_gamma_C_Trho(4,i)
                    end do
                end do
            end do

        else

            !Calculate help variables
            r_i = 0.D0
            q_i = 0.D0
            Do i = 1, gl%ncomp
                Do k = 1, gl%nr_of_groups_i(i)
                    r_i(i) = r_i(i) + gl%v_ik(i, k) * gl%R_ik(i, k)
                    q_i(i) = q_i(i) + gl%v_ik(i, k) * gl%Q_ik(i, k)
                End do
            End do

            sum_rjxj = 0.D0
            sum_qjxj = 0.D0
            Do j = 1, gl%ncomp
                sum_rjxj = sum_rjxj + r_i(j) * x(j)
                sum_qjxj = sum_qjxj + q_i(j) * x(j)
            End do
            Do i = 1, gl%ncomp
                phi_i(i) = r_i(i) / sum_rjxj
                Theta_i(i) = q_i(i) / sum_qjxj
            End do

            !Calculate the derivativs of phi_i and Theta_i with respect to xa. These derivatives are also needed for all tau and delta derivatives
            Do a = 1, gl%ncomp-1
                Do i = 1, gl%ncomp
                    dphi_i_dxa(i,a) = -r_i(i) * (r_i(a) - r_i(gl%ncomp)) / sum_rjxj**2
                    dTheta_i_dxa(i,a) = -q_i(i) * (q_i(a) - q_i(gl%ncomp)) / sum_qjxj**2
                End do
            End do

            !Calculate the derivativs of phi_i and Theta_i with respect to xa and xb.
            Do b = 1, gl%ncomp -1
                Do a = 1, gl%ncomp-1
                    Do i = 1, gl%ncomp
                        d2phi_i_dxadxb(i,a,b) = 2.D0 * r_i(i) * (r_i(a) - r_i(gl%ncomp)) * (r_i(b) - r_i(gl%ncomp)) / sum_rjxj**3
                        d2Theta_i_dxadxb(i,a,b) = 2.D0 * q_i(i) * (q_i(a) - q_i(gl%ncomp)) * (q_i(b) - q_i(gl%ncomp)) / sum_qjxj**3
                    End do
                End do
            end do

            !Compute the mole fraction derivatives of all activity coefficients of the combinatorial part
            do b = 1, gl%ncomp-1
                do a = 1, gl%ncomp-1
                    Do i = 1, gl%ncomp
                        ln_gamma_C_Trho_dxadxb(1,i,a,b) = (1.D0 / phi_i(i) - 1.D0) * d2phi_i_dxadxb(i,a,b) &
                            & - dphi_i_dxa(i,a) * dphi_i_dxa(i,b) / phi_i(i)**2 &
                            & - 5.D0 * q_i(i) &
                            & * ((Theta_i(i) / phi_i(i) - 1.D0) * (d2phi_i_dxadxb(i,a,b) / Theta_i(i)  &
                            & - phi_i(i) * d2Theta_i_dxadxb(i,a,b) / Theta_i(i)**2 &
                            & - dTheta_i_dxa(i,a) * dphi_i_dxa(i,b) / Theta_i(i)**2 &
                            & - dTheta_i_dxa(i,b) * dphi_i_dxa(i,a) / Theta_i(i)**2 &
                            & + 2.D0 * phi_i(i) * dTheta_i_dxa(i,a) * dTheta_i_dxa(i,b) / Theta_i(i)**3) &
                            & + (Theta_i(i) * dphi_i_dxa(i,a) - phi_i(i) * dTheta_i_dxa(i,a)) &
                            & * (phi_i(i) * dTheta_i_dxa(i,b) - Theta_i(i) * dphi_i_dxa(i,b)) &
                            & / Theta_i(i)**2 / phi_i(i)**2)
                    End do
                End do
            end do

            !Transformation of variables from T,rho,x to tau,delta,x
            do b = 1, gl%ncomp-1
                do a = 1, gl%ncomp-1
                    do i = 1, gl%ncomp
                        gl%ge%ln_gamma_C_dxadxb(1,i,a,b) = ln_gamma_C_Trho_dxadxb(1,i,a,b) &
                            & + dTred_dxa(a) * dTred_dxa(b) / tau**2 * ln_gamma_C_Trho(5,i) &
                            & + dTred_dxa(a) / tau * ln_gamma_C_Trho_dxa(4,i,b) &
                            & + dTred_dxa(b) / tau * ln_gamma_C_Trho_dxa(4,i,a) &
                            & + d2Tred_dxadxb(a,b) / tau * ln_gamma_C_Trho(4,i)
                    end do
                end do
            end do

        end if

        if (GETDER(1) .eq. 1) then

            if ((gl%ge%GETDER_dxadxb_prev(1) == 1) .and. (recalc .eqv. .false.)) then

                gl%ge%gE_C_dxadxb(1,:,:) = gl%ge%gE_C_dxadxb_prev(1,:,:)

                !Transformation of variables from tau,delta, and x to T,rho,x
                do b = 1, gl%ncomp-1
                    do a = 1, gl%ncomp-1
                        gE_C_Trho_dxadxb(1,a,b) = gl%ge%gE_C_dxadxb(1,a,b) &
                            & - dTred_dxa(a) * dTred_dxa(b) / tau**2 * gE_C_Trho(5) &
                            & - dTred_dxa(a) / tau * gE_C_Trho_dxa(4,b) &
                            & - dTred_dxa(b) / tau * gE_C_Trho_dxa(4,a) &
                            & - d2Tred_dxadxb(a,b) / tau * gE_C_Trho(4)
                    end do
                end do

            else

                Do b = 1, gl%ncomp-1
                    Do a = 1,gl%ncomp-1
                        Do i = 1, gl%ncomp
                            gE_C_Trho_dxadxb(1,a,b) = gE_C_Trho_dxadxb(1,a,b) + R_const * Temp * x(i) * ln_gamma_C_Trho_dxadxb(1,i,a,b)
                        end do
                        gE_C_Trho_dxadxb(1,a,b) = gE_C_Trho_dxadxb(1,a,b) &
                            & + R_const * Temp * (ln_gamma_C_Trho_dxa(1,a,b) - ln_gamma_C_Trho_dxa(1,gl%ncomp,b)) &
                            & + R_const * Temp * (ln_gamma_C_Trho_dxa(1,b,a) - ln_gamma_C_Trho_dxa(1,gl%ncomp,a))
                    End do
                End do

                !Transformation of variables from T,rho,x to tau,delta,x
                do b = 1, gl%ncomp-1
                    do a = 1, gl%ncomp-1
                        gl%ge%gE_C_dxadxb(1,a,b) = gE_C_Trho_dxadxb(1,a,b) &
                            & + dTred_dxa(a) * dTred_dxa(b) / tau**2 * gE_C_Trho(5) &
                            & + dTred_dxa(a) / tau * gE_C_Trho_dxa(4,b) &
                            & + dTred_dxa(b) / tau * gE_C_Trho_dxa(4,a) &
                            & + d2Tred_dxadxb(a,b) / tau * gE_C_Trho(4)
                    end do
                end do

            end if
        end if

        !All other derivatives of the combinatorial part are 0 (except for the tau-only derivatives of gE), because the combinatorial part is a function of composition only
        !Derivative wrt delta
        if (GETDER(2) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(2,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(2,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(2,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(2,:,:) = 0.D0
        end if
        !Second derivative wrt delta
        if (GETDER(3) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(3,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(3,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(3,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(3,:,:) = 0.D0
        end if

        !Derivative wrt tau
        if (GETDER(4) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(4,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(4,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(4,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(4,:,:) = 0.D0
        end if

        !Second derivative wrt tau
        if (GETDER(5) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(5,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(5,:,:) = 0.D0
            gl%ge%ln_gamma_C_dxadxb(5,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(5,:,:) = 0.D0
        end if
        !Second derivative wrt delta and tau
        if (GETDER(6) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(6,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(6,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(6,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(6,:,:) = 0.D0
        end if
        !Third derivative wrt delta, tau, and tau
        if (GETDER(7) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(7,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(7,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(7,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(7,:,:) = 0.D0
        end if
        !Third derivative wrt delta
        if (GETDER(8) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(8,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(8,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(8,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(8,:,:) = 0.D0
        end if
        !Third derivative wrt tau
        if (GETDER(9) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(9,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(9,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(9,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(9,:,:) = 0.D0
        end if
        !Third derivative wrt tau, delta, and delta
        if (GETDER(10) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(10,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(10,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(10,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(10,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta
        if (GETDER(11) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(11,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(11,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(11,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(11,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, delta and tau
        if (GETDER(12) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(12,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(12,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(12,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(12,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, tau, tau
        if (GETDER(13) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(13,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(13,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(13,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(13,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, tau, tau, tau
        if (GETDER(14) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(14,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(14,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(14,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(14,:,:) = 0.D0
        end if
        !Fourth derivative wrt tau
        if (GETDER(15) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(15,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(15,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(15,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(15,:,:) = 0.D0
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end if


    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 2)) then
        !Calculate the residual part of the activity coefficients
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !ln_gamma_R and all other quantities are AT THE MOMENT NOT always needed for the following calculations
        !As soon as higher derivatives are required, this may change!
        if ((gl%ge%GETDER_dxadxb_prev(1) == 1) .and. (recalc .eqv. .false.)) then

            gl%ge%ln_gamma_R_dxadxb(1,:,:,:) = gl%ge%ln_gamma_R_dxadxb_prev(1,:,:,:)

            !Transformation of variables from tau,delta, and x to T,rho,x
            do b = 1, gl%ncomp-1
                do a = 1, gl%ncomp-1
                    do i = 1, gl%ncomp
                        ln_gamma_R_Trho_dxadxb(1,i,a,b) = gl%ge%ln_gamma_R_dxadxb(1,i,a,b) &
                            & - dTred_dxa(a) * dTred_dxa(b) / tau**2 * ln_gamma_R_Trho(5,i) &
                            & - dTred_dxa(a) / tau * ln_gamma_R_Trho_dxa(4,i,b) &
                            & - dTred_dxa(b) / tau * ln_gamma_R_Trho_dxa(4,i,a) &
                            & - d2Tred_dxadxb(a,b) / tau * ln_gamma_R_Trho(4,i)
                    end do
                end do
            end do

        else

            !----------------------------------------------------------------------
            !Calculate the mole fractions of all groups in the mixture
            !----------------------------------------------------------------------
            sumsum_vnxj = 0.D0
            Do j = 1, gl%ncomp                     !sum over all components in the mixture
                Do n = 1, gl%nr_of_groups_i(j)     !Sum over all groups of component j
                    sumsum_vnxj = sumsum_vnxj + gl%v_ik(j, n) * x(j)
                End do
            End do

            m = 1
            Do j = 1, gl%ncomp                     !Sum over all components in the mixture
                Do k = 1, gl%nr_of_groups_i(j)     !Sum over all groups of component j
                    X_m(m) = gl%v_ik(j, k) * x(j) / sumsum_vnxj
                    m = m + 1
                End do
            End do
            !----------------------------------------------------------------------


            !----------------------------------------------------------------------
            !Calculate Theta_m for all groups in the mixture
            !----------------------------------------------------------------------
            m = 1
            sum_QnXn = 0.D0
            Do j = 1, gl%ncomp                      !Sum over all components in the mixture
                Do k = 1, gl%nr_of_groups_i(j)      !Sum over all groups of component j
                    sum_QnXn = sum_QnXn + gl%Q_ik(j, k) * X_m(m)
                    m = m + 1
                End do
            End do

            m = 1
            Do j = 1, gl%ncomp                      !Sum over all components in the mixture
                Do k = 1, gl%nr_of_groups_i(j)      !Sum over all groups of component j
                    Theta_m(m) = gl%Q_ik(j, k) * X_m(m) / sum_QnXn
                    !Set help variable Q_m to list all group areas in one vector in order to not have to map them back to Q_ik(i, k) for derivative calculation
                    Q_m(m) = gl%Q_ik(j, k)
                    m = m + 1
                End do
            End do
            !----------------------------------------------------------------------


            !----------------------------------------------------------------------
            !Calculate the derivative of the group mole fractions Xm with respect
            !to xa. This Derivative is needed for all following derivatives
            Do a = 1, gl%ncomp-1
                Do m = 1, gl%nr_of_groups_mix
                    !Three cases have to be distinguished:
                    !Case 1: Group m is part of molecule a
                    !Case 2: Group m is part of molecule N (Last component)
                    !Case 3: Group m is neither part of molecule a nor N

                    !Case 1
                    if (m >= gl%group_start_index(a) .and. m < gl%group_start_index(a+1)) then
                        !The 0.D0 in the following equation stands for the number of groups of type m in molecule ncomp, which is 0, because here m is a group of molecule a
                        !v_ik(a, m+1-group_start_index(a)) refers to the number of groups of type m in molecule of type a
                        dX_m_dxa(m,a) = ((gl%v_ik(a, m+1-gl%group_start_index(a)) - 0.D0)*sumsum_vnxj &
                            & - (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(a) * gl%v_ik(a, m+1-gl%group_start_index(a))) &
                            &  / sumsum_vnxj**2

                        !Case 2
                    elseif (m >= gl%group_start_index(gl%ncomp)) then
                        !The 0.D0 in the following equation stands for the number of groups of type m in molecule a, which is 0, because here m is a group of molecule N
                        !v_ik(ncomp, m+1-group_start_index(ncomp)) refers to the number of groups of type m in molecule of type ncomp
                        dX_m_dxa(m,a) = ((0.D0 - gl%v_ik(gl%ncomp, m+1-gl%group_start_index(gl%ncomp))) * sumsum_vnxj &
                            & - (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(gl%ncomp) * gl%v_ik(gl%ncomp, m+1-gl%group_start_index(gl%ncomp))) &
                            &  / sumsum_vnxj**2

                        !Case 3
                    else
                        !The first 0.D0 stands for the number of groups of type m in molecule a and the second 0.D0 stands for the number of groups of type m in molecule N.
                        !Both are 0, because the group belongs to another molecule of type k which needs to be determined
                        !Determine which molecule k the group m belongs to
                        do k=1,gl%ncomp-1
                            if (m >= gl%group_start_index(k) .and. m < gl%group_start_index(k+1)) then
                                exit
                            end if
                        end do
                        dX_m_dxa(m,a) = ((0.D0 - 0.D0) * sumsum_vnxj &
                            & - (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(k) * gl%v_ik(k, m+1-gl%group_start_index(k))) &
                            &  / sumsum_vnxj**2

                    end if

                end do
            end do
            !----------------------------------------------------------------------


            !----------------------------------------------------------------------
            !Calculate the derivative of Theta_m with respect to xa
            !This Derivative is needed for all following derivatives
            sum_QndXn_dxa = 0.D0
            do a = 1,gl%ncomp-1
                Do m = 1, gl%nr_of_groups_mix
                    sum_QndXn_dxa(a) = sum_QndXn_dxa(a) + Q_m(m) * dX_m_dxa(m,a)
                end do
            end do

            Do a = 1, gl%ncomp-1
                Do m = 1, gl%nr_of_groups_mix
                    dTheta_m_dxa(m,a) = (Q_m(m) * dX_m_dxa(m,a) * sum_QnXn - Q_m(m) * X_m(m) * sum_QndXn_dxa(a)) / sum_QnXn**2
                end do
            end do
            !----------------------------------------------------------------------




            !----------------------------------------------------------------------
            !Calculate the derivative of the group mole fractions Xm with respect
            !to xa and xb. This Derivative is needed for all following derivatives
            Do b = 1, gl%ncomp-1
                Do a = 1, gl%ncomp-1
                    Do m = 1, gl%nr_of_groups_mix
                        !Four cases have to be distinguished:
                        !Case 1: Group m is part of molecule a and b (i.e. a and b are the SAME molecule, a=b)
                        !Case 2: Group m is part of molecule a
                        !Case 3: Group m is part of molecule b
                        !Case 4: Group m is part of molecule N (Last component)
                        !Case 5: Group m is neither part of molecule a nor N

                        !Case 1
                        if ((a .eq. b) .and. (m >= gl%group_start_index(a) .and. m < gl%group_start_index(a+1))) then
                            !The 0.D0s in the following equation stand for the number of groups of type m in molecule ncomp and groups of type m in molecule b, which are 0, because here m is a group of molecule a
                            !v_ik(a, m+1-group_start_index(a)) refers to the number of groups of type m in molecule of type a
                            d2X_m_dxadxb(m,a,b) = (-(gl%v_ik(a, m+1-gl%group_start_index(a)) - 0.D0) * (gl%nr_of_groups_i(b) - gl%nr_of_groups_i(gl%ncomp)) &
                                & - (gl%v_ik(b, m+1-gl%group_start_index(b)) - 0.D0) * (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp))) &
                                & / sumsum_vnxj**2 &
                                & + 2.D0 * (gl%nr_of_groups_i(b) - gl%nr_of_groups_i(gl%ncomp)) &
                                & * (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(a) * gl%v_ik(a, m+1-gl%group_start_index(a)) &
                                & / sumsum_vnxj**3

                            !Case 2
                        elseif (m >= gl%group_start_index(a) .and. m < gl%group_start_index(a+1)) then
                            !The 0.D0s in the following equation stand for the number of groups of type m in molecule ncomp and groups of type m in molecule b, which are 0, because here m is a group of molecule a
                            !v_ik(a, m+1-group_start_index(a)) refers to the number of groups of type m in molecule of type a
                            d2X_m_dxadxb(m,a,b) = (-(gl%v_ik(a, m+1-gl%group_start_index(a)) - 0.D0) * (gl%nr_of_groups_i(b) - gl%nr_of_groups_i(gl%ncomp)) &
                                & - (0.D0 - 0.D0) * (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp))) &
                                & / sumsum_vnxj**2 &
                                & + 2.D0 * (gl%nr_of_groups_i(b) - gl%nr_of_groups_i(gl%ncomp)) &
                                & * (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(a) * gl%v_ik(a, m+1-gl%group_start_index(a)) &
                                & / sumsum_vnxj**3

                            !Case 3
                        elseif (m >= gl%group_start_index(b) .and. m < gl%group_start_index(b+1)) then
                            !The 0.D0s in the following equation stand for the number of groups of type m in molecule ncomp and groups of type m in molecule a, which are 0, because here m is a group of molecule b
                            !v_ik(b, m+1-group_start_index(b)) refers to the number of groups of type m in molecule of type b
                            d2X_m_dxadxb(m,a,b) = (-(0.D0 - 0.D0) * (gl%nr_of_groups_i(b) - gl%nr_of_groups_i(gl%ncomp)) &
                                & - (gl%v_ik(b, m+1-gl%group_start_index(b)) - 0.D0) * (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp))) &
                                & / sumsum_vnxj**2 &
                                & + 2.D0 * (gl%nr_of_groups_i(b) - gl%nr_of_groups_i(gl%ncomp)) &
                                & * (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(b) * gl%v_ik(b, m+1-gl%group_start_index(b)) &
                                & / sumsum_vnxj**3


                            !Case 4
                        elseif (m >= gl%group_start_index(gl%ncomp)) then
                            !The 0.D0s in the following equation stand for the number of groups of type m in molecule a and groups of type m in molecule b, which are 0, because here m is a group of molecule ncomp
                            !v_ik(ncomp, m+1-group_start_index(ncomp)) refers to the number of groups of type m in molecule of type ncomp
                            d2X_m_dxadxb(m,a,b) = (-(0.D0 - gl%v_ik(gl%ncomp, m+1-gl%group_start_index(gl%ncomp))) * (gl%nr_of_groups_i(b) - gl%nr_of_groups_i(gl%ncomp)) &
                                & - (0.D0 - gl%v_ik(gl%ncomp, m+1-gl%group_start_index(gl%ncomp))) * (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp))) &
                                & / sumsum_vnxj**2 &
                                & + 2.D0 * (gl%nr_of_groups_i(b) - gl%nr_of_groups_i(gl%ncomp)) &
                                & * (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(gl%ncomp) * gl%v_ik(gl%ncomp, m+1-gl%group_start_index(gl%ncomp)) &
                                & / sumsum_vnxj**3

                            !Case 5
                        else
                            !The 0.D0s in the following equation stand for the number of groups of type m in molecule a and groups of type m in molecule b and groups of type m in molecule ncomp, which are 0, because here m is a group of molecule k
                            !molecule of type k needs to be determined
                            !Determine which molecule k the group m belongs to
                            do k=1,gl%ncomp-1
                                if (m >= gl%group_start_index(k) .and. m < gl%group_start_index(k+1)) then
                                    exit
                                end if
                            end do
                            d2X_m_dxadxb(m,a,b) = (-(0.D0 - 0.D0) * (gl%nr_of_groups_i(b) - gl%nr_of_groups_i(gl%ncomp)) &
                                & - (0.D0 - 0.D0) * (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp))) &
                                & / sumsum_vnxj**2 &
                                & + 2.D0 * (gl%nr_of_groups_i(b) - gl%nr_of_groups_i(gl%ncomp)) &
                                & * (gl%nr_of_groups_i(a) - gl%nr_of_groups_i(gl%ncomp)) * x(k) * gl%v_ik(k, m+1-gl%group_start_index(k)) &
                                & / sumsum_vnxj**3

                        end if

                    end do
                end do
            end do
            !----------------------------------------------------------------------

            !----------------------------------------------------------------------
            !Calculate the derivative of Theta_m with respect to xa and xb
            !This Derivative is needed for all following derivatives
            sum_QndXn_dxadxb = 0.D0
            do b = 1,gl%ncomp-1
                do a = 1,gl%ncomp-1
                    Do m = 1, gl%nr_of_groups_mix
                        sum_QndXn_dxadxb(a,b) = sum_QndXn_dxadxb(a,b) + Q_m(m) * d2X_m_dxadxb(m,a,b)
                    end do
                end do
            end do

            do b = 1, gl%ncomp-1
                Do a = 1, gl%ncomp-1
                    Do m = 1, gl%nr_of_groups_mix
                        d2Theta_m_dxadxb(m,a,b) = (Q_m(m) * d2X_m_dxadxb(m,a,b) * sum_QnXn - Q_m(m) * dX_m_dxa(m,a) * sum_QndXn_dxa(b) &
                            & - Q_m(m) * dX_m_dxa(m,b) * sum_QndXn_dxa(a) - Q_m(m) * X_m(m) * sum_QndXn_dxadxb(a,b) ) / sum_QnXn**2 &
                            & + (2.D0 *  sum_QndXn_dxa(b) * Q_m(m) * X_m(m) * sum_QndXn_dxa(a)) / sum_QnXn**3
                    end do
                end do
            end do
            !----------------------------------------------------------------------



            !----------------------------------------------------------------------
            !Calculate all group interaction parameters Psi_nm
            !----------------------------------------------------------------------
            Do n = 1, gl%nr_of_groups_mix
                Do m = 1, gl%nr_of_groups_mix
                    Psi_nm(n, m) = dexp(-(gl%a_nm(n, m) + gl%b_nm(n, m) * Temp + gl%c_nm(n, m) * Temp**2) / Temp)
                End do
            End do
            !----------------------------------------------------------------------


            !Loop over all molefactions xa and xb
            do b = 1, gl%ncomp-1
                Do a = 1, gl%ncomp-1
                    help_k = 0     !Help index, because for Psi_nm a different index for k is needed (k runs through the groups of each molecule, help_k runs through ALL groups)
                    !Loop over all components in the mixture
                    Do i = 1, gl%ncomp

                        !----------------------------------------------------------------------
                        !Calculate dlnGamma_k_dxa
                        !Also needed for the following calculations
                        !----------------------------------------------------------------------
                        Do k = 1, gl%nr_of_groups_i(i)
                            sum_ThetamPsimk = 0.D0
                            sum_dThetam_dxa_Psimk = 0.D0
                            sum_dThetam_dxb_Psimk = 0.D0
                            sum_d2Thetam_dxadxb_Psimk = 0.D0
                            sumsum_ThetamPsikm_dxadxb = 0.D0
                            help_k = help_k + 1
                            Do m = 1, gl%nr_of_groups_mix
                                sum_ThetamPsimk = sum_ThetamPsimk + Theta_m(m) * Psi_nm(m, help_k)
                                sum_dThetam_dxa_Psimk = sum_dThetam_dxa_Psimk + dTheta_m_dxa(m,a) * Psi_nm(m, help_k)
                                sum_dThetam_dxb_Psimk = sum_dThetam_dxb_Psimk + dTheta_m_dxa(m,b) * Psi_nm(m, help_k)
                                sum_d2Thetam_dxadxb_Psimk = sum_d2Thetam_dxadxb_Psimk + d2Theta_m_dxadxb(m,a,b) * Psi_nm(m, help_k)
                                sum_ThetanPsinm = 0.D0
                                sum_dThetan_dxa_Psinm = 0.D0
                                sum_dThetan_dxb_Psinm = 0.D0
                                sum_d2Thetan_dxadxb_Psinm = 0.D0
                                Do n = 1, gl%nr_of_groups_mix
                                    sum_ThetanPsinm = sum_ThetanPsinm + Theta_m(n) * Psi_nm(n, m)
                                    sum_dThetan_dxa_Psinm = sum_dThetan_dxa_Psinm + dTheta_m_dxa(n,a) * Psi_nm(n, m)
                                    sum_dThetan_dxb_Psinm = sum_dThetan_dxb_Psinm + dTheta_m_dxa(n,b) * Psi_nm(n, m)
                                    sum_d2Thetan_dxadxb_Psinm = sum_d2Thetan_dxadxb_Psinm + d2Theta_m_dxadxb(n,a,b) * Psi_nm(n, m)
                                End do
                                sumsum_ThetamPsikm_dxadxb = sumsum_ThetamPsikm_dxadxb &
                                    & - ( (d2Theta_m_dxadxb(m,a,b) * Psi_nm(help_k, m) * sum_ThetanPsinm &
                                    & - dTheta_m_dxa(m,a) * Psi_nm(help_k, m) * sum_dThetan_dxb_Psinm &
                                    & - dTheta_m_dxa(m,b) * Psi_nm(help_k, m) * sum_dThetan_dxa_Psinm &
                                    & - Theta_m(m) * Psi_nm(help_k, m) * sum_d2Thetan_dxadxb_Psinm) &
                                    & / sum_ThetanPsinm**2 &
                                    & + 2.D0 * Theta_m(m) * Psi_nm(help_k, m) &
                                    & * sum_dThetan_dxa_Psinm * sum_dThetan_dxb_Psinm / sum_ThetanPsinm**3 )
                            End do
                            d2lnGamma_k_dxadxb(k,a,b) = gl%Q_ik(i, k) * sum_dThetam_dxa_Psimk * sum_dThetam_dxb_Psimk / sum_ThetamPsimk**2 &
                                & - gl%Q_ik(i, k) * sum_d2Thetam_dxadxb_Psimk / sum_ThetamPsimk &
                                & + gl%Q_ik(i, k) * sumsum_ThetamPsikm_dxadxb
                        End do
                        !----------------------------------------------------------------------


                        !Calculate the derivative of the logarithm of the residual activity coefficients with respect to xa
                        !!!NOTE THAT LNGAMMA_ki OF PURE COMPONENT i IS NOT A FUNCTION OF COMPOSITION, HENCE dlnGamma_ki_dxa(k,i,:) = 0!!! This is indicated with the 0.D0 in the following equation
                        Do k = 1, gl%nr_of_groups_i(i)
                            ln_gamma_R_Trho_dxadxb(1,i,a,b) = ln_gamma_R_Trho_dxadxb(1,i,a,b) + gl%v_ik(i, k) * (d2lnGamma_k_dxadxb(k,a,b) - 0.D0)
                        End do


                    End do
                end do
            end do

            !Transformation of variables from T,rho,x to tau,delta,x
            do b = 1, gl%ncomp-1
                do a = 1, gl%ncomp-1
                    do i = 1, gl%ncomp
                        gl%ge%ln_gamma_R_dxadxb(1,i,a,b) = ln_gamma_R_Trho_dxadxb(1,i,a,b) &
                            & + dTred_dxa(a) * dTred_dxa(b) / tau**2 * ln_gamma_R_Trho(5,i) &
                            & + dTred_dxa(a) / tau * ln_gamma_R_Trho_dxa(4,i,b) &
                            & + dTred_dxa(b) / tau * ln_gamma_R_Trho_dxa(4,i,a) &
                            & + d2Tred_dxadxb(a,b) / tau * ln_gamma_R_Trho(4,i)
                    end do
                end do
            end do

        end if

        !Note in the following that all derivatives with respect to delta at constant tau and x are 0, because the residual part is not a function of density
        if (GETDER(1) .eq. 1) then

            if ((gl%ge%GETDER_dxadxb_prev(1) == 1) .and. (recalc .eqv. .false.)) then

                gl%ge%gE_R_dxadxb(1,:,:) = gl%ge%gE_R_dxadxb_prev(1,:,:)

                !Transformation of variables from tau,delta, and x to T,rho,x
                do b = 1, gl%ncomp-1
                    do a = 1, gl%ncomp-1
                        gE_R_Trho_dxadxb(1,a,b) = gl%ge%gE_R_dxadxb(1,a,b) &
                            & - dTred_dxa(a) * dTred_dxa(b) / tau**2 * gE_R_Trho(5) &
                            & - dTred_dxa(a) / tau * gE_R_Trho_dxa(4,b) &
                            & - dTred_dxa(b) / tau * gE_R_Trho_dxa(4,a) &
                            & - d2Tred_dxadxb(a,b) / tau * gE_R_Trho(4)
                    end do
                end do

            else

                Do b = 1, gl%ncomp-1
                    Do a = 1,gl%ncomp-1
                        Do i = 1, gl%ncomp
                            gE_R_Trho_dxadxb(1,a,b) = gE_R_Trho_dxadxb(1,a,b) + R_const * Temp * x(i) * ln_gamma_R_Trho_dxadxb(1,i,a,b)
                        end do
                        gE_R_Trho_dxadxb(1,a,b) = gE_R_Trho_dxadxb(1,a,b) &
                            & + R_const * Temp * (ln_gamma_R_Trho_dxa(1,a,b) - ln_gamma_R_Trho_dxa(1,gl%ncomp,b)) &
                            & + R_const * Temp * (ln_gamma_R_Trho_dxa(1,b,a) - ln_gamma_R_Trho_dxa(1,gl%ncomp,a))
                    End do
                End do

                !Transformation of variables from T,rho,x to tau,delta,x
                do b = 1, gl%ncomp-1
                    do a = 1, gl%ncomp-1
                        gl%ge%gE_R_dxadxb(1,a,b) = gE_R_Trho_dxadxb(1,a,b) &
                            & + dTred_dxa(a) * dTred_dxa(b) / tau**2 * gE_R_Trho(5) &
                            & + dTred_dxa(a) / tau * gE_R_Trho_dxa(4,b) &
                            & + dTred_dxa(b) / tau * gE_R_Trho_dxa(4,a) &
                            & + d2Tred_dxadxb(a,b) / tau * gE_R_Trho(4)
                    end do
                end do

            end if

        end if

        !Second derivative wrt delta and tau
        if (GETDER(2) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(2,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(2,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(2,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(2,:,:) = 0.D0
        end if

        !Third derivative wrt delta, tau, and tau
        if (GETDER(3) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(3,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(3,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(3,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(3,:,:) = 0.D0
        end if

        !Derivative wrt tau
        if (GETDER(4) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(4,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(4,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(4,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(4,:,:) = 0.D0
        end if

        !Second derivative wrt tau
        if (GETDER(5) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(5,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(5,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(5,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(5,:,:) = 0.D0
        end if
        !Second derivative wrt delta and tau
        if (GETDER(6) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(6,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(6,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(6,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(6,:,:) = 0.D0
        end if
        !Third derivative wrt delta, tau, and tau
        if (GETDER(7) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(7,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(7,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(7,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(7,:,:) = 0.D0
        end if
        !Third derivative wrt delta
        if (GETDER(8) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(8,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(8,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(8,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(8,:,:) = 0.D0
        end if
        !Third derivative wrt tau
        if (GETDER(9) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(9,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(9,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(9,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(9,:,:) = 0.D0
        end if
        !Third derivative wrt tau, delta, and delta
        if (GETDER(10) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(10,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(10,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(10,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(10,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta
        if (GETDER(11) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(11,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(11,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(11,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(11,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, delta and tau
        if (GETDER(12) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(12,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(12,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(12,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(12,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, tau, tau
        if (GETDER(13) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(13,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(13,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(13,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(13,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, tau, tau, tau
        if (GETDER(14) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(14,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(14,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(14,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(14,:,:) = 0.D0
        end if
        !Fourth derivative wrt tau
        if (GETDER(15) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(15,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(15,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(15,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(15,:,:) = 0.D0
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    End if


    !Write new values in save variables
    !gl%ge%Temp_dxadxb_prev = Temp
    !gl%ge%GETDER_dxadxb_prev = GETDER
    !gl%ge%molfrac_dxadxb_prev = gl%molfractions
    !gl%ge%Eq_type_dxadxb_prev = gl%Eq_type
    !gl%ge%mixtype_dxadxb_prev = gl%mix_type
    !gl%ge%C_or_R_dxadxb_prev = C_or_R
    !gl%ge%gE_C_dxadxb_prev = gl%ge%gE_C_dxadxb
    !gl%ge%gE_R_dxadxb_prev = gl%ge%gE_R_dxadxb
    !gl%ge%ln_gamma_C_dxadxb_prev = gl%ge%ln_gamma_C_dxadxb
    !gl%ge%ln_gamma_R_dxadxb_prev = gl%ge%ln_gamma_R_dxadxb


    end subroutine
    !--------------------------------------------------------------------------------------



    !Erik, July 2018
    !Derivatives of COSMO-SAC with respect to mole frations
    !Needed for phase equilibrium routines for new Helmholtz+gE model
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine gE_COSMO_SAC_MIXDERIVS_dxa(gl,Temp, GETDER, C_or_R, errval)
    !------------------------------------------------------------------------------------------------------------------------------
    !Subroutine to calculate the excess Gibbs energy and the activity coefficients of the
    !components in the mixture according to the COSMO-SAC model as well as tau- and delta-
    !derivatives of the model.
    !All quantities are given as derivatives with respect to mole fractions xa for all components
    !in the mixture. The convention that the "last" mole fraction xN of component N is replaced with
    !1-sum(xk) is used here.


    implicit none

    type(type_gl) :: gl


    double precision:: temp
    integer, dimension(nderivs):: GETDER                !array specifier to indicate, which derivative is needed
    integer:: C_or_R
    integer:: errval

    double precision:: delta_x                          !Step in x for numerical derivative
    double precision, dimension(30):: molfractions_orig !Original value of molfractions

    double precision :: R_const, R_const_cal            !universal gas constant, needs to be in calories for cosmo version 1

    !Derivatives of the combinatorial part with respect to "natural" variables
    double precision, dimension(nderivs, 30) :: ln_gamma_C_Trho       !Combinatorial activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs, 30) :: ln_gamma_R_Trho       !Residual activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs) :: gE_C_Trho        !Combinatorial part of gE and derivatives with respect to T and rho
    double precision, dimension(nderivs) :: gE_R_Trho        !Residual part of gE and derivatives with respect to T and rho
    double precision, allocatable:: ln_gamma_C_Trho_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to rho and T, and xa
    double precision, allocatable:: ln_gamma_R_Trho_dxa(:,:,:)       !Residual activity coefficients and derivatives with respect to rho and T, and xa
    double precision, dimension(nderivs, 30) :: gE_C_Trho_dxa        !Combinatorial part of gE and derivatives with respect to T and rho and xa
    double precision, dimension(nderivs, 30) :: gE_R_Trho_dxa        !Residual part of gE and derivatives with respect to T and rho and xa

    !Required derivatives of activities and excess Gibbs energy only with respect to tau and delta (not xa, however, needed for the calculation of derivatives with respect to xa)
    integer, dimension(nderivs):: GETDER_no_xa                !array specifier to indicate, which derivative is needed

    !Help variables
    double precision, dimension(30) :: x
    integer :: i, j, k, m, n, o, a
    double precision, dimension(51,30) :: dseggammadxa_mix, dseggammadxa_mix_gauss!, dseggammadxa_mix_test
    double precision, dimension (51,30) :: Fi_xa
    double precision :: Det_seggamma, Det_xa
    !double precision, dimension(60,60):: MatrixA, MatrixB, dMatrixAdxa, dMatrixBdxa, MatrixC, MatrixD!, MatrixA_p, MatrixA_m, MatrixB_p, MatrixB_m
    double precision, dimension(:,:), allocatable:: MatrixA
    double precision, dimension(:,:), allocatable:: MatrixB
    double precision, dimension(:,:), allocatable:: MatrixC
    double precision, dimension(:,:), allocatable:: MatrixD
    !double precision, dimension(51,60,60):: MatrixA_mix_p, MatrixA_mix_m
    !November 2018, Erik, allocatable
    !double precision, dimension(60,60,51,30) :: MatrixA_mix
    double precision, dimension(:,:,:,:), allocatable :: MatrixA_mix

    !double precision, dimension(60,60) ::MatrixB_mix!, MatrixB_mix_p, MatrixB_mix_m
    double precision, dimension(:,:), allocatable ::MatrixB_mix
    !double precision, dimension(51,30,60,60) :: MatrixA_pure, MatrixA_pure_p, MatrixA_pure_m
    !double precision, dimension(30,60,60) :: MatrixB_pure, MatrixB_pure_p, MatrixB_pure_m
    integer:: rankA
    double precision, dimension(60):: vectorb, vectorb_gauss
    double precision, dimension(60):: vectorx
    integer:: errorflag
    double precision, dimension(15,30) :: ln_gamma_R_test
    double precision, dimension(51,30) :: dlnseggammadxa_mix, dlnseggammadxa_tau_mix!, dlnseggammadtau_mix
    !double precision, dimension(51,30) :: dlnseggammadxa_pure, dlnseggammadxa_tau_pure!, dlnseggammadtau_pure
    integer, dimension(51) :: pos

    double precision :: sum_xiAi
    double precision, dimension(51) :: sum_xiAisigma
    double precision, dimension(15,30) :: gE_R_dxa_test
    !double precision, dimension(51,30) :: DetMatrixA_pure! , DetMatrixA_pure_p, DetMatrixA_pure_m
    !double precision, dimension(30) :: DetMatrixB_pure!, DetMatrixB_pure_p, DetMatrixB_pure_m
    double precision, dimension(51,30) :: DetMatrixA_mix!, DetMatrixA_mix_test!, DetMatrixA_mix_p, DetMatrixA_mix_m
    double precision :: DetMatrixB_mix

    !help variables mixed derivation dxa*dT
    double precision, dimension(51) :: d2Fi_dxa_dT_mix
    double precision, dimension(51) :: sum_d2Fi_dseggammadT
    !November 2018, Erik, allocatable
    !double precision, dimension(51,51) :: d2Fi_dseggammadT_mix
    double precision, dimension(:,:), allocatable :: d2Fi_dseggammadT_mix
    double precision, dimension(51,30) :: dDetAdT_mix, d2seggamma_dxa_dT_mix, d2lnseggamma_dxa_dT_mix
    !double precision, dimension(60,60):: adj_A, adj_B, dMatrixBdT, dMatrixAdT
    double precision, dimension(:,:), allocatable:: adj_A
    double precision, dimension(:,:), allocatable:: adj_B
    double precision, dimension(:,:), allocatable:: dMatrixBdT
    double precision, dimension(:,:), allocatable:: dMatrixAdT

    double precision :: dDetBdT_mix
    double precision, dimension(51,30) :: dseggammadxa_mix_num

    double precision, dimension(30):: dTred_dxa, drhored_dxa    !Variables for the derivatives of the reducing functions with respect to molfractions xa

    !Corresponding states relevant variables
    double precision:: tau, del

    !Variable that indicates if derivatives need to be recalculated
    logical:: recalc

    !COSMO-SAC variables
    double precision, dimension(51) :: sigma_profile_mix, seggamma

    !double precision, dimension(51,51) :: dFi_dseggamma
    double precision, dimension(:,:), allocatable :: dFi_dseggamma


    !numerical test of d2seggamma_dxa_dT_mix
    double precision :: Temp_p, Temp_m
    !double precision, dimension(51,30) :: dseggamma_dxa_p, dseggamma_dxa_m
    !double precision, dimension(51,51,30) :: Fi_seggamma_mix_p, Fi_seggamma_mix_m, d2Fdseggammadxa_mix
    double precision, dimension(51) :: d2FdT2_mix
    double precision, dimension(51,30) :: seggamma_p, seggamma_m!, seggamma_pure_p, seggamma_pure_m, d2FdT2_pure, Fi_T_pure_p, Fi_T_pure_m,  d2Fdxadxb_mix
    !double precision, dimension(51,51,30) :: Fi_seggamma_pure_p, Fi_seggamma_pure_m, d2Fi_dseggammadT_pure, d2FdTdseggamma_pure
    !double precision, dimension(51,30) :: dDetAdT_mix_num, dDetBdT_mix_num
    !double precision, dimension(51,30) :: DetMatrixA_mix_p, DetMatrixA_mix_m, d2seggammadxadT_mix_num, dseggammadxa_mix_m, dseggammadxa_mix_p
    !double precision :: DetMatrixB_mix_p, DetMatrixB_mix_m
    !double precision, dimension(51,30) :: d2lnseggammadxadT_mix_num, dlnseggammadxa_mix_p, dlnseggammadxa_mix_m, Fi_xa_p, Fi_xa_m
    !double precision, dimension(15,30,30) :: ln_gamma_R_Trho_dxa_num, ln_gamma_R_Trho_dxa_p, ln_gamma_R_Trho_dxa_m

    double precision :: delta_tau

    !test variables
    double precision :: tredmix_orig, rhoredmix_orig, Temperature
    integer, dimension(nderivs):: getder_dep
    double precision, dimension(51,30) :: d2segamma_dxadxb_num
    double precision, dimension(51,30) :: dseggamma_dxa
    double precision, dimension(30,30) :: x_p, x_m
    double precision, dimension(30) :: sum_xiAi_p, sum_xiAi_m
    double precision, dimension(51,30) :: sum_xiAisigma_p, sum_xiAisigma_m
    double precision, dimension(nderivs) :: gE_R_p, gE_R_m
    double precision, dimension(15,30) :: ln_gamma_R_p, ln_gamma_R_m

    !allocation of large arrays
    if(.not. allocated(gl%ge%ln_gamma_C_dxa)) allocate(gl%ge%ln_gamma_C_dxa(nderivs,30,30))
    if(.not. allocated(gl%ge%ln_gamma_R_dxa)) allocate(gl%ge%ln_gamma_R_dxa(nderivs,30,30))
    if(.not. allocated(ln_gamma_C_Trho_dxa)) allocate(ln_gamma_C_Trho_dxa(nderivs,30,30))
    if(.not. allocated(ln_gamma_R_Trho_dxa)) allocate(ln_gamma_R_Trho_dxa(nderivs,30,30))

    !November 2018, Erik, allocatable
    if(.not. allocated(MatrixA_mix)) allocate(MatrixA_mix(60,60,51,30))
    if(.not. allocated(dFi_dseggamma)) allocate(dFi_dseggamma(51,51))
    if(.not. allocated(MatrixB_mix)) allocate(MatrixB_mix(60,60))
    if(.not. allocated(d2Fi_dseggammadT_mix)) allocate(d2Fi_dseggammadT_mix(51,51))
    if(.not. allocated(adj_A)) allocate(adj_A(60,60))
    if(.not. allocated(adj_B)) allocate(adj_B(60,60))
    if(.not. allocated(dMatrixBdT)) allocate(dMatrixBdT(60,60))
    if(.not. allocated(dMatrixAdT)) allocate(dMatrixAdT(60,60))
    if(.not. allocated(MatrixA)) allocate(MatrixA(60,60))
    if(.not. allocated(MatrixB)) allocate(MatrixB(60,60))
    if(.not. allocated(MatrixC)) allocate(MatrixC(60,60))
    if(.not. allocated(MatrixD)) allocate(MatrixD(60,60))



    errval = 0

    !Check for wrong input to C_or_R
    if ((C_or_R < 0) .or. (C_or_R > 2)) then
        errval = -1111
        gl%ge%gE_C_dxa = errval
        gl%ge%gE_R_dxa = errval
        gl%ge%ln_gamma_C_dxa = errval
        gl%ge%ln_gamma_R_dxa = errval
        return
    end if

    recalc = .true.

    !Check whether calculations have to be redone or not.
    !If the temperature, (density), composition, equation of state and, the required part of COSMO-SAC stayed the same, no need to recalculate
    if(dabs(gl%ge%Temp_dxa_prev - Temp) > 1.D-16) then
        recalc = .true.
    end if
    do i=1,gl.ncomp
        if(dabs(gl%ge%molfrac_dxa_prev(i) - gl.molfractions(i)) > 1.D-16) then
            recalc = .true.
        end if
        if(gl%ge%Eq_type_dxa_prev(i) .ne. gl.Eq_type(i)) then
            recalc = .true.
        end if
    end do
    if (gl%ge%mixtype_dxa_prev .ne. gl.mix_type) then
        recalc = .true.
    end if
    if (gl%ge%C_or_R_dxa_prev .ne. C_or_R) then
        recalc = .true.
    end if


    !Initialize variables
    ln_gamma_C_Trho = 0.D0
    ln_gamma_R_Trho = 0.D0
    gE_C_Trho = 0.D0
    gE_R_Trho = 0.d0
    ln_gamma_C_Trho_dxa = 0.D0
    ln_gamma_R_Trho_dxa = 0.d0
    gE_C_Trho_dxa = 0.D0
    gE_R_Trho_dxa = 0.D0
    gl%ge%ln_gamma_C_dxa = 0.D0
    gl%ge%ln_gamma_R_dxa = 0.D0
    gl%ge%gE_C_dxa = 0.D0
    gl%ge%gE_R_dxa = 0.D0

    !Value for the ideal gas constant
    R_const = 8.3144598D0
    R_const_cal = 8.3144598D0 / 4184.D0

    !Get mole fractions of the mixture from module variable
    x = gl.molfractions

    !Calculate tau
    tau = gl.tredmix / Temp

    !Set dummy value for delta (Not needed because UNIFAC is not a function of density. However, already implemented here for possible later modifications)
    del = 1.D0

    !Derivatives of tredmix with respect to xa are required for the following derivatives
    call dYr_dxi(gl,dTred_dxa, drhored_dxa)


    !Get the required derivatives of activities and excess Gibbs energies with respect tau (and not with respect to xa, however, needed for the derivatives here)
    GETDER_no_xa = 0
    if (GETDER(1) .eq. 1)  then
        GETDER_no_xa(1) = 1
        GETDER_no_xa(4) = 1
    end if
    if (GETDER(4) .eq. 1)  then
        GETDER_no_xa(1) = 1
        GETDER_no_xa(4) = 1
        GETDER_no_xa(5) = 1
    end if

    call gE_COSMO_SAC_MIXDERIVS(gl,Temp, GETDER_no_xa, C_or_R, errval)

    if (GETDER(1) .eq. 1) then
        !Transform the tau derivatives into the needed temperature derivatives
        ln_gamma_C_Trho(1,:) = gl%ge%ln_gamma_C(1,:)
        gE_C_Trho(1) = gl%ge%gE_C(1)
        ln_gamma_R_Trho(1,:) = gl%ge%ln_gamma_R(1,:)
        gE_R_Trho(1) = gl%ge%gE_R(1)
        ln_gamma_C_Trho(4,:) = - tau**2 / gl.tredmix * gl%ge%ln_gamma_C(4,:)
        gE_C_Trho(4) = - tau**2 / gl.tredmix * gl%ge%gE_C(4)
        ln_gamma_R_Trho(4,:) = - tau**2 / gl.tredmix * gl%ge%ln_gamma_R(4,:)
        gE_R_Trho(4) = - tau**2 / gl.tredmix * gl%ge%gE_R(4)

        !COSMO-SAC variables needed here
        sigma_profile_mix = gl.cosmo.sigma_profile_mix_gl
        dFi_dseggamma = gl.cosmo.dFi_mix_v1_gl
        seggamma = gl.cosmo.seggamma_gl
    end if
    if (GETDER(4) .eq. 1) then
        !Transform the tau derivatives into the needed temperature derivatives
        ln_gamma_C_Trho(1,:) = gl%ge%ln_gamma_C(1,:)
        gE_C_Trho(1) = gl%ge%gE_C(1)
        ln_gamma_R_Trho(1,:) = gl%ge%ln_gamma_R(1,:)
        gE_R_Trho(1) = gl%ge%gE_R(1)
        ln_gamma_C_Trho(4,:) = - tau**2 / gl.tredmix * gl%ge%ln_gamma_C(4,:)
        gE_C_Trho(4) = - tau**2 / gl.tredmix * gl%ge%gE_C(4)
        ln_gamma_R_Trho(4,:) = - tau**2 / gl.tredmix * gl%ge%ln_gamma_R(4,:)
        gE_R_Trho(4) = - tau**2 / gl.tredmix * gl%ge%gE_R(4)
        ln_gamma_C_Trho(5,:) = tau**4 / gl.tredmix**2 * gl%ge%ln_gamma_C(5,:) - 2.D0 * tau / gl.tredmix * ln_gamma_C_Trho(4,:)
        gE_C_Trho(5) = tau**4 / gl.tredmix**2 * gl%ge%gE_C(5) - 2.D0 * tau / gl.tredmix * gE_C_Trho(4)
        ln_gamma_R_Trho(5,:) = tau**4 / gl.tredmix**2 * gl%ge%ln_gamma_R(5,:) - 2.D0 * tau / gl.tredmix * ln_gamma_R_Trho(4,:)
        gE_R_Trho(5) = tau**4 / gl.tredmix**2 * gl%ge%gE_R(5) - 2.D0 * tau / gl.tredmix * gE_R_Trho(4)
    end if
    !------------------------------------------------------------------------------------------------------
    !!numerical test second derivatives in respect to xa and T, remove after testing
    !
    !    delta_tau = 1.d-4
    !    Temp_p = gl.tredmix / (tau + delta_tau)
    !    Temp_m = gl.tredmix / (tau - delta_tau)
    !    call gE_COSMO_SAC_CALC(gl, Temp_p, C_or_R, errval)
    !    !gE_C_p = gl.gE_C
    !    !gE_R_p = gl.gE_R
    !    !ln_gamma_C_p = gl.ln_gamma_C
    !    !ln_gamma_R_p = gl.ln_gamma_R
    !    !seggamma_pure_p = gl.cosmo.seggamma_pure_gl
    !    seggamma_p = gl.cosmo.seggamma_gl
    !    !Fi_seggamma_pure_p = gl.cosmo.dFi_pure_v1_gl
    !    Fi_seggamma_mix_p = gl.cosmo.dFi_mix_v1_gl
    !    call gE_COSMO_SAC_CALC(gl, Temp_m, C_or_R, errval)
    !    !gE_C_m = gl.gE_C
    !    !gE_R_m = gl.gE_R
    !    !ln_gamma_C_m = gl.ln_gamma_C
    !    !ln_gamma_R_m = gl.ln_gamma_R
    !    !seggamma_pure_m = gl.cosmo.seggamma_pure_gl
    !    seggamma_m = gl.cosmo.seggamma_gl
    !    !Fi_seggamma_pure_m = gl.cosmo.dFi_pure_v1_gl
    !    Fi_seggamma_mix_m = gl.cosmo.dFi_mix_v1_gl
    !
    !Numerical derivative of the excess based departure function with respect to xi at constant tau and del
    delta_x = 1.D-4
    molfractions_orig = gl.molfractions
    tredmix_orig = gl.tredmix
    rhoredmix_orig = gl.rhoredmix
    getder_dep = 1
    Temperature = Temp
    do i =1, gl.ncomp-1
        !Increase xi
        gl.molfractions(i) = molfractions_orig(i) + delta_x
        gl.molfractions(gl.ncomp) = molfractions_orig(gl.ncomp) - delta_x
        x_p(:,i) = gl.molfractions
        call reduced_parameters_calc(gl, Temperature)
        Temp_p = gl.tredmix/tredmix_orig * Temperature
        Temp_p = Temp
        !Dens_p = rhoredmix/rhoredmix_orig * Density
        call gE_COSMO_SAC_MIXDERIVS(gl, Temp_p, getder_dep, C_or_R, errval)
        !gE_C_p = gl.gE_C
        gE_R_p = gl%ge%gE_R
        !ln_gamma_C_p = gl.ln_gamma_C
        ln_gamma_R_p = gl%ge%ln_gamma_R
        !test
        !dseggamma_dxa_p= gl.cosmo.dseggamma_dxa_mix
        seggamma_p(:,i) = gl.cosmo.seggamma_gl
        !Fi_seggamma_mix_p(:,:,i) = gl.cosmo.dFi_mix_v1_gl

        !    alpha_dep_p = der_dep(1)

        !Decrease xi
        gl.molfractions(i) = molfractions_orig(i) - delta_x
        gl.molfractions(gl.ncomp) = molfractions_orig(gl.ncomp) + delta_x
        x_m(:,i) = gl.molfractions
        call reduced_parameters_calc(gl, Temperature)
        Temp_m = gl.tredmix/tredmix_orig * Temperature
        Temp_m = Temp
        !    Dens_m = rhoredmix/rhoredmix_orig * Density
        call gE_COSMO_SAC_MIXDERIVS(gl, Temp_m, getder_dep, C_or_R, errval)
        !gE_C_m = gl.gE_C
        gE_R_m = gl%ge%gE_R
        !ln_gamma_C_m = gl.ln_gamma_C
        ln_gamma_R_m = gl%ge%ln_gamma_R
        !test
        !dseggamma_dxa_m= gl.cosmo.dseggamma_dxa_mix
        seggamma_m(:,i) = gl.cosmo.seggamma_gl
        !Fi_seggamma_mix_m(:,:,i) = gl.cosmo.dFi_mix_v1_gl

        !    alpha_dep_m = der_dep(1)
        !Calculate the numerical derivative
        !gl.gE_C_dxa(1,i) = (gE_C_p(1) - gE_C_m(1)) / (2.D0 * delta_x)
        gl%ge%gE_R_dxa(1,i) = (gE_R_p(1) - gE_R_m(1)) / (2.D0 * delta_x)
        !test
        dseggamma_dxa(:,i) = (seggamma_p(:,i) - seggamma_m(:,i)) / (2.D0 * delta_x)
        !d2segamma_dxadxb_num = (dseggamma_dxa_p - dseggamma_dxa_m) / (2.D0 * delta_x)
        gl.molfractions = molfractions_orig
    end do

    call reduced_parameters_calc(gl, Temperature)
    !------------------------------------------------------------------------------------------------------


    !!set all COSMO-SAC values to orig
    !!can be removed after successful testing of analytical derivs
    call gE_COSMO_SAC_CALC(gl, temp, C_or_R, errval)
    call gE_COSMO_SAC_MIXDERIVS(gl, Temp, getder_dep, C_or_R, errval)
    !getder_dep = 1
    !call gE_COSMO_SAC_MIXDERIVS(gl, Temp, GETDER_no_xa, C_or_R, errval)
    !skipped all the combinatorial contributions for now as they are not needed
    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 2)) then
        !Calculate the residual part of the activity coefficients
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !ln_gamma_R is always needed for the following calculations, thus it is always calculated

        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        !mixture
        !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        ! Calculation Determinante of Matrix B (only needs to be calculated ones)
        !Matrix partial derivatives in respect to seggamma only
        m = 0
        do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
            n = 0
            if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                cycle
            end if
            m = m + 1
            do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                    cycle
                end if
                n = n + 1
                MatrixB(m,n) = dFi_dseggamma(j,k)
                MatrixB_mix(j,k) = MatrixB(m,n)
            end do
        end do
        rankA = m
        call Gauss_algorithm(gl,MatrixB, rankA, vectorb, vectorx, Det_seggamma, errval)
        DetMatrixB_mix = Det_seggamma

        !-------------------------------------------------------------------------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!test, REMOVE AFTER TESTING
        ! ! Calculation Determinante of Matrix B (only needs to be calculated ones)
        !!Matrix partial derivatives in respect to seggamma only
        !do a = 1, gl.ncomp
        !    m = 0
        !    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
        !        n = 0
        !        if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
        !            cycle
        !        end if
        !        m = m + 1
        !        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
        !            if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
        !                cycle
        !            end if
        !            n = n + 1
        !            !MatrixB_mix_p(j,k) = Fi_seggamma_mix_p(j,k)
        !            MatrixB_p(m,n) = Fi_seggamma_mix_p(j,k,a)
        !            MatrixB_m(m,n) = Fi_seggamma_mix_m(j,k,a)
        !            !MatrixB_mix_m(j,k) = Fi_seggamma_mix_m(j,k)
        !        end do
        !    end do
        !    rankA = m
        !    !-------------------------------------------------------------------------------------
        !    !nummerical derivatives test
        !    call Gauss_algorithm(gl,MatrixB_p(:,:,a), rankA, vectorb, vectorx, Det_seggamma, errval)
        !    DetMatrixB_mix_p(a) = Det_seggamma
        !    !nummerical derivatives test
        !    call Gauss_algorithm(gl,MatrixB_m(:,:,a), rankA, vectorb, vectorx, Det_seggamma, errval)
        !    DetMatrixB_mix_m(a) = Det_seggamma
        !
        !    dDetBdT_mix_num(:,a) = (DetMatrixB_mix_p(a) - DetMatrixB_mix_m(a)) / (2.D0 * delta_x)
        !    gl.cosmo.dDetBdt_mix_numeric(:,a) = dDetBdT_mix_num(:,a)
        !
        !    !numerical derivative in respect to xa of Matrix B
        !    gl.cosmo.dMatrixB_dxa(:,:,a) = (MatrixB_p - MatrixB_m) / (2.D0 * delta_x)
        !    !-------------------------------------------------------------------------------------
        !end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        sum_xiAi = 0.D0
        sum_xiAisigma = 0.D0
        do i = 1, gl.ncomp
            sum_xiAi = sum_xiAi + x(i) * gl.cosmo.Acosmo(i)
        end do
        do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
            do i = 1, gl.ncomp
                sum_xiAisigma(j) = sum_xiAisigma(j) + x(i) * gl.cosmo.sigma(j,i)
            end do
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!test, REMOVE AFTER TESTING
        !sum_xiAi_p = 0.D0
        !sum_xiAi_m = 0.D0
        !sum_xiAisigma_p = 0.D0
        !sum_xiAisigma_m = 0.D0
        !do a = 1, gl.ncomp -1
        !    do i = 1, gl.ncomp
        !        !test
        !        sum_xiAi_p(a) = sum_xiAi_p(a) + x_p(i,a) * gl.cosmo.Acosmo(i)
        !        sum_xiAi_m(a) = sum_xiAi_m(a) + x_m(i,a) * gl.cosmo.Acosmo(i)
        !    end do
        !    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
        !        do i = 1, gl.ncomp
        !            !test
        !            sum_xiAisigma_p(j,a) = sum_xiAisigma_p(j,a) + x_p(i,a) * gl.cosmo.sigma(j,i)
        !            sum_xiAisigma_m(j,a) = sum_xiAisigma_m(j,a) + x_m(i,a) * gl.cosmo.sigma(j,i)
        !        end do
        !    end do
        !end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !help variable Fi_xa
        Fi_xa = 0.D0
        !Fi_xa_p = 0.D0
        !Fi_xa_m = 0.D0
        do a = 1, gl.ncomp-1 !xa
            do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !all Fi_xa's
                !if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                !    cycle
                !end if
                do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !sum for Fi_T
                    !if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                    !    cycle
                    !end if
                    !partial derivatives in respect to xa
                    Fi_xa(j,a) = Fi_xa(j,a) + ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k)) / (sum_xiAi) ** 2.D0 * seggamma(k) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*temp))
                    !Fi_xa(j) = Fi_xa(j) + ((gl.cosmo.sigma(k,i)) * sum_xiAi - (gl.cosmo.Acosmo(i)) * sum_xiAisigma(k)) / (sum_xiAi) ** 2.D0 * seggamma(k) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*temp))

                    !test for second derivative of Fi_xa in respect to Temperature, can be removed after test
                    !Fi_xa_p(j,a) = Fi_xa_p(j,a) + ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi_p(a) - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma_p(k,a)) / sum_xiAi_p(a) ** 2.D0 * seggamma_p(k,a) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp_p))
                    !Fi_xa_m(j,a) = Fi_xa_m(j,a) + ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi_m(a) - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma_m(k,a)) / sum_xiAi_m(a) ** 2.D0 * seggamma_m(k,a) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp_m))
                    !gemische Ableitung dxadT
                    !d2FdTdseggamma_mix(j,k) = (Fi_seggamma_mix_p(j,k) - Fi_seggamma_mix_m(j,k)) / (Temp_p - Temp_m)
                    !gemischte Ableitung dsegammadxa
                    !d2Fdseggammadxa_mix(j,k,a) = (Fi_seggamma_mix_p(j,k,a) - Fi_seggamma_mix_m(j,k,a)) / (2.D0 * delta_x)
                end do
                Fi_xa(j,a) = seggamma(j) * Fi_xa(j,a)
                !test for second derivative of Fi_xa in respect to Temperature, can be removed after test
                !Fi_xa_p(j,a) = seggamma_p(j,a) * Fi_xa_p(j,a)
                !Fi_xa_m(j,a) = seggamma_m(j,a) * Fi_xa_m(j,a)
                !gemische Ableitung dxadT
                !d2FdxadT_mix(j) = (Fi_xa_p(j) - Fi_xa_m(j)) / (Temp_p - Temp_m)
                !gemischte Ableitung dxadxb
                !d2Fdxadxb_mix(j,a) = (Fi_xa_p(j,a) - Fi_xa_m(j,a)) / (2.D0 * delta_x)
            end do

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !test, REMOVE AFTER TESTING
            !gl.cosmo.d2F_dseggammadxa_num(:,:,a) = d2Fdseggammadxa_mix(:,:,a)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end do
        !save for use in second derivative in respect to xa and xb in ge_COSMO_SAC_MIXDERIVS_dxadxb
        gl.cosmo.dFi_dxa = Fi_xa

        !-------------------------------------------------------------------------------------
        pos = 0
        do a = 1, gl.ncomp-1 !xa
            !do i = 1, gl.ncomp
            do o = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !all derivatives of seggamma in respect to xa
                if (abs(gl.cosmo.sigma_profile_mix_gl(o)) < 1.D-14) then
                    cycle
                end if
                !save position of occupied intervalls of sigma profile to match dseggammadxa_mix to the respective intervall
                !as this information is lost using gauss algorithm to solve dseggammadxa_mix
                pos(o) = 1
                m = 0
                do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Zeile MatrixA
                    if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                        cycle
                    end if
                    m = m + 1
                    n = 0
                    do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Spalte (Ableitung) MatrixA
                        if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                            cycle
                        end if
                        n = n + 1
                        if (k == o) then
                            !MatrixA(m,n) = Fi_xa(j,a)
                            MatrixA_mix(j,k,o,a) = Fi_xa(j,a)
                            vectorb_gauss(m) = - Fi_xa(j,a)
                            !nummerical derivatives test
                            !!MatrixA_mix_p(o,j,k) = Fi_T_p(j)
                            !!MatrixA_mix_m(o,j,k) = Fi_T_m(j)
                            !MatrixA_p(m,n) = Fi_xa_p(j)
                            !MatrixA_m(m,n) = Fi_xa_m(j)
                        else
                            !partial derivatives in respect to seggamma
                            !MatrixA(m,n) = dFi_dseggamma(j,k)
                            MatrixA_mix(j,k,o,a) = dFi_dseggamma(j,k)
                            !nummerical derivatives test
                            !!MatrixA_mix_p(o,j,k) = Fi_seggamma_mix_p(j,k)
                            !!MatrixA_mix_m(o,j,k) = Fi_seggamma_mix_m(j,k)
                            !MatrixA_p(m,n) = Fi_seggamma_mix_p(j,k)
                            !MatrixA_m(m,n) = Fi_seggamma_mix_m(j,k)
                        end if
                    end do
                end do

                rankA = n
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !replaced with gauss algorithm
                !vectorb = 0.D0
                !
                !call Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_xa, errval)
                !DetMatrixA_mix(o,a) = Det_xa
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !numerical test
                !gl.cosmo.detMatrixA_mix_num(o,a) = Det_xa
                !!nummerical derivatives test
                !call Gauss_algorithm(gl,MatrixA_p, rankA, vectorb, vectorx, Det_xa, errval)
                !DetMatrixA_mix_p(o,a) = Det_xa
                !!nummerical derivatives test
                !call Gauss_algorithm(gl,MatrixA_m, rankA, vectorb, vectorx, Det_xa, errval)
                !DetMatrixA_mix_m(o,a) = Det_xa

                !nummerical derivatives test
                !dDetAdT_mix_num(o) = (DetMatrixA_mix_p(o,a) - DetMatrixA_mix_m(o,a)) / (Temp_p - Temp_m)

                !gl.cosmo.dDetAdxb_mix_numeric(o,a) = (DetMatrixA_mix_p(o,a) - DetMatrixA_mix_m(o,a)) / (2.D0 * delta_x)

                !dseggammadxa_mix_num(o,a) = - DetMatrixA_mix_num(o,a) / DetMatrixB_mix_num

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !replaced with gauss algorithm
                !calculate dseggammadxa_mix
                !dseggammadxa_mix(o,a) = - DetMatrixA_mix(o,a) / DetMatrixB_mix
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !!numerical test
                !dseggammadxa_mix_p(o,a) = - DetMatrixA_mix_p(o,a) / DetMatrixB_mix_p
                !dseggammadxa_mix_m(o,a) = - DetMatrixA_mix_m(o,a) / DetMatrixB_mix_m
                !!d2seggammadxadT_mix_num(o,a) = (dseggammadxa_mix_p(o,a) - dseggammadxa_mix_m(o,a)) / (Temp_p - Temp_m)
                !gl.cosmo.d2seggamma_dxa_dxb_numeric(o,a) = (dseggammadxa_mix_p(o,a) - dseggammadxa_mix_m(o,a)) / (2.D0 * delta_x)

                !!numerical derivatives of Matrix A in respect to xa
                !do j = 1, 60
                !    do k = 1, 60
                !        gl.cosmo.dMatrixA_dxa_num(j,k,o,a) = (MatrixA_p(j,k) - MatrixA_m(j,k)) / (2.D0 * delta_x)
                !    end do
                !end do

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !test, remove after test
                !gl.cosmo.MatrixA(:,:,o,a) = MatrixA(:,:)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            end do
            !Gauss Algorithm is faster than Cramer's rule to calculate derivatives of system of equations
            call Gauss_algorithm(gl,MatrixB, rankA, vectorb_gauss, vectorx, Det_xa, errval)
            dseggammadxa_mix_gauss(:,a) = vectorx(1:51)

            m = 0
            do o = 1, 51
                If (pos(o) == 0) then
                else
                    m = m + 1
                    dseggammadxa_mix(o,a) = dseggammadxa_mix_gauss(m,a)
                end if
            end do
        end do

        !for Derivative wrt tau
        DetMatrixA_mix = - dseggammadxa_mix * DetMatrixB_mix

        !save for second derivative of help function Fi in respect to xa and xb
        gl.cosmo.dseggamma_dxa_mix = dseggammadxa_mix

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !test, REMOVE AFTER TEST
        ! dseggammadxa_mix = dseggamma_dxa
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Derivations are dy/dx, need to be transformed to dlny/dx; dlny/dx = 1/y dy/dx
        do a = 1, gl.ncomp - 1
            do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                !dlnseggammadxa_pure(j,a) = 0.D0
                dlnseggammadxa_mix(j,a) = dseggammadxa_mix(j,a) / seggamma(j)

                !numerical test
                !dlnseggammadxa_mix_p(j,a) = dseggammadxa_mix_p(j,a) / seggamma_p(j)
                !dlnseggammadxa_mix_m(j,a) = dseggammadxa_mix_m(j,a) / seggamma_m(j)
                !!d2lnseggammadxadT_mix_num(j,a) = (dlnseggammadxa_mix_p(j,a) - dlnseggammadxa_mix_m(j,a)) / (Temp_p - Temp_m)
                !d2lnseggammadxadT_mix_num(j,a) = (dlnseggammadxa_mix_p(j,a) - dlnseggammadxa_mix_m(j,a)) / (2.D0 * delta_x)
                !gl.cosmo.d2lnseggammadxadT_mix_numeric(j,a) = d2lnseggammadxadT_mix_num(j,a)

            end do
        end do

        !test
        !ln_gamma_R_Trho_dxa_p = 0.D0
        !ln_gamma_R_Trho_dxa_m = 0.D0
        !ln_gamma_R_Trho_dxa_num = 0.D0
        do a = 1, gl.ncomp - 1    !xa
            !do k = 1, gl.ncomp - 1
            do i = 1, gl.ncomp  !component
                do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    !!!NOTE THAT segment gammas OF PURE COMPONENT i IS NOT A FUNCTION OF COMPOSITION, HENCE dlnseggammadxa_pure(j,i) = 0!!! This is indicated with the 0.D0 in the following equation
                    ln_gamma_R_Trho_dxa(1,i,a) = ln_gamma_R_Trho_dxa(1,i,a) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (dlnseggammadxa_mix(j,a) - 0.D0)
                    !!numerical test
                    !ln_gamma_R_Trho_dxa_p(1,i,a) = ln_gamma_R_Trho_dxa_p(1,i,a) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (dlnseggammadxa_mix_p(j,a) - 0.D0)
                    !ln_gamma_R_Trho_dxa_m(1,i,a) = ln_gamma_R_Trho_dxa_m(1,i,a) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (dlnseggammadxa_mix_m(j,a) - 0.D0)

                end do
                !numerical test
                !!ln_gamma_R_Trho_dxa_num(4,i,a) = (ln_gamma_R_Trho_dxa_p(1,i,a) - ln_gamma_R_Trho_dxa_m(1,i,a)) / (Temp_p - Temp_m)
                !ln_gamma_R_Trho_dxa_num(4,i,a) = (ln_gamma_R_Trho_dxa_p(1,i,a) - ln_gamma_R_Trho_dxa_m(1,i,a)) / (2.D0 * delta_x)
                !gl.cosmo.ln_gamma_R_Trho_dxa_numeric(4,i,a) = ln_gamma_R_Trho_dxa_num(4,i,a)
            end do
        end do

        !Transformation of variables from T,rho,x, to tau,delta, and x
        do a = 1, gl.ncomp - 1
            do i = 1, gl.ncomp
                gl%ge%ln_gamma_R_dxa(1,i,a) =  dTred_dxa(a) / tau * ln_gamma_R_Trho(4,i) + ln_gamma_R_Trho_dxa(1,i,a)
            end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !test, REMOVE AFTER TESTING
        gl.cosmo.ln_gamma_R_Trho_dxa_numeric = ln_gamma_R_Trho_dxa
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Note in the following that all derivatives with respect to delta at constant tau and x are 0, because the residual part is not a function of density
        if ((GETDER(1) .eq. 1) .or. (GETDER(4) .eq. 1) .or. (GETDER(5) .eq. 1)) then

            if ((gl%ge%GETDER_dxa_prev(1) == 1) .and. (recalc .eqv. .false.)) then

                gl%ge%gE_R_dxa(1,:) = gl%ge%gE_R_dxa_prev(1,:)
                !Transformation of variables from tau,delta,x to T,rho,x
                do a=1,gl.ncomp-1
                    gE_R_Trho_dxa(1,a) = gl%ge%gE_R_dxa(1,a) - dTred_dxa(a) / tau * gE_R_Trho(4)
                end do

            else

                Do a = 1, gl.ncomp - 1
                    Do i = 1, gl.ncomp
                        gE_R_Trho_dxa(1,a) = gE_R_Trho_dxa(1,a) + R_const * Temp * x(i) * ln_gamma_R_Trho_dxa(1,i,a)
                    end do
                    gE_R_Trho_dxa(1,a) = gE_R_Trho_dxa(1,a) + R_const * Temp * (ln_gamma_R_Trho(1,a) - ln_gamma_R_Trho(1,gl.ncomp))
                End do

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !test, REMOVE AFTER TESTING
                !gl.cosmo.gE_R_Trho_dxa_num = gE_R_Trho_dxa
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !Calculate the derivative of gE_R with respect to temperature T at constant rho and x. Also needed for derivatives with respect to tau and xa
                !Transformation of variables from T,rho,x, to tau,delta, and x
                do a = 1, gl.ncomp - 1
                    gl%ge%gE_R_dxa(1,a) =  dTred_dxa(a) / tau * gE_R_Trho(4) + gE_R_Trho_dxa(1,a)
                end do

            end if

        end if


        !Second derivative wrt delta and tau
        if (GETDER(2) .eq. 1) then
            ln_gamma_R_Trho_dxa(2,:,:) = 0.D0
            gE_R_Trho_dxa(2,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(2,:,:) = 0.d0
            gl%ge%gE_R_dxa(2,:) = 0.D0
        end if

        !Third derivative wrt delta, tau, and tau
        if (GETDER(3) .eq. 1) then
            ln_gamma_R_Trho_dxa(3,:,:) = 0.D0
            gE_R_Trho_dxa(3,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(3,:,:) = 0.d0
            gl%ge%gE_R_dxa(3,:) = 0.D0
        end if

        !Derivative wrt tau, times tau
        if ((GETDER(4) .eq. 1) .or. (GETDER(5) .eq. 1)) then
            !---------------------------------------------------------------------------------------
            !Calculate the first temperature derivative of all segment activity coefficients
            !Do this always because it might be needed for the following derivatives, too
            !---------------------------------------------------------------------------------------

            !!set all COSMO-SAC values to orig
            !can be removed after successful testing of analytical derivs
            !call gE_COSMO_SAC_CALC(gl, temp, C_or_R, errval)

            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            !calculate Derivative of Determinant of Matrix B
            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            !Matrix dB/dT, second partial derivatives of B in respect to T
            m = 0
            do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                n = 0
                if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                    cycle
                end if
                m = m + 1
                do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                        cycle
                    end if
                    n = n + 1
                    dMatrixBdT(m,n) = gl.cosmo.d2Fi_dseggamma_dT_mix(j,k)
                end do
            end do

            !Matrix B is independent of derivative in respect to xa.
            !Note: needs MatrixB_mix from first derivative in respect to xa
            MatrixB = 0.D0
            m = 0
            do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                n = 0
                if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                    cycle
                end if
                m = m + 1
                do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                        cycle
                    end if
                    n = n + 1
                    MatrixB(m,n) = MatrixB_mix(j,k)
                end do
            end do

            rankA = m

            call Adjugate(gl,MatrixB, rankA, adj_B, errval)
            call Mat_mult(gl,adj_B, rankA, rankA, dMatrixBdT, rankA, rankA, MatrixD, errval)
            call Trace(gl,MatrixD, rankA, dDetBdT_mix, errval)

            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            !calculate Derivative of Determinant of Matrix A
            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            !sum_d2Fi_dT2 = 0.D0
            !sum2_d2Fi_dT2 = 0.D0
            !sum_d2Fi_dseggammadT = 0.D0
            !d2Fi_dT2_pure_test = 0.D0
            !sum3_d2Fi_dT2 = 0.D0
            sum_xiAi = 0.D0
            sum_xiAisigma = 0.D0
            do i = 1, gl.ncomp
                sum_xiAi = sum_xiAi + x(i) * gl.cosmo.Acosmo(i)
            end do

            do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                do i = 1, gl.ncomp
                    sum_xiAisigma(j) = sum_xiAisigma(j) + x(i) * gl.cosmo.sigma(j,i)
                end do
            end do
            do a = 1, gl.ncomp - 1  !for all first derivatives in respect to xa
                !d2Fi_dxa_dT_mix for all derivatives of Fi
                d2Fi_dxa_dT_mix = 0.D0
                do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    !if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                    !    cycle
                    !end if
                    do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                        !    cycle
                        !end if
                        !wrong derivation
                        !d2Fi_dxa_dT_mix(j) = d2Fi_dxa_dT_mix(j) + ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k)) / (sum_xiAi) ** 2.D0 * seggamma(k) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*temp)) *(gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0))
                        d2Fi_dxa_dT_mix(j) = d2Fi_dxa_dT_mix(j) + (gl.cosmo.dseggamma_dT_mix(j) * seggamma(k) + seggamma(j) * seggamma(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)) + seggamma(j) * gl.cosmo.dseggamma_dT_mix(k)) * ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*temp))
                        !Sum for d2Fi_dseggammadT
                        !sum_d2Fi_dseggammadT(j) = sum_d2Fi_dseggammadT(j) + gl.cosmo.sigma_profile_mix_gl(k) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggammadt_mix(k) + seggamma(k) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)))

                    end do
                    !kann vielleicht aus der Mixderivs genommen werden, zumindest kann es vor der a-Schleife platziert werden
                    !do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    !    if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                    !        cycle
                    !    end if
                    !    if (j == k) then
                    !        d2Fi_dseggammadT_mix(j,k) = sum_d2Fi_dseggammadT(j) + gl.cosmo.sigma_profile_mix_gl(j) * exp(-gl.cosmo.delta_w_gl(j,j) / (R_const_cal*Temp)) * (dseggammadt_mix(j) + seggamma(j) * (gl.cosmo.delta_w_gl(j,j) / (R_const_cal*Temp ** 2.D0)))
                    !    else
                    !        d2Fi_dseggammadT_mix(j,k) = gl.cosmo.sigma_profile_mix_gl(k) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggammadt_mix(j) + seggamma(j) * (gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp ** 2.D0)))
                    !    end if
                    !end do
                end do
                d2Fi_dxa_dT_mix = d2Fi_dxa_dT_mix / (sum_xiAi) ** 2.D0
                !d2Fi_dxa_dT_mix = d2FdxadT_mix

                !matrix dA/dT, partial derivatives of A in respect to T
                do o = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)   !all derivatuves of seggamma
                    if (abs(gl.cosmo.sigma_profile_mix_gl(o)) < 1.D-14) then
                        cycle
                    end if
                    m = 0
                    dMatrixAdT = dMatrixBdT
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Zeile
                        if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        n = 0
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Spalte (Ableitung)
                            if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            if (k == o) then
                                dMatrixAdT(m,n) = d2Fi_dxa_dT_mix(j)
                            else
                            end if
                        end do
                    end do

                    !Adjugates of MatrixA
                    MatrixA = 0.D0
                    m = 0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        n = 0
                        if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            MatrixA(m,n) = MatrixA_mix(j,k,o,a)
                        end do
                    end do

                    !m = 0
                    !do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Zeile MatrixA
                    !    if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                    !        cycle
                    !    end if
                    !    m = m + 1
                    !    n = 0
                    !    do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Spalte (Ableitung) MatrixA
                    !        if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                    !            cycle
                    !        end if
                    !        n = n + 1
                    !        if (k == o) then
                    !            MatrixA(m,n) = Fi_xa(j,a)
                    !        else
                    !            !partial derivatives in respect to seggamma
                    !            MatrixA(m,n) = dFi_dseggamma(j,k)
                    !        end if
                    !    end do
                    !end do
                    rankA = m

                    call Adjugate(gl,MatrixA, rankA, adj_A, errval)
                    !kann außerhalb der Schleife berechnet werden
                    !call Adjugate(gl,MatrixB, rankA, adj_B, errval)

                    !Multiply adj(MatrixA) * dMatrixAdT
                    call Mat_mult(gl,adj_A, rankA, rankA, dMatrixAdT, rankA, rankA, MatrixC, errval)
                    !vor der a-Schleife
                    !call Mat_mult(gl,adj_B, rankA, rankA, dMatrixBdT, rankA, rankA, MatrixD, errval)

                    !calculate dDetA/dT = Trace(Adjugate(A) dMatrixA/dT)
                    call Trace(gl,MatrixC, rankA, dDetAdT_mix(o,a), errval)
                    !vor der a-Schleife
                    !call Trace(gl,MatrixD, rankA, dDetBdT_mix(o), errval)

                end do
            end do

            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            !Calculate second derivative of segment gamma wrt T for mixture
            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            !Note detMatrixA_mix and detMatrixB_mix are taken from first derivative of gE_R in respect to xa

            do a = 1, gl.ncomp - 1
                do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    d2seggamma_dxa_dT_mix(j,a) = - ( dDetAdT_mix(j,a) * detMatrixB_mix - dDetBdT_mix * detMatrixA_mix(j,a)) / (detMatrixB_mix ** 2.D0)
                end do
            end do

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !test
            !gl.cosmo.d2seggamma_dxa_dT_mix_test = d2seggamma_dxa_dT_mix
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            !Derivatives need to be for tau not Temp
            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            do a = 1, gl.ncomp - 1
                do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    !if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                    !    cycle
                    !end if
                    d2lnseggamma_dxa_dT_mix(j,a) =  - 1.D0 / seggamma(j) ** 2.D0 * gl.cosmo.dseggamma_dT_mix(j) * dseggammadxa_mix(j,a) + d2seggamma_dxa_dT_mix(j,a) / seggamma(j)
                end do
            end do

            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            !Calculate the derivative of the logarithm of the residual activity coefficients with respect to xa and temperature
            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            !!!NOTE THAT LNGAMMA_ki OF PURE COMPONENT i IS NOT A FUNCTION OF COMPOSITION, HENCE d2lnseggamma_dxa_dT_pure(j,a,i) = 0!!! This is indicated with the 0.D0 in the following equation
            do a = 1, gl.ncomp - 1
                do i = 1, gl.ncomp
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                        !    cycle
                        !end if
                        ln_gamma_R_Trho_dxa(4,i,a) = ln_gamma_R_Trho_dxa(4,i,a) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (d2lnseggamma_dxa_dT_mix(j,a) - 0.D0)
                    end do
                end do
            end do
            !ln_gamma_R_Trho_dxa(4,:,:) = ln_gamma_R_Trho_dxa_num(4,:,:)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !test
            !gl.cosmo.ln_gamma_R_Trho_dxa_test = ln_gamma_R_Trho_dxa
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            !Compute the derivative of gE_R with respect to xa and T (natural variables)
            Do a = 1,gl.ncomp-1
                Do i = 1, gl.ncomp
                    gE_R_Trho_dxa(4,a) = gE_R_Trho_dxa(4,a) + R_const * x(i) * ln_gamma_R_Trho_dxa(1,i,a) &
                        & +  R_const * Temp * x(i) * ln_gamma_R_Trho_dxa(4,i,a)
                end do
                gE_R_Trho_dxa(4,a) = gE_R_Trho_dxa(4,a) + R_const * (ln_gamma_R_Trho(1,a) - ln_gamma_R_Trho(1,gl.ncomp)) &
                    &  + R_const * Temp * (ln_gamma_R_Trho(4,a) - ln_gamma_R_Trho(4,gl.ncomp))
            End do

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !test
            gl.cosmo.gE_R_Trho_dxa_test = gE_R_Trho_dxa
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Transformation of variables from T,rho,x, to tau,delta, and x
            do a=1,gl.ncomp-1
                Do i = 1, gl.ncomp
                    gl%ge%ln_gamma_R_dxa(4,i,a) = - dTred_dxa(a) / tau**2 * ln_gamma_R_Trho(4,i) &
                        & - gl.tredmix / tau**3 * dTred_dxa(a) * ln_gamma_R_Trho(5,i) &
                        & - gl.tredmix / tau**2 * ln_gamma_R_Trho_dxa(4,i,a)
                end do
                gl%ge%gE_R_dxa(4,a) = - dTred_dxa(a) / tau**2 * gE_R_Trho(4) &
                    & - gl.tredmix / tau**3 * dTred_dxa(a) * gE_R_Trho(5) &
                    & - gl.tredmix / tau**2 * gE_R_Trho_dxa(4,a)
            end do
            !
            !end if

        end if


        !Second derivative wrt tau NOT YET IMPLEMENTED
        if (GETDER(5) .eq. 1) then
            ln_gamma_R_Trho_dxa(5,:,:) = 0.D0
            gE_R_Trho_dxa(5,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(5,:,:) = 0.d0
            gl%ge%gE_R_dxa(5,:) = 0.D0
        end if
        !Second derivative wrt delta and tau
        if (GETDER(6) .eq. 1) then
            ln_gamma_R_Trho_dxa(6,:,:) = 0.D0
            gE_R_Trho_dxa(6,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(6,:,:) = 0.d0
            gl%ge%gE_R_dxa(6,:) = 0.D0
        end if
        !Third derivative wrt delta, tau, and tau
        if (GETDER(7) .eq. 1) then
            ln_gamma_R_Trho_dxa(7,:,:) = 0.D0
            gE_R_Trho_dxa(7,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(7,:,:) = 0.d0
            gl%ge%gE_R_dxa(7,:) = 0.D0
        end if
        !Third derivative wrt delta
        if (GETDER(8) .eq. 1) then
            ln_gamma_R_Trho_dxa(8,:,:) = 0.D0
            gE_R_Trho_dxa(8,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(8,:,:) = 0.d0
            gl%ge%gE_R_dxa(8,:) = 0.D0
        end if
        !Third derivative wrt tau
        if (GETDER(9) .eq. 1) then
            ln_gamma_R_Trho_dxa(9,:,:) = 0.D0
            gE_R_Trho_dxa(9,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(9,:,:) = 0.d0
            gl%ge%gE_R_dxa(9,:) = 0.D0
        end if
        !Third derivative wrt tau, delta, and delta
        if (GETDER(10) .eq. 1) then
            ln_gamma_R_Trho_dxa(10,:,:) = 0.D0
            gE_R_Trho_dxa(10,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(10,:,:) = 0.d0
            gl%ge%gE_R_dxa(10,:) = 0.D0
        end if
        !Fourth derivative wrt delta
        if (GETDER(11) .eq. 1) then
            ln_gamma_R_Trho_dxa(11,:,:) = 0.D0
            gE_R_Trho_dxa(11,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(11,:,:) = 0.d0
            gl%ge%gE_R_dxa(11,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, delta and tau
        if (GETDER(12) .eq. 1) then
            ln_gamma_R_Trho_dxa(12,:,:) = 0.D0
            gE_R_Trho_dxa(12,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(12,:,:) = 0.d0
            gl%ge%gE_R_dxa(12,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, tau, tau
        if (GETDER(13) .eq. 1) then
            ln_gamma_R_Trho_dxa(13,:,:) = 0.D0
            gE_R_Trho_dxa(13,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(13,:,:) = 0.d0
            gl%ge%gE_R_dxa(13,:) = 0.D0
        end if
        !Fourth derivative wrt delta, tau, tau, tau
        if (GETDER(14) .eq. 1) then
            ln_gamma_R_Trho_dxa(14,:,:) = 0.D0
            gE_R_Trho_dxa(14,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(14,:,:) = 0.d0
            gl%ge%gE_R_dxa(14,:) = 0.D0
        end if
        !Fourth derivative wrt tau
        if (GETDER(15) .eq. 1) then
            ln_gamma_R_Trho_dxa(15,:,:) = 0.D0
            gE_R_Trho_dxa(15,:) = 0.d0
            gl%ge%ln_gamma_R_dxa(15,:,:) = 0.d0
            gl%ge%gE_R_dxa(15,:) = 0.D0
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !End if
    end if

    end subroutine

    !Erik, August, 2018
    !Derivatives of COSMO-SAC with respect to mole frations xa and xb
    !Needed for phase equilibrium routines for PSRK(gl,and VTPR and new Helmholtz+gE) model
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine gE_COSMO_SAC_MIXDERIVS_dxadxb(gl,Temp, GETDER, C_or_R, errval)
    !------------------------------------------------------------------------------------------------------------------------------
    !Subroutine to calculate the excess Gibbs energy and the activity coefficients of the
    !components in the mixture according to the COSMO-SAC model as well as tau- and delta-
    !derivatives of the model.
    !All quantities are given as derivatives with respect to mole fractions xa and xb for all components
    !in the mixture. The convention that the "last" mole fraction xN of component N is replaced with
    !1-sum(xk) is used here.
    !
    !INPUT
    !Temp       -    Temperature in K
    !Note that the model makes use of the global variable "molfractions" defined in the module
    !"module_fluid_parameters"
    ! GETDER      - AN ARRAY WITH nderivs ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED. !!!ALL DERIVATIVES ARE ADDITIONALLY DERIVATIVES WITH RESPECT TO xa and xb!!!:
    !                1. gE from UNIFAC as a function of temperature and composition
    !                2. 1ST DERIVATIVE OF gE WITH RESPECT TO DEL AT CONSTANT TAU
    !                3. 2ND DERIVATIVE OF gE WITH RESPECT TO DEL AT CONSTANT TAU
    !                4. 1ST DERIVATIVE OF gE WITH RESPECT TO TAU AT CONSTANT DEL
    !                5: 2ND DERIVATIVE OF gE WITH RESPECT TO TAU AT CONSTANT DEL
    !                6: 2ND MIXED DERIVATIVE OF gE WITH RESPECT TO DEL AND TAU
    !                7: 3RD MIXED DERIVATIVE OF gE WITH RESPECT TO DEL, TAU, AND TAU
    !                8: 3RD DERIVATIVE OF gE WITH RESPECT TO DEL
    !                9: 3RD DERIVATIVE OF gE WITH RESPECT TO TAU
    !               10: 3RD MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, DEL, AND DEL
    !               11: 4TH DERIVATIVE OF gE WITH RESPECT TO DEL
    !               12: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, DEL, DEL AND DEL
    !               13: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, TAU, DEL, AND DEL
    !               14: 4TH MIXED DERIVATIVE OF gE WITH RESPECT TO TAU, TAU, TAU AND DEL
    !               15: 4TH DERIVATIVE OF gE WITH RESPECT TO TAU
    !C_or_R       - Integer that specifies if the combinatorial, the residual or both parts shall be calculated
    !               C_or_R = 0: Calculate both parts
    !               C_or_R = 1: Calculate the combinatorial part only
    !               C_or_R = 2: Calculate the the residual part only
    !
    !OUTPUT
    !gE_C_dxadxb         -    Matrix with molar excess Gibbs energy of the combinatorial part (C) in J/molK and derivatives with respect to tau and delta and xa
    !gE_R_dxadxb         -    Matrix with molar excess Gibbs energy of the residual part (R) in J/molK and derivatives with respect to tau and delta and xa
    !ln_gamma_C_dxadxb   -    Matrix containing the natural logarithm of the activity coefficients of the combinatorial part for each component and derivatives with respect to tau and delta and xa
    !ln_gamma_R_dxadxb   -    Matrix containing the natural logarithm of the activity coefficients of the residual part for each component and derivatives with respect to tau and delta and xa
    !errval              -    Indicates if an error occurred during calculations
    !
    !Andreas Jäger, April 2017
    !------------------------------------------------------------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl


    double precision :: Temp
    integer, dimension(nderivs):: GETDER                                    !array specifier to indicate, which derivative is needed
    !double precision, allocatable :: gl.gE_C_dxadxb(:,:,:)                 !Combinatorial part of gE and derivatives with respect to tau and delta and xa and xb
    !double precision, allocatable :: gl.gE_R_dxadxb(:,:,:)                 !Residual part of gE and derivatives with respect to tau and delta and xa and xb
    !double precision, allocatable :: gl.ln_gamma_C_dxadxb(:,:,:,:)     !Combinatorial activity coefficients and derivatives with respect to delta and tau and xa and xb
    !double precision, allocatable :: gl.ln_gamma_R_dxadxb(:,:,:,:)     !Residual activity coefficients and derivatives with respect to delta and tau and xa and xb
    integer:: C_or_R
    integer:: errval

    !double precision, dimension(30) :: r_i              !Molecular vdW volume of components i
    !double precision, dimension(30) :: q_i              !Molecular vdW surface area of components i

    double precision :: R_const, R_const_cal            !universal gas constant

    !Help variables
    double precision, dimension(30) :: x
    !double precision, dimension(30) :: phi_i
    !double precision, dimension(30) :: Theta_i
    !!double precision, dimension(100) :: lnGamma_k
    !!double precision, dimension(100,30) :: lnGamma_ki   !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    !double precision, dimension(100) :: Theta_m
    !double precision, dimension(100,30) :: Theta_mi     !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    !double precision, dimension(100) :: X_m             !Mole fraction of groups in the mixture
    !double precision, dimension(100,30) :: X_mi         !Mole fraction of groups in pure component i. !ADDED INDEX FOR PURE FLUID i. VARIABLE IS NEEDED FOR CALCULATION OF DERIVATIVES, DOES NOT CHANGE WITH COMPOSITION AND IS THUS SAFED. Andreas Jäger, April 2017
    !double precision, dimension(100,100) :: Psi_nm
    !double precision, dimension(100):: Q_m              !Help variable to list all group areas in one vector in order to not have to map them back to Q_ik(i, k) for derivative calculation

    !More help variables
    !double precision:: sum_rjxj
    !double precision:: sum_qjxj
    !double precision:: sum_QnXn
    !double precision, dimension(30):: sum_QndXn_dxa
    !double precision, dimension(30,30):: sum_QndXn_dxadxb
    !double precision:: sumsum_vnxj
    !!double precision:: sum_xjLj
    !double precision:: sum_ThetamPsimk
    !!double precision:: sumsum_ThetamPsikm
    !double precision:: sum_ThetanPsinm

    !Derivatives of the combinatorial part with respect to "natural" variables
    double precision, dimension(nderivs, 30) :: ln_gamma_C_Trho       !Combinatorial activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs, 30) :: ln_gamma_R_Trho       !Residual activity coefficients and derivatives with respect to rho and T
    double precision, dimension(nderivs) :: gE_C_Trho        !Combinatorial part of gE and derivatives with respect to T and rho
    double precision, dimension(nderivs) :: gE_R_Trho        !Residual part of gE and derivatives with respect to T and rho
    double precision, allocatable  :: ln_gamma_C_Trho_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to rho and T, and xa
    double precision, allocatable  :: ln_gamma_R_Trho_dxa(:,:,:)       !Residual activity coefficients and derivatives with respect to rho and T, and xa
    double precision, dimension(nderivs, 30) :: gE_C_Trho_dxa        !Combinatorial part of gE and derivatives with respect to T and rho and xa
    double precision, dimension(nderivs, 30) :: gE_R_Trho_dxa        !Residual part of gE and derivatives with respect to T and rho and xa
    double precision, allocatable :: gE_C_Trho_dxadxb(:,:,:)                 !Combinatorial part of gE and derivatives with respect to T and rho and xa and xb
    double precision, allocatable :: gE_R_Trho_dxadxb(:,:,:)                 !Residual part of gE and derivatives with respect to T and rho and xa and xb
    double precision, allocatable :: ln_gamma_C_Trho_dxadxb(:,:,:,:)       !Combinatorial activity coefficients and derivatives with respect to rho and T and xa and xb
    double precision, allocatable :: ln_gamma_R_Trho_dxadxb(:,:,:,:)       !Residual activity coefficients and derivatives with respect to rho and T and xa and xb

    !Derivatives of the residual part with respect to "natural" variables
    !For combinatorial part
    !double precision, dimension(30,30):: dphi_i_dxa
    !double precision, dimension(30,30):: dTheta_i_dxa
    !double precision, allocatable:: d2phi_i_dxadxb(:,:,:)
    !double precision, allocatable:: d2Theta_i_dxadxb(:,:,:)
    !!For residual part
    !!double precision, dimension(100,30):: dlnGamma_k_dxa
    !!double precision, dimension(100,30):: d2lnGamma_k_dxadT
    !double precision, allocatable:: d2lnGamma_k_dxadxb(:,:,:)
    !!double precision, dimension(100,100) :: dPsi_nm_dT
    !double precision, dimension(100,30):: dTheta_m_dxa
    !double precision, allocatable:: d2Theta_m_dxadxb(:,:,:)
    !double precision, dimension(100,30):: dX_m_dxa
    !double precision, allocatable:: d2X_m_dxadxb(:,:,:)
    !double precision:: sum_dThetam_dxa_Psimk
    !double precision:: sum_dThetam_dxb_Psimk
    !double precision:: sum_d2Thetam_dxadxb_Psimk
    !double precision:: sum_dThetan_dxa_Psinm
    !double precision:: sum_dThetan_dxb_Psinm
    !double precision:: sum_d2Thetan_dxadxb_Psinm
    !double precision:: sumsum_ThetamPsikm_dxadxb          !Variable for double sum in xa and xb derivative
    !double precision:: sum_Thetam_dPsimk_dT
    !double precision:: sum_Thetan_dPsinm_dT
    !double precision:: sum_dThetam_dxa_dPsimk_dT
    !double precision:: sum_dThetan_dxa_dPsinm_dT
    !double precision:: sumsum_ThetamPsikm_dxadT        !Variable for double sum in xa and temperature derivative

    double precision, dimension(30):: dTred_dxa, drhored_dxa            !Variables for the derivatives of the reducing functions with respect to molfractions xa
    double precision, dimension(30):: d2Tred_dxa2, drhored_dxa2         !Variables for the second derivatives of the reducing functions with respect to molfractions xa
    double precision, dimension(30,30):: d2Tred_dxadxb, drhored_dxadxb  !Variables for the second derivatives of the reducing functions with respect to molfractions xa and xb

    !Required derivatives of activities and excess Gibbs energy only with respect to tau and delta (not xa and xb, however, needed for the calculation of derivatives with respect to xa)
    integer, dimension(nderivs):: GETDER_no_xaxb        !array specifier to indicate, which derivative is needed
    !double precision, dimension(nderivs) :: gl.gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs) :: gl.gE_R        !Residual part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs, 30) :: gl.ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, dimension(nderivs, 30) :: gl.ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau

    !Required derivatives of activities and excess Gibbs energy only with respect to tau and delta and xa (not xb, however, needed for the calculation of derivatives with respect to xa and xb)
    integer, dimension(nderivs):: GETDER_no_xb                      !array specifier to indicate, which derivative is needed
    !double precision, dimension(nderivs,30) :: gl.gE_C_dxa                !Combinatorial part of gE and derivatives with respect to tau and delta
    !double precision, dimension(nderivs,30) :: gl.gE_R_dxa                !Residual part of gE and derivatives with respect to tau and delta
    !double precision, allocatable:: gl.ln_gamma_C_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to delta and tau
    !double precision, allocatable:: gl.ln_gamma_R_dxa(:,:,:)      !Residual activity coefficients and derivatives with respect to delta and tau

    !COSMO-SAC variables
    double precision, dimension(51) :: sigma_profile_mix, seggamma
    !double precision, dimension(51,51) :: dFi_dseggamma
    double precision, dimension(:,:), allocatable :: dFi_dseggamma
    double precision, dimension(51,30) :: dseggamma_dxa, dFi_dxa


    !help variables
    double precision :: sum_xiAi
    double precision, dimension(51) :: sum_xiAisigma, sum_d2F_dxadxb, sum_d2F_dxadxb2
    !double precision, dimension(51,30,30) :: d2F_dxadxb!, d2F_dxadxb2    !d2F_dxadxb2 kann nach test gelöscht werden
    double precision, dimension(:,:,:), allocatable :: d2F_dxadxb
    double precision, dimension(51) :: sum_d2Fi_dseggammadxb
    !double precision, dimension(51,51) :: d2Fi_dseggammadxb_mix
    double precision, dimension(:,:), allocatable :: d2Fi_dseggammadxb_mix
    !double precision, dimension(60,60) :: dMatrixA_dxb, dMatrixB_dxb, MatrixA, MatrixB, adj_A, adj_B, MatrixC, MatrixD
    double precision, dimension(:,:), allocatable:: dMatrixA_dxb
    double precision, dimension(:,:), allocatable:: dMatrixB_dxb
    double precision, dimension(:,:), allocatable:: MatrixA
    double precision, dimension(:,:), allocatable:: MatrixB
    double precision, dimension(:,:), allocatable:: MatrixC
    double precision, dimension(:,:), allocatable:: MatrixD
    double precision, dimension(:,:), allocatable:: adj_A
    double precision, dimension(:,:), allocatable:: adj_B

    double precision, dimension(51) :: dDetA_dxb_mix
    double precision :: dDetB_dxb_mix
    !double precision, dimension(51,30,30) :: d2seggamma_dxadxb_mix, dlnseggamma_dxadxb_mix
    double precision, dimension(:,:,:), allocatable:: d2seggamma_dxadxb_mix
    double precision, dimension(:,:,:), allocatable:: dlnseggamma_dxadxb_mix

    double precision, dimension(51,30) :: detMatrixA_mix
    double precision :: DetMatrixB_mix, Det_seggamma, Det_xa
    !November 2018, Erik, allocatable
    !double precision, dimension(60,60,51,30) :: MatrixA_mix
    double precision, dimension(:,:,:,:), allocatable :: MatrixA_mix
    !double precision, dimension(51,51) :: MatrixB_mix
    double precision, dimension(:,:), allocatable :: MatrixB_mix

    !Gauss algorithm
    integer:: rankA
    double precision, dimension(60):: vectorb
    double precision, dimension(60):: vectorx

    !Corresponding states relevant variables
    double precision:: tau, del

    !Summation variables
    integer:: a, b, i, j, k, m, n, o, help_k

    !Variable that indicates if derivatives need to be recalculated
    logical:: recalc

    !--------------------------------------------------------------------------------------------------
    !test variables for numerical derivatives, remove after testing
    !double precision :: delta_x, Temp_p, Temp_m, tredmix_orig, rhoredmix_orig, Temperature
    !double precision, dimension(30) :: molfractions_orig
    !double precision, dimension(51) :: seggamma_p, seggamma_m, seggamma_dxa
    !integer, dimension(nderivs):: getder_dep
    !double precision, dimension(51,30) :: dseggamma_dxa_m, dseggamma_dxa_p, dFi_dxa_p, dFi_dxa_m, detMatrixA_mix_num_p, detMatrixA_mix_num_m
    !double precision, dimension(51,30,30) :: d2segamma_dxadxb_num, d2Fi_dxadxb_num, ddetMatrixA_dxb_mix_num
    !double precision, dimension(15,30,30) :: ln_gamma_R_dxa_p, ln_gamma_R_dxa_m, gE_R_Trho_dxadxb_num, gE_R_dxadxb_num
    !double precision, dimension(15,30,30,30) :: ln_gamma_R_dxadxb_num
    !double precision, dimension(15,30) :: gE_R_Trho_dxa_p, gE_R_Trho_dxa_m, gE_R_dxa_p, gE_R_dxa_m
    !double precision, dimension(51) :: sum_d2F_dxadxb3, sum_d2F_dxadxb4, sum_d2F_dxadxb5
    !double precision, dimension(51,51,30) :: Fi_seggamma_mix_p, Fi_seggamma_mix_m
    !double precision, dimension(51,51,30) :: d2F_dseggammadxb_num
    !double precision, dimension(60,60,51,30) :: MatrixA_p, MatrixA_m
    !double precision, dimension(60,60,51,30,30) :: dMatrixA_dxb_num
    !double precision, dimension(60,60,30,30) :: dMatrixA_dxb_num2
    !--------------------------------------------------------------------------------------------------

    !allocation of large arrays
    if(.not. allocated(gl%ge%gE_C_dxadxb)) allocate(gl%ge%gE_C_dxadxb(nderivs, 30, 30))
    if(.not. allocated(gl%ge%gE_R_dxadxb)) allocate(gl%ge%gE_R_dxadxb(nderivs, 30, 30))
    if(.not. allocated(gl%ge%ln_gamma_C_dxadxb)) allocate(gl%ge%ln_gamma_C_dxadxb(nderivs,30, 30, 30))
    if(.not. allocated(gl%ge%ln_gamma_R_dxadxb)) allocate(gl%ge%ln_gamma_R_dxadxb(nderivs,30, 30, 30))
    if(.not. allocated(ln_gamma_C_Trho_dxa)) allocate(ln_gamma_C_Trho_dxa(nderivs, 30, 30))
    if(.not. allocated(ln_gamma_R_Trho_dxa)) allocate(ln_gamma_R_Trho_dxa(nderivs, 30, 30))
    if(.not. allocated(gE_C_Trho_dxadxb)) allocate(gE_C_Trho_dxadxb(nderivs, 30, 30))
    if(.not. allocated(gE_R_Trho_dxadxb)) allocate(gE_R_Trho_dxadxb(nderivs, 30, 30))
    if(.not. allocated(ln_gamma_C_Trho_dxadxb)) allocate(ln_gamma_C_Trho_dxadxb(nderivs,30, 30, 30))
    if(.not. allocated(ln_gamma_R_Trho_dxadxb)) allocate(ln_gamma_R_Trho_dxadxb(nderivs,30, 30, 30))
    !if(.not. allocated(d2phi_i_dxadxb)) allocate(d2phi_i_dxadxb(30,30,30))
    !if(.not. allocated(d2Theta_i_dxadxb)) allocate(d2Theta_i_dxadxb(30,30,30))
    !
    !if(.not. allocated(d2lnGamma_k_dxadxb)) allocate(d2lnGamma_k_dxadxb(100,30,30))
    !if(.not. allocated(d2Theta_m_dxadxb)) allocate(d2Theta_m_dxadxb(100,30,30))
    !if(.not. allocated(d2X_m_dxadxb)) allocate(d2X_m_dxadxb(100,30,30))

    if(.not. allocated(gl%ge%ln_gamma_C_dxa)) allocate(gl%ge%ln_gamma_C_dxa(nderivs, 30, 30))
    if(.not. allocated(gl%ge%ln_gamma_R_dxa)) allocate(gl%ge%ln_gamma_R_dxa(nderivs, 30, 30))

    if(.not. allocated(MatrixA_mix)) allocate(MatrixA_mix(60, 60, 51, 30))
    if(.not. allocated(dFi_dseggamma)) allocate(dFi_dseggamma(51, 51))
    if(.not. allocated(d2F_dxadxb)) allocate(d2F_dxadxb(51, 30, 30))
    if(.not. allocated(d2Fi_dseggammadxb_mix)) allocate(d2Fi_dseggammadxb_mix(51, 51))
    if(.not. allocated(adj_A)) allocate(adj_A(60,60))
    if(.not. allocated(adj_B)) allocate(adj_B(60,60))
    if(.not. allocated(dMatrixA_dxb)) allocate(dMatrixA_dxb(60,60))
    if(.not. allocated(dMatrixB_dxb)) allocate(dMatrixB_dxb(60,60))
    if(.not. allocated(MatrixA)) allocate(MatrixA(60,60))
    if(.not. allocated(MatrixB)) allocate(MatrixB(60,60))
    if(.not. allocated(MatrixC)) allocate(MatrixC(60,60))
    if(.not. allocated(MatrixD)) allocate(MatrixD(60,60))
    if(.not. allocated(d2seggamma_dxadxb_mix)) allocate(d2seggamma_dxadxb_mix(51,30,30))
    if(.not. allocated(dlnseggamma_dxadxb_mix)) allocate(dlnseggamma_dxadxb_mix(51,30,30))
    if(.not. allocated(MatrixB_mix)) allocate(MatrixB_mix(51, 51))

    !end of allocation


    errval = 0

    !Check for wrong input to C_or_R
    if ((C_or_R < 0) .or. (C_or_R > 2)) then
        errval = -1111
        gl%ge%gE_C_dxadxb = errval
        gl%ge%gE_R_dxadxb = errval
        gl%ge%ln_gamma_C_dxadxb = errval
        gl%ge%ln_gamma_R_dxadxb = errval
        return
    end if

    recalc = .true.

    !Check whether calculations have to be redone or not.
    !If the temperature, (density), composition, equation of state and, the required part of UNIFAC stayed the same, no need to recalculate
    if(dabs(gl%ge%Temp_dxadxb_prev - Temp) > 1.D-16) then
        recalc = .true.
    end if
    do i=1,gl.ncomp
        if(dabs(gl%ge%molfrac_dxadxb_prev(i) - gl.molfractions(i)) > 1.D-16) then
            recalc = .true.
        end if
        if(gl%ge%Eq_type_dxadxb_prev(i) .ne. gl.Eq_type(i)) then
            recalc = .true.
        end if
    end do
    if (gl%ge%mixtype_dxadxb_prev .ne. gl.mix_type) then
        recalc = .true.
    end if
    if (gl%ge%C_or_R_dxadxb_prev .ne. C_or_R) then
        recalc = .true.
    end if


    !Initialize variables
    ln_gamma_C_Trho = 0.D0
    ln_gamma_R_Trho = 0.D0
    gE_C_Trho = 0.D0
    gE_R_Trho = 0.d0
    ln_gamma_C_Trho_dxa = 0.D0
    ln_gamma_R_Trho_dxa = 0.d0
    gE_C_Trho_dxa = 0.D0
    gE_R_Trho_dxa = 0.D0
    gl%ge%ln_gamma_C_dxa = 0.D0
    gl%ge%ln_gamma_R_dxa = 0.D0
    gl%ge%gE_C_dxa = 0.D0
    gl%ge%gE_R_dxa = 0.D0
    gl%ge%gE_C_dxadxb = 0.D0
    gl%ge%gE_R_dxadxb = 0.D0
    gl%ge%ln_gamma_C_dxadxb = 0.D0
    gl%ge%ln_gamma_R_dxadxb = 0.D0
    gE_C_Trho_dxadxb = 0.D0
    gE_R_Trho_dxadxb = 0.D0
    ln_gamma_C_Trho_dxadxb = 0.D0
    ln_gamma_R_Trho_dxadxb = 0.D0

    !Value for the ideal gas constant
    R_const = 8.3144598D0
    R_const_cal = 8.3144598D0 / 4184.D0

    !Get mole fractions of the mixture from module variable
    x = gl.molfractions

    !Calculate tau
    tau = gl.tredmix / Temp

    !Set dummy value for delta (Not needed because UNIFAC is not a function of density. However, already implemented here for possible later modifications)
    del = 1.D0

    !Derivatives of tredmix with respect to xa are required for the following derivatives
    call dYr_dxi(gl,dTred_dxa, drhored_dxa)

    !Derivatives of tredmix with respect to xa and xb are required for the following derivatives
    call d2Yr_dxi2(gl,d2Tred_dxa2, drhored_dxa2)
    call d2Yr_dxidxj(gl,d2Tred_dxadxb, drhored_dxadxb)
    !write the derivatives with respect to two times the same mole fraction also in the variable d2Tred_dxadxb
    do i = 1, gl.ncomp-1
        d2Tred_dxadxb(i,i) = d2Tred_dxa2(i)
    end do

    !Get the required derivatives of activities and excess Gibbs energies with respect tau (and not with respect to xa and xb, however, needed for the derivatives here)
    GETDER_no_xaxb = 0
    if (GETDER(1) .eq. 1)  then
        GETDER_no_xaxb(1) = 1
        GETDER_no_xaxb(4) = 1
        GETDER_no_xaxb(5) = 1
    end if
    call gE_COSMO_SAC_MIXDERIVS(gl,Temp, GETDER_no_xaxb, C_or_R, errval)
    if (GETDER(1) .eq. 1) then
        !Transform the tau derivatives into the needed temperature derivatives
        ln_gamma_C_Trho(1,:) = gl%ge%ln_gamma_C(1,:)
        gE_C_Trho(1) = gl%ge%gE_C(1)
        ln_gamma_R_Trho(1,:) = gl%ge%ln_gamma_R(1,:)
        gE_R_Trho(1) = gl%ge%gE_R(1)
        ln_gamma_C_Trho(4,:) = - tau**2 / gl.tredmix * gl%ge%ln_gamma_C(4,:)
        gE_C_Trho(4) = - tau**2 / gl.tredmix * gl%ge%gE_C(4)
        ln_gamma_R_Trho(4,:) = - tau**2 / gl.tredmix * gl%ge%ln_gamma_R(4,:)
        gE_R_Trho(4) = - tau**2 / gl.tredmix * gl%ge%gE_R(4)
        ln_gamma_C_Trho(5,:) = tau**4 / gl.tredmix**2 * gl%ge%ln_gamma_C(5,:) - 2.D0 * tau / gl.tredmix * ln_gamma_C_Trho(4,:)
        gE_C_Trho(5) = tau**4 / gl.tredmix**2 * gl%ge%gE_C(5) - 2.D0 * tau / gl.tredmix * gE_C_Trho(4)
        ln_gamma_R_Trho(5,:) = tau**4 / gl.tredmix**2 * gl%ge%ln_gamma_R(5,:) - 2.D0 * tau / gl.tredmix * ln_gamma_R_Trho(4,:)
        gE_R_Trho(5) = tau**4 / gl.tredmix**2 * gl%ge%gE_R(5) - 2.D0 * tau / gl.tredmix * gE_R_Trho(4)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !test, rewmove after test
        !gE_R_Trho(5) = gl.cosmo.gE_R_Trho_num(5)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !COSMO-SAC variables needed here
        sigma_profile_mix = gl.cosmo.sigma_profile_mix_gl
        dFi_dseggamma = gl.cosmo.dFi_mix_v1_gl
        seggamma = gl.cosmo.seggamma_gl

    end if



    !Get the required derivatives of activities and excess Gibbs energies with respect tau (and not with respect to xa and xb, however, needed for the derivatives here)
    GETDER_no_xb = 0
    if (GETDER(1) .eq. 1)  then
        GETDER_no_xb(1) = 0
        GETDER_no_xb(4) = 1
        !GETDER_no_xb(5) = 1
    end if
    call gE_COSMO_SAC_MIXDERIVS_dxa(gl,Temp, GETDER_no_xb, C_or_R, errval)
    if (GETDER(1) .eq. 1) then
        !Transform the tau derivatives into the needed temperature derivatives
        do a = 1, gl.ncomp-1
            do i = 1, gl.ncomp
                !ln_gamma_C_Trho_dxa(1,i,a) =  gl.ln_gamma_C_dxa(1,i,a) - dTred_dxa(a) / tau * ln_gamma_C_Trho(4,i)
                ln_gamma_R_Trho_dxa(1,i,a) =  gl%ge%ln_gamma_R_dxa(1,i,a) - dTred_dxa(a) / tau * ln_gamma_R_Trho(4,i)
                !ln_gamma_C_Trho_dxa(4,i,a) = - %ge%tau**2 / gl.tredmix * gl.ln_gamma_C_dxa(4,i,a) &
                !& - dTred_dxa(a) / gl.tredmix * ln_gamma_C_Trho(4,i) &
                !& - dTred_dxa(a) / tau * ln_gamma_C_Trho(5,i)
                ln_gamma_R_Trho_dxa(4,i,a) = - tau**2 / gl.tredmix * gl%ge%ln_gamma_R_dxa(4,i,a) &
                    & - dTred_dxa(a) / gl.tredmix * ln_gamma_R_Trho(4,i) &
                    & - dTred_dxa(a) / tau * ln_gamma_R_Trho(5,i)
            end do
            !gE_C_Trho_dxa(1,a) = gl.gE_C_dxa(1,a) - dTred_dxa(a) / tau * gE_C_Trho(4)
            gE_R_Trho_dxa(1,a) = gl%ge%gE_R_dxa(1,a) - dTred_dxa(a) / tau * gE_R_Trho(4)
            !gE_C_Trho_dxa(4,a) = - tau**2 / gl.tredmix * gl.gE_C_dxa(4,a) &
            !& - dTred_dxa(a) / gl.tredmix * gE_C_Trho(4) &
            !& - dTred_dxa(a) / tau * gE_C_Trho(5)
            gE_R_Trho_dxa(4,a) = - tau**2 / gl.tredmix * gl%ge%gE_R_dxa(4,a) &
                & - dTred_dxa(a) / gl.tredmix * gE_R_Trho(4) &
                & - dTred_dxa(a) / tau * gE_R_Trho(5)
        end do

        !COSMO-SAC variables
        dFi_dxa = gl.cosmo.dFi_dxa
        dseggamma_dxa = gl.cosmo.dseggamma_dxa_mix
    end if



    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 1)) then
        !Calculate the combinatorial part of the activity coefficients
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !ln_gamma_C is always needed for the following calculations, thus it is always calculated

        !if ((gl%ge%GETDER_dxadxb_prev(1) == 1) .and. (recalc .eqv. .false.)) then
        !
        !    gl.ln_gamma_C_dxadxb(1,:,:,:) = gl%ge%ln_gamma_C_dxadxb_prev(1,:,:,:)
        !
        !    !Transformation of variables from tau,delta, and x to T,rho,x
        !    do b = 1, gl.ncomp-1
        !        do a = 1, gl.ncomp-1
        !            do i = 1, gl.ncomp
        !                ln_gamma_C_Trho_dxadxb(1,i,a,b) = gl.ln_gamma_C_dxadxb(1,i,a,b) &
        !                                             & - dTred_dxa(a) * dTred_dxa(b) / tau**2 * ln_gamma_C_Trho(5,i) &
        !                                             & - dTred_dxa(a) / tau * ln_gamma_C_Trho_dxa(4,i,b) &
        !                                             & - dTred_dxa(b) / tau * ln_gamma_C_Trho_dxa(4,i,a) &
        !                                             & - d2Tred_dxadxb(a,b) / tau * ln_gamma_C_Trho(4,i)
        !            end do
        !        end do
        !    end do
        !
        !else
        !
        !    !Calculate help variables
        !    r_i = 0.D0
        !    q_i = 0.D0
        !    Do i = 1, gl.ncomp
        !        Do k = 1, gl.nr_of_groups_i(i)
        !            r_i(i) = r_i(i) + gl.v_ik(i, k) * gl.R_ik(i, k)
        !            q_i(i) = q_i(i) + gl.v_ik(i, k) * gl.Q_ik(i, k)
        !        End do
        !    End do
        !
        !    sum_rjxj = 0.D0
        !    sum_qjxj = 0.D0
        !    Do j = 1, gl.ncomp
        !        sum_rjxj = sum_rjxj + r_i(j) * x(j)
        !        sum_qjxj = sum_qjxj + q_i(j) * x(j)
        !    End do
        !    Do i = 1, gl.ncomp
        !        phi_i(i) = r_i(i) / sum_rjxj
        !        Theta_i(i) = q_i(i) / sum_qjxj
        !    End do
        !
        !    !Calculate the derivativs of phi_i and Theta_i with respect to xa. These derivatives are also needed for all tau and delta derivatives
        !    Do a = 1, gl.ncomp-1
        !        Do i = 1, gl.ncomp
        !            dphi_i_dxa(i,a) = -r_i(i) * (r_i(a) - r_i(gl.ncomp)) / sum_rjxj**2
        !            dTheta_i_dxa(i,a) = -q_i(i) * (q_i(a) - q_i(gl.ncomp)) / sum_qjxj**2
        !        End do
        !    End do
        !
        !    !Calculate the derivativs of phi_i and Theta_i with respect to xa and xb.
        !    Do b = 1, gl.ncomp -1
        !        Do a = 1, gl.ncomp-1
        !            Do i = 1, gl.ncomp
        !                d2phi_i_dxadxb(i,a,b) = 2.D0 * r_i(i) * (r_i(a) - r_i(gl.ncomp)) * (r_i(b) - r_i(gl.ncomp)) / sum_rjxj**3
        !                d2Theta_i_dxadxb(i,a,b) = 2.D0 * q_i(i) * (q_i(a) - q_i(gl.ncomp)) * (q_i(b) - q_i(gl.ncomp)) / sum_qjxj**3
        !            End do
        !        End do
        !    end do
        !
        !    !Compute the mole fraction derivatives of all activity coefficients of the combinatorial part
        !    do b = 1, gl.ncomp-1
        !        do a = 1, gl.ncomp-1
        !            Do i = 1, gl.ncomp
        !                ln_gamma_C_Trho_dxadxb(1,i,a,b) = (1.D0 / phi_i(i) - 1.D0) * d2phi_i_dxadxb(i,a,b) &
        !                                              & - dphi_i_dxa(i,a) * dphi_i_dxa(i,b) / phi_i(i)**2 &
        !                                              & - 5.D0 * q_i(i) &
        !                                              & * ((Theta_i(i) / phi_i(i) - 1.D0) * (d2phi_i_dxadxb(i,a,b) / Theta_i(i)  &
        !                                              & - phi_i(i) * d2Theta_i_dxadxb(i,a,b) / Theta_i(i)**2 &
        !                                              & - dTheta_i_dxa(i,a) * dphi_i_dxa(i,b) / Theta_i(i)**2 &
        !                                              & - dTheta_i_dxa(i,b) * dphi_i_dxa(i,a) / Theta_i(i)**2 &
        !                                              & + 2.D0 * phi_i(i) * dTheta_i_dxa(i,a) * dTheta_i_dxa(i,b) / Theta_i(i)**3) &
        !                                              & + (Theta_i(i) * dphi_i_dxa(i,a) - phi_i(i) * dTheta_i_dxa(i,a)) &
        !                                              & * (phi_i(i) * dTheta_i_dxa(i,b) - Theta_i(i) * dphi_i_dxa(i,b)) &
        !                                              & / Theta_i(i)**2 / phi_i(i)**2)
        !            End do
        !        End do
        !    end do
        !
        !    !Transformation of variables from T,rho,x to tau,delta,x
        !    do b = 1, gl.ncomp-1
        !        do a = 1, gl.ncomp-1
        !            do i = 1, gl.ncomp
        !                gl.ln_gamma_C_dxadxb(1,i,a,b) = ln_gamma_C_Trho_dxadxb(1,i,a,b) &
        !                                             & + dTred_dxa(a) * dTred_dxa(b) / tau**2 * ln_gamma_C_Trho(5,i) &
        !                                             & + dTred_dxa(a) / tau * ln_gamma_C_Trho_dxa(4,i,b) &
        !                                             & + dTred_dxa(b) / tau * ln_gamma_C_Trho_dxa(4,i,a) &
        !                                             & + d2Tred_dxadxb(a,b) / tau * ln_gamma_C_Trho(4,i)
        !            end do
        !        end do
        !    end do
        !
        !end if
        !
        !if (GETDER(1) .eq. 1) then
        !
        !    if ((gl%ge%GETDER_dxadxb_prev(1) == 1) .and. (recalc .eqv. .false.)) then
        !
        !       gl.gE_C_dxadxb(1,:,:) = gl%ge%gE_C_dxadxb_prev(1,:,:)
        !
        !        !Transformation of variables from tau,delta, and x to T,rho,x
        !        do b = 1, gl.ncomp-1
        !            do a = 1, gl.ncomp-1
        !                    gE_C_Trho_dxadxb(1,a,b) = gl.gE_C_dxadxb(1,a,b) &
        !                                            & - dTred_dxa(a) * dTred_dxa(b) / tau**2 * gE_C_Trho(5) &
        !                                            & - dTred_dxa(a) / tau * gE_C_Trho_dxa(4,b) &
        !                                            & - dTred_dxa(b) / tau * gE_C_Trho_dxa(4,a) &
        !                                            & - d2Tred_dxadxb(a,b) / tau * gE_C_Trho(4)
        !            end do
        !        end do
        !
        !    else
        !
        !        Do b = 1, gl.ncomp-1
        !            Do a = 1,gl.ncomp-1
        !                Do i = 1, gl.ncomp
        !                    gE_C_Trho_dxadxb(1,a,b) = gE_C_Trho_dxadxb(1,a,b) + R_const * Temp * x(i) * ln_gamma_C_Trho_dxadxb(1,i,a,b)
        !                end do
        !                gE_C_Trho_dxadxb(1,a,b) = gE_C_Trho_dxadxb(1,a,b) &
        !                                        & + R_const * Temp * (ln_gamma_C_Trho_dxa(1,a,b) - ln_gamma_C_Trho_dxa(1,gl.ncomp,b)) &
        !                                        & + R_const * Temp * (ln_gamma_C_Trho_dxa(1,b,a) - ln_gamma_C_Trho_dxa(1,gl.ncomp,a))
        !            End do
        !        End do
        !
        !        !Transformation of variables from T,rho,x to tau,delta,x
        !        do b = 1, gl.ncomp-1
        !            do a = 1, gl.ncomp-1
        !                    gl.gE_C_dxadxb(1,a,b) = gE_C_Trho_dxadxb(1,a,b) &
        !                                            & + dTred_dxa(a) * dTred_dxa(b) / tau**2 * gE_C_Trho(5) &
        !                                            & + dTred_dxa(a) / tau * gE_C_Trho_dxa(4,b) &
        !                                            & + dTred_dxa(b) / tau * gE_C_Trho_dxa(4,a) &
        !                                            & + d2Tred_dxadxb(a,b) / tau * gE_C_Trho(4)
        !            end do
        !        end do
        !
        !    end if
        !end if

        !All other derivatives of the combinatorial part are 0 (except for the tau-only derivatives of gE), because the combinatorial part is a function of composition only
        !Derivative wrt delta
        if (GETDER(2) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(2,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(2,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(2,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(2,:,:) = 0.D0
        end if
        !Second derivative wrt delta
        if (GETDER(3) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(3,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(3,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(3,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(3,:,:) = 0.D0
        end if

        !Derivative wrt tau
        if (GETDER(4) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(4,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(4,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(4,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(4,:,:) = 0.D0
        end if

        !Second derivative wrt tau
        if (GETDER(5) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(5,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(5,:,:) = 0.D0
            gl%ge%ln_gamma_C_dxadxb(5,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(5,:,:) = 0.D0
        end if
        !Second derivative wrt delta and tau
        if (GETDER(6) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(6,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(6,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(6,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(6,:,:) = 0.D0
        end if
        !Third derivative wrt delta, tau, and tau
        if (GETDER(7) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(7,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(7,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(7,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(7,:,:) = 0.D0
        end if
        !Third derivative wrt delta
        if (GETDER(8) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(8,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(8,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(8,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(8,:,:) = 0.D0
        end if
        !Third derivative wrt tau
        if (GETDER(9) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(9,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(9,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(9,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(9,:,:) = 0.D0
        end if
        !Third derivative wrt tau, delta, and delta
        if (GETDER(10) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(10,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(10,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(10,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(10,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta
        if (GETDER(11) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(11,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(11,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(11,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(11,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, delta and tau
        if (GETDER(12) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(12,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(12,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(12,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(12,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, tau, tau
        if (GETDER(13) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(13,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(13,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(13,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(13,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, tau, tau, tau
        if (GETDER(14) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(14,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(14,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(14,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(14,:,:) = 0.D0
        end if
        !Fourth derivative wrt tau
        if (GETDER(15) .eq. 1) then
            ln_gamma_C_Trho_dxadxb(15,:,:,:) = 0.D0
            gE_C_Trho_dxadxb(15,:,:) = 0.d0
            gl%ge%ln_gamma_C_dxadxb(15,:,:,:) = 0.d0
            gl%ge%gE_C_dxadxb(15,:,:) = 0.D0
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end if


    if ((C_or_R .eq. 0) .or. (C_or_R .eq. 2)) then
        !Calculate the residual part of the activity coefficients
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !ln_gamma_R and all other quantities are AT THE MOMENT NOT always needed for the following calculations
        !As soon as higher derivatives are required, this may change!
        if ((gl%ge%GETDER_dxadxb_prev(1) == 1) .and. (recalc .eqv. .false.)) then

            gl%ge%ln_gamma_R_dxadxb(1,:,:,:) = gl%ge%ln_gamma_R_dxadxb_prev(1,:,:,:)

            !Transformation of variables from tau,delta, and x to T,rho,x
            do b = 1, gl.ncomp-1
                do a = 1, gl.ncomp-1
                    do i = 1, gl.ncomp
                        ln_gamma_R_Trho_dxadxb(1,i,a,b) = gl%ge%ln_gamma_R_dxadxb(1,i,a,b) &
                            & - dTred_dxa(a) * dTred_dxa(b) / tau**2 * ln_gamma_R_Trho(5,i) &
                            & - dTred_dxa(a) / tau * ln_gamma_R_Trho_dxa(4,i,b) &
                            & - dTred_dxa(b) / tau * ln_gamma_R_Trho_dxa(4,i,a) &
                            & - d2Tred_dxadxb(a,b) / tau * ln_gamma_R_Trho(4,i)
                    end do
                end do
            end do

        else
            !------------------------------------------------------------------------------------------------------
            !test
            !------------------------------------------------------------------------------------------------------

            !!!Numerical derivative of the excess based departure function with respect to xi at constant tau and del
            !delta_x = 1.D-4
            !molfractions_orig = gl.molfractions
            !tredmix_orig = gl.tredmix
            !rhoredmix_orig = gl.rhoredmix
            !getder_dep = 1
            !Temperature = Temp
            !
            !    do b = 1, gl.ncomp-1
            !        !Increase xi
            !        gl.molfractions(b) = molfractions_orig(b) + delta_x
            !        gl.molfractions(gl.ncomp) = molfractions_orig(gl.ncomp) - delta_x
            !        call reduced_parameters_calc(gl, Temperature)
            !        Temp_p = gl.tredmix/tredmix_orig * Temperature
            !        !Temp_p = Temp
            !        call gE_COSMO_SAC_MIXDERIVS_dxa(gl, Temp_p, getder_dep, C_or_R, errval)
            !        dFi_dxa_p = gl.cosmo.dFi_dxa
            !        Fi_seggamma_mix_p(:,:,b) = gl.cosmo.dFi_mix_v1_gl
            !        dseggamma_dxa_p = gl.cosmo.dseggamma_dxa_mix
            !        ln_gamma_R_dxa_p(1,:,:) = gl.ln_gamma_R_dxa(1,:,:)
            !        gE_R_Trho_dxa_p(1,:) = gl.cosmo.gE_R_Trho_dxa_num(1,:)
            !        gE_R_dxa_p(1,:) = gl.gE_R_dxa(1,:)
            !        detMatrixA_mix_num_p = gl.cosmo.detMatrixA_mix_num
            !        MatrixA_p = gl.cosmo.MatrixA
            !
            !        !Decrease xi
            !        gl.molfractions(b) = molfractions_orig(b) - delta_x
            !        gl.molfractions(gl.ncomp) = molfractions_orig(gl.ncomp) + delta_x
            !        call reduced_parameters_calc(gl, Temperature)
            !        Temp_m = gl.tredmix/tredmix_orig * Temperature
            !        !Temp_m = Temp
            !        call gE_COSMO_SAC_MIXDERIVS_dxa(gl, Temp_m, getder_dep, C_or_R, errval)
            !        dFi_dxa_m = gl.cosmo.dFi_dxa
            !        Fi_seggamma_mix_m(:,:,b) = gl.cosmo.dFi_mix_v1_gl
            !        dseggamma_dxa_m = gl.cosmo.dseggamma_dxa_mix
            !        ln_gamma_R_dxa_m(1,:,:) = gl.ln_gamma_R_dxa(1,:,:)
            !        gE_R_Trho_dxa_m(1,:) = gl.cosmo.gE_R_Trho_dxa_num(1,:)
            !        gE_R_dxa_m(1,:) = gl.gE_R_dxa(1,:)
            !        detMatrixA_mix_num_m = gl.cosmo.detMatrixA_mix_num
            !        MatrixA_m = gl.cosmo.MatrixA
            !
            !        !Calculate the numerical derivative
            !        d2Fi_dxadxb_num(:,:,b) = (dFi_dxa_p - dFi_dxa_m) / (2.D0 * delta_x)
            !        d2F_dseggammadxb_num(:,:,b) = (Fi_seggamma_mix_p(:,:,b) - Fi_seggamma_mix_m(:,:,b)) / (2.D0 * delta_x)
            !        ddetMatrixA_dxb_mix_num(:,:,b) = (detMatrixA_mix_num_p - detMatrixA_mix_num_m) / (2.D0 * delta_x)
            !        dMatrixA_dxb_num(:,:,:,:,b) = (MatrixA_p - MatrixA_m ) / (2.D0 * delta_x)
            !        d2segamma_dxadxb_num(:,:,b) = (dseggamma_dxa_p - dseggamma_dxa_m) / (2.D0 * delta_x)
            !        ln_gamma_R_dxadxb_num(1,:,:,b) = (ln_gamma_R_dxa_p(1,:,:) - ln_gamma_R_dxa_m(1,:,:)) / (2.D0 * delta_x)
            !        gE_R_Trho_dxadxb_num(1,:,b) = (gE_R_Trho_dxa_p(1,:) - gE_R_Trho_dxa_m(1,:)) / (2.D0 * delta_x)
            !        gE_R_dxadxb_num(1,:,b) = (gE_R_dxa_p(1,:) - gE_R_dxa_m(1,:)) / (2.D0 * delta_x)
            !        gl.molfractions = molfractions_orig
            !    end do
            !        dMatrixA_dxb_num2(:,:,:,:) = dMatrixA_dxb_num(:,:,22,:,:)
            !
            !call reduced_parameters_calc(gl, Temperature)
            !call gE_COSMO_SAC_MIXDERIVS_dxa(gl,Temp, getder_dep, C_or_R, errval)
            !------------------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            !mixture
            !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            ! Calculation Determinante of Matrix B (only needs to be calculated ones)
            !Matrix partial derivatives in respect to seggamma only
            m = 0
            do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                n = 0
                if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                    cycle
                end if
                m = m + 1
                do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                    if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                        cycle
                    end if
                    n = n + 1
                    MatrixB(m,n) = dFi_dseggamma(j,k)
                    MatrixB_mix(j,k) = MatrixB(m,n)
                end do
            end do
            rankA = m
            call Gauss_algorithm(gl,MatrixB, rankA, vectorb, vectorx, Det_seggamma, errval)
            DetMatrixB_mix = Det_seggamma

            ! Calculation Determinante of Matrix A
            do a = 1, gl.ncomp-1
                do o = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !all derivatives of seggamma in respect to xa
                    if (abs(gl.cosmo.sigma_profile_mix_gl(o)) < 1.D-14) then
                        cycle
                    end if

                    m = 0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Zeile MatrixA
                        if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        n = 0
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Spalte (Ableitung) MatrixA
                            if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            if (k == o) then
                                MatrixA(m,n) = gl.cosmo.dFi_dxa(j,a)
                                MatrixA_mix(j,k,o,a) = MatrixA(m,n)
                            else
                                !partial derivatives in respect to seggamma
                                MatrixA(m,n) = dFi_dseggamma(j,k)
                                MatrixA_mix(j,k,o,a) = MatrixA(m,n)
                            end if
                        end do
                    end do

                    rankA = n
                    vectorb = 0.D0

                    call Gauss_algorithm(gl,MatrixA, rankA, vectorb, vectorx, Det_xa, errval)
                    DetMatrixA_mix(o,a) = Det_xa
                end do
            end do

            !sums for second derivative of help function Fi in respect to xa and xb
            sum_xiAi = 0.D0
            sum_xiAisigma = 0.D0
            do i = 1, gl.ncomp
                sum_xiAi = sum_xiAi + x(i) * gl.cosmo.Acosmo(i)
            end do
            do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                do i = 1, gl.ncomp
                    sum_xiAisigma(j) = sum_xiAisigma(j) + x(i) * gl.cosmo.sigma(j,i)
                end do
            end do

            do b = 1, gl.ncomp-1
                do a = 1, gl.ncomp-1
                    !calculate second derivative of help function Fi in respect to xa and xb
                    sum_d2F_dxadxb = 0.D0
                    sum_d2F_dxadxb2 = 0.D0
                    !sum_d2F_dxadxb3 = 0.D0
                    !sum_d2F_dxadxb4 = 0.D0
                    !sum_d2F_dxadxb5 = 0.D0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                        !    cycle
                        !end if
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            !if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                            !    cycle
                            !end if
                            !Sums for d2F_dxadxb
                            sum_d2F_dxadxb(j) = sum_d2F_dxadxb(j) + seggamma(k) * ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k)) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                            !sum_d2F_dxadxb2(j) = sum_d2F_dxadxb2(j) + exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (gl.cosmo.dseggamma_dxa_mix(k,b) * ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k)) + &
                            sum_d2F_dxadxb2(j) = sum_d2F_dxadxb2(j) + exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggamma_dxa(k,b) * ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k)) + &
                                & seggamma(k) * ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * (gl.cosmo.Acosmo(b) - gl.cosmo.Acosmo(gl.ncomp)) - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * (gl.cosmo.sigma(k,b) - gl.cosmo.sigma(k,gl.ncomp))))

                            !sum_d2F_dxadxb3(j) = sum_d2F_dxadxb3(j) + seggamma(k) * ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k)) / ( sum_xiAi ** 2.D0) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                            !sum_d2F_dxadxb4(j) = sum_d2F_dxadxb4(j) + dseggamma_dxa(k,b) * ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k)) / ( sum_xiAi ** 2.D0) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                            !sum_d2F_dxadxb5(j) = sum_d2F_dxadxb5(j) + seggamma(k) * (((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * (gl.cosmo.Acosmo(b) - gl.cosmo.Acosmo(gl.ncomp)) - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * (gl.cosmo.sigma(k,b) - gl.cosmo.sigma(k,gl.ncomp))) * sum_xiAi - &
                            !                                        & 2.D0 * (gl.cosmo.Acosmo(b) - gl.cosmo.Acosmo(gl.ncomp)) * ((gl.cosmo.sigma(k,a) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(a) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k))) / ( sum_xiAi ** 3.D0) * exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp))
                        end do
                        !d2F_dxadxb(j,a,b) = sum_d2F_dxadxb(j) * (gl.cosmo.dseggamma_dxa_mix(j,b) * sum_xiAi - 2.D0 * (gl.cosmo.Acosmo(b) - gl.cosmo.Acosmo(gl.ncomp)) * seggamma(j)) / sum_xiAi ** 3.D0 + seggamma(j) / sum_xiAi ** 2.D0 * sum_d2F_dxadxb2(j)
                        d2F_dxadxb(j,a,b) = sum_d2F_dxadxb(j) * (dseggamma_dxa(j,b) * sum_xiAi - 2.D0 * (gl.cosmo.Acosmo(b) - gl.cosmo.Acosmo(gl.ncomp)) * seggamma(j)) / sum_xiAi ** 3.D0 + seggamma(j) / sum_xiAi ** 2.D0 * sum_d2F_dxadxb2(j)
                        !d2F_dxadxb2(j,a,b) = dseggamma_dxa(j,b) * sum_d2F_dxadxb3(j) + seggamma(j) * ( sum_d2F_dxadxb4(j) + sum_d2F_dxadxb5(j))
                    end do
                    !calculate second derivative of help function Fi in respect to seggamma and xb
                    !sum for d2Fi_dseggammadxb_mix
                    !könnte außerhalb der a Schleife sein
                    sum_d2Fi_dseggammadxb = 0.D0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                        !    cycle
                        !end if
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            !if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                            !    cycle
                            !end if
                            !sum_d2Fi_dseggammadxb(j) = sum_d2Fi_dseggammadxb(j) + exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (gl.cosmo.dseggamma_dxa_mix(k,b) * gl.cosmo.sigma_profile_mix_gl(k) + &
                            sum_d2Fi_dseggammadxb(j) = sum_d2Fi_dseggammadxb(j) + exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggamma_dxa(k,b) * gl.cosmo.sigma_profile_mix_gl(k) + &
                                & seggamma(k) / sum_xiAi ** 2.D0 * ((gl.cosmo.sigma(k,b) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(b) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k)))
                        end do
                    end do

                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                        !    cycle
                        !end if
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            !if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                            !    cycle
                            !end if
                            if (j == k) then
                                !d2Fi_dseggammadxb_mix(j,k) = sum_d2Fi_dseggammadxb(j) + exp(-gl.cosmo.delta_w_gl(j,j) / (R_const_cal*Temp)) * (gl.cosmo.dseggamma_dxa_mix(j,b) * gl.cosmo.sigma_profile_mix_gl(j) + &
                                d2Fi_dseggammadxb_mix(j,k) = sum_d2Fi_dseggammadxb(j) + exp(-gl.cosmo.delta_w_gl(j,j) / (R_const_cal*Temp)) * (dseggamma_dxa(j,b) * gl.cosmo.sigma_profile_mix_gl(j) + &
                                    & seggamma(j) / sum_xiAi ** 2.D0 * ((gl.cosmo.sigma(j,b) - gl.cosmo.sigma(j,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(b) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(j)))
                            else
                                !d2Fi_dseggammadxb_mix(j,k) = exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (gl.cosmo.dseggamma_dxa_mix(j,b) * gl.cosmo.sigma_profile_mix_gl(k) + &
                                d2Fi_dseggammadxb_mix(j,k) = exp(-gl.cosmo.delta_w_gl(j,k) / (R_const_cal*Temp)) * (dseggamma_dxa(j,b) * gl.cosmo.sigma_profile_mix_gl(k) + &
                                    & seggamma(j) / sum_xiAi ** 2.D0 * ((gl.cosmo.sigma(k,b) - gl.cosmo.sigma(k,gl.ncomp)) * sum_xiAi - (gl.cosmo.Acosmo(b) - gl.cosmo.Acosmo(gl.ncomp)) * sum_xiAisigma(k)))
                            end if
                        end do
                    end do
                    !numerical test,REMOVE AFTER TESTING!!!!!!
                    !d2F_dxadxb = gl.cosmo.d2F_dxadxb_num(:,b)
                    !d2Fi_dseggammadxb_mix = gl.cosmo.d2F_dseggammadxa_num(:,:,b)

                    !Matrix dB/dxb, second partial derivatives of B in respect to xb
                    m = 0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        n = 0
                        if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            dMatrixB_dxb(m,n) = d2Fi_dseggammadxb_mix(j,k)
                        end do
                    end do

                    !Adjugates of MatrixA and MatrixB
                    MatrixB = 0.D0
                    m = 0
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        n = 0
                        if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                            cycle
                        end if
                        m = m + 1
                        do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                                cycle
                            end if
                            n = n + 1
                            MatrixB(m,n) = MatrixB_mix(j,k)
                        end do
                    end do

                    rankA = m

                    call Adjugate(gl,MatrixB, rankA, adj_B, errval)

                    !Multiply adj(MatrixA) * dMatrixAdxb
                    call Mat_mult(gl,adj_B, rankA, rankA, dMatrixB_dxb, rankA, rankA, MatrixD, errval)

                    !calculate dDetA_dxb = Trace(Adjugate(A) dMatrixA/dxb)
                    call Trace(gl,MatrixD, rankA, dDetB_dxb_mix, errval)

                    !matrix dA/dxb, partial derivatives of A in respect to xb
                    do o = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)   !all derivatuves of seggamma
                        if (abs(gl.cosmo.sigma_profile_mix_gl(o)) < 1.D-14) then
                            cycle
                        end if
                        m = 0
                        dMatrixA_dxb = dMatrixB_dxb
                        do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Zeile
                            if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                                cycle
                            end if
                            m = m + 1
                            n = 0
                            do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)      !Spalte (Ableitung)
                                if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                                    cycle
                                end if
                                n = n + 1
                                if (k == o) then
                                    dMatrixA_dxb(m,n) = d2F_dxadxb(j,a,b)
                                else
                                end if
                            end do
                        end do

                        !Adjugates of MatrixA and MatrixB
                        MatrixA = 0.D0
                        !MatrixB = 0.D0
                        m = 0
                        do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            n = 0
                            if (abs(gl.cosmo.sigma_profile_mix_gl(j)) < 1.D-14) then
                                cycle
                            end if
                            m = m + 1
                            do k = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                                if (abs(gl.cosmo.sigma_profile_mix_gl(k)) < 1.D-14) then
                                    cycle
                                end if
                                n = n + 1
                                MatrixA(m,n) = MatrixA_mix(j,k,o,a)
                                !MatrixB(m,n) = MatrixB_mix(j,k)
                            end do
                        end do

                        rankA = m

                        call Adjugate(gl,MatrixA, rankA, adj_A, errval)
                        !call Adjugate(gl,MatrixB, rankA, adj_B, errval)

                        !Multiply adj(MatrixA) * dMatrixAdxb
                        call Mat_mult(gl,adj_A, rankA, rankA, dMatrixA_dxb, rankA, rankA, MatrixC, errval)
                        !call Mat_mult(gl,adj_B, rankA, rankA, dMatrixB_dxb, rankA, rankA, MatrixD, errval)

                        !calculate dDetA_dxb = Trace(Adjugate(A) dMatrixA/dxb)
                        call Trace(gl,MatrixC, rankA, dDetA_dxb_mix(o), errval)
                        !call Trace(gl,MatrixD, rankA, dDetB_dxb_mix(o), errval)
                    end do
                    !dDetA_dxb_mix = gl.cosmo.dDetAdxb_mix_numeric(:,a)
                    !dDetB_dxb_mix = gl.cosmo.dDetBdt_mix_numeric
                    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    !Calculate second derivative of segment gamma wrt T for pure substance and mixture
                    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                    !mixture
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        d2seggamma_dxadxb_mix(j,a,b) = - ( dDetA_dxb_mix(j) * detMatrixB_mix - dDetB_dxb_mix * detMatrixA_mix(j,a)) / (detMatrixB_mix ** 2.D0)
                    end do
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !test, REMOVE AFTER TESTING
                    !d2seggamma_dxadxb_mix = d2segamma_dxadxb_num
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    !Derivations are dy/dx, need to be transformed to dlny/dx;
                    do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                        !dlnseggamma_dxadxb_mix(j) = -1.D0 / seggamma(j) ** 2.D0 * gl.cosmo.dseggamma_dxa_mix(j,b) * gl.cosmo.dseggamma_dxa_mix(j,a) + d2seggamma_dxadxb_mix(j) / seggamma(j)
                        dlnseggamma_dxadxb_mix(j,a,b) = -1.D0 / seggamma(j) ** 2.D0 * dseggamma_dxa(j,b) * dseggamma_dxa(j,a) + d2seggamma_dxadxb_mix(j,a,b) / seggamma(j)
                    end do

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !test, REMOVE AFTER TESTING
                    !dlnseggamma_dxadxb_mix(:) = gl.cosmo.d2lnseggammadxadT_mix_numeric(:,b)
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    !----------------------------------------------------------------------


                    !Calculate the derivative of the logarithm of the residual activity coefficients with respect to xa
                    !!!NOTE THAT LNGAMMA_ki OF PURE COMPONENT i IS NOT A FUNCTION OF COMPOSITION, HENCE dlnGamma_ki_dxa(k,i,:) = 0!!! This is indicated with the 0.D0 in the following equation
                    do i = 1, gl.ncomp
                        do j = gl.cosmo.start_sigma(1), gl.cosmo.end_sigma(1)
                            ln_gamma_R_Trho_dxadxb(1,i,a,b) = ln_gamma_R_Trho_dxadxb(1,i,a,b) + gl.cosmo.sigma(j,i) / gl.cosmo.aeff_gl * (dlnseggamma_dxadxb_mix(j,a,b) - 0.D0)
                        end do
                    end do
                end do
            end do

            !Transformation of variables from T,rho,x to tau,delta,x
            do b = 1, gl.ncomp-1
                do a = 1, gl.ncomp-1
                    do i = 1, gl.ncomp
                        gl%ge%ln_gamma_R_dxadxb(1,i,a,b) = ln_gamma_R_Trho_dxadxb(1,i,a,b) &
                            & + dTred_dxa(a) * dTred_dxa(b) / tau**2 * ln_gamma_R_Trho(5,i) &
                            & + dTred_dxa(a) / tau * ln_gamma_R_Trho_dxa(4,i,b) &
                            & + dTred_dxa(b) / tau * ln_gamma_R_Trho_dxa(4,i,a) &
                            & + d2Tred_dxadxb(a,b) / tau * ln_gamma_R_Trho(4,i)
                    end do
                end do
            end do

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !test, REMOVE AFTER TESTING
            !gl.ln_gamma_R_dxadxb = ln_gamma_R_dxadxb_num
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end if

        !Note in the following that all derivatives with respect to delta at constant tau and x are 0, because the residual part is not a function of density
        if (GETDER(1) .eq. 1) then

            if ((gl%ge%GETDER_dxadxb_prev(1) == 1) .and. (recalc .eqv. .false.)) then

                gl%ge%gE_R_dxadxb(1,:,:) = gl%ge%gE_R_dxadxb_prev(1,:,:)

                !Transformation of variables from tau,delta, and x to T,rho,x
                do b = 1, gl.ncomp-1
                    do a = 1, gl.ncomp-1
                        gE_R_Trho_dxadxb(1,a,b) = gl%ge%gE_R_dxadxb(1,a,b) &
                            & - dTred_dxa(a) * dTred_dxa(b) / tau**2 * gE_R_Trho(5) &
                            & - dTred_dxa(a) / tau * gE_R_Trho_dxa(4,b) &
                            & - dTred_dxa(b) / tau * gE_R_Trho_dxa(4,a) &
                            & - d2Tred_dxadxb(a,b) / tau * gE_R_Trho(4)
                    end do
                end do

            else

                Do b = 1, gl.ncomp-1
                    Do a = 1,gl.ncomp-1
                        Do i = 1, gl.ncomp
                            gE_R_Trho_dxadxb(1,a,b) = gE_R_Trho_dxadxb(1,a,b) + R_const * Temp * x(i) * ln_gamma_R_Trho_dxadxb(1,i,a,b)
                        end do
                        gE_R_Trho_dxadxb(1,a,b) = gE_R_Trho_dxadxb(1,a,b) &
                            & + R_const * Temp * (ln_gamma_R_Trho_dxa(1,a,b) - ln_gamma_R_Trho_dxa(1,gl.ncomp,b)) &
                            & + R_const * Temp * (ln_gamma_R_Trho_dxa(1,b,a) - ln_gamma_R_Trho_dxa(1,gl.ncomp,a))
                    End do
                End do

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !test; REMOVE AFTER TESTING
                !gE_R_Trho_dxadxb = gE_R_Trho_dxadxb_num
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                !Transformation of variables from T,rho,x to tau,delta,x
                do b = 1, gl.ncomp-1
                    do a = 1, gl.ncomp-1
                        gl%ge%gE_R_dxadxb(1,a,b) = gE_R_Trho_dxadxb(1,a,b) &
                            & + dTred_dxa(a) * dTred_dxa(b) / tau**2 * gE_R_Trho(5) &
                            & + dTred_dxa(a) / tau * gE_R_Trho_dxa(4,b) &
                            & + dTred_dxa(b) / tau * gE_R_Trho_dxa(4,a) &
                            & + d2Tred_dxadxb(a,b) / tau * gE_R_Trho(4)
                        !gl.gE_R_dxadxb(1,a,b) = gE_R_Trho_dxadxb(1,a,b) + dTred_dxa(a) * dTred_dxa(b) / tau**2 * gE_R_Trho(5) + dTred_dxa(a) / tau * gE_R_Trho_dxa(4,b) + dTred_dxa(b) / tau * gE_R_Trho_dxa(4,a) + d2Tred_dxadxb(a,b) / tau * gE_R_Trho(4)
                    end do
                end do
                !test, REMOVE AFTER TESTING
                !gl.gE_R_dxadxb(1,1,:) = gl.cosmo.gE_R_dxadxb_num(:)
            end if

        end if

        !Second derivative wrt delta and tau
        if (GETDER(2) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(2,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(2,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(2,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(2,:,:) = 0.D0
        end if

        !Third derivative wrt delta, tau, and tau
        if (GETDER(3) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(3,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(3,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(3,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(3,:,:) = 0.D0
        end if

        !Derivative wrt tau
        if (GETDER(4) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(4,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(4,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(4,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(4,:,:) = 0.D0
        end if

        !Second derivative wrt tau
        if (GETDER(5) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(5,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(5,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(5,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(5,:,:) = 0.D0
        end if
        !Second derivative wrt delta and tau
        if (GETDER(6) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(6,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(6,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(6,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(6,:,:) = 0.D0
        end if
        !Third derivative wrt delta, tau, and tau
        if (GETDER(7) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(7,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(7,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(7,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(7,:,:) = 0.D0
        end if
        !Third derivative wrt delta
        if (GETDER(8) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(8,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(8,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(8,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(8,:,:) = 0.D0
        end if
        !Third derivative wrt tau
        if (GETDER(9) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(9,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(9,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(9,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(9,:,:) = 0.D0
        end if
        !Third derivative wrt tau, delta, and delta
        if (GETDER(10) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(10,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(10,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(10,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(10,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta
        if (GETDER(11) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(11,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(11,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(11,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(11,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, delta and tau
        if (GETDER(12) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(12,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(12,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(12,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(12,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, delta, tau, tau
        if (GETDER(13) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(13,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(13,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(13,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(13,:,:) = 0.D0
        end if
        !Fourth derivative wrt delta, tau, tau, tau
        if (GETDER(14) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(14,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(14,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(14,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(14,:,:) = 0.D0
        end if
        !Fourth derivative wrt tau
        if (GETDER(15) .eq. 1) then
            ln_gamma_R_Trho_dxadxb(15,:,:,:) = 0.D0
            gE_R_Trho_dxadxb(15,:,:) = 0.d0
            gl%ge%ln_gamma_R_dxadxb(15,:,:,:) = 0.d0
            gl%ge%gE_R_dxadxb(15,:,:) = 0.D0
        end if
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    End if


    !Write new values in save variables
    gl%ge%Temp_dxadxb_prev = Temp
    gl%ge%GETDER_dxadxb_prev = GETDER
    gl%ge%molfrac_dxadxb_prev = gl.molfractions
    gl%ge%Eq_type_dxadxb_prev = gl.Eq_type
    gl%ge%mixtype_dxadxb_prev = gl.mix_type
    gl%ge%C_or_R_dxadxb_prev = C_or_R
    gl%ge%gE_C_dxadxb_prev = gl%ge%gE_C_dxadxb
    gl%ge%gE_R_dxadxb_prev = gl%ge%gE_R_dxadxb
    gl%ge%ln_gamma_C_dxadxb_prev = gl%ge%ln_gamma_C_dxadxb
    gl%ge%ln_gamma_R_dxadxb_prev = gl%ge%ln_gamma_R_dxadxb


    end subroutine
    !--------------------------------------------------------------------------------------
    !**************************************************************************
    subroutine db_HelmhE_dxi(gl,dbHelmgEdxi)
    !**************************************************************************
    !This subroutine calculates the partial derivative of the mixed parameter b
    !of the covolume for the gE based departure function with respect to the compositions xi
    !
    !OUTPUT: vector dbHelmgEdxi, containing the derivatives db / dxi at const. del, tau, xk





    implicit none

    type(type_gl) :: gl


    double precision, dimension(30):: dbHelmgEdxi
    integer:: i, j

    dbHelmgEdxi = 0.D0

    Do i = 1, gl%ncomp - 1
        dbHelmgEdxi(i) = gl%bi_HelmgE(i) - gl%bi_HelmgE(gl%ncomp)
    End do

    end subroutine db_HelmhE_dxi
    !**************************************************************************

    subroutine d2b_HelmhE_dxidxj(gl,d2bHelmgEdxidxj)
    !**************************************************************************
    !This subroutine calculates the second partial derivative of the mixed parameter b
    !of the covolume for the gE based departure function with respect to the compositions xi and xj
    !
    !OUTPUT: vector d2bHelmgEdxidxj, containing the derivatives d2b / dxidxj at const. del, tau, xk





    implicit none

    type(type_gl) :: gl


    double precision, dimension(30,30):: d2bHelmgEdxidxj
    integer:: i, j

    d2bHelmgEdxidxj = 0.D0  !All second derivates 0 at the moment because b is not a linear function of the mole fractions

    end subroutine d2b_HelmhE_dxidxj
    !**************************************************************************



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !**************************************************************************
    !   Andreas, November 2013
    !   Derivative of ln(fug_coef) w.r.t ni at constant T,p,nk implemented
    !   Not needed for our phase equilibrium routines, implemented for SINTEF
    !**************************************************************************
    subroutine ndlnphii_dnj(gl,Temperature, Density, dphidnj)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE FUGACITY COEFFICIENT
    ! OF COMPONENT i W.R.T. nj AT CONSTANT T AND P!!!
    !                            --- n*d(ln(phi_i)/d(nj)) ---
    ! A. Jäger , J. Gernert, Aug. 2010
    !--------------------------------------------------------------------------------------------------
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 118 EQ. 7.31)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! dphidnj     - A MATRIX OF(30 x 30)IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------




    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30,30):: dphidnj

    double precision, dimension(30):: dPdn
    double precision, dimension(30,30)::DERIVSANINJ
    double precision:: dPdV, Rmix, T
    integer::i, j

    T = Temperature
    ! call the derivative n*d(d(n*ar)/d(ni))dnj
    call nddna_dnidnj(gl,T, DENSITY, DERIVSANINJ)

    ! call n*d(P)/d(ni) for each component i
    call ndP_dni_TV (gl,T, DENSITY, dPdn)

    call R_mix_calc(gl,Rmix)
    ! call n*d(P)/d(V)
    call ndP_dV (gl,T, DENSITY, dPdV)

    do j = 1, gl%ncomp
        do i = j, gl%ncomp
            dphidnj(j, i) = 0.D0
            dphidnj(j, i) = DERIVSANINJ(j, i) + 1.D0 + dPdn(j)*dPdn(i)/(Rmix*T*dPdV)
            dphidnj(i, j) = dphidnj(j, i)   ! the matrix is symmetric ...
        end do
    end do

    end subroutine ndlnphii_dnj

    !**************************************************************************


    !**************************************************************************
    !   Andreas, November 2013
    !   Additional derivatives needed for the derivative ln(fug_coef) w.r.t ni
    !   at constant T,p,nk
    !**************************************************************************

    subroutine nddna_dnidnj(gl,TEMPERATURE, DENSITY, DERIVSANINJ)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF alpha_r W.R.T. nj AND
    ! ni - n*d(d(n*ar)/d(ni))dnj AT CONSTANT T and V
    ! A. Jäger , J. Gernert, Aug. 2010
    !--------------------------------------------------------------------------------------------------
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 120 EQ. 7.46)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVSANINJ     - A MATRIX OF(30 x 30)IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30,30):: DERIVSANINJ
    double precision, dimension(30)::DANIDEL, DANITAU, ndTred_dni, ndrhored_dni, CHEMPOT
    double precision, dimension(30,30)::dndadndx
    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivs):: getder_res
    double precision:: xk, nda_dnj, dna_dnj, ar
    integer::i, j, k, OIR

    DERIVSANINJ = 0.D0

    !Andreas, February 2016:
    !Old derivative calls changed to the actual ones. There is no difference, but the old ones are not maintained anymore
    !This derivative is not needed for flash calculations but will be needed for the critical point routines

    ! get the derivatives del*d(n* d(ar)/d(ni))/d(del) for all ni
    !call dndar_dniddel_old(TEMPERATURE, DENSITY, DANIDEL)
    call dndar_dniddel(gl,TEMPERATURE, DENSITY, DANIDEL)

    ! get the derivatives tau*d(n* d(ar)/d(ni))/d(tau) for all ni
    !call dndar_dnidtau_old (TEMPERATURE, DENSITY, DANITAU)
    call dndar_dnidtau (gl,TEMPERATURE, DENSITY, DANITAU)

    ! get the derivatives n*d(Tred)/d(ni) and n*d(rhored)/d(ni) for all ni
    !call ndYr_dni_old(ndTred_dni, ndrhored_dni)
    call ndYr_dni(gl,ndTred_dni, ndrhored_dni)

    ! get the derivatives d(n*d(ar)/d(ni))dxj for all xj, ni
    !call dnda_dnidxj_old(TEMPERATURE, DENSITY, dndadndx)
    call dnda_dnidxj(gl,TEMPERATURE, DENSITY, dndadndx)

    ! get the residual part of the cemical potential
    OIR = 2
    call dna_dni (gl,TEMPERATURE, DENSITY, CHEMPOT, OIR) !use the new function

    ! get alpha_r
    getder_res = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    Call MIXDERIVSFNR (gl,Temperature, Density, getder_res, der_res)
    ar = der_res(1)

    ! calculate the derivative n*d(d(n*ar)/d(ni))dnj [GERG-Monograph, p.120, Equ. 7.46]
    ! the triple loop includes the second part of Equ. 7.46, which is given in Equ. 7.47
    do j = 1, gl%ncomp
        dna_dnj = CHEMPOT(j)
        nda_dnj = dna_dnj - ar  !calculate n*d(ar)/d(nj)
        do i = 1, gl%ncomp
            DERIVSANINJ(j, i) = DANIDEL(i)*(1.D0 - ndrhored_dni(j)/gl%rhoredmix) &
                & + DANITAU(i)*ndTred_dni(j)/gl%tredmix + dndadndx(j, i) + nda_dnj
            do k = 1, gl%ncomp
                xk = gl%molfractions(k)
                DERIVSANINJ(j, i) = DERIVSANINJ(j, i) - xk*dndadndx(k, i)
            end do
        end do
    end do

    end subroutine nddna_dnidnj
    !**************************************************************************
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine dndar_dniddel_old(gl,TEMPERATURE, DENSITY, DERIVS)
    !**************************************************************************
    !subroutine for the calculation of the derivative of n*d(ar)/d(ni) with
    !respect to del, multiplied with del - del*d(n*d(ar)/d(ni))/d(del)
    !J. Gernert, Aug. 2010
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 121, Equ. 7.50)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    ! OUTPUT PARAMETERS:
    ! DERIVS      - A VECTOR (30)IN WHICH THE DERIVATIVE WITH RESPECT TO THE AMOUNT OF SUBSTANCE OF
    !               COMPONENT i AND DELTA IS STORED FOR EACH COMPONENT OF THE MIXTURE.
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVS
    double precision, dimension(30):: ndTred_dni, ndrhored_dni, DERIVFXDEL
    double precision, dimension(nderivs):: FNRDER
    integer, dimension(nderivs):: GETDER
    double precision:: dadd, d2add2, d2adddt, d2addeldxi, d2addeldxk, xk
    integer:: i, k

    ! get derivatives of alpha_r w.r.t. del and tau
    GETDER = (/0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR(gl,TEMPERATURE, DENSITY, GETDER, FNRDER)
    ! 1st der. of alpha_r w.r.t del, multiplied with del - del*d(ar)/d(del)
    dadd = FNRDER(2)
    ! 2nd der. of alpha_r w.r.t del, multiplied with del^2 - del^2*d^2(ar)/d(del)^2
    d2add2 = FNRDER(3)
    ! 1st mixed der. of alpha_r w.r.t. del and tau, multiplied with del and tau - del*tau*d^2(ar)/d(del)d(tau)
    d2adddt = FNRDER(6)

    ! get derivatives of the reducing parameters Tred and rhored w.r.t all ni
    call ndYr_dni_old(gl,ndTred_dni, ndrhored_dni)

    !get derivatives of alpha_r w.r.t del and all xi
    call d2ar_dxiddel_old (gl,TEMPERATURE, DENSITY, DERIVFXDEL)

    ! calculate the derivative
    do i = 1, gl%ncomp
        DERIVS(i)= 0.D0
        d2addeldxi = DERIVFXDEL(i)
        DERIVS(i) = ((dadd + d2add2)*(1.D0 - ndrhored_dni(i)/gl%rhoredmix) &
            & + d2adddt*ndTred_dni(i)/gl%tredmix + d2addeldxi)
        ! inner summation over all components ...
        do k = 1, gl%ncomp
            xk = gl%MOLFRACTIONS(k)
            d2addeldxk = DERIVFXDEL(k)
            DERIVS(i) = DERIVS(i) - xk*d2addeldxk
        end do
    end do

    end subroutine dndar_dniddel_old
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine dndar_dnidtau_old (gl,TEMPERATURE, DENSITY, DANITAU)
    !**************************************************************************
    ! A. Jäger, Aug. 2010
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF n* d(ar)/d(ni) WITH RESPECT TO tau, MULTIPLIED BY tau
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 121 EQ. 7.51)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A VECTOR (30)IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DANITAU

    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivs):: getder_res
    integer:: i,j
    double precision, dimension(30):: drhoredmix_dni, dtredmix_dni, dda_dxdt
    !dda_dtdd = second deriv of a with respect to del and tau
    !da_dt = first deriv of a with respect to tau
    !dda_dtdt = second deriv of a with respect to tau
    !dda_dxdt = second deriv of a with respect to tau and xi
    double precision:: dda_dtdd, da_dt, dda_dtdt, dda_dxidt, dda_dxjdt

    getder_res = (/0,0,0,1,1,1,0,0,0,0,0,0,0,0,0/)
    der_res = 0.d0

    Do i = 1, 30
        DANITAU(i) = 0.D0               !Initialize return vector
    End do

    dda_dtdd = 0.D0                     !Initialize
    da_dt  = 0.D0                       !Initialize
    dda_dtdt= 0.D0                      !Initialize

    call ndYr_dni_old(gl,dtredmix_dni, drhoredmix_dni)                     !Get derivatives of reduced parameters with respect to ni
    Call d2ar_dxidtau_old(gl,temperature, density, dda_dxdt)               !Get derivatives of a with respect to tau and xi
    Call MIXDERIVSFNR (gl,temperature, density, getder_res, der_res)   !Get derivatives of ar
    dda_dtdd = der_res(6)
    da_dt  = der_res(4)
    dda_dtdt= der_res(5)

    Do i = 1, gl%NCOMP                     !Loop over all components
        dda_dxidt = dda_dxdt(i)
        DANITAU(i) = (dda_dtdd*(1.d0-drhoredmix_dni(i)/gl%rhoredmix)+(da_dt+dda_dtdt)*dtredmix_dni(i)/gl%tredmix+dda_dxidt)
        Do j = 1, gl%NCOMP
            dda_dxjdt = dda_dxdt(j)
            DANITAU(i) = DANITAU(i) - gl%molfractions(j) * dda_dxjdt
        end do
    end do

    end subroutine dndar_dnidtau_old
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine d2ar_dxiddel_old (gl,TEMPERATURE, DENSITY, DERIVFXDEL)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE COMPOSITION Xi AND DELTA FOR EACH COMPONENT
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 111 table 7.5)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A MATRIX (30)IN WHICH THE DERIVATIVE WITH RESPECT TO COMPOSITION AND DEL FOR
    !              EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVFXDEL

    double precision, dimension(nderivs):: der_res, der_dep
    integer, dimension(nderivs):: getder_res, getder_dep
    integer:: i, j
    double precision:: d_alphaijr_d_delta, d_F0ir_d_delta
    double precision, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij
    logical, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij_calculated


    getder_res = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    getder_dep = (/0,1,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.d0
    der_dep = 0.d0

    Do i = 1, 30
        DERIVFXDEL(i) = 0.D0                        !Initialize return matrix DERIVFX
    End do
    d_alphaijr_d_delta = 0.D0                       !Initialize departure function derivative
    d_F0ir_d_delta  = 0.D0                          !Initialize residual helmholtz pure derivative

    DEPFUNCFNR_ij_calculated = .false.

    !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP and delta
    Do i = 1, gl%NCOMP                                 !Loop over all components
        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
        d_F0ir_d_delta  = der_res(2)                !Derivative with respect to delta of the residual part of F for component i
        DERIVFXDEL(i) = d_F0ir_d_delta
        Do j = 1, gl%NCOMP                             !Sum over all components with departure function
            if(gl%Fij(i,j)  /=  0.d0) then
                Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, j, 2, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                d_alphaijr_d_delta = der_dep(2)     !Derivative of the residual part of the departure function with respect to delta
                DERIVFXDEL(i) = DERIVFXDEL(i)+ gl%molfractions(j)* gl%Fij(i,j) * d_alphaijr_d_delta
            end if
        end do
    end do

    end subroutine d2ar_dxiddel_old
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine nd2ar_dnidtau (gl,TEMPERATURE, DENSITY, danitau)
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE RESIDUAL HELMHOLTZ ENERGY
    ! WITH RESPECT TO THE REDUCED TEMPERATURE TAU AND THE MOLAR NUMBER OF COMPONENT i IN THE
    ! MIXTURE. THE DERIVATIVE IS GIVEN TIMES THE TOTAL NUMBER OF MOLES n AND TAU. THE
    ! CALCULATION IS BASED ON THE FORMULATION OF THE HELMHOLTZ FREE ENERGY AS PUBLISHED BY:
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 123 EQ. (7.65))
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! Temperature
    ! Density
    !
    ! OUTPUT PARAMETERS:
    ! danitau -  vector (30) the derivative n * delta * d(Fr)_d(ni)_d(del) is stored in.
    !-------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: danitau

    integer, dimension(nderivs):: getder_res
    double precision, dimension(nderivs):: der_res
    double precision, dimension(30):: dda_dxdt

    !derivatives of the reducingfunctions with respect to ni multiplied by n
    double precision, dimension(30):: drhoredmix_dni, dtredmix_dni
    !second derivatives of Fr with respect to delta, tau, x and combinations of those are needed
    double precision:: dda_dtdd, dda_dtdt, dda_dtdxi, dda_dtdxj
    integer:: i, j

    getder_res = (/0,0,0,0,1,1,0,0,0,0,0,0,0,0,0/)                                !Derivatives with respect to del^2 and to del and tau are needed
    der_res = 0.d0
    dda_dxdt = 0.D0
    danitau = 0.D0

    call ndYr_dni_old(gl,dtredmix_dni, drhoredmix_dni)                             !Get derivatives of reduced parameters with respect to ni
    Call MIXDERIVSFNR (gl,Temperature, Density, getder_res, der_res)           !Second (mixed) derivatives of Fr with respect to delta and tau
    Call d2ar_dxidtau_old (gl,Temperature, Density, dda_dxdt)                      !Second Derivative of F_res with respect to xi and delta
    dda_dtdd = der_res(6)
    dda_dtdt = der_res(5)

    Do i = 1, gl%NCOMP
        dda_dtdxi = dda_dxdt(i)
        danitau(i) = dda_dtdd*(1.d0-1.d0/gl%rhoredmix*drhoredmix_dni(i))+ &
            & dda_dtdt/gl%tredmix*dtredmix_dni(i) + dda_dtdxi
        Do j = 1, gl%NCOMP
            dda_dtdxj = dda_dxdt(j)
            danitau(i) = danitau(i) - gl%molfractions(j) * dda_dtdxj
        end do
    end do

    end subroutine nd2ar_dnidtau
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine d2ar_dxidtau_old (gl,TEMPERATURE, DENSITY, DERIVFXTAU)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE COMPOSITION Xi AND TAU FOR EACH COMPONENT
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 111 table 7.5)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A MATRIX (30)IN WHICH THE DERIVATIVE WITH RESPECT TO COMPOSITION AND TAU FOR
    !              EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVFXTAU

    double precision, dimension(nderivs):: der_res, der_dep
    integer, dimension(nderivs):: getder_res, getder_dep
    integer:: i, j
    double precision:: d_alphaijr_d_tau, d_F0ir_d_tau
    double precision, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij
    logical, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij_calculated


    getder_res = (/0,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)
    getder_dep = (/0,0,0,1,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.d0
    der_dep = 0.d0

    Do i = 1, 30
        DERIVFXTAU(i) = 0.D0                !Initialize return matrix DERIVFX
    End do
    d_alphaijr_d_tau = 0.D0                 !Initialize departure function
    d_F0ir_d_tau  = 0.D0                    !Initialize residual helmholtz pure

    DEPFUNCFNR_ij_calculated = .false.

    !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP and tau
    Do i = 1, gl%NCOMP                        !Loop over all components
        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
        d_F0ir_d_tau  = der_res(4)          !Derivative with respect to tau of the residual part of F for component i
        DERIVFXTAU(i) = d_F0ir_d_tau
        Do j = 1, gl%NCOMP                         !Sum over all components with departure function
            if(gl%Fij(i,j)  /=  0.d0) then
                Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, j, 4, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                d_alphaijr_d_tau = der_dep(4)     !Derivative of the residual part of the departure function with respect to tau
                DERIVFXTAU(i) = DERIVFXTAU(i)+ gl%molfractions(j)* gl%Fij(i,j) * d_alphaijr_d_tau
            end if
        end do
    end do

    end subroutine d2ar_dxidtau_old
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine dnda_dnidxj_old(gl,TEMPERATURE, DENSITY, dndadndx)
    !**************************************************************************
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF alpha_r W.R.T. xj AND
    ! ni - d(n*d(ar)/d(ni))dxj
    ! A. Jäger , J. Gernert, Aug. 2010
    !--------------------------------------------------------------------------------------------------
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 121 EQ. 7.52)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! dndadndx     - A MATRIX OF(30 x 30)IN WHICH THE DERIVATIVE IS STORED
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30,30):: dndadndx
    integer, dimension(nderivs)::GETDER
    double precision, dimension(nderivs)::MIXDERFNR
    double precision, dimension(30)::DERIVFXDEL, DERIVFXTAU, ndTred_dni, ndrhored_dni
    double precision, dimension(30)::dTred_dxi, drhored_dxi
    double precision, dimension(30,30)::dndTr_dnidxj, dndrhor_dnidxj
    double precision, dimension(30):: DERIVFX
    double precision, dimension(30,30):: DERIVFXiXj
    double precision:: da_dd, da_dt, dda_dxjdd, dda_dxjdt, da_dxj, dda_dxidxj, xm, dda_dxjdxm!, dummy
    integer:: i, j, m

    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            dndadndx(i,j) = 0.D0
        end do
    end do

    !get the derivatives tau*d(ar)/d(tau) and del*d(ar)/d(del)
    GETDER = (/0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    call MIXDERIVSFNR (gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)

    ! get the derivatives del*d^2(ar)/d(xi)d(del) for all xi
    call d2ar_dxiddel_old (gl,TEMPERATURE, DENSITY, DERIVFXDEL)

    ! get the derivatives tau*d^2(ar)/d(xi)d(tau) for all xi
    call d2ar_dxidtau_old (gl,TEMPERATURE, DENSITY, DERIVFXTAU)

    ! get the derivatives n*d(Tred)/d(ni) and n*d(rhored)/d(ni) for all ni
    call ndYr_dni_old(gl,ndTred_dni, ndrhored_dni)

    ! get the derivatives d(Tred)/d(xi) and d(rhored)/d(xi) for all xi
    call dYr_dxi_old(gl,dTred_dxi, drhored_dxi)

    ! get the derivatives d(n*d(Tred)/d(ni))dxj and d(n*d(rhored)/d(ni))dxj for all ni, xj
    call dndYr_dnidxj_old(gl,dndTr_dnidxj, dndrhor_dnidxj)

    ! get the derivatives d(ar)/d(xi) for all xi
    call dar_dxi_old (gl,TEMPERATURE, DENSITY, DERIVFX)

    ! get the derivatives d^2(ar)/d(xi)d(xj) for all xi, xj
    call d2ar_dxidxj_old (gl,TEMPERATURE, DENSITY, DERIVFXiXj)

    da_dd = MIXDERFNR(2)
    da_dt = MIXDERFNR(4)
    do j = 1, gl%ncomp
        do i = 1, gl%ncomp
            dda_dxjdd = DERIVFXDEL(j)
            dda_dxjdt = DERIVFXTAU(j)
            da_dxj = DERIVFX(j)

            dndadndx(j, i) = dda_dxjdd*(1.D0 - ndrhored_dni(i)/gl%rhoredmix) + dda_dxjdt*ndTred_dni(i)/gl%tredmix &
                & - da_dd/gl%rhoredmix*(dndrhor_dnidxj(j, i) - ndrhored_dni(i)/gl%rhoredmix*drhored_dxi(j)) &
                & + da_dt/gl%tredmix*(dndtr_dnidxj(j, i) - ndtred_dni(i)/gl%tredmix*dtred_dxi(j)) &
                & - da_dxj
            if (j  /=  i) then
                dda_dxidxj = DERIVFXiXj(j, i)
                dndadndx(j, i) = dndadndx(j, i) + dda_dxidxj
            end if
            !dummy = 0.D0
            do m = 1, gl%ncomp
                xm = gl%molfractions(m)
                if (j  /=  m) then
                    dda_dxjdxm = DERIVFXiXj(j, m)
                    !dummy = dummy + xm*dda_dxjdxm
                    dndadndx(j, i) = dndadndx(j, i) - xm*dda_dxjdxm
                end if
            end do
        end do
    end do

    end subroutine dnda_dnidxj_old
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine dar_dxi_old (gl,TEMPERATURE, DENSITY, DERIVFX)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE MOLEFRACTION Xi FOR EACH COMPONENT.
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 111 table 7.5)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A MATRIX (30)IN WHICH THE DERIVATIVE WITH RESPECT TO COMPOSITION FOR
    !               EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------







    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVFX

    double precision, dimension(nderivs):: der_res, der_dep
    integer, dimension(nderivs):: getder_res, getder_dep
    integer:: i, j
    double precision:: alphaij_r, F0i_r
    double precision, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij
    logical, dimension (gl%ncomp,gl%ncomp) :: DEPFUNCFNR_ij_calculated


    getder_res = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    getder_dep = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    der_res = 0.d0
    der_dep = 0.d0
    DERIVFX = 0.D0
    alphaij_r = 0.D0                    !Initialize departure function
    F0i_r = 0.D0                        !Initialize residual helmholtz pure

    DEPFUNCFNR_ij_calculated = .false.

    !Derivative of the residual Helmholtz energy with respect to all molefractions x1, x2, x3,...,xNCOMP
    Do i = 1, gl%NCOMP                             !Loop over all components
        Call FNRDERIVS (gl,temperature, density, getder_res, der_res, i)
        F0i_r = der_res(1)
        DERIVFX(i) = F0i_r                      !residual part of component i: F0i_r
        Do j = 1, gl%NCOMP                         !Sum over all components with departure function
            if(gl%Fij(i,j)  /=  0.d0) then
                Call DEPFUNCFNR_caller (gl,temperature, density, getder_dep, der_dep, i, j, 1, DEPFUNCFNR_ij, DEPFUNCFNR_ij_calculated)
                alphaij_r = der_dep(1)
                DERIVFX(i) = DERIVFX(i)+ gl%molfractions(j)*gl%Fij(i,j) * alphaij_r
            end if
        end do
    end do

    end subroutine dar_dxi_old
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine d2ar_dxi2_old (gl,TEMPERATURE, DENSITY, DERIVFX2)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE MOLEFRACTION Xi FOR EACH COMPONENT.
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 111 table 7.5)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX2    - A MATRIX (30)IN WHICH THE SECOND DERIVATIVE WITH RESPECT TO COMPOSITION FOR
    !               EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------


    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: DERIVFX2

    integer::i

    Do i = 1, 30
        DERIVFX2(i) = 0.D0               !The second derivative of Fr with respect to composition is 0
    End do

    end subroutine d2ar_dxi2_old
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger
    !**************************************************************************
    subroutine d2ar_dxidxj_old (gl,TEMPERATURE, DENSITY, DERIVFXiXj)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE MIXED SECOND DERIVATIVE OF THE RESIDUAL PART OF THE
    ! HELMHOLTZ FREE ENERGY WITH RESPECT TO THE MOLEFRACTION Xi and Xj FOR EACH COMPONENT.
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 111 table 7.5)
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   mol/m^3
    !
    ! OUTPUT PARAMETERS:
    ! DERIVFX     - A MATRIX (30,30)IN WHICH THE MIXED DERIVATIVE WITH RESPECT TO MOLFRACTION i and j FOR
    !              EACH COMPONENT OF THE MIXTURE IS STORED.
    !--------------------------------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30,30):: DERIVFXiXj

    double precision, dimension(nderivs):: der_dep
    integer, dimension(nderivs):: getder_dep
    integer:: i, j
    double precision:: Fij_r

    getder_dep = (/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
    der_dep = 0.D0

    Do i = 1, 30
        Do j = 1, 30
            DERIVFXiXj(i,j) = 0.D0      !Initialize return matrix DERIVFXiXj
        End do
    End do
    Fij_r = 0.D0                        !Initialize departure function

    Do i = 1, (gl%NCOMP - 1)
        Do j = (i+1), gl%NCOMP
            if (gl%Fij(i,j)  /=  0.d0) then
                Call DEPFUNCFNR (gl,temperature, density, getder_dep, der_dep, i,j)
                Fij_r = der_dep(1)
                DERIVFXiXj(i,j) = gl%Fij(i,j) * Fij_r
                DERIVFXiXj(j,i) = gl%Fij(j,i) * Fij_r
            end if
        end do
    end do


    end subroutine d2ar_dxidxj_old
    !**************************************************************************


    !**************************************************************************
    !   Andreas Jäger (new function is used)
    !**************************************************************************
    subroutine dna_dni_old (gl,TEMPERATURE, DENSITY, CHEMPOT, OIR)
    ! SUBROUTINE FOR THE CALCULATION OF THE REDUCED CHEMICAL POTENTIALS. THE CALCULATION
    ! IS BASED ON THE FORMULATION OF THE HELMHOLTZ FREE ENERGY AS PUBLISHED BY:
    !                   Kunz, O. et al.
    !                   The Gerg-2004 Wide-Range Equation of State for Natural Gases and Other Mixtures
    !                   GERG TM15, 2007
    !                   (page 58 Eq. 5.36, page 109 Table 7.4)
    !-------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! Temperature
    ! Density
    ! OIR       - 0: Overall Chemical potential is given    d(nF)/d(ni)
    !           - 1: Ideal chemical potential is given      d(nFi)/d(ni)
    !           - 2: Residual chemical potential is given   d(nFr)/d(ni)
    !
    ! OUTPUT PARAMETERS:
    ! CHEMPOT - n x 1 vector where n is the number of components in the mixture
    !-------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision:: TEMPERATURE, DENSITY
    double precision, dimension(30):: CHEMPOT
    integer:: OIR

    double precision, dimension(30):: chempot_ideal, chempot_res
    integer, dimension(nderivs):: getder_res
    double precision, dimension(nderivs):: der_res
    integer, dimension(nderivsi):: getder_ideal
    double precision, dimension(nderivsi):: der_ideal
    double precision, dimension(30):: der_FRX

    !derivatives of the reducingfunctions with respect to ni multiplied by n
    double precision, dimension(30):: d_rhoredmix_d_ni, d_tredmix_d_ni
    double precision:: F0i_0, F_res, F_res_d_delta, F_res_d_tau, F_res_d_xi, F_res_d_x
    integer:: i, j


    getder_res = (/1,1,0,1,0,0,0,0,0,0,0,0,0,0,0/)                                !F_res, d_F_res_d_delta and d_F_res_d_tau are needed
    getder_ideal = (/1,0,0,0,0,0,0,0,0,0/)                                  !F_0 for each pure component is needed

    der_res = 0.d0
    der_ideal = 0.d0

    der_FRX = 0.D0
    CHEMPOT = 0.D0
    chempot_ideal = 0.D0
    chempot_res = 0.D0

    if (gl%NCOMP  ==  1) Then                                                    !PURE COMPONENT
        Call FNIDERIVS (gl,Temperature, Density, getder_ideal, der_ideal,1)    !Get F of of the pure component
        F0i_0 = der_ideal(1)
        chempot_ideal(1) = 1.D0 + F0i_0                                        !Ideal chemical potential for pure component
        Call FNRDERIVS (gl,Temperature, Density, getder_res, der_res, 1)
        F_res = der_res(1)
        F_res_d_delta = der_res(2)
        chempot_res(1) = F_res + F_res_d_delta                              !Residual chemical potential for pure component

        select case (OIR)
        case (0)
            chempot(1) = chempot_ideal(1) + chempot_res(1)              !Overall chemical potential
        case (1)
            chempot(1) = chempot_ideal(1)                               !Ideal chemical potential
        case (2)
            chempot(1) = chempot_res(1)                                 !Residual chemical potential
        end select
    else                                                                    !MIXTURE
        Do i = 1, gl%NCOMP
            Call FNIDERIVS (gl,Temperature, Density, getder_ideal, der_ideal,i)    !Get F of the pure component
            F0i_0 = der_ideal(1)
            if (gl%molfractions(i)  >  0.D0) then
                chempot_ideal(i) = 1.d0 + F0i_0 + dlog(gl%molfractions(i))             !Ideal chem.pot of comp i = F0i + ln(xi)
            end if
        end do

        call ndYr_dni(gl,d_tredmix_d_ni, d_rhoredmix_d_ni)                       !Get derivatives of reduced parameters with respect to ni
        Call MIXDERIVSFNR (gl,Temperature, Density, getder_res, der_res)           !F_res and derivatives with respect to delta and tau
        Call dar_dxi_old (gl,Temperature, Density, der_FRX)                             !Derivative of F_res with respect to xi
        F_res = der_res(1)
        F_res_d_delta = der_res(2)
        F_res_d_tau = der_res(4)

        Do i = 1, gl%NCOMP
            F_res_d_xi = der_FRX(i)
            chempot_res(i) = F_res + F_res_d_delta*(1.d0-1.d0/gl%rhoredmix*d_rhoredmix_d_ni(i))+ &
                & F_res_d_tau/gl%tredmix*d_tredmix_d_ni(i) + F_res_d_xi
            Do j = 1, gl%NCOMP
                F_res_d_x = der_FRX(j)
                chempot_res(i) = chempot_res(i) - gl%molfractions(j) * F_res_d_x
            end do
        end do

        select case (OIR)
        case (0)
            Do i = 1, gl%NCOMP
                chempot(i) = chempot_ideal(i) + chempot_res(i)                  !Overall chemical potential
            end do
        case (1)
            Do i = 1, gl%NCOMP
                chempot(i) = chempot_ideal(i)                                   !Ideal chemical potential
            end do
        case (2)
            Do i = 1, gl%NCOMP
                chempot(i) = chempot_res(i)                                     !Residual chemical potential
            end do
        end select

    end if

    end subroutine dna_dni_old
    !**************************************************************************

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    !**************************************************************************
    subroutine dP_dXi_taudel(gl,TEMPERATURE, DENSITY, dPdXi)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE PRESSURE WITH RESPECT TO THE MOL FRACTION
    ! Xi AT CONSTANT TAU AND DELTA
    !                                           d(p)/d(xi) at const. TAU, DELTA, Xj
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOLE/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dPdXi       - Vector(30) the derivatives for each component are stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision, intent(in):: TEMPERATURE, DENSITY
    double precision, dimension(30), intent(out):: dPdXi
    double precision, dimension(30):: d2ar__dxiddel, dTreddxi, drhoreddxi
    double precision, dimension(15):: MIXDERFNR
    double precision:: Rmix, dar_ddel, tau, delta
    integer:: i
    integer, dimension(15):: GETDER

    ! get the 1st derivative of alpha-r w.r.t. delta
    GETDER = 0
    GETDER(2) = 1
    call MIXDERIVSFNR(gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)
    dar_ddel = MIXDERFNR(2)

    ! get the derivative del*d^2(ar)/d(xi)d(del) at const. TAU, DELTA
    call d2ar_dxiddel(gl,TEMPERATURE, DENSITY, d2ar__dxiddel)

    ! get the 1st derivative of the reducing functions w.r.t xi
    call dYr_dxi(gl,dTreddxi, drhoreddxi)

    ! get the mixture Molar Gas Constant
    call R_mix_calc(gl,Rmix)

    tau = gl%tredmix/TEMPERATURE
    delta = DENSITY/gl%rhoredmix

    do i = 1, gl%ncomp-1
        dPdXi(i) = (delta*Rmix)/tau*(drhoreddxi(i)*gl%tredmix+gl%rhoredmix*dTreddxi(i))+(delta*Rmix)/tau*(drhoreddxi(i)*gl%tredmix*dar_ddel+gl%rhoredmix*dTreddxi(i)*dar_ddel+gl%rhoredmix*gl%tredmix*d2ar__dxiddel(i))
    end do

    end subroutine

    !**************************************************************************
    subroutine dP_dXi_Trho(gl,TEMPERATURE, DENSITY, dPdXi)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE PRESSURE WITH RESPECT TO THE MOL FRACTION
    ! Xi AT CONSTANT TEMPERATURE AND DENSITY
    !                                           d(p)/d(xi) at const. T, rho, Xj
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOLE/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dPdXi       - Vector(30) the derivatives for each component are stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision, intent(in):: TEMPERATURE, DENSITY
    double precision, dimension(30), intent(out):: dPdXi
    double precision, dimension(15):: MIXDERFNR
    double precision, dimension(30):: dTreddxi, drhoreddxi, dP_dXi__taudel
    double precision:: Rmix, dar_ddel, d2ar_d2del, d2ar_ddeldtau, tau, delta
    integer:: i
    integer, dimension(15):: GETDER

    ! get the derivatives of alpha-r
    GETDER = 0
    GETDER(2) = 1
    GETDER(3) = 1
    GETDER(6) = 1
    call MIXDERIVSFNR(gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)
    dar_ddel = MIXDERFNR(2)
    d2ar_d2del = MIXDERFNR(3)
    d2ar_ddeldtau = MIXDERFNR(6)

    ! get the 1st derivative of the reducing functions w.r.t xi
    call dYr_dxi(gl,dTreddxi, drhoreddxi)

    ! get the mixture Molar Gas Constant
    call R_mix_calc(gl,Rmix)

    ! get the 1st derivative of the pressure p w.r.t xi at const. TAU, DELTA, xj
    call dP_dXi_taudel(gl,TEMPERATURE, DENSITY, dP_dXi__taudel)

    tau = gl%tredmix/TEMPERATURE
    delta = DENSITY/gl%rhoredmix

    do i = 1, gl%ncomp-1
        dPdXi(i) = DENSITY*Rmix*dTreddxi(i)*((d2ar_ddeldtau/tau)-1.D0/tau*(1.D0+dar_ddel)) - Rmix*TEMPERATURE*delta*drhoreddxi(i)*(1.D0+2.D0*dar_ddel+d2ar_d2del) + dP_dXi__taudel(i)
    end do

    end subroutine

    !**************************************************************************
    subroutine dlnfi_drho_TX(gl,TEMPERATURE, DENSITY, dlnfidrho)

    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF ln(fi) WITH RESPECT TO THE DENSITY rho
    ! AT CONSTANT TEMPERATURE AND MOL FRACTION X
    !                                           d(ln(fi))/d(rho) at const. T, X
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOLE/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dPdXi       - Vector(30) the derivatives for each component are stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision, intent(in):: TEMPERATURE, DENSITY
    double precision, dimension(30), intent(out):: dlnfidrho
    double precision, dimension(15):: MIXDERFNR
    double precision, dimension(30):: DANIDEL
    double precision:: Rmix, dar_ddel, tau, delta
    integer:: i
    integer, dimension(15):: GETDER

    ! get the 1st derivative of alpha-r w.r.t. delta
    GETDER = 0
    GETDER(2) = 1
    call MIXDERIVSFNR(gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)
    dar_ddel = MIXDERFNR(2)

    ! get the mixture Molar Gas Constant
    call R_mix_calc(gl,Rmix)

    ! get del*d(n*d(ar)/d(ni))/d(del)
    call dndar_dniddel(gl,TEMPERATURE, DENSITY, DANIDEL)

    tau = gl%tredmix/TEMPERATURE
    delta = DENSITY/gl%rhoredmix

    do i = 1, gl%ncomp
        dlnfidrho(i) = 1/DENSITY * (1 + dar_ddel + DANIDEL(i))
    end do

    end subroutine

    !**************************************************************************
    subroutine dlnfi_dT_rhoX(gl,TEMPERATURE, DENSITY, dlnfidT)

    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF ln(fi) WITH RESPECT TO THE TEMPERATURE T
    ! AT CONSTANT DENSITY rho AND MOL FRACTION X
    !                                           d(ln(fi))/d(T) at const. rho, X
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOLE/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dPdXi       - Vector(30) the derivatives for each component are stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision, intent(in):: TEMPERATURE, DENSITY
    double precision, dimension(30), intent(out):: dlnfidT
    double precision, dimension(15):: MIXDERFNR
    double precision, dimension(30):: DANITAU
    double precision:: Rmix, dar_dtau, tau, delta
    integer:: i
    integer, dimension(15):: GETDER

    ! get the 1st derivative of alpha-r w.r.t. tau
    GETDER = 0
    GETDER(4) = 1
    call MIXDERIVSFNR(gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)
    dar_dtau = MIXDERFNR(4)

    ! get the mixture Molar Gas Constant
    call R_mix_calc(gl,Rmix)

    ! get d/dtau*(n*dalpha-r/dni) at const. delta, X
    call dndar_dnidtau (gl,TEMPERATURE, DENSITY, DANITAU)

    do i = 1, gl%ncomp
        dlnfidT(i) = 1.D0/TEMPERATURE * (1.D0 - dar_dtau - DANITAU(i))
    end do

    end subroutine

    !**************************************************************************
    subroutine dlnfi_dXi_TrhoXj(gl,TEMPERATURE, DENSITY, dlnfidXi)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF ln(fi) WITH RESPECT TO THE MOLE FRACTION X
    ! AT CONSTANT TEMPERATURE AND DENSITY rho
    !                                           d(ln(fi))/d(Xi) at const. T, rho, xj
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOLE/M^3
    !
    ! OUTPUT PARAMETERS:
    ! dPdXi       - Vector(30) the derivatives for each component are stored in
    !--------------------------------------------------------------------------------------------------





    implicit none

    type(type_gl) :: gl


    double precision, intent(in):: TEMPERATURE, DENSITY
    double precision, dimension(30,30), intent(out):: dlnfidXi
    double precision, dimension(15):: MIXDERFNR
    double precision, dimension(30):: DANITAU, DERIVFX, DANIDEL, dTreddxi, drhoreddxi
    double precision, dimension(30,30):: dndadndx
    double precision:: Rmix, dar_dtau, tau, delta, dar_ddel
    integer:: i, j
    integer, dimension(15):: GETDER

    ! get del*d(n*d(ar)/d(ni))/d(del)
    call dndar_dniddel(gl,TEMPERATURE, DENSITY, DANIDEL)

    ! get the 1st derivative of the reducing functions w.r.t xi
    call dYr_dxi(gl,dTreddxi, drhoreddxi)

    ! get the derivative of alpha-r
    GETDER = 0
    GETDER(2) = 1
    GETDER(4) = 1
    call MIXDERIVSFNR(gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)
    dar_ddel = MIXDERFNR(2)
    dar_dtau = MIXDERFNR(4)

    ! get d/dtau*(n*dalpha-r/dni) at const. delta, X
    call dndar_dnidtau (gl,TEMPERATURE, DENSITY, DANITAU)

    ! get d(n*dar/dni)dXi at const. tau, delta
    call dnda_dnidxj(gl,TEMPERATURE, DENSITY, dndadndx)

    ! get dardxi at const. tau, delta
    call dar_dxi (gl,TEMPERATURE, DENSITY, DERIVFX)

    tau = gl%tredmix/TEMPERATURE
    delta = DENSITY/gl%rhoredmix

    do j = 1, gl%ncomp-1
        do i = 1, gl%ncomp-1
            if (i == j) then
                dlnfidXi(j, i) = (-1/tau+dar_dtau/tau+DANITAU(i)/tau)*1/TEMPERATURE*dTreddxi(j) + (1 + dar_ddel + DANIDEL(i))*((-1/gl%rhoredmix)*drhoreddxi(j)) + 1/(gl%rhoredmix*gl%tredmix) * (gl%rhoredmix*gl%tredmix/gl%molfractions(j)+drhoreddxi(j)*gl%tredmix+gl%rhoredmix*dTreddxi(j)) + DERIVFX(j) + dndadndx(j,i)
            else
                dlnfidXi(j, i) = (-1/tau+dar_dtau/tau+DANITAU(i)/tau)*1/TEMPERATURE*dTreddxi(j) + (1 + dar_ddel + DANIDEL(i))*((-1/gl%rhoredmix)*drhoreddxi(j)) + 1/(gl%rhoredmix*gl%tredmix) * (drhoreddxi(j)*gl%tredmix+gl%rhoredmix*dTreddxi(j)) + DERIVFX(j) + dndadndx(j,i)
            end if
        end do
        dlnfidXi(j, gl%ncomp) = (-1/tau+dar_dtau/tau+DANITAU(gl%ncomp)/tau)*1/TEMPERATURE*dTreddxi(j) + (1 + dar_ddel + DANIDEL(gl%ncomp))*((-1/gl%rhoredmix)*drhoreddxi(j)) + 1/(gl%rhoredmix*gl%tredmix) * (-gl%rhoredmix*gl%tredmix/gl%molfractions(gl%ncomp)+drhoreddxi(j)*gl%tredmix+gl%rhoredmix*dTreddxi(j)) + DERIVFX(j) + dndadndx(j,gl%ncomp)
    end do

    end subroutine
    !***********************************************************************************


    !************************************************************************************
    !************************************************************************************
    subroutine chempot_reac_calc(gl, temp, dens, p, chempot, phase)
    !BS 09/2018, seawater


    implicit none

    type(type_gl) :: gl


    double precision, intent(in):: temp, dens, p
    double precision, dimension(30), intent(out):: chempot
    double precision, dimension(30) ::molfractions, molfractions_orig, a_up,a_low
    double precision:: Rmix, delta_num, tredmix_orig, rhoredmix_orig, factorx, xfactor_up, xfactor_low, dens_est, dens1, dens2, molup, mollow
    integer:: i, j, k,phase, oir !Phase=1=liquid, phase=2=vapor, necessary for seawater/reactive/gE additive calc flac

    !Save the original composition and reducing parameters
    molfractions_orig = gl%molfractions
    tredmix_orig = gl%tredmix
    rhoredmix_orig = gl%rhoredmix

    !a_up = 0.d0
    !a_low = 0.d0
    !chempot = 0.d0
    !molfractions = 0.d0
    !dens_est = dens

    call chempot_calc(gl, temp, dens, chempot, oir)       ! gives back the reactive part of the chemical potential occording to the reative model. Anything futher will happen in that routine, e.g. choosing the right models

    if(gl%seacalc) then
        chempot = chempot
        !else
        !    call chempot_reaction(gl, temp, d, p, phase)
        !    do i=1,gl%ncomp
        !        do j=1,gl%nreac
        !            if(i .eq. reacpos(j)) then
        !                chempot(i) = chempot(i) + chempot_reac(i)
        !            else
        !                chempot(i) = chempot(i)
        !
        !            end if
        !        end do
        !    end do
    end if



    !Write the original composition back
    !gl%molfractions = molfractions_orig
    !gl%tredmix = tredmix_orig
    !gl%rhoredmix = rhoredmix_orig

    end subroutine chempot_reac_calc
    !*******************************************************************************************
    !*******************************************************************************************


    !*******************************************************************************************
    !*******************************************************************************************
    subroutine dchempot_dxj_reac_tp(gl, temp, Dens_est, p, dchempot_dxi, phase)
    !BS, 09/2018, seawater

    implicit none

    type(type_gl) :: gl


    double precision, intent(in):: temp, dens_est, p
    double precision, dimension(30):: chempot1, chempot2
    double precision, dimension(30,30), intent(out) :: dchempot_dxi
    double precision, dimension(30,30) ::  dchempotdx, dchemp_reac_dx
    double precision, dimension(30) ::molfractions, molfractions_orig
    double precision:: Rmix, delta_num, tredmix_orig, rhoredmix_orig, factorx, xfactor_up, xfactor_low, dens1, dens2
    integer:: i, j, k,phase, reacpos !Phase=1=liquid, phase=2=vapor, necessary for seawater calc flac

    !Save the original composition and reducing parameters
    molfractions_orig = gl%molfractions
    tredmix_orig = gl%tredmix
    rhoredmix_orig = gl%rhoredmix

    !It is:
    !dlnfidXj_num = (lnfi(T,p,xj_up) - lnfi(T,p,xj_low)) / (xj_up - xj_low)
    factorx = 1.d-7!0.00001D0
    xfactor_up = 1.D0 + factorx
    xfactor_low = 1.D0 - factorx
    chempot1 = 0.d0
    dchempotdx= 0.d0
    dchempot_dxi = 0.d0
    molfractions = 0.d0
    dchemp_reac_dx = 0.d0

    call R_mix_calc(gl,Rmix)
    !call molswitch(gl, gl%molfraction, xp, xm)

    !starting numdiff
    !    open(unit=57, file='C:\Users\semrau\Documents\work\seawater\chempot.txt',status='old',position='append',  action='write')
    !2369 format(e21.14, 6(e21.14))

    call d2na_dnidxj_PT(gl, temp, dens_est, dchempotdx, 0)

    if((gl%seacalc) .and. (phase .eq.1)) then
        call d2na_dnidxj_PT_reac(gl, temp,p, dens_est, dchemp_reac_dx, reacpos)
        dchempot_dxi =( dchempotdx * Rmix * Temp / gl%wm(1) + dchemp_reac_dx)  *gl%sea%wm_sea
    else
        do i = 1,gl%ncomp-1
            dchempot_dxi(:,i) = dchempotdx(:,i) * Rmix * Temp !/ gl%sea%wm_sea
        end do
    end if

    ! write(57,2369) dchempotdx(1,1), dchemp_reac_dx(1,1),  dchempotdx(2,1), dchemp_reac_dx(2,1), dchempot_dxi(1,1)
    ! close(57)
    !do j = 1, gl%ncomp-1            !copied form dlnfi_dxj_TP_num
    !
    !    ! increase xj, xN stays constant
    !    molfractions(j) = molfractions_orig(j) * xfactor_up
    !    !Calculate xN
    !    molfractions(gl%ncomp) = 1.D0
    !    do k = 1, gl%ncomp - 1
    !        molfractions(gl%ncomp) = molfractions(gl%ncomp) - molfractions(k)
    !    end do
    !    gl%molfractions = molfractions
    !    call reduced_parameters_calc(gl, Temp)
    !    dens1 = rhomix_calc(gl, Temp, p, Dens_est, 0, 0)
    !    !call chempot_num_calc(gl, temp, dens1, p, chempot1, phase)       !brauche ich hier p?!
    !
    !
    !    ! decrease xj, xN stays constant
    !    molfractions(j) = molfractions_orig(j) * xfactor_low
    !    !Calculate xN
    !    molfractions(gl%ncomp) = 1.D0
    !    do k = 1, gl%ncomp - 1
    !        molfractions(gl%ncomp) = molfractions(gl%ncomp) - molfractions(k)
    !    end do
    !    gl%molfractions = molfractions
    !    call reduced_parameters_calc(gl,Temp)
    !    dens2 = rhomix_calc(gl, Temp, p, Dens_est, 0, 0)
    !    ! call chempot_num_calc(gl, temp, dens2, p, chempot2, phase)       !brauche ich hier p?!
    !
    !    do i = 1,gl%ncomp
    !        dchempot_dxi(j,i) = (chempot1(i) - chempot2(i)) / (molfractions_orig(j)*2.D0*factorx)
    !    end do
    !
    !    molfractions(j) = molfractions_orig(j)
    !    molfractions(gl%ncomp) = molfractions_orig(gl%ncomp)
    !
    !end do

    !Write the original composition back
    gl%molfractions = molfractions_orig
    gl%tredmix = tredmix_orig
    gl%rhoredmix = rhoredmix_orig

    end subroutine dchempot_dxj_reac_tp!calc
    !*******************************************************************************************
    !*******************************************************************************************





    !  !*******************************************************************************************
    !  !*******************************************************************************************
    !  subroutine dchempot_dxi_num_calc_false(gl, temp, Dens_est, p, dchempot_dxi, phase)
    !  !BS, 09/2018, seawater


    !  implicit none
    !
    !  type(type_gl) :: gl
    !
    !
    !  double precision, intent(in):: temp, dens_est, p
    !  double precision, dimension(30):: chempot1, chempot2
    !  double precision, dimension(30,30), intent(out) :: dchempot_dxi
    !  double precision, dimension(30) ::molfractions, molfractions_orig, aup1, aup2, alow1, alow2, a_norm
    !  double precision:: Rmix, delta_num, tredmix_orig, rhoredmix_orig, factorx, xfactor_up, xfactor_low, dens1, dens2, dens12, dens22, mollowdiff, molupdiff, molup2, mollow2, molup, mol
    !  integer:: i, j, k,phase !Phase=1=liquid, phase=2=vapor, necessary for seawater calc flac
    !
    !  !Save the original composition and reducing parameters
    !  molfractions_orig = gl%molfractions
    !  tredmix_orig = gl%tredmix
    !  rhoredmix_orig = gl%rhoredmix
    !
    !  !It is:
    !  !dlnfidXj_num = (lnfi(T,p,xj_up) - lnfi(T,p,xj_low)) / (xj_up - xj_low)
    !  factorx = 1.d-6!0.00001D0
    !  xfactor_up = 1.D0 + factorx
    !  xfactor_low = 1.D0 - factorx
    !  chempot1 = 0.d0
    !  chempot2= 0.d0
    !  dchempot_dxi = 0.d0
    !  molfractions = 0.d0
    !  aup1 =0.d0
    !  aup2=0.d0
    !  alow=0.d0
    !  alow2=0.d0
    !  a_norm=0.d0
    !
    !  !call molswitch(gl, gl%molfraction, xp, xm)
    !
    !  !starting numdiff
    !
    !  do j = 1, gl%ncomp-1            !copied form dlnfi_dxj_TP_num
    !
    !      a_norm(j) = a_calc(gl, temp, dens_est, 0)
    !
    !      ! increase xj, xN stays constant
    !      molfractions(j) = molfractions_orig(j) * xfactor_up
    !      molupdiff = molfractions(j)
    !      !Calculate xN
    !      !molfractions(gl%ncomp) = 1.D0
    !      do k = 1, gl%ncomp !- 1
    !          !molfractions(gl%ncomp) = molfractions(gl%ncomp) - molfractions(k)
    !
    !          if(k .ne. j) then
    !              if(gl%ncomp .eq. 2)then
    !                  molfractions(k) = 1.d0 - molfractions(j)
    !              else
    !                  !molfractions(k) = molfractions(k) - ( ( molup - molfractions_orig(j) ) / ( 1.d0 - molfractions_orig(j) ) * molfractions_orig(k) )
    !              end if
    !          end if
    !      end do
    !      gl%molfractions = molfractions
    !      call reduced_parameters_calc(gl, Temp)
    !      dens1 = rhomix_calc(gl, Temp, p, Dens_est, 0, 0)
    !
    !      !if(phase==1) gl%seacalc = .true.
    !      aup1(j) = A_CALC(gl, Temp, dens1, 0)!/gl%req(j) / Temp
    !      gl%seacalc = .false.
    !
    !      ! decrease xj, xN stays constant
    !      molfractions(j) = molfractions_orig(j) * xfactor_low
    !      mollowdiff = molfractions(j)
    !      !Calculate xN
    !      !molfractions(gl%ncomp) = 1.D0
    !
    !      do k = 1, gl%ncomp !- 1
    !          !molfractions(gl%ncomp) = molfractions(gl%ncomp) - molfractions(k)
    !          if(k .ne. j) then
    !              if(gl%ncomp .eq. 2)then
    !                  molfractions(k) = 1.d0 - molfractions(j)
    !              else
    !                  !molfractions(k) = molfractions(k) - ( ( molup - molfractions_orig(j) ) / ( 1.d0 - molfractions_orig(j) ) * molfractions_orig(k) )
    !              end if
    !          end if
    !          !end do
    !      end do
    !      gl%molfractions = molfractions
    !      call reduced_parameters_calc(gl, temp)
    !      dens2 = rhomix_calc(gl, temp, p, Dens_est, 0, 0)
    !
    !      !if(phase==1) gl%seacalc = .true.
    !      alow1(j) = a_calc(gl, temp, dens2, 0)!/gl%req(j) / Temp
    !      gl%seacalc = .false.
    !
    !      !-------------------------------------------- lijke first oder
    !
    !
    !      molfractions(j) = molfractions_orig(j) + (2.d0*molupdiff)
    !      molup2 = molfractions(j)
    !      !Calculate xN
    !      !molfractions(gl%ncomp) = 1.D0
    !      do k = 1, gl%ncomp !- 1
    !          !molfractions(gl%ncomp) = molfractions(gl%ncomp) - molfractions(k)
    !
    !          if(k .ne. j) then
    !              if(gl%ncomp .eq. 2)then
    !                  molfractions(k) = 1.d0 - molfractions(j)
    !              else
    !                  molfractions(k) = molfractions(k) - ( ( molup2 - molfractions_orig(j) ) / ( 1.d0 - molfractions_orig(j) ) * molfractions_orig(k) )
    !              end if
    !          end if
    !      end do
    !      gl%molfractions = molfractions
    !      call reduced_parameters_calc(gl, Temp)
    !      dens12 = rhomix_calc(gl, Temp, p, Dens_est, 0, 0)
    !
    !      !if(phase==1) gl%seacalc = .true.
    !      aup2(j) = A_CALC(gl, Temp, dens12, 0)!/gl%req(j) / Temp
    !      gl%seacalc = .false.
    !
    !      ! decrease xj, xN stays constant
    !      molfractions(j) = molfractions_orig(j) - (2.d0*mollowdiff)
    !      mollow2 = molfractions(j)
    !      !Calculate xN
    !      !molfractions(gl%ncomp) = 1.D0
    !
    !      do k = 1, gl%ncomp !- 1
    !          !molfractions(gl%ncomp) = molfractions(gl%ncomp) - molfractions(k)
    !          if(k .ne. j) then
    !              if(gl%ncomp .eq. 2)then
    !                  molfractions(k) = 1.d0 - molfractions(j)
    !              else
    !                  molfractions(k) = molfractions(k) - ( ( molup - molfractions_orig(j) ) / ( 1.d0 - molfractions_orig(j) ) * molfractions_orig(k) )
    !              end if
    !          end if
    !          !end do
    !      end do
    !      gl%molfractions = molfractions
    !      call reduced_parameters_calc(gl, temp)
    !      dens22 = rhomix_calc(gl, temp, p, Dens_est, 0, 0)
    !
    !      !if(phase==1) gl%seacalc = .true.
    !      alow2(j) = a_calc(gl, temp, dens22, 0)!/gl%req(j) / Temp
    !      gl%seacalc = .false.
    !
    !
    !      !---------------------------------------------second deriv stencil method
    !      do i=1, gl%ncomp
    !      dchempot_dxi(j,i) = ( -aup2(j) + 16.d0*aup1(j) - 30.d0*a_norm(j) + 16.d0*alow1(j) - alow2(j) ) / (12.d0 * ((molupdiff - mollowdiff)**2.d0))
    !      end do
    !      !chempot(j) = ((a_up(j) - a_low(j)))/(molup-mollow)!(molfractions_orig(j)*2.D0*factorx) !/gl%req(j) / temp
    !
    !      molfractions(j) = molfractions_orig(j)
    !      molfractions(gl%ncomp) = molfractions_orig(gl%ncomp)
    !
    !  end do
    !
    !  !Write the original composition back
    !  gl%molfractions = molfractions_orig
    !  gl%tredmix = tredmix_orig
    !  gl%rhoredmix = rhoredmix_orig
    !
    !  end subroutine dchempot_dxi_num_calc_false
    !  !*******************************************************************************************
    !  !*******************************************************************************************




    end module vle_derivs_module
