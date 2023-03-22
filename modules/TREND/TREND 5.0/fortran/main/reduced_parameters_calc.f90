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

    ! module for file reduced_parameters_calc.f90
    submodule (reduced_parameters_calc_module) impl
    !global use inclusion
    use module_all_types
    use cubic_eos_module


    contains




!New subroutine for PURE SUBSTANCES. This routine updates temperature dependent parameters of other EOS than Helmholtz EOS
!e.g. cubic EOS like the SRK. THe parameter a of the SRK is temperature dependent and thus the parameter has to be recalculated, whenever the 
!temperature changes
!Andreas March 2013
module subroutine update_parameters(gl,Temp, nrsubst)






implicit none

    type(type_gl) :: gl


double precision:: Temp, T
integer :: nrsubst

if (nrsubst < 1) return !Error, wrong nrsubst

Select Case (gl%Eq_type(nrsubst))
    Case(2) !SRK
        T = Temp
        call init_SRK(gl,T, nrsubst)
    Case(3) !PR
        T = Temp
        call init_PR(gl,T, nrsubst)        
End select
    
end subroutine update_parameters
    
    
    
module subroutine reduced_parameters_calc(gl,Temp)
!DEC$ ATTRIBUTES DLLEXPORT :: reduced_parameters_calc
!-------------------------------------------------------------
! Reducing function for density and temperature
! Based on quadratic mixing rules
! Hailong Li, Sep 2009
! modified: Johannes Gernert, July 2010
!-------------------------------------------------------------








implicit none

    type(type_gl) :: gl


double precision :: Temp
double precision:: T

double precision:: xi, xj, rhoci, rhocj, tci, tcj
double precision:: betatij, betavij, gammatij, gammavij
double precision:: vredmix!, wmix
double precision:: vstar1, vstar2, vstar3, vij      !auxiliary variable for Costald mixing rules
integer:: i, j


gl%tredmix = 0.D0
gl%rhoredmix = 0.D0
vredmix = 0.d0

if ((gl%Mix_type == 1) .or. (gl%Mix_type == 11) .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) then  !Lorentz Berthelot or modified mixing rules used Andreas March 2012. 
                                                                                !Mixtype = 11 is for quadratic mixing rules for the reduced residual Helmholtz energy (Andreas, February 2016)
                                                                                !Mixtype = 12 is for excess based departure function mixing rules (Andreas, July 2016)

    do i = 1, gl%ncomp     !v over the critical parameters of the pure fluid
        gl%tredmix = gl%tredmix + gl%molfractions(i)**2*gl%tc(i)
        vredmix = vredmix + gl%molfractions(i)**2/gl%rhored(i)!*wm(i)
    end do

    do i = 1, gl%ncomp - 1          ! sum over the binary reducing functions
      do j = i+1, gl%ncomp
        xi = gl%MOLFRACTIONS(i)
        xj = gl%MOLFRACTIONS(j)
        rhoci = gl%rhored(i)!/wm(i)
        rhocj = gl%rhored(j)!/wm(j)
        tci = gl%tc(i)
        tcj = gl%tc(j)
        betatij = gl%RFBETAT(i, j)
        betavij = gl%RFBETARHO(i,j)
        gammatij = gl%RFGAMMAT(i, j)
        gammavij = gl%RFGAMMARHO(i, j)
        vredmix = vredmix + 2.D0*xi*xj*betavij*gammavij*(xi + xj)/(betavij**2*xi + xj) &
                  & *0.125d0*(1.d0/(rhoci)**(1.d0/3.d0) + 1.d0/(rhocj)**(1.d0/3.d0))**3
        gl%tredmix = gl%tredmix + 2.D0*xi*xj*betatij*gammatij*(xi + xj)/(betatij**2*xi + xj)*(tci*tcj)**0.5d0
      end do
    end do
    !!switched back to molar units! J.Gernert, Aug. 2010
    !call wm_mix_calc(wmix)
    gl%rhoredmix = 1.D0/vredmix !*wmix
    
    !Calculate the covolume of the mixture in case of
    if ((gl%Mix_type == 12) .or. (gl%Mix_type == 13)) then
        gl%b_HelmgE = 0.d0
        do i = 1, gl%ncomp
            !Andreas, Erik, Dezember 2018
            !gl%bi_HelmgE(i) = 0.08664D0*R_HelmgE* gl%tc(i) / (gl%pc(i)*1.D6/gl%factorpress)
            !------
            !!EXCEPTION FOR WATER, TEST, MAYBE DELETE LATER OR ADJUST FOR ALL FLUIDS, Andreas Jäger, January 2018
            !if ((dabs(gl%tc(i) - 647.096D0) < 0.01D0) .and. (dabs(gl%pc(i) - 22.064D0) < 0.01D0)) then
            !    !If the substance is water, adjust the covolume b of water such that u is 1.17 for water (it really is 0.889757593595087)
            !    gl%bi_HelmgE(i) = gl%bi_HelmgE(i) / gl%u_pack * 0.889757593595087D0
            !end if 
            !!EXCEPTION FOR ACETONE, TEST, MAYBE DELETE LATER OR ADJUST FOR ALL FLUIDS, Andreas Jäger, January 2018
            !if ((dabs(gl%tc(i) - 508.1) < 0.01D0) .and. (dabs(gl%pc(i) - 4.69241592505793D0) < 0.01D0)) then
            !    !If the substance is acetone, adjust the covolume b of acetone such that u is 1.17 for acetone (it really is 0.994174158860908)
            !    gl%bi_HelmgE(i) = gl%bi_HelmgE(i) / gl%u_pack * 0.994174158860908D0
            !end if
            !!EXCEPTION FOR CO2, TEST, MAYBE DELETE LATER OR ADJUST FOR ALL FLUIDS, Andreas Jäger, May 2018
            !if ((dabs(gl%tc(i) - 304.1282D0) < 0.01D0) .and. (dabs(gl%pc(i) - 7.3773D0) < 0.01D0)) then
            !    !If the substance is co2, adjust the covolume b of co2 such that u is 1.17 for co2 (it really is 1.25770722213069)
            !    gl%bi_HelmgE(i) = gl%bi_HelmgE(i) / gl%u_pack * 1.25770722213069D0
            !    !write(*,*) "CO2_dens_adj"
            !end if
            !!EXCEPTION FOR Ethanol, TEST, MAYBE DELETE LATER OR ADJUST FOR ALL FLUIDS, Andreas Jäger, May 2018
            !if ((dabs(gl%tc(i) - 514.71D0) < 0.01D0) .and. (dabs(gl%pc(i) - 6.2680D0) < 0.01D0)) then
            !    !If the substance is ethanol, adjust the covolume b of ethanol such that u is 1.17 for ethanol (it really is 1.05753963623426)
            !    gl%bi_HelmgE(i) = gl%bi_HelmgE(i) / gl%u_pack * 1.05753963623426D0
            !    !write(*,*) "Ethanol_dens_adj"
            !end if
            !!EXCEPTION FOR Propene, TEST, MAYBE DELETE LATER OR ADJUST FOR ALL FLUIDS, Andreas Jäger, May 2018
            !if ((dabs(gl%tc(i) - 364.211D0) < 0.01D0) .and. (dabs(gl%pc(i) - 4.555D0) < 0.01D0)) then
            !    !If the substance is propene, adjust the covolume b of propene such that u is 1.17 for propene (it really is 1.19752490734425)
            !    gl%bi_HelmgE(i) = gl%bi_HelmgE(i) / gl%u_pack * 1.19752490734425D0
            !    !write(*,*) "Propene_dens_adj"
            !end if 
            !!EXCEPTION FOR Benzene, TEST, MAYBE DELETE LATER OR ADJUST FOR ALL FLUIDS, Andreas Jäger, May 2018
            !if ((dabs(gl%tc(i) - 562.02D0) < 0.01D0) .and. (dabs(gl%pc(i) - 4.907277D0) < 0.01D0)) then
            !    !If the substance is Benzene, adjust the covolume b of Benzene such that u is 1.17 for Benzene (it really is 1.16328536377193)
            !    gl%bi_HelmgE(i) = gl%bi_HelmgE(i) / gl%u_pack * 1.16328536377193D0
            !    !write(*,*) "Benzene_dens_adj"
            !end if
            !!EXCEPTION FOR Ethane, TEST, MAYBE DELETE LATER OR ADJUST FOR ALL FLUIDS, Andreas Jäger, May 2018
            !if ((dabs(gl%tc(i) - 305.322D0) < 0.01D0) .and. (dabs(gl%pc(i) - 4.8722D0) < 0.01D0)) then
            !    !If the substance is Ethane, adjust the covolume b of Ethane such that u is 1.17 for Ethane (it really is 1.22482003184118)
            !    gl%bi_HelmgE(i) = gl%bi_HelmgE(i) / gl%u_pack * 1.22482003184118D0
            !    !write(*,*) "Ethane_dens_adj"
            !end if 
            !!EXCEPTION FOR SO2, TEST, MAYBE DELETE LATER OR ADJUST FOR ALL FLUIDS, Andreas Jäger, August 2018
            !if ((dabs(gl%tc(i) - 430.64D0) < 0.01D0) .and. (dabs(gl%pc(i) - 7.88657889302166D0) < 0.01D0)) then
            !    !If the substance is SO2, adjust the covolume b of SO2 such that u is 1.17 for SO2 (it really is 1.11308150289714)
            !    gl%bi_HelmgE(i) = gl%bi_HelmgE(i) / gl%u_pack * 1.11308150289714D0
            !    !write(*,*) "SO2_dens_adj"
            !end if 
            !!EXCEPTION FOR H2S, TEST, MAYBE DELETE LATER OR ADJUST FOR ALL FLUIDS, Andreas Jäger, August 2018
            !if ((dabs(gl%tc(i) -373.1D0) < 0.01D0) .and. (dabs(gl%pc(i) - 8.99873466564311D0) < 0.01D0)) then
            !    !If the substance is H2S, adjust the covolume b of H2S such that u is 1.17 for H2S (it really is 1.11308150289714)
            !    gl%bi_HelmgE(i) = gl%bi_HelmgE(i) / gl%u_pack * 1.11308150289714D0
            !    !write(*,*) "H2S_dens_adj"
            !end if 			
            !------
            gl%b_HelmgE = gl%b_HelmgE + gl%molfractions(i) * gl%bi_HelmgE(i)
        end do
               
        !Calculate the density of the mixture with the SRK-mixture covolume and the constant packing fraction
        gl%rho_mix_ref = 1.D0 / gl%u_pack / gl%b_HelmgE
        !Calculate the densities of the pure fluids with the SRK covolumes and the constant packing fraction
        !Andreas, Erik, Dezember 2018
        !do i = 1, gl%ncomp
        !    gl%rho_i_ref(i) = 1.D0 / gl%u_pack / gl%bi_HelmgE(i)
        !end do
    end if

elseif ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then    !SRK with SRK mixing rules used

    gl%rhoredmix = 1.D0
    gl%tredmix = 1.D0
    T = Temp
    
    call init_SRK(gl,T, 0)

elseif ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then    !PR with PR mixing rules used

    gl%rhoredmix = 1.D0
    gl%tredmix = 1.D0
    T = Temp
    call init_PR(gl,T, 0)

elseif (gl%Mix_type == 4) then    !LKP with LKP mixing rules used
    
    !U. Plöcker, H. Knapp, J. Prausnitz, Ind. Eng. Chem. Process Des. Dev., 17(3): 32- (1978)
    
    gl%accenLKPmix= 0.d0
    
    !p. 325, Eq. (11)
    do i=1, gl%ncomp 
        gl%accenLKPmix = gl%accenLKPMix + gl%molfractions(i) * gl%accen(i)
    end do
    
    !p. 325, Eq. (14)    
    if (gl%ncomp > 1) then
        gl%zcLKP=0.2905d0 - 0.085d0 * gl%accenLKPmix
    end if
    
    !Andreas, October 2014
    !It can be dangerous to check if the module variables are non zero at this point. 
    !Since the variables are only initialized once in the new TREND version, this will lead to problems when the composition is changed
    !If ((DABS(tc_ij_LKP(1,1)) < 1.d-12) .OR. (DABS(vc_ij_LKP(1,1)) < 1.d-12)) Then 
    
        !p. 325, Eq. (12); Eq. (13)
        do i = 1, gl%ncomp
            do j = 1, gl%ncomp
                gl%tc_ij_LKP(i,j) = gl%kij_LKP(i, j) * (gl%tc(i) * gl%tc(j)) ** (0.5d0)
                gl%vc_ij_LKP(i,j) = 0.125d0 * (gl%vc_LKP(i) ** (1.d0 / 3.d0) + gl%vc_LKP(j) ** (1.d0 / 3.d0)) ** 3
            end do
        end do
    !End if
    
    !p. 325, Eq. (10); Eq. (9)
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            vredmix = vredmix + gl%molfractions(i) * gl%molfractions(j) * gl%vc_ij_LKP(i,j)
            gl%tredmix = gl%tredmix + gl%molfractions(i) * gl%molfractions(j) * gl%vc_ij_LKP(i,j) ** (0.25d0) * gl%tc_ij_LKP(i,j)
        end do
    end do

    gl%rhoredmix=1.d0/vredmix
    
    gl%tredmix = gl%tredmix * vredmix ** (-0.25d0)

elseif (gl%mix_type == 6) then
    !do i = 1, ncomp     !v over the critical parameters of the pure fluid
    !    tredmix = tredmix + molfractions(i)*tc(i)
    !    vredmix = vredmix + molfractions(i)/rhored(i)!*wm(i)
    !end do
    !rhoredmix = 1.d0 / vredmix
    gl%tredmix = 1.d0
    gl%rhoredmix = 1.d0
    
elseif (gl%mix_type == 9) then     !Costald with Costald mixing rules used
!-------------------------------------------------------------------------------
!Taken from
!G. H. Thomson, K. R. Brobst, R. W. Hankinson, AIChE Journal, 28(4):671-676 (1982)
!-------------------------------------------------------------------------------

    gl%omegaCOS = 0.d0

    do i=1, gl%ncomp
        !Eq. (16) acentric faktor mixture
        gl%omegaCOS = gl%omegaCOS + gl%molfractions(i) * gl%omega(i)
    end do
    
    vstar1 = 0.d0
    vstar2 = 0.d0
    vstar3 = 0.d0
    
    do i=1, gl%ncomp
        vstar1 = vstar1 + gl%molfractions(i) * gl%vstar(i)
        vstar2 = vstar2 + gl%molfractions(i) * gl%vstar(i) ** (2.d0 / 3.d0)
        vstar3 = vstar3 + gl%molfractions(i) * gl%vstar(i) ** (1.d0 / 3.d0)
    end do
    
    !Eq. (15) critical volumen mixture
    gl%vstarCOS = 0.25d0 * (vstar1 + 3.d0 * vstar2 * vstar3)
    
    vij = 0.d0
    
    do i=1, gl%ncomp
        do j=1, gl%ncomp
            !Numerator of Eq. (13) with Eq. (14) 
            vij = vij +  gl%molfractions(i) * gl%molfractions(j) * ((gl%vstar(i) * gl%tc(i) * gl%vstar(j) * gl%tc(j)) ** (0.5d0))
        end do
    end do
    
    !Eq. (13) critical temperature mixture
    gl%TcCOS = vij / gl%vstarCOS       

elseif (gl%mix_type == 19) then     !One-fluid mixing rules for multiparameter EOS, no reducing parameters necessary. Andreas Jäger, July 2016
    
    gl%rhoredmix = 1.D0
    gl%tredmix = 1.D0    
    
end if



end subroutine reduced_parameters_calc
!**************************************************************************

! This set of routines is a adaptation of the routines found in reduced_parameters_calc_old.f95
! The changes are based on the algorithm of K. Hall et al., where the function g(T,p,x) for 
! N components has only N-1 + 2 variables, namely T, p and x_1, x_2, ... x_N-1. The composition 
! of the last component is replaced by: x_N = 1 - (x_1 + x_2 + ...+ x_N-1)
! All derivatives with respect to xi or ni are changed according to the above definition.
!
! J. Gernert, A. Jäger, Nov. 2010

!**************************************************************************
module subroutine ndYr_dni(gl,ndTred_dni, ndrhored_dni)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE FIRST DERIVATIVE OF THE REDUCING 
!FUNCTIONS FOR DENSITY AND TEMPERATURE WITH REGARD TO THE AMOUNT OF SUBSTANCE
! OF THE COMPONENT i: n*(d rhored/d ni)_nj, n*(d tred/d ni)_nj
!**************************************************************************
! Johannes Gernert, Bochum, Nov. 2010
!--------------------------------------------------------------------------
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! ndTred_dni     -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF tred WITH RESPECT TO ALL AMOUNTS OF SUBSTANCE
! ndrhored_dni   -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF rhored WITH RESPECT TO ALL  AMOUNTS OF SUBSTANCE
!--------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl


double precision, dimension(30) ::ndTred_dni, ndrhored_dni   ! Return vectors, maximum 30 components
double precision, dimension(30)::dTred_dxi, drhored_dxi   
integer :: i, k
double precision:: xk

ndTred_dni = 0.D0
ndrhored_dni = 0.D0

!Andreas March 2012
!If SRK with SRK mixing rules is used, all derivatives of the reducing parameters wrt n are 0
if ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then
    ndTred_dni = 0.D0
    ndrhored_dni = 0.D0   
    return
End if

!Andreas November 2015
!If PR with PR mixing rules is used, all derivatives of the reducing parameters wrt n are 0
if ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then
    ndTred_dni = 0.D0
    ndrhored_dni = 0.D0   
    return
End if

!Theresa February 2018
!If PC-SAFT, all derivatives of the reducing parameters wrt n are 0
if (gl%Mix_type == 6) then
    ndTred_dni = 0.D0
    ndrhored_dni = 0.D0   
    return
End if

!Andreas Jäger, July 2016
!One-fluid mixing rules for multiparameter EOS, no reducing parameters necessary. Andreas Jäger, July 2016
if (gl%mix_type == 19) then
    ndTred_dni = 0.D0
    ndrhored_dni = 0.D0   
    return
End if

!Do the following for Mix_type = 1 and Mix_type = 11. Mixtype = 11 is for quadratic mixing rules for the reduced residual Helmholtz energy (Andreas, February 2016)

call dYr_dxi(gl,dTred_dxi, drhored_dxi)

! calculate n*(drhored/dni) = drhored/dxi - SUM(xk*(drhored/dxk))
do i = 1, gl%NCOMP
    if (i < gl%ncomp) then
        ndTred_dni(i) = dTred_dxi(i)
        ndrhored_dni(i) = drhored_dxi(i)
    end if
    do k = 1, gl%NCOMP-1
        xk = gl%MOLFRACTIONS(k)
        ndTred_dni(i) = ndTred_dni(i) - xk*dTred_dxi(k)
        ndrhored_dni(i) = ndrhored_dni(i) - xk*drhored_dxi(k)
    end do
end do

end subroutine ndYr_dni
!**************************************************************************

!**************************************************************************
module subroutine dYr_dxi(gl,dTred_dxi, drhored_dxi)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE FIRST DERIVATIVE OF THE REDUCING 
!FUNCTIONS FOR DENSITY AND TEMPERATURE WITH REGARD TO THE MOLE FRACTION
! OF THE COMPONENT i: d(rhored)/d(xi)_xj, d(tred)/d(xi)_xj
!**************************************************************************
! Johannes Gernert, Nov. 2010
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! dTred_dxi      -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVATIVES
!                   OF tred WITH RESPECT TO ALL MOLAR FRACTIONS
! drhored_dxi    -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF rhored WITH RESPECT TO ALL  MOLAR FRACTIONS
!--------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl


double precision, dimension(30)::dTred_dxi, drhored_dxi
double precision, dimension(30)::dvred_dxi
integer::i, j, k
double precision:: xi, xj, xk, xn, rhoci, rhocj, rhock, rhocn, tci, tcj, tck, tcn
! parameters for the 1st loop
double precision:: betatji, betarhoji, gammatji, gammarhoji
double precision:: ct_ji, crho_ji
double precision:: dftji_dxi, dfrhoji_dxi
! parameters for the 2nd loop
double precision:: betatij, betarhoij, gammatij, gammarhoij
double precision:: ct_ij, crho_ij
double precision:: dftij_dxi, dfrhoij_dxi
! parameters for the combination i, n
double precision:: betatin, betarhoin, gammatin, gammarhoin
double precision:: ct_in, crho_in
double precision:: dftin_dxi, dfrhoin_dxi
! parameters for the 3rd loop
double precision:: betatkn, betarhokn, gammatkn, gammarhokn
double precision:: ct_kn, crho_kn
double precision:: dftkn_dxi, dfrhokn_dxi
! LKP variables
double precision:: vredmix
double precision, dimension(30) :: help_da, help_db
double precision :: help_a, help_b

double precision :: b_1_3

do i = 1, gl%ncomp
    dTred_dxi(i) = 0.D0
    drhored_dxi(i) = 0.D0
    dvred_dxi(i) = 0.D0
end do

b_1_3 = 1.d0/3.d0

!Andreas March 2012
If ((gl%Mix_type == 1) .or. (gl%Mix_type == 11) .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) Then     !Mixtype = 11 is for quadratic mixing rules for the reduced residual Helmholtz energy (Andreas, February 2016),  Mixtype = 12 is for excess based departure function mixing rules (Andreas, July 2016)
    xn = gl%molfractions(gl%ncomp)
    rhocn = gl%rhored(gl%ncomp)
    tcn = gl%tc(gl%ncomp)

    do i = 1, gl%NCOMP - 1
        rhoci = gl%rhored(i)
        tci = gl%tc(i)
        xi = gl%MOLFRACTIONS(i)
        dTred_dxi(i) = 2.D0*(xi*tci - xn*tcn)
        dvred_dxi(i) = 2.D0*(xi/rhoci - xn/rhocn)
        ! 1st loop for j < i
        do j = 1, i-1
            xj = gl%MOLFRACTIONS(j)
            rhocj = gl%rhored(j)
            tcj = gl%tc(j)
            betatji = gl%RFBETAT(j, i)
            betarhoji = gl%RFBETARHO(j, i)
            gammatji = gl%RFGAMMAT(j, i)
            gammarhoji = gl%RFGAMMARHO(j, i)
            ct_ji = 2.D0*betatji*gammatji*(tci*tcj)**0.5D0      !critical term, temperature
            crho_ji = 2.D0*betarhoji*gammarhoji*0.125d0*(1.D0/(rhoci)**b_1_3 + 1.D0/(rhocj)**b_1_3)** 3  !critical term, density
            !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
            dftji_dxi = xj*(xj + xi)/(betatji**2*xj + xi) + &
                        & xj*xi/(betatji**2*xj + xi)*(1.D0 - (xj + xi)/(betatji**2*xj + xi))
            !derivative of the composition-dependent part of the reducing function for the density w.r.t. xi:  
            dfrhoji_dxi = xj*(xj + xi)/(betarhoji**2*xj + xi) + &
                        & xj*xi/(betarhoji**2*xj + xi)*(1.D0 - (xj + xi)/(betarhoji**2*xj + xi))
            ! 1st derivative of Tred w.r.t xi - 1st summation for j = 1 to i-1
            dTred_dxi(i) = dTred_dxi(i) +  ct_ji * dftji_dxi
            ! 1st derivative of rhored w.r.t xi - 1st summation for j = 1 to i-1
            dvred_dxi(i) = dvred_dxi(i) +  crho_ji * dfrhoji_dxi
        end do
        !2nd loop for j > i
        do j = i + 1, gl%NCOMP-1
            xj = gl%MOLFRACTIONS(j)
            rhocj = gl%rhored(j)
            tcj = gl%tc(j)        
            betatij = gl%RFBETAT(i, j)
            betarhoij = gl%RFBETARHO(i, j)
            gammatij = gl%RFGAMMAT(i, j)
            gammarhoij = gl%RFGAMMARHO(i, j)
            ct_ij = 2.D0*betatij*gammatij*(tci*tcj)**0.5D0      !critical term, temperature
            crho_ij = 2.D0*betarhoij*gammarhoij*0.125d0*(1.D0/(rhoci)**b_1_3 + 1.D0/(rhocj)**b_1_3)** 3  !critical term, density
            !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
            dftij_dxi = xj*(xi + xj)/(betatij**2*xi + xj) + &
                     & xi*xj/(betatij**2*xi + xj)*(1d0 - betatij**2*(xi + xj)/(betatij**2*xi + xj))
            !derivative of the composition-dependent part of the reducing function for the density w.r.t. xi: 
            dfrhoij_dxi = xj*(xi + xj)/(betarhoij**2*xi + xj) + &
                       & xi*xj/(betarhoij**2*xi + xj)*(1d0 - betarhoij**2*(xi + xj)/(betarhoij**2*xi + xj))
            ! 1st derivative of Tred w.r.t xi - 2nd summation for j = i+1 to ncomp
            dTred_dxi(i) = dTred_dxi(i) +  ct_ij * dftij_dxi
            ! 1st derivative of rhored w.r.t xi - 2nd summation for j = i+1 to ncomp
            dvred_dxi(i) = dvred_dxi(i) +  crho_ij * dfrhoij_dxi
        end do
        ! Now the additional part with the  combinations of component i with component N       
        betatin = gl%RFBETAT(i, gl%ncomp)
        betarhoin = gl%RFBETARHO(i, gl%ncomp)
        gammatin = gl%RFGAMMAT(i, gl%ncomp)
        gammarhoin = gl%RFGAMMARHO(i, gl%ncomp)
        ct_in = 2.D0*betatin*gammatin*(tci*tcn)**0.5D0      !critical term, temperature
        crho_in = 2.D0*betarhoin*gammarhoin*0.125d0*(1.D0/(rhoci)**b_1_3 + 1.D0/(rhocn)**b_1_3)** 3  !critical term, density
        dftin_dxi = xn*(xi + xn)/(betatin**2*xi + xn) + (1.D0 - betatin**2)*xi*xn*xn/(betatin**2*xi + xn)**2
        dfrhoin_dxi = xn*(xi + xn)/(betarhoin**2*xi + xn) + (1.D0 - betarhoin**2)*xi*xn*xn/(betarhoin**2*xi + xn)**2
        dTred_dxi(i) = dTred_dxi(i) + ct_in*dftin_dxi
        dvred_dxi(i) = dvred_dxi(i) +  crho_in * dfrhoin_dxi
        ! the following sumation takes into account that xn = 1 - (x1 + x2 + ... + xn-1) creates derivatives with respect to xi
        ! even if the outer sum is now at another component xk
        do k = 1, gl%ncomp - 1
            xk = gl%MOLFRACTIONS(k)
            rhock = gl%rhored(k)
            tck = gl%tc(k)
            betatkn = gl%RFBETAT(k, gl%ncomp)
            betarhokn = gl%RFBETARHO(k, gl%ncomp)
            gammatkn = gl%RFGAMMAT(k, gl%ncomp)
            gammarhokn = gl%RFGAMMARHO(k, gl%ncomp)
            ct_kn = 2.D0*betatkn*gammatkn*(tck*tcn)**0.5D0      !critical term, temperature
            crho_kn = 2.D0*betarhokn*gammarhokn*0.125d0*(1.D0/(rhock)**b_1_3 + 1.D0/(rhocn)**b_1_3)** 3  !critical term, density    
            dftkn_dxi = - xk*(xk + xn)/(betatkn**2*xk + xn) + (1.D0 - betatkn**2)*xk*xk*xn/(betatkn**2*xk + xn)**2
            dfrhokn_dxi = - xk*(xk + xn)/(betarhokn**2*xk + xn) + (1.D0 - betarhokn**2)*xk*xk*xn/(betarhokn**2*xk + xn)**2
            dTred_dxi(i) = dTred_dxi(i) + ct_kn*dftkn_dxi
            dvred_dxi(i) = dvred_dxi(i) +  crho_kn * dfrhokn_dxi
        end do    
        drhored_dxi(i) = - gl%rhoredmix**2*dvred_dxi(i) ! drhored/dxi = -rhored^2*(d(1/rhored)dxi) - the routine above actually calculates dvred/dxi ...
    end do

!If SRK with SRK mixing rules is used or PR with PR mixing rules, all derivatives of the reducing parameters wrt x are 0
elseif ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22) .OR. (gl%Mix_type == 3) .OR. (gl%Mix_type == 31)) then
    dTred_dxi = 0.D0
    drhored_dxi = 0.D0  
    return

elseif(gl%Mix_type == 4) then
    
    help_a=0.d0
    help_b=0.d0
    help_Da=0.d0
    help_Db=0.d0
    
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            help_b = help_b + gl%molfractions(i) * gl%molfractions(j) * (gl%vc_ij_LKP(i, j) ** 0.25d0) * gl%tc_ij_LKP(i, j)
        end do
    end do
    
    do i = 1, gl%ncomp - 1
        do j = 1, gl%ncomp - 1
            dvred_dxi(i) = dvred_dxi(i) + gl%molfractions(j) * (gl%vc_ij_LKP(j, i) + gl%vc_ij_LKP(i, j) - gl%vc_ij_LKP(j, gl%ncomp) - gl%vc_ij_LKP(gl%ncomp, j))
        end do
        dvred_dxi(i) = dvred_dxi(i) - 2.d0 * gl%molfractions(gl%ncomp) * gl%vc_ij_LKP(gl%ncomp, gl%ncomp) + gl%molfractions(gl%ncomp) * (gl%vc_ij_LKP(gl%ncomp, i) + gl%vc_ij_LKP(i, gl%ncomp))
    end do 
    
    drhored_dxi = - gl%rhoredmix**2*dvred_dxi ! drhored/dxi = -rhored^2*(d(1/rhored)dxi) - the routine above actually calculates dvred/dxi ...
   
    vredmix = 1.d0/gl%rhoredmix
    help_a = vredmix ** (-0.25d0)

    do i = 1, gl%ncomp - 1
        help_Da(i) = -0.25d0 * (vredmix **(-1.25d0)) * dvred_dxi(i)
    end do
    do i = 1, gl%ncomp - 1
        do j = 1, gl%ncomp - 1
            help_Db(i) = help_Db(i) + gl%molfractions(j) * ((gl%vc_ij_LKP(j, i) ** 0.25d0) * gl%tc_ij_LKP(j, i) + (gl%vc_ij_LKP(i, j) ** 0.25d0) * gl%tc_ij_LKP(i, j) &
                     & - (gl%vc_ij_LKP(j, gl%ncomp) ** 0.25d0) * gl%tc_ij_LKP(j, gl%ncomp) - (gl%vc_ij_LKP(gl%ncomp, j) ** 0.25d0) * gl%tc_ij_LKP(gl%ncomp, j))
        end do
        help_Db(i) = help_Db(i) - 2.d0 * gl%molfractions(gl%ncomp) * gl%vc_ij_LKP(gl%ncomp, gl%ncomp) ** 0.25d0 * gl%tc_ij_LKP(gl%ncomp, gl%ncomp) &
                 & + gl%molfractions(gl%ncomp) * ((gl%vc_ij_LKP(gl%ncomp, i) ** 0.25d0) * gl%tc_ij_LKP(gl%ncomp, i) + (gl%vc_ij_LKP(i, gl%ncomp) ** 0.25d0) * gl%tc_ij_LKP(i, gl%ncomp))
        dTred_dxi(i) = help_a * help_Db(i) + help_Da(i) * help_b
    end do

    
elseif (gl%mix_type ==6) then     !PC-SAFT, no reducing parameters necessary. Theresa, Frebruary 2018
    dTred_dxi = 0.D0
    drhored_dxi = 0.D0  
    return    

    
    elseif (gl%mix_type == 19) then     !One-fluid mixing rules for multiparameter EOS, no reducing parameters necessary. Andreas Jäger, July 2016
    dTred_dxi = 0.D0
    drhored_dxi = 0.D0  
    return    
    
End If

end subroutine dYr_dxi
!**************************************************************************
!
!**************************************************************************
module subroutine d2Yr_dxi2(gl,d2Tred_dxi2, d2rhored_dxi2)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE REDUCING 
!FUNCTIONS FOR DENSITY AND TEMPERATURE WITH REGARD TO THE COMPOSITION OF THE 
!COMPONENT i (d^2 rhored/d xi^2)_xj, (d^2 tred/d xi^2)_xj
!**************************************************************************
! Johannes Gernert, Nov. 2010
!--------------------------------------------------------------------------
!
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! d2Tred_dxi2     -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF tred WITH RESPECT TO ALL COMPOSITIONS
! d2rhored_dxi2   -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF rhored WITH RESPECT TO ALL COMPOSITIONS
!--------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl



double precision, dimension(30)::d2Tred_dxi2, d2rhored_dxi2   ! Return vectors, maximum 30 components
double precision, dimension(30)::dTred_dxi, drhored_dxi         ! Return vectors, 1st derivative, maximum 30 components
double precision, dimension(30)::d2vred_dxi2                    ! internal vector for vred
integer::i, j, k
double precision:: xi, xj, xk, xn, rhoci, rhocj, rhock, rhocn, tci, tcj, tck, tcn
! parameters for the 1st loop
double precision:: betatji, betarhoji, gammatji, gammarhoji
double precision:: ct_ji, crho_ji
double precision:: d2ftji_dxi2, d2frhoji_dxi2
! parameters for the 2nd loop
double precision:: betatij, betarhoij, gammatij, gammarhoij
double precision:: ct_ij, crho_ij
double precision:: d2ftij_dxi2, d2frhoij_dxi2
! parameters for the i-n part
double precision:: betatin, betarhoin, gammatin, gammarhoin
double precision:: ct_in, crho_in
double precision:: d2ftin_dxi2, d2frhoin_dxi2
! parameters for the 3rd loop
double precision:: betatkn, betarhokn, gammatkn, gammarhokn
double precision:: ct_kn, crho_kn
double precision:: d2ftkn_dxi2, d2frhokn_dxi2

double precision :: b_1_3

d2Tred_dxi2 = 0.D0
d2vred_dxi2 = 0.D0
d2rhored_dxi2 = 0.D0
dTred_dxi = 0.D0
drhored_dxi = 0.D0

b_1_3 = 1.d0/3.d0


If ((gl%Mix_type == 1) .or. (gl%Mix_type == 11) .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) then !Mixtype = 11 is for quadratic mixing rules for the reduced residual Helmholtz energy (Andreas, February 2016), Mixtype = 12 is for excess based departure function mixing rules (Andreas, July 2016)

    ! the first derivative of rhored w.r.t. xi is needed to calculate d^2rhored/dxi^2
    call dYr_dxi(gl,dTred_dxi, drhored_dxi)

    rhocn = gl%rhored(gl%ncomp)
    tcn = gl%tc(gl%ncomp)
    xn = gl%MOLFRACTIONS(gl%ncomp)

    do i = 1, gl%NCOMP-1
        rhoci = gl%rhored(i)
        tci = gl%tc(i)
        xi = gl%MOLFRACTIONS(i)
        d2Tred_dxi2(i) = 2.D0*(tci + tcn)
        d2vred_dxi2(i) = 2.D0*(1.D0/rhoci + 1.D0/rhocn)
        ! 1st loop for j < i
        do j = 1, i-1
            xj = gl%MOLFRACTIONS(j)
            rhocj = gl%rhored(j)
            tcj = gl%tc(j)
            betatji = gl%RFBETAT(j, i)
            betarhoji = gl%RFBETARHO(j, i)
            gammatji = gl%RFGAMMAT(j, i)
            gammarhoji = gl%RFGAMMARHO(j, i)
            ct_ji = 2.D0*betatji*gammatji*(tci*tcj)**0.5D0      !critical term, temperature
            crho_ji = 2.D0*betarhoji*gammarhoji*0.125d0*(1.D0/(rhoci)**b_1_3 + 1.D0/(rhocj)**b_1_3)** 3  !critical term, density
            !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
            d2ftji_dxi2 = (1.D0 - (xj + xi)/(betatji**2*xj + xi)) / (betatji**2*xj + xi) * &
                      & (2.D0*xj - xj*xi*2.D0/(betatji**2*xj + xi)) 
            !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi:  
            d2frhoji_dxi2 = (1.D0 - (xj + xi)/(betarhoji**2*xj + xi)) / (betarhoji**2*xj + xi) * &
                      & (2.D0*xj - xj*xi*2.D0/(betarhoji**2*xj + xi)) 
            ! 2nd derivative of Tred w.r.t xi - 1st summation for j = 1 to i-1
            d2Tred_dxi2(i) = d2Tred_dxi2(i) +  ct_ji * d2ftji_dxi2
            ! 2nd derivative of vred w.r.t xi - 1st summation for j = 1 to i-1
            d2vred_dxi2(i) = d2vred_dxi2(i) +  crho_ji * d2frhoji_dxi2
        end do
        !2nd loop for j > i
        do j = i + 1, gl%NCOMP-1
            xj = gl%MOLFRACTIONS(j)
            rhocj = gl%rhored(j)
            tcj = gl%tc(j)        
            betatij = gl%RFBETAT(i, j)
            betarhoij = gl%RFBETARHO(i, j)
            gammatij = gl%RFGAMMAT(i, j)
            gammarhoij = gl%RFGAMMARHO(i, j)
            ct_ij = 2.D0*betatij*gammatij*(tci*tcj)**0.5D0      !critical term, temperature
            crho_ij = 2.D0*betarhoij*gammarhoij*0.125d0*(1.D0/(rhoci)**b_1_3 + 1.D0/(rhocj)**b_1_3)** 3  !critical term, density
            !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
            d2ftij_dxi2 = (1.D0 - betatij**2*(xi + xj)/(betatij**2*xi + xj)) / (betatij**2*xi + xj) * &
                       & (2.D0*xj - xi*xj*2.D0*betatij**2/(betatij**2*xi + xj))
            !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi: 
            d2frhoij_dxi2 = (1.D0 - betarhoij**2*(xi + xj)/(betarhoij**2*xi + xj)) / (betarhoij**2*xi + xj) * &
                       & (2.D0*xj - xi*xj*2.D0*betarhoij**2/(betarhoij**2*xi + xj))
            ! 2nd derivative of Tred w.r.t xi - 2nd summation for j = i+1 to ncomp
            d2Tred_dxi2(i) = d2Tred_dxi2(i) +  ct_ij * d2ftij_dxi2
            ! 2nd derivative of vred w.r.t xi - 2nd summation for j = i+1 to ncomp
            d2vred_dxi2(i) = d2vred_dxi2(i) +  crho_ij * d2frhoij_dxi2
        end do
        ! now the binary mixture between component i and the last component
        betatin = gl%RFBETAT(i, gl%ncomp)
        betarhoin = gl%RFBETARHO(i, gl%ncomp)
        gammatin = gl%RFGAMMAT(i, gl%ncomp)
        gammarhoin = gl%RFGAMMARHO(i, gl%ncomp)   
        ct_in = 2.D0*betatin*gammatin*(tci*tcn)**0.5D0      !critical term, temperature
        crho_in = 2.D0*betarhoin*gammarhoin*0.125d0*(1.D0/(rhoci)**b_1_3 + 1.D0/(rhocn)**b_1_3)** 3  !critical term, density    
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
        d2ftin_dxi2 = 2.D0/(betatin**2*xi + xn)*(-(xi + xn) &
                    & + (1.D0 - betatin**2)/(betatin**2*xi + xn)*(xn**2 & 
                    & + ((1.D0 - betatin**2)*xi*xn**2 - betatin**2*xi**2*xn)/(betatin**2*xi + xn)))   
        !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi:
        d2frhoin_dxi2 = 2.D0/(betarhoin**2*xi + xn)*(-(xi + xn) &
                    & + (1.D0 - betarhoin**2)/(betarhoin**2*xi + xn)*(xn**2 & 
                    & + ((1.D0 - betarhoin**2)*xi*xn**2 - betarhoin**2*xi**2*xn)/(betarhoin**2*xi + xn)))     
        ! 2nd derivative of Tred w.r.t xi - 2nd summation for j = i+1 to ncomp
        d2Tred_dxi2(i) = d2Tred_dxi2(i) +  ct_in * d2ftin_dxi2
        ! 2nd derivative of vred w.r.t xi - 2nd summation for j = i+1 to ncomp
        d2vred_dxi2(i) = d2vred_dxi2(i) +  crho_in * d2frhoin_dxi2    
        ! the following summation takes into account that xn = 1 - (x1 + x2 + ... + xn-1) creates derivatives with respect to xi
        ! even if the outer sum is now at another component xk    
        do k = 1, gl%ncomp - 1
            xk = gl%MOLFRACTIONS(k)
            rhock = gl%rhored(k)
            tck = gl%tc(k) 
            betatkn = gl%RFBETAT(k, gl%ncomp)
            betarhokn = gl%RFBETARHO(k, gl%ncomp)
            gammatkn = gl%RFGAMMAT(k, gl%ncomp)
            gammarhokn = gl%RFGAMMARHO(k, gl%ncomp)
            ct_kn = 2.D0*betatkn*gammatkn*(tck*tcn)**0.5D0      !critical term, temperature
            crho_kn = 2.D0*betarhokn*gammarhokn*0.125d0*(1.D0/(rhock)**b_1_3 + 1.D0/(rhocn)**b_1_3)** 3  !critical term, density 
            !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
            d2ftkn_dxi2 = 2.D0*xk**2*(1.D0 - betatkn**2)/((betatkn**2*xk + xn)**2)*(xn/(betatkn**2*xk + xn) - 1.D0)    
            !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi:
            d2frhokn_dxi2 = 2.D0*xk**2*(1.D0 - betarhokn**2)/((betarhokn**2*xk + xn)**2)*(xn/(betarhokn**2*xk + xn) - 1.D0)
            ! 2nd derivative of Tred w.r.t xi - 2nd summation for j = i+1 to ncomp
            d2Tred_dxi2(i) = d2Tred_dxi2(i) +  ct_kn * d2ftkn_dxi2
            ! 2nd derivative of vred w.r.t xi - 2nd summation for j = i+1 to ncomp
            d2vred_dxi2(i) = d2vred_dxi2(i) +  crho_kn * d2frhokn_dxi2              
        end do
        ! the routine above actually calculates d^2vred/dxi^2, therefore must be transformed to:
        ! d^2rhored/dxi^2 = 2/rhored*(drhored/dxi)^2 - rhored^2*(d^2rhored/dxi^2)
        d2rhored_dxi2(i) = 2.D0/gl%rhoredmix*(drhored_dxi(i))**2 - gl%rhoredmix**2*d2vred_dxi2(i) 
    end do

!Andreas March 2012
!If SRK with SRK mixing rules or PR with PR mixing rules is used, all derivatives of the reducing parameters wrt x are 0
Else if ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22) .OR. (gl%Mix_type == 3) .OR. (gl%Mix_type == 31)) then
    d2Tred_dxi2 = 0.D0
    d2rhored_dxi2 = 0.D0
    return

!Stefan & Andreas July 2014
!This is not correct and has to be updated - Preliminary fast and easy solution!!!
!Correct derivatives d2Tred_dxi2 and d2rhored_dxi2 are included in d2Yr_dxidxj!!!
Else if (gl%Mix_type == 4) then
    d2Tred_dxi2 = 0.D0
    d2rhored_dxi2 = 0.D0
    return 
    
elseif (gl%mix_type == 6) then     !PC-SAFT, no reducing parameters necessary. Theresa, Feb 2018
    d2Tred_dxi2 = 0.D0
    d2rhored_dxi2 = 0.D0
    return 

elseif (gl%mix_type == 19) then     !One-fluid mixing rules for multiparameter EOS, no reducing parameters necessary. Andreas Jäger, July 2016
        d2Tred_dxi2 = 0.D0
    d2rhored_dxi2 = 0.D0
    return 

end if
    
end subroutine d2Yr_dxi2
!
!**************************************************************************

!**************************************************************************
module subroutine d2Yr_dxidxj(gl,d2Tred_dxidxj, d2rhored_dxidxj)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE REDUCING 
!FUNCTIONS FOR DENSITY AND TEMPERATURE WITH REGARD TO THE COMPOSITION OF THE 
!COMPONENTS i AND j (d^2(rhored)/d(xi)d(xj)), (d^2(tred)/d(xi)d(xj))
!**************************************************************************
! Johannes Gernert, Nov. 2010
!--------------------------------------------------------------------------
!
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! d2Tred_dxidxj  -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF tred WITH RESPECT TO ALL COMPOSITIONS
! d2rhored_dxdxj -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF rhored WITH RESPECT TO ALL COMPOSITIONS
!--------------------------------------------------------------------------







implicit none

    type(type_gl) :: gl


double precision, dimension(30,30) :: d2Tred_dxidxj, d2rhored_dxidxj
double precision, dimension(30,30):: d2vred_dxidxj
double precision, dimension(30)::dTred_dxi, drhored_dxi         ! 1st derivative, maximum 30 components

double precision::d2fvij_dxidxj, d2fTij_dxidxj, cv_ij, cT_ij
double precision::betatij, betarhoij, gammatij, gammarhoij

double precision::d2fvin_dxidxj, d2fTin_dxidxj, cv_in, cT_in
double precision::betatin, betarhoin, gammatin, gammarhoin

double precision::d2fvjn_dxidxj, d2fTjn_dxidxj, cv_jn, cT_jn
double precision::betatjn, betarhojn, gammatjn, gammarhojn

double precision::d2fvkn_dxidxj, d2fTkn_dxidxj, cv_kn, cT_kn
double precision::betatkn, betarhokn, gammatkn, gammarhokn

double precision:: xi, xj, xk, xn, rhoci, rhocj, rhock, rhocn, tci, tcj, tck, tcn, help
integer::i, j, k

! LKP variables
double precision:: vredmix
double precision, dimension(30) :: help_Da, help_Db
double precision, dimension(30,30) :: help_D2a, help_D2b
double precision :: help_a, help_b

double precision :: b_1_3


d2Tred_dxidxj = 0.D0
d2vred_dxidxj = 0.D0
d2rhored_dxidxj = 0.D0

b_1_3 = 1.d0/3.d0


!Andreas March 2012
If ((gl%Mix_type == 1) .or. (gl%Mix_type == 11) .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) Then !Mixtype = 11 is for quadratic mixing rules for the reduced residual Helmholtz energy (Andreas, February 2016), Mixtype = 12 is for excess based departure function mixing rules (Andreas, July 2016)
    
    ! the first derivative of rhored w.r.t. xi is needed to calculate d^2(rhored)/d(xi)d(xj)
    call dYr_dxi(gl,dTred_dxi, drhored_dxi)  

    xn = gl%MOLFRACTIONS(gl%ncomp)
    tcn = gl%tc(gl%ncomp)
    rhocn = gl%rhored(gl%ncomp)

    do i = 1, (gl%ncomp - 2)
        do j = (i + 1), gl%ncomp - 1
            betatij = gl%RFBETAT(i,j)
            betarhoij = gl%RFBETARHO(i,j)
            gammatij = gl%RFGAMMAT(i,j)
            gammarhoij = gl%RFGAMMARHO(i,j)
            xi = gl%MOLFRACTIONS(i)
            xj = gl%MOLFRACTIONS(j)
            tci = gl%tc(i)
            tcj = gl%tc(j)
            rhoci = gl%rhored(i)
            rhocj = gl%rhored(j)
            !critical term, temperature
            cT_ij = 2.D0*betatij*gammatij*(tci*tcj)**0.5D0
            !critical term, volume
            cv_ij = 2.D0*betarhoij*gammarhoij*0.125d0*(1.D0/(rhoci**b_1_3) + 1.D0/(rhocj**b_1_3))**3
            !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi and xj:
            d2fTij_dxidxj = (xi + xj)/(betatij**2*xi + xj) + xj/(betatij**2*xi + xj)*(1.D0 - (xi + xj)/(betatij**2*xi + xj)) &
                        & + xi/(betatij**2*xi + xj)*(1.D0 - betatij**2*(xi + xj)/(betatij**2*xi + xj)) &
                        & - xi*xj/((betatij**2*xi + xj)**2)*(1.D0 + betatij**2 - 2.D0*betatij**2*(xi + xj)/(betatij**2*xi + xj))
            !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi and xj:
            help = 1.d0 / (betarhoij**2*xi + xj)
            d2fvij_dxidxj = (xi + xj)*help &
                        & + xj*help*(1.D0 - (xi + xj)*help) &
                        & + xi*help*(1.D0 - betarhoij**2*(xi + xj)*help) &
                        & - xi*xj/((betarhoij**2*xi + xj)**2)*(1.D0 + betarhoij**2 &
                        & - 2.D0*betarhoij**2*(xi + xj)*help)
            !2nd derivative of Tred w.r.t xi and xj
            d2Tred_dxidxj(i,j) = 2.D0* tcn + cT_ij*d2fTij_dxidxj
            !2nd derivative of vred w.r.t xi and xj
            d2vred_dxidxj(i,j) = 2.D0/rhocn + cv_ij*d2fvij_dxidxj
            !--------------------------------------------------------------------------------
            ! now the part with the binary combination of component i with the last component
            !--------------------------------------------------------------------------------
            betatin = gl%RFBETAT(i,gl%ncomp)
            betarhoin = gl%RFBETARHO(i,gl%ncomp)
            gammatin = gl%RFGAMMAT(i,gl%ncomp)
            gammarhoin = gl%RFGAMMARHO(i,gl%ncomp)        
            !critical term, temperature
            cT_in = 2.D0*betatin*gammatin*(tci*tcn)**0.5D0
            !critical term, volume
            cv_in = 2.D0*betarhoin*gammarhoin*0.125d0*(1.D0/(rhoci**b_1_3) + 1.D0/(rhocn**b_1_3))**3
            !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi and xj:
            d2fTin_dxidxj = (1.D0 - betatin**2)*2.D0*xi*xn**2/((betatin**2*xi + xn)**3) &
                        & - (1.D0 - betatin**2)*xi*xn/((betatin**2*xi + xn)**2) &
                        & - (xi + xn)/(betatin**2*xi + xn)
            !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi and xj:
            d2fvin_dxidxj = (1.D0 - betarhoin**2)*2.D0*xi*xn**2/((betarhoin**2*xi + xn)**3) &
                        & - (1.D0 - betarhoin**2)*xi*xn/((betarhoin**2*xi + xn)**2) &
                        & - (xi + xn)/(betarhoin**2*xi + xn)
            !2nd derivative of Tred w.r.t xi and xj
            d2Tred_dxidxj(i,j) = d2Tred_dxidxj(i,j) + cT_in*d2fTin_dxidxj      
            !2nd derivative of vred w.r.t xi and xj
            d2vred_dxidxj(i,j) = d2vred_dxidxj(i,j) + cv_in*d2fvin_dxidxj
            !--------------------------------------------------------------------------------
            ! now the part with the binary combination of component j with the last component
            !--------------------------------------------------------------------------------
            betatjn = gl%RFBETAT(j,gl%ncomp)
            betarhojn = gl%RFBETARHO(j,gl%ncomp)
            gammatjn = gl%RFGAMMAT(j,gl%ncomp)
            gammarhojn = gl%RFGAMMARHO(j,gl%ncomp)        
            !critical term, temperature
            cT_jn = 2.D0*betatjn*gammatjn*(tcj*tcn)**0.5D0
            !critical term, volume
            cv_jn = 2.D0*betarhojn*gammarhojn*0.125d0*(1.D0/(rhocj**b_1_3) + 1.D0/(rhocn**b_1_3))**3        
            !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi and xj:
            d2fTjn_dxidxj = - (1.D0 - betatjn**2)*2.D0*xj**2*xn*betatjn**2/((betatjn**2*xj + xn)**3) &
                          & + (1.D0 - betatjn**2)*xj*xn/((betatjn**2*xj + xn)**2) &
                          & - (xj + xn)/(betatjn**2*xj + xn)
            !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi and xj:
            d2fvjn_dxidxj = - (1.D0 - betarhojn**2)*2.D0*xj**2*xn*betarhojn**2/((betarhojn**2*xj + xn)**3) &
                          & + (1.D0 - betarhojn**2)*xj*xn/((betarhojn**2*xj + xn)**2) &
                          & - (xj + xn)/(betarhojn**2*xj + xn)
            !2nd derivative of Tred w.r.t xi and xj
            d2Tred_dxidxj(i,j) = d2Tred_dxidxj(i,j) + cT_jn*d2fTjn_dxidxj        
            !2nd derivative of vred w.r.t xi and xj
            d2vred_dxidxj(i,j) = d2vred_dxidxj(i,j) + cv_jn*d2fvjn_dxidxj                                                  
            ! calculate d^2(rhored)/d(xi)d(xj) from d^2(1/rhored)/d(xi)d(xj)
            do k = 1, gl%ncomp - 1
                betatkn = gl%RFBETAT(k,gl%ncomp)
                betarhokn = gl%RFBETARHO(k,gl%ncomp)
                gammatkn = gl%RFGAMMAT(k,gl%ncomp)
                gammarhokn = gl%RFGAMMARHO(k,gl%ncomp)
                xk = gl%MOLFRACTIONS(k)
                tck = gl%tc(k)
                rhock = gl%rhored(k)
                ct_kn = 2.D0*betatkn*gammatkn*(tck*tcn)**0.5D0      !critical term, temperature
                cv_kn = 2.D0*betarhokn*gammarhokn*0.125D0*(1.D0/(rhock)**b_1_3 + 1.D0/(rhocn)**b_1_3)** 3  !critical term, density 
                !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
                d2fTkn_dxidxj = 2.D0*xk**2*(1.D0 - betatkn**2)/((betatkn**2*xk + xn)**2)*(xn/(betatkn**2*xk + xn) - 1.D0)    
                !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi:
                d2fvkn_dxidxj = 2.D0*xk**2*(1.D0 - betarhokn**2)/((betarhokn**2*xk + xn)**2)*(xn/(betarhokn**2*xk + xn) - 1.D0)
                !2nd derivative of Tred w.r.t xi and xj
                d2Tred_dxidxj(i,j) = d2Tred_dxidxj(i,j) + cT_kn*d2fTkn_dxidxj
                !2nd derivative of vred w.r.t xi and xj
                d2vred_dxidxj(i,j) = d2vred_dxidxj(i,j) + cv_kn*d2fvkn_dxidxj     
            end do
            d2rhored_dxidxj(i,j) = 2.D0/gl%rhoredmix*drhored_dxi(i)*drhored_dxi(j) - gl%rhoredmix**2*d2vred_dxidxj(i,j)
            d2Tred_dxidxj(j,i) = d2Tred_dxidxj(i,j)
            d2rhored_dxidxj(j,i) = d2rhored_dxidxj(i,j)
        end do
    end do
    
!If SRK with SRK mixing rules or PR with PR mixing rules is used, all derivatives of the reducing parameters wrt x are 0
Else if ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22) .OR. (gl%Mix_type == 3) .OR. (gl%Mix_type == 31)) then
    d2Tred_dxidxj = 0.D0
    d2rhored_dxidxj = 0.D0
    return
    
!LKP goes up, baby!
Else if (gl%Mix_type == 4) then
    
    help_a = 0.d0
    help_b = 0.d0
    help_Da = 0.d0
    help_Db = 0.d0
    help_D2a = 0.d0
    help_D2b = 0.d0
    
    ! the first derivative of rhored w.r.t. xi is needed to calculate d^2(rhored)/d(xi)d(xj)
    call dYr_dxi(gl,dTred_dxi, drhored_dxi)
        
    vredmix=1.d0/gl%rhoredmix
    
    do i = 1, gl%ncomp - 1
        do j = 1, gl%ncomp - 1
            d2vred_dxidxj(i, j) = gl%vc_ij_LKP(j, i) + gl%vc_ij_LKP(i, j) - gl%vc_ij_LKP(j, gl%ncomp) - gl%vc_ij_LKP(gl%ncomp, j) + 2.d0 * gl%vc_ij_LKP(gl%ncomp, gl%ncomp) &
                               & - gl%vc_ij_LKP(gl%ncomp, i) - gl%vc_ij_LKP(i, gl%ncomp)
            d2rhored_dxidxj(i,j) = 2.D0*vredmix*drhored_dxi(i)*drhored_dxi(j) - (gl%rhoredmix**2)*d2vred_dxidxj(i,j)
        end do
    end do
      
    
    help_a = vredmix **(-0.25d0)
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            help_b = help_b + gl%molfractions(i) * gl%molfractions(j) * (gl%vc_ij_LKP(i, j) ** 0.25d0) * gl%tc_ij_LKP(i, j)
        end do
    end do
    do i = 1, gl%ncomp - 1
        help_Da(i) = -0.25d0 * (vredmix ** (-1.25d0)) * (-drhored_dxi(i))*vredmix**2
    end do
    do i = 1, gl%ncomp - 1
        do j = 1, gl%ncomp - 1
            help_Db(i) = help_Db(i) + gl%molfractions(j) * ((gl%vc_ij_LKP(j, i) ** 0.25d0) * gl%tc_ij_LKP(j, i) + (gl%vc_ij_LKP(i, j) ** 0.25d0) * gl%tc_ij_LKP(i, j) &
                     & - (gl%vc_ij_LKP(j, gl%ncomp) ** 0.25d0) * gl%tc_ij_LKP(j, gl%ncomp) - (gl%vc_ij_LKP(gl%ncomp, j) ** 0.25d0) * gl%tc_ij_LKP(gl%ncomp, j))
            help_D2a(i, j) = -0.25d0 * ((vredmix ** (-1.25d0)) * d2vred_dxidxj(i, j) - 1.25d0 * (vredmix ** (-2.25d0)) * (-drhored_dxi(i))*(vredmix**2) * (-drhored_dxi(j))/(gl%rhoredmix**2))
            help_D2b(i, j) = (gl%vc_ij_LKP(j, i) ** 0.25d0) * gl%tc_ij_LKP(j, i) + (gl%vc_ij_LKP(i, j) ** 0.25d0) * gl%tc_ij_LKP(i, j) &
                         & - (gl%vc_ij_LKP(j, gl%ncomp) ** 0.25d0) * gl%tc_ij_LKP(j, gl%ncomp) - (gl%vc_ij_LKP(gl%ncomp, j) ** 0.25d0) * gl%tc_ij_LKP(gl%ncomp, j) &
                         & + 2.d0 * (gl%vc_ij_LKP(gl%ncomp, gl%ncomp) ** 0.25d0) * gl%tc_ij_LKP(gl%ncomp, gl%ncomp) &
                         & - (gl%vc_ij_LKP(gl%ncomp, i) ** 0.25d0) * gl%tc_ij_LKP(gl%ncomp, i) - (gl%vc_ij_LKP(i, gl%ncomp) ** 0.25d0) * gl%tc_ij_LKP(i, gl%ncomp)
        end do
        help_Db(i) = help_Db(i) - 2.d0 * gl%molfractions(gl%ncomp) * (gl%vc_ij_LKP(gl%ncomp, gl%ncomp) ** 0.25d0) * gl%tc_ij_LKP(gl%ncomp, gl%ncomp) &
                 & + gl%molfractions(gl%ncomp) * ((gl%vc_ij_LKP(gl%ncomp, i) ** 0.25d0) * gl%tc_ij_LKP(gl%ncomp, i) + (gl%vc_ij_LKP(i, gl%ncomp) ** 0.25d0) * gl%tc_ij_LKP(i, gl%ncomp))
    end do
    do i = 1, gl%ncomp - 1
        do j = 1, gl%ncomp - 1
            d2Tred_dxidxj(i, j) = help_a * help_D2b(i, j) + help_Da(j) * help_Db(i) + help_Da(i) * help_Db(j) + help_D2a(i, j) * help_b
        end do
    end do 
   
    
elseif (gl%mix_type == 6) then     !PC-SAFT, no reducing parameters necessary. Theresa, Feb 2016
    d2Tred_dxidxj = 0.D0
    d2rhored_dxidxj = 0.D0
    return

elseif (gl%mix_type == 19) then     !One-fluid mixing rules for multiparameter EOS, no reducing parameters necessary. Andreas Jäger, July 2016
        d2Tred_dxidxj = 0.D0
    d2rhored_dxidxj = 0.D0
    return

end if
    
end subroutine d2Yr_dxidxj
!**************************************************************************

!**************************************************************************
module subroutine dndYr_dnidxj(gl,dndTr_dnidxj, dndrhor_dnidxj)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE REDUCING FUNCTIONS 
!W.R.T. ni AND xj - d(n*d(Tr)/d(ni))d(xj) AND d(n*d(rhor)/d(ni))d(xj) 
!**************************************************************************
! Johannes Gernert, Nov. 2010
!--------------------------------------------------------------------------
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! dndTr_dnidxj     -  RETURN VECTOR OF THE SIZE 30 x 30 WHICH HOLDS THE DERIVATIVES
!                     OF tred WITH RESPECT TO ALL ni AND xj
! dndrhor_dnidxj   -  RETURN VECTOR OF THE SIZE 30 x 30 WHICH HOLDS THE DERIVATIVES
!                     OF rhored WITH RESPECT TO ALL ni AND xj
!--------------------------------------------------------------------------




implicit none

    type(type_gl) :: gl


double precision, dimension(30, 30)::dndTr_dnidxj, dndrhor_dnidxj   ! Return matrices, maximum 30 components
double precision, dimension(30):: dTred_dxi, drhored_dxi 
double precision, dimension(30):: d2Tred_dxi2, d2rhored_dxi2
double precision, dimension(30, 30):: d2Tred_dxidxj, d2rhored_dxidxj
integer :: i, j, k
double precision:: xk


dndTr_dnidxj = 0.D0
dndrhor_dnidxj = 0.D0

!Andreas March 2012
!If SRK with SRK mixing rules is used, all derivatives of the reducing parameters wrt x and n are 0
if ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then
    dndTr_dnidxj = 0.D0
    dndrhor_dnidxj = 0.D0
    return
End if

!Andreas November 2015
!If PR with PR mixing rules is used, all derivatives of the reducing parameters wrt x and n are 0
if ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then
    dndTr_dnidxj = 0.D0
    dndrhor_dnidxj = 0.D0
    return
End if

if (gl%mix_type == 6) then
    dndTr_dnidxj = 0.D0
    dndrhor_dnidxj = 0.D0   
    return
End if

!Andreas Jäger, July 2016
!One-fluid mixing rules for multiparameter EOS, no reducing parameters necessary. Andreas Jäger, July 2016
if (gl%mix_type == 19) then
    dndTr_dnidxj = 0.D0
    dndrhor_dnidxj = 0.D0   
    return
End if

!Do the following for Mix_type = 1 and Mix_type = 11 and Mix_type = 12. Mixtype = 11 is for quadratic mixing rules for the reduced residual Helmholtz energy (Andreas, February 2016)

! get 1st der. of Tr and rhor w.r.t. all xi
call dYr_dxi(gl,dTred_dxi, drhored_dxi)
! get 2nd der. of Tr and rhor w.r.t. all xi
call d2Yr_dxi2(gl,d2Tred_dxi2, d2rhored_dxi2)
! get the mixed der. of Tred and rhored w.r.t.  all xi and xj 
call d2Yr_dxidxj(gl,d2Tred_dxidxj, d2rhored_dxidxj)


!Stefan & Andreas July 2014
!This is not correct and has to be updated - Preliminary fast and easy solution!!!
!Correct derivatives d2Tred_dxi2 and d2rhored_dxi2 are included in d2Yr_dxidxj!!!
if (gl%Mix_type == 4) then
    Do i=1, gl%ncomp
        d2Tred_dxi2(i) = d2Tred_dxidxj(i,i)
        d2rhored_dxi2(i) = d2rhored_dxidxj(i,i)
    End do
End if

! calculate d(n*(dTred/dni))/d(xj) and d(n*(drhored/dni))/d(xj)
! outer loop runs through all xj
do j = 1, gl%NCOMP-1
    ! inner loop runs through all ni
    do i = 1, gl%NCOMP
        if (j == i) then    !the mixed deriv. becomes the 2nd der w.r.t. xj
            dndTr_dnidxj(j, i) = d2Tred_dxi2(j) - dTred_dxi(j)
            dndrhor_dnidxj(j, i) = d2rhored_dxi2(j) - drhored_dxi(j)
        else if (i < gl%ncomp) then
            dndTr_dnidxj(j, i) = d2Tred_dxidxj(j, i) - dTred_dxi(j)
            dndrhor_dnidxj(j, i) = d2rhored_dxidxj(j, i) - drhored_dxi(j)
        else
            dndTr_dnidxj(j, i) =  - dTred_dxi(j)
            dndrhor_dnidxj(j, i) =  - drhored_dxi(j)        
        end if            
        do k = 1, gl%ncomp-1
            xk = gl%MOLFRACTIONS(k)
            if (j == k) then    !the mixed deriv. becomes the 2nd der w.r.t. xj
                dndTr_dnidxj(j, i) = dndTr_dnidxj(j, i) - xk*d2Tred_dxi2(k)
                dndrhor_dnidxj(j, i) = dndrhor_dnidxj(j, i) - xk*d2rhored_dxi2(k)
            else
                dndTr_dnidxj(j, i) = dndTr_dnidxj(j, i) - xk*d2Tred_dxidxj(j,k)
                dndrhor_dnidxj(j, i) = dndrhor_dnidxj(j, i) - xk*d2rhored_dxidxj(j,k)                
            end if
        end do
    end do
end do

end subroutine dndYr_dnidxj
!**************************************************************************





!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!NEW DERIVATIVES NEEDED FOR CRITICAL POINT CALCULATION OF A MIXTURE
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!**************************************************************************
module subroutine d3Yr_dxidxjdxk(gl,d3Tred_dxidxjdxk, d3rhored_dxidxjdxk, d2Tred_dxidxj, d2rhored_dxidxj, dTred_dxi, drhored_dxi)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE THIRD DERIVATIVE OF THE REDUCING 
!FUNCTIONS FOR DENSITY AND TEMPERATURE WITH REGARD TO THE COMPOSITION OF THE 
!COMPONENTS i, j, and k 
!(d^3(rhored)/d(xi)d(xj)d(xk)), (d^3(tred)/d(xi)d(xj)d(xk))
!All derivatives are programed according to the article
! Bell and Jäger (2016)
!Eqs. (65) - (74) of the supplementary material
!**************************************************************************
! Andreas Jäger, February 2016
!--------------------------------------------------------------------------
!
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! d3Tred_dxidxjdxk  -  RETURN MATRIX OF THE SIZE ncomp X ncomp X ncomp WHICH  
!                      HOLDS THE DERIVATIVES OF tred WITH RESPECT TO ALL COMPOSITIONS
! d3rhored_dxdxjdxk -  RETURN MATRIX OF THE SIZE ncomp X ncomp X ncomp WHICH  
!                      HOLDS THE DERIVATIVES OF tred WITH RESPECT TO ALL COMPOSITIONS
!--------------------------------------------------------------------------







implicit none

    type(type_gl) :: gl


double precision, dimension(30,30,30) :: d3Tred_dxidxjdxk, d3rhored_dxidxjdxk

double precision, dimension(30,30,30):: d3vred_dxidxjdxk
double precision, dimension(30):: dvred_dxi
double precision, dimension(30,30):: d2vred_dxidxj

!For testing purposes
double precision, dimension(30):: dTred_dxi
double precision, dimension(30,30)::d2Tred_dxidxj
double precision, dimension(30):: drhored_dxi
double precision, dimension(30,30)::d2rhored_dxidxj

!Help variables
double precision, dimension(30):: x
double precision, dimension(30):: vci
double precision, dimension(30):: vci_13
double precision, dimension(30,30):: cv_ij, cT_ij

double precision, dimension(30,30):: dfT_ij_dxi, dfT_ji_dxi, dfv_ij_dxi, dfv_ji_dxi 
double precision, dimension(30,30):: d2fT_ij_dxidxj, d2fv_ij_dxidxj
double precision, dimension(30,30):: d2fT_ij_dxi2, d2fT_ji_dxi2, d2fv_ij_dxi2, d2fv_ji_dxi2 

double precision, dimension(30,30):: d3fT_ij_dxi3, d3fT_ji_dxi3, d3fv_ij_dxi3, d3fv_ji_dxi3 
double precision, dimension(30,30):: d3fT_ij_dxidxj2, d3fv_ij_dxidxj2 
double precision, dimension(30,30):: d3fT_ij_dxi2dxj, d3fv_ij_dxi2dxj 

double precision:: help

integer::i, j, k, m

d3Tred_dxidxjdxk = 0.D0
d3rhored_dxidxjdxk = 0.D0

d3vred_dxidxjdxk = 0.D0
dvred_dxi = 0.D0
d2vred_dxidxj = 0.D0

!Write molfractions in variable xi (for shorter code)
x = gl%molfractions
!Calculate the critical volumes from the critical densities and calculate the critical volumes to the power 1/3
Do i = 1, gl%ncomp
    vci(i) = 1.D0 / gl%rhoc(i)
    vci_13(i) = vci(i)**(1.D0/3.D0)
end do
!Calculate the help variables cij
do i = 1, gl%ncomp
    do j = 1, gl%ncomp
        cv_ij(i,j) = 2.D0 * gl%rfbetarho(i,j) * gl%rfgammarho(i,j) * (vci_13(i) + vci_13(j))**3 / 8.D0  
        cT_ij(i,j) = 2.D0 * gl%rfbetat(i,j) * gl%rfgammat(i,j) * (gl%tc(i) * gl%tc(j))**0.5  
    end do
end do

!Calculate the derivatives of the help function ft and fv with respect to the molfractions for the case 
!THAT ALL x ARE INDEPENDENT!!! (This will be corrected later)
do i = 1, gl%ncomp
    do j = 1, gl%ncomp
        
        !FIRST DERIVATIVE OF THE HELP FUNCTIONS fT and frho WITH RESPECT TO xi
        !--------------------------------------------------------------------------------------------------------------------------
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
        dfT_ji_dxi(j,i) = x(j)*(x(j) + x(i))/(gl%rfbetat(j,i)**2*x(j) + x(i)) + &
            & x(j)*x(i)/(gl%rfbetat(j,i)**2*x(j) + x(i))*(1.D0 - (x(j) + x(i))/(gl%rfbetat(j,i)**2*x(j) + x(i)))
        !derivative of the composition-dependent part of the reducing function for the density w.r.t. xi:  
        dfv_ji_dxi(j,i) = x(j)*(x(j) + x(i))/(gl%rfbetarho(j,i)**2*x(j) + x(i)) + &
            & x(j)*x(i)/(gl%rfbetarho(j,i)**2*x(j) + x(i))*(1.D0 - (x(j) + x(i))/(gl%rfbetarho(j,i)**2*x(j) + x(i)))
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
        dfT_ij_dxi(i,j) = x(j)*(x(i) + x(j))/(gl%rfbetat(i,j)**2*x(i) + x(j)) + &
            & x(i)*x(j)/(gl%rfbetat(i,j)**2*x(i) + x(j))*(1.D0 - gl%rfbetat(i,j)**2*(x(i) + x(j))/(gl%rfbetat(i,j)**2*x(i) + x(j)))
        !derivative of the composition-dependent part of the reducing function for the density w.r.t. xi: 
        dfv_ij_dxi(i,j) = x(j)*(x(i) + x(j))/(gl%rfbetarho(i,j)**2*x(i) + x(j)) + &
            & x(i)*x(j)/(gl%rfbetarho(i,j)**2*x(i) + x(j))*(1.D0 - gl%rfbetarho(i,j)**2*(x(i) + x(j))/(gl%rfbetarho(i,j)**2*x(i) + x(j)))
        !--------------------------------------------------------------------------------------------------------------------------    

        
        !SECOND DERIVATIVE OF THE HELP FUNCTIONS fT and frho WITH RESPECT TO xi and xj
        !--------------------------------------------------------------------------------------------------------------------------        
        !xi = xj
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
        d2fT_ji_dxi2(j,i) = (1.D0 - (x(j) + x(i))/(gl%rfbetat(j,i)**2*x(j) + x(i))) / (gl%rfbetat(j,i)**2*x(j) + x(i)) * &
                    & (2.D0*x(j) - x(j)*x(i)*2.D0/(gl%rfbetat(j,i)**2*x(j) + x(i))) 
        !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi:  
        d2fv_ji_dxi2(j,i) = (1.D0 - (x(j) + x(i))/(gl%rfbetarho(j,i)**2*x(j) + x(i))) / (gl%rfbetarho(j,i)**2*x(j) + x(i)) * &
                    & (2.D0*x(j) - x(j)*x(i)*2.D0/(gl%rfbetarho(j,i)**2*x(j) + x(i)))     
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
        d2fT_ij_dxi2(i,j) = (1.D0 - gl%rfbetat(i,j)**2*(x(i) + x(j))/(gl%rfbetat(i,j)**2*x(i) + x(j))) / (gl%rfbetat(i,j)**2*x(i) + x(j)) * &
                    & (2.D0*x(j) - x(i)*x(j)*2.D0*gl%rfbetat(i,j)**2/(gl%rfbetat(i,j)**2*x(i) + x(j)))
        !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi: 
        d2fv_ij_dxi2(i,j) = (1.D0 - gl%rfbetarho(i,j)**2*(x(i) + x(j))/(gl%rfbetarho(i,j)**2*x(i) + x(j))) / (gl%rfbetarho(i,j)**2*x(i) + x(j)) * &
                    & (2.D0*x(j) - x(i)*x(j)*2.D0*gl%rfbetarho(i,j)**2/(gl%rfbetarho(i,j)**2*x(i) + x(j)))
        
        !xi /= xj
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi and xj:
        d2fT_ij_dxidxj(i,j) = (x(i) + x(j))/(gl%rfbetat(i,j)**2*x(i) + x(j)) + x(j)/(gl%rfbetat(i,j)**2*x(i) + x(j))*(1.D0 - (x(i) + x(j))/(gl%rfbetat(i,j)**2*x(i) + x(j))) &
                    & + x(i)/(gl%rfbetat(i,j)**2*x(i) + x(j))*(1.D0 - gl%rfbetat(i,j)**2*(x(i) + x(j))/(gl%rfbetat(i,j)**2*x(i) + x(j))) &
                    & - x(i)*x(j)/((gl%rfbetat(i,j)**2*x(i) + x(j))**2)*(1.D0 + gl%rfbetat(i,j)**2 - 2.D0*gl%rfbetat(i,j)**2*(x(i) + x(j))/(gl%rfbetat(i,j)**2*x(i) + x(j)))
        !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi and xj:       
        help = 1.d0 / (gl%rfbetarho(i,j)**2*x(i) + x(j))
        d2fv_ij_dxidxj(i,j) = (x(i) + x(j))*help &
                    & + x(j)*help*(1.D0 - (x(i) + x(j))*help) &
                    & + x(i)*help*(1.D0 - gl%rfbetarho(i,j)**2*(x(i) + x(j))*help) &
                    & - x(i)*x(j)/((gl%rfbetarho(i,j)**2*x(i) + x(j))**2)*(1.D0 + gl%rfbetarho(i,j)**2 &
                    & - 2.D0*gl%rfbetarho(i,j)**2*(x(i) + x(j))*help)  
        !--------------------------------------------------------------------------------------------------------------------------  
        

        !THIRD DERIVATIVE OF THE HELP FUNCTIONS fT and frho WITH RESPECT TO xi and xj
        !-------------------------------------------------------------------------------------------------------------------------- 
        !d3f_dxi3
        d3fT_ij_dxi3(i,j) = 6.0D0 * gl%rfbetat(i,j)**2 * x(j)**3 * (gl%rfbetat(i,j)**2 - 1.D0) / (gl%rfbetat(i,j) * x(i) + x(j))**4
        d3fT_ji_dxi3(j,i) = -6.0D0 * gl%rfbetat(j,i)**2 * x(j)**3 * (gl%rfbetat(j,i)**2 - 1.D0) / (gl%rfbetat(j,i) * x(j) + x(i))**4
        
        d3fv_ij_dxi3(i,j) = 6.0D0 * gl%rfbetarho(i,j)**2 * x(j)**3 * (gl%rfbetarho(i,j)**2 - 1.D0) / (gl%rfbetarho(i,j) * x(i) + x(j))**4
        d3fv_ji_dxi3(j,i) = -6.0D0 * gl%rfbetarho(j,i)**2 * x(j)**3 * (gl%rfbetarho(j,i)**2 - 1.D0) / (gl%rfbetarho(j,i) * x(j) + x(i))**4
        
        d3fT_ij_dxidxj2(i,j) = 6.0D0 * gl%rfbetat(i,j)**2 * x(i)**2 * x(j) * (gl%rfbetat(i,j)**2 - 1.D0) / (gl%rfbetat(i,j) * x(i) + x(j))**4
        d3fv_ij_dxidxj2(i,j) = 6.0D0 * gl%rfbetarho(i,j)**2 * x(i)**2 * x(j) * (gl%rfbetarho(i,j)**2 - 1.D0) / (gl%rfbetarho(i,j) * x(i) + x(j))**4
        
        d3fT_ij_dxi2dxj(i,j) = -6.0D0 * gl%rfbetat(i,j)**2 * x(i) * x(j)**2 * (gl%rfbetat(i,j)**2 - 1.D0) / (gl%rfbetat(i,j) * x(i) + x(j))**4
        d3fv_ij_dxi2dxj(i,j) = -6.0D0 * gl%rfbetarho(i,j)**2 * x(i) * x(j)**2 * (gl%rfbetarho(i,j)**2 - 1.D0) / (gl%rfbetarho(i,j) * x(i) + x(j))**4
        !--------------------------------------------------------------------------------------------------------------------------  

    end do
end do

!Andreas March 2012
If ((gl%Mix_type == 1) .or. (gl%Mix_type == 11) .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) Then !Mixtype = 11 is for quadratic mixing rules for the reduced residual Helmholtz energy (Andreas, February 2016),  Mixtype = 12 is for excess based departure function mixing rules (Andreas, July 2016)
    
    !---
    !Calculate dvred_dxi first, as it will be needed in for the calculation of d3rhored_dxidxjdxk
    !Also calculate dtred_dxi for a test
    do i = 1, gl%ncomp-1       !OUTER LOOP, goes through the derivatives with respect to xi (dvred_dxi)
        dvred_dxi(i) = 2.D0 * x(i) * vci(i) - 2.D0 * x(gl%ncomp) * vci(gl%ncomp)      
        dTred_dxi(i) = 2.D0 * x(i) * gl%tc(i) - 2.D0 * x(gl%ncomp) * gl%tc(gl%ncomp)        !FOR TESTING PURPOSES, CAN BE DELETED LATERON
        Do k = 1, gl%ncomp - 1
            if (k < i) then
                dvred_dxi(i) = dvred_dxi(i) + cv_ij(k,i) * dfv_ji_dxi(k,i)
                dTred_dxi(i) = dTred_dxi(i) + cT_ij(k,i) * dfT_ji_dxi(k,i)
            elseif (k > i) then
                dvred_dxi(i) = dvred_dxi(i) + cv_ij(i,k) * dfv_ij_dxi(i,k)
                dTred_dxi(i) = dTred_dxi(i) + cT_ij(i,k) * dfT_ij_dxi(i,k) 
            end if
            dvred_dxi(i) = dvred_dxi(i) - cv_ij(k,gl%ncomp) * dfv_ji_dxi(k,gl%ncomp)
            dTred_dxi(i) = dTred_dxi(i) - cT_ij(k,gl%ncomp) * dfT_ji_dxi(k,gl%ncomp)
        end do
        dvred_dxi(i) = dvred_dxi(i) + cv_ij(i,gl%ncomp) * dfv_ij_dxi(i,gl%ncomp)
        dTred_dxi(i) = dTred_dxi(i) + cT_ij(i,gl%ncomp) * dfT_ij_dxi(i,gl%ncomp)
        
        !TEST
        drhored_dxi(i) = -gl%rhoredmix**2 * dvred_dxi(i)
    end do
    !---
    
    !---
    !Calculate d2vred_dxi2, as it will also be needed in for the calculation of d3rhored_dxidxjdxk
    !Also calculate d2tred_dxi2 for a test
    do i = 1, gl%ncomp-1       !OUTER LOOP, goes through the derivatives with respect to xi (dvred_dxi)
        do j = 1, gl%ncomp-1       !OUTER LOOP, goes through the derivatives with respect to xj (dvred_dxj) 
 
            if (i == j) then
                
                d2vred_dxidxj(i,j) = 2.D0 * vci(i) + 2.D0 * vci(gl%ncomp) 
                d2Tred_dxidxj(i,j) = 2.D0 * gl%tc(i) + 2.D0 * gl%tc(gl%ncomp) 
                Do k = 1, gl%ncomp - 1
                    if (k < i) then
                        d2vred_dxidxj(i,j) = d2vred_dxidxj(i,j) + cv_ij(k,i) * d2fv_ji_dxi2(k,i)
                        d2Tred_dxidxj(i,j) = d2Tred_dxidxj(i,j) + cT_ij(k,i) * d2fT_ji_dxi2(k,i)
                    elseif (k > i) then
                        d2vred_dxidxj(i,j) = d2vred_dxidxj(i,j) + cv_ij(i,k) * d2fv_ij_dxi2(i,k)
                        d2Tred_dxidxj(i,j) = d2Tred_dxidxj(i,j) + cT_ij(i,k) * d2fT_ij_dxi2(i,k) 
                    end if
                    d2vred_dxidxj(i,j) = d2vred_dxidxj(i,j) + cv_ij(k,gl%ncomp) * d2fv_ji_dxi2(k,gl%ncomp)
                    d2Tred_dxidxj(i,j) = d2Tred_dxidxj(i,j) + cT_ij(k,gl%ncomp) * d2fT_ji_dxi2(k,gl%ncomp)
                end do
                d2vred_dxidxj(i,j) = d2vred_dxidxj(i,j) + cv_ij(i,gl%ncomp) * (-2.D0 * d2fv_ij_dxidxj(i,gl%ncomp) + d2fv_ij_dxi2(i,gl%ncomp)) 
                d2Tred_dxidxj(i,j) = d2Tred_dxidxj(i,j) + cT_ij(i,gl%ncomp) * (-2.D0 * d2fT_ij_dxidxj(i,gl%ncomp) + d2fT_ij_dxi2(i,gl%ncomp))
                
            else
                
                d2vred_dxidxj(i,j) = 2.D0 * vci(gl%ncomp) 
                d2Tred_dxidxj(i,j) = 2.D0 * gl%tc(gl%ncomp) 
                Do k = 1, gl%ncomp - 1
                    d2vred_dxidxj(i,j) = d2vred_dxidxj(i,j) + cv_ij(k,gl%ncomp) * d2fv_ji_dxi2(k,gl%ncomp)
                    d2Tred_dxidxj(i,j) = d2Tred_dxidxj(i,j) + cT_ij(k,gl%ncomp) * d2fT_ji_dxi2(k,gl%ncomp)
                end do
                d2vred_dxidxj(i,j) = d2vred_dxidxj(i,j) + cv_ij(i,j) * d2fv_ij_dxidxj(i,j) - cv_ij(j,gl%ncomp) * d2fv_ij_dxidxj(j,gl%ncomp) - cv_ij(i,gl%ncomp) * d2fv_ij_dxidxj(i,gl%ncomp)  
                d2Tred_dxidxj(i,j) = d2Tred_dxidxj(i,j) + cT_ij(i,j) * d2fT_ij_dxidxj(i,j) - cT_ij(j,gl%ncomp) * d2fT_ij_dxidxj(j,gl%ncomp) - cT_ij(i,gl%ncomp) * d2fT_ij_dxidxj(i,gl%ncomp) 
                
            end if
            
            !TEST
            d2rhored_dxidxj(i,j) = 2.D0/gl%rhoredmix*drhored_dxi(i)*drhored_dxi(j) - gl%rhoredmix**2*d2vred_dxidxj(i,j)
        end do
    end do
    !---
    
    
    do i = 1, gl%ncomp-1
        do j = 1, gl%ncomp-1
            do k = 1, gl%ncomp-1
                
                if ((i == k) .and. (j == k)) then
                    
                elseif ((i /= k) .and. (j /= k) .and. (i /= j)) then
                    
                elseif ((j == k) .and. (i /= k)) then
                    
                elseif ((i == k) .and. (j /= k)) then
                
                elseif ((i == j) .and. (i /= k)) then
                    !d3vred_dxidxjdxk(i,j,k)
                end if
                
                
                !Calculate d3rhored_dxidxjdxk from d3vred_dxidxjdxk 
                d3rhored_dxidxjdxk(i,j,k) = -gl%rhoredmix**2 * d3vred_dxidxjdxk(i,j,k)  & 
                & + 2.D0 * gl%rhoredmix**3 * (dvred_dxi(k) * d2vred_dxidxj(i,j) &
                & + dvred_dxi(j) * d2vred_dxidxj(i,k) + dvred_dxi(i) * d2vred_dxidxj(j,k) ) &
                & - 6.D0 * gl%rhoredmix**4 * dvred_dxi(i) * dvred_dxi(j) * dvred_dxi(k)
            end do
        end do
    end do
    
    
!If SRK with SRK mixing rules or PR with PR mixing rules is used, all derivatives of the reducing parameters wrt x are 0
Else if ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22) .OR. (gl%Mix_type == 3) .OR. (gl%Mix_type == 31)) then
    d3Tred_dxidxjdxk = 0.D0
    d3rhored_dxidxjdxk = 0.D0
    return
    
!LKP 
Else if (gl%Mix_type == 4) then
    
    !NOT YET IMPLEMENTED!!
    !Calculate d3rhored_dxidxjdxk from d3vred_dxidxjdxk 
    do i = 1, gl%ncomp-1
        do j = 1, gl%ncomp-1
            do k = 1, gl%ncomp-1
                d3rhored_dxidxjdxk(i,j,k) = -gl%rhoredmix**2 * d3vred_dxidxjdxk(i,j,k)  & 
                & + 2.D0 * gl%rhoredmix**3 * (dvred_dxi(k) * d2vred_dxidxj(i,j) &
                & + dvred_dxi(j) * d2vred_dxidxj(i,k) + dvred_dxi(i) * d2vred_dxidxj(j,k) ) &
                & - 6.D0 * gl%rhoredmix**4 * dvred_dxi(i) * dvred_dxi(j) * dvred_dxi(k)
            end do
        end do
    end do
    
    
elseif (gl%mix_type == 6) then     
    d3Tred_dxidxjdxk = 0.D0
    d3rhored_dxidxjdxk = 0.D0
    return
    
!One-fluid mixing rules for multiparameter EOS, no reducing parameters necessary. Andreas Jäger, July 2016   
elseif (gl%mix_type == 19) then     
    d3Tred_dxidxjdxk = 0.D0
    d3rhored_dxidxjdxk = 0.D0
    return
    
end if
    
end subroutine d3Yr_dxidxjdxk
!**************************************************************************

!**************************************************************************
module subroutine d2ndYr_dnidxjdxk(gl,d2ndTr_dnidxjdxk, d2ndrhor_dnidxjdxk)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE REDUCING FUNCTIONS 
!W.R.T. ni, xj, and xk - 
! d2(n*d(Tr)/d(ni))d(xj)d(xk) AND d2(n*d(rhor)/d(ni))d(xj)d(xk) 
!**************************************************************************
! Andreas Jäger, February 2016
!--------------------------------------------------------------------------
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! d2ndTr_dnidxjdxk     -  RETURN VECTOR OF THE SIZE 30 x 30 x 30 WHICH 
!                       HOLDS THE DERIVATIVES
!                       OF tred WITH RESPECT TO ALL ni AND xj
! d2ndrhor_dnidxjdxk   -  RETURN VECTOR OF THE SIZE 30 x 30 x 30 WHICH 
!                       HOLDS THE DERIVATIVES
!                       OF rhored WITH RESPECT TO ALL ni AND xj
!--------------------------------------------------------------------------




implicit none

    type(type_gl) :: gl


double precision, dimension(30, 30, 30)::d2ndTr_dnidxjdxk, d2ndrhor_dnidxjdxk   ! Return matrices, maximum 30 components


double precision, dimension(30):: d2Tred_dxi2, d2rhored_dxi2
double precision, dimension(30, 30):: d2Tred_dxidxj, d2rhored_dxidxj
double precision, dimension(30,30,30):: d3Tred_dxidxjdxk, d3rhored_dxidxjdxk
integer :: i, j, k, m

d2ndTr_dnidxjdxk = 0.D0
d2ndrhor_dnidxjdxk = 0.D0

!Andreas March 2012
!If SRK with SRK mixing rules is used, all derivatives of the reducing parameters wrt x and n are 0
if ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then
    d2ndTr_dnidxjdxk = 0.D0
    d2ndrhor_dnidxjdxk = 0.D0
    return
End if

!Andreas November 2015
!If PR with PR mixing rules is used, all derivatives of the reducing parameters wrt x and n are 0
if ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then
    d2ndTr_dnidxjdxk = 0.D0
    d2ndrhor_dnidxjdxk = 0.D0
    return
End if

if (gl%mix_type == 6) then
    d2ndTr_dnidxjdxk = 0.D0
    d2ndrhor_dnidxjdxk = 0.D0   
    return
end if

!Andreas Jäger, July 2016
!One-fluid mixing rules for multiparameter EOS, no reducing parameters necessary. Andreas Jäger, July 2016
if (gl%mix_type == 19) then
    d2ndTr_dnidxjdxk = 0.D0
    d2ndrhor_dnidxjdxk = 0.D0   
    return
End if

!Do the following for Mix_type = 1 and Mix_type = 11. Mixtype = 11 is for quadratic mixing rules for the reduced residual Helmholtz energy (Andreas, February 2016)


! get the mixed der. of Tred and rhored w.r.t.  all xi and xj 
call d2Yr_dxidxj(gl,d2Tred_dxidxj, d2rhored_dxidxj)
! Get the mixed der. of Tred and rhored w.r.t. all xi, xj, and xk
!call d3Yr_dxidxjdxk(d3Tred_dxidxjdxk, d3rhored_dxidxjdxk)

!For Helmholtz mixtures write the diagonal elements of the derivative matrix on the other matrix
if ((gl%Mix_type == 1) .or. (gl%Mix_type == 11)  .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) then
    ! get 2nd der. of Tr and rhor w.r.t. all xi
    call d2Yr_dxi2(gl,d2Tred_dxi2, d2rhored_dxi2)
    Do i=1, gl%ncomp
        d2Tred_dxidxj(i,i) = d2Tred_dxi2(i) 
        d2rhored_dxidxj(i,i) = d2rhored_dxi2(i)
    End do
End if

do i = 1, gl%ncomp !Loop over ni
    do j = 1, gl%ncomp-1 !Loop over xj
        do k = 1, gl%ncomp-1 !Loop over xk
            d2ndTr_dnidxjdxk(i,j,k) = d3Tred_dxidxjdxk(i,j,k) - 2.D0 * d2Tred_dxidxj(j,k)
            d2ndrhor_dnidxjdxk(i,j,k) = d3rhored_dxidxjdxk(i,j,k) - 2.D0 * d2rhored_dxidxj(j,k)
            do m = 1, gl%ncomp-1 !Loop over mk
                d2ndTr_dnidxjdxk(i,j,k) = d2ndTr_dnidxjdxk(i,j,k) - gl%molfractions(m) * d3Tred_dxidxjdxk(j,k,m)       
                d2ndrhor_dnidxjdxk(i,j,k) = d2ndrhor_dnidxjdxk(i,j,k) - gl%molfractions(m) * d3rhored_dxidxjdxk(j,k,m)  
            end do
        end do
    end do
end do

end subroutine d2ndYr_dnidxjdxk
!**************************************************************************


!**************************************************************************
module subroutine Help_PSI_Derivs(gl)
!**************************************************************************
!Routine that calculates the help variables PSI_rho and PSI_T according
!to the paper of Bell and Jäger (2016). The equations can be found in
!the supplementary material of the article, page4, Eqs. (26 - 33)
!Andreas Jäger, February 2016
!--------------------------------------------------------------------------
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! PSI_Y         - Help variables for Tred and rhored  
! dPSI_Y_dxj    - Derivatives of help variables wrt xj 
! dPSI_Y_dxjdxk - Derivatives of help variables wrt xj and xk
!--------------------------------------------------------------------------





implicit none

    type(type_gl) :: gl



!Derivatives of reducing functions wrt to ni and xj (and xk) which are needed for the help variables
double precision, dimension(30):: ndTr_dni, ndrhor_dni
double precision, dimension(:,:),allocatable:: dndTr_dnidxj, dndrhor_dnidxj
double precision, dimension(:,:,:),allocatable:: d2ndTr_dnidxjdxk, d2ndrhor_dnidxjdxk

!Derivatives of reducing functions wrt to xi and xj and xk which are needed for the help variables
double precision, dimension(30):: dTr_dxi, drhor_dxi, d2Tr_dxi2, d2rhor_dxi2
double precision, dimension(:,:),allocatable:: d2Tr_dxidxj, d2rhor_dxidxj
double precision, dimension(:,:,:),allocatable:: d3Tr_dxidxjdxk, d3rhor_dxidxjdxk

integer:: i,j,k

if(.not. allocated(d2ndTr_dnidxjdxk)) then
    allocate(d2ndTr_dnidxjdxk(30,30,30))
    allocate(d2ndrhor_dnidxjdxk,d3Tr_dxidxjdxk,d3rhor_dxidxjdxk,mold=d2ndTr_dnidxjdxk)
    allocate(dndTr_dnidxj(30,30))
    allocate(dndrhor_dnidxj,d2Tr_dxidxj, d2rhor_dxidxj,mold=dndTr_dnidxj)
endif
!Get the necessary derivatives
!--
call ndYr_dni(gl,ndTr_dni, ndrhor_dni)
call dndYr_dnidxj(gl,dndTr_dnidxj, dndrhor_dnidxj)
call d2ndYr_dnidxjdxk(gl,d2ndTr_dnidxjdxk, d2ndrhor_dnidxjdxk)

call dYr_dxi(gl,dTr_dxi, drhor_dxi)
call d2Yr_dxidxj(gl,d2Tr_dxidxj, d2rhor_dxidxj)
!call d3Yr_dxidxjdxk(d3Tr_dxidxjdxk, d3rhor_dxidxjdxk)
!--

if ((gl%mix_type == 1) .or. (gl%mix_type == 11) .or. (gl%Mix_type == 12) .or. (gl%mix_type == 13)) then
    call d2Yr_dxi2(gl,d2Tr_dxi2, d2rhor_dxi2)
    do i = 1, gl%ncomp-1
        d2Tr_dxidxj(i,i) = d2Tr_dxi2(i)
        d2rhor_dxidxj(i,i) = d2rhor_dxi2(i)
    end do  
end if
gl%PSI_Y(1,:) = 1.D0 / gl%tredmix * ndTr_dni
gl%PSI_Y(2,:) = 1.D0 - 1.D0 / gl%rhoredmix * ndrhor_dni

do i = 1, gl%ncomp
    do j = 1, gl%ncomp-1
        gl%dPSI_Y_dxj(1,i,j) = -1.D0 / gl%tredmix * (dndTr_dnidxj(i,j) - dTr_dxi(j) * gl%PSI_Y(1,i))
        gl%dPSI_Y_dxj(2,i,j) = -1.D0 / gl%rhoredmix * (dndrhor_dnidxj(i,j) + drhor_dxi(j) * (1.D0 - gl%PSI_Y(2,i)))
        do k = 1, gl%ncomp-1
            gl%d2PSI_Y_dxjdxk(1,i,j,k) = 1.D0 / gl%tredmix * d2ndTr_dnidxjdxk(i,j,k) &
                        & - 1.D0 / gl%tredmix**2 * dTr_dxi(k) * dndTr_dnidxj(i,j) &
                        & - 1.D0 / gl%tredmix * dTr_dxi(j) * gl%dPSI_Y_dxj(1,i,k) &
                        & - (1.D0 / gl%tredmix * d2Tr_dxidxj(j,k) - 1.D0 / gl%tredmix**2 * dTr_dxi(k) * dTr_dxi(j)) * gl%PSI_Y(1,i)
            gl%d2PSI_Y_dxjdxk(2,i,j,k) = 1.D0 / gl%rhoredmix * d2ndrhor_dnidxjdxk(i,j,k) &
                        & + 1.D0 / gl%rhoredmix**2 * drhor_dxi(k) * dndrhor_dnidxj(i,j) &
                        & - 1.D0 / gl%rhoredmix * drhor_dxi(j) * gl%dPSI_Y_dxj(2,i,k) &
                        & + (1.D0 / gl%rhoredmix * d2rhor_dxidxj(j,k) &
                        & - 1.D0 / gl%rhoredmix**2 * drhor_dxi(k) * drhor_dxi(j)) * (1.D0 - gl%PSI_Y(2,i))
        end do
    end do
end do

end subroutine Help_PSI_Derivs        
!**************************************************************************


!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!END OF NEW DERIVATIVES NEEDED FOR CRITICAL POINT CALCULATION OF A MIXTURE
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   Andreas, November 2013
!   Additional derivatives implemented that are needed for the derivative 
!   of ln(fug_coef) w.r.t ni at constant T,p,nk.
!   Not needed for our phase equilibrium routines, implemented for SINTEF

!**************************************************************************
module subroutine ndYr_dni_old(gl,ndTred_dni, ndrhored_dni)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE FIRST DERIVATIVE OF THE REDUCING 
!FUNCTIONS FOR DENSITY AND TEMPERATURE WITH REGARD TO THE AMOUNT OF SUBSTANCE
! OF THE COMPONENT i: n*(d rhored/d ni)_nj, n*(d tred/d ni)_nj
!**************************************************************************
! Johannes Gernert, Bochum, July 2010
!--------------------------------------------------------------------------
!THE CALCULATION IS BASED ON THE EQUATIONS PUBLISHED IN THE GERG-2004 
!TECHNICAL MONOGRAPH 15, p. 116:
!                ****************************************
!                *   Kunz, O. et al.,                   *
!                *   The GERG-2004 Wide-Range Equation  * 
!                *   of State for Natural Gases and     * 
!                *   Other Mixtures                     *
!                *   GERG TM 15, 2007                   *
!                ****************************************
!--------------------------------------------------------------------------
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! ndTred_dni     -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF tred WITH RESPECT TO ALL AMOUNTS OF SUBSTANCE
! ndrhored_dni   -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF rhored WITH RESPECT TO ALL  AMOUNTS OF SUBSTANCE
!--------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl


double precision, dimension(30)::ndTred_dni, ndrhored_dni   ! Return vectors, maximum 30 components
double precision, dimension(30)::dTred_dxi, drhored_dxi   
integer :: i, k
double precision:: xk

do i = 1, 30
    ndTred_dni(i) = 0.D0
    ndrhored_dni(i) = 0.D0
end do

call dYr_dxi_old(gl,dTred_dxi, drhored_dxi)

! calculate n*(drhored/dni) = drhored/dxi - SUM(xk*(drhored/dxk))
do i = 1, gl%NCOMP
    ndTred_dni(i) = dTred_dxi(i)
    ndrhored_dni(i) = drhored_dxi(i)
    do k = 1, gl%NCOMP
        xk = gl%MOLFRACTIONS(k)
        ndTred_dni(i) = ndTred_dni(i) - xk*dTred_dxi(k)
        ndrhored_dni(i) = ndrhored_dni(i) - xk*drhored_dxi(k)
    end do
end do

end subroutine ndYr_dni_old
!**************************************************************************

!**************************************************************************
module subroutine dYr_dxi_old(gl,dTred_dxi, drhored_dxi)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE FIRST DERIVATIVE OF THE REDUCING 
!FUNCTIONS FOR DENSITY AND TEMPERATURE WITH REGARD TO THE MOLE FRACTION
! OF THE COMPONENT i: d(rhored)/d(xi)_xj, d(tred)/d(xi)_xj
!**************************************************************************
! Johannes Gernert, Boulder, Aug. 2010
!--------------------------------------------------------------------------
!THE CALCULATION IS BASED ON THE EQUATIONS PUBLISHED IN THE GERG-2004 
!TECHNICAL MONOGRAPH 15, p. 116:
!                ****************************************
!                *   Kunz, O. et al.,                   *
!                *   The GERG-2004 Wide-Range Equation  * 
!                *   of State for Natural Gases and     * 
!                *   Other Mixtures                     *
!                *   GERG TM 15, 2007                   *
!                ****************************************
!--------------------------------------------------------------------------
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! dTred_dxi      -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVATIVES
!                   OF tred WITH RESPECT TO ALL MOLAR FRACTIONS
! drhored_dxi    -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF rhored WITH RESPECT TO ALL  MOLAR FRACTIONS
!--------------------------------------------------------------------------





implicit none

    type(type_gl) :: gl


double precision, dimension(30)::dTred_dxi, drhored_dxi  
double precision, dimension(30)::dvred_dxi
integer::i, j
double precision:: xi, xj, rhoci, rhocj, tci, tcj
! parameters for the 1st loop
double precision:: betatji, betarhoji, gammatji, gammarhoji
double precision:: ct_ji, crho_ji
double precision:: dftji_dxi, dfrhoji_dxi
! parameters for the 2nd loop
double precision:: betatij, betarhoij, gammatij, gammarhoij
double precision:: ct_ij, crho_ij
double precision:: dftij_dxi, dfrhoij_dxi

do i = 1, gl%ncomp
    dTred_dxi(i) = 0.D0
    drhored_dxi(i) = 0.D0
    dvred_dxi(i) = 0.D0
end do

do i = 1, gl%NCOMP
    rhoci = gl%rhoc(i)
    tci = gl%tc(i)
    xi = gl%MOLFRACTIONS(i)
    dTred_dxi(i) = 2.D0*xi*tci
    dvred_dxi(i) = 2.D0*xi/rhoci
    ! 1st loop for j < i
    do j = 1, i-1
        xj = gl%MOLFRACTIONS(j)
        rhocj = gl%rhoc(j)
        tcj = gl%tc(j)
        betatji = gl%RFBETAT(j, i)
        betarhoji = gl%RFBETARHO(j, i)
        gammatji = gl%RFGAMMAT(j, i)
        gammarhoji = gl%RFGAMMARHO(j, i)
        ct_ji = 2.D0*betatji*gammatji*(tci*tcj)**0.5D0      !critical term, temperature
        crho_ji = 2.D0*betarhoji*gammarhoji*0.125d0*(1.D0/(rhoci)**(1.D0/3.D0) + 1.D0/(rhocj)**(1.D0/3.D0))**3  !critical term, density
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
        dftji_dxi = xj*(xj + xi)/(betatji**2*xj + xi) + &
                    & xj*xi/(betatji**2*xj + xi)*(1.D0 - (xj + xi)/(betatji**2*xj + xi))
        !derivative of the composition-dependent part of the reducing function for the density w.r.t. xi:  
        dfrhoji_dxi = xj*(xj + xi)/(betarhoji**2*xj + xi) + &
                    & xj*xi/(betarhoji**2*xj + xi)*(1.D0 - (xj + xi)/(betarhoji**2*xj + xi))
        ! 1st derivative of Tred w.r.t xi - 1st summation for j = 1 to i-1
        dTred_dxi(i) = dTred_dxi(i) +  ct_ji * dftji_dxi
        ! 1st derivative of rhored w.r.t xi - 1st summation for j = 1 to i-1
        dvred_dxi(i) = dvred_dxi(i) +  crho_ji * dfrhoji_dxi
    end do
    !2nd loop for j > i
    do j = i + 1, gl%NCOMP
        xj = gl%MOLFRACTIONS(j)
        rhocj = gl%rhoc(j)
        tcj = gl%tc(j)        
        betatij = gl%RFBETAT(i, j)
        betarhoij = gl%RFBETARHO(i, j)
        gammatij = gl%RFGAMMAT(i, j)
        gammarhoij = gl%RFGAMMARHO(i, j)
        ct_ij = 2.D0*betatij*gammatij*(tci*tcj)**0.5D0      !critical term, temperature
        crho_ij = 2.D0*betarhoij*gammarhoij*0.125d0*(1.D0/(rhoci)**(1.D0/3.D0) + 1.D0/(rhocj)**(1.D0/3.D0))**3  !critical term, density
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
        dftij_dxi = xj*(xi + xj)/(betatij**2*xi + xj) + &
                 & xi*xj/(betatij**2*xi + xj)*(1 - betatij**2*(xi + xj)/(betatij**2*xi + xj))
        !derivative of the composition-dependent part of the reducing function for the density w.r.t. xi: 
        dfrhoij_dxi = xj*(xi + xj)/(betarhoij**2*xi + xj) + &
                   & xi*xj/(betarhoij**2*xi + xj)*(1 - betarhoij**2*(xi + xj)/(betarhoij**2*xi + xj))
        ! 1st derivative of Tred w.r.t xi - 2nd summation for j = i+1 to ncomp
        dTred_dxi(i) = dTred_dxi(i) +  ct_ij * dftij_dxi
        ! 1st derivative of rhored w.r.t xi - 2nd summation for j = i+1 to ncomp
        dvred_dxi(i) = dvred_dxi(i) +  crho_ij * dfrhoij_dxi
    end do
    drhored_dxi(i) = - gl%rhoredmix**2*dvred_dxi(i) ! drhored/dxi = -rhored^2*(d(1/rhored)dxi) - the routine above actually calculates dvred/dxi ...
end do

end subroutine dYr_dxi_old
!**************************************************************************

!**************************************************************************
module subroutine d2Yr_dxi2_old(gl,d2Tred_dxi2, d2rhored_dxi2)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE REDUCING 
!FUNCTIONS FOR DENSITY AND TEMPERATURE WITH REGARD TO THE COMPOSITION OF THE 
!COMPONENT i (d^2 rhored/d xi^2)_xj, (d^2 tred/d xi^2)_xj
!**************************************************************************
! Johannes Gernert, Bochum, July 2010
!--------------------------------------------------------------------------
!THE CALCULATION IS BASED ON THE EQUATIONS PUBLISHED IN THE GERG-2004 
!TECHNICAL MONOGRAPH 15, p. 116:
!                ****************************************
!                *   Kunz, O. et al.,                   *
!                *   The GERG-2004 Wide-Range Equation  * 
!                *   of State for Natural Gases and     * 
!                *   Other Mixtures                     *
!                *   GERG TM 15, 2007                   *
!                ****************************************
!--------------------------------------------------------------------------
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! d2Tred_dxi2     -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF tred WITH RESPECT TO ALL COMPOSITIONS
! d2rhored_dxi2   -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF rhored WITH RESPECT TO ALL COMPOSITIONS
!--------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl



double precision, dimension(30)::d2Tred_dxi2, d2rhored_dxi2   ! Return vectors, maximum 30 components
double precision, dimension(30)::dTred_dxi, drhored_dxi         ! Return vectors, 1st derivative, maximum 30 components
double precision, dimension(30)::d2vred_dxi2                    ! internal vector for vred
integer::i, j
double precision:: xi, xj, rhoci, rhocj, tci, tcj
! parameters for the 1st loop
double precision:: betatji, betarhoji, gammatji, gammarhoji
double precision:: ct_ji, crho_ji
double precision:: d2ftji_dxi2, d2frhoji_dxi2
! parameters for the 2nd loop
double precision:: betatij, betarhoij, gammatij, gammarhoij
double precision:: ct_ij, crho_ij
double precision:: d2ftij_dxi2, d2frhoij_dxi2

do i = 1, gl%ncomp
    d2Tred_dxi2(i) = 0.D0
    d2vred_dxi2(i) = 0.D0
    d2rhored_dxi2(i) = 0.D0
    dTred_dxi(i) = 0.D0
    drhored_dxi(i) = 0.D0
end do

! the first derivative of rhored w.r.t. xi is needed to calculate d^2rhored/dxi^2
call dYr_dxi_old(gl,dTred_dxi, drhored_dxi)

do i = 1, gl%NCOMP
    rhoci = gl%rhoc(i)
    tci = gl%tc(i)
    xi = gl%MOLFRACTIONS(i)
    d2Tred_dxi2(i) = 2.D0*tci
    d2vred_dxi2(i) = 2.D0/rhoci
    ! 1st loop for j < i
    do j = 1, i-1
        xj = gl%MOLFRACTIONS(j)
        rhocj = gl%rhoc(j)
        tcj = gl%tc(j)
        betatji = gl%RFBETAT(j, i)
        betarhoji = gl%RFBETARHO(j, i)
        gammatji = gl%RFGAMMAT(j, i)
        gammarhoji = gl%RFGAMMARHO(j, i)
        ct_ji = 2.D0*betatji*gammatji*(tci*tcj)**0.5D0      !critical term, temperature
        crho_ji = 2.D0*betarhoji*gammarhoji*0.125d0*(1.D0/(rhoci)**(1.D0/3.D0) + 1.D0/(rhocj)**(1.D0/3.D0))**3  !critical term, density
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
        d2ftji_dxi2 = (1.D0 - (xj + xi)/(betatji**2*xj + xi)) / (betatji**2*xj + xi) * &
                  & (2.D0*xj - xj*xi*2.D0/(betatji**2*xj + xi)) 
        !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi:  
        d2frhoji_dxi2 = (1.D0 - (xj + xi)/(betarhoji**2*xj + xi)) / (betarhoji**2*xj + xi) * &
                  & (2.D0*xj - xj*xi*2.D0/(betarhoji**2*xj + xi)) 
        ! 2nd derivative of Tred w.r.t xi - 1st summation for j = 1 to i-1
        d2Tred_dxi2(i) = d2Tred_dxi2(i) +  ct_ji * d2ftji_dxi2
        ! 2nd derivative of vred w.r.t xi - 1st summation for j = 1 to i-1
        d2vred_dxi2(i) = d2vred_dxi2(i) +  crho_ji * d2frhoji_dxi2
    end do
    !2nd loop for j > i
    do j = i + 1, gl%NCOMP
        xj = gl%MOLFRACTIONS(j)
        rhocj = gl%rhoc(j)
        tcj = gl%tc(j)        
        betatij = gl%RFBETAT(i, j)
        betarhoij = gl%RFBETARHO(i, j)
        gammatij = gl%RFGAMMAT(i, j)
        gammarhoij = gl%RFGAMMARHO(i, j)
        ct_ij = 2.D0*betatij*gammatij*(tci*tcj)**0.5D0      !critical term, temperature
        crho_ij = 2.D0*betarhoij*gammarhoij*0.125d0*(1.D0/(rhoci)**(1.D0/3.D0) + 1.D0/(rhocj)**(1.D0/3.D0))**3  !critical term, density
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi:
        d2ftij_dxi2 = (1.D0 - betatij**2*(xi + xj)/(betatij**2*xi + xj)) / (betatij**2*xi + xj) * &
                   & (2.D0*xj - xi*xj*2.D0*betatij**2/(betatij**2*xi + xj))
        !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi: 
        d2frhoij_dxi2 = (1.D0 - betarhoij**2*(xi + xj)/(betarhoij**2*xi + xj)) / (betarhoij**2*xi + xj) * &
                   & (2.D0*xj - xi*xj*2.D0*betarhoij**2/(betarhoij**2*xi + xj))
        ! 2nd derivative of Tred w.r.t xi - 2nd summation for j = i+1 to ncomp
        d2Tred_dxi2(i) = d2Tred_dxi2(i) +  ct_ij * d2ftij_dxi2
        ! 2nd derivative of vred w.r.t xi - 2nd summation for j = i+1 to ncomp
        d2vred_dxi2(i) = d2vred_dxi2(i) +  crho_ij * d2frhoij_dxi2
    end do
    ! the routine above actually calculates d^2vred/dxi^2, therefore must be transformed to:
    ! d^2rhored/dxi^2 = 2/rhored*(drhored/dxi)^2 - rhored^2*(d^2rhored/dxi^2)
    d2rhored_dxi2(i) = 2.D0/gl%rhoredmix*(drhored_dxi(i))**2 - gl%rhoredmix**2*d2vred_dxi2(i) 
end do

! calculation of 
end subroutine d2Yr_dxi2_old
!**************************************************************************

!**************************************************************************
module subroutine d2Yr_dxidxj_old(gl,d2Tred_dxidxj, d2rhored_dxidxj)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE SECOND DERIVATIVE OF THE REDUCING 
!FUNCTIONS FOR DENSITY AND TEMPERATURE WITH REGARD TO THE COMPOSITION OF THE 
!COMPONENTS i AND j (d^2(rhored)/d(xi)d(xj)), (d^2(tred)/d(xi)d(xj))
!**************************************************************************
! Johannes Gernert, Boulder, Aug. 2010
!--------------------------------------------------------------------------
!THE CALCULATION IS BASED ON THE EQUATIONS PUBLISHED IN THE GERG-2004 
!TECHNICAL MONOGRAPH 15, p. 116:
!                ****************************************
!                *   Kunz, O. et al.,                   *
!                *   The GERG-2004 Wide-Range Equation  * 
!                *   of State for Natural Gases and     * 
!                *   Other Mixtures                     *
!                *   GERG TM 15, 2007                   *
!                ****************************************
!--------------------------------------------------------------------------
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! d2Tred_dxidxj  -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF tred WITH RESPECT TO ALL COMPOSITIONS
! d2rhored_dxdxj -  RETURN VECTOR OF THE SIZE ncomp WHICH HOLDS THE DERIVETIVES
!                   OF rhored WITH RESPECT TO ALL COMPOSITIONS
!--------------------------------------------------------------------------






implicit none

    type(type_gl) :: gl


double precision, dimension(30,30):: d2Tred_dxidxj, d2rhored_dxidxj
double precision, dimension(30,30):: d2vred_dxidxj
double precision, dimension(30)::dTred_dxi, drhored_dxi         ! 1st derivative, maximum 30 components
double precision::d2fvij_dxidxj, d2fTij_dxidxj, cv_ij, cT_ij
double precision::betatij, betarhoij, gammatij, gammarhoij
double precision:: xi, xj, rhoci, rhocj, tci, tcj, help
integer::i, j


do i = 1, gl%ncomp
    do j = 1, gl%ncomp
        d2Tred_dxidxj(i, j) = 0.D0
        d2vred_dxidxj(i, j) = 0.D0
        d2rhored_dxidxj(i, j) = 0.D0
    end do
end do

! the first derivative of rhored w.r.t. xi is needed to calculate d^2(rhored)/d(xi)d(xj)
call dYr_dxi_old(gl,dTred_dxi, drhored_dxi)

do i = 1, (gl%ncomp - 1)
    do j = (i + 1), gl%ncomp
        betatij = gl%RFBETAT(i,j)
        betarhoij = gl%RFBETARHO(i,j)
        gammatij = gl%RFGAMMAT(i,j)
        gammarhoij = gl%RFGAMMARHO(i,j)
        xi = gl%MOLFRACTIONS(i)
        xj = gl%MOLFRACTIONS(j)
        tci = gl%tc(i)
        tcj = gl%tc(j)
        rhoci = gl%rhoc(i)
        rhocj = gl%rhoc(j)
        !critical term, temperature
        cT_ij = 2.D0*betatij*gammatij*(tci*tcj)**0.5D0
        !critical term, volume
        cv_ij = 2.D0*betarhoij*gammarhoij*0.125d0*(1.D0/(rhoci**(1.D0/3.D0)) + 1.D0/(rhocj**(1.D0/3.D0)))**3
        !derivative of the composition-dependent part of the reducing function for the temperature w.r.t. xi and xj:
        d2fTij_dxidxj = (xi + xj)/(betatij**2*xi + xj) + xj/(betatij**2*xi + xj)*(1.D0 - (xi + xj)/(betatij**2*xi + xj)) &
                    & + xi/(betatij**2*xi + xj)*(1.D0 - betatij**2*(xi + xj)/(betatij**2*xi + xj)) &
                    & - xi*xj/((betatij**2*xi + xj)**2)*(1.D0 + betatij**2 - 2.D0*betatij**2*(xi + xj)/(betatij**2*xi + xj))
        !derivative of the composition-dependent part of the reducing function for the volume w.r.t. xi and xj:
        help = 1.d0/(betarhoij**2*xi + xj)
        d2fvij_dxidxj = (xi + xj)*help &
                    & + xj*help*(1.D0 - (xi + xj)*help) &
                    & + xi*help*(1.D0 - betarhoij**2*(xi + xj)*help) &
                    & - xi*xj/((betarhoij**2*xi + xj)**2)*(1.D0 + betarhoij**2 &
                    & - 2.D0*betarhoij**2*(xi + xj)*help)
        !2nd derivative of Tred w.r.t xi and xj
        d2Tred_dxidxj(i,j) = cT_ij*d2fTij_dxidxj
        d2Tred_dxidxj(j,i) = d2Tred_dxidxj(i,j)
        !2nd derivative of vred w.r.t xi and xj
        d2vred_dxidxj(i,j) = cv_ij*d2fvij_dxidxj
        d2vred_dxidxj(j,i) = d2vred_dxidxj(i,j)
        ! calculate d^2(rhored)/d(xi)d(xj) from d^2(1/rhored)/d(xi)d(xj)
        d2rhored_dxidxj(i,j) = 2.D0/gl%rhoredmix*drhored_dxi(i)*drhored_dxi(j) - gl%rhoredmix**2*d2vred_dxidxj(i,j)
        d2rhored_dxidxj(j,i) = d2rhored_dxidxj(i,j)
    end do
end do

end subroutine d2Yr_dxidxj_old
!**************************************************************************

!**************************************************************************
module subroutine dndYr_dnidxj_old(gl,dndTr_dnidxj, dndrhor_dnidxj)
!**************************************************************************
!SUBROUTINE FOR THE CALCULATION OF THE DERIVATIVE OF THE REDUCING FUNCTIONS 
!W.R.T. ni AND xj - d(n*d(Tr)/d(ni))d(xj) AND d(n*d(rhor)/d(ni))d(xj) 
!**************************************************************************
! Johannes Gernert, Aug. 2010
!--------------------------------------------------------------------------
!THE CALCULATION IS BASED ON THE EQUATIONS PUBLISHED IN THE GERG-2004 
!TECHNICAL MONOGRAPH 15, p. 122, EQU. 7.55 AND 7.56 
!                ****************************************
!                *   Kunz, O. et al.,                   *
!                *   The GERG-2004 Wide-Range Equation  * 
!                *   of State for Natural Gases and     * 
!                *   Other Mixtures                     *
!                *   GERG TM 15, 2007                   *
!                ****************************************
!--------------------------------------------------------------------------
! INPUT PARAMETERS:
! - NONE
!
! OUTPUT PARAMETERS
! dndTr_dnidxj     -  RETURN VECTOR OF THE SIZE 30 x 30 WHICH HOLDS THE DERIVATIVES
!                     OF tred WITH RESPECT TO ALL ni AND xj
! dndrhor_dnidxj   -  RETURN VECTOR OF THE SIZE 30 x 30 WHICH HOLDS THE DERIVATIVES
!                     OF rhored WITH RESPECT TO ALL ni AND xj
!--------------------------------------------------------------------------




implicit none

    type(type_gl) :: gl


double precision, dimension(30, 30)::dndTr_dnidxj, dndrhor_dnidxj   ! Return matrices, maximum 30 components
double precision, dimension(30):: dTred_dxi, drhored_dxi 
double precision, dimension(30):: d2Tred_dxi2, d2rhored_dxi2
double precision, dimension(30, 30):: d2Tred_dxidxj, d2rhored_dxidxj
integer :: i, j, k
double precision:: xk

do i = 1, gl%ncomp
    do j = 1, gl%ncomp
        dndTr_dnidxj(i, j) = 0.D0
        dndrhor_dnidxj(i, j) = 0.D0
    end do
end do

! get 1st der. of Tr and rhor w.r.t. all xi
call dYr_dxi_old(gl,dTred_dxi, drhored_dxi)
! get 2nd der. of Tr and rhor w.r.t. all xi
call d2Yr_dxi2_old(gl,d2Tred_dxi2, d2rhored_dxi2)
! get the mixed der. of Tred and rhored w.r.t.  all xi and xj 
call d2Yr_dxidxj_old(gl,d2Tred_dxidxj, d2rhored_dxidxj)

! calculate d(n*(dTred/dni))/d(xj) and d(n*(drhored/dni))/d(xj)
! outer loop runs through all xj
do j = 1, gl%NCOMP
    ! inner loop runs through all ni
    do i = 1, gl%NCOMP
        if (j == i) then    !the mixed deriv. becomes the 2nd der w.r.t. xj
            dndTr_dnidxj(j, i) = d2Tred_dxi2(j) - dTred_dxi(j)
            dndrhor_dnidxj(j, i) = d2rhored_dxi2(j) - drhored_dxi(j)
        else
            dndTr_dnidxj(j, i) = d2Tred_dxidxj(j, i) - dTred_dxi(j)
            dndrhor_dnidxj(j, i) = d2rhored_dxidxj(j, i) - drhored_dxi(j)
        end if            
        do k = 1, gl%ncomp
            xk = gl%MOLFRACTIONS(k)
            if (j == k) then    !the mixed deriv. becomes the 2nd der w.r.t. xj
                dndTr_dnidxj(j, i) = dndTr_dnidxj(j, i) - xk*d2Tred_dxi2(k)
                dndrhor_dnidxj(j, i) = dndrhor_dnidxj(j, i) - xk*d2rhored_dxi2(k)
            else
                dndTr_dnidxj(j, i) = dndTr_dnidxj(j, i) - xk*d2Tred_dxidxj(j,k)
                dndrhor_dnidxj(j, i) = dndrhor_dnidxj(j, i) - xk*d2rhored_dxidxj(j,k)                
            end if
        end do
    end do
end do

end subroutine dndYr_dnidxj_old
!**************************************************************************

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    end submodule impl
