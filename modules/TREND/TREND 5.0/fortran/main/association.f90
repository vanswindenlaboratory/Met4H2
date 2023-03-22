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

    ! module for file association.f90
    module association_module
    !global use inclusion
    use module_all_types


    contains

	




!**************************************************************************
!           --------------------------------------------------
!           Routines for the calculation of the derivatives of 
!           the association grade ASSOGRADE with helpfunctions
!
!           S. Hielscher, Bochum, 04.2011 
!           --------------------------------------------------
!**************************************************************************

subroutine ASSOGRADEDERIVS (gl,T, D, GETASSOGRADEDER, SETASSOGRADEDER, Fluid)








implicit none

    type(type_gl) :: gl


double precision:: T, D  
integer, dimension(8):: GETASSOGRADEDER            ! array specifier to indicate, which derivative is needed 
double precision, dimension(8)::SETASSOGRADEDER     ! array with the computed values for the derivatives of the associationgrade
integer:: Fluid                                    ! fluidnumber for which fluid the assogradederivs have to be calculated 
double precision:: del, tau                                      ! reduced temperature and density

double precision:: ASSOGRADE, ASSOGRADED, ASSOGRADEDD            ! associationgrade and its
double precision:: ASSOGRADET, ASSOGRADETT, ASSOGRADEDT          ! derivatives with respect to
double precision:: ASSOGRADEDTT, ASSOGRADEDDD                    ! density and temperature

double precision:: ASTRENGTH                                     ! Constants to run helpfunctions
double precision:: AALFA, ABETA, AGAMMA                          ! only once
double precision:: AZETA, ATHETA, AMY, AXI                       ! 

double precision:: DXDDEL, DXDDDEL, DXDDDDEL, DXDDDGAMMA                        ! 
double precision:: DXDDDALFA, DXDDDX, DXDASSOSTRENGTH, DXDDDASSOSTRENGTH        ! 
double precision:: DXDTDX, DXDTDASSOSTRENGTH, DXDTDALFA, DXDTDBETA, DXDTDTHETA  ! Derivatives and sub-totals needed in the 3rd Model
double precision:: group12dd, group3dd, group4dd                                ! 
double precision:: group0dddn, group12dddn, group3dddn, group4dddn              ! 


! computes del and tau from reducing parameters
del = D / gl%rhored(fluid)      ! reduced density
tau = gl%tred(fluid) / T        ! inverse reduced temperature

    ASTRENGTH = 0.D0
    AALFA = 0.D0
    ABETA = 0.D0
    AGAMMA = 0.D0
    ATHETA = 0.D0
    AZETA = 0.D0
    AXI = 0.D0
    AMY = 0.D0

! ASTRENGTH is needed for every Derivative of the Associationgrade    
    ASTRENGTH = ASSOSTRENGTH(gl,Tau, Del, Fluid)

! Values of the Helpfunctions are transferd to constants to run each function only once
! not every Derivative of the Associationgrade needs each Helpfunction
! Helpfunctions can be found in ASSOMODEL Helpfunctions.f95
if ((GETASSOGRADEDER(2) >= 1).or.(GETASSOGRADEDER(3) >= 1).or.(GETASSOGRADEDER(5) >= 1).or.(GETASSOGRADEDER(6) >= 1).or.&
    &(GETASSOGRADEDER(7) >= 1).or.(GETASSOGRADEDER(8) >= 1)) then
    AALFA = ASSOALFA(gl,tau, del, Fluid)
endif
if ((GETASSOGRADEDER(4) >= 1).or.(GETASSOGRADEDER(5) >= 1).or.(GETASSOGRADEDER(6) >= 1).or.(GETASSOGRADEDER(7) >= 1)) then
    ABETA = ASSOBETA(gl,tau, del, Fluid)
endif
if ((GETASSOGRADEDER(3) >= 1).or.(GETASSOGRADEDER(8) >= 1)) then
    AGAMMA = ASSOGAMMA(gl,tau, del, Fluid)
endif
if ((GETASSOGRADEDER(6) >= 1).or.(GETASSOGRADEDER(7) >= 1)) then
    ATHETA = ASSOTHETA(gl,tau, del, Fluid)
endif
if ((GETASSOGRADEDER(5) >= 1).or.(GETASSOGRADEDER(7) >= 1)) then
    AZETA = ASSOZETA(gl,tau, del, Fluid)
endif
if (GETASSOGRADEDER(7) >= 1) then
    AXI = ASSOXI(gl,tau, del, Fluid)
endif
if (GETASSOGRADEDER(8) >= 1) then
    AMY = ASSOMY(gl,tau, del, Fluid)
endif

SETASSOGRADEDER = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)   ! initialize the return vector of the associationgrade

!***********************************************************
!   ASSOGRADE      CALCULATION OF THE ASSOCIATION GRADE [-]
!***********************************************************
    ASSOGRADE = 0.D0
    if ((gl%eos_coeff%nlna(fluid) == 1) .or. (gl%eos_coeff%nlna(fluid) == 2)) then
    ASSOGRADE = 2.D0/((1.D0 + 4.D0*ASTRENGTH*del)**0.5D0 + 1.D0)
    
    elseif (gl%eos_coeff%nlna(fluid) == 3) then
    ASSOGRADE = 2.D0/(1 - ASTRENGTH*del)
    end if
    
SETASSOGRADEDER(1) = ASSOGRADE

!***********************************************************
!   ASSOGRADED     CALCULATION OF THE FIRST DERIVATIVE OF THE ASSOCIATION GRADE
!                  WITH RESPECT TO DELTA [-]
!***********************************************************
if (GETASSOGRADEDER(2) >= 1) then
    
    ASSOGRADED = 0.D0
    if ((gl%eos_coeff%nlna(fluid) == 1) .or. (gl%eos_coeff%nlna(fluid) == 2)) then
    ASSOGRADED = -(ASSOGRADE**2)/(2.D0*ASTRENGTH*del*ASSOGRADE+1.D0)*(ASTRENGTH + del*AALFA)
    
    elseif (gl%eos_coeff%nlna(fluid) == 3) then
    ASSOGRADED = 0.5D0*ASSOGRADE**2*((ASTRENGTH + del*AALFA) - (ASTRENGTH**2*del + 3.D0*ASTRENGTH + (ASTRENGTH*&
    &del**2 + 3.D0*del)*AALFA)/(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**0.5D0)
    end if

SETASSOGRADEDER(2) = ASSOGRADED
end if
!***********************************************************
!   ASSOGRADEDD    CALCULATION OF THE SECOND DERIVATIVE OF THE ASSOCIATION GRADE
!                  WITH RESPECT TO DELTA [-]
!   Model 3 Derivs:
!       DXDDDEL     1st Derivative of ASSOGRADED(=XD) with respect to Del at constant ASSOGRADE(=X), ASSOSTRENGTH(gl,=DELTAAN), ASSOALFA therefore its NOT equal to ASSOGRADEDD
!       group12dd   sub-total
!       group3dd    sub-total
!       group4dd    sub-total
!***********************************************************
if (GETASSOGRADEDER(3) >= 1) then
    
    ASSOGRADEDD = 0.D0
    if ((gl%eos_coeff%nlna(fluid) == 1) .or. (gl%eos_coeff%nlna(fluid) == 2)) then
    ASSOGRADEDD = (ASSOGRADE**2*(2.D0*ASTRENGTH**2*ASSOGRADE - AALFA))/((2.D0*ASTRENGTH*del*ASSOGRADE + 1.D0)**2) - &
    &(ASTRENGTH + del*AALFA)*(2.D0*(ASTRENGTH*del*ASSOGRADE**2 + ASSOGRADE))/((2.D0*ASTRENGTH*del*ASSOGRADE + 1.D0)**2)*&
    &(-(ASTRENGTH*ASSOGRADE**2)/(2.D0*ASTRENGTH*del*ASSOGRADE+1.D0)) - (ASTRENGTH+del*AALFA)*(2.D0*(ASTRENGTH*del*&
    &ASSOGRADE**2 + ASSOGRADE))/((2.D0*ASTRENGTH*del*ASSOGRADE + 1.D0)**2)*(-(del*ASSOGRADE**2)/&
    &(2.D0*ASTRENGTH*del*ASSOGRADE+1.D0))*AALFA + ((ASSOGRADE**2*(2.D0*del**2*ASSOGRADE*AALFA - 1.D0))/&
    &((2.D0*ASTRENGTH*del*ASSOGRADE + 1.D0)**2))*AALFA-((del*ASSOGRADE**2)/(2.D0*ASTRENGTH*del*ASSOGRADE+1.D0))*AGAMMA
    
    elseif (gl%eos_coeff%nlna(fluid) == 3) then
        DXDDDEL = 0.D0
        DXDDDEL = 0.5D0*ASSOGRADE**2*(AALFA - (ASTRENGTH**2 + AALFA*(2.D0*ASTRENGTH*DEL + 3.D0))/(1.D0 + &
        &ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**0.5D0 + 0.5D0*(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*&
        &del)**(-1.5D0)*(ASTRENGTH**2*del + 3.D0*ASTRENGTH + (ASTRENGTH*del**2 + 3.D0*del)*AALFA)*(2.D0*&
        &ASTRENGTH**2*del + 6.D0*ASTRENGTH))
        
        group12dd = 0.D0
        group12dd = 0.5D0*(ASTRENGTH + AALFA*del)*ASSOGRADE**3*((ASTRENGTH + del*AALFA) - (ASTRENGTH**2*del + &
        &3.D0*ASTRENGTH + (ASTRENGTH*del**2 + 3.D0*del)*AALFA)/(1.D0 + ASTRENGTH**2*del**2 + 6.D0*&
        &ASTRENGTH*DEL)**0.5D0)*(1.D0 - (ASTRENGTH*del + 3.D0)/(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*&
        &del)**0.5D0)
        
        group3dd = 0.D0
        group3dd = 0.5D0*AALFA*ASSOGRADE**2*(1.D0 - (2.D0*ASTRENGTH*del + 3.D0 + AALFA*del**2)/(1.D0 + &
        &ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**0.5D0 + 0.5D0*(ASTRENGTH**2*del + 3.D0*ASTRENGTH + &
        &(ASTRENGTH*del**2 + 3.D0*del)*AALFA)*(2.D0*ASTRENGTH*del**2 + 6.D0*del)/(1.D0 + ASTRENGTH**2*del**2 + &
        &6.D0*ASTRENGTH*del)**0.5D0)
        
        group4dd = 0.D0
        group4dd = 0.5D0*ASSOGRADE**2*(del - (ASTRENGTH*del**2 + 3.D0*del)/(1.D0 + ASTRENGTH**2*del**2 + &
        &6.D0*ASTRENGTH*del)**0.5D0)*AGAMMA
        
    ASSOGRADEDD = DXDDDEL + group12dd + group3dd + group4dd
    end if
    
SETASSOGRADEDER(3) = ASSOGRADEDD
end if
!***********************************************************
!   ASSOGRADET     CALCULATION OF THE FIRST DERIVATIVE OF THE ASSOCIATION GRADE
!                  WITH RESPECT TO TAU [-]
!***********************************************************
if (GETASSOGRADEDER(4) >= 1) then
    
    ASSOGRADET = 0.D0
    if ((gl%eos_coeff%nlna(fluid) == 1) .or. (gl%eos_coeff%nlna(fluid) == 2)) then
    ASSOGRADET = -(del*ASSOGRADE**2)/(2.D0*ASTRENGTH*del*ASSOGRADE+1.D0)*ABETA
    
    elseif (gl%eos_coeff%nlna(fluid) == 3) then
    ASSOGRADET = 0.5D0*ASSOGRADE**2*del*(1.D0 - (ASTRENGTH*del + 3.D0)/(1.D0 + ASTRENGTH**2*del**2 + &
    &6.D0*ASTRENGTH*del)**0.5D0)*ABETA
    end if
    
SETASSOGRADEDER(4) = ASSOGRADET
end if
!***********************************************************
!   ASSOGRADETT    CALCULATION OF THE SECOND DERIVATIVE OF THE ASSOCIATION GRADE
!                  WITH RESPECT TO TAU [-]
!***********************************************************
if (GETASSOGRADEDER(5) >= 1) then
    ASSOGRADETT = 0.D0
    if ((gl%eos_coeff%nlna(fluid) == 1) .or. (gl%eos_coeff%nlna(fluid) == 2)) then
    ASSOGRADETT = -2.D0*ABETA*del*ASSOGRADE*(ASTRENGTH*del*ASSOGRADE+1.D0)/((1.D0 + 2.D0*ASTRENGTH*del*ASSOGRADE)**2)*&
    &(-(del*ASSOGRADE**2)/(2.D0*ASTRENGTH*del*ASSOGRADE+1.D0))*ABETA + (2.D0*del**2*ASSOGRADE**3)/&
    &((1.D0 + 2.D0*ASTRENGTH*del*ASSOGRADE)**2)*ABETA**2-(del*ASSOGRADE**2)/(2*ASTRENGTH*del*ASSOGRADE+1.D0)*&
    &AZETA
    
    elseif (gl%eos_coeff%nlna(fluid) == 3) then
    ASSOGRADETT = ABETA**2*0.5D0*ASSOGRADE**3*del**2*(1.D0 - (ASTRENGTH*del + 3.D0)/(1.D0 + ASTRENGTH**2*&
    &del**2 + 6.D0*ASTRENGTH*del)**0.5D0)**2 + 4.D0*del**2*ASSOGRADE**2*ABETA**2/(1.D0 + ASTRENGTH**2*&
    &del**2 + 6.D0*ASTRENGTH*del)**1.5D0 + 0.5D0*ASSOGRADE**2*del*(1.D0 - (ASTRENGTH*del + 3.D0)/(1.D0 + &
    &ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**0.5D0)*AZETA
    end if
    
SETASSOGRADEDER(5) = ASSOGRADETT
end if
!***********************************************************
!   ASSOGRADEDT    CALCULATION OF THE FIRST MIXED DERIVATIVE OF THE ASSOCIATION GRADE
!                  WITH RESPECT TO DELTA AND TAU [-]
!***********************************************************
if (GETASSOGRADEDER(6) >= 1) then
    ASSOGRADEDT = 0.D0
    if ((gl%eos_coeff%nlna(fluid) == 1) .or. (gl%eos_coeff%nlna(fluid) == 2)) then
    ASSOGRADEDT = -(ASTRENGTH + del*AALFA)*(2.D0*(ASTRENGTH*del*ASSOGRADE**2+ASSOGRADE)/((2.D0*ASTRENGTH*del*ASSOGRADE&
    &+ 1.D0)**2))*(-(del*ASSOGRADE**2)/(2.D0*ASTRENGTH*del*ASSOGRADE+1.D0))*ABETA + (ASSOGRADE**2*&
    &(2.D0*del**2*ASSOGRADE*AALFA - 1.D0))/((2.D0*ASTRENGTH*del*ASSOGRADE + 1.D0)**2)*ABETA - (del*ASSOGRADE**2)/&
    &(2.D0*ASTRENGTH*del*ASSOGRADE+1.D0)*ATHETA
    
    elseif (gl%eos_coeff%nlna(fluid) == 3) then
    ASSOGRADEDT = 0.5D0*ASSOGRADE**3*del*((ASTRENGTH + del*AALFA) - (ASTRENGTH**2*del + 3.D0*ASTRENGTH + (ASTRENGTH*&
    &del**2 + 3.D0*del)*AALFA/(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**0.5D0)*(1.D0 - (ASTRENGTH*del + &
    &3.D0)/(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**0.5D0))*ABETA + 0.5D0*ASSOGRADE**2*(1.D0 - &
    &(2.D0*ASTRENGTH*del + 3.D0 + AALFA*del**2)/(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**0.5D0 + &
    &0.5D0*(ASTRENGTH**2*del + 3.D0*ASTRENGTH + (ASTRENGTH*del**2 + 3.D0*del)*AALFA)*(2.D0*ASTRENGTH*del**2 + &
    &6.D0*del)/(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**1.5D0)*ABETA + 0.5D0*ASSOGRADE**2*del*(1.D0 - &
    &(ASTRENGTH*del + 3.D0)/(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**0.5D0)*ATHETA
    end if
    
SETASSOGRADEDER(6) = ASSOGRADEDT
end if

! Warning Model 3
! Third Derivatives are not verified!
!***********************************************************
!   ASSOGRADEDTT    CALCULATION OF THE SECOND DERIVATIVE OF THE ASSOCIATION GRADE
!                   WITH RESPECT TO DELTA [-]
!   used Derivatives:
!       DXDTDX              Derivative of ASSOGRADEDT with respect to Assograde(=X)
!       DXDTDASSOSTRENGTH   Derivative of ASSOGRADEDT with respect to ASSOSTRENGTH(gl,=Deltaan)
!       DXDTDALFA           Derivative of ASSOGRADEDT with respect to ASSOALFA
!       DXDTDBETA           Derivative of ASSOGRADEDT with respect to ASSOBETA
!       DXDTDTHETA          Derivative of ASSOGRADEDT with respect to ASSOTHETA
!       DXDASSOSTRENGTH     Derivative of ASSOGRADE with respect to ASSOSTRENGTH(gl,=Deltaan)
!
!***********************************************************
if (GETASSOGRADEDER(7) >= 1) then
    ASSOGRADEDTT = 0.D0
    if ((gl%eos_coeff%nlna(fluid) == 1) .or. (gl%eos_coeff%nlna(fluid) == 2) .or. (gl%eos_coeff%nlna(fluid) == 3)) then
        DXDTDX = 0.D0
        DXDTDX = (2.D0*ASSOGRADE*(-ABETA + ABETA*ASTRENGTH*DEL*ASSOGRADE + 4.D0*ABETA*&
        &ASTRENGTH**2*DEL**2*ASSOGRADE**2 + 6.D0*DEL**2*ASSOGRADE*ABETA*AALFA - 5.D0*&
        &DEL**2*ATHETA*ASTRENGTH*ASSOGRADE - 8.D0*DEL**3*ATHETA*ASTRENGTH**2*&
        &ASSOGRADE**2 - DEL*ATHETA + 12.D0*DEL**3*ASSOGRADE**2*ABETA*AALFA*&
        &ASTRENGTH + 2.D0*DEL**3*ASSOGRADE**3*ABETA*ASTRENGTH**3 - 4.D0*DEL**4*&
        &ATHETA*ASTRENGTH**3*ASSOGRADE**3 + 6.D0*DEL**4*ASSOGRADE**3*ABETA*&
        &AALFA*ASTRENGTH**2))/(1.D0 + 2.D0*ASTRENGTH*del*ASSOGRADE)**4

        DXDTDASSOSTRENGTH = 0.D0
        DXDTDASSOSTRENGTH = -(2.D0*ASSOGRADE**3*DEL*(-2.D0*ABETA*ASTRENGTH*DEL*ASSOGRADE + 2.D0*&
        &ABETA*ASTRENGTH**2*DEL**2*ASSOGRADE**2 - 3.D0*ABETA + 9.D0*DEL**2*ASSOGRADE*&
        &ABETA*AALFA + 12.D0*DEL**3*ASSOGRADE**2*ABETA*AALFA*ASTRENGTH - &
        &DEL*ATHETA - 4.D0*DEL**2*ATHETA*ASTRENGTH*ASSOGRADE - 4.D0*DEL**3*ATHETA*&
        &ASTRENGTH**2*ASSOGRADE**2))/(1.D0 + 2.D0*ASTRENGTH*DEL*ASSOGRADE)**4
    
        DXDTDALFA = 0.D0
        DXDTDALFA = 2.D0*ASSOGRADE**3*ABETA*DEL**2*(3.D0*ASTRENGTH*DEL*ASSOGRADE + 2.D0)/&
        &((1.D0 + 2.D0*ASTRENGTH*DEL*ASSOGRADE)**3)
    
        DXDTDBETA = 0.D0
        DXDTDBETA = (ASSOGRADE**2*(2.D0*ASTRENGTH**2*DEL**2*ASSOGRADE**2 + 6.D0*DEL**3*ASSOGRADE**2*&
        &AALFA*ASTRENGTH + 4.D0*DEL**2*ASSOGRADE*AALFA - 1.D0))/(1.D0 + 2.D0*ASTRENGTH*DEL*ASSOGRADE)**3
    
        DXDTDTHETA = 0.D0
        DXDTDTHETA = -DEL*ASSOGRADE**2/(1.D0 + 2.D0*ASTRENGTH*DEL*ASSOGRADE)
    
        DXDASSOSTRENGTH = 0.D0
        DXDASSOSTRENGTH = -DEL*ASSOGRADE**2/(2.D0*ASTRENGTH*DEL*ASSOGRADE + 1.D0)    
    
    ASSOGRADEDTT = DXDTDX*DXDASSOSTRENGTH*ABETA + DXDTDASSOSTRENGTH*ABETA + DXDTDALFA*ATHETA + DXDTDBETA*AZETA + &
    &DXDTDTHETA*AXI
    end if

SETASSOGRADEDER(7) = ASSOGRADEDTT
end if

! Warning Model 3
! Third Derivatives are not verified!
!***********************************************************
!   ASSOGRADEDDD    CALCULATION OF THE THIRD DERIVATIVE OF THE ASSOCIATION GRADE
!                   WITH RESPECT TO DELTA [-]
!   used Derivatives:
!       DXDDDDEL            First Derivative of ASSOGRADEDD with respect to Del at constant ASSOGRADE(=X), ASSOSTRENGTH(gl,=DELTAAN), ASSOALFA, ASSOGAMMA therefore its NOT equal to ASSOGRADEDDD
!       DXDDDX              First Derivative of ASSOGRADEDD with respect to the Associationgrade(=X)
!       DXDDEL              Derivative of ASSOGRADE(=X) with respect to Del(=Delta)
!       DXDASSOSTRENGTH     Derivative of ASSOGRADE(=X) with respect to ASSOSTRENGTH(gl,=Deltaan)
!       DXDDDASSOSTRENGTH   Derivative of ASSOGRADEDD with respect to ASSOSTRENGTH(gl,=Deltaan)
!       DXDDDALFA           Derivative of ASSOGRADEDD with respect to ASSOALFA
!       DXDDDGAMMA          Derivative of ASSOGRADEDD with respect to ASSOGAMMA
!
!       group12dddn         sub-total
!       group3dddn          sub-total
!       group4dddn          sub-total
!***********************************************************
 if (GETASSOGRADEDER(8) >= 1) then
    !First partial Derivative of the second derived Associationgrade X_dd with respect to Delta   
    if ((gl%eos_coeff%nlna(fluid) == 1) .or. (gl%eos_coeff%nlna(fluid) == 2)) then
    DXDDDDEL = 0.D0
    DXDDDDEL = (ASSOGRADE**2*(-AGAMMA + 12.D0*ASSOGRADE*ASTRENGTH*AALFA - 24.D0*&
    &ASSOGRADE**3*ASTRENGTH**4*DEL + 8.D0*ASSOGRADE**2*ASTRENGTH**2*AALFA*DEL&
    &-8.D0*ASSOGRADE**3*ASTRENGTH**3*AALFA*DEL**2 + 10.D0*AALFA**2*&
    &ASSOGRADE**2*DEL**2*ASTRENGTH - 18.D0*ASSOGRADE**2*ASTRENGTH**3 + 8.D0*&
    &AALFA**2*ASSOGRADE*DEL - 4.D0*AGAMMA*ASTRENGTH**2*DEL**2*ASSOGRADE**2 - 4.D0*&
    &AGAMMA*ASTRENGTH*DEL*ASSOGRADE))/((2.D0*ASTRENGTH*DEL*ASSOGRADE + 1.D0)**4)
    
    elseif (gl%eos_coeff%nlna(fluid) == 3) then
        group0dddn = -4.D0*ASSOGRADE**2*ASTRENGTH*(-2.D0*AALFA + 9.D0*ASTRENGTH**2 - 3.D0*AALFA*ASTRENGTH*del + &
        &AALFA*ASTRENGTH**2*del**2+3.D0*ASTRENGTH**3*del)*(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*&
        &del)**(-2.5D0)
        
        group12dddn = 0.5D0*AALFA*ASSOGRADE**3*(ASTRENGTH + del*AALFA - (ASTRENGTH**2*del + 3.D0*ASTRENGTH + &
        &AALFA*(ASTRENGTH*del**2 + 3.D0*del))*(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**(-0.5D0))*&
        &(1.D0 - (3.D0 + ASTRENGTH*del)*(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**(-0.5D0))+0.5D0*&
        &(ASTRENGTH + del*AALFA)*(ASSOGRADE**3)*((AALFA) - ((ASTRENGTH**2 + AALFA*(2.D0*ASTRENGTH*del + 3.D0))*&
        &(1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**(-0.5D0)) + 0.5D0*(ASTRENGTH**2*del + 3.D0*&
        &ASTRENGTH + AALFA*(ASTRENGTH*del**2 + 3.D0*del))*((1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*&
        &del)**(-1.5D0))*(2.D0*ASTRENGTH**2*del + 6.D0*ASTRENGTH))*(1.D0 - (3.D0 + ASTRENGTH*del)*(1.D0 + &
        &ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**(-0.5D0)) + 0.5D0*(ASTRENGTH + del*AALFA)*(ASSOGRADE**3)*&
        &(ASTRENGTH+del*AALFA - (ASTRENGTH**2*del + 3.D0*ASTRENGTH + AALFA*(ASTRENGTH*del**2 + 3.D0*del))*(1.D0 +&
        &ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*del)**(-0.5D0))*(-(ASTRENGTH*(1.D0 + ASTRENGTH**2*del**2 + &
        &6.D0*ASTRENGTH*del)**(-0.5D0))+0.5D0*(3.D0 + ASTRENGTH*del)*((1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*&
        &del)**(-1.5D0))*(2.D0*ASTRENGTH**2*del+6.D0*ASTRENGTH))
        
        group3dddn = -4.0D0*AALFA*ASSOGRADE**2*(-2.0D0*ASTRENGTH - 3.0D0*ASTRENGTH**2*del - 2.0D0*AALFA*del - 3.0D0*&
        &AALFA*ASTRENGTH*del**2 + ASTRENGTH**3*del**2 + AALFA*ASTRENGTH**2*del**3)*(1.0D0 + &
        &ASTRENGTH**2*del**2 + 6.0D0*ASTRENGTH*del)**(-2.5D0)
        
        group4dddn = 0.5D0*ASSOGRADE**2*(1.D0 - (ASTRENGTH*del + 3.D0)*(1.D0 + ASTRENGTH**2*del**2 + &
        &6.D0*ASTRENGTH*del)**(-0.5D0))*AGAMMA + 0.5D0*ASSOGRADE**2*del*(-(ASTRENGTH*(1.D0 + ASTRENGTH**2*del**2 + &
        &6.D0*ASTRENGTH*del)**(-0.5D0)) + 0.5D0*(ASTRENGTH*del + 3.D0)*((1.D0 + ASTRENGTH**2*del**2 + 6.D0*ASTRENGTH*&
        &del)**(-1.5D0))*(2*ASTRENGTH**2*del + 6.D0*ASTRENGTH))*AGAMMA
    
    DXDDDDEL = 0.D0
    DXDDDDEL = group0dddn+group12dddn+group3dddn+group4dddn
    end if
    
    DXDDEL = 0.D0
    DXDDEL = - ASTRENGTH*ASSOGRADE**2/(2.D0*ASTRENGTH*DEL*ASSOGRADE + 1.D0)

    DXDDDGAMMA = 0.D0
    DXDDDGAMMA = - DEL*ASSOGRADE**2/(2.D0*ASTRENGTH*DEL*ASSOGRADE + 1.D0)
        
    DXDDDALFA = 0.D0
    DXDDDALFA = (2.D0*ASSOGRADE**2*(-1.D0 + 2.D0*ASTRENGTH**2*DEL**2*ASSOGRADE**2 + 6.D0*&
    &DEL**3*ASSOGRADE**2*AALFA*ASTRENGTH + 4.D0*DEL**2*ASSOGRADE*AALFA))/&
    &(2.D0*ASTRENGTH*DEL*ASSOGRADE + 1.D0)**3
        
    DXDDDX = 0.D0
    DXDDDX = 1/((2.D0*ASTRENGTH*DEL*ASSOGRADE + 1.D0)**4)*(2.D0*ASSOGRADE*(-2.D0*AALFA + 2.D0*AALFA*&
    &ASTRENGTH*DEL*ASSOGRADE + 8.D0*ASSOGRADE**2*ASTRENGTH**2*DEL**2*AALFA + 12.D0&
    &*ASSOGRADE**2*DEL**3*AALFA**2*ASTRENGTH + 6.D0*ASTRENGTH**2*ASSOGRADE + &
    &12.D0*ASTRENGTH**3*ASSOGRADE**2*DEL + 6.D0*ASSOGRADE*DEL**2*&
    &AALFA**2 - DEL*AGAMMA - 8.D0*DEL**3*AGAMMA*ASTRENGTH**2*ASSOGRADE**2 - &
    &5.D0*DEL**2*AGAMMA*ASTRENGTH*ASSOGRADE + 6.D0*ASTRENGTH**4*ASSOGRADE**3*&
    &DEL**2 + 4.D0*ASSOGRADE**3*ASTRENGTH**3*DEL**3*AALFA + 6.D0*ASSOGRADE**3*DEL**4*&
    &AALFA**2*ASTRENGTH**2 - 4.D0*DEL**4*AGAMMA*ASTRENGTH**3*ASSOGRADE**3))
        
    DXDASSOSTRENGTH = 0.D0
    DXDASSOSTRENGTH = -DEL*ASSOGRADE**2/(2.D0*ASTRENGTH*DEL*ASSOGRADE + 1.D0)
      
    DXDDDASSOSTRENGTH = 0.D0
    DXDDDASSOSTRENGTH = -1.D0/((2.D0*ASTRENGTH*DEL*ASSOGRADE + 1.D0)**4)*(2.D0*ASSOGRADE**3*(-4.D0*&
    &ASTRENGTH - 4.D0*DEL**2*AALFA*ASTRENGTH*ASSOGRADE - 6.D0*DEL*AALFA - &
    &5.D0*ASTRENGTH**2*DEL*ASSOGRADE + 4.D0*ASTRENGTH**2*DEL**3*ASSOGRADE**2*&
    &AALFA + 12.D0*DEL**4*ASSOGRADE**2*AALFA**2*ASTRENGTH + 9.D0*DEL**3*&
    &ASSOGRADE*AALFA**2 - DEL**2*AGAMMA - 4.D0*DEL**4*AGAMMA*&
    &ASTRENGTH**2*ASSOGRADE**2 - 4.D0*DEL**3*AGAMMA*ASTRENGTH*ASSOGRADE))

ASSOGRADEDDD = 0.D0
            
ASSOGRADEDDD = DXDDDDEL + DXDDDX * DXDDEL + DXDDDX * DXDASSOSTRENGTH * AALFA + DXDDDASSOSTRENGTH *&
&AALFA + DXDDDALFA * AGAMMA + DXDDDGAMMA * AMY
   
SETASSOGRADEDER(8) = ASSOGRADEDDD
end if
end subroutine ASSOGRADEDERIVS


!**************************************************************************
!           --------------------------------------------------
!           Routines for the calculation of the functions of
!           the association grade (and if available its derivatives)
!
!           S. Hielscher, Bochum, 04.2011 
!           --------------------------------------------------
!           every function needs reduced Density and/or Temperature
!**************************************************************************
double precision function ASSODISTFUNCT(gl,del,Fluid)




implicit none

    type(type_gl) :: gl


    double precision:: del,ASSODENRED
    integer:: Fluid
    
    ASSODENRED = 0.d0
    ASSODISTFUNCT = 0.d0
    
    ASSODENRED = gl%assovolred(Fluid) * del
    !assovolred from module_asso
    
    ASSODISTFUNCT = 0.5D0*((2 - ASSODENRED)/((1-ASSODENRED)**3))
    

end function

double precision function ASSODISTFUNCTD_ETA(gl,del,Fluid)   !ASSODISTFUNCT Derivative _with respect to ETA


    

implicit none

    type(type_gl) :: gl

    
    double precision:: del,ASSODENRED
    integer:: Fluid
    
    ASSODENRED = 0.d0
    ASSODISTFUNCTD_ETA = 0.d0
    
    ASSODENRED = gl%assovolred(Fluid) * del
    !assovolred from module_asso
    
    ASSODISTFUNCTD_ETA = -0.5D0*(1/((1-ASSODENRED)**3))+ 1.5D0*(2-ASSODENRED)/((1-ASSODENRED)**4)
    

end function

double precision function ASSODISTFUNCTDD_ETA(gl,del, Fluid)    !ASSODISTFUNCT 2nd Derivative _with respect to ETA




implicit none

    type(type_gl) :: gl


    double precision:: del,ASSODENRED
    integer:: Fluid
    
    ASSODENRED = 0.d0
    ASSODISTFUNCTDD_ETA = 0.d0
    
    ASSODENRED = gl%assovolred(Fluid) * del
    !assovolred from module_asso
    
    ASSODISTFUNCTDD_ETA = -3.D0/((1-ASSODENRED)**4)+(6*(2-ASSODENRED))/((1-ASSODENRED)**5)
    

end function

double precision function ASSODISTFUNCTDDD_ETA(gl,del, Fluid)    !ASSODISTFUNCT 3rd Derivative _with respect to ETA




implicit none

    type(type_gl) :: gl


    double precision:: del,ASSODENRED
    integer:: Fluid

    
    ASSODENRED = 0.d0
    ASSODISTFUNCTDDD_ETA = 0.d0
    
    ASSODENRED = gl%assovolred(Fluid) * del
    !assovolred from module_asso
    
    ASSODISTFUNCTDDD_ETA = -18.D0/((1-ASSODENRED)**5)+(30*(2-ASSODENRED))/((1-ASSODENRED)**6)
    

end function

double precision function ASSOSTRENGTH(gl,tau, del, Fluid)




implicit none

    type(type_gl) :: gl


    double precision:: tau,del
    integer:: Fluid

    ASSOSTRENGTH = 0.D0
    
    ASSOSTRENGTH = ASSODISTFUNCT(gl,del,Fluid)*(dexp(gl%ASSOENERGY(Fluid)*TAU)-1)*gl%ASSOVOLINT(Fluid)
    
   
end function

double precision function ASSOALFA(gl,tau, del, Fluid)




implicit none

    type(type_gl) :: gl


    double precision:: tau,del
    integer:: Fluid
    
    ASSOALFA = 0.D0
    
    ASSOALFA = ASSODISTFUNCTD_ETA(gl,del,Fluid)*(dexp(gl%ASSOENERGY(Fluid)*tau)-1)*gl%ASSOVOLINT(Fluid)*gl%ASSOVOLRED(Fluid)
    
   
end function

double precision function ASSOBETA(gl,tau, del, Fluid)




implicit none

    type(type_gl) :: gl


    double precision:: tau,del
    integer:: Fluid
    
    ASSOBETA = 0.D0
    
    ASSOBETA = ASSODISTFUNCT(gl,del,Fluid)*gl%ASSOVOLINT(Fluid)*dexp(gl%ASSOENERGY(Fluid)*tau)*gl%ASSOENERGY(Fluid)
    
    
end function

double precision function ASSOGAMMA(gl,tau, del, Fluid)




implicit none

    type(type_gl) :: gl


    double precision:: tau,del
    integer:: Fluid
    
    ASSOGAMMA = 0.D0
    
    ASSOGAMMA = ASSODISTFUNCTDD_ETA(gl,del,Fluid)*(dexp(gl%ASSOENERGY(Fluid)*tau)-1)*gl%ASSOVOLINT(Fluid)*gl%ASSOVOLRED(Fluid)**2
    
  
end function

double precision function ASSOTHETA(gl,tau, del, Fluid)




implicit none

    type(type_gl) :: gl


    double precision:: tau,del
    integer:: Fluid
   
    ASSOTHETA = 0.D0
    
    ASSOTHETA = ASSODISTFUNCTD_ETA(gl,del,Fluid)*dexp(gl%ASSOENERGY(Fluid)*tau)*gl%ASSOENERGY(Fluid)*gl%ASSOVOLINT(Fluid)* &
    & gl%ASSOVOLRED(Fluid)
    
    
end function

double precision function ASSOZETA(gl,tau, del, Fluid)




implicit none

    type(type_gl) :: gl


    double precision:: tau,del
    integer:: Fluid
   
    ASSOZETA = 0.D0
    
    ASSOZETA = ASSODISTFUNCT(gl,del,Fluid)*gl%ASSOVOLINT(Fluid)*dexp(gl%ASSOENERGY(Fluid)*tau)*gl%ASSOENERGY(Fluid)**2
    
    
end function

double precision function ASSOMY(gl,tau, del, FLuid)




implicit none

    type(type_gl) :: gl


    double precision:: tau,del
    integer:: Fluid
   
    ASSOMY = 0.D0
    
    ASSOMY = ASSODISTFUNCTDDD_ETA(gl,del,Fluid)*(dexp(gl%ASSOENERGY(Fluid)*tau)-1)*gl%ASSOVOLINT(Fluid)*gl%ASSOVOLRED(Fluid)**3
    
    
end function

double precision function ASSOXI(gl,tau, del, Fluid)




implicit none

    type(type_gl) :: gl


    double precision:: tau,del
    integer:: Fluid
    
    ASSOXI = 0.D0
    
    ASSOXI = ASSODISTFUNCTD_ETA(gl,del,Fluid)*dexp(gl%ASSOENERGY(Fluid)*tau)*gl%ASSOENERGY(Fluid)**2*gl%ASSOVOLINT(Fluid)* &
    & gl%ASSOVOLRED(Fluid)
    
    
end function



!**************************************************************************
!           --------------------------------------------------
!           Routines for the calculation of the derivatives of 
!           the association part of the helmholz energy
!
!           S. Hielscher, Bochum, 04.2011 
!           --------------------------------------------------
!**************************************************************************
subroutine FNRASSODERIVS (gl,Temperature, Density, GETASSODER, SETFNRASSODER,Fluid)
!**************************************************************************
! SUBROUTINE FOR THE CALCULATION OF ALL DERIVATIVES OF THE ASSOCIATION PART
! OF THE HELMHOLTZ FREE ENERGY
! THE CALCULATION IS BASED ON THE FORMULATION OF THE HELMHOLTZ FREE ENERGY
! AS PUBLISHED BY:         
!                   Piazza, L; Span, R
!                   An equation of state for acetic acid including the
!                   association term of SAFT
!                   Fluid Phase Equilibria, 2011
!-------------------------------------------------------------------------
! INPUT PARAMETERS:     
! TEMPERATURE - T   K
! DENSITY     - D   KG/M^3
! GETASSODER  - AN ARRAY WITH 8 ENTRIES WITH VALUES EITHER "1" OR "0", 
!                INDICATING WHICH DERIVATIVES ARE NEEDED:
!                1. NORMALIZED RESIDUAL ASSOCIATION HELMHOLTZ ENERGY F AS A FUNCTION OF D AND T
!                2. 1ST DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY del*ASSOC
!                3. 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY del^2*ASSOC
!                4. 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU*ASSOC
!                5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2*ASSOC
!                6: 1ST MIXED DERIVATIVE OF F WITH RESPECT TO D AND T, MULTIPLIED BY TAU*del*ASSOC
!                7: 2ND MIXED DERIVATIVE OF F WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*del*ASSOC
!                8: 3RD DERIVATIVE OF F WITH RESPECT TO D, MULTIPLIED BY del^3*ASSOC
! OUTPUT PARAMETERS: 
! SETFNRASSODER - AN ARRAY WITH 8 ENTRIES WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
!                AS INDICATED IN "GETASSODER"
!-------------------------------------------------------------------------








implicit none

    type(type_gl) :: gl


double precision:: Temperature, Density  
integer, dimension(8):: GETASSODER             ! array specifier to indicate, which derivative of the residual association helmholtzenergy is needed 
integer, dimension(8):: GETASSOGRADEDER                    ! array specifier to indicate, which derivative of the associationgrade is needed
integer:: Fluid                                ! fluidnumber for which fluid fnrassoderivs have to be calculate for
double precision, dimension(8)::SETFNRASSODER   ! array with the computed values for the derivatives of the association helmholtzenergy
double precision, dimension(8):: SETASSOGRADEDER             ! array with the computed values for the derivatives of the associationgrade
double precision:: del, tau                                  ! reduced temperature and density

double precision:: M_ARASSO, MD_ARASSOD, MDD_ARASSODD        ! residual association helmholtz energy and its derivatives with
double precision:: MT_ARASSOT, MTT_ARASSOTT, MDT_ARASSODT    ! respect to density and temperature multiplied with the 
double precision:: MDTT_ARASSODTT, MDDD_ARASSODDD            ! muliplicative coefficient "m"


double precision:: ASSOGRADE, ASSOGRADED, ASSOGRADEDD        ! associationgrade and its
double precision:: ASSOGRADET, ASSOGRADETT, ASSOGRADEDT      ! derivatives with respect to
double precision:: ASSOGRADEDTT, ASSOGRADEDDD                ! density and temperature


! computes del and tau from reducing parameters
del = Density / gl%rhored(fluid)      ! reduced density
tau = gl%tred(fluid) / Temperature    ! inverse reduced temperature
if (gl%eos_coeff%nlna(fluid) /= 0) then


!WARNING SAFT Model 3 is not verified
!if (nna(fluid) == 3) then
!    if(nrwarn <= 20) then
!        nrwarn = nrwarn + 1
!        warn(nrwarn) = 10113
!    end if
!end if
!End WARNING

SETFNRASSODER = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)   ! initialize the return vector
GETASSOGRADEDER = (/0, 0, 0, 0, 0, 0, 0, 0/)   ! initialize the associationgrade get vector

if (GETASSODER(1) == 1) then
GETASSOGRADEDER = GETASSOGRADEDER + (/1,0,0,0,0,0,0,0/) ! Vector tells the subroutine which derivative is neccessary to calculate the property (0 or  > 1)
end if                                                  ! 0 means derivative is not needed,  > 1 means derivative is needed
if (GETASSODER(2) == 1) then
GETASSOGRADEDER = GETASSOGRADEDER + (/1,1,0,0,0,0,0,0/) ! for the first derivative of the associationgrade with respect to del the associationgrade itself is needed
end if
if (GETASSODER(3) == 1) then
GETASSOGRADEDER = GETASSOGRADEDER + (/1,1,1,0,0,0,0,0/) ! for the second derivative of the associationgrade with respect to del^2 the firt derivative and the associationgrade is needed
end if
if (GETASSODER(4) == 1) then
GETASSOGRADEDER = GETASSOGRADEDER + (/1,0,0,1,0,0,0,0/) ! and so on
end if
if (GETASSODER(5) == 1) then
GETASSOGRADEDER = GETASSOGRADEDER + (/1,0,0,1,1,0,0,0/)
end if
if (GETASSODER(6) == 1) then
GETASSOGRADEDER = GETASSOGRADEDER + (/1,1,0,1,0,1,0,0/) 
end if
if (GETASSODER(7) == 1) then
GETASSOGRADEDER = GETASSOGRADEDER + (/1,1,0,1,1,1,1,0/)
end if
if (GETASSODER(8) == 1) then
GETASSOGRADEDER = GETASSOGRADEDER + (/1,1,1,0,0,0,0,1/) 
end if

CALL ASSOGRADEDERIVS(gl,Temperature, Density, GETASSOGRADEDER, SETASSOGRADEDER,Fluid) ! subroutine for the calculation of the associationgrade in "ASSOGRADEDER.f95"

ASSOGRADE = SETASSOGRADEDER(1)                         ! Variables are assigned to Vectorrow
ASSOGRADED = SETASSOGRADEDER(2)
ASSOGRADEDD = SETASSOGRADEDER(3)
ASSOGRADET = SETASSOGRADEDER(4)
ASSOGRADETT = SETASSOGRADEDER(5)
ASSOGRADEDT = SETASSOGRADEDER(6)
ASSOGRADEDTT = SETASSOGRADEDER(7)
ASSOGRADEDDD = SETASSOGRADEDER(8)


!*****************************************************
!   M_ARASSO         RESIDUAL ASSOCIATION HELMHOLTZ FREE ENERGY F multiplied
!                    with the multiplicative coefficient ASSOC [-]
!*****************************************************
M_ARASSO = 0.D0     ! residual association helmholtzenergy is always needed

    if (gl%eos_coeff%nlna(fluid) == 1) then                                                 ! nna represents the Associationmodel (1,2B or 3), there are different calculations for each model
        M_ARASSO = (dlog(ASSOGRADE)- ASSOGRADE/2.D0 + 0.5D0)*gl%ASSOC(Fluid)         ! Calculates the first derivative of the asso-part 
    elseif (gl%eos_coeff%nlna(fluid) == 2) then
        M_ARASSO = (dlog(ASSOGRADE)- ASSOGRADE/2.D0 + 0.5D0)*2.D0*gl%ASSOC(Fluid)
    elseif (gl%eos_coeff%nlna(fluid) == 3) then
        M_ARASSO = (2.D0*(dlog(ASSOGRADE)- ASSOGRADE/2.D0) * (DLOG(2.D0*ASSOGRADE - 1.D0) - (2.D0*ASSOGRADE - 1.D0)/2.D0)&
        & + 1.5D0) * gl%ASSOC(Fluid)
    end if
SETFNRASSODER(1) = M_ARASSO      !writes the result into the outputvector

    !*****************************************************
    !   MD_ARASSOD       1ST DERIVATIVE OF F WITH RESPECT TO 
    !                    del MULTIPLIED WITH del*ASSOC [-]
    !*****************************************************
    if (GETASSODER(2) == 1) then    ! saving time if derivative is not needed
    MD_ARASSOD = 0.D0

        if (gl%eos_coeff%nlna(fluid) == 1) then
        MD_ARASSOD = ((1.D0/ASSOGRADE - 0.5D0)*ASSOGRADED)*del*gl%ASSOC(Fluid)    
        elseif (gl%eos_coeff%nlna(fluid) == 2) then                                              
        MD_ARASSOD = (2.D0*(1.D0/ASSOGRADE - 0.5D0)*ASSOGRADED)*del*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 3) then                                                         
        MD_ARASSOD = ((1.D0/ASSOGRADE - 1.D0 + 1.D0/(2.D0*ASSOGRADE - 1.D0))*ASSOGRADED)*2.D0*del*gl%ASSOC(Fluid)
        end if
    SETFNRASSODER(2) = MD_ARASSOD   !writes the result into the outputvector

    end if
    !*****************************************************
    !   MDD_ARASSODD    2ND DERIVATIVE OF F WITH RESPECT TO 
    !                   del MULTIPLIED WITH del^2*ASSOC [-]
    !*****************************************************
    if (GETASSODER(3) == 1) then
    MDD_ARASSODD = 0.D0

        if (gl%eos_coeff%nlna(fluid) == 1) then
        MDD_ARASSODD = ((1.D0/ASSOGRADE - 0.5D0)*ASSOGRADEDD - (ASSOGRADED**2)/(ASSOGRADE**2))*del**2*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 2) then
        MDD_ARASSODD = ((1.D0/ASSOGRADE - 0.5D0)*ASSOGRADEDD - (ASSOGRADED**2)/(ASSOGRADE**2))*2.D0*del**2*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 3) then 
        MDD_ARASSODD = ((1.D0/ASSOGRADE - 1.D0 + 1.D0/(2.D0*ASSOGRADE - 1.D0))*ASSOGRADEDD + (-1.D0/ASSOGRADE**2 - &
        &2.D0/((2.D0*ASSOGRADE - 1.D0)**2))*ASSOGRADED**2)*2.D0*del**2*gl%ASSOC(Fluid)
        end if
    SETFNRASSODER(3) = MDD_ARASSODD

    end if
    !*****************************************************
    !   MT_ARASSOT     1ST DERIVATIVE OF F WITH RESPECT TO 
    !                  tau MULTIPLIED WITH tau*ASSOC [-]
    !*****************************************************
    if (GETASSODER(4) == 1) then
    MT_ARASSOT = 0.D0

        if (gl%eos_coeff%nlna(fluid) == 1) then
        MT_ARASSOT = ((1.D0/ASSOGRADE - 0.5D0)*ASSOGRADET)*tau*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 2) then
        MT_ARASSOT = ((1.D0/ASSOGRADE - 0.5D0)*ASSOGRADET)*2.D0*tau*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 3) then
        MT_ARASSOT = ((1.D0/ASSOGRADE - 1.D0 + 1.D0/(2.D0*ASSOGRADE - 1.D0))*ASSOGRADET)*2.D0*tau*gl%ASSOC(Fluid)
        endif
    SETFNRASSODER(4) = MT_ARASSOT

    end if
    !*****************************************************
    !   MTT_ARASSOTT     2ND DERIVATIVE OF F WITH RESPECT TO 
    !                    tau MULTIPLIED WITH tau^2*ASSOC [-]
    !*****************************************************
    if (GETASSODER(5) == 1) then
    MTT_ARASSOTT = 0.D0

        if (gl%eos_coeff%nlna(fluid) == 1) then
        MTT_ARASSOTT = ((1.D0/ASSOGRADE - 0.5D0)*ASSOGRADETT - (ASSOGRADET**2)/(ASSOGRADE**2))*tau**2*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 2) then
        MTT_ARASSOTT = ((1.D0/ASSOGRADE - 0.5D0)*ASSOGRADETT - (ASSOGRADET**2)/(ASSOGRADE**2))*2.D0*tau**2*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 3) then
        MTT_ARASSOTT = ((1.D0/ASSOGRADE - 1.D0 + 1.D0/(2.D0*ASSOGRADE - 1.D0))*ASSOGRADETT + (-1.D0/ASSOGRADE**2 - &
        &2.D0/((2.D0*ASSOGRADE - 1.D0)**2))*ASSOGRADET**2)*2.D0*tau**2*gl%ASSOC(Fluid)
        endif
    SETFNRASSODER(5) = MTT_ARASSOTT

    end if
    !*****************************************************
    !   MDT_ARASSODT    1ST MIXED DERIVATIVE OF F WITH RESPECT TO 
    !                   del AND tau MULTIPLIED WITH del*tau*ASSOC [-]
    !*****************************************************
    if (GETASSODER(6) == 1) then
    MDT_ARASSODT = 0.D0

        if (gl%eos_coeff%nlna(fluid) == 1) then
        MDT_ARASSODT = ((-ASSOGRADET/(ASSOGRADE**2))*ASSOGRADED + ASSOGRADEDT*(1/ASSOGRADE - 0.5D0))*del*tau*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 2) then
        MDT_ARASSODT = ((-ASSOGRADET/(ASSOGRADE**2))*ASSOGRADED + ASSOGRADEDT*(1/ASSOGRADE - 0.5D0))*2.D0*del*tau*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 3) then
        MDT_ARASSODT = ((-1.D0/ASSOGRADE**2 - 2.D0/(2.D0*ASSOGRADE - 1.D0)**2)*ASSOGRADET*ASSOGRADED + &
        &(1.D0/ASSOGRADE - 1.D0 + 1.D0/(2.D0*ASSOGRADE - 1.D0))*ASSOGRADEDT)*2.D0*del*tau*gl%ASSOC(Fluid)
        endif
    SETFNRASSODER(6) = MDT_ARASSODT

    end if
    !*****************************************************
    !   MDTT_ARASSODTT   2ND MIXED DERIVATIVE OF F WITH RESPECT TO 
    !                    del AND tau MULTIPLIED WITH del*tau^2*ASSOC [-]
    !*****************************************************
    if (GETASSODER(7) == 1) then
    MDTT_ARASSODTT = 0.D0

        if (gl%eos_coeff%nlna(fluid) == 1) then
        MDTT_ARASSODTT = ((-ASSOGRADET/(ASSOGRADE**2))*ASSOGRADEDT + &
        & (-(ASSOGRADETT*ASSOGRADE**2 - 2*ASSOGRADE*ASSOGRADET**2)/(ASSOGRADE**4))*ASSOGRADED + &
        & ASSOGRADEDTT*(1.D0/ASSOGRADE - 0.5D0) + ASSOGRADEDT*(-ASSOGRADET/(ASSOGRADE**2)))*del*tau**2*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 2) then
        MDTT_ARASSODTT = ((-ASSOGRADET/(ASSOGRADE**2))*ASSOGRADEDT + &
        & (-(ASSOGRADETT*ASSOGRADE**2 - 2*ASSOGRADE*ASSOGRADET**2)/(ASSOGRADE**4))*ASSOGRADED + &
        & ASSOGRADEDTT*(1.D0/ASSOGRADE - 0.5D0) + ASSOGRADEDT*(-ASSOGRADET/(ASSOGRADE**2)))*del*tau**2*gl%ASSOC(Fluid)*2.D0
        elseif (gl%eos_coeff%nlna(fluid) == 3) then
        MDTT_ARASSODTT = ((1.D0/ASSOGRADE**3 + 4.D0/(2.D0*ASSOGRADE - 1.D0)**3)*2.D0*ASSOGRADET**2*ASSOGRADED - &
        & (1.D0/ASSOGRADE**2 + 2.D0/(2.D0*ASSOGRADE - 1.D0)**2)*(ASSOGRADETT*ASSOGRADED + 2.D0*ASSOGRADET*ASSOGRADEDT) + &
        & (1.D0/ASSOGRADE - 1.D0 + 1.D0/(2.D0*ASSOGRADE - 1.D0))*ASSOGRADEDTT)*del*tau**2*gl%ASSOC(Fluid)*2.D0
        endif
    SETFNRASSODER(7) = MDTT_ARASSODTT

    end if
    !*****************************************************
    !   MDDD_ARASSODDD   3RD DERIVATIVE OF F WITH RESPECT TO 
    !                    del MULTIPLIED WITH del^3*ASSOC [-]
    !*****************************************************
    if (GETASSODER(8) == 1) then
    MDDD_ARASSODDD = 0.D0

        if (gl%eos_coeff%nlna(fluid) == 1) then
        MDDD_ARASSODDD = ((1.D0/ASSOGRADE - 0.5D0)*ASSOGRADEDDD - ASSOGRADED/ASSOGRADE**2*ASSOGRADEDD - &
        & 2.D0*ASSOGRADED/ASSOGRADE**3*(-ASSOGRADED**2 + ASSOGRADEDD*ASSOGRADE))*del**3*gl%ASSOC(Fluid)
        elseif (gl%eos_coeff%nlna(fluid) == 2) then
        MDDD_ARASSODDD = ((1.D0/ASSOGRADE - 0.5D0)*ASSOGRADEDDD - ASSOGRADED/ASSOGRADE**2*ASSOGRADEDD - &
        & 2.D0*ASSOGRADED/ASSOGRADE**3*(-ASSOGRADED**2 + ASSOGRADEDD*ASSOGRADE))*2.D0*del**3*gl%ASSOC(Fluid)  
        elseif (gl%eos_coeff%nlna(fluid) == 3) then
        MDDD_ARASSODDD = ((-1.D0/ASSOGRADE**2 - 2.D0/(2.D0*ASSOGRADE - 1.D0)**2)*3.D0*ASSOGRADED*ASSOGRADEDD + &
        & (1.D0/ASSOGRADE - 1.D0 + 1.D0/(2*ASSOGRADE - 1.D0))*ASSOGRADEDDD + (2.D0/ASSOGRADE**3 + &
        & 8.D0/(2.D0*ASSOGRADE - 1.D0)**3)*ASSOGRADED**3)*2.D0*del**3*gl%ASSOC (Fluid)
        endif
    SETFNRASSODER(8) = MDDD_ARASSODDD

    end if

elseif (gl%eos_coeff%nlna(fluid) == 0) then
    SETFNRASSODER = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
end if
end subroutine FNRASSODERIVS

!
!!Dummy routines
!
!subroutine ARDERIVS(T, D, GETDERAR, DERAR)
!

!
!    implicit none
!    
!    !Declarations
!    double precision :: T, D
!    integer, dimension (nderivs) :: GETDERAR
!    double precision, dimension (nderivs) :: DERAR
!    
!    !Initializations
!    DERAR = 0.d0
!        
!end subroutine 


    end module association_module
