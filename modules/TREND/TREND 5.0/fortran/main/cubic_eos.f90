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
    submodule (cubic_eos_module) impl
    !global use inclusion
    use module_all_types
    use ancillary_equations_mix_module
    use setup_module
    use vle_derivs_module
    use reduced_parameters_calc_module

    
    
    contains

	






double precision module Function SRK(gl,T,P0_lc,iphase, nrsubst) 
!-----------------------------------------------------------
!Solves SRK_equation
!P=R*T/(V-b) - a/(V*V+b*V)
!= >  v^3+(-r*t/p)*v^2+(-b^2-r*t*b/p+a/p)v+(-a*b/p) = 0
!Input: P(Pa), T(K), and phase status (1-liquid; 2-vapor)
!Output: density kg/m3
!Algorithm for liquid density:
!1 SRK(gl,T,P)- > RhoL0; 2 ancillary equs (T) - >  Ps and RohLs
!3 SRK(gl,T,Ps)- > RhoLs'; 4 RhoL=RhoL0+(RohLs-RhoLs') 
!Hailong Li, 2009 Sep
!Denmark
!
!PS: if srk is negative, it means solving cubic equation is failed
!-----------------------------------------------------------



implicit none

    type(type_gl) :: gl

	real *8 :: t,p,P0_lc ! input T and P
    integer :: iphase !input iphase: 1-liquid; 2-vapor
    real *8 :: a,b,aa,tr,msrk,alfa_lc ! parameters of SRK
    !real *8 :: tc,pc,accen,wm ! properties: critical T and P, ascentric factor, and mole mass
    real *8 :: cc(3),roots(3)
    !real *8 :: req
    real *8 :: vsrk,srk1,vvpe,dv,ps
    integer :: nrsubst
	
	!CO2 properties (for test)
	!pc=7.39d6 !pa!
	!tc=304.1 !K
	!accen=0.225 
	!wm=44.01 !g/mol
	!req=8.314

	! read properties (done in the main program)
	!call read_file('co2.fld',i)
    p=P0_lc*1.d6    ! MPa- > Pa
    
	! parameters of SRK EOS
    	tr=t/gl%tc(nrsubst)
    
	aa=.42748d0*gl%req(nrsubst)*gl%req(nrsubst)*gl%tc(nrsubst)*gl%tc(nrsubst)/gl%pc(nrsubst)*1.d-6
	b=.08664d0*gl%req(nrsubst)*gl%tc(nrsubst)/gl%pc(nrsubst)*1.d-6
	msrk=0.48d0+1.574d0*gl%accen(nrsubst)-0.176d0*gl%accen(nrsubst)*gl%accen(nrsubst)
	alfa_lc=(1.d0+msrk*(1.d0-tr**0.5))**2
	a=aa*alfa_lc
    
	! parameter of cubic function
    !cc(1) = 1.d0
    cc(1) = -gl%req(nrsubst)*t/p
    cc(2) = -b**2-gl%req(nrsubst)*t*b/p+a/p
    cc(3) = -a*b/p
	
    ! Solves cc(1)*z^3 + cc(2)*z^2+cc(3)*z+cc(4) = 0
    call cubic_nt(gl,cc,roots) 
    
	if (roots(1) < -10.d0)then 
    	srk=-1.d0 !failed solving cubic equation
    else
      if (roots(1) < -0.d0)then! one real root
    	!vSRK=roots(2) !no matter liquid or vapor since only one real root
        VSRK=1.d0/(roots(2))!/wm(nrsubst) !kg/m3 - in mol/m^3, J.Gernert, Aug.2010
      else ! three real root
        if (iphase == 1)then
          !vSRK=roots(3) !liquid volume m3/mol
          VSRK=1.d0/(roots(3))!/wm(nrsubst) !mol/m3
        else
          !vSRK=roots(2) !vapor volume
          VSRK=1.d0/(roots(2))!/wm(nrsubst) !mol/m3
        endif
      endif
    endif
    srk=vsrk
    
!--------------------------------------------------------------------------------------0
	! modify the liquid density given by SRK based on the correction of saturated liquid density
	   
	if (iphase == 1)then 
		ps=vp_eq(gl,t, nrsubst)*1.d6 ! calculate saturated pressure with vapor pressure equation

		! calculate saturated liquid density with SRK---------------------1
    	!cc(1) = 1.d0
    	cc(1) = -gl%req(nrsubst)*t/ps
    	cc(2) = -b**2-gl%req(nrsubst)*t*b/ps+a/ps
    	cc(3) = -a*b/ps
    
    	call cubic_nt(gl,cc,roots) 
    
		if (roots(1) < 0)then ! one real root
    		!SRK1=roots(2) !no matter liquid or vapor since only one real root
        	SRK1=1.d0/(roots(2))!/wm(nrsubst) !mol/m3
    	else ! three real root
       	!SRK1=roots(3) !liquid volume m3/mol
          	SRK1=1.d0/(roots(3))!/wm(nrsubst) !mol/m3
    	endif
    
   	!------------------------------------------------------------------1
   
		vvpe=dl_eq(gl,t,nrsubst) ! calculate saturated liquid density with vapor pressure equation mol/l
        !vvpe=vvpe*wm
		dv=vvpe-srk1 ! calculate the difference 
  		srk=vsrk+dv !modify the corresponding liquid density given by SRK
	endif
!-------------------------------------------------------------------------------------------0
!	srk=1.d0/srk*wm !dm3/mol
End Function


module subroutine init_SRK(gl,Temp, nrsubst)







implicit none

    type(type_gl) :: gl


double precision:: temp
integer :: i, j, nrsubst
double precision :: tau

!Variables for PSRK
double precision :: gE_unifc
integer, dimension(nderivs):: GETDER                !array specifier to indicate, which derivative is needed 
!double precision, dimension(nderivs) :: gl%ge%gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
!double precision, dimension(nderivs) :: gl%ge%gE_R        !Residual part of gE and derivatives with respect to tau and delta
!double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
!double precision, dimension(nderivs, 30) :: gl%ge%ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau
integer:: C_or_R                !Specify which part of UNIFAC is needed (0:both, 1:only combinatorial, 2: only residual)
integer:: errval

!It is very important to initialize these variables here (set to 0)
gl%a_SRK = 0.D0
gl%b_SRK = 0.D0

!If whole mixture properties shall be calculated -->  nrsubst = 0
if (nrsubst == 0) then  
    Do i = 1, gl%ncomp        
        gl%ac_SRK(i) = 0.42747D0*gl%R_SRK*gl%R_SRK*gl%tc(i)*gl%tc(i)/(gl%pc(i)*1.D6/gl%factorpress)
        !TEST_AJ: Function for ac used by Stradi et al (2001). Used temporalily to check critical points of mixtures. Andreas May 2016
        !ac_SRK(i) = 0.42748D0*R_SRK*R_SRK*tc(i)*tc(i)/(pc(i)*1.D6/factorpress)
        gl%m_SRK(i) = 0.48D0+1.574D0*gl%accen(i)-0.176D0*gl%accen(i)*gl%accen(i)
    End do

    do i = 1, gl%ncomp
        !New version: Mathias and Copeman (Fluid Phase Equilib., 13, 91-108, 1983) temperature dependence for the parameter ai
        if ((dabs(gl%Cji_cubic(1,i)) < 1.D-12) .and. (dabs(gl%Cji_cubic(2,i)) < 1.D-12) .and. (dabs(gl%Cji_cubic(3,i)) < 1.D-12)) then   !If no Mathias Copeman parameters are given, then use classical formulation by setting C1i = mi
            gl%Cji_cubic(1,i) = gl%m_SRK(i)        
        end if
        gl%ai_SRK(i) = gl%ac_SRK(i) * (1.D0 + gl%Cji_cubic(1,i)*(1.D0-(Temp/gl%tc(i))**0.5D0) & 
                                    & + gl%Cji_cubic(2,i)*(1.D0-(Temp/gl%tc(i))**0.5D0)**2 & 
                                    & + gl%Cji_cubic(3,i)*(1.D0-(Temp/gl%tc(i))**0.5D0)**3 )**2
        !Old "classical" calculation of the parameter ai
        !ai_SRK(i) = ac_SRK(i) *(1.D0+m_SRK(i)*(1.D0-(Temp/tc(i))**0.5D0))**2
        gl%bi_SRK(i) = 0.08664D0*gl%R_SRK* gl%tc(i) / (gl%pc(i)*1.D6/gl%factorpress)
    end do
    
    if ((gl%mix_type == 2) .or. (gl%mix_type == 21)) then !quadratic mixing rules for parameter "a" and either linear or quadratic mixing rules for "b"  
        do j = 1, gl%ncomp
            do i = 1, gl%ncomp
                gl%aij_SRK(i,j) = (gl%ai_SRK(i) * gl%ai_SRK(j))**0.5D0 * (1.D0 - gl%kij_SRK(i,j)) 
                gl%a_SRK = gl%a_SRK + gl%molfractions(i) * gl%molfractions(j) * gl%aij_SRK(i,j)
            end do
            if (gl%mix_type == 2) then
                gl%b_SRK = gl%b_SRK +  gl%molfractions(j) * gl%bi_SRK(j)
            end if
        end do
 
        !Andreas, February 2016. Quadratic mixing rules for b
        if (gl%mix_type == 21) then
            do j = 1, gl%ncomp
                do i = 1, gl%ncomp
                    gl%bij_SRK(i,j) = (gl%bi_SRK(i) * gl%bi_SRK(j))**0.5D0 * (1.D0 - gl%lij_SRK(i,j)) 
                    gl%b_SRK = gl%b_SRK + gl%molfractions(i) * gl%molfractions(j) * gl%bij_SRK(i,j)
                end do
            end do    
        end if
    end if
    
    !Mixing rules for PSRK added, Andreas Jäger, January 2017
    if (gl%mix_type == 22) then
        !Linear mixing rules for "b"
        do i = 1,gl%ncomp
            gl%b_SRK = gl%b_SRK +  gl%molfractions(i) * gl%bi_SRK(i)   
        end do
        
        !PSRK mixing rules for "a"
        !Get the excess Gibbs energy from the UNIFAC model
        GETDER = 0
        GETDER(1) = 1   !Only the Gibbs energy itself is needed, no derivatives with respect to tau or delta
        C_or_R = 0      !Both, the combinatorial and the residual part, are needed
        !call gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER, gl%ge%gE_C, gl%ge%gE_R, gl%ge%ln_gamma_C, gl%ge%ln_gamma_R, C_or_R, errval)
        !Andreas Jäger, February 2019. Only call this routine for mixtures (gl%ncomp > 1). For pure components set gE to 0.
        If (gl%ncomp > 1) then
            call gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER, C_or_R, errval)
        else
            gl%ge%gE_C(1) = 0.D0 
            gl%ge%gE_R(1) = 0.D0   
        end if
        gE_unifc = gl%ge%gE_C(1) + gl%ge%gE_R(1)
        gl%a_SRK = gl%b_SRK * gE_unifc / A1_PSRK
        do i = 1,gl%ncomp
            gl%a_SRK = gl%a_SRK + gl%b_SRK * gl%molfractions(i) * gl%ai_SRK(i) / gl%bi_SRK(i) + gl%b_SRK * gl%R_SRK * Temp / A1_PSRK * gl%molfractions(i) * dlog(gl%b_SRK / gl%bi_SRK(i))         
        end do    
    end if
    
    !Important, Tred_SRK must be set to 1 here, when full SRK model is used
    !Andreas Jäger, Dec. 2013
    gl%Tred_SRK = 1.D0
else    
    !If a pure component in a mixture with the SRK is calculated, only the parameters of this pure component need to be calculated
    
    !-----------------------------------------------------------------------------------------------
    !Andreas Nov 2013. Implemented SRK in Helmholtz mixtures
    !Belege Tred_SRK
    if ((gl%Mix_Type == 1) .or. (gl%Mix_Type == 12)  .or. (gl%mix_type == 13)) then    !SRK als Helmholtz im Gemisch, GERG-Mischungsregeln werden verwendet
        gl%Tred_SRK = gl%TredMix
    else    !SRK Mischungsregeln werden verwendet
        !Tred_SRK = Tc(1)    !This was wrong, nrsubst should be standing here. Important if a pure component is to be calculated in a mixture. Corrected July 2016, Andreas Jäger
        gl%Tred_SRK = gl%Tc(nrsubst)
    end if

    !Berechne Tau
    tau = gl%Tred_SRK / Temp
    !-----------------------------------------------------------------------------------------------
       
    gl%ac_SRK(nrsubst) = 0.42747D0*gl%R_SRK*gl%R_SRK*gl%tc(nrsubst)*gl%tc(nrsubst)/(gl%pc(nrsubst)*1.D6/gl%factorpress)
    !TEST_AJ: Function for ac used by Stradi et al (2001). Used temporalily to check critical points of mixtures. Andreas May 2016
    !ac_SRK(nrsubst) = 0.42748D0*R_SRK*R_SRK*tc(nrsubst)*tc(nrsubst)/(pc(nrsubst)*1.D6/factorpress)
    gl%m_SRK(nrsubst) = 0.48D0+1.574D0*gl%accen(nrsubst)-0.176D0*gl%accen(nrsubst)*gl%accen(nrsubst)
    !New version: Mathias and Copeman (Fluid Phase Equilib., 13, 91-108, 1983) temperature dependence for the parameter ai
    if ((dabs(gl%Cji_cubic(1,nrsubst)) < 1.D-12) .and. (dabs(gl%Cji_cubic(2,nrsubst)) < 1.D-12) .and. (dabs(gl%Cji_cubic(3,nrsubst)) < 1.D-12)) then   !If no Mathias Copeman parameters are given, then use classical formulation by setting C1i = mi
        gl%Cji_cubic(1,nrsubst) = gl%m_SRK(nrsubst)        
    end if
    gl%ai_SRK(nrsubst) = gl%ac_SRK(nrsubst) * (1.D0 + gl%Cji_cubic(1,nrsubst)*(1.D0-(1.d0/tau)**0.5D0) & 
                                & + gl%Cji_cubic(2,nrsubst)*(1.D0-(1.d0/tau)**0.5D0)**2 & 
                                & + gl%Cji_cubic(3,nrsubst)*(1.D0-(1.d0/tau)**0.5D0)**3 )**2
    !Old "classical" calculation of the parameter ai
    !!Andreas, November 2015, below line ist the same, no need to devide tc by tc
    !!ai_SRK(nrsubst) = ac_SRK(nrsubst) *(1.D0+m_SRK(nrsubst)*(1.D0-((tc(nrsubst)/tau)/tc(nrsubst))**0.5D0))**2
    !gl%ai_SRK(nrsubst) = gl%ac_SRK(nrsubst) *(1.D0+gl%m_SRK(nrsubst)*(1.D0-(1.d0/tau)**0.5D0))**2
    gl%bi_SRK(nrsubst) = 0.08664D0*gl%R_SRK* gl%tc(nrsubst) / (gl%pc(nrsubst)*1.D6/gl%factorpress)
    gl%a_SRK = gl%ai_SRK(nrsubst)
    gl%b_SRK = gl%bi_SRK(nrsubst)
End if

! Andreas March 2013
!Save the temperature the SRK was last initialized with
gl%Temp_init = Temp 

end subroutine init_SRK




!Test of SRK_EOS transformed to the residual part of the Helmholtz-Energy
!Andreas March 2012

double precision module function rho_SRK (gl,Temp, press_in, Iphase, nrsubst)
    




implicit none

    type(type_gl) :: gl

    
    double precision:: cc(3),roots(3)                 ! 3 parameters of cEOS, 3 roots
    double precision:: Temp, press, press_in
    integer:: Iphase, nrsubst
    
    call init_SRK(gl,Temp, nrsubst)
        
    roots = 0.D0
    press = press_in * 1.D6/gl%factorpress    !Mpa -- >  Pa
    
    cc(1) = -gl%R_SRK*Temp/press
    cc(2) = -gl%b_SRK**2-gl%R_SRK*temp*gl%b_SRK/press+gl%a_SRK/press
    cc(3) = -gl%a_SRK*gl%b_SRK/press
    
    call cubic_nt(gl,cc,roots)                 ! Solves z^3 + cc(1)*z^2+cc(2)*z+cc(3) = 0
    
    !--handling roots----------------------------------------------------------------
	if (roots(1) < -10.d0)then 
    	rho_SRK = 0.d0                       ! failed solving cubic equation
    	
    else
      if (roots(1) < -0.d0)then               ! one real root
    	                                    ! vSRK=roots(2) !no matter liquid or vapor since only one real root
        rho_SRK=1.d0/(roots(2))!/wmm)          !mol/m^3 J.G. ! kg/m3
        
      else                                  ! three real roots
        if (iphase == 1)then                  ! vSRK=roots(3) !liquid volume m3/mol
          rho_SRK=1.d0/(roots(3))!/wmm)          !mol/m^3 J.G. ! kg/m3
          
        else                                ! vSRK=roots(2) !vapor volume
          rho_SRK=1.d0/(roots(2))!/wmm)          !mol/m^3 J.G. ! kg/m3
        endif
        
      endif
      
    endif
    
    if (rho_SRK > 1.d0/gl%b_SRK) then
        rho_SRK = 0.d0              ! solution of root is larger than pole
    endif
    
end function rho_SRK



!subroutine MIXDERIVSFNR_SRK(gl,TEMP, DENS, GETDER, MIXDERFNR)
!!**************************************************************************
!!--------------------------------------------------------------------------------------------------
!! SUBROUTINE FOR THE CALCULATION OF DERIVATIVES OF THE RESIDUAL PART OF THE 
!! HELMHOLTZ FREE ENERGY FOR MIXTURES
!! THE ROUTINE CALLS THE ROUTINE 'FNRDERIVS' FOR THE SINGLE FLUIDS AND THE ROUTINE 'DEPFUNCFNR' FOR 
!! THE BINARY DEPARTURE FUNCTION AT THE REDUCING PARAMETERS FOR T_RED AND RHO_RED. IT RETURNS THE 
!! RESIDUAL MIXTURE HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES
!!--------------------------------------------------------------------------------------------------
!! INPUT PARAMETERS:     
!! TEMPERATURE - T   K
!! DENSITY     - D   mol/m^3
!! GETDER      - AN ARRAY WITH 8 ENTRIES WITH VALUES EITHER "1" OR "0", 
!!                INDICATING WHICH DERIVATIVES ARE NEEDED:
!!                1. NORMALIZED RESIDUAL MIXTURE HELMHOLTZ ENERGY F AS A FUNCTION OF DEL AND TAU AND X
!!                2. 1ST DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
!!                3. 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
!!                4. 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
!!                5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
!!                6: 1ST MIXED DERIVATIVE OF F WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
!!                7: 2ND MIXED DERIVATIVE OF F WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
!!                8: 3RD DERIVATIVE OF F WITH RESPECT TO D, MULTIPLIED BY DEL^3
!! 
!! OUTPUT PARAMETERS: 
!! MIXDERFNR   - AN ARRAY WITH 8 ENTRIES WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
!!                AS INDECATED IN "GETDER"
!!--------------------------------------------------------------------------------------------------


!

!implicit none
!
!    type(type_gl) :: gl
!
!
!double precision:: TEMP, DENS 
!integer, dimension(8):: GETDER            ! array specifier to indicate, which derivative is needed 
!double precision,  dimension(8)::MIXDERFNR  ! array with the computed values for the derivatives
!
!double precision:: delta, tau
!
!double precision:: daSRK_dtau, da_SRK_dtau ,dT_dtau
!
!MIXDERFNR = 0.d0
!
!delta = Dens / gl%rhored_SRK
!tau = gl%Tred_SRK / Temp
!
!
!if (GETDER(1) == 1) then
!    MIXDERFNR(1) = -dlog(1.D0-gl%b_SRK*gl%rhored_SRK*delta)-gl%a_SRK*tau/gl%R_SRK/gl%tred_SRK/gl%b_SRK*dlog(1.D0+gl%b_SRK*gl%rhored_SRK*delta)
!end if
!
!if (GETDER(2) == 1) then
!    MIXDERFNR(2) = (gl%b_SRK * gl%rhored_SRK / (1.D0 - gl%b_SRK * gl%rhored_SRK * delta) - &
!                 & gl%a_SRK * gl%rhored_SRK * tau / (gl%R_SRK * gl%tred_SRK*(1.D0+gl%b_SRK*gl%rhored_SRK*delta)))*delta
!end if
!if (GETDER(3) == 1) then
!    MIXDERFNR(3) = delta*delta*(((gl%b_SRK*gl%rhored_SRK)/(1-gl%b_SRK*gl%rhored_SRK*delta))**2 + &
!                 & gl%a_SRK * tau / gl%R_SRK / gl%tred_SRK * gl%b_SRK * gl%rhored_SRK**2 / (1.d0+gl%b_SRK*gl%rhored_SRK*delta)**2)
!end if
!if (GETDER(4) == 1) then
!    daSRK_dtau = da_SRK_dtau(gl,Temp, 0)
!    dT_dtau = - Temp / tau
!    MIXDERFNR(4) = -(1.D0 / gl%R_SRK / gl%b_SRK * dlog(1.D0 + gl%b_SRK*gl%rhored_SRK*delta)* &
!                & (Temp*daSRK_dtau-dT_dtau*gl%a_SRK) / Temp**2) * tau                 
!end if
!if (GETDER(5) == 1) then
!    !Not yet implemented
!    MIXDERFNR(5) = 0.D0                
!end if
!if (GETDER(6) == 1) then
!    daSRK_dtau = da_SRK_dtau(gl,Temp, 0)
!    dT_dtau = - Temp / tau
!    MIXDERFNR(6) = - (gl%rhored_SRK / gl%R_SRK / (1.D0 + gl%b_SRK*gl%rhored_SRK*delta) * &
!                & (Temp*daSRK_dtau-dT_dtau*gl%a_SRK) / Temp**2) * delta * tau            
!end if
!if (GETDER(7) == 1) then
!    !Not yet implemented
!    MIXDERFNR(7) = 0.D0                
!end if
!if (GETDER(8) == 1) then
!    !Not yet implemented
!    MIXDERFNR(8) = 0.D0                
!end if
!
!
!end subroutine MIXDERIVSFNR_SRK


!**************************************************************************
module subroutine da_SRK_dxi(gl,daSRKdxi)
!**************************************************************************
!This subroutine calculates the partial derivative of the mixed parameter a
!of the SRK equation of state with respect to the compositions xi
!
!OUTPUT: vector daSRKdxi, containing the derivatives da / dxi at const. del, tau, xk





implicit none

    type(type_gl) :: gl


double precision, dimension(30):: daSRKdxi
integer:: i, j

daSRKdxi = 0.D0

Do i = 1, gl%ncomp - 1
    Do j = 1, gl%ncomp - 1
        daSRKdxi(i) = daSRKdxi(i) + gl%molfractions(j)*(gl%aij_SRK(j,i) + gl%aij_SRK(i,j) - gl%aij_SRK(j,gl%ncomp) - gl%aij_SRK(gl%ncomp,j))
    End do
    daSRKdxi(i) = daSRKdxi(i) - 2.D0 * gl%molfractions(gl%ncomp)*gl%aij_SRK(gl%ncomp,gl%ncomp) + & 
                & gl%molfractions(gl%ncomp) * (gl%aij_SRK(i,gl%ncomp) + gl%aij_SRK(gl%ncomp,i))
End do

end subroutine da_SRK_dxi
!**************************************************************************


!OLD VERSION OF THE DERIVATIVE OF THE MIXED PARAMETER a WRT the inverse
!reduced temperature tau
!This version may be deleted if the new version (see below) proofs to work
!Andreas March 2013
!!**************************************************************************
!double precision function da_SRK_dtau(gl,Temp)
!!**************************************************************************
!!This subroutine calculates the partial derivative of the mixed parameter a
!!of the SRK equation of state with respect to the inversed reduces temperature tau
!!
!!OUTPUT: scalar daSRKdT, containing the derivative da / dtau at const. del, xk
!


!
!implicit none
!      
!double precision, dimension(30):: daiSRKdT      !Temperature derivative of the parameter a of fluid i
!double precision, dimension(30,30)::daijSRKdT   !Temperature derivative of the binary parameter aij
!integer:: i, j
!double precision:: tau, temp
!
!da_SRK_dtau = 0.D0
!daiSRKdT = 0.D0
!daijSRKdT = 0.D0
!
!tau = Tred_SRK / Temp
!
!Do i = 1, ncomp
!    daiSRKdT(i) = -ac_SRK(i) *(1.D0+m_SRK(i)*(1.D0-(Temp/tc(i))**0.5))*m_SRK(i)/tc(i)/(Temp/tc(i))**0.5
!End do
!
!Do i = 1, ncomp
!    Do j = 1, ncomp
!        daijSRKdT(i,j) = (1.D0 - kij_SRK(i,j))/2.D0/(ai_SRK(i)*ai_SRK(j))**0.5*(ai_SRK(i)*daiSRKdT(j)+ai_SRK(j)*daiSRKdT(i))
!        da_SRK_dtau = da_SRK_dtau + molfractions(i)*molfractions(j)*daijSRKdT(i,j)
!    End do
!End do
!
!da_SRK_dtau = -da_SRK_dtau * Temp / tau
!
!end function da_SRK_dtau
!!**************************************************************************



!**************************************************************************
double precision module function da_SRK_dtau(gl,Temp, nrsubst)
!**************************************************************************
!This subroutine calculates the partial derivative of the mixed parameter a
!of the SRK equation of state with respect to the inversed reduces temperature tau
!
!OUTPUT: scalar daSRKdT, containing the derivative da / dtau at const. del, xk






implicit none

    type(type_gl) :: gl

      
double precision, dimension(30):: daiSRKdT      !Temperature derivative of the parameter a of fluid i
double precision, dimension(30,30)::daijSRKdT   !Temperature derivative of the binary parameter aij
double precision :: da_SRK_dT
integer:: i, j
double precision:: tau, temp

integer:: nrsubst

da_SRK_dtau = 0.D0
daiSRKdT = 0.D0
daijSRKdT = 0.D0
da_SRK_dT = 0.d0


!Belege Tred_SRK, Added mix_type 12 and 13 for new models, Andreas Jäger, January 2018
if ((gl%Mix_Type == 1) .or. (gl%Mix_Type == 12) .or. (gl%Mix_Type == 13)) then    !SRK als Helmholtz im Gemisch, GERG-Mischungsregeln werden verwendet
    gl%Tred_SRK = gl%TredMix
else    !SRK Mischungsregeln werden verwendet
    if (gl%ncomp == 1) then
        !Tred_SRK = Tc(1)    !Andreas Jäger, July 2016
        gl%Tred_SRK = gl%Tc(nrsubst)
    else
        gl%Tred_SRK = 1.d0
    end if
end if

!Berechne Tau
tau = gl%Tred_SRK / Temp
    
!!Reinstoffe    
!if (nrsubst /= 0) then    
!    da_SRK_dT = -ac_SRK(nrsubst) * (1.d0 + m_SRK(nrsubst) * (1.d0 - (Temp / tc(nrsubst)) ** 0.5d0)) &
!                & * m_SRK(nrsubst) / tc(nrsubst) / (Temp / tc(nrsubst)) ** 0.5d0
!    da_SRK_dtau = -da_SRK_dT * Temp / tau
!    
!!Gemische
!else
!    do i = 1, ncomp
!        daiSRKdT(i) = -ac_SRK(i) * (1.d0 + m_SRK(i) * (1.d0 - (Temp / tc(i)) ** 0.5d0)) * m_SRK(i) / tc(i) / (Temp / tc(i)) ** 0.5d0
!    end do
!    do i = 1, ncomp
!        do j = 1, ncomp
!            daijSRKdT(i, j) = (1.d0 - kij_SRK(i, j)) / 2.d0 / (ai_SRK(i) * ai_SRK(j)) **0.5d0 * &
!                            & (ai_SRK(i) * daiSRKdT(j) + ai_SRK(j) * daiSRKdT(i))
!            da_SRK_dT = da_SRK_dT + molfractions(i) * molfractions(j)  *daijSRKdT(i, j)
!        end do
!    end do
!    da_SRK_dtau = -da_SRK_dT * Temp / tau
!end if

!Andreas, November 2013
!Reinstoffe    
if (nrsubst /= 0) then    
    !da_SRK_dT = -ac_SRK(nrsubst) * (1.d0 + m_SRK(nrsubst) * (1.d0 - ((tc(nrsubst)/tau) / tc(nrsubst)) ** 0.5d0)) &
    !            & * m_SRK(nrsubst) / tc(nrsubst) / ((tc(nrsubst)/tau) / tc(nrsubst)) ** 0.5d0
    !da_SRK_dtau = -da_SRK_dT * (tc(nrsubst)/tau) / tau   
    !Alternative, Andreas Nov 2015
    da_SRK_dtau = gl%ac_SRK(nrsubst) * (1.d0 + gl%m_SRK(nrsubst) * (1.d0 - (1.D0/tau) ** 0.5d0)) &
                & * gl%m_SRK(nrsubst) * (1.D0/tau)**(-0.5D0) / tau**2
!Gemische
else    !SRK Mischungsregeln werden verwendet
    do i = 1, gl%ncomp
        !daiSRKdT(i) = -ac_SRK(i) * (1.d0 + m_SRK(i) * (1.d0 - ((Tred_SRK/tau) / tc(i)) ** 0.5d0)) * m_SRK(i) / tc(i) / ((Tred_SRK/tau) / tc(i)) ** 0.5d0
        !Andreas, December 2015
        daiSRKdT(i) = -gl%ac_SRK(i) * (1.d0 + gl%m_SRK(i) * (1.d0 - (Temp / gl%tc(i)) ** 0.5d0)) * gl%m_SRK(i) / gl%tc(i) / (Temp / gl%tc(i)) ** 0.5d0
    end do
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            daijSRKdT(i, j) = (1.d0 - gl%kij_SRK(i, j)) / 2.d0 / (gl%ai_SRK(i) * gl%ai_SRK(j)) **0.5d0 * &
                            & (gl%ai_SRK(i) * daiSRKdT(j) + gl%ai_SRK(j) * daiSRKdT(i))
            da_SRK_dT = da_SRK_dT + gl%molfractions(i) * gl%molfractions(j)  *daijSRKdT(i, j)
        end do
    end do
    da_SRK_dtau = -da_SRK_dT * (gl%Tred_SRK/tau) / tau
end if

end function da_SRK_dtau
!**************************************************************************
    
 
!**************************************************************************
double precision module function d2a_SRK_dtau2(gl,Temp, nrsubst)
!**************************************************************************
!Lars Hüttermann
!15.02.2013
!This function calculates the second partial derivate of the mixed parameter a of SRK with respect to the inversed reduces temperature tau.
!**************************************************************************








implicit none

    type(type_gl) :: gl



double precision :: Temp
double precision :: Tau
double precision :: da_dT
double precision :: d2a_dT2
double precision, dimension(30, 30) :: daij_dT
double precision, dimension(30, 30) :: d2aij_dT2
double precision, dimension(30) :: dai_dT
double precision, dimension(30) :: d2ai_dT2
integer :: i
integer :: j

integer:: nrsubst


Tau = 0.d0
da_dT = 0.d0
d2a_dT2 = 0.d0
daij_dT = 0.d0
d2aij_dT2 = 0.d0
dai_dT = 0.d0
d2ai_dT2 = 0.d0
i = 0
j = 0


!Belege Tred_SRK, Added mix_type 12 and 13 for new models, Andreas Jäger, January 2018
if ((gl%Mix_Type == 1) .or. (gl%Mix_Type == 12) .or. (gl%Mix_Type == 13)) then    !SRK als Helmholtz im Gemisch, GERG-Mischungsregeln werden verwendet
    gl%Tred_SRK = gl%TredMix
else    !SRK Mischungsregeln werden verwendet
    if (gl%ncomp == 1) then
        !Tred_SRK = Tc(1)    !Andreas Jäger, July 2016
        gl%Tred_SRK = gl%Tc(nrsubst)
    else
        gl%Tred_SRK = 1.d0
    end if
end if

!Berechne Tau
Tau = gl%Tred_SRK / Temp

!!Reinstoffe    
!if (nrsubst /= 0) then 
!    da_dT = -ac_SRK(nrsubst) * (1.d0 + m_SRK(nrsubst) * (1.d0 - (Temp / tc(nrsubst)) ** 0.5d0)) & 
!          & * m_SRK(nrsubst) / tc(nrsubst) / (Temp / tc(nrsubst)) ** 0.5d0
!    d2a_dT2 = -ac_SRK(nrsubst) * m_SRK(nrsubst) / tc(nrsubst) ** 0.5d0 * & 
!            & (-0.5d0 * Temp ** (-3.d0 / 2.d0) - 0.5d0 * m_SRK(nrsubst) * Temp ** (-3.d0 / 2.d0))
!    d2a_SRK_dtau2 = d2a_dT2 * Temp ** 2.d0 / Tau ** 2.d0 + 2.d0 * da_dT * Temp / Tau ** 2.d0
!
!!Gemische
!else
!    do i = 1, ncomp
!        dai_dT(i) = -ac_SRK(i) * (1.d0 + m_SRK(i) * (1.d0 - (Temp / tc(i)) ** 0.5d0)) * m_SRK(i) / tc(i) / (Temp / tc(i)) ** 0.5d0
!        d2ai_dT2(i) = -ac_SRK(i) * m_SRK(i) / tc(i) ** 0.5d0 * (-0.5d0 * Temp ** (-3.d0 / 2.d0) - 0.5d0 * m_SRK(i) * Temp ** (-3.d0 / 2.d0))
!    end do
!    do i = 1, ncomp
!        do j = 1, ncomp
!            daij_dT(i, j) = 0.5d0 * (1.d0 - kij_SRK(i, j)) * (ai_SRK(i) * ai_SRK(j)) ** (-0.5d0) * &
!                          & (ai_SRK(i) * dai_dT(j) + ai_SRK(j) * dai_dT(i))
!            d2aij_dT2(i, j) = 0.5d0 * (1.d0 - kij_SRK(i, j)) * &
!                            & (-0.5d0 * (ai_SRK(i) * ai_SRK(j)) ** (-3.d0 / 2.d0) * &
!                            & (ai_SRK(i) * dai_dT(j) + dai_dT(i) * ai_SRK(j)) ** 2.d0 + &
!                            & (ai_SRK(i) * ai_SRK(j)) ** (-0.5d0) * &
!                            & (d2ai_dT2(i) * ai_SRK(j) + d2ai_dT2(j) * ai_SRK(i) + &
!                            & 2.d0 * dai_dT(i) * dai_dT(j)))
!            da_dT = da_dT + molfractions(i) * molfractions(j) * daij_dT(i, j)
!            d2a_dT2 = d2a_dT2 + molfractions(i) * molfractions(j) * d2aij_dT2(i, j)
!        end do
!    end do
!    d2a_SRK_dtau2 = d2a_dT2 * Temp ** 2.d0 / Tau ** 2.d0 + 2.d0 * da_dT * Temp / Tau ** 2.d0
!end if

!Andreas, November 2013
!Reinstoffe    
if (nrsubst /= 0) then 
    da_dT = -gl%ac_SRK(nrsubst) * (1.d0 + gl%m_SRK(nrsubst) * (1.d0 - ((gl%tc(nrsubst)/tau) / gl%tc(nrsubst)) ** 0.5d0)) & 
          & * gl%m_SRK(nrsubst) / gl%tc(nrsubst) / ((gl%tc(nrsubst)/tau) / gl%tc(nrsubst)) ** 0.5d0
    d2a_dT2 = -gl%ac_SRK(nrsubst) * gl%m_SRK(nrsubst) / gl%tc(nrsubst) ** 0.5d0 * & 
            & (-0.5d0 * (gl%tc(nrsubst)/tau) ** (-3.d0 / 2.d0) - 0.5d0 * gl%m_SRK(nrsubst) * (gl%tc(nrsubst)/tau) ** (-3.d0 / 2.d0))
    d2a_SRK_dtau2 = d2a_dT2 * (gl%tc(nrsubst)/tau) ** 2 / Tau ** 2 + 2.d0 * da_dT * (gl%tc(nrsubst)/tau) / Tau ** 2

!Gemische
else    !SRK Mischungsregeln werden verwendet
    
    do i = 1, gl%ncomp
        dai_dT(i) = -gl%ac_SRK(i) * (1.d0 + gl%m_SRK(i) * (1.d0 - ((gl%Tred_SRK/tau) / gl%tc(i)) ** 0.5d0)) * gl%m_SRK(i) / gl%tc(i) / ((gl%Tred_SRK/tau) / gl%tc(i)) ** 0.5d0
        d2ai_dT2(i) = -gl%ac_SRK(i) * gl%m_SRK(i) / gl%tc(i) ** 0.5d0 * (-0.5d0 * (gl%Tred_SRK/tau) ** (-3.d0 / 2.d0) - 0.5d0 * gl%m_SRK(i) * (gl%Tred_SRK/tau) ** (-3.d0 / 2.d0))
    end do
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            daij_dT(i, j) = 0.5d0 * (1.d0 - gl%kij_SRK(i, j)) * (gl%ai_SRK(i) * gl%ai_SRK(j)) ** (-0.5d0) * &
                            & (gl%ai_SRK(i) * dai_dT(j) + gl%ai_SRK(j) * dai_dT(i))
            d2aij_dT2(i, j) = 0.5d0 * (1.d0 - gl%kij_SRK(i, j)) * &
                            & (-0.5d0 * (gl%ai_SRK(i) * gl%ai_SRK(j)) ** (-3.d0 / 2.d0) * &
                            & (gl%ai_SRK(i) * dai_dT(j) + dai_dT(i) * gl%ai_SRK(j)) ** 2 + &
                            & (gl%ai_SRK(i) * gl%ai_SRK(j)) ** (-0.5d0) * &
                            & (d2ai_dT2(i) * gl%ai_SRK(j) + d2ai_dT2(j) * gl%ai_SRK(i) + &
                            & 2.d0 * dai_dT(i) * dai_dT(j)))
            da_dT = da_dT + gl%molfractions(i) * gl%molfractions(j) * daij_dT(i, j)
            d2a_dT2 = d2a_dT2 + gl%molfractions(i) * gl%molfractions(j) * d2aij_dT2(i, j)
        end do
    end do
    d2a_SRK_dtau2 = d2a_dT2 * (gl%Tred_SRK/tau) ** 2 / Tau ** 2 + 2.d0 * da_dT * (gl%Tred_SRK/tau) / Tau ** 2

end if

end function d2a_SRK_dtau2


    
!**************************************************************************
double precision module function d3a_SRK_dtau3(gl,Temp, nrsubst)
!**************************************************************************
!Lars Hüttermann
!15.02.2013
!This function calculates the third partial derivate of the mixed parameter a of SRK with respect to the inversed reduces temperature tau.
!**************************************************************************






implicit none

    type(type_gl) :: gl



double precision :: Temp
double precision :: Tau
double precision :: da_dT
double precision :: d2a_dT2
double precision :: d3a_dT3
double precision, dimension(30, 30) :: daij_dT
double precision, dimension(30, 30) :: d2aij_dT2
double precision, dimension(30, 30) :: d3aij_dT3
double precision, dimension(30) :: dai_dT
double precision, dimension(30) :: d2ai_dT2
double precision, dimension(30) :: d3ai_dT3
integer :: i
integer :: j

integer:: nrsubst

Tau = 0.d0
da_dT = 0.d0
d2a_dT2 = 0.d0
d3a_dT3 = 0.d0
daij_dT = 0.d0
d2aij_dT2 = 0.d0
d3aij_dT3 = 0.d0
dai_dT = 0.d0
d2ai_dT2 = 0.d0
d3ai_dT3 = 0.d0
i = 0
j = 0


!Belege Tred_SRK, Added mix_type 12 and 13 for new models, Andreas Jäger, January 2018
if ((gl%Mix_Type == 1) .or. (gl%Mix_Type == 12) .or. (gl%Mix_Type == 13)) then    !SRK als Helmholtz im Gemisch, GERG-Mischungsregeln werden verwendet
    gl%Tred_SRK = gl%TredMix
else    !SRK Mischungsregeln werden verwendet
    if (gl%ncomp == 1) then
        !Tred_SRK = Tc(1)    !Andreas Jäger, July 2016
        gl%Tred_SRK = gl%Tc(nrsubst)
    else
        gl%Tred_SRK = 1.d0
    end if
end if

!Berechne Tau
Tau = gl%Tred_SRK / Temp

!!Reinstoffe    
!if (nrsubst /= 1) then 
!    da_dT = -ac_SRK(nrsubst) * (1.d0 + m_SRK(nrsubst) * (1.d0 - (Temp / tc(nrsubst)) ** 0.5d0)) &
!          & * m_SRK(nrsubst) / tc(nrsubst) / (Temp / tc(nrsubst)) ** 0.5d0
!    d2a_dT2 = -ac_SRK(nrsubst) * m_SRK(nrsubst) / tc(nrsubst) ** 0.5d0 * & 
!            & (-0.5d0 * Temp ** (-3.d0 / 2.d0) - 0.5d0 * m_SRK(nrsubst) * Temp ** (-3.d0 / 2.d0))
!    d3a_dT3 = -ac_SRK(nrsubst) * m_SRK(nrsubst) / tc(nrsubst) ** 0.5d0 * &
!            & (3.d0 / 4.d0 * Temp ** (-5.d0 / 2.d0) + 3.d0 / 4.d0 * m_SRK(nrsubst) * &
!            & Temp ** (-5.d0 / 2.d0))
!    d3a_SRK_dtau3 = -d3a_dT3 * Temp ** 3.d0 / Tau ** 3.d0 - 6.d0 * d2a_dT2 * Temp ** 2.d0 / Tau ** 3.d0 - &
!                  & 6.d0 * da_dT * Temp / Tau ** 3.d0
!
!!Gemische
!else
!    do i = 1, ncomp
!        dai_dT(i) = -ac_SRK(i) * (1.d0 + m_SRK(i) * (1.d0 - (Temp / tc(i)) ** 0.5d0)) * m_SRK(i) / tc(i) / (Temp / tc(i)) ** 0.5d0
!        d2ai_dT2(i) = -ac_SRK(i) * m_SRK(i) / tc(i) ** 0.5d0 * (-0.5d0 * Temp ** (-3.d0 / 2.d0) - 0.5d0 * m_SRK(i) * Temp ** (-3.d0 / 2.d0))
!        d3ai_dT3(i) = -ac_SRK(i) * m_SRK(i) / tc(i) ** 0.5d0 * (3.d0 / 4.d0 * Temp ** (-5.d0 / 2.d0) + 3.d0 / 4.d0 * m_SRK(i) * &
!                    & Temp ** (-5.d0 / 2.d0))
!    end do
!    do i = 1, ncomp
!        do j = 1, ncomp
!            daij_dT(i, j) = 0.5d0 * (1.d0 - kij_SRK(i, j)) * (ai_SRK(i) * ai_SRK(j)) ** (-0.5d0) * &
!                          & (ai_SRK(i) * dai_dT(j) + ai_SRK(j) * dai_dT(i))
!            d2aij_dT2(i, j) = 0.5d0 * (1.d0 - kij_SRK(i, j)) * &
!                            & (-0.5d0 * (ai_SRK(i) * ai_SRK(j)) ** (-3.d0 / 2.d0) * &
!                            & (ai_SRK(i) * dai_dT(j) + dai_dT(i) * ai_SRK(j)) ** 2.d0 + &
!                            & (ai_SRK(i) * ai_SRK(j)) ** (-0.5d0) * &
!                            & (d2ai_dT2(i) * ai_SRK(j) + d2ai_dT2(j) * ai_SRK(i) + &
!                            & 2.d0 * dai_dT(i) * dai_dT(j)))
!            d3aij_dT3(i, j) = 1.d0 / 2.d0 * (1.d0 - kij_SRK(i, j)) * &
!                            & (3.d0 / 4.d0 * (ai_SRK(i) * ai_SRK(j)) ** (-5.d0 / 2.d0) * &
!                            & (ai_SRK(i) * dai_dT(j) + dai_dT(i) * ai_SRK(j)) ** 3.d0 - &
!                            & (ai_SRK(i) * ai_SRK(j)) ** (-3.d0 / 2.d0) * &
!                            & (ai_SRK(i) * dai_dT(j) + dai_dT(i) * ai_SRK(j)) * &
!                            & (d2ai_dT2(i) * ai_SRK(j) + ai_SRK(i) * d2ai_dT2(j) + 2.d0 * dai_dT(i) * dai_dT(j)) + &
!                            & (ai_SRK(i) * ai_SRK(j)) ** (-1.d0 / 2.d0) * &
!                            & (d3ai_dT3(i) * ai_SRK(j) + ai_SRK(i) * d3ai_dT3(j) + &
!                            & 3.d0 * d2ai_dT2(i) * dai_dT(j) + 3.d0 * dai_dT(i) * d2ai_dT2(j)) - &
!                            & 1.d0 / 2.d0 * (ai_SRK(i) * dai_dT(j) + dai_dT(i) * ai_SRK(j)) * &
!                            & (ai_SRK(i) * ai_SRK(j)) ** (-3.d0 / 2.d0) * &
!                            & (d2ai_dT2(i) * ai_SRK(j) + ai_SRK(i) * d2ai_dT2(j) + 2.d0 * dai_dT(i) * dai_dT(j)))
!            da_dT = da_dT + molfractions(i) * molfractions(j) * daij_dT(i, j)
!            d2a_dT2 = d2a_dT2 + molfractions(i) * molfractions(j) * d2aij_dT2(i, j)
!            d3a_dT3 = d3a_dT3 + molfractions(i) * molfractions(j) * d3aij_dT3(i, j)
!        end do
!    end do
!    d3a_SRK_dtau3 = -d3a_dT3 * Temp ** 3.d0 / Tau ** 3.d0 - 6.d0 * d2a_dT2 * Temp ** 2.d0 / Tau ** 3.d0 - &
!                  & 6.d0 * da_dT * Temp / Tau ** 3.d0
!end if


!Reinstoffe    
if (nrsubst /= 0) then 
    da_dT = -gl%ac_SRK(nrsubst) * (1.d0 + gl%m_SRK(nrsubst) * (1.d0 - ((gl%tc(nrsubst)/tau) / gl%tc(nrsubst)) ** 0.5d0)) &
          & * gl%m_SRK(nrsubst) / gl%tc(nrsubst) / ((gl%tc(nrsubst)/tau) / gl%tc(nrsubst)) ** 0.5d0
    d2a_dT2 = -gl%ac_SRK(nrsubst) * gl%m_SRK(nrsubst) / gl%tc(nrsubst) ** 0.5d0 * & 
            & (-0.5d0 * (gl%tc(nrsubst)/tau) ** (-3.d0 / 2.d0) - 0.5d0 * gl%m_SRK(nrsubst) * (gl%tc(nrsubst)/tau) ** (-3.d0 / 2.d0))
    d3a_dT3 = -gl%ac_SRK(nrsubst) * gl%m_SRK(nrsubst) / gl%tc(nrsubst) ** 0.5d0 * &
            & (3.d0 / 4.d0 * (gl%tc(nrsubst)/tau) ** (-5.d0 / 2.d0) + 3.d0 / 4.d0 * gl%m_SRK(nrsubst) * &
            & (gl%tc(nrsubst)/tau) ** (-5.d0 / 2.d0))
    d3a_SRK_dtau3 = -d3a_dT3 * (gl%tc(nrsubst)/tau) ** 3 / Tau ** 3 - 6.d0 * d2a_dT2 * (gl%tc(nrsubst)/tau) ** 2 / Tau ** 3 - &
                  & 6.d0 * da_dT * (gl%tc(nrsubst)/tau) / Tau ** 3
!Gemische
else    !SRK als Helmholtz im Gemisch, GERG-Mischungsregeln werden verwendet
            
    do i = 1, gl%ncomp
        dai_dT(i) = -gl%ac_SRK(i) * (1.d0 + gl%m_SRK(i) * (1.d0 - ((gl%Tred_SRK/tau) / gl%tc(i)) ** 0.5d0)) * gl%m_SRK(i) / gl%tc(i) / ((gl%Tred_SRK/tau) / gl%tc(i)) ** 0.5d0
        d2ai_dT2(i) = -gl%ac_SRK(i) * gl%m_SRK(i) / gl%tc(i) ** 0.5d0 * (-0.5d0 * (gl%Tred_SRK/tau) ** (-3.d0 / 2.d0) - 0.5d0 * gl%m_SRK(i) * (gl%Tred_SRK/tau) ** (-3.d0 / 2.d0))
        d3ai_dT3(i) = -gl%ac_SRK(i) * gl%m_SRK(i) / gl%tc(i) ** 0.5d0 * (3.d0 / 4.d0 * (gl%Tred_SRK/tau) ** (-5.d0 / 2.d0) + 3.d0 / 4.d0 * gl%m_SRK(i) * &
                    & (gl%Tred_SRK/tau) ** (-5.d0 / 2.d0))
    end do
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            daij_dT(i, j) = 0.5d0 * (1.d0 - gl%kij_SRK(i, j)) * (gl%ai_SRK(i) * gl%ai_SRK(j)) ** (-0.5d0) * &
                            & (gl%ai_SRK(i) * dai_dT(j) + gl%ai_SRK(j) * dai_dT(i))
            d2aij_dT2(i, j) = 0.5d0 * (1.d0 - gl%kij_SRK(i, j)) * &
                            & (-0.5d0 * (gl%ai_SRK(i) * gl%ai_SRK(j)) ** (-3.d0 / 2.d0) * &
                            & (gl%ai_SRK(i) * dai_dT(j) + dai_dT(i) * gl%ai_SRK(j)) ** 2 + &
                            & (gl%ai_SRK(i) * gl%ai_SRK(j)) ** (-0.5d0) * &
                            & (d2ai_dT2(i) * gl%ai_SRK(j) + d2ai_dT2(j) * gl%ai_SRK(i) + &
                            & 2.d0 * dai_dT(i) * dai_dT(j)))
            d3aij_dT3(i, j) = 1.d0 / 2.d0 * (1.d0 - gl%kij_SRK(i, j)) * &
                            & (3.d0 / 4.d0 * (gl%ai_SRK(i) * gl%ai_SRK(j)) ** (-5.d0 / 2.d0) * &
                            & (gl%ai_SRK(i) * dai_dT(j) + dai_dT(i) * gl%ai_SRK(j)) ** 3 - &
                            & (gl%ai_SRK(i) * gl%ai_SRK(j)) ** (-3.d0 / 2.d0) * &
                            & (gl%ai_SRK(i) * dai_dT(j) + dai_dT(i) * gl%ai_SRK(j)) * &
                            & (d2ai_dT2(i) * gl%ai_SRK(j) + gl%ai_SRK(i) * d2ai_dT2(j) + 2.d0 * dai_dT(i) * dai_dT(j)) + &
                            & (gl%ai_SRK(i) * gl%ai_SRK(j)) ** (-1.d0 / 2.d0) * &
                            & (d3ai_dT3(i) * gl%ai_SRK(j) + gl%ai_SRK(i) * d3ai_dT3(j) + &
                            & 3.d0 * d2ai_dT2(i) * dai_dT(j) + 3.d0 * dai_dT(i) * d2ai_dT2(j)) - &
                            & 1.d0 / 2.d0 * (gl%ai_SRK(i) * dai_dT(j) + dai_dT(i) * gl%ai_SRK(j)) * &
                            & (gl%ai_SRK(i) * gl%ai_SRK(j)) ** (-3.d0 / 2.d0) * &
                            & (d2ai_dT2(i) * gl%ai_SRK(j) + gl%ai_SRK(i) * d2ai_dT2(j) + 2.d0 * dai_dT(i) * dai_dT(j)))
            da_dT = da_dT + gl%molfractions(i) * gl%molfractions(j) * daij_dT(i, j)
            d2a_dT2 = d2a_dT2 + gl%molfractions(i) * gl%molfractions(j) * d2aij_dT2(i, j)
            d3a_dT3 = d3a_dT3 + gl%molfractions(i) * gl%molfractions(j) * d3aij_dT3(i, j)
        end do
    end do
    d3a_SRK_dtau3 = -d3a_dT3 * (gl%Tred_SRK/tau) ** 3 / Tau ** 3 - 6.d0 * d2a_dT2 * (gl%Tred_SRK/tau) ** 2 / Tau ** 3 - &
                    & 6.d0 * da_dT * (gl%Tred_SRK/tau) / Tau ** 3        
   
end if


end function d3a_SRK_dtau3
!**************************************************************************


!**************************************************************************
module subroutine d2a_SRK_dxidtau(gl,Temp,d2aSRKdxidtau)
!**************************************************************************
!Second partial derivative of the mixed parameter a
!of the SRK equation of state with respect to the compositions xi and tau
!
!OUTPUT: vector d2aSRKdxidtau, containing the derivatives d2a / dxi dtau 





implicit none

    type(type_gl) :: gl


double precision, dimension(30):: d2aSRKdxidtau
double precision:: temp
double precision, dimension(30,30)::daijSRKdT   !Temperature derivative of the binary parameter aij
double precision, dimension(30):: daiSRKdT      !Temperature derivative of the parameter a of fluid i
integer:: i, j
double precision:: tau, dT_dtau

d2aSRKdxidtau = 0.D0

tau = gl%Tred_SRK / temp
dT_dtau = - Temp / tau

!Get the temperature derivative of the parameter a for all pure fluid equations
Do i = 1, gl%ncomp
    daiSRKdT(i) = -gl%ac_SRK(i) *(1.D0+gl%m_SRK(i)*(1.D0-(Temp/gl%tc(i))**0.5d0))*gl%m_SRK(i)/gl%tc(i)/(Temp/gl%tc(i))**0.5d0
End do

!Get the temperature derivative of the binary parameter aij
Do i = 1, gl%ncomp
    Do j = 1, gl%ncomp
         daijSRKdT(i,j) = (1.D0 - gl%kij_SRK(i,j))/2.D0/(gl%ai_SRK(i)*gl%ai_SRK(j))**0.5d0*(gl%ai_SRK(i)*daiSRKdT(j)+gl%ai_SRK(j)*daiSRKdT(i))
    End do
End do


Do i = 1, gl%ncomp - 1
    Do j = 1, gl%ncomp - 1        
        d2aSRKdxidtau(i) = d2aSRKdxidtau(i) + &
                         & gl%molfractions(j)*(daijSRKdT(j,i) + daijSRKdT(i,j) - daijSRKdT(j,gl%ncomp) - daijSRKdT(gl%ncomp,j))
    End do
    d2aSRKdxidtau(i) = d2aSRKdxidtau(i) - 2.D0 * gl%molfractions(gl%ncomp)*daijSRKdT(gl%ncomp,gl%ncomp) + & 
                     & gl%molfractions(gl%ncomp) * (daijSRKdT(i,gl%ncomp) + daijSRKdT(gl%ncomp,i))
    d2aSRKdxidtau(i) = d2aSRKdxidtau(i) * dT_dtau
End do

end subroutine d2a_SRK_dxidtau
!**************************************************************************

!**************************************************************************
module subroutine d2a_SRK_dxidxj(gl,d2aSRKdxidxj)
!**************************************************************************
!This subroutine calculates the second derivative of the mixed parameter a
!of the SRK equation of state with respect to the compositions xi and xj
!
!OUTPUT: matrix daSRKdxi, containing the derivatives da / dxi at const. del, tau, xk





implicit none

    type(type_gl) :: gl


double precision, dimension(30, 30):: d2aSRKdxidxj
integer:: i, j

d2aSRKdxidxj = 0.D0

Do i = 1, gl%ncomp - 1
    Do j = 1, gl%ncomp - 1
        d2aSRKdxidxj(i,j) = gl%aij_SRK(j,i) + gl%aij_SRK(i,j) - gl%aij_SRK(j,gl%ncomp) - gl%aij_SRK(gl%ncomp,j) &
                        & + 2.D0 * gl%aij_SRK(gl%ncomp,gl%ncomp) - gl%aij_SRK(gl%ncomp,i) - gl%aij_SRK(i,gl%ncomp)
    End do
End do

end subroutine d2a_SRK_dxidxj







module subroutine db_SRK_dxi(gl,dbSRKdxi)
!**************************************************************************
!This subroutine calculates the partial derivative of the mixed parameter b
!of the SRK equation of state with respect to the compositions xi
!
!OUTPUT: vector dbSRKdxi, containing the derivatives db / dxi at const. del, tau, xk





implicit none

    type(type_gl) :: gl


double precision, dimension(30):: dbSRKdxi
integer:: i, j

dbSRKdxi = 0.D0

if ((gl%mix_type == 2) .or. (gl%mix_type == 22)) then
    Do i = 1, gl%ncomp - 1
        dbSRKdxi(i) = gl%bi_SRK(i) - gl%bi_SRK(gl%ncomp)
    End do
!Andreas, Feb.2016. Quadratic mixing rules for b -> mix_type = 21
elseif (gl%mix_type == 21) then
    Do i = 1, gl%ncomp - 1
        Do j = 1, gl%ncomp - 1
            dbSRKdxi(i) = dbSRKdxi(i) + gl%molfractions(j)*(gl%bij_SRK(j,i) + gl%bij_SRK(i,j) - gl%bij_SRK(j,gl%ncomp) - gl%bij_SRK(gl%ncomp,j))
        End do
        dbSRKdxi(i) = dbSRKdxi(i) - 2.D0 * gl%molfractions(gl%ncomp)*gl%bij_SRK(gl%ncomp,gl%ncomp) + & 
                    & gl%molfractions(gl%ncomp) * (gl%bij_SRK(i,gl%ncomp) + gl%bij_SRK(gl%ncomp,i))
    End do
end if

end subroutine db_SRK_dxi
!**************************************************************************


!Andreas, January 2016
!**************************************************************************
module subroutine d2b_SRK_dxidxj(gl,d2bSRKdxidxj)
!**************************************************************************
!This subroutine calculates the second derivative of the mixed parameter b
!of the SRK equation of state with respect to the compositions xi and xj
!
!OUTPUT: matrix d2bSRKdxidxj, containing the derivatives d2b / dxidxj at const. del, tau, xk





implicit none

    type(type_gl) :: gl


double precision, dimension(30, 30):: d2bSRKdxidxj
integer:: i, j

d2bSRKdxidxj = 0.D0

!Andreas, Feb.2016. Quadratic mixing rules for b -> mix_type = 21
if (gl%mix_type == 21) then
    Do i = 1, gl%ncomp - 1
        Do j = 1, gl%ncomp - 1
            d2bSRKdxidxj(i,j) = gl%bij_SRK(j,i) + gl%bij_SRK(i,j) - gl%bij_SRK(j,gl%ncomp) - gl%bij_SRK(gl%ncomp,j) &
                            & + 2.D0 * gl%bij_SRK(gl%ncomp,gl%ncomp) - gl%bij_SRK(gl%ncomp,i) - gl%bij_SRK(i,gl%ncomp)
        End do
    End do
end if

end subroutine d2b_SRK_dxidxj
!**************************************************************************


!!**************************************************************************
!subroutine dar_dxi_SRK(Temp, dens, dardxi_SRK)
!!**************************************************************************
!!Partial derivative of the reduced Helmholtz energy (using the SRK!) with 
!!respect to xi at constant tau, del and xk
!!
!!OUTPUT: vector dardxi_SRK, containing the derivatives dar / dxi at const. del, tau, xk
!


!

!implicit none
!
!    type(type_gl) :: gl
!
!
!double precision:: TEMP, DENS 
!double precision, dimension(30):: dardxi_SRK
!double precision, dimension(30)::dbSRKdxi, daSRKdxi
!integer:: i
!double precision:: delta, tau
!
!
!delta = Dens / gl%rhored_SRK
!tau = gl%Tred_SRK / Temp
!
!dardxi_SRK = 0.D0
!call da_SRK_dxi(gl,daSRKdxi)
!call db_SRK_dxi(gl,dbSRKdxi)
!
!Do i = 1, gl%ncomp - 1
!    dardxi_SRK(i) = gl%rhored_SRK*delta / (1.D0 - gl%b_SRK * gl%rhored_SRK * delta) * dbSRKdxi(i) - &
!                & tau / gl%R_SRK / gl%tred_SRK * (dlog(1.D0+gl%b_SRK*gl%rhored_SRK*delta)*(gl%b_SRK*daSRKdxi(i)-gl%a_SRK*dbSRKdxi(i))/gl%b_SRK**2 + &
!                & gl%a_SRK / gl%b_SRK *gl%rhored_SRK * delta * dbSRKdxi(i) / (1.D0 + gl%b_SRK * gl%rhored_SRK * delta))
!End do
!
!end subroutine dar_dxi_SRK
!!**************************************************************************
!
!!**************************************************************************
!subroutine d2ar_dxidxj_SRK(Temp, dens, d2ardxidxj_SRK)
!!**************************************************************************
!!Second derivative of the reduced Helmholtz energy (using the SRK!) with 
!!respect to xi and xj
!!
!!OUTPUT: matrix d2ardxidxj_SRK, containing the derivatives d2ar / dxi dxj
!


!

!implicit none
!
!    type(type_gl) :: gl
!
!
!double precision:: TEMP, DENS 
!double precision, dimension(30,30):: d2ardxidxj_SRK
!double precision, dimension(30)::dbSRKdxi, daSRKdxi
!double precision, dimension(30,30):: d2aSRKdxidxj
!integer:: i, j
!double precision:: delta, tau
!
!
!delta = Dens / gl%rhored_SRK
!tau = gl%Tred_SRK / Temp
!
!d2ardxidxj_SRK = 0.D0
!call da_SRK_dxi(gl,daSRKdxi)
!call db_SRK_dxi(gl,dbSRKdxi)
!call d2a_SRK_dxidxj(gl,d2aSRKdxidxj)
!
!Do i = 1, gl%ncomp - 1
!    Do j = 1, gl%ncomp - 1
!        d2ardxidxj_SRK(i,j) = gl%rhored_SRK**2*delta**2/(1.D0-gl%b_SRK*gl%rhored_SRK*delta)**2*dbSRKdxi(i)*dbSRKdxi(j) - &
!                            & 1.D0/gl%R_SRK/Temp*(gl%rhored_SRK*delta/(1.D0+gl%b_SRK*gl%rhored_SRK*delta)*dbSRKdxi(j)* &
!                            & (gl%b_SRK*daSRKdxi(i)-gl%a_SRK*dbSRKdxi(i))/gl%b_SRK**2 + &
!                            & dlog(1.D0 + gl%b_SRK * gl%rhored_SRK* delta) * &
!                            & ((gl%b_SRK*d2aSRKdxidxj(i,j)-daSRKdxi(i)*dbSRKdxi(j))/gl%b_SRK**2 - &
!                            &  dbSRKdxi(i) * (gl%b_SRK**2*daSRKdxi(j)-2.D0*gl%b_SRK*dbSRKdxi(j)*gl%a_SRK)/gl%b_SRK**4) + &
!                            & gl%rhored_SRK*delta*dbSRKdxi(i)*(-gl%rhored_SRK*delta/(1.D0+ gl%b_SRK*gl%rhored_SRK*delta)**2 * &
!                            & dbSRKdxi(j) * gl%a_SRK / gl%b_SRK + 1.D0 / (1.D0 + gl%b_SRK * gl%rhored_SRK * delta) * &
!                            & (gl%b_SRK*daSRKdxi(j)-gl%a_SRK*dbSRKdxi(j))/gl%b_SRK**2))
!    End do  
!End do
!
!end subroutine d2ar_dxidxj_SRK
!**************************************************************************

!**************************************************************************
!subroutine d2ar_dxiddel_SRK(Temp, dens, d2ardxiddel_SRK)
!!**************************************************************************
!!Second derivative of the reduced Helmholtz energy (using the SRK!) with 
!!respect to xi and del
!!
!!OUTPUT: vector d2ardxiddel_SRK, containing the derivatives d2ar / d(xi)d(del) 
!!                                MULTIPLIED WITH DELTA!
!


!

!implicit none
!
!    type(type_gl) :: gl
!
!
!double precision:: TEMP, DENS 
!double precision, dimension(30):: d2ardxiddel_SRK
!double precision, dimension(30)::dbSRKdxi, daSRKdxi
!integer:: i
!double precision:: delta, tau
!
!
!delta = Dens / gl%rhored_SRK
!tau = gl%Tred_SRK / Temp
!
!d2ardxiddel_SRK = 0.D0
!call da_SRK_dxi(gl,daSRKdxi)
!call db_SRK_dxi(gl,dbSRKdxi)
!
!Do i = 1, gl%ncomp - 1
!    !d2ardxiddel_SRK(i) = dbSRKdxi(i) * rhored_SRK / (1.D0-b_SRK*rhored_SRK*delta)**2 - Tred_SRK / tau / R_SRK * &
!!                       & (b_SRK*rhored_SRK/(1.D0+b_SRK*rhored_SRK*delta)*(b_SRK*daSRKdxi(i)-a_SRK*dbSRKdxi(i)) / b_SRK**2 + &
!!                       & dbSRKdxi(i) * rhored_SRK / (1.D0+b_SRK*rhored_SRK*delta)**2 *a_SRK / b_SRK)
!    d2ardxiddel_SRK(i) = (dbSRKdxi(i) * gl%rhored_SRK / (1.D0-gl%b_SRK*gl%rhored_SRK*delta)**2 - tau / gl%tred_SRK / gl%R_SRK * &
!                       & ((1.D0+gl%b_SRK*gl%rhored_SRK*delta)*daSRKdxi(i)- &
!                       &  gl%a_SRK*gl%rhored_SRK*delta*dbSRKdxi(i))/(1.D0+gl%b_SRK*gl%rhored_SRK*delta)**2)*delta
!End do
!
!end subroutine d2ar_dxiddel_SRK
!**************************************************************************

!**************************************************************************
!subroutine d2ar_dxidtau_SRK(Temp, dens, d2ardxidtau_SRK)
!!**************************************************************************
!!Second derivative of the reduced Helmholtz energy (using the SRK!) with 
!!respect to xi and tau
!!
!!OUTPUT: vector d2ardxidtau_SRK, containing the derivatives d2ar / d(xi)d(tau) 
!


!

!implicit none
!
!    type(type_gl) :: gl
!
!
!double precision:: TEMP, DENS 
!double precision, dimension(30):: d2ardxidtau_SRK
!double precision, dimension(30)::dbSRKdxi, daSRKdxi, d2a_dxidtau
!double precision:: daSRK_dtau, da_SRK_dtau, dT_dtau
!integer:: i
!double precision:: delta, tau
!
!
!delta = Dens / gl%rhored_SRK
!tau = gl%Tred_SRK / Temp
!
!d2ardxidtau_SRK = 0.D0
!call da_SRK_dxi(gl,daSRKdxi)
!call db_SRK_dxi(gl,dbSRKdxi)
!call d2a_SRK_dxidtau(gl,Temp, d2a_dxidtau)
!
!daSRK_dtau = da_SRK_dtau(gl,Temp, 0)
!dT_dtau = - Temp / tau
!
!Do i = 1, gl%ncomp - 1
!    d2ardxidtau_SRK(i) = (-1.D0 / gl%R_SRK / Temp * & 
!                      & (log(1.D0 + gl%b_SRK*gl%rhored_SRK*delta) *(d2a_dxidtau(i)/gl%b_SRK-daSRK_dtau*dbSRKdxi(i)/gl%b_SRK**2) + &
!                      & daSRK_dtau *gl%rhored_SRK * delta / gl%b_SRK / (1.D0 + gl%b_SRK * gl%rhored_SRK*delta) * dbSRKdxi(i)) + &
!                      & dT_dtau / gl%R_SRK / Temp**2 * (dlog(1.D0+gl%b_SRK*gl%rhored_SRK*delta)*(gl%b_SRK*daSRKdxi(i)- & 
!                      & gl%a_SRK*dbSRKdxi(i))/gl%b_SRK**2 + gl%a_SRK / gl%b_SRK *gl%rhored_SRK * delta * dbSRKdxi(i) / & 
!                      & (1.D0 + gl%b_SRK * gl%rhored_SRK * delta))) * tau
!End do
!
!end subroutine d2ar_dxidtau_SRK
!**************************************************************************


!!!THIS ROUTINES IS NOT NEEDED, Andreas, Nov 2015
!!**************************************************************************
!subroutine dnar_dni_SRK(Temp, dens, dnardni_SRK)
!!**************************************************************************
!!Partial derivative of the reduced Helmholtz energy (using the SRK!) times
!!the total moles with respect to ni at constant T, V and nk
!!
!!OUTPUT: vector dnardni_SRK, containing the derivatives dnar / dni at const. T, V, nk
!


!
!implicit none
!
!double precision:: TEMP, DENS 
!double precision, dimension(30):: dnardni_SRK
!integer:: i, k
!double precision:: delta, tau, ar, dar_ddel
!double precision, dimension(30):: dardxi_SRK
!double precision, dimension(8):: MIXDERFNR
!integer, dimension(8):: GETDER
!
!
!delta = Dens / rhored_SRK
!tau = Tred_SRK / Temp
!
!dnardni_SRK = 0.D0
!
!GETDER = (/1,1,0,0,0,0,0,0/)
!call MIXDERIVSFNR_SRK(TEMP, DENS, GETDER, MIXDERFNR)
!ar = MIXDERFNR(1)
!dar_ddel = MIXDERFNR(2)
!
!call dar_dxi_SRK(Temp, dens, dardxi_SRK)
!
!
!Do i = 1, ncomp
!    dnardni_SRK(i) = ar + dar_ddel 
!    if (i < ncomp) dnardni_SRK(i) = dnardni_SRK(i)+ dardxi_SRK(i) 
!    Do k = 1, ncomp-1
!        dnardni_SRK(i) = dnardni_SRK(i) - molfractions(k) * dardxi_SRK(k)
!    End do
!End do
!
!end subroutine dnar_dni_SRK
!!**************************************************************************

!!**************************************************************************
!subroutine d2nar_dnidT_Vn_SRK(Temp, dens, d2nardnidT_SRK)
!!**************************************************************************
!!Second derivative of the reduced Helmholtz energy (using the SRK!) times
!!the total moles with respect to ni and T at constant V and n
!!
!!OUTPUT: vector d2nardnidT_SRK, containing the derivatives dnar / dni dT 
!


!
!implicit none
!
!double precision, intent(in):: TEMP, DENS 
!double precision, dimension(30), intent(out):: d2nardnidT_SRK
!integer:: i, k
!double precision:: delta, tau, ar, dar_ddel
!double precision, dimension(30):: dardxi_SRK
!double precision, dimension(8):: MIXDERFNR
!integer, dimension(8):: GETDER
!
!
!delta = Dens / rhored_SRK
!tau = Tred_SRK / Temp
!
!dnardni_SRK = 0.D0
!
!GETDER = (/1,1,0,0,0,0,0,0/)
!call MIXDERIVSFNR_SRK(TEMP, DENS, GETDER, MIXDERFNR)
!ar = MIXDERFNR(1)
!dar_ddel = MIXDERFNR(2)
!
!call dar_dxi_SRK(Temp, dens, dardxi_SRK)
!
!
!Do i = 1, ncomp
!    dnardni_SRK(i) = ar + dar_ddel 
!    if (i < ncomp) dnardni_SRK(i) = dnardni_SRK(i)+ dardxi_SRK(i) 
!    Do k = 1, ncomp-1
!        dnardni_SRK(i) = dnardni_SRK(i) - molfractions(k) * dardxi_SRK(k)
!    End do
!End do
!
!end subroutine d2nar_dnidT_Vn_SRK
!!**************************************************************************


!!THIS ROUTINES IS NOT NEEDED, Andreas, Nov 2015
!!**************************************************************************
!subroutine lnf_SRK(Temp, dens, lnfSRK)
!!**************************************************************************
!!The logarithm of the fugacities calculated from the SRK EOS
!!
!!OUTPUT: vector lnf_SRK, containing the fugacities of each component
!


!
!implicit none
!
!double precision:: TEMP, DENS 
!double precision, dimension(30):: lnfSRK
!integer:: i
!double precision:: delta, tau
!double precision, dimension(30):: dnardni
!
!lnfSRK = 0.D0
!
!delta = Dens / rhored_SRK
!tau = Tred_SRK / Temp
!
!call dnar_dni_SRK(Temp, dens, dnardni)
!
!Do i = 1, ncomp
!    lnfSRK(gl,i) = dlog(molfractions(i)*dens*R_SRK*Temp*1.D-6) + dnardni(i)
!End do
!
!    end subroutine lnf_SRK
!!**************************************************************************
    
    
    

Double Precision module Function PSRK(gl,T,P0_lc,iphase, nrsubst) 
!------------------------------------------------------------------------------------
!Solves SRK_equation
! P=R*T/(v-b) - a/(v*v+b*v)
! = >  v^3+(-r*t/p)*v^2+(-b^2-r*t*b/p+a/p)v+(-a*b/p) = 0
!Input: P(MPa), T(K), and phase status (1-liquid; 2-vapor)
!Output: density mol/m3 
!
! Theresa Wiens
! April 2010
! Bochum
!
! if psrk is negative, it means solving cubic equation is failed
!------------------------------------------------------------------------------------


	! double precision tc, pc, accen, wm, req(1) : crit.T+P, acentric factor, mole mass, Req(1)
	! double precision dimension (30) molfractions

	! double precision, dimension (30,3) ccoeff: Mathias Copeman
	
	

implicit none

    type(type_gl) :: gl

	
    !--Variable declaration--
    !--------------------------------------------------------------------------------
	integer :: iphase                                ! input iphase: 1-liquid; 2-vapor
	integer:: nrsubst                             ! if 0, then mixture, otherwise pure fluid
	integer::nc                                   ! number of components
	double precision :: T,P0_lc,p                         ! input: T,p0  transfered to p
	integer :: i,j,k                                   ! counter parameter
    double precision, dimension (30):: ai,bi,ci,aa  ! help parameters of cEOS
    double precision, dimension (30):: msrk,alfa_lc,tr ! help parameters of cEOS
    double precision :: a,a1,b, ex                     ! help parameters of cEOS
    double precision, dimension (30,30):: bij       ! help parameters of cEOS
    double precision :: ge, ge_res                     ! variable and function allocating mixture ge for cEOS
    double precision :: cc(3),roots(3)                 ! 3 parameters of cEOS, 3 roots
    !double precision :: wmm                            ! mixture molar weight    
    double precision :: vsrk,vfl,c                     ! parameters to correct liquid densities
    double precision :: Rmix                           ! mixture gas constant
    !  ------------------------------------------------------------------------------

    !--default values--
    !--------------------------------------------------------------------------------
    p = P0_lc*1.d6 / gl%factorpress   ! transfering MPa to Pa
    nc = gl%ncomp
    vfl = 0.d0
    a1 = 0.d0
    a = 0.d0
    b = 0.d0
    c = 0.d0
    ex = 4.d0 / 3.d0
    !wmm = 0.d0
    ge_res = 0.d0
    
    
    ai = 0.d0
    bi = 0.d0
    ci = 0.d0
    aa = 0.d0
    tr = 0.d0
    msrk = 0.d0
    alfa_lc = 0.d0    
    roots = 0.d0
    cc = 0.d0  
    bij= 0.d0
    

    !-- Program Start--
    !-------------------------------------------------------------------------------------------------------------           
    k = 1
    !--PSRK/RKS---loop start-----for each component in mixture-----------------------
    do while (k <= nc)
        ! check if calculation for pure fluid or mixture - J.G. 09.2010
        if (nrsubst /= 0) then
            k = nrsubst     ! DO Loop is only computed once, all module values are taken from component no. "nrsubst"
            nc = nrsubst
        end if
        tr(k)=T/gl%tc(k)                                               !reduced temperature
        
	    ! help parameters of cEOS
        aa(k)=.42748d0*gl%req(k)**2*gl%tc(k)**2/gl%pc(k)*1.d-6*gl%factorpress            ! /1.d6:  pc(i) [MPa] ->  [Pa]
        
	    if (.NOT.gl%psrk_ok .OR. tr(k) > 1 .OR. gl%ccoeff(k,1) == 0) then     ! if no PSRK parameters are available the RKS cEOS is used
	        msrk(k)=0.48d0+1.574d0*gl%accen(k)-0.176d0*gl%accen(k)**2
	        alfa_lc(k)=(1.d0+msrk(k)*(1.d0-tr(k)**0.5d0))**2
	        
	    else  !PSRK: Mathias Copemann coefficients
	        alfa_lc(k) = (1.D0 + gl%ccoeff(k,1)*(1.D0-tr(k)**0.5d0) + gl%ccoeff(k,2)*(1.D0-tr(k)**0.5d0)**2 &
	                & + gl%ccoeff(k,3)*(1.D0-tr(k)**0.5d0)**3)**2
	    endif
    	
	    ai(k)=aa(k)*alfa_lc(k)
	    bi(k)=.08664d0*gl%req(k)*gl%tc(k)/gl%pc(k)*1.d-6*gl%factorpress 
	    ci(k)=.40768d0*gl%req(k)*gl%tc(k)/gl%pc(k)*1.d-6 *gl%factorpress * (.29441d0-gl%pc(k)*1.d6/gl%factorpress /gl%rhoc(k)/gl%req(k)/gl%tc(k))!*wm(i) in mol/m^3 J.G. ! Peneloux correction in m³/kg
    	
    	
    	
	    ! mixing rules-----------------------------------------
	    if (.NOT.gl%psrk_ok) then  
	        a1 = a1 + gl%molfractions(k)* ai(k)**.5d0
    	else  !PSRK              
    	    a1 = a1 + gl%molfractions(k)*ai(k)/bi(k)
    	end if
    	
    	
    	!original b mixing rule still used for SRK
    	if (.NOT.gl%psrk_ok) then 
	        b = b + gl%molfractions(k)*bi(k)
	    else
	            
	    end if
	    	    
	    c = c + gl%molfractions(k)*ci(k)   ! Peneloux correction
        ! end of mixing rules-----------------------------------
        
	    ! calculating mixture molar weight
	    !wmm = wmm + wm(i)*molfractions(i)
	    k = k + 1
	end do	
    !--loop end-------------------
        
    
    !--mixing rules--
    !if (.NOT.psrk_ok) then    
        a = a1**2
        
    !else   !PSRK
    !    
    !    !Li mixing rule for parameter b
    !    do i=1, ncomp
    !        do j=1, ncomp
    !            bij(i,j) = (0.5d0*(bi(i)**.75d0 + bi(j)**.75d0))**ex
    !            b = b + molfractions(i)*molfractions(j)*bij(i,j)
    !        end do
    !    end do
    !    
    !    ge_res = ge(T)      ! calling ge function
    !    a = (a1 + ge_res/(-0.64663d0)) * b
    !end if
    
    if (nrsubst /= 0) then  ! pure fluid SRK - J.G. 09.2010
        a = ai(nrsubst)
        b = bi(nrsubst)	        
	    c = ci(nrsubst)    
    end if
    !--PSRK/RKS---end----------------------------------------------------------------
    

    
	! parameters of cubic function
	if (nrsubst > 0) then
        cc(1) = -gl%req(nrsubst)*t/p
        cc(2) = -b**2-gl%req(nrsubst)*t*b/p+a/p
        cc(3) = -a*b/p
	else
	    call R_mix_calc(gl,Rmix)
	    cc(1) = -Rmix*t/p
        cc(2) = - b**2 - Rmix*t*b/p+a/p
        cc(3) = -a*b/p
    end if
	!write (6,996) a,b!, cc(1), cc(2), cc(3)
    !996 format(/,f12.9,2x,f12.9)!,/,f15.6,2x,f15.6,2x,f15.6)
    
    call cubic_nt(gl,cc,roots)                 ! Solves z^3 + cc(1)*z^2+cc(2)*z+cc(3) = 0
    
    
    !--handling roots----------------------------------------------------------------
	if (roots(1) < -10.d0)then 
    	psrk=-1.d0                           ! failed solving cubic equation
    	
    else
      if (roots(1) < -0.d0)then               ! one real root
    	                                    ! vSRK=roots(2) !no matter liquid or vapor since only one real root
        VSRK=1.d0/(roots(2))!/wmm)          !mol/m^3 J.G. ! kg/m3
        
      else                                  ! three real roots
        if (iphase == 1)then                  ! vSRK=roots(3) !liquid volume m3/mol
          VSRK=1.d0/(roots(3))!/wmm)          !mol/m^3 J.G. ! kg/m3
          
        else                                ! vSRK=roots(2) !vapor volume
          VSRK=1.d0/(roots(2))!/wmm)          !mol/m^3 J.G. ! kg/m3
        endif
        
      endif
      
    endif
    
    
    !--correction for liquid density (Péneloux)
    if (iphase == 1) then!.AND.psrk_ok) then 
        vfl = 1/vsrk    ! in m^3/mol   m³/kg
        c = c!/wmm
        vfl = vfl - c 
        vsrk = 1/vfl    ! transforming volume into density again
    endif
    
    psrk=vsrk
End Function

    
    

module subroutine CUBIC_NT(gl,para,volume)
!----------------------------------------------------------
! Algorithm: arccos method
! Theresa Wiens
! April 2010
! Bochum
!-----------------------------------------------------------


implicit none

    type(type_gl) :: gl

	
    ! Variable declaration:
    !  ------------------------------------------------------------------------------
	double precision :: para(3)                ! a,b,c input parameters
    double precision :: volume(3)              ! output [m³/mol]
    ! volume(1)- the number of real roots: < 0: one real root, > 0: three
    ! volume(2)- the biggest root
    ! volume(3)- the smallest root. if the number of real roots is 1, then volume(3)=0
	double precision :: a,b,c,qq,pp,nn
	double precision :: r(3)                   ! roots of cubic eqn
    double precision :: pie,aa,bb
	integer :: nr                            ! number of real roots
    !  ------------------------------------------------------------------------------
    	
	pie=3.14159265358979d0
	
	!x^3+ax^2+bx+c=0
	a=para(1)
	b=para(2)
	c=para(3)

	pp=(a*a-3.d0*b)/9.d0    ! = -p/3
    qq=(2.d0*a**3-9.d0*a*b+27.d0*c)/54.d0   ! = q/2



    !--case differentiation regarding D-----------D = -pp^3 + qq^2-----------------------------------------------
    
    if (pp**3 > qq**2) then ! D < 0: casus irreducibilis and pp > 0 !!!------------
      nn=dacos(-qq/pp**1.5d0)
      !nn=dacos(-qq/sqrt(pp**3))
      nr=3   ! there are 3 real solutions
      ! roots
      r(1)=2.d0*pp**0.5d0*dcos(nn/3.d0)-a/3.d0
      r(2)=2.d0*pp**0.5d0*dcos((nn+2.d0*pie)/3.d0)-a/3.d0
      r(3)=2.d0*pp**0.5d0*dcos((nn+4.d0*pie)/3.d0)-a/3.d0
        
      
    else ! D >= 0 ---------------------------------------------------------------
      nr=-1  ! there is only one real root, perhaps additional dual root which is ignored
     
     ! u or v:
      aa=-sign(1.d0,qq)*(abs(qq)+(qq**2-pp**3)**0.5d0)**(1.d0/3.d0)
      
      if (aa == 0.d0)then ! triple root
        bb=0.d0 !v
      else
        bb=pp/aa !v
      endif 
       
         
      r(1)=(aa+bb)-a/3.d0
      r(2)=0.d0
      r(3)=0.d0
	endif
	!--end of case differentiation concerning D-------------------------------------------------
    
    
    !--creating output vector-------------------------------------------------------------------  
    if (nr < 0d0) then                            ! there is only one real root  WAS IST MIT D=0 ?
        volume(1)=-1d0
        volume(2)=r(1)
    	volume(3)=0d0   
    	 
	else	! D < 0                             ! three real roots
		volume(1)=3d0
		volume(2)=max(r(2),r(3),r(1))           ! vapor volume is the biggest root
        
        ! in case some roots are negative   needed since minimum > 0 is searched for!
          qq=1.d10
          pp=1.d10
          nn=1.d10
          if (r(1) > 0) then
            qq=r(1)
          endif
          if (r(2) > 0) then
            pp=r(2)
          endif
          if (r(3) > 0) then
            nn=r(3)
          endif
          volume(3)=min(qq,pp,nn)        
	endif
	!--end creating output vector---------------------------------------------------------------
    end subroutine

    
module subroutine init_PR(gl,Temp, nrsubst)






implicit none

    type(type_gl) :: gl


double precision:: temp
integer :: i, j, nrsubst
double precision :: tau

!It is very important to initialize these variables here (set to 0)
gl%a_PR = 0.D0
gl%b_PR = 0.D0

!If whole mixture properties shall be calculated -->  nrsubst = 0
if (nrsubst == 0) then 
    Do i = 1, gl%ncomp
        gl%ac_PR(i) = 0.45724D0 * gl%R_PR**2 * gl%tc(i)**2 / (gl%pc(i)*1.D6)
        !Added missing different calculation of kappai for an acentric factor acceni > 0.491. Andreas Jäger, January 2017
        If (gl%accen(i) .le. 0.491D0) then
            gl%kappa(i) = 0.37464D0 + 1.54226D0 * gl%accen(i) - 0.26992D0 * gl%accen(i)**2
        else
            gl%kappa(i) = 0.379642D0 + 1.48503D0 * gl%accen(i) - 0.164423D0 * gl%accen(i)**2 + 0.016666D0 * gl%accen(i)**3
        end if
    End do

    do i = 1, gl%ncomp
        !New version: Mathias and Copeman (Fluid Phase Equilib., 13, 91-108, 1983) temperature dependence for the parameter ai
        if ((dabs(gl%Cji_cubic(1,i)) < 1.D-12) .and. (dabs(gl%Cji_cubic(2,i)) < 1.D-12) .and. (dabs(gl%Cji_cubic(3,i)) < 1.D-12)) then   !If no Mathias Copeman parameters are given, then use classical formulation by setting C1i = kappai
            gl%Cji_cubic(1,i) = gl%kappa(i)        
        end if
        gl%ai_PR(i) = gl%ac_PR(i) * (1.D0 + gl%Cji_cubic(1,i)*(1.D0-(Temp/gl%tc(i))**0.5D0) & 
                                  & + gl%Cji_cubic(2,i)*(1.D0-(Temp/gl%tc(i))**0.5D0)**2 & 
                                  & + gl%Cji_cubic(3,i)*(1.D0-(Temp/gl%tc(i))**0.5D0)**3 )**2
        !Old "classical" calculation of the parameter ai        
        !ai_PR(i) = ac_PR(i) * (1.D0 + kappa(i) * (1.D0 - (Temp / tc(i))**0.5D0))**2
        gl%bi_PR(i) = 0.0778D0*gl%R_PR * gl%tc(i) / (gl%pc(i)*1.D6)
    end do
    
    do j = 1, gl%ncomp
        do i = 1, gl%ncomp
            gl%aij_PR(i,j) = (1.D0 - gl%kij_PR(i,j)) * (gl%ai_PR(i) * gl%ai_PR(j))**0.5D0
            gl%a_PR = gl%a_PR + gl%molfractions(i) * gl%molfractions(j) * gl%aij_PR(i,j)
        end do
        if (gl%mix_type == 3) then
            gl%b_PR = gl%b_PR + gl%molfractions(j) * gl%bi_PR(j)
        end if
    end do

    !Andreas, February 2016. Quadratic mixing rules for b
    if (gl%mix_type == 31) then
        do j = 1, gl%ncomp
            do i = 1, gl%ncomp
                gl%bij_PR(i,j) = (gl%bi_PR(i) * gl%bi_PR(j))**0.5D0 * (1.D0 - gl%lij_PR(i,j)) 
                gl%b_PR = gl%b_PR + gl%molfractions(i) * gl%molfractions(j) * gl%bij_PR(i,j)
            end do
        end do    
    end if    
    
    !Important, Tred_PR must be set to 1 here, when full PR model is used
    !Andreas Jäger, Nov. 2015
    gl%Tred_PR = 1.D0
else 
    !-----------------------------------------------------------------------------------------------
    !Andreas Nov 2015. 
    
    if ((gl%Mix_Type == 1) .or. (gl%Mix_Type == 12)  .or. (gl%mix_type == 13)) then   !PR in Helmholtz mixture
        gl%Tred_PR = gl%TredMix
    else                        !PR with PR mixing rules
        !Tred_PR = Tc(1)         !Andreas Jäger, July 2016
        gl%Tred_PR = gl%Tc(nrsubst) 
    end if
    
    tau = gl%Tred_PR / Temp
    !-----------------------------------------------------------------------------------------------
    
    
    !If a pure component in a mixture with the PR is calculated, only the parameters of this pure component need to be calculated
    gl%ac_PR(nrsubst) = 0.45724D0 * gl%R_PR**2 * gl%tc(nrsubst)**2 / (gl%pc(nrsubst)*1.D6)
    !Added missing different calculation of kappai for an acentric factor acceni > 0.491. Andreas Jäger, January 2017
    If (gl%accen(nrsubst) .le. 0.491D0) then
        gl%kappa(nrsubst) = 0.37464D0 + 1.54226D0 * gl%accen(nrsubst) - 0.26992D0 * gl%accen(nrsubst)**2
    else
        gl%kappa(nrsubst) = 0.379642D0 + 1.48503D0 * gl%accen(nrsubst) - 0.164423D0 * gl%accen(nrsubst)**2 + 0.016666D0 * gl%accen(nrsubst)**3    
    end if
    !New version: Mathias and Copeman (Fluid Phase Equilib., 13, 91-108, 1983) temperature dependence for the parameter ai
    if ((dabs(gl%Cji_cubic(1,nrsubst)) < 1.D-12) .and. (dabs(gl%Cji_cubic(2,nrsubst)) < 1.D-12) .and. (dabs(gl%Cji_cubic(3,nrsubst)) < 1.D-12)) then   !If no Mathias Copeman parameters are given, then use classical formulation by setting C1i = mi
        gl%Cji_cubic(1,nrsubst) = gl%kappa(nrsubst)        
    end if
    gl%ai_PR(nrsubst) = gl%ac_PR(nrsubst) * (1.D0 + gl%Cji_cubic(1,nrsubst)*(1.D0-(1.d0/tau)**0.5D0) & 
                                & + gl%Cji_cubic(2,nrsubst)*(1.D0-(1.d0/tau)**0.5D0)**2 & 
                                & + gl%Cji_cubic(3,nrsubst)*(1.D0-(1.d0/tau)**0.5D0)**3 )**2
    !Old "classical" calculation of the parameter ai
    !ai_PR(nrsubst) = ac_PR(nrsubst) * (1.D0 + kappa(nrsubst) * (1.D0 - (1.D0 / tau)**0.5D0))**2
    gl%bi_PR(nrsubst) = 0.0778D0 * gl%R_PR * gl%tc(nrsubst) / (gl%pc(nrsubst)*1.D6) 
    gl%a_PR = gl%ai_PR(nrsubst)
    gl%b_PR = gl%bi_PR(nrsubst)
End if


!Save the temperature the PR was last initialized with
gl%Temp_init = Temp 

end subroutine init_PR
    
    
double precision module function rho_PR (gl,Temp, press_in, Iphase, nrsubst)
    




implicit none

    type(type_gl) :: gl

    
double precision:: cc(3),roots(3)                 ! 3 parameters of cEOS, 3 roots
double precision:: Temp, press, press_in
integer:: Iphase, nrsubst
    
call init_PR(gl,Temp, nrsubst)
        
roots = 0.D0
press = press_in * 1.D6    !Mpa -- >  Pa
    
cc(1) = -gl%R_PR*Temp/press+gl%b_PR
cc(2) = -3.d0*gl%b_PR**2-2.d0*gl%R_PR*temp*gl%b_PR/press+gl%a_PR/press   
cc(3) = (-gl%a_PR*gl%b_PR+gl%b_PR**2*gl%R_PR*temp)/press+gl%b_PR**3
    
call cubic_nt(gl,cc,roots)                 ! Solves z^3 + cc(1)*z^2+cc(2)*z+cc(3) = 0
    
!--handling roots----------------------------------------------------------------
if (roots(1) < -10.d0)then 
    rho_PR = 0.d0                       ! failed solving cubic equation
    	
else
    if (roots(1) < -0.d0)then               ! one real root
    	                                ! vPR=roots(2) !no matter liquid or vapor since only one real root
    rho_PR=1.d0/(roots(2))!/wmm)          !mol/m^3 J.G. ! kg/m3
        
    else                                  ! three real roots
    if (iphase == 1)then                  ! vPR=roots(3) !liquid volume m3/mol
        rho_PR=1.d0/(roots(3))!/wmm)          !mol/m^3 J.G. ! kg/m3
          
    else                                ! vPR=roots(2) !vapor volume
        rho_PR=1.d0/(roots(2))!/wmm)          !mol/m^3 J.G. ! kg/m3
    endif
        
    endif
      
endif

if (rho_PR > 1.d0/gl%b_PR) then
    rho_PR = 0.d0              ! solution of root is larger than pole
endif

end function rho_PR

    
double precision module function da_PR_dtau (gl,Temp, nrsubst)
!**************************************************************************
!This subroutine calculates the partial derivative of the mixed parameter a
!of the PR equation of state with respect to the inversed reduces temperature tau
!
!OUTPUT: scalar daPRdT, containing the derivative da / dtau at const. del, xk
!**************************************************************************






implicit none

    type(type_gl) :: gl

      
double precision, dimension(30):: daiPR_dT      !Temperature derivative of the parameter a of fluid i
double precision, dimension(30,30)::daijPR_dT   !Temperature derivative of the binary parameter aij
double precision :: daPR_dT 
integer:: i, j
double precision:: tau, temp

integer:: nrsubst

Tau = 0.d0
da_PR_dtau = 0.D0
daiPR_dT = 0.D0
daijPR_dT = 0.D0
daPR_dT = 0.d0


!Belege Tred_PR, added new models, Andreas Jäger, January 2018
if ((gl%Mix_Type == 1) .or. (gl%Mix_Type == 12) .or. (gl%Mix_Type == 13)) then    !PR als Helmholtz im Gemisch, GERG-Mischungsregeln werden verwendet
    gl%Tred_PR = gl%TredMix
else    !PR Mischungsregeln werden verwendet
    if (gl%ncomp == 1) then
        !Tred_PR = Tc(1)         !Andreas Jäger, July 2016
        gl%Tred_PR = gl%Tc(nrsubst) 
    else
        gl%Tred_PR = 1.d0
    end if
end if

!Berechne Tau
tau = gl%Tred_PR / Temp
    
!Reinstoffe    
if (nrsubst /= 0) then    
    !daPR_dT = -ac_PR(nrsubst) * (1.d0 + kappa(nrsubst) * (1.d0 - (Temp / tc(nrsubst))**0.5d0)) &
    !        & * kappa(nrsubst) / tc(nrsubst) / (Temp / tc(nrsubst))**0.5d0
    !da_PR_dtau = -daPR_dT * Temp / tau
    !Correction for PR in Helmholtz mixtures, Andreas Nov 2015
    da_PR_dtau = gl%ac_PR(nrsubst) * (1.d0 + gl%kappa(nrsubst) * (1.d0 - (1.D0/tau) ** 0.5d0)) &
                & * gl%kappa(nrsubst) * (1.D0/tau)**(-0.5D0) / tau**2    
    
!Gemische
else
    do i = 1, gl%ncomp
        daiPR_dT(i) = -gl%ac_PR(i) * (1.d0 + gl%kappa(i) * (1.d0 - (Temp / gl%tc(i))**0.5d0)) &
                & * gl%kappa(i) / gl%tc(i) / (Temp / gl%tc(i))**0.5d0
    end do
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            daijPR_dT(i, j) = (1.D0 - gl%kij_PR(i, j)) / 2.D0 / (gl%ai_PR(i) * gl%ai_PR(j))**0.5d0 * &
                            & (gl%ai_PR(i) * daiPR_dT(j) + gl%ai_PR(j) * daiPR_dT(i))
            daPR_dT = daPR_dT + gl%molfractions(i) * gl%molfractions(j) * daijPR_dT(i, j)
        end do
    end do
    da_PR_dtau = -daPR_dT * (gl%Tred_PR/tau) / tau
end if

end function da_PR_dtau 

    
double precision module function d2a_PR_dtau2(gl,Temp, nrsubst)
!**************************************************************************
!This function calculates the second partial derivate of the mixed parameter a of PR with respect to the inversed reduces temperature tau.
!**************************************************************************








implicit none

    type(type_gl) :: gl



double precision :: Temp
double precision :: Tau
double precision :: daPR_dT
double precision :: d2aPR_dT2
double precision, dimension(30, 30) :: daijPR_dT
double precision, dimension(30, 30) :: d2aijPR_dT2
double precision, dimension(30) :: daiPR_dT
double precision, dimension(30) :: d2aiPR_dT2
integer :: i
integer :: j

integer:: nrsubst

Tau = 0.d0
d2a_PR_dtau2 = 0.d0
daPR_dT = 0.d0
d2aPR_dT2 = 0.d0
daijPR_dT = 0.d0
d2aijPR_dT2 = 0.d0
daiPR_dT = 0.d0
d2aiPR_dT2 = 0.d0
i = 0
j = 0


!Belege Tred_PR, added mix types 12 and 13, Andreas Jäger, January 2018
if ((gl%Mix_Type == 1) .or. (gl%Mix_Type == 12) .or. (gl%Mix_Type == 13)) then    !PR als Helmholtz im Gemisch, GERG-Mischungsregeln werden verwendet
    gl%Tred_PR = gl%TredMix
else    !PR Mischungsregeln werden verwendet
    if (gl%ncomp == 1) then
        !Tred_PR = Tc(1)         !Andreas Jäger, July 2016
        gl%Tred_PR = gl%Tc(nrsubst)  
    else
        gl%Tred_PR = 1.d0
    end if
end if

!Berechne Tau
Tau = gl%Tred_PR / Temp

!Pure substance (either really pure or in Helmholtz mixture)    
if (nrsubst /= 0) then 
    !daPR_dT = -ac_PR(nrsubst) * (1.d0 + kappa(nrsubst) * (1.d0 - (Temp / tc(nrsubst)) ** 0.5d0)) & 
    !        & * kappa(nrsubst) / tc(nrsubst) / (Temp / tc(nrsubst)) ** 0.5d0
    !d2aPR_dT2 = -ac_PR(nrsubst) * kappa(nrsubst) / tc(nrsubst) ** 0.5d0 * & 
    !        & (-0.5d0 / Temp ** (1.5d0) - 0.5d0 * kappa(nrsubst) / Temp ** (1.5d0))
    !d2a_PR_dtau2 = d2aPR_dT2 * Temp ** 2 / Tau ** 2 + 2.d0 * daPR_dT * Temp / Tau ** 2
    !Correction for PR in Helmholtz mixtures, Andreas Nov 2015
    daPR_dT = -gl%ac_PR(nrsubst) * (1.d0 + gl%kappa(nrsubst) * (1.d0 - (1.D0/tau) ** 0.5d0)) & 
            & * gl%kappa(nrsubst) / gl%tc(nrsubst) / (1.D0/tau) ** 0.5d0
    d2aPR_dT2 = -gl%ac_PR(nrsubst) * gl%kappa(nrsubst) / gl%tc(nrsubst) ** 0.5d0 * & 
            & (-0.5d0 / (gl%tc(nrsubst)/tau) ** (1.5d0) - 0.5d0 * gl%kappa(nrsubst) / (gl%tc(nrsubst)/tau) ** (1.5d0))
    d2a_PR_dtau2 = d2aPR_dT2 * (gl%tc(nrsubst)/tau) ** 2 / Tau ** 2 + 2.d0 * daPR_dT * (gl%tc(nrsubst)/tau) / Tau ** 2
    
!Full PR mixture model
else
    do i = 1, gl%ncomp
    daiPR_dT(i) = -gl%ac_PR(i) * (1.d0 + gl%kappa(i) * (1.d0 - (Temp / gl%tc(i)) ** 0.5d0)) & 
               & * gl%kappa(i) / gl%tc(i) / (Temp / gl%tc(i)) ** 0.5d0
    d2aiPR_dT2(i) = -gl%ac_PR(i) * gl%kappa(i) / gl%tc(i) ** 0.5d0 * & 
                 & (-0.5d0 / Temp ** (1.5d0) - 0.5d0 * gl%kappa(i) / Temp ** (1.5d0))
    end do
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            daijPR_dT(i, j) = (1.D0 - gl%kij_PR(i, j)) / (2.D0 * (gl%ai_PR(i) * gl%ai_PR(j))**0.5d0) * &
                            & (gl%ai_PR(i) * daiPR_dT(j) + gl%ai_PR(j) * daiPR_dT(i))
            d2aijPR_dT2(i, j) = ((gl%kij_PR(i, j) - 1.D0) * (-2.D0 * gl%ai_PR(i) * gl%ai_PR(j) * (daiPR_dT(i) * &
                              & daiPR_dT(j) + gl%ai_PR(i) * d2aiPR_dT2(j)) + (gl%ai_PR(j)**2 * ((daiPR_dT(i))**2 - &
                              & 2.D0 * gl%ai_PR(i) * d2aiPR_dT2(i)) + gl%ai_PR(i)**2 * (daiPR_dT(j))**2))) / &
                              & (4.D0 * (gl%ai_PR(i) * gl%ai_PR(j))**(1.5d0))
            daPR_dT = daPR_dT + gl%molfractions(i) * gl%molfractions(j) * daijPR_dT(i, j)
            d2aPR_dT2 = d2aPR_dT2 + gl%molfractions(i) * gl%molfractions(j) * d2aijPR_dT2(i, j)
        end do
    end do
    d2a_PR_dtau2 = d2aPR_dT2 * Temp ** 2 / Tau ** 2 + 2.d0 * daPR_dT * Temp / Tau ** 2
end if

end function d2a_PR_dtau2

    
double precision module function d3a_PR_dtau3(gl,Temp, nrsubst)
!**************************************************************************
!This function calculates the third partial derivate of the mixed parameter a of PR with respect to the inversed reduces temperature tau.
!**************************************************************************






implicit none

    type(type_gl) :: gl



double precision :: Temp
double precision :: Tau
double precision :: daPR_dT
double precision :: d2aPR_dT2
double precision :: d3aPR_dT3
double precision, dimension(30, 30) :: daijPR_dT
double precision, dimension(30, 30) :: d2aijPR_dT2
double precision, dimension(30, 30) :: d3aijPR_dT3
double precision, dimension(30) :: daiPR_dT
double precision, dimension(30) :: d2aiPR_dT2
double precision, dimension(30) :: d3aiPR_dT3
integer :: i
integer :: j

integer:: nrsubst



Tau = 0.d0
d3a_PR_dtau3 = 0.d0
daPR_dT = 0.d0
d2aPR_dT2 = 0.d0
d3aPR_dT3 = 0.d0
daijPR_dT = 0.d0
d2aijPR_dT2 = 0.d0
d3aijPR_dT3 = 0.d0
daiPR_dT = 0.d0
d2aiPR_dT2 = 0.d0
d3aiPR_dT3 = 0.d0
i = 0
j = 0



!Belege Tred_PR, added mix types 12 and 13, Andreas Jäger, January 2018
if ((gl%Mix_Type == 1) .or. (gl%Mix_Type == 12) .or. (gl%Mix_Type == 13)) then    !PR als Helmholtz im Gemisch, GERG-Mischungsregeln werden verwendet
    gl%Tred_PR = gl%TredMix
else    !PR Mischungsregeln werden verwendet
    if (gl%ncomp == 1) then
        !Tred_PR = Tc(1)         !Andreas Jäger, July 2016
        gl%Tred_PR = gl%Tc(nrsubst) 
    else
        gl%Tred_PR = 1.d0
    end if
end if

!Berechne Tau
Tau = gl%Tred_PR / Temp

!Pure substance (either really pure or in Helmholtz mixture)    
if (nrsubst /= 0) then 
    daPR_dT = -gl%ac_PR(nrsubst) * (1.d0 + gl%kappa(nrsubst) * (1.d0 - (1.D0/tau) ** 0.5d0)) &
            & * gl%kappa(nrsubst) / gl%tc(nrsubst) / (1.D0/tau) ** 0.5d0
    d2aPR_dT2 = -gl%ac_PR(nrsubst) * gl%kappa(nrsubst) / gl%tc(nrsubst) ** 0.5d0 * & 
              & (-0.5d0 / (gl%tc(nrsubst)/tau) ** (1.5d0) - 0.5d0 * gl%kappa(nrsubst) / (gl%tc(nrsubst)/tau) ** (1.5d0))
    d3aPR_dT3 = -gl%ac_PR(nrsubst) * gl%kappa(nrsubst) / gl%tc(nrsubst) ** 0.5d0 * &
              & (0.75d0 / (gl%tc(nrsubst)/tau) ** (2.5d0) + 0.75d0 * gl%kappa(nrsubst) / &
              & (gl%tc(nrsubst)/tau) ** (2.5d0))
    d3a_PR_dtau3 = -d3aPR_dT3 * (gl%tc(nrsubst)/tau) ** 3 / Tau ** 3 - 6.d0 * d2aPR_dT2 * (gl%tc(nrsubst)/tau) ** 2 / Tau ** 3 - &
                 & 6.d0 * daPR_dT * (gl%tc(nrsubst)/tau) / Tau ** 3
    
!Full PR mixture model
else
    do i = 1, gl%ncomp
        daiPR_dT(i) = -gl%ac_PR(i) * (1.d0 + gl%kappa(i) * (1.d0 - (Temp / gl%tc(i)) ** 0.5d0)) * gl%kappa(i) / gl%tc(i) / (Temp / gl%tc(i)) ** 0.5d0
        d2aiPR_dT2(i) = -gl%ac_PR(i) * gl%kappa(i) / gl%tc(i) ** 0.5d0 * &
                      & (-0.5d0 / Temp ** (1.5d0) - 0.5d0 * gl%kappa(i) / Temp ** (1.5d0))
        d3aiPR_dT3(i) = -gl%ac_PR(i) * gl%kappa(i) / gl%tc(i) ** 0.5d0 * (0.75d0 / Temp ** (2.5d0) + 0.75d0 * gl%kappa(i) / &
                    & Temp ** (2.5d0))
    end do
    do i = 1, gl%ncomp
        do j = 1, gl%ncomp
            daijPR_dT(i, j) = 0.5d0 * (1.d0 - gl%kij_PR(i, j)) * (gl%ai_PR(i) * gl%ai_PR(j)) ** (-0.5d0) * &
                            & (gl%ai_PR(i) * daiPR_dT(j) + gl%ai_PR(j) * daiPR_dT(i))
            d2aijPR_dT2(i, j) = 0.5d0 * (1.d0 - gl%kij_PR(i, j)) * &
                            & (-0.5d0 / (gl%ai_PR(i) * gl%ai_PR(j)) ** (1.5d0) * &
                            & (gl%ai_PR(i) * daiPR_dT(j) + daiPR_dT(i) * gl%ai_PR(j)) ** 2 + &
                            & (gl%ai_PR(i) * gl%ai_PR(j)) ** (-0.5d0) * &
                            & (d2aiPR_dT2(i) * gl%ai_PR(j) + d2aiPR_dT2(j) * gl%ai_PR(i) + &
                            & 2.d0 * daiPR_dT(i) * daiPR_dT(j)))
            d3aijPR_dT3(i, j) = 0.5d0 * (1.d0 - gl%kij_PR(i, j)) * &
                            & (0.75d0 / (gl%ai_PR(i) * gl%ai_PR(j)) ** (2.5d0) * &
                            & (gl%ai_PR(i) * daiPR_dT(j) + daiPR_dT(i) * gl%ai_PR(j)) ** 3 - &
                            & 1.d0/(gl%ai_PR(i) * gl%ai_PR(j)) ** (1.5d0) * &
                            & (gl%ai_PR(i) * daiPR_dT(j) + daiPR_dT(i) * gl%ai_PR(j)) * &
                            & (d2aiPR_dT2(i) * gl%ai_PR(j) + gl%ai_PR(i) * d2aiPR_dT2(j) + 2.d0 * daiPR_dT(i) * daiPR_dT(j)) + &
                            & 1.d0/(gl%ai_PR(i) * gl%ai_PR(j)) ** (0.5d0) * &
                            & (d3aiPR_dT3(i) * gl%ai_PR(j) + gl%ai_PR(i) * d3aiPR_dT3(j) + &
                            & 3.d0 * d2aiPR_dT2(i) * daiPR_dT(j) + 3.d0 * daiPR_dT(i) * d2aiPR_dT2(j)) - &
                            & 0.5d0 * (gl%ai_PR(i) * daiPR_dT(j) + daiPR_dT(i) * gl%ai_PR(j)) / &
                            & (gl%ai_PR(i) * gl%ai_PR(j)) ** (1.5d0) * &
                            & (d2aiPR_dT2(i) * gl%ai_PR(j) + gl%ai_PR(i) * d2aiPR_dT2(j) + 2.d0 * daiPR_dT(i) * daiPR_dT(j)))
            daPR_dT = daPR_dT + gl%molfractions(i) * gl%molfractions(j) * daijPR_dT(i, j)
            d2aPR_dT2 = d2aPR_dT2 + gl%molfractions(i) * gl%molfractions(j) * d2aijPR_dT2(i, j)
            d3aPR_dT3 = d3aPR_dT3 + gl%molfractions(i) * gl%molfractions(j) * d3aijPR_dT3(i, j)
        end do
    end do
    d3a_PR_dtau3 = -d3aPR_dT3 * Temp ** 3 / Tau ** 3 - 6.d0 * d2aPR_dT2 * Temp ** 2 / Tau ** 3 - &
                  & 6.d0 * daPR_dT * Temp / Tau ** 3
end if

end function d3a_PR_dtau3



!Andreas, November 2015
!**************************************************************************
module subroutine da_PR_dxi(gl,daPRdxi)
!**************************************************************************
!This subroutine calculates the partial derivative of the mixed parameter a
!of the PR equation of state with respect to the compositions xi
!
!OUTPUT: vector daPRdxi, containing the derivatives da / dxi at const. del, tau, xk





implicit none

    type(type_gl) :: gl


double precision, dimension(30):: daPRdxi
integer:: i, j

daPRdxi = 0.D0

Do i = 1, gl%ncomp - 1
    Do j = 1, gl%ncomp - 1
        daPRdxi(i) = daPRdxi(i) + gl%molfractions(j)*(gl%aij_PR(j,i) + gl%aij_PR(i,j) - gl%aij_PR(j,gl%ncomp) - gl%aij_PR(gl%ncomp,j))
    End do
    daPRdxi(i) = daPRdxi(i) - 2.D0 * gl%molfractions(gl%ncomp)*gl%aij_PR(gl%ncomp,gl%ncomp) + & 
                & gl%molfractions(gl%ncomp) * (gl%aij_PR(i,gl%ncomp) + gl%aij_PR(gl%ncomp,i))
End do

end subroutine da_PR_dxi
!**************************************************************************


!Andreas, November 2015
!**************************************************************************
module subroutine d2a_PR_dxidxj(gl,d2aPRdxidxj)
!**************************************************************************
!This subroutine calculates the second derivative of the mixed parameter a
!of the PR equation of state with respect to the compositions xi and xj
!
!OUTPUT: matrix d2aPRdxidxj, containing the derivatives d2a / dxidxj at const. del, tau, xk





implicit none

    type(type_gl) :: gl


double precision, dimension(30, 30):: d2aPRdxidxj
integer:: i, j

d2aPRdxidxj = 0.D0

Do i = 1, gl%ncomp - 1
    Do j = 1, gl%ncomp - 1
        d2aPRdxidxj(i,j) = gl%aij_PR(j,i) + gl%aij_PR(i,j) - gl%aij_PR(j,gl%ncomp) - gl%aij_PR(gl%ncomp,j) &
                        & + 2.D0 * gl%aij_PR(gl%ncomp,gl%ncomp) - gl%aij_PR(gl%ncomp,i) - gl%aij_PR(i,gl%ncomp)
    End do
End do

end subroutine d2a_PR_dxidxj
!**************************************************************************


!Andreas, November 2015
!**************************************************************************
module subroutine d2a_PR_dxidtau(gl,Temp,d2aPRdxidtau)
!**************************************************************************
!Second partial derivative of the mixed parameter a
!of the PR equation of state with respect to the compositions xi and tau
!
!OUTPUT: vector d2aPRdxidtau, containing the derivatives d2a / dxi dtau 





implicit none

    type(type_gl) :: gl


double precision, dimension(30):: d2aPRdxidtau
double precision:: temp
double precision, dimension(30,30)::daijPRdT   !Temperature derivative of the binary parameter aij
double precision, dimension(30):: daiPRdT      !Temperature derivative of the parameter a of fluid i
integer:: i, j
double precision:: tau, dT_dtau

d2aPRdxidtau = 0.D0

tau = gl%Tred_PR / temp
dT_dtau = - Temp / tau

!Get the temperature derivative of the parameter a for all pure fluid equations
Do i = 1, gl%ncomp
    daiPRdT(i) = -gl%ac_PR(i) *(1.D0+gl%kappa(i)*(1.D0-(Temp/gl%tc(i))**0.5d0))*gl%kappa(i)/gl%tc(i)/(Temp/gl%tc(i))**0.5d0
End do

!Get the temperature derivative of the binary parameter aij
Do i = 1, gl%ncomp
    Do j = 1, gl%ncomp
         daijPRdT(i,j) = (1.D0 - gl%kij_PR(i,j))/2.D0/(gl%ai_PR(i)*gl%ai_PR(j))**0.5d0 * (gl%ai_PR(i)*daiPRdT(j)+gl%ai_PR(j)*daiPRdT(i))
    End do
End do


Do i = 1, gl%ncomp - 1
    Do j = 1, gl%ncomp - 1        
        d2aPRdxidtau(i) = d2aPRdxidtau(i) + &
                         & gl%molfractions(j)*(daijPRdT(j,i) + daijPRdT(i,j) - daijPRdT(j,gl%ncomp) - daijPRdT(gl%ncomp,j))
    End do
    d2aPRdxidtau(i) = d2aPRdxidtau(i) - 2.D0 * gl%molfractions(gl%ncomp)*daijPRdT(gl%ncomp,gl%ncomp) + & 
                     & gl%molfractions(gl%ncomp) * (daijPRdT(i,gl%ncomp) + daijPRdT(gl%ncomp,i))
    d2aPRdxidtau(i) = d2aPRdxidtau(i) * dT_dtau
End do

end subroutine d2a_PR_dxidtau
!**************************************************************************


!Andreas, November 2015
!**************************************************************************
module subroutine db_PR_dxi(gl,dbPRdxi)
!**************************************************************************
!This subroutine calculates the partial derivative of the mixed parameter b
!of the PR equation of state with respect to the compositions xi
!
!OUTPUT: vector dbPRdxi, containing the derivatives db / dxi at const. del, tau, xk





implicit none

    type(type_gl) :: gl


double precision, dimension(30):: dbPRdxi
integer:: i, j

dbPRdxi = 0.D0

!Andreas, Feb.2016. Quadratic mixing rules for b -> mix_type = 31
if (gl%Mix_type == 3) then
    
    Do i = 1, gl%ncomp - 1
        dbPRdxi(i) = gl%bi_PR(i) - gl%bi_PR(gl%ncomp)
    End do
    
elseif (gl%Mix_type == 31) then
    
    Do i = 1, gl%ncomp - 1
        Do j = 1, gl%ncomp - 1
            dbPRdxi(i) = dbPRdxi(i) + gl%molfractions(j)*(gl%bij_PR(j,i) + gl%bij_PR(i,j) - gl%bij_PR(j,gl%ncomp) - gl%bij_PR(gl%ncomp,j))
        End do
        dbPRdxi(i) = dbPRdxi(i) - 2.D0 * gl%molfractions(gl%ncomp)*gl%bij_PR(gl%ncomp,gl%ncomp) + & 
                    & gl%molfractions(gl%ncomp) * (gl%bij_PR(i,gl%ncomp) + gl%bij_PR(gl%ncomp,i))
    End do   
    
end if

end subroutine db_PR_dxi
!**************************************************************************


!Andreas, December 2015
!**************************************************************************
module subroutine d2b_PR_dxidxj(gl,d2bPRdxidxj)
!**************************************************************************
!This subroutine calculates the second derivative of the mixed parameter b
!of the PR equation of state with respect to the compositions xi and xj
!
!OUTPUT: matrix d2bPRdxidxj, containing the derivatives d2b / dxidxj at const. del, tau, xk





implicit none

    type(type_gl) :: gl


double precision, dimension(30, 30):: d2bPRdxidxj
integer:: i, j

d2bPRdxidxj = 0.D0

!Andreas, Feb.2016. Quadratic mixing rules for b -> mix_type = 21
if (gl%Mix_type == 31) then
    Do i = 1, gl%ncomp - 1
        Do j = 1, gl%ncomp - 1
            d2bPRdxidxj(i,j) = gl%bij_PR(j,i) + gl%bij_PR(i,j) - gl%bij_PR(j,gl%ncomp) - gl%bij_PR(gl%ncomp,j) &
                            & + 2.D0 * gl%bij_PR(gl%ncomp,gl%ncomp) - gl%bij_PR(gl%ncomp,i) - gl%bij_PR(i,gl%ncomp)
        End do
    End do
end if

end subroutine d2b_PR_dxidxj
!**************************************************************************


!Andreas, November 2015
!**************************************************************************
module subroutine dar_dxi_PR(gl,Temp, dens, dardxi_PR)
!**************************************************************************
!This subroutine calculates the partial derivatives of the reduced Helmholtz 
!energy formulation of the PR equation with respect to all mole fractions xi
!at constant tau, del, and xk
!
!OUTPUT: vector dardxi_PR, containing the derivatives db / dxi at const. del, tau, xk





implicit none

    type(type_gl) :: gl


double precision:: TEMP, DENS 
double precision, dimension(30):: dardxi_PR
double precision, dimension(30)::dbPRdxi, daPRdxi
integer:: i
double precision:: delta, tau

dardxi_PR = 0.D0

delta = Dens / gl%rhored_PR
tau = gl%Tred_PR / Temp

call da_PR_dxi(gl,daPRdxi)
call db_PR_dxi(gl,dbPRdxi)

Do i = 1, gl%ncomp - 1
    dardxi_PR(i) = gl%rhored_PR*delta / (1.D0 - gl%b_PR * gl%rhored_PR * delta) * dbPRdxi(i) - &
                & ( tau / (2.D0 * sqrt_2 * gl%R_PR * gl%Tred_PR * gl%b_PR**2) * (gl%b_PR * daPRdxi(i) - gl%a_PR * dbPRdxi(i)) * &
                & dlog((1.D0 + consA_PR * gl%b_PR * gl%rhored_PR * delta)/(1.D0 - consB_PR * gl%b_PR * gl%rhored_PR * delta)) + &
                & gl%a_PR * tau / (2.D0 * sqrt_2 * gl%R_PR * gl%Tred_PR * gl%b_PR) * &
                & ((consA_PR * gl%rhored_PR * delta * dbPRdxi(i))/(1.D0 + consA_PR * gl%rhored_PR * delta * gl%b_PR) + &
                & (consB_PR * gl%rhored_PR * delta * dbPRdxi(i))/(1.D0 - consB_PR * gl%rhored_PR * delta * gl%b_PR)) )
End do

end subroutine dar_dxi_PR
!**************************************************************************

!Andreas, December 2015
!**************************************************************************
module subroutine d2ar_dxidxj_PR(gl,Temp, dens, d2ardxidxj_PR)
!**************************************************************************
!Second derivative of the reduced Helmholtz energy (using the PR!) with 
!respect to xi and xj
!
!OUTPUT: matrix d2ardxidxj_PR, containing the derivatives d2ar / dxi dxj





implicit none

    type(type_gl) :: gl


double precision:: TEMP, DENS 
double precision, dimension(30,30) :: d2ardxidxj_PR
double precision, dimension(30) :: dbPRdxi, daPRdxi
Double precision, dimension(30,30) :: d2aPRdxidxj, d2bPRdxidxj
integer:: i, j
double precision:: delta, tau

d2ardxidxj_PR = 0.D0

delta = Dens / gl%rhored_PR
tau = gl%Tred_PR / Temp

call da_PR_dxi(gl,daPRdxi)
call db_PR_dxi(gl,dbPRdxi)
call d2a_PR_dxidxj(gl,d2aPRdxidxj)
call d2b_PR_dxidxj(gl,d2bPRdxidxj)

Do i = 1, gl%ncomp - 1
    do j = 1, gl%ncomp - 1
        d2ardxidxj_PR(i,j) = gl%rhored_PR*delta / (1.D0 - gl%b_PR * gl%rhored_PR * delta) * d2bPRdxidxj(i,j) &
            & + (gl%rhored_PR*delta / (1.D0 - gl%b_PR * gl%rhored_PR * delta))**2 * dbPRdxi(i) * dbPRdxi(j) &
            & - tau / 2.D0 / sqrt_2 / gl%R_PR / gl%Tred_PR &
            & * (   ( d2aPRdxidxj(i,j) / gl%b_PR - daPRdxi(i) * dbPRdxi(j) / gl%b_PR**2 - daPRdxi(j) * dbPRdxi(i) / gl%b_PR**2 & 
            & - gl%a_PR * d2bPRdxidxj(i,j) / gl%b_PR**2 + 2.D0 * gl%a_PR * dbPRdxi(i) * dbPRdxi(j) / gl%b_PR**3 ) &
            & * dlog((1.D0 + consA_PR * gl%b_PR * gl%rhored_PR * delta)/(1.D0 - consB_PR * gl%b_PR * gl%rhored_PR * delta)) &
            & + (daPRdxi(i) / gl%b_PR - gl%a_PR * dbPRdxi(i) / gl%b_PR**2) &
            & * ((consA_PR * gl%rhored_PR * delta * dbPRdxi(j))/(1.D0 + consA_PR * gl%rhored_PR * delta * gl%b_PR) &
            & + (consB_PR * gl%rhored_PR * delta * dbPRdxi(j))/(1.D0 - consB_PR * gl%rhored_PR * delta * gl%b_PR)) &        
            & + (daPRdxi(j) / gl%b_PR - gl%a_PR * dbPRdxi(j) / gl%b_PR**2) &
            & * ((consA_PR * gl%rhored_PR * delta * dbPRdxi(i))/(1.D0 + consA_PR * gl%rhored_PR * delta * gl%b_PR) &
            & + (consB_PR * gl%rhored_PR * delta * dbPRdxi(i))/(1.D0 - consB_PR * gl%rhored_PR * delta * gl%b_PR)) &   
            & + gl%a_PR / gl%b_PR * ( ((1.D0 + consA_PR * gl%rhored_PR * delta * gl%b_PR) * consA_PR * delta * d2bPRdxidxj(i,j) & 
            & - consA_PR**2 * gl%rhored_PR**2 * delta**2 * dbPRdxi(i) * dbPRdxi(j)) &
            & / (1.D0 + consA_PR * gl%rhored_PR * delta * gl%b_PR)**2 &
            & + ((1.D0 - consB_PR * gl%rhored_PR * delta * gl%b_PR) * consB_PR * delta * d2bPRdxidxj(i,j) & 
            & + consB_PR**2 * gl%rhored_PR**2 * delta**2 * dbPRdxi(i) * dbPRdxi(j)) &
            & / (1.D0 - consB_PR * gl%rhored_PR * delta * gl%b_PR)**2 )   ) 
    End do
End do

end subroutine d2ar_dxidxj_PR
!**************************************************************************


!Andreas, December 2015
!**************************************************************************
module subroutine d2ar_dxidtau_PR(gl,Temp, dens, d2ardxidtau_PR)
!**************************************************************************
!Second derivative of the reduced Helmholtz energy (using the PR!) with 
!respect to xi and tau
!
!OUTPUT: vector d2ardxidtau_PR, containing the derivatives d2ar / d(xi)d(tau) 
!                               MULTIPLIED WITH TAU!





implicit none

    type(type_gl) :: gl


double precision:: TEMP, DENS 
double precision, dimension(30)::d2ardxidtau_PR
double precision, dimension(30)::dbPRdxi, daPRdxi, d2a_dxidtau
double precision :: daPR_dtau
integer:: i
double precision:: delta, tau

d2ardxidtau_PR = 0.D0

delta = Dens / gl%rhored_PR
tau = gl%Tred_PR / Temp

call da_PR_dxi(gl,daPRdxi)
call db_PR_dxi(gl,dbPRdxi)
call d2a_PR_dxidtau(gl,Temp, d2a_dxidtau)
daPR_dtau = da_PR_dtau(gl,Temp, 0)

Do i = 1, gl%ncomp - 1
    d2ardxidtau_PR(i) =  - 0.5D0 / sqrt_2 / gl%R_PR / gl%Tred_PR * ((daPRdxi(i) / gl%b_PR - gl%a_PR / gl%b_PR**2 * dbPRdxi(i) &
                      &  + tau * (d2a_dxidtau(i) / gl%b_PR - daPR_dtau * dbPrdxi(i) / gl%b_PR**2)) &
                      &  * dlog((1.D0 + consA_PR * gl%b_PR * gl%rhored_PR * delta)/(1.D0 - consB_PR * gl%b_PR * gl%rhored_PR * delta)) &
                      &  + (gl%a_PR / gl%b_PR + tau * daPR_dtau / gl%b_PR) & 
                      &  * ((consA_PR * gl%rhored_PR * delta * dbPRdxi(i))/(1.D0 + consA_PR * gl%rhored_PR * delta * gl%b_PR) &
                      &  + (consB_PR * gl%rhored_PR * delta * dbPRdxi(i))/(1.D0 - consB_PR * gl%rhored_PR * delta * gl%b_PR)) )
End do

d2ardxidtau_PR = d2ardxidtau_PR * tau


end subroutine d2ar_dxidtau_PR
!**************************************************************************


!Andreas, November 2015
!**************************************************************************
module subroutine d2ar_dxiddel_PR(gl,Temp, dens, d2ardxiddel_PR)
!**************************************************************************
!Second derivative of the reduced Helmholtz energy (using the PR!) with 
!respect to xi and del
!
!OUTPUT: vector d2ardxiddel_PR, containing the derivatives d2ar / d(xi)d(del) 
!                               MULTIPLIED WITH DELTA!





implicit none

    type(type_gl) :: gl


double precision:: TEMP, DENS 
double precision, dimension(30)::d2ardxiddel_PR
double precision, dimension(30)::dbPRdxi, daPRdxi
integer:: i
double precision:: delta, tau

d2ardxiddel_PR = 0.D0

delta = Dens / gl%rhored_PR
tau = gl%Tred_PR / Temp

call da_PR_dxi(gl,daPRdxi)
call db_PR_dxi(gl,dbPRdxi)

Do i = 1, gl%ncomp - 1
    !d2ardxiddel_PR(i) =  dbPRdxi(i) * rhored_PR / (1.D0-b_PR*rhored_PR*delta)**2 - tau / 2.D0 / sqrt_2 / tred_PR / R_PR / b_PR * &
    !                   & ( (b_PR * daPRdxi(i) - a_PR * dbPRdxi(i)) * (consA_PR * rhored_PR / (1.D0 + consA_PR * b_PR * rhored_PR * delta) + &
    !                   & consB_PR * rhored_PR / (1.D0 - consB_PR * b_PR * rhored_PR * delta) ) + & 
    !                   & a_PR * (consA_PR * rhored_PR * dbPRdxi(i) / (1.D0 + consA_PR * b_PR * rhored_PR * delta)**2 + &
    !                   & consB_PR * rhored_PR * dbPRdxi(i) / (1.D0 - consB_PR * b_PR * rhored_PR * delta)**2) )
    d2ardxiddel_PR(i) =  dbPRdxi(i) * gl%rhored_PR / (1.D0-gl%b_PR*gl%rhored_PR*delta)**2 - tau / 2.D0 / sqrt_2 / gl%tred_PR / gl%R_PR  * &
                       & ( daPRdxi(i) * (consA_PR * gl%rhored_PR / (1.D0 + consA_PR * gl%b_PR * gl%rhored_PR * delta) + &
                       & consB_PR * gl%rhored_PR / (1.D0 - consB_PR * gl%b_PR * gl%rhored_PR * delta) ) + & 
                       & gl%a_PR * (-consA_PR**2 * gl%rhored_PR**2 * delta * dbPRdxi(i) / (1.D0 + consA_PR * gl%b_PR * gl%rhored_PR * delta)**2 + &
                       & consB_PR**2 * gl%rhored_PR**2 * dbPRdxi(i) * delta / (1.D0 - consB_PR * gl%b_PR * gl%rhored_PR * delta)**2) )
End do

d2ardxiddel_PR = d2ardxiddel_PR * delta

end subroutine d2ar_dxiddel_PR
!**************************************************************************



!The following higher order derivatives are mainly needed for the calculation of the critical point of mixtures
!Andreas, December 2015

double precision module function d4a_SRK_dtau4(gl,Temp, nrsubst)
!**************************************************************************
!This function calculates the fourth partial derivate of the mixed 
!parameter a of SRK with respect to the inversed reduced temperature tau.
!Andreas, Dec 2015
!**************************************************************************





implicit none

    type(type_gl) :: gl


Double precision:: temp
integer :: nrsubst

d4a_SRK_dtau4 = 0.D0

end function d4a_SRK_dtau4


double precision module function d4a_PR_dtau4(gl,Temp, nrsubst)
!**************************************************************************
!This function calculates the third partial derivate of the mixed 
!parameter a of PR with respect to the inversed reduced temperature tau.
!Andreas, Dec 2015
!**************************************************************************





implicit none

    type(type_gl) :: gl


Double precision:: temp
integer :: nrsubst

d4a_PR_dtau4 = 0.D0

end function d4a_PR_dtau4


!**************************************************************************
module subroutine MIXDERIVSFNR_HIGHER_CUBIC (gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR)
!**************************************************************************
! Andreas, December 2015
!--------------------------------------------------------------------------------------------------
! SUBROUTINE FOR THE CALCULATION OF DERIVATIVES OF THE RESIDUAL PART OF THE 
! HELMHOLTZ FREE ENERGY FOR MIXTURES CALCULATED WITH CUBIC EOS (SRK AND PR)
!--------------------------------------------------------------------------------------------------
! INPUT PARAMETERS:     
! TEMPERATURE - T   K
! DENSITY     - D   mol/m^3
! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", 
!               INDICATING WHICH DERIVATIVES ARE NEEDED:
!                1. NORMALIZED RESIDUAL MIXTURE HELMHOLTZ ENERGY F AS A FUNCTION OF DEL AND TAU AND X)
!                2. 1ST DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
!                3. 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
!                4. 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
!                5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
!                6: 1ST MIXED DERIVATIVE OF F WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
!                7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
!                8: 3RD DERIVATIVE OF F WITH RESPECT TO D, MULTIPLIED BY DEL^3
!                9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
!               10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
!               11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^4, MULTIPLIED BY DEL^4
!               12: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^3 and TAU, MULTIPLIED BY DEL^3*TAU
!               13: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^2 and TAU^2, MULTIPLIED BY DEL^2*TAU^2
!               14: 4TH DERIVATIVE OF F WITH RESPECT TO DEL and TAU^3, MULTIPLIED BY DEL*TAU^3
!               15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU^4, MULTIPLIED BY TAU^4
! 
! OUTPUT PARAMETERS: 
! MIXDERFNR   - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
!                AS INDECATED IN "GETDER"
!--------------------------------------------------------------------------------------------------








implicit none

    type(type_gl) :: gl


double precision:: TEMPERATURE, DENSITY  
integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed 
double precision, dimension(15)::MIXDERFNR  ! array with the computed values for the derivatives

double precision:: tau, delta

double precision:: D1_cubic, D2_cubic

!Supporting variables
double precision:: psi_m, psi_p
double precision:: dpsi_m_ddel, d2psi_m_ddel2, d3psi_m_ddel3, d4psi_m_ddel4 
double precision:: dpsi_p_ddel, d2psi_p_ddel2, d3psi_p_ddel3, d4psi_p_ddel4

Double precision:: Pi_12, dPi12_ddel, d2Pi12_ddel2, d3Pi12_ddel3
Double precision:: Pi_12_2, Pi_12_3, Pi_12_4

Double precision:: da_dtau, d2a_dtau2, d3a_dtau3, d4a_dtau4
Double precision, dimension(15):: dnadtaun
double precision, dimension(15,30)::dnadxidtaum       
double precision, dimension(15,30,30)::dnadxi2dtaum       
integer :: get_a_deriv

Double precision:: dtaua_dtau, d2taua_dtau2, d3taua_dtau3, d4taua_dtau4

Double precision:: a_cubic, b_cubic, R_cubic
Double precision:: rhored_cubic, tred_cubic

Double precision:: tau2, tau3, tau4, del2, del3, del4        !Powers of delta and tau

integer:: i,j

!All CUBIC equations can all be written in the same form:
!See Bell et al. (2016), "Helmholtz energy translations for common cubic equations of state for use in one-fluid and multi-fluid mixture models", to be submitted to ??
!$$ p = \frac{R*T}{v-b} - \frac{a}{(v+D1*b)*(v+D2*b)} $$
!For the different cubics, it is:
!   vdW: D1 = 0,        D2 = 0
!   SRK: D1 = 1,        D2 = 0
!   PR:  D1 = 1+2**0.5  D2 = 1 - 2**0.5
!The Helmholtz translation of the cubics can then also be written in a universal way:
!$$ \alpha_r = \psi^(-) - \frac{\tau*a}{R*T_r} * \psi^(+) $$
!Where:
!$$ \psi^(-) = -ln(1-b * \delta * \rho_r) $$
!$$ \psi^(+) = \frac{ln(\frac{D1*b*\rho_r*\delta}{D2*b*\rho_r*\delta})}{b*(D1-D2)} 
!The variable \psi^(-) is called psi_m in the code, this variable is only a function of delta
!The variable \psi^(+) is called psi_p in the code, this variable is only a function of delta

If ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then
    !SRK
    D1_cubic = 1.D0
    D2_cubic = 0.D0
    rhored_cubic = gl%rhored_SRK
    tred_cubic = gl%tred_SRK
    a_cubic = gl%a_SRK
    b_cubic = gl%b_SRK
    R_cubic = gl%R_SRK
elseif ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then
    !PR
    D1_cubic = consA_PR
    D2_cubic = -consB_PR    !THE MINUS IS VERY IMPORTANT!! NOTE THE DIFFERENT DEFINITIONS OF THE CONSTANTS IN THE PAPER AND IN THE CODE!!!
    rhored_cubic = gl%rhored_PR
    tred_cubic = gl%tred_PR
    a_cubic = gl%a_PR
    b_cubic = gl%b_PR
    R_cubic = gl%R_PR
end if

delta = Density / rhored_cubic
tau = Tred_cubic / Temperature

!Variables for speeding up computations
tau2 = tau**2
tau3 = tau*tau2
tau4 = tau2*tau2
del2 = delta**2
del3 = del2 * delta
del4 = del2 * del2

!do i = 1, ncomp
!    sqrt_tr_tci(i) = (Tred_cubic / tc(i))**0.5
!end do

!Calculate the derivatives of psi^(-) wrt delta
psi_m = -dlog(1.D0-b_cubic*delta*rhored_cubic)
dpsi_m_ddel = b_cubic * rhored_cubic / (1.D0 - b_cubic * delta * rhored_cubic)
d2psi_m_ddel2 = dpsi_m_ddel**2
d3psi_m_ddel3 = 2.D0 * d2psi_m_ddel2 * dpsi_m_ddel  
d4psi_m_ddel4 = 6.D0 * d2psi_m_ddel2**2


!Calculate the derivatives of psi^(+) wrt delta
!First calculate supporting variable Pi_12 and its derivatives wrt delta
Pi_12 = (1.D0 + D1_cubic * b_cubic * rhored_cubic * delta) * (1.D0 + D2_cubic * b_cubic * rhored_cubic * delta)
Pi_12_2 = Pi_12 * Pi_12
Pi_12_3 = Pi_12 * Pi_12_2
Pi_12_4 = Pi_12_2 * Pi_12_2
dPi12_ddel = b_cubic * rhored_cubic * (2.D0 * D1_cubic * D2_cubic * b_cubic * delta * rhored_cubic + D1_cubic + D2_cubic)
d2Pi12_ddel2 = 2.D0 * D1_cubic * D2_cubic * b_cubic**2 * rhored_cubic**2
d3Pi12_ddel3 = 0.D0

psi_p = dlog((D1_cubic*b_cubic*rhored_cubic*delta+1.D0)/(D2_cubic*b_cubic*rhored_cubic*delta+1.D0))/(b_cubic*(D1_cubic - D2_cubic))
dpsi_p_ddel = rhored_cubic / Pi_12
d2psi_p_ddel2 = -rhored_cubic * dPi12_ddel / Pi_12_2
d3psi_p_ddel3 = rhored_cubic * (-Pi_12 * d2Pi12_ddel2 + 2.D0 * dPi12_ddel**2) / Pi_12_3
d4psi_p_ddel4 = rhored_cubic * (-Pi_12_2 * d3Pi12_ddel3 + 6.D0 * Pi_12 * dPi12_ddel * d2Pi12_ddel2 - 6.D0 * dPi12_ddel**3) / Pi_12_4


! Get the derivatives of the parameter a with respect to tau
get_a_deriv = 4 !Get all derivatives of a with respect to tau
call da_dxi_dtau_cubic_all(gl,Temperature, density, get_a_deriv, dnadtaun, dnadxidtaum, dnadxi2dtaum)
da_dtau = dnadtaun(1)
d2a_dtau2 = dnadtaun(2)
d3a_dtau3 = dnadtaun(3)
d4a_dtau4 = dnadtaun(4)

!Calculate the derivatives of the product of (a * tau) with respect to tau
dtaua_dtau = tau * da_dtau + a_cubic
d2taua_dtau2 = tau * d2a_dtau2 + 2.D0 * da_dtau
d3taua_dtau3 = tau * d3a_dtau3 + 3.D0 * d2a_dtau2
d4taua_dtau4 = tau * d4a_dtau4 + 4.D0 * d3a_dtau3


!alphar
if (GETDER(1) == 1) then
    MIXDERFNR(1) = psi_m - tau * a_cubic / R_cubic / Tred_cubic * psi_p
end if

!d_alphar_d_delta * delta
if (GETDER(2) == 1) then
    MIXDERFNR(2) = dpsi_m_ddel - tau * a_cubic / R_cubic / Tred_cubic * dpsi_p_ddel
    MIXDERFNR(2) = MIXDERFNR(2) * delta
end if

!d2_alphar_d_delta2 * delta^2
if (GETDER(3) == 1) then
    MIXDERFNR(3) = d2psi_m_ddel2 - tau * a_cubic / R_cubic / Tred_cubic * d2psi_p_ddel2
    MIXDERFNR(3) = MIXDERFNR(3) * del2
end if

!d_alphar_d_tau * tau
if (GETDER(4) == 1) then
    MIXDERFNR(4) = - dtaua_dtau / R_cubic / Tred_cubic * psi_p
    MIXDERFNR(4) = MIXDERFNR(4) * tau
end if

!d2_alphar_d_tau2 * tau^2
if (GETDER(5) == 1) then
    MIXDERFNR(5) = - d2taua_dtau2 / R_cubic / Tred_cubic * psi_p
    MIXDERFNR(5) = MIXDERFNR(5) * tau2
end if

!d2_alphar_d_delta_d_tau * delta * tau
if (GETDER(6) == 1) then
    MIXDERFNR(6) = - dtaua_dtau / R_cubic / Tred_cubic * dpsi_p_ddel
    MIXDERFNR(6) = MIXDERFNR(6) * tau * delta
end if

!d3_alphar_d_delta_d_tau2 * delta * tau^2
if (GETDER(7) == 1) then
    MIXDERFNR(7) = - d2taua_dtau2 / R_cubic / Tred_cubic * dpsi_p_ddel
    MIXDERFNR(7) = MIXDERFNR(7) * tau2 * delta
end if

!d3_alphar_d_delta3 * delta^3
if (GETDER(8) == 1) then
    MIXDERFNR(8) = d3psi_m_ddel3 - tau * a_cubic / R_cubic / Tred_cubic * d3psi_p_ddel3
    MIXDERFNR(8) = MIXDERFNR(8) * del3
end if

!d3_alphar_d_tau3 * tau^3
if (GETDER(9) == 1) then
    MIXDERFNR(9) = - d3taua_dtau3 / R_cubic / Tred_cubic * psi_p
    MIXDERFNR(9) = MIXDERFNR(9) * tau3
end if

!d3_alphar_d_delta2_dtau * delta^2 * tau
if (GETDER(10) == 1) then
    MIXDERFNR(10) = - dtaua_dtau / R_cubic / Tred_cubic * d2psi_p_ddel2
    MIXDERFNR(10) = MIXDERFNR(10) * del2 * tau
end if

!d4_alphar_d_delta4 * delta^4
if (GETDER(11) == 1) then
    MIXDERFNR(11) = d4psi_m_ddel4 - tau * a_cubic / R_cubic / Tred_cubic * d4psi_p_ddel4
    MIXDERFNR(11) = MIXDERFNR(11) * del4
end if

!d4_alphar_d_delta3_dtau * delta^3*tau
if (GETDER(12) == 1) then
    MIXDERFNR(12) = - dtaua_dtau / R_cubic / Tred_cubic * d3psi_p_ddel3
    MIXDERFNR(12) = MIXDERFNR(12) * del3 * tau    
end if

!d4_alphar_d_delta2_dtau2 * delta^2*tau^2
if (GETDER(13) == 1) then
    MIXDERFNR(13) = - d2taua_dtau2 / R_cubic / Tred_cubic * d2psi_p_ddel2
    MIXDERFNR(13) = MIXDERFNR(13) * del2 * tau2
end if

!d4_alphar_d_delta_dtau3 * delta*tau^3
if (GETDER(14) == 1) then
    MIXDERFNR(14) = - d3taua_dtau3 / R_cubic / Tred_cubic * dpsi_p_ddel
    MIXDERFNR(14) = MIXDERFNR(14) * delta * tau3
end if

!d4_alphar_d_tau4 * tau^4
if (GETDER(15) == 1) then
    MIXDERFNR(15) = - d4taua_dtau4 / R_cubic / Tred_cubic * psi_p
    MIXDERFNR(15) = MIXDERFNR(15) * tau4
end if


End subroutine MIXDERIVSFNR_HIGHER_CUBIC





!Andreas, December 2015
!**************************************************************************
module subroutine da_dxi_dtau_cubic_all (gl,Temp, dens, get_a_deriv, dnadtaun, dnadxidtaum, dnadxi2dtaum)
!This subroutine calculates the nth derivatives of the parameter a with 
!respect to tau
!INPUT: get_a_deriv: Integer that indicates up to which derivatives are calculated
!                    get_a_deriv = 1 --> calculate da_dtau (and d2a_dxidtau and d3a_dxi2dtau)
!                    get_a_deriv = 2 --> calculate da_dtau and d2a_dtau2 (and d3a_dxidtau2 and d4a_dxi2dtau2)
!                    get_a_deriv = 3 --> calculate da_dtau, d2a_dtau2 and d3a_dtau3 (and d4a_dxidtau3 and d5a_dxi2dtau3)
!                    get_a_deriv = 4 --> calculate up to d4a_dtau4 (and d5a_dxidtau4 and d6a_dxi2dtau4)

!OUTPUT:    dnadtaun(1) : da_dtau
!           dnadtaun(2) : d2a_dtau2
!           dnadtaun(3) : d3a_dtau3
!           dnadtaun(4) : d4a_dtau4

!           dnadxidtaum(1,i) : da_dxi
!           dnadxidtaum(2,i) : d2a_dxidtau
!           dnadxidtaum(3,i) : d3a_dxidtau2
!           dnadxidtaum(4,i) : d4a_dxidtau3
!           dnadxidtaum(5,i) : d5a_dxidtau4

!           dnadxidtaum(1,i,j) : da_dxidxj
!           dnadxidtaum(2,i,j) : d2a_dxidxjdtau
!           dnadxidtaum(3,i,j) : d3a_dxidxjdtau2
!           dnadxidtaum(4,i,j) : d4a_dxidxjdtau3
!           dnadxidtaum(5,i,j) : d5a_dxidxjdtau4











implicit none

    type(type_gl) :: gl


double precision:: TEMP, DENS
integer get_a_deriv                             ! array specifier to indicate, which derivatives are needed 
double precision, dimension(15)::dnadtaun       ! array with the computed values for the derivatives
double precision, dimension(15,30)::dnadxidtaum       ! array with the computed values for the derivatives
double precision, dimension(15,30,30) :: dnadxi2dtaum       ! array with the computed values for the derivatives

double precision:: tau, delta

Double precision:: a_cubic, b_cubic, R_cubic
Double precision:: rhored_cubic, tred_cubic
Double precision, dimension(30):: ai_cubic, dai_dtau, d2ai_dtau2, d3ai_dtau3, d4ai_dtau4, bi_cubic
double precision, dimension(30):: ac_cubic
Double precision, dimension(30):: mi_cubic
double precision, dimension(:,:), allocatable:: aij_cubic, daij_dtau, d2aij_dtau2, d3aij_dtau3, d4aij_dtau4
double precision, dimension(:,:), allocatable:: kij_cubic
double precision, dimension(:,:), allocatable:: term_uij_cubic, duij_dtau, d2uij_dtau2, d3uij_dtau3, d4uij_dtau4
Double precision, dimension(30):: term_B_cubic

Double precision:: tau2, tau3, tau4, tau5, del2, del3, del4, tau05        !Powers of delta and tau, tau05 = sqrt(tau)
Double precision, dimension(30):: sqrt_tr_tci
double precision, dimension(:,:), allocatable:: term_uij_cubic05   !Square root of term uij

!New variables for terms of Mathias and Copeman temperature dependence of ai, Andreas Jäger, January 2017
Double precision, dimension(30) :: term_D_cubic 
Double precision, dimension(30) :: dterm_D_cubic_dtau, d2term_D_cubic_dtau2, d3term_D_cubic_dtau3, d4term_D_cubic_dtau4
Double precision, dimension(30) :: dterm_B_cubic_dtau, d2term_B_cubic_dtau2, d3term_B_cubic_dtau3, d4term_B_cubic_dtau4

integer:: i, j, k

!Variables for numerical and analytical derivatives for PSRK, Andreas Jäger, January 2017
double precision :: a_PSRK_orig, a_PSRK_p, a_PSRK_m
double precision :: a_PSRK_iptaup, a_PSRK_imtaup, a_PSRK_iptaum, a_PSRK_imtaum
double precision :: a_PSRK_ipjp, a_PSRK_imjp, a_PSRK_ipjm, a_PSRK_imjm
double precision :: tau_p, tau_m, Temp_p, Temp_m, Temp_pm
double precision, dimension(30) :: molfractions_orig
double precision :: delta_x, delta_tau
integer, dimension(nderivs):: GETDER                !array specifier to indicate, which derivative is needed 

integer:: C_or_R                !Specify which part of UNIFAC is needed (0:both, 1:only combinatorial, 2: only residual)
integer:: errval
double precision, dimension(30) :: db_dxi

if (.not. allocated(aij_cubic))         allocate(aij_cubic(30,30))                        ! allocate(aij_cubic(30,30))
if (.not. allocated(daij_dtau))         allocate(daij_dtau(30,30))            ! allocate(daij_dtau(gl%ncomp,gl%ncomp))
if (.not. allocated(d2aij_dtau2))       allocate(d2aij_dtau2(30,30))          ! allocate(d2aij_dtau2(gl%ncomp,gl%ncomp))
if (.not. allocated(d3aij_dtau3))       allocate(d3aij_dtau3(30,30))          ! allocate(d3aij_dtau3(gl%ncomp,gl%ncomp))
if (.not. allocated(d4aij_dtau4))       allocate(d4aij_dtau4(30,30))          ! allocate(d4aij_dtau4(gl%ncomp,gl%ncomp))
if (.not. allocated(kij_cubic))         allocate(kij_cubic(30,30))                        ! allocate(kij_cubic(30,30))
if (.not. allocated(term_uij_cubic))    allocate(term_uij_cubic(30,30))       ! allocate(term_uij_cubic(gl%ncomp,gl%ncomp))
if (.not. allocated(duij_dtau))         allocate(duij_dtau(30,30))            ! allocate(duij_dtau(gl%ncomp,gl%ncomp))
if (.not. allocated(d2uij_dtau2))       allocate(d2uij_dtau2(30,30))          ! allocate(d2uij_dtau2(gl%ncomp,gl%ncomp))
if (.not. allocated(d3uij_dtau3))       allocate(d3uij_dtau3(30,30))          ! allocate(d3uij_dtau3(gl%ncomp,gl%ncomp))
if (.not. allocated(d4uij_dtau4))       allocate(d4uij_dtau4(30,30))          ! allocate(d4uij_dtau4(gl%ncomp,gl%ncomp))
if (.not. allocated(term_uij_cubic05))  allocate(term_uij_cubic05(30,30))     ! allocate(term_uij_cubic05(gl%ncomp,gl%ncomp))

If ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then
    !SRK
    rhored_cubic = gl%rhored_SRK
    tred_cubic = gl%tred_SRK
    a_cubic = gl%a_SRK
    b_cubic = gl%b_SRK
    R_cubic = gl%R_SRK
    ai_cubic = gl%ai_SRK   !ac*(1+m*(1-(T/Tc)**0.5))**2
    bi_cubic = gl%bi_SRK
    ac_cubic = gl%ac_SRK   
    aij_cubic = gl%aij_SRK
    kij_cubic = gl%kij_SRK
    mi_cubic = gl%m_SRK
elseif ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then
    !PR
    rhored_cubic = gl%rhored_PR
    tred_cubic = gl%tred_PR
    a_cubic = gl%a_PR
    b_cubic = gl%b_PR
    bi_cubic = gl%bi_PR
    R_cubic = gl%R_PR
    ac_cubic = gl%ac_PR
    ai_cubic = gl%ai_PR   !ac*(1+m*(1-(T/Tc)**0.5))**2
    aij_cubic = gl%aij_PR
    kij_cubic = gl%kij_PR
    mi_cubic = gl%kappa
end if

delta = Dens / rhored_cubic
tau = Tred_cubic / Temp


!If (Mix_type /= 22) then
    
!Variables for speeding up computations
tau05 = tau**0.5
tau2 = tau**2
tau3 = tau*tau2
tau4 = tau2*tau2
tau5 = tau4*tau

del2 = delta**2
del3 = del2 * delta
del4 = del2 * del2

do i = 1, gl%ncomp
    sqrt_tr_tci(i) = (Tred_cubic / gl%tc(i))**0.5
    !New variable needed for Mathias and Copeman function
    term_D_cubic(i) = 1.D0 - sqrt_tr_tci(i) / tau05
    dterm_D_cubic_dtau(i) = 0.5D0 * sqrt_tr_tci(i) / (tau * tau05)
    if (get_a_deriv > 1) then
        d2term_D_cubic_dtau2(i) = -0.75D0 * sqrt_tr_tci(i) / (tau2 * tau05)
    end if
    if (get_a_deriv > 2) then
        d3term_D_cubic_dtau3(i) = 1.875D0 * sqrt_tr_tci(i) / (tau3 * tau05)
    end if
    if (get_a_deriv > 3) then
        d4term_D_cubic_dtau4(i) = -6.5625D0 * sqrt_tr_tci(i) / (tau4 * tau05)
    end if
end do


!--------------
!Calculate tau derivatives of the mixed parameter a_cubic of the cubic equation
!Calculate the supporting term Bi and derivatives of aii wrt tau
!OLD VERSION, BEFORE MATHIAS COPEMAN CHANGES, Andreas Jäger, January 2017
!do i = 1, ncomp
!    term_B_cubic(i) = 1.D0 + mi_cubic(i)*(1.D0 - sqrt_tr_tci(i) / tau05)
!    dai_dtau(i) = ac_cubic(i) * mi_cubic(i) * term_B_cubic(i) * sqrt_tr_tci(i) / (tau * tau05)    
!    if (get_a_deriv > 1) then
!    d2ai_dtau2(i) = 0.5D0 * ac_cubic(i) * mi_cubic(i) * (mi_cubic(i) * sqrt_tr_tci(i)**2 / tau3 & 
!                    & - 3.D0 * term_B_cubic(i) * sqrt_tr_tci(i) / (tau2 * tau05))
!    end if
!    if (get_a_deriv > 2) then
!    d3ai_dtau3(i) = 0.75D0 * ac_cubic(i) * mi_cubic(i) * (-3.D0 * mi_cubic(i) * sqrt_tr_tci(i)**2 / tau4 & 
!                    & + 5.D0 * term_B_cubic(i) * sqrt_tr_tci(i) / (tau3 * tau05))
!    end if
!    if (get_a_deriv > 3) then
!    d4ai_dtau4(i) = 0.375D0 * ac_cubic(i) * mi_cubic(i) * (29.D0 * mi_cubic(i) * sqrt_tr_tci(i)**2 / tau5 & 
!                    & - 35.D0 * term_B_cubic(i) * sqrt_tr_tci(i) / (tau4 * tau05))
!    end if
!end do

!NEW VERSION WITH MATHIAS COPEMAN CHANGES
!Calculate the help variable B and its tau derivatives
!Contrary to the paper of Bell and Jäger (2016) (JRES 121, 238 - 263, 2016) the formulation with the sum over n is not used here in order to avoid problems like 0^0, especially when calculating pure substances with SRK or PR at the critical temperature 
do i = 1, gl%ncomp
    term_B_cubic(i) = 1.D0 + gl%Cji_cubic(1,i) * term_D_cubic(i) + gl%Cji_cubic(2,i) * term_D_cubic(i)**2 &
                  & + gl%Cji_cubic(3,i) * term_D_cubic(i)**3 
       
    dterm_B_cubic_dtau(i) =  dterm_D_cubic_dtau(i) * (1.D0 * gl%Cji_cubic(1,i) &
                        & + 2.D0 * gl%Cji_cubic(2,i) * term_D_cubic(i) &
                        & + 3.D0 * gl%Cji_cubic(3,i) * term_D_cubic(i)**2) 

    dai_dtau(i) = 2.D0 * ac_cubic(i) * term_B_cubic(i) * dterm_B_cubic_dtau(i)  
    
    if (get_a_deriv > 1) then
        d2term_B_cubic_dtau2(i) = gl%Cji_cubic(1,i) * d2term_D_cubic_dtau2(i) &
                              & + 2.D0 * gl%Cji_cubic(2,i) * (dterm_D_cubic_dtau(i)**2 &
                              & + term_D_cubic(i) * d2term_D_cubic_dtau2(i)) &
                              & + 3.D0 * gl%Cji_cubic(3,i) * (2.D0 * dterm_D_cubic_dtau(i)**2 &
                              & + term_D_cubic(i) * d2term_D_cubic_dtau2(i)) * term_D_cubic(i)
        
        d2ai_dtau2(i) = 2.D0 * ac_cubic(i) * (term_B_cubic(i) * d2term_B_cubic_dtau2(i) + dterm_B_cubic_dtau(i)**2)     
    end if   
    
    if (get_a_deriv > 2) then
        d3term_B_cubic_dtau3(i) = gl%Cji_cubic(1,i) * d3term_D_cubic_dtau3(i) &
                              & + 2.D0 * gl%Cji_cubic(2,i) * (3.D0 * dterm_D_cubic_dtau(i) * d2term_D_cubic_dtau2(i) &
                              & + term_D_cubic(i) * d3term_D_cubic_dtau3(i)) &
                              & + 3.D0 * gl%Cji_cubic(3,i) * (6.D0 * term_D_cubic(i) * dterm_D_cubic_dtau(i) &
                              & * d2term_D_cubic_dtau2(i) + 2.D0 * dterm_D_cubic_dtau(i)**3 & 
                              & + term_D_cubic(i)**2 * d3term_D_cubic_dtau3(i))
        
        d3ai_dtau3(i) = 2.D0 * ac_cubic(i) * (term_B_cubic(i) * d3term_B_cubic_dtau3(i) &
                    & + 3.D0 * dterm_B_cubic_dtau(i) * d2term_B_cubic_dtau2(i))  
    end if   
    
    if (get_a_deriv > 3) then
        d4term_B_cubic_dtau4(i) = gl%Cji_cubic(1,i) * d4term_D_cubic_dtau4(i) &
                              & + 2.D0 * gl%Cji_cubic(2,i) * (4.D0 * dterm_D_cubic_dtau(i) * d3term_D_cubic_dtau3(i) & 
                              & + 3.D0 * d2term_D_cubic_dtau2(i)**2 + term_D_cubic(i) * d4term_D_cubic_dtau4(i)) &
                              & + 3.D0 * gl%Cji_cubic(3,i) * (12.D0 * dterm_D_cubic_dtau(i)**2 * d2term_D_cubic_dtau2(i) & 
                              & + (8.D0 * dterm_D_cubic_dtau(i) * d3term_D_cubic_dtau3(i) & 
                              & + 6.D0 * d2term_D_cubic_dtau2(i)**2) * term_D_cubic(i) & 
                              & + term_D_cubic(i)**2 * d4term_D_cubic_dtau4(i))
        d4ai_dtau4(i) = 2.D0 * ac_cubic(i) * (term_B_cubic(i) * d4term_B_cubic_dtau4(i) &
                      & + 4.D0 * dterm_B_cubic_dtau(i) * d3term_B_cubic_dtau3(i) + 3.D0 * d2term_B_cubic_dtau2(i)**2)  
    end if
    
end do

If (gl%Mix_type /= 22) then        !In this part, derivatives of a with respect to tau and x for quadratic mixing rules for a are calulated
    
!Calculate uij and its derivatives
do i = 1, gl%ncomp
    do j= 1, gl%ncomp
        if (i <= j) then     !matrix is symmetric, so only one half and the diagonal must be calculated
            term_uij_cubic(i,j) = ai_cubic(i) * ai_cubic(j)
            term_uij_cubic05(i,j) = term_uij_cubic(i,j)**0.5
            duij_dtau(i,j) = ai_cubic(i) * dai_dtau(j) + ai_cubic(j) * dai_dtau(i)
            if (get_a_deriv > 1) then
                d2uij_dtau2(i,j) = ai_cubic(i) * d2ai_dtau2(j) + ai_cubic(j) * d2ai_dtau2(i) + &
                                &   2.D0 * dai_dtau(i) * dai_dtau(j)
            end if
            if (get_a_deriv > 2) then
                d3uij_dtau3(i,j) = ai_cubic(i) * d3ai_dtau3(j) + ai_cubic(j) * d3ai_dtau3(i) & 
                            & + 3.D0 * dai_dtau(i) * d2ai_dtau2(j) + 3.D0 * dai_dtau(j) * d2ai_dtau2(i)
            end if
            if (get_a_deriv > 3) then
                d4uij_dtau4(i,j) = ai_cubic(i) * d4ai_dtau4(j) + ai_cubic(j) * d4ai_dtau4(i) & 
                            & + 4.D0 * dai_dtau(i) * d3ai_dtau3(j) + 4.D0 * dai_dtau(j) * d3ai_dtau3(i) &
                            & + 6.D0 * d2ai_dtau2(i) * d2ai_dtau2(j)           
            end if
            daij_dtau(i,j) = 0.5D0 * (1.D0 - kij_cubic(i,j)) / term_uij_cubic05(i,j) * duij_dtau(i,j)       
            if (get_a_deriv > 1) then
                d2aij_dtau2(i,j) = 0.25D0 * (1.D0 - kij_cubic(i,j)) / (term_uij_cubic05(i,j) * term_uij_cubic(i,j)) * &
                             & (2.D0 * term_uij_cubic(i,j) * d2uij_dtau2(i,j) - duij_dtau(i,j)**2) 
            end if
            if (get_a_deriv > 2) then
                d3aij_dtau3(i,j) = 0.125D0 * (1.D0 - kij_cubic(i,j)) / (term_uij_cubic05(i,j) * term_uij_cubic(i,j)**2) * &
                             & (4.D0 * term_uij_cubic(i,j)**2 * d3uij_dtau3(i,j) - & 
                             & 6.D0 * term_uij_cubic(i,j) * duij_dtau(i,j) * d2uij_dtau2(i,j) + 3.D0 * duij_dtau(i,j)**3)     
            end if
            if (get_a_deriv > 3) then
                d4aij_dtau4(i,j) = 0.0625D0 * (1.D0 - kij_cubic(i,j)) / (term_uij_cubic05(i,j) * term_uij_cubic(i,j)**3) * &
                             & (-4.D0 * term_uij_cubic(i,j)**2 * &
                             & (4.D0 * duij_dtau(i,j) * d3uij_dtau3(i,j) + 3.D0 * d2uij_dtau2(i,j)**2) + &
                             & 8.D0 * term_uij_cubic(i,j)**3 * d4uij_dtau4(i,j) & 
                             & + 36.D0 * term_uij_cubic(i,j) * duij_dtau(i,j)**2 * d2uij_dtau2(i,j) &
                             & - 15.D0 *duij_dtau(i,j)**4)
            end if
            
        else            !This is for speeding up the code. Get elements from other half of the matrix       
            daij_dtau(i,j) = daij_dtau(j,i)        
            if (get_a_deriv > 1) then
                d2aij_dtau2(i,j) = d2aij_dtau2(j,i)
            end if
            if (get_a_deriv > 2) then
                d3aij_dtau3(i,j) = d3aij_dtau3(j,i)
            end if
            if (get_a_deriv > 3) then
                d4aij_dtau4(i,j) = d4aij_dtau4(j,i)
            end if
        end if
    end do
end do

!Get the tau derivatives of a_cubic
dnadtaun = 0.D0
do i = 1, gl%ncomp
    do j= i, gl%ncomp  !matrix is symmetric, so only one half and the diagonal must be calculated
        if (i < j) then     
            dnadtaun(1) = dnadtaun(1) + 2.D0 * gl%molfractions(i) * gl%molfractions(j) * daij_dtau(i,j)
            if (get_a_deriv > 1) then
                dnadtaun(2) = dnadtaun(2) + 2.D0 * gl%molfractions(i) * gl%molfractions(j) * d2aij_dtau2(i,j)
            end if
            if (get_a_deriv > 2) then
                dnadtaun(3) = dnadtaun(3) + 2.D0 * gl%molfractions(i) * gl%molfractions(j) * d3aij_dtau3(i,j)
            end if
            if (get_a_deriv > 3) then
                dnadtaun(4) = dnadtaun(4) + 2.D0 * gl%molfractions(i) * gl%molfractions(j) * d4aij_dtau4(i,j)
            end if
        end if
        if (i == j) then
            dnadtaun(1) = dnadtaun(1) +  gl%molfractions(i) * gl%molfractions(j) * daij_dtau(i,j)
            if (get_a_deriv > 1) then
                dnadtaun(2) = dnadtaun(2) + gl%molfractions(i) * gl%molfractions(j) * d2aij_dtau2(i,j)
            end if
            if (get_a_deriv > 2) then
                dnadtaun(3) = dnadtaun(3) + gl%molfractions(i) * gl%molfractions(j) * d3aij_dtau3(i,j)
            end if
            if (get_a_deriv > 3) then
                dnadtaun(4) = dnadtaun(4) + gl%molfractions(i) * gl%molfractions(j) * d4aij_dtau4(i,j)
            end if
        end if
    end do
end do


!Andreas, December 2015
!Derivatives of a with respect to xi and tau
!----------------------------------------------------------------------------------------------------------------------------------------
dnadxidtaum = 0.D0
Do i = 1, gl%ncomp - 1
    Do j = 1, gl%ncomp - 1
        !*** da_dxi
        dnadxidtaum(1,i) = dnadxidtaum(1,i) + gl%molfractions(j)*(aij_cubic(j,i) + aij_cubic(i,j) - aij_cubic(j,gl%ncomp) - aij_cubic(gl%ncomp,j))
        !*** d2a_dxidtau
        dnadxidtaum(2,i) = dnadxidtaum(2,i) + gl%molfractions(j)*(daij_dtau(j,i) + daij_dtau(i,j) - daij_dtau(j,gl%ncomp) - daij_dtau(gl%ncomp,j))
        !*** d3a_dxidtau2
        if (get_a_deriv > 1) then    
            dnadxidtaum(3,i) = dnadxidtaum(3,i) + gl%molfractions(j)*(d2aij_dtau2(j,i) + d2aij_dtau2(i,j) - d2aij_dtau2(j,gl%ncomp) - d2aij_dtau2(gl%ncomp,j))
        end if
        !*** d4a_dxidtau3
        if (get_a_deriv > 2) then
            dnadxidtaum(4,i) = dnadxidtaum(4,i) + gl%molfractions(j)*(d3aij_dtau3(j,i) + d3aij_dtau3(i,j) - d3aij_dtau3(j,gl%ncomp) - d3aij_dtau3(gl%ncomp,j))   
        end if
        !*** d4a_dxidtau4
        if (get_a_deriv > 3) then
            dnadxidtaum(5,i) = dnadxidtaum(5,i) + gl%molfractions(j)*(d4aij_dtau4(j,i) + d4aij_dtau4(i,j) - d4aij_dtau4(j,gl%ncomp) - d4aij_dtau4(gl%ncomp,j))    
        end if
    End do
    !*** da_dxi
    dnadxidtaum(1,i) = dnadxidtaum(1,i) - 2.D0 * gl%molfractions(gl%ncomp)*aij_cubic(gl%ncomp,gl%ncomp) + & 
                & gl%molfractions(gl%ncomp) * (aij_cubic(i,gl%ncomp) + aij_cubic(gl%ncomp,i))
    !*** d2a_dxidtau
    dnadxidtaum(2,i) = dnadxidtaum(2,i) - 2.D0 * gl%molfractions(gl%ncomp)*daij_dtau(gl%ncomp,gl%ncomp) + & 
            & gl%molfractions(gl%ncomp) * (daij_dtau(i,gl%ncomp) + daij_dtau(gl%ncomp,i))
    if (get_a_deriv > 1) then
        !*** d3a_dxidtau2       
        dnadxidtaum(3,i) = dnadxidtaum(3,i) - 2.D0 * gl%molfractions(gl%ncomp)*d2aij_dtau2(gl%ncomp,gl%ncomp) + & 
                    & gl%molfractions(gl%ncomp) * (d2aij_dtau2(i,gl%ncomp) + d2aij_dtau2(gl%ncomp,i))    
    end if
    !*** d4a_dxidtau3
    if (get_a_deriv > 2) then
        dnadxidtaum(4,i) = dnadxidtaum(4,i) - 2.D0 * gl%molfractions(gl%ncomp)*d3aij_dtau3(gl%ncomp,gl%ncomp) + & 
                & gl%molfractions(gl%ncomp) * (d3aij_dtau3(i,gl%ncomp) + d3aij_dtau3(gl%ncomp,i))        
    end if
    !*** d4a_dxidtau4
    if (get_a_deriv > 3) then
        dnadxidtaum(5,i) = dnadxidtaum(5,i) - 2.D0 * gl%molfractions(gl%ncomp)*d4aij_dtau4(gl%ncomp,gl%ncomp) + & 
                    & gl%molfractions(gl%ncomp) * (d4aij_dtau4(i,gl%ncomp) + d4aij_dtau4(gl%ncomp,i))            
    end if
End do
!----------------------------------------------------------------------------------------------------------------------------------------


!Andreas, December 2015
!Derivatives of a with respect to xi, xj, and tau
!----------------------------------------------------------------------------------------------------------------------------------------
dnadxi2dtaum = 0.D0
Do i = 1, gl%ncomp - 1
    Do j = 1, gl%ncomp - 1
        !*** d2a_dxidxj
        dnadxi2dtaum(1,i,j) = aij_cubic(j,i) + aij_cubic(i,j) - aij_cubic(j,gl%ncomp) - aij_cubic(gl%ncomp,j) &
                        & + 2.D0 * aij_cubic(gl%ncomp,gl%ncomp) - aij_cubic(gl%ncomp,i) - aij_cubic(i,gl%ncomp)
        !*** d3a_dxidxjdtau
        dnadxi2dtaum(2,i,j) = daij_dtau(j,i) + daij_dtau(i,j) - daij_dtau(j,gl%ncomp) - daij_dtau(gl%ncomp,j) &
                        & + 2.D0 * daij_dtau(gl%ncomp,gl%ncomp) - daij_dtau(gl%ncomp,i) - daij_dtau(i,gl%ncomp)
        !*** d4a_dxidxjdtau2
        if (get_a_deriv > 1) then
            dnadxi2dtaum(3,i,j) = d2aij_dtau2(j,i) + d2aij_dtau2(i,j) - d2aij_dtau2(j,gl%ncomp) - d2aij_dtau2(gl%ncomp,j) &
                & + 2.D0 * d2aij_dtau2(gl%ncomp,gl%ncomp) - d2aij_dtau2(gl%ncomp,i) - d2aij_dtau2(i,gl%ncomp)
        end if
        !*** d5a_dxidxjdtau3
        if (get_a_deriv > 2) then
            dnadxi2dtaum(4,i,j) = d3aij_dtau3(j,i) + d3aij_dtau3(i,j) - d3aij_dtau3(j,gl%ncomp) - d3aij_dtau3(gl%ncomp,j) &
                & + 2.D0 * d3aij_dtau3(gl%ncomp,gl%ncomp) - d3aij_dtau3(gl%ncomp,i) - d3aij_dtau3(i,gl%ncomp)
        end if
        !*** d6a_dxidxjdtau4
        if (get_a_deriv > 3) then
            dnadxi2dtaum(5,i,j) = d4aij_dtau4(j,i) + d4aij_dtau4(i,j) - d4aij_dtau4(j,gl%ncomp) - d4aij_dtau4(gl%ncomp,j) &
                & + 2.D0 * d4aij_dtau4(gl%ncomp,gl%ncomp) - d4aij_dtau4(gl%ncomp,i) - d4aij_dtau4(i,gl%ncomp)
        end if
    End do
End do
!----------------------------------------------------------------------------------------------------------------------------------------

else    !Mix_type = 22, PSRK mixing rules!! These derivatives are done numerically at the moment. Andreas Jäger, January 2017 
        !UPDATE: Numerical erivatives replaced with analytical derivatives, Andreas Jäger, May 2017
    
    dnadtaun = 0.D0
    dnadxidtaum = 0.D0
    dnadxi2dtaum = 0.D0
    

    !!FOR NUMERICAL DERIVATIVES   
    !!-------------------------------------------------------------------------------------------------
    !!Save original values of parameter a and of mole fractions
    !a_PSRK_orig = a_cubic
    !molfractions_orig = molfractions
    !
    !!Specify delta_tau
    !delta_tau = 1.D-5
    !
    !!increase tau and calculate the corresponding temperature (at constant delta and mole fractions)
    !tau_p = tau + delta_tau
    !Temp_p = Tred_cubic / tau_p
    !!decrease tau and calculate the corresponding temperature (at constant delta and mole fractions)
    !tau_m = tau - delta_tau
    !Temp_m = Tred_cubic / tau_m
    !
    !!Calculate parameter a at increased tau
    !call reduced_parameters_calc(gl,Temp_p)
    !a_PSRK_p = a_SRK
    !
    !!Calculate parameter a at decreased tau
    !call reduced_parameters_calc(gl,Temp_m)
    !a_PSRK_m = a_SRK    
    !
    !!Set back to original parameter
    !call reduced_parameters_calc(gl,Temp)  
    !
    !!First derivative of parameter a with respect to tau: da/dtau NUMERICAL
    !dnadtaun(1) = (a_PSRK_p - a_PSRK_m) / (2.D0 * delta_tau)
    !!Second derivative of parameter a with respect to tau: d2a/dtau2 NUMERICAL
    !dnadtaun(2) = (a_PSRK_p - 2.D0 * a_PSRK_orig + a_PSRK_m) / (delta_tau**2)
    !
    !a_SRK = a_PSRK_orig
    !!-------------------------------------------------------------------------------------------------
    
    !Compute the required derivatives of a with respect to tau  
    GETDER = 0
    GETDER(1) = 1   !Needed for derivatives with respect to xa and xb!
    GETDER(4) = 1   !Get derivative of gE with respect to tau
    if (get_a_deriv > 1) then
        GETDER(5) = 1   !Get second derivative of gE with respect to tau    
    end if
    C_or_R = 0      !Both, the combinatorial and the residual part, are needed
    !call gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER, gl%ge%gE_C, gl%ge%gE_R, gl%ge%ln_gamma_C, gl%ge%ln_gamma_R, C_or_R, errval)
    call gE_UNIFAC_MIXDERIVS(gl,Temp, GETDER, C_or_R, errval)

    
    !First derivative of parameter a with respect to tau: da/dtau ANALYTICAL 
    dnadtaun(1) = b_cubic * (gl%ge%gE_C(4) + gl%ge%gE_R(4)) / A1_PSRK
    do i = 1,gl%ncomp
        dnadtaun(1) = dnadtaun(1) + b_cubic * gl%molfractions(i) * dai_dtau(i) / bi_cubic(i) - b_cubic * R_cubic * tred_cubic / tau**2 / A1_PSRK * gl%molfractions(i) * dlog(b_cubic / bi_cubic(i))         
    end do 

    !Second derivative of parameter a with respect to tau: d2a/dtau2 ANALYTICAL
    if (get_a_deriv > 1) then
        dnadtaun(2) = b_cubic * (gl%ge%gE_C(5) + gl%ge%gE_R(5)) / A1_PSRK
        do i = 1,gl%ncomp
            dnadtaun(2) = dnadtaun(2) + b_cubic * gl%molfractions(i) * d2ai_dtau2(i) / bi_cubic(i) + 2.D0 * b_cubic * R_cubic * tred_cubic / tau**3 / A1_PSRK * gl%molfractions(i) * dlog(b_cubic / bi_cubic(i))         
        end do    
    end if

    !Third derivative of parameter a with respect to tau: d3a/dtau3
    dnadtaun(3) = 0.D0  !NOT IMPLEMENTED YET

    !Fourth derivative of parameter a with respect to tau: d4a/dtau4
    dnadtaun(4) = 0.D0  !NOT IMPLEMENTED YET
    
    
    
    if (gl%ncomp > 1) then !Only for mixtures
        
        !!FOR NUMERICAL DERIVATIVES   
        !!----------------------------------------------------------------------------------------------------------------------------
        !!First derivative of parameter a with respect to xi NUMERICAL
        !!---------------------------------------------------------------
        !delta_x = 1.D-7   
        !molfractions_orig = molfractions
        !do i =1, ncomp-1
        !    !Increase xi
        !    molfractions(i) = molfractions_orig(i) + delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !    call reduced_parameters_calc(gl,Temp)
        !    a_PSRK_p = a_SRK
        !    
        !    !Decrease xi
        !    molfractions(i) = molfractions_orig(i) - delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !    call reduced_parameters_calc(gl,Temp)
        !    a_PSRK_m = a_SRK
        !    
        !    !Calculate the numerical derivative
        !    dnadxidtaum(1,i) = (a_PSRK_p - a_PSRK_m) / (2.D0 * delta_x)
        !    molfractions = molfractions_orig
        !end do
        !
        !call reduced_parameters_calc(gl,Temp)
        !a_SRK = a_PSRK_orig
        !!---------------------------------------------------------------
        !
        !
        !
        !!Second derivative of parameter a with respect to tau and xi
        !!---------------------------------------------------------------
        !Do i = 1, ncomp-1
        !    !Increase xi, increase tau
        !    molfractions(i) = molfractions_orig(i) + delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !    tau_p = tau + delta_tau
        !    Temp_pm = tredmix/tau_p
        !    call reduced_parameters_calc(gl,Temp_pm)
        !    a_PSRK_iptaup = a_SRK
        !    
        !    !Increase xi, decrease tau
        !    molfractions(i) = molfractions_orig(i) + delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !    tau_p = tau - delta_tau
        !    Temp_pm = tredmix/tau_p
        !    call reduced_parameters_calc(gl,Temp_pm)
        !    a_PSRK_iptaum = a_SRK
        !
        !    !Decrease xi, increase tau
        !    molfractions(i) = molfractions_orig(i) - delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !    tau_p = tau + delta_tau
        !    Temp_pm = tredmix/tau_p
        !    call reduced_parameters_calc(gl,Temp_pm)
        !    a_PSRK_imtaup = a_SRK
        !
        !     !decrease xi, decrease tau
        !    molfractions(i) = molfractions_orig(i) - delta_x
        !    molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !    tau_p = tau - delta_tau
        !    Temp_pm = tredmix/tau_p
        !    call reduced_parameters_calc(gl,Temp_pm)
        !    a_PSRK_imtaum = a_SRK           
        !    
        !    dnadxidtaum(2,i) = (a_PSRK_iptaup + a_PSRK_imtaum - a_PSRK_iptaum - a_PSRK_imtaup) / 4.D0 / delta_x / delta_tau
        !    
        !    molfractions = molfractions_orig
        !    call reduced_parameters_calc(gl,Temp)
        !    a_SRK = a_PSRK_orig
        !end do    
        !
        !call reduced_parameters_calc(gl,Temp)
        !a_SRK = a_PSRK_orig
        !!---------------------------------------------------------------
        !!----------------------------------------------------------------------------------------------------------------------------
        
        !Compute the required derivatives of a with respect to tau and xa 
        GETDER = 0
        GETDER(1) = 1   !Get derivative of gE with respect to xa
        GETDER(4) = 1   !Get derivative of gE with respect to tau and xa
        C_or_R = 0      !Both, the combinatorial and the residual part, are needed       
        call gE_UNIFAC_MIXDERIVS_dxa(gl,Temp, GETDER, C_or_R, errval)
        
        !Compute required derivatives of b
        call db_SRK_dxi(gl,db_dxi)
        
        !First derivative of parameter a with respect to xi ANALYTICAL
        dnadxidtaum(1,:) = b_cubic * (gl%ge%gE_C_dxa(1,:) + gl%ge%gE_R_dxa(1,:)) / A1_PSRK + db_dxi(:) * (gl%ge%gE_C(1) + gl%ge%gE_R(1)) / A1_PSRK 
        do j = 1, gl%ncomp - 1
            dnadxidtaum(1,j) = dnadxidtaum(1,j) + b_cubic * (ai_cubic(j) / bi_cubic(j) - ai_cubic(gl%ncomp) / bi_cubic(gl%ncomp)) &
                           & + R_cubic / tau / A1_PSRK * tred_cubic * b_cubic * dlog(bi_cubic(gl%ncomp)/bi_cubic(j))
            do i = 1, gl%ncomp
                dnadxidtaum(1,j) = dnadxidtaum(1,j) + db_dxi(j) * gl%molfractions(i) * ai_cubic(i) / bi_cubic(i) &
                               & + R_cubic / tau / A1_PSRK & 
                               & * (tred_cubic * db_dxi(j) * gl%molfractions(i) * dlog(b_cubic/bi_cubic(i)) & 
                               & +  tred_cubic * gl%molfractions(i) * db_dxi(j))
            end do
        end do
        
        !Second derivative of parameter a with respect to tau and xi ANALYTICAL
        dnadxidtaum(2,:) = b_cubic * (gl%ge%gE_C_dxa(4,:) + gl%ge%gE_R_dxa(4,:)) / A1_PSRK + db_dxi(:) * (gl%ge%gE_C(4) + gl%ge%gE_R(4)) / A1_PSRK 
        do j = 1, gl%ncomp - 1
            dnadxidtaum(2,j) = dnadxidtaum(2,j) + b_cubic * (dai_dtau(j) / bi_cubic(j) - dai_dtau(gl%ncomp) / bi_cubic(gl%ncomp)) &
                           & - R_cubic / tau**2 / A1_PSRK * tred_cubic * b_cubic * dlog(bi_cubic(gl%ncomp)/bi_cubic(j))
            do i = 1, gl%ncomp
                dnadxidtaum(2,j) = dnadxidtaum(2,j) + db_dxi(j) * gl%molfractions(i) * dai_dtau(i) / bi_cubic(i) &
                               & - R_cubic / tau**2 / A1_PSRK & 
                               & * (tred_cubic * db_dxi(j) * gl%molfractions(i) * dlog(b_cubic/bi_cubic(i)) & 
                               & +  tred_cubic * gl%molfractions(i) * db_dxi(j))
            end do
        end do
        
    
        !Third derivative of parameter a with respect to tau, tau and xi
        dnadxidtaum(3,:) = 0.D0    !NOT IMPLEMENTED YET
    
        !Fourth derivative of parameter a with respect to tau, tau, tau and xi
        dnadxidtaum(4,:) = 0.D0    !NOT IMPLEMENTED YET
    
        !Fifth derivative of parameter a with respect to tau, tau, tau, tau and xi
        dnadxidtaum(5,:) = 0.D0    !NOT IMPLEMENTED YET
    
    
        !!FOR NUMERICAL DERIVATIVES   
        !!----------------------------------------------------------------------------------------------------------------------------
        !!Second derivative of parameter a with respect to xi, and xj NUMERICAL
        !!---------------------------------------------------------------
        !delta_x = 1.D-4
        !Do i = 1, ncomp-1
        !    do j = i, ncomp-1
        !        !Increase xi, increase xj
        !        molfractions(i) = molfractions_orig(i) + delta_x
        !        molfractions(j) = molfractions_orig(j) + delta_x
        !        if (i==j) then
        !            molfractions(ncomp) = molfractions_orig(ncomp) - delta_x
        !        else
        !            molfractions(ncomp) = molfractions_orig(ncomp) - 2.D0 * delta_x
        !        end if
        !        call reduced_parameters_calc(gl,temp)
        !        a_PSRK_ipjp = a_SRK  
        !    
        !        if (i /= j) then
        !            !Increase xi, decrease xj
        !            molfractions(i) = molfractions_orig(i) + delta_x
        !            molfractions(j) = molfractions_orig(j) - delta_x
        !            molfractions(ncomp) = molfractions_orig(ncomp)
        !            call reduced_parameters_calc(gl,temp)
        !            a_PSRK_ipjm = a_SRK
        !
        !            !Decrease xi, increase xj
        !            molfractions(i) = molfractions_orig(i) - delta_x
        !            molfractions(j) = molfractions_orig(j) + delta_x
        !            molfractions(ncomp) = molfractions_orig(ncomp)
        !            call reduced_parameters_calc(gl,temp)
        !            a_PSRK_imjp = a_SRK
        !        end if
        !        
        !        !decrease xi, decrease xj
        !        molfractions(i) = molfractions_orig(i) - delta_x
        !        molfractions(j) = molfractions_orig(j) - delta_x
        !        if (i == j) then
        !            molfractions(ncomp) = molfractions_orig(ncomp) + delta_x
        !        else
        !            molfractions(ncomp) = molfractions_orig(ncomp) + 2.D0 * delta_x
        !        end if
        !        call reduced_parameters_calc(gl,temp)
        !        a_PSRK_imjm = a_SRK         
        !        
        !        if (i == j) then
        !            dnadxi2dtaum(1,i,j) = (a_PSRK_ipjp - 2.D0 * a_PSRK_orig + a_PSRK_imjm) / (delta_x**2)    
        !        else
        !            dnadxi2dtaum(1,i,j) = (a_PSRK_ipjp + a_PSRK_imjm - a_PSRK_ipjm - a_PSRK_imjp) / 4.D0 / delta_x**2
        !            dnadxi2dtaum(1,j,i) = dnadxi2dtaum(1,i,j)
        !        end if
        !        molfractions = molfractions_orig
        !    end do
        !end do
        !
        !call reduced_parameters_calc(gl,Temp)
        !a_SRK = a_PSRK_orig
        !!---------------------------------------------------------------
        !!----------------------------------------------------------------------------------------------------------------------------
        
        !Compute the required derivatives of a with respect to xa and xb  
        GETDER = 0
        GETDER(1) = 1   !Get derivative of gE with respect to xa and xb
        C_or_R = 0      !Both, the combinatorial and the residual part, are needed       
        call gE_UNIFAC_MIXDERIVS_dxadxb(gl,Temp, GETDER, C_or_R, errval)
        
        
        !Second derivative of parameter a with respect to xi, and xj ANALYTICAL
        do k = 1, gl%ncomp - 1
            do j = 1, gl%ncomp - 1
                dnadxi2dtaum(1,j,k) = b_cubic * (gl%ge%gE_C_dxadxb(1,j,k) + gl%ge%gE_R_dxadxb(1,j,k)) / A1_PSRK &
                                  & + db_dxi(j) * (gl%ge%gE_C_dxa(1,k) + gl%ge%gE_R_dxa(1,k)) / A1_PSRK &
                                  & + db_dxi(k) * (gl%ge%gE_C_dxa(1,j) + gl%ge%gE_R_dxa(1,j)) / A1_PSRK &
                                  & + db_dxi(j) * (ai_cubic(k) / bi_cubic(k) - ai_cubic(gl%ncomp) / bi_cubic(gl%ncomp)) &
                                  & + db_dxi(k) * (ai_cubic(j) / bi_cubic(j) - ai_cubic(gl%ncomp) / bi_cubic(gl%ncomp)) &
                                  & + R_cubic / tau / A1_PSRK &
                                  & * (tred_cubic * db_dxi(j) * dlog(bi_cubic(gl%ncomp)/bi_cubic(k)) &
                                  & + tred_cubic * db_dxi(k) * dlog(bi_cubic(gl%ncomp)/bi_cubic(j)) &
                                  & + tred_cubic / b_cubic * db_dxi(j) * db_dxi(k))
            end do
        end do
        
    
        !Third derivative of parameter a with respect to tau, xi, and xj
        dnadxi2dtaum(2,:,:) = 0.D0      !NOT IMPLEMENTED YET
    
        !Fourth derivative of parameter a with respect to tau, tau and xi, and xj
        dnadxi2dtaum(3,:,:) = 0.D0      !NOT IMPLEMENTED YET
    
        !Fifth derivative of parameter a with respect to tau, tau, tau and xi, and xj
        dnadxi2dtaum(4,:,:) = 0.D0      !NOT IMPLEMENTED YET
    
        !Sixth derivative of parameter a with respect to tau, tau, tau, tau and xi, and xj
        dnadxi2dtaum(5,:,:) = 0.D0      !NOT IMPLEMENTED YET
    end if
    
end if

 deallocate(aij_cubic)
 deallocate(daij_dtau)
 deallocate(d2aij_dtau2)
 deallocate(d3aij_dtau3)
 deallocate(d4aij_dtau4)
 deallocate(kij_cubic)
 deallocate(term_uij_cubic)  
 deallocate(duij_dtau)
 deallocate(d2uij_dtau2)
 deallocate(d3uij_dtau3)
 deallocate(d4uij_dtau4)
 deallocate(term_uij_cubic05)

end subroutine da_dxi_dtau_cubic_all
!**************************************************************************




!**************************************************************************
module subroutine MIXDERIVSFNR_dxi_CUBIC (gl,TEMPERATURE, DENSITY, GETDER, MIXDERFNR_dxi)
!**************************************************************************
! Andreas, December 2015
!--------------------------------------------------------------------------------------------------
! SUBROUTINE FOR THE CALCULATION OF DERIVATIVES OF THE RESIDUAL PART OF THE 
! HELMHOLTZ FREE ENERGY FOR MIXTURES CALCULATED WITH CUBIC EOS (SRK AND PR)
! ALL FIRST COMPOSITION DERIVATIVES ARE CALCULATED IN THIS ROUTINE
!--------------------------------------------------------------------------------------------------
! INPUT PARAMETERS:     
! TEMPERATURE - T   K
! DENSITY     - D   mol/m^3
! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
!               ALL THE LISTED DERIVATIVES ARE DERIVATED WITH RESPECT TO xi
!                1. NORMALIZED RESIDUAL MIXTURE HELMHOLTZ ENERGY F AS A FUNCTION OF DEL AND TAU AND X 
!                2. 1ST DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
!                3. 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
!                4. 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
!                5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
!                6: 1ST MIXED DERIVATIVE OF F WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
!                7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
!                8: 3RD DERIVATIVE OF F WITH RESPECT TO D, MULTIPLIED BY DEL^3
!                9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
!               10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
!               11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^4, MULTIPLIED BY DEL^4
!               12: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^3 and TAU, MULTIPLIED BY DEL^3*TAU
!               13: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^2 and TAU^2, MULTIPLIED BY DEL^2*TAU^2
!               14: 4TH DERIVATIVE OF F WITH RESPECT TO DEL and TAU^3, MULTIPLIED BY DEL*TAU^3
!               15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU^4, MULTIPLIED BY TAU^4
! 
! OUTPUT PARAMETERS: 
! MIXDERFNR_dxi   - A MATRIX WITH 15 ROWS AND 30 COLUMNS WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
!                   AS INDICATED IN "GETDER"
!                   THE COLUMNS ARE FOR THE DIFFERENT COMPOSITIONS
!                   E.g.    MIXDERFNR_dxi(1,1) --> dalphar_dx1
!                           MIXDERFNR_dxi(1,2) --> dalphar_dx2
!                           MIXDERFNR_dxi(2,1) --> d2alphar_ddel_dx1 * del
!                           MIXDERFNR_dxi(2,3) --> d2alphar_ddel_dx3 * del
!                           MIXDERFNR_dxi(10,4) --> d4alphar_ddel2_dtau_dx4 * del^2 * tau
!--------------------------------------------------------------------------------------------------








implicit none

    type(type_gl) :: gl


double precision:: TEMPERATURE, DENSITY  
integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed 
double precision, dimension(15,30)::MIXDERFNR_dxi  ! array with the computed values for the derivatives

double precision:: tau, delta

double precision:: D1_cubic, D2_cubic

!Supporting variables
double precision, dimension(30):: dpsi_m_dxi, d2psi_m_dxiddel, d3psi_m_dxiddel2, d4psi_m_dxiddel3, d5psi_m_dxiddel4
double precision:: psi_p 
double precision:: dpsi_p_ddel, d2psi_p_ddel2, d3psi_p_ddel3, d4psi_p_ddel4
double precision, dimension(30):: dpsi_p_dxi, d2psi_p_dxiddel, d3psi_p_dxiddel2, d4psi_p_dxiddel3, d5psi_p_dxiddel4

Double precision:: Pi_12, dPi12_ddel, d2Pi12_ddel2!, d3Pi12_ddel3
Double precision:: Pi_12_2, Pi_12_3, Pi_12_4

Double precision, dimension(30):: dPI12_dxi, d2PI12_dxiddel, d3PI12_dxiddel2!, d4PI12_dxiddel3, d5PI12_dxiddel4 

Double precision:: help_A
Double precision , dimension(30):: dhelpA_dxi

!Variables for the derivatives of parameter a wrt tau and xi
Double precision:: da_dtau, d2a_dtau2, d3a_dtau3, d4a_dtau4
Double precision, dimension(15):: dnadtaun
double precision, dimension(15,30)::dnadxidtaum       
double precision, dimension(15,30,30)::dnadxi2dtaum
integer :: get_a_deriv

!Variables for the derivatives of the product tau*a wrt tau
Double precision:: dtaua_dtau, d2taua_dtau2, d3taua_dtau3, d4taua_dtau4

Double precision:: a_cubic, b_cubic, R_cubic
Double precision:: rhored_cubic, tred_cubic

Double precision:: tau2, tau3, tau4, del2, del3, del4        !Powers of delta and tau

double precision, dimension(30):: dbdxi_cubic
double precision, dimension(30):: da_dxi, d2a_dxidtau, d3a_dxidtau2, d4a_dxidtau3, d5a_dxidtau4
double precision, dimension(30):: dtaua_dxi, d2taua_dxidtau, d3taua_dxidtau2, d4taua_dxidtau3, d5taua_dxidtau4

!integer:: i,j

!All CUBIC equations can be written in the same form:
!See Bell et al. (2016), "Helmholtz energy translations for common cubic equations of state for use in one-fluid and multi-fluid mixture models", to be submitted to ??
!$$ p = \frac{R*T}{v-b} - \frac{a}{(v+D1*b)*(v+D2*b)} $$
!For the different cubics, it is:
!   vdW: D1 = 0,        D2 = 0
!   SRK: D1 = 1,        D2 = 0
!   PR:  D1 = 1+2**0.5  D2 = 1 - 2**0.5
!The Helmholtz translation of the cubics can then also be written in a universal way:
!$$ \alpha_r = \psi^(-) - \frac{\tau*a}{R*T_r} * \psi^(+) $$
!Where:
!$$ \psi^(-) = -ln(1-b * \delta * \rho_r) $$
!$$ \psi^(+) = \frac{ln(\frac{D1*b*\rho_r*\delta}{D2*b*\rho_r*\delta})}{b*(D1-D2)} 
!The variable \psi^(-) is called psi_m in the code, this variable is only a function of delta
!The variable \psi^(+) is called psi_p in the code, this variable is only a function of delta

If ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then
    !SRK
    D1_cubic = 1.D0
    D2_cubic = 0.D0
    rhored_cubic = gl%rhored_SRK
    tred_cubic = gl%tred_SRK
    a_cubic = gl%a_SRK
    b_cubic = gl%b_SRK
    R_cubic = gl%R_SRK
    call db_SRK_dxi(gl,dbdxi_cubic)
elseif ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then
    !PR
    D1_cubic = consA_PR
    D2_cubic = -consB_PR    !THE MINUS IS VERY IMPORTANT!! NOTE THE DIFFERENT DEFINITIONS OF THE CONSTANTS IN THE PAPER AND IN THE CODE!!!
    rhored_cubic = gl%rhored_PR
    tred_cubic = gl%tred_PR
    a_cubic = gl%a_PR
    b_cubic = gl%b_PR
    R_cubic = gl%R_PR
    call db_PR_dxi(gl,dbdxi_cubic)
end if

delta = Density / rhored_cubic
tau = Tred_cubic / Temperature

!Variables for speeding up computations
tau2 = tau**2
tau3 = tau*tau2
tau4 = tau2*tau2
del2 = delta**2
del3 = del2 * delta
del4 = del2 * del2

!Calculate the derivatives of psi^(-) wrt xi
!psi_m = -dlog(1.D0-b_cubic*delta*rhored_cubic)
dpsi_m_dxi = delta * dbdxi_cubic * rhored_cubic / (1.D0 - b_cubic * delta * rhored_cubic)
d2psi_m_dxiddel = dbdxi_cubic * rhored_cubic / (1.D0 - b_cubic * delta * rhored_cubic)**2
d3psi_m_dxiddel2 = 2.D0 * b_cubic* dbdxi_cubic * rhored_cubic**2 / (1.D0 - b_cubic * delta * rhored_cubic)**3  
d4psi_m_dxiddel3 = 6.D0 * b_cubic**2 * dbdxi_cubic * rhored_cubic**3 / (1.D0 - b_cubic * delta * rhored_cubic)**4
d5psi_m_dxiddel4 = 24.D0 * b_cubic**3 * dbdxi_cubic * rhored_cubic**4 / (1.D0 - b_cubic * delta * rhored_cubic)**5



!Calculate the derivatives of psi^(+) wrt delta and xi
!Calculate supporting variable Pi_12 and its derivatives wrt delta
!---
Pi_12 = (1.D0 + D1_cubic * b_cubic * rhored_cubic * delta) * (1.D0 + D2_cubic * b_cubic * rhored_cubic * delta)
Pi_12_2 = Pi_12 * Pi_12
Pi_12_3 = Pi_12 * Pi_12_2
Pi_12_4 = Pi_12_2 * Pi_12_2
dPi12_ddel = b_cubic * rhored_cubic * (2.D0 * D1_cubic * D2_cubic * b_cubic * delta * rhored_cubic + D1_cubic + D2_cubic)
d2Pi12_ddel2 = 2.D0 * D1_cubic * D2_cubic * b_cubic**2 * rhored_cubic**2
!d3Pi12_ddel3 = 0.D0
!---
!Calculate supporting variable A
!---
help_A = dlog((delta * rhored_cubic * b_cubic * D1_cubic + 1.D0) / (delta * rhored_cubic * b_cubic * D2_cubic + 1.D0))
dhelpA_dxi = delta * rhored_cubic * dbdxi_cubic * (D1_cubic - D2_cubic) / Pi_12
!---
!Calculate derivatives of supporting variable Pi_12 wrt xi and delta
!---
dPI12_dxi = delta * rhored_cubic * dbdxi_cubic * (D1_cubic * (D2_cubic * delta * rhored_cubic * b_cubic + 1.D0) & 
          & + D2_cubic * (D1_cubic * delta * rhored_cubic * b_cubic + 1.D0))
d2PI12_dxiddel = rhored_cubic * dbdxi_cubic * (D1_cubic * (D2_cubic * delta * rhored_cubic * b_cubic + 1.D0) & 
          & + D2_cubic * (D1_cubic * delta * rhored_cubic * b_cubic + 1.D0) & 
          & + 2.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * b_cubic)
d3PI12_dxiddel2 = 4.D0 * D1_cubic * D2_cubic * rhored_cubic**2 * b_cubic * dbdxi_cubic 
!---
!Calculate the delta derivatives
psi_p = dlog((D1_cubic*b_cubic*rhored_cubic*delta+1.D0)/(D2_cubic*b_cubic*rhored_cubic*delta+1.D0)) &
      & /(b_cubic*(D1_cubic - D2_cubic))
dpsi_p_ddel = rhored_cubic / Pi_12
d2psi_p_ddel2 = -rhored_cubic * dPi12_ddel / Pi_12_2
d3psi_p_ddel3 = rhored_cubic * (-Pi_12 * d2Pi12_ddel2 + 2.D0 * dPi12_ddel**2) / Pi_12_3
d4psi_p_ddel4 = rhored_cubic * (6.D0 * Pi_12 * dPi12_ddel * d2Pi12_ddel2 - 6.D0 * dPi12_ddel**3) / Pi_12_4
!Calculate the xi and delta derivatives
dpsi_p_dxi = (b_cubic * dhelpA_dxi - help_A * dbdxi_cubic) / ((D1_cubic - D2_cubic) * b_cubic**2) 
d2psi_p_dxiddel = -rhored_cubic / Pi_12_2 * dPI12_dxi
d3psi_p_dxiddel2 = -rhored_cubic / Pi_12_2 * (d2PI12_dxiddel + 2.D0 * Pi_12 / rhored_cubic * dPi12_ddel * d2psi_p_dxiddel)
d4psi_p_dxiddel3 = -1.D0 / Pi_12_2 * (rhored_cubic * d3PI12_dxiddel2 &
                 & + 2.D0 * d2psi_p_dxiddel * (dPi12_ddel**2 + Pi_12 * d2Pi12_ddel2) + &
                 & 4.D0 * Pi_12 * dPi12_ddel * d3psi_p_dxiddel2)
d5psi_p_dxiddel4 = -1.D0 / Pi_12_2 * (2.D0 * d2psi_p_dxiddel * 3.D0 * dPi12_ddel * d2Pi12_ddel2 + &
                 & 6.D0 * d3psi_p_dxiddel2 * (dPi12_ddel**2 + Pi_12 * d2Pi12_ddel2) + &
                 & 6.D0 * Pi_12 * dPi12_ddel * d4psi_p_dxiddel3)



get_a_deriv = 4 !Get all derivatives of a with respect to tau
call da_dxi_dtau_cubic_all(gl,Temperature, density, get_a_deriv, dnadtaun, dnadxidtaum, dnadxi2dtaum)
!Derivatives of a with respect to tau
da_dtau = dnadtaun(1)
d2a_dtau2 = dnadtaun(2)
d3a_dtau3 = dnadtaun(3)
d4a_dtau4 = dnadtaun(4)
!Derivatives of a with respect to tau and xi
da_dxi = dnadxidtaum(1,:)
d2a_dxidtau = dnadxidtaum(2,:)
d3a_dxidtau2 = dnadxidtaum(3,:)
d4a_dxidtau3 = dnadxidtaum(4,:)
d5a_dxidtau4 = dnadxidtaum(5,:)


!Calculate the derivatives of the product of (a * tau) with respect to tau
dtaua_dtau = tau * da_dtau + a_cubic
d2taua_dtau2 = tau * d2a_dtau2 + 2.D0 * da_dtau
d3taua_dtau3 = tau * d3a_dtau3 + 3.D0 * d2a_dtau2
d4taua_dtau4 = tau * d4a_dtau4 + 4.D0 * d3a_dtau3

!Calculate the derivatives of the product of (a * tau) with respect to tau and xi
dtaua_dxi = tau * da_dxi 
d2taua_dxidtau = tau * d2a_dxidtau + da_dxi
d3taua_dxidtau2 = tau * d3a_dxidtau2 + 2.D0 * d2a_dxidtau
d4taua_dxidtau3 = tau * d4a_dxidtau3 + 3.D0 * d3a_dxidtau2
d5taua_dxidtau4 = tau * d5a_dxidtau4 + 4.D0 * d4a_dxidtau3


!ALL DERIVATIVES ALSO WITH RESPECT TO xi
!alphar
if (GETDER(1) == 1) then
    MIXDERFNR_dxi(1,:) = dpsi_m_dxi - 1.D0 / R_cubic / Tred_cubic * (dtaua_dxi * psi_p + tau * a_cubic * dpsi_p_dxi)
end if

!d_alphar_d_delta * delta
if (GETDER(2) == 1) then
    MIXDERFNR_dxi(2,:) = d2psi_m_dxiddel - &
                       & 1.D0 / R_cubic / Tred_cubic * (dtaua_dxi * dpsi_p_ddel + tau * a_cubic * d2psi_p_dxiddel)
    MIXDERFNR_dxi(2,:) = MIXDERFNR_dxi(2,:) * delta
end if

!d2_alphar_d_delta2 * delta^2
if (GETDER(3) == 1) then
    MIXDERFNR_dxi(3,:) = d3psi_m_dxiddel2 - &
                       & 1.D0 / R_cubic / Tred_cubic * (dtaua_dxi * d2psi_p_ddel2 + tau * a_cubic * d3psi_p_dxiddel2)
    MIXDERFNR_dxi(3,:) = MIXDERFNR_dxi(3,:) * del2
end if

!d_alphar_d_tau * tau
if (GETDER(4) == 1) then
    MIXDERFNR_dxi(4,:) =  -1.D0 / R_cubic / Tred_cubic * (d2taua_dxidtau * psi_p + dtaua_dtau * dpsi_p_dxi)
    MIXDERFNR_dxi(4,:) = MIXDERFNR_dxi(4,:) * tau
end if

!d2_alphar_d_tau2 * tau^2
if (GETDER(5) == 1) then
    MIXDERFNR_dxi(5,:) = -1.D0 / R_cubic / Tred_cubic * (d3taua_dxidtau2 * psi_p + d2taua_dtau2 * dpsi_p_dxi)
    MIXDERFNR_dxi(5,:) = MIXDERFNR_dxi(5,:) * tau2
end if

!d2_alphar_d_delta_d_tau * delta * tau
if (GETDER(6) == 1) then
    MIXDERFNR_dxi(6,:) = -1.D0 / R_cubic / Tred_cubic * (d2taua_dxidtau * dpsi_p_ddel + dtaua_dtau * d2psi_p_dxiddel)
    MIXDERFNR_dxi(6,:) = MIXDERFNR_dxi(6,:) * tau * delta
end if

!d3_alphar_d_delta_d_tau2 * delta * tau^2
if (GETDER(7) == 1) then
    MIXDERFNR_dxi(7,:) = -1.D0 / R_cubic / Tred_cubic * (d3taua_dxidtau2 * dpsi_p_ddel + d2taua_dtau2 * d2psi_p_dxiddel)
    MIXDERFNR_dxi(7,:) = MIXDERFNR_dxi(7,:) * tau2 * delta
end if

!d3_alphar_d_delta3 * delta^3
if (GETDER(8) == 1) then
    MIXDERFNR_dxi(8,:) = d4psi_m_dxiddel3 - &
                       & 1.D0 / R_cubic / Tred_cubic * (dtaua_dxi * d3psi_p_ddel3 + tau * a_cubic * d4psi_p_dxiddel3)
    MIXDERFNR_dxi(8,:) = MIXDERFNR_dxi(8,:) * del3
end if

!d3_alphar_d_tau3 * tau^3
if (GETDER(9) == 1) then
    MIXDERFNR_dxi(9,:) = -1.D0 / R_cubic / Tred_cubic * (d4taua_dxidtau3 * psi_p + d3taua_dtau3 * dpsi_p_dxi)
    MIXDERFNR_dxi(9,:) = MIXDERFNR_dxi(9,:) * tau3
end if

!d3_alphar_d_delta2_dtau * delta^2 * tau
if (GETDER(10) == 1) then
    MIXDERFNR_dxi(10,:) = -1.D0 / R_cubic / Tred_cubic * (d2taua_dxidtau * d2psi_p_ddel2 + dtaua_dtau * d3psi_p_dxiddel2)
    MIXDERFNR_dxi(10,:) = MIXDERFNR_dxi(10,:) * del2 * tau
end if

!d4_alphar_d_delta4 * delta^4
if (GETDER(11) == 1) then
    MIXDERFNR_dxi(11,:) = d5psi_m_dxiddel4 - &
                       & 1.D0 / R_cubic / Tred_cubic * (dtaua_dxi * d4psi_p_ddel4 + tau * a_cubic * d5psi_p_dxiddel4)
    MIXDERFNR_dxi(11,:) = MIXDERFNR_dxi(11,:) * del4
end if

!d4_alphar_d_delta3_dtau * delta^3*tau
if (GETDER(12) == 1) then
    MIXDERFNR_dxi(12,:) = -1.D0 / R_cubic / Tred_cubic * (d2taua_dxidtau * d3psi_p_ddel3 + dtaua_dtau * d4psi_p_dxiddel3)
    MIXDERFNR_dxi(12,:) = MIXDERFNR_dxi(12,:) * del3 * tau    
end if

!d4_alphar_d_delta2_dtau2 * delta^2*tau^2
if (GETDER(13) == 1) then
    MIXDERFNR_dxi(13,:) = -1.D0 / R_cubic / Tred_cubic * (d3taua_dxidtau2 * d2psi_p_ddel2 + d2taua_dtau2 * d3psi_p_dxiddel2)
    MIXDERFNR_dxi(13,:) = MIXDERFNR_dxi(13,:) * del2 * tau2
end if

!d4_alphar_d_delta_dtau3 * delta*tau^3
if (GETDER(14) == 1) then
    MIXDERFNR_dxi(14,:) = -1.D0 / R_cubic / Tred_cubic * (d4taua_dxidtau3 * dpsi_p_ddel + d3taua_dtau3 * d2psi_p_dxiddel)
    MIXDERFNR_dxi(14,:) = MIXDERFNR_dxi(14,:) * delta * tau3
end if

!d4_alphar_d_tau4 * tau^4
if (GETDER(15) == 1) then
    MIXDERFNR_dxi(15,:) = -1.D0 / R_cubic / Tred_cubic * (d5taua_dxidtau4 * psi_p + d4taua_dtau4 * dpsi_p_dxi)
    MIXDERFNR_dxi(15,:) = MIXDERFNR_dxi(15,:) * tau4
end if


End subroutine MIXDERIVSFNR_dxi_CUBIC




!************************************************************************************
module subroutine MIXDERIVSFNR_dxidxj_CUBIC (gl,TEMPERATURE, DENSITY, GETDER)
!************************************************************************************
! Andreas, Januar 2016
!--------------------------------------------------------------------------------------------------
! SUBROUTINE FOR THE CALCULATION OF DERIVATIVES OF THE RESIDUAL PART OF THE 
! HELMHOLTZ FREE ENERGY FOR MIXTURES CALCULATED WITH CUBIC EOS (SRK AND PR)
! ALL SECOND COMPOSITION DERIVATIVES ARE CALCULATED IN THIS ROUTINE
!--------------------------------------------------------------------------------------------------
! INPUT PARAMETERS:     
! TEMPERATURE - T   K
! DENSITY     - D   mol/m^3
! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
!               ALL THE LISTED DERIVATIVES ARE DERIVATED WITH RESPECT TO xi
!                1. NORMALIZED RESIDUAL MIXTURE HELMHOLTZ ENERGY F AS A FUNCTION OF DEL AND TAU AND X 
!                2. 1ST DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
!                3. 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
!                4. 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
!                5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
!                6: 1ST MIXED DERIVATIVE OF F WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
!                7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
!                8: 3RD DERIVATIVE OF F WITH RESPECT TO D, MULTIPLIED BY DEL^3
!                9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
!               10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
!               11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^4, MULTIPLIED BY DEL^4
!               12: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^3 and TAU, MULTIPLIED BY DEL^3*TAU
!               13: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^2 and TAU^2, MULTIPLIED BY DEL^2*TAU^2
!               14: 4TH DERIVATIVE OF F WITH RESPECT TO DEL and TAU^3, MULTIPLIED BY DEL*TAU^3
!               15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU^4, MULTIPLIED BY TAU^4
! 
! OUTPUT PARAMETERS: 
! gl%MIXDERIVFNR_dxidxj   - A FIELD WITH 15 X 30 X 30 VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
!                       AS INDICATED IN "GETDER"
!                       THE ENTRIES IN THE MATRIX ARE BUILT UP AS FOLLOWS:
!                       E.g.    MIXDERFNR_dxi(1,1,1) --> d2alphar_dx1dx1
!                               MIXDERFNR_dxi(1,2,1) --> d2alphar_dx2dx1
!                               MIXDERFNR_dxi(2,1,1) --> d3alphar_ddel_dx1dx1 * del
!                               MIXDERFNR_dxi(2,3,4) --> d3alphar_ddel_dx3dx4 * del
!                               MIXDERFNR_dxi(10,3,5) --> d5alphar_ddel2_dtau_dx3dx5 * del^2 * tau
!--------------------------------------------------------------------------------------------------








implicit none

    type(type_gl) :: gl


double precision:: TEMPERATURE, DENSITY  
integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed 
!SH variable is part of the gl type
!double precision, dimension(15,30,30) ::gl%MIXDERIVFNR_dxidxj  ! array with the computed values for the derivatives

double precision:: tau, delta

double precision:: D1_cubic, D2_cubic

!Supporting variables
double precision, dimension(30):: dpsi_m_dxi, d2psi_m_dxiddel, d3psi_m_dxiddel2, d4psi_m_dxiddel3, d5psi_m_dxiddel4
double precision, dimension(:,:), allocatable:: d2psi_m_dxidxj, d3psi_m_dxidxjddel, d4psi_m_dxidxjddel2, d5psi_m_dxidxjddel3
double precision, dimension(:,:), allocatable:: d6psi_m_dxidxjddel4
double precision:: psi_p 
double precision:: dpsi_p_ddel, d2psi_p_ddel2, d3psi_p_ddel3, d4psi_p_ddel4
double precision, dimension(30):: dpsi_p_dxi, d2psi_p_dxiddel, d3psi_p_dxiddel2, d4psi_p_dxiddel3, d5psi_p_dxiddel4
double precision, dimension(:,:), allocatable:: d2psi_p_dxidxj, d3psi_p_dxidxjddel, d4psi_p_dxidxjddel2, d5psi_p_dxidxjddel3
double precision, dimension(:,:), allocatable:: d6psi_p_dxidxjddel4

Double precision:: Pi_12, dPi12_ddel, d2Pi12_ddel2!, d3Pi12_ddel3
Double precision:: Pi_12_2, Pi_12_3, Pi_12_4

Double precision, dimension(30):: dPI12_dxi, d2PI12_dxiddel, d3PI12_dxiddel2!, d4PI12_dxiddel3, d5PI12_dxiddel4 
Double precision, dimension(:,:), allocatable:: d2PI12_dxidxj, d3PI12_dxidxjddel, d4PI12_dxidxjddel2

Double precision:: help_A
Double precision, dimension(30):: dhelpA_dxi
Double precision, dimension(:,:), allocatable:: d2helpA_dxidxj

!Variables for the derivatives of parameter a wrt tau and xi
Double precision:: da_dtau, d2a_dtau2, d3a_dtau3, d4a_dtau4
Double precision, dimension(nderivs):: dnadtaun
double precision, dimension(:,:), allocatable::dnadxidtaum       
double precision, dimension(15,30,30) ::dnadxi2dtaum
integer :: get_a_deriv

!Variables for the derivatives of the product tau*a wrt tau
Double precision:: dtaua_dtau, d2taua_dtau2, d3taua_dtau3, d4taua_dtau4

Double precision:: a_cubic, b_cubic, R_cubic
Double precision:: rhored_cubic, tred_cubic

Double precision:: tau2, tau3, tau4, del2, del3, del4        !Powers of delta and tau

double precision, dimension(30):: dbdxi_cubic
double precision, dimension(:,:), allocatable:: d2b_dxidxj
double precision, dimension(30):: da_dxi, d2a_dxidtau, d3a_dxidtau2, d4a_dxidtau3, d5a_dxidtau4
double precision, dimension(:,:), allocatable:: d2a_dxidxj, d3a_dxidxjdtau, d4a_dxidxjdtau2, d5a_dxidxjdtau3, d4a_dxidxjdtau4
double precision, dimension(:,:), allocatable:: d6a_dxidxjdtau4
double precision, dimension(30):: dtaua_dxi, d2taua_dxidtau, d3taua_dxidtau2, d4taua_dxidtau3, d5taua_dxidtau4
double precision, dimension(:,:), allocatable:: d2taua_dxidxj, d3taua_dxidxjdtau, d4taua_dxidxjdtau2, d5taua_dxidxjdtau3
double precision, dimension(:,:), allocatable:: d6taua_dxidxjdtau4

integer:: i,j

!Supportive variables to speed up calculations
double precision:: denom_psim, denom_psim_2, denom_psim_3, denom_psim_4, denom_psim_5, denom_psim_6

    if (.not.(allocated(gl%MIXDERIVFNR_dxidxj)))         allocate(gl%MIXDERIVFNR_dxidxj(15,30,30))
    if (.not.(allocated(d2psi_m_dxidxj)))           allocate(d2psi_m_dxidxj(30,30))                          !allocate(d2psi_m_dxidxj(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d3psi_m_dxidxjddel)))       allocate(d3psi_m_dxidxjddel(30,30))                      !allocate(d3psi_m_dxidxjddel(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d4psi_m_dxidxjddel2)))      allocate(d4psi_m_dxidxjddel2(30,30))                     !allocate(d4psi_m_dxidxjddel2(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d5psi_m_dxidxjddel3)))      allocate(d5psi_m_dxidxjddel3(30,30))                     !allocate(d5psi_m_dxidxjddel3(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d6psi_m_dxidxjddel4)))      allocate(d6psi_m_dxidxjddel4(30,30))                     !allocate(d6psi_m_dxidxjddel4(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d2psi_p_dxidxj)))           allocate(d2psi_p_dxidxj(30,30))                          !allocate(d2psi_p_dxidxj(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d3psi_p_dxidxjddel)))       allocate(d3psi_p_dxidxjddel(30,30))                      !allocate(d3psi_p_dxidxjddel(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d4psi_p_dxidxjddel2)))      allocate(d4psi_p_dxidxjddel2(30,30))                     !allocate(d4psi_p_dxidxjddel2(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d5psi_p_dxidxjddel3)))      allocate(d5psi_p_dxidxjddel3(30,30))                     !allocate(d5psi_p_dxidxjddel3(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d6psi_p_dxidxjddel4)))      allocate(d6psi_p_dxidxjddel4(30,30))                     !allocate(d6psi_p_dxidxjddel4(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d2PI12_dxidxj)))            allocate(d2PI12_dxidxj(30,30))                           !allocate(d2PI12_dxidxj(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d3PI12_dxidxjddel)))        allocate(d3PI12_dxidxjddel(30,30))                       !allocate(d3PI12_dxidxjddel(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d4PI12_dxidxjddel2)))       allocate(d4PI12_dxidxjddel2(30,30))                      !allocate(d4PI12_dxidxjddel2(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d2helpA_dxidxj)))           allocate(d2helpA_dxidxj(30,30))                          !allocate(d2helpA_dxidxj(gl%ncomp,gl%ncomp))
    if (.not.(allocated(dnadxidtaum)))              allocate(dnadxidtaum(nderivs,30))                              !allocate(dnadxidtaum(nderivs,gl%ncomp))
    !if (.not.(allocated(dnadxi2dtaum)))             allocate(dnadxi2dtaum(nderivs,gl%ncomp,gl%ncomp))                   !allocate(dnadxi2dtaum(nderivs,gl%ncomp,gl%ncomp))
                                                                                                                         !
    if (.not.(allocated(d2b_dxidxj)))               allocate(d2b_dxidxj(30,30))                              !allocate(d2b_dxidxj(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d2a_dxidxj)))               allocate(d2a_dxidxj(30,30))                              !allocate(d2a_dxidxj(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d3a_dxidxjdtau)))           allocate(d3a_dxidxjdtau(30,30))                          !allocate(d3a_dxidxjdtau(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d4a_dxidxjdtau2)))          allocate(d4a_dxidxjdtau2(30,30))                         !allocate(d4a_dxidxjdtau2(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d5a_dxidxjdtau3)))          allocate(d5a_dxidxjdtau3(30,30))                         !allocate(d5a_dxidxjdtau3(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d4a_dxidxjdtau4)))          allocate(d4a_dxidxjdtau4(30,30))                         !allocate(d4a_dxidxjdtau4(gl%ncomp,gl%ncomp))
    if (.not.(allocated(d6a_dxidxjdtau4)))          allocate(d6a_dxidxjdtau4(30,30))                         !allocate(d6a_dxidxjdtau4(gl%ncomp,gl%ncomp))
                                                                                                                         !
     if (.not.(allocated(d2taua_dxidxj)))           allocate(d2taua_dxidxj(30,30))                           !allocate(d2taua_dxidxj(gl%ncomp,gl%ncomp))
     if (.not.(allocated(d3taua_dxidxjdtau)))       allocate(d3taua_dxidxjdtau(30,30))                       !allocate(d3taua_dxidxjdtau(gl%ncomp,gl%ncomp))
     if (.not.(allocated(d4taua_dxidxjdtau2)))      allocate(d4taua_dxidxjdtau2(30,30))                      !allocate(d4taua_dxidxjdtau2(gl%ncomp,gl%ncomp))
     if (.not.(allocated(d5taua_dxidxjdtau3)))      allocate(d5taua_dxidxjdtau3(30,30))                      !allocate(d5taua_dxidxjdtau3(gl%ncomp,gl%ncomp))
     if (.not.(allocated(d6taua_dxidxjdtau4)))      allocate(d6taua_dxidxjdtau4(30,30))                      !allocate(d6taua_dxidxjdtau4(gl%ncomp,gl%ncomp))

!All CUBIC equations can be written in the same form:
!See Bell et al. (2016), "Helmholtz energy translations for common cubic equations of state for use in one-fluid and multi-fluid mixture models", to be submitted to ??
!$$ p = \frac{R*T}{v-b} - \frac{a}{(v+D1*b)*(v+D2*b)} $$
!For the different cubics, it is:
!   vdW: D1 = 0,        D2 = 0
!   SRK: D1 = 1,        D2 = 0
!   PR:  D1 = 1+2**0.5  D2 = 1 - 2**0.5
!The Helmholtz translation of the cubics can then also be written in a universal way:
!$$ \alpha_r = \psi^(-) - \frac{\tau*a}{R*T_r} * \psi^(+) $$
!Where:
!$$ \psi^(-) = -ln(1-b * \delta * \rho_r) $$
!$$ \psi^(+) = \frac{ln(\frac{D1*b*\rho_r*\delta}{D2*b*\rho_r*\delta})}{b*(D1-D2)} 
!The variable \psi^(-) is called psi_m in the code, this variable is only a function of delta
!The variable \psi^(+) is called psi_p in the code, this variable is only a function of delta

     gl%MIXDERIVFNR_dxidxj = 0.d0
     
If ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then
    !SRK
    D1_cubic = 1.D0
    D2_cubic = 0.D0
    rhored_cubic = gl%rhored_SRK
    tred_cubic = gl%tred_SRK
    a_cubic = gl%a_SRK
    b_cubic = gl%b_SRK
    R_cubic = gl%R_SRK
    call db_SRK_dxi(gl,dbdxi_cubic)
    call d2b_SRK_dxidxj(gl,d2b_dxidxj)
elseif ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then
    !PR
    D1_cubic = consA_PR
    D2_cubic = -consB_PR    !THE MINUS IS VERY IMPORTANT!! NOTE THE DIFFERENT DEFINITIONS OF THE CONSTANTS IN THE PAPER AND IN THE CODE!!!
    rhored_cubic = gl%rhored_PR
    tred_cubic = gl%tred_PR
    a_cubic = gl%a_PR
    b_cubic = gl%b_PR
    R_cubic = gl%R_PR
    call db_PR_dxi(gl,dbdxi_cubic)
    call d2b_PR_dxidxj(gl,d2b_dxidxj)
end if

delta = Density / rhored_cubic
tau = Tred_cubic / Temperature

!Variables for speeding up computations
tau2 = tau**2
tau3 = tau*tau2
tau4 = tau2*tau2
del2 = delta**2
del3 = del2 * delta
del4 = del2 * del2

!Calculate the derivatives of psi^(-) wrt xi
!psi_m = -dlog(1.D0-b_cubic*delta*rhored_cubic)
!Help variables:
denom_psim = (1.D0 - b_cubic * delta * rhored_cubic)
denom_psim_2 = denom_psim**2
denom_psim_3 = denom_psim_2 * denom_psim
denom_psim_4 = denom_psim_2 * denom_psim_2
denom_psim_5 = denom_psim_2 * denom_psim_3
denom_psim_6 = denom_psim_3 * denom_psim_3
dpsi_m_dxi = delta * dbdxi_cubic * rhored_cubic / denom_psim
d2psi_m_dxiddel = dbdxi_cubic * rhored_cubic / denom_psim_2
d3psi_m_dxiddel2 = 2.D0 * b_cubic* dbdxi_cubic * rhored_cubic**2 / denom_psim_3  
d4psi_m_dxiddel3 = 6.D0 * b_cubic**2 * dbdxi_cubic * rhored_cubic**3 / denom_psim_4
d5psi_m_dxiddel4 = 24.D0 * b_cubic**3 * dbdxi_cubic * rhored_cubic**4 / denom_psim_5

!Calculate the derivatives of psi^(-) wrt xi and xj
do i=1, gl%ncomp-1
    do j=1, gl%ncomp-1        
        if (i <= j) then!Speed up calculations, use symmetry
            d2psi_m_dxidxj(i,j) = delta * rhored_cubic * d2b_dxidxj(i,j) / denom_psim + &
                            & delta**2 * rhored_cubic**2 * dbdxi_cubic(i) * dbdxi_cubic(j) / denom_psim_2 
            d3psi_m_dxidxjddel(i,j) = rhored_cubic * d2b_dxidxj(i,j) / denom_psim_2 + &
                            & 2.D0 * delta * rhored_cubic**2 * dbdxi_cubic(i) * dbdxi_cubic(j) / denom_psim_3       
            d4psi_m_dxidxjddel2(i,j) = 2.D0 * rhored_cubic**2 * d2b_dxidxj(i,j) / denom_psim_3 + &
                            & 2.D0 * rhored_cubic**2 * dbdxi_cubic(i) * dbdxi_cubic(j) / denom_psim_4 * &
                            & (2.D0 * delta * rhored_cubic * b_cubic + 1.D0)
            d5psi_m_dxidxjddel3(i,j) = 6.D0 * rhored_cubic**3 * d2b_dxidxj(i,j) / denom_psim_4 + &
                            & 6.D0 * rhored_cubic**3 * b_cubic * dbdxi_cubic(i) * dbdxi_cubic(j) / denom_psim_5 * &
                            & (2.D0 * delta * rhored_cubic * b_cubic + 2.D0)
            d6psi_m_dxidxjddel4(i,j) = 24.D0 * rhored_cubic**4 * d2b_dxidxj(i,j) / denom_psim_5 + &
                            & 24.D0 * rhored_cubic**4 * b_cubic**2 * dbdxi_cubic(i) * dbdxi_cubic(j) / denom_psim_6 * &
                            & (2.D0 * delta * rhored_cubic * b_cubic + 3.D0)
        else 
            d2psi_m_dxidxj(i,j) = d2psi_m_dxidxj(j,i) 
            d3psi_m_dxidxjddel(i,j) = d3psi_m_dxidxjddel(j,i)
            d4psi_m_dxidxjddel2(i,j) = d4psi_m_dxidxjddel2(j,i) 
            d5psi_m_dxidxjddel3(i,j) = d5psi_m_dxidxjddel3(j,i)
            d6psi_m_dxidxjddel4(i,j) = d6psi_m_dxidxjddel4(j,i)
        end if
       end do
end do

!Calculate the derivatives of psi^(+) wrt delta, xi, and xj
!Calculate supporting variable Pi_12 and its derivatives wrt delta
!---
Pi_12 = (1.D0 + D1_cubic * b_cubic * rhored_cubic * delta) * (1.D0 + D2_cubic * b_cubic * rhored_cubic * delta)
Pi_12_2 = Pi_12 * Pi_12
Pi_12_3 = Pi_12 * Pi_12_2
Pi_12_4 = Pi_12_2 * Pi_12_2
dPi12_ddel = b_cubic * rhored_cubic * (2.D0 * D1_cubic * D2_cubic * b_cubic * delta * rhored_cubic + D1_cubic + D2_cubic)
d2Pi12_ddel2 = 2.D0 * D1_cubic * D2_cubic * b_cubic**2 * rhored_cubic**2
!d3Pi12_ddel3 = 0.D0
!---
!Calculate derivatives of supporting variable Pi_12 wrt xi and delta
!---
dPI12_dxi = delta * rhored_cubic * dbdxi_cubic * (D1_cubic * (D2_cubic * delta * rhored_cubic * b_cubic + 1.D0) & 
          & + D2_cubic * (D1_cubic * delta * rhored_cubic * b_cubic + 1.D0))
d2PI12_dxiddel = rhored_cubic * dbdxi_cubic * (D1_cubic * (D2_cubic * delta * rhored_cubic * b_cubic + 1.D0) & 
          & + D2_cubic * (D1_cubic * delta * rhored_cubic * b_cubic + 1.D0) & 
          & + 2.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * b_cubic)
d3PI12_dxiddel2 = 4.D0 * D1_cubic * D2_cubic * rhored_cubic**2 * b_cubic * dbdxi_cubic 
!---
!Calculate derivatives of supporting variable Pi_12 wrt xi, xj, and delta
!---
do i=1,gl%ncomp-1
    do j=1,gl%ncomp-1
        if (i <= j) then!Speed up calculations, use symmetry
            d2PI12_dxidxj(i,j) = delta * rhored_cubic * & 
                      & (2.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * dbdxi_cubic(i) * dbdxi_cubic(j) + &
                      & (2.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * b_cubic + D1_cubic + D2_cubic) * d2b_dxidxj(i,j))
            d3PI12_dxidxjddel(i,j) = rhored_cubic * & 
                      & (4.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * dbdxi_cubic(i) * dbdxi_cubic(j) + &
                      & (4.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * b_cubic + D1_cubic + D2_cubic) * d2b_dxidxj(i,j))
            d4PI12_dxidxjddel2(i,j) = 4.D0 * D1_cubic * D2_cubic * rhored_cubic**2 * &
                                & (b_cubic * d2b_dxidxj(i,j) + dbdxi_cubic(i) * dbdxi_cubic(j)) 
        else
            d2PI12_dxidxj(i,j) = d2PI12_dxidxj(j,i)
            d3PI12_dxidxjddel(i,j) = d3PI12_dxidxjddel(j,i)
            d4PI12_dxidxjddel2(i,j) = d4PI12_dxidxjddel2(j,i)
        end if
    end do
end do
!---
!Calculate supporting variable A
!---
help_A = dlog((delta * rhored_cubic * b_cubic * D1_cubic + 1.D0) / (delta * rhored_cubic * b_cubic * D2_cubic + 1.D0))
dhelpA_dxi = delta * rhored_cubic * dbdxi_cubic * (D1_cubic - D2_cubic) / Pi_12
do i = 1, gl%ncomp-1
    do j=1,gl%ncomp-1
        if (i <= j) then!Speed up calculations, use symmetry
            d2helpA_dxidxj(i,j) = delta * rhored_cubic * (D1_cubic - D2_cubic) / Pi_12_2 * &
                            & (Pi_12 * d2b_dxidxj(i,j) - dPI12_dxi(j) * dbdxi_cubic(i))
        else
            d2helpA_dxidxj(i,j) = d2helpA_dxidxj(j,i)
        end if
    end do
end do
!---
!Calculate the delta derivatives
psi_p = dlog((D1_cubic*b_cubic*rhored_cubic*delta+1.D0)/(D2_cubic*b_cubic*rhored_cubic*delta+1.D0)) &
      & /(b_cubic*(D1_cubic - D2_cubic))
dpsi_p_ddel = rhored_cubic / Pi_12
d2psi_p_ddel2 = -rhored_cubic * dPi12_ddel / Pi_12_2
d3psi_p_ddel3 = rhored_cubic * (-Pi_12 * d2Pi12_ddel2 + 2.D0 * dPi12_ddel**2) / Pi_12_3
d4psi_p_ddel4 = rhored_cubic * (6.D0 * Pi_12 * dPi12_ddel * d2Pi12_ddel2 - 6.D0 * dPi12_ddel**3) / Pi_12_4
!Calculate the xi and delta derivatives
dpsi_p_dxi = (b_cubic * dhelpA_dxi - help_A * dbdxi_cubic) / ((D1_cubic - D2_cubic) * b_cubic**2) 
d2psi_p_dxiddel = -rhored_cubic / Pi_12_2 * dPI12_dxi
d3psi_p_dxiddel2 = -rhored_cubic / Pi_12_2 * (d2PI12_dxiddel + 2.D0 * Pi_12 / rhored_cubic * dPi12_ddel * d2psi_p_dxiddel)
d4psi_p_dxiddel3 = -1.D0 / Pi_12_2 * (rhored_cubic * d3PI12_dxiddel2 &
                 & + 2.D0 * d2psi_p_dxiddel * (dPi12_ddel**2 + Pi_12 * d2Pi12_ddel2) + &
                 & 4.D0 * Pi_12 * dPi12_ddel * d3psi_p_dxiddel2)
d5psi_p_dxiddel4 = -1.D0 / Pi_12_2 * (2.D0 * d2psi_p_dxiddel * 3.D0 * dPi12_ddel * d2Pi12_ddel2 + &
                 & 6.D0 * d3psi_p_dxiddel2 * (dPi12_ddel**2 + Pi_12 * d2Pi12_ddel2) + &
                 & 6.D0 * Pi_12 * dPi12_ddel * d4psi_p_dxiddel3)
!Calculate the xi, xj, and delta derivatives
do i = 1, gl%ncomp-1
    do j=1,gl%ncomp-1
        if (i <= j) then!Speed up calculations, use symmetry
            d2psi_p_dxidxj(i,j) = (-b_cubic * (help_A * d2b_dxidxj(i,j) + dhelpA_dxi(i) * dbdxi_cubic(j) + & 
                            & dhelpA_dxi(j) * dbdxi_cubic(i)) + 2.D0 * help_A * dbdxi_cubic(i) * dbdxi_cubic(j) + &
                            & b_cubic**2 * d2helpA_dxidxj(i,j)) / ((D1_cubic - D2_cubic) * b_cubic**3)
            d3psi_p_dxidxjddel(i,j) = -1.D0 / Pi_12_2 * &
                                & (rhored_cubic * d2PI12_dxidxj(i,j) + 2.D0 * Pi_12 * dPI12_dxi(j) * d2psi_p_dxiddel(i))
            d4psi_p_dxidxjddel2(i,j) = -1.D0 / Pi_12_2 * (rhored_cubic * d3PI12_dxidxjddel(i,j) + &
                                 & 2.D0 * d2psi_p_dxiddel(i) * &
                                 & (Pi_12 * d2PI12_dxiddel(j) + dPi12_ddel * dPI12_dxi(j)) + 2.D0 * Pi_12 * dPi12_ddel * &
                                 & d3psi_p_dxidxjddel(i,j) + 2.D0 * Pi_12 * dPI12_dxi(j) * d3psi_p_dxiddel2(i))
            d5psi_p_dxidxjddel3(i,j) = -1.D0 / Pi_12_2 * (rhored_cubic * d4PI12_dxidxjddel2(i,j) + &
                                 & 2.d0 * d3psi_p_dxidxjddel(i,j) * (Pi_12 * d2Pi12_ddel2 + dPi12_ddel**2) + &
                                 & 4.D0 * d3psi_p_dxiddel2(i) * (Pi_12 * d2PI12_dxiddel(j) + dPI12_ddel * dPI12_dxi(j)) + &
                                 & 2.D0 * d2psi_p_dxiddel(i) * (Pi_12 * d3PI12_dxiddel2(j) + &
                                 & 2.D0 * dPi12_ddel * d2PI12_dxiddel(j) + dPI12_dxi(j) * d2Pi12_ddel2) + &
                                 & 4.D0 * Pi_12 * dPi12_ddel * d4psi_p_dxidxjddel2(i,j) + &
                                 & 2.D0 * Pi_12 * dPI12_dxi(j) * d4psi_p_dxiddel3(i))
            d6psi_p_dxidxjddel4(i,j) = -1.D0 / Pi_12_2 * &
                                 & (6.D0 * d4psi_p_dxidxjddel2(i,j) * (Pi_12 * d2Pi12_ddel2 + dPi12_ddel**2) + & 
                                 & 2.D0 * d3psi_p_dxidxjddel(i,j) * 3.D0 * dPi12_ddel * d2Pi12_ddel2 + &
                                 & 6.D0 * d4psi_p_dxiddel3(i) * (Pi_12 * d2PI12_dxiddel(j) + dPi12_ddel * dPi12_dxi(j)) +&
                                 & 6.D0 * d3psi_p_dxiddel2(i) * (Pi_12 * d3PI12_dxiddel2(j) + & 
                                 & 2.D0 * dPi12_ddel * d2PI12_dxiddel(j) + dPi12_dxi(j) * d2Pi12_ddel2) + &
                                 & 2.D0 * d2psi_p_dxiddel(i) * (3.D0 * dPi12_ddel * d3PI12_dxiddel2(j) + & 
                                 & 3.D0 * d2Pi12_ddel2 * d2PI12_dxiddel(j)) + &
                                 & 6.D0 * Pi_12 * dPi12_ddel * d5psi_p_dxidxjddel3(i,j) + &
                                 & 2.D0 * Pi_12 * dPi12_dxi(j) * d5psi_p_dxiddel4(i)) 
        else
            d2psi_p_dxidxj(i,j) = d2psi_p_dxidxj(j,i)
            d3psi_p_dxidxjddel(i,j) = d3psi_p_dxidxjddel(j,i) 
            d4psi_p_dxidxjddel2(i,j) = d4psi_p_dxidxjddel2(j,i)
            d5psi_p_dxidxjddel3(i,j) = d5psi_p_dxidxjddel3(j,i)
            d6psi_p_dxidxjddel4(i,j) = d6psi_p_dxidxjddel4(j,i)
        end if
    end do
end do


get_a_deriv = 4 !Get all derivatives of a with respect to tau
call da_dxi_dtau_cubic_all(gl,Temperature, density, get_a_deriv, dnadtaun, dnadxidtaum, dnadxi2dtaum)
!Derivatives of a with respect to tau
da_dtau = dnadtaun(1)
d2a_dtau2 = dnadtaun(2)
d3a_dtau3 = dnadtaun(3)
d4a_dtau4 = dnadtaun(4)
!Derivatives of a with respect to tau and xi
da_dxi = dnadxidtaum(1,:)
d2a_dxidtau = dnadxidtaum(2,:)
d3a_dxidtau2 = dnadxidtaum(3,:)
d4a_dxidtau3 = dnadxidtaum(4,:)
d5a_dxidtau4 = dnadxidtaum(5,:)
!Derivatives of a with respect to tau, xi, and xj
do i = 1, gl%ncomp-1
    do j = 1, gl%ncomp-1
        d2a_dxidxj(i,j) = dnadxi2dtaum(1,i,j)
        d3a_dxidxjdtau(i,j) = dnadxi2dtaum(2,i,j)
        d4a_dxidxjdtau2(i,j) = dnadxi2dtaum(3,i,j)
        d5a_dxidxjdtau3(i,j) = dnadxi2dtaum(4,i,j)
        d6a_dxidxjdtau4(i,j) = dnadxi2dtaum(5,i,j)
    end do
end do


!Calculate the derivatives of the product of (a * tau) with respect to tau
dtaua_dtau = tau * da_dtau + a_cubic
d2taua_dtau2 = tau * d2a_dtau2 + 2.D0 * da_dtau
d3taua_dtau3 = tau * d3a_dtau3 + 3.D0 * d2a_dtau2
d4taua_dtau4 = tau * d4a_dtau4 + 4.D0 * d3a_dtau3

!Calculate the derivatives of the product of (a * tau) with respect to tau and xi
dtaua_dxi = tau * da_dxi 
d2taua_dxidtau = tau * d2a_dxidtau + da_dxi
d3taua_dxidtau2 = tau * d3a_dxidtau2 + 2.D0 * d2a_dxidtau
d4taua_dxidtau3 = tau * d4a_dxidtau3 + 3.D0 * d3a_dxidtau2
d5taua_dxidtau4 = tau * d5a_dxidtau4 + 4.D0 * d4a_dxidtau3

!Calculate the derivatives of the product of (a * tau) with respect to tau, xi, and xj
do i = 1, gl%ncomp-1
    do j = 1, gl%ncomp-1
        d2taua_dxidxj(i,j) = tau * d2a_dxidxj(i,j) 
        d3taua_dxidxjdtau(i,j) = tau * d3a_dxidxjdtau(i,j) + d2a_dxidxj(i,j)
        d4taua_dxidxjdtau2(i,j) = tau * d4a_dxidxjdtau2(i,j) + 2.D0 * d3a_dxidxjdtau(i,j)
        d5taua_dxidxjdtau3(i,j) = tau * d5a_dxidxjdtau3(i,j) + 3.D0 * d4a_dxidxjdtau2(i,j)
        d6taua_dxidxjdtau4(i,j) = tau * d6a_dxidxjdtau4(i,j) + 4.D0 * d5a_dxidxjdtau3(i,j)
    end do
end do

!ALL DERIVATIVES ALSO WITH RESPECT TO xi and xj
do i =1, gl%ncomp - 1
    do j =1, gl%ncomp -1
        !alphar
        if (GETDER(1) == 1) then
            gl%MIXDERIVFNR_dxidxj(1,i,j) = d2psi_m_dxidxj(i,j) - 1.D0 / R_cubic / Tred_cubic * &
                                    & (d2taua_dxidxj(i,j) * psi_p + dtaua_dxi(i) * dpsi_p_dxi(j) + &
                                    & dtaua_dxi(j) * dpsi_p_dxi(i) + tau * a_cubic * d2psi_p_dxidxj(i,j)) 
        end if

        !d_alphar_d_delta * delta
        if (GETDER(2) == 1) then
            gl%MIXDERIVFNR_dxidxj(2,i,j) = d3psi_m_dxidxjddel(i,j) - 1.D0 / R_cubic / Tred_cubic * &
                                    & (d2taua_dxidxj(i,j) * dpsi_p_ddel + dtaua_dxi(i) * d2psi_p_dxiddel(j) + &
                                    & dtaua_dxi(j) * d2psi_p_dxiddel(i) + tau * a_cubic * d3psi_p_dxidxjddel(i,j))
            gl%MIXDERIVFNR_dxidxj(2,i,j) = gl%MIXDERIVFNR_dxidxj(2,i,j) * delta
        end if

        !d2_alphar_d_delta2 * delta^2
        if (GETDER(3) == 1) then
            gl%MIXDERIVFNR_dxidxj(3,i,j) = d4psi_m_dxidxjddel2(i,j) - 1.D0 / R_cubic / Tred_cubic * &
                                    & (d2taua_dxidxj(i,j) * d2psi_p_ddel2 + dtaua_dxi(i) * d3psi_p_dxiddel2(j) + &
                                    & dtaua_dxi(j) * d3psi_p_dxiddel2(i) + tau * a_cubic * d4psi_p_dxidxjddel2(i,j)) 
            gl%MIXDERIVFNR_dxidxj(3,i,j) = gl%MIXDERIVFNR_dxidxj(3,i,j) * del2
        end if

        !d_alphar_d_tau * tau
        if (GETDER(4) == 1) then
            gl%MIXDERIVFNR_dxidxj(4,i,j) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d3taua_dxidxjdtau(i,j) * psi_p + d2taua_dxidtau(i) * dpsi_p_dxi(j) + &
                                    & d2taua_dxidtau(j) * dpsi_p_dxi(i) + dtaua_dtau * d2psi_p_dxidxj(i,j))
            gl%MIXDERIVFNR_dxidxj(4,i,j) = gl%MIXDERIVFNR_dxidxj(4,i,j) * tau
        end if

        !d2_alphar_d_tau2 * tau^2
        if (GETDER(5) == 1) then
            gl%MIXDERIVFNR_dxidxj(5,i,j) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d4taua_dxidxjdtau2(i,j) * psi_p + d3taua_dxidtau2(i) * dpsi_p_dxi(j) + &
                                    & d3taua_dxidtau2(j) * dpsi_p_dxi(i) + d2taua_dtau2 * d2psi_p_dxidxj(i,j))
            gl%MIXDERIVFNR_dxidxj(5,i,j) = gl%MIXDERIVFNR_dxidxj(5,i,j) * tau2
        end if

        !d2_alphar_d_delta_d_tau * delta * tau
        if (GETDER(6) == 1) then
            gl%MIXDERIVFNR_dxidxj(6,i,j) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d3taua_dxidxjdtau(i,j) * dpsi_p_ddel + d2taua_dxidtau(i) * d2psi_p_dxiddel(j) + &
                                    & d2taua_dxidtau(j) * d2psi_p_dxiddel(i) + dtaua_dtau * d3psi_p_dxidxjddel(i,j))
            gl%MIXDERIVFNR_dxidxj(6,i,j) = gl%MIXDERIVFNR_dxidxj(6,i,j) * tau * delta
        end if

        !d3_alphar_d_delta_d_tau2 * delta * tau^2
        if (GETDER(7) == 1) then
            gl%MIXDERIVFNR_dxidxj(7,i,j) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d4taua_dxidxjdtau2(i,j) * dpsi_p_ddel + d3taua_dxidtau2(i) * d2psi_p_dxiddel(j) + &
                                    & d3taua_dxidtau2(j) * d2psi_p_dxiddel(i) + d2taua_dtau2 * d3psi_p_dxidxjddel(i,j))
            gl%MIXDERIVFNR_dxidxj(7,i,j) = gl%MIXDERIVFNR_dxidxj(7,i,j) * tau2 * delta
        end if

        !d3_alphar_d_delta3 * delta^3
        if (GETDER(8) == 1) then
            gl%MIXDERIVFNR_dxidxj(8,i,j) = d5psi_m_dxidxjddel3(i,j) - 1.D0 / R_cubic / Tred_cubic * &
                                    & (d2taua_dxidxj(i,j) * d3psi_p_ddel3 + dtaua_dxi(i) * d4psi_p_dxiddel3(j) + &
                                    & dtaua_dxi(j) * d4psi_p_dxiddel3(i) + tau * a_cubic * d5psi_p_dxidxjddel3(i,j))
            gl%MIXDERIVFNR_dxidxj(8,i,j) = gl%MIXDERIVFNR_dxidxj(8,i,j) * del3
        end if

        !d3_alphar_d_tau3 * tau^3
        if (GETDER(9) == 1) then
            gl%MIXDERIVFNR_dxidxj(9,i,j) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d5taua_dxidxjdtau3(i,j) * psi_p + d4taua_dxidtau3(i) * dpsi_p_dxi(j) + &
                                    & d4taua_dxidtau3(j) * dpsi_p_dxi(i) + d3taua_dtau3 * d2psi_p_dxidxj(i,j))
            gl%MIXDERIVFNR_dxidxj(9,i,j) = gl%MIXDERIVFNR_dxidxj(9,i,j) * tau3
        end if

        !d3_alphar_d_delta2_dtau * delta^2 * tau
        if (GETDER(10) == 1) then
            gl%MIXDERIVFNR_dxidxj(10,i,j) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d3taua_dxidxjdtau(i,j) * d2psi_p_ddel2 + d2taua_dxidtau(i) * d3psi_p_dxiddel2(j) + &
                                    & d2taua_dxidtau(j) * d3psi_p_dxiddel2(i) + dtaua_dtau * d4psi_p_dxidxjddel2(i,j))
            gl%MIXDERIVFNR_dxidxj(10,i,j) = gl%MIXDERIVFNR_dxidxj(10,i,j) * del2 * tau
        end if

        !d4_alphar_d_delta4 * delta^4
        if (GETDER(11) == 1) then
            gl%MIXDERIVFNR_dxidxj(11,i,j) = d6psi_m_dxidxjddel4(i,j) - 1.D0 / R_cubic / Tred_cubic * &
                                    & (d2taua_dxidxj(i,j) * d4psi_p_ddel4 + dtaua_dxi(i) * d5psi_p_dxiddel4(j) + &
                                    & dtaua_dxi(j) * d5psi_p_dxiddel4(i) + tau * a_cubic * d6psi_p_dxidxjddel4(i,j))
            gl%MIXDERIVFNR_dxidxj(11,i,j) = gl%MIXDERIVFNR_dxidxj(11,i,j) * del4
        end if

        !d4_alphar_d_delta3_dtau * delta^3*tau
        if (GETDER(12) == 1) then
            gl%MIXDERIVFNR_dxidxj(12,i,j) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d3taua_dxidxjdtau(i,j) * d3psi_p_ddel3 + d2taua_dxidtau(i) * d4psi_p_dxiddel3(j) + &
                                    & d2taua_dxidtau(j) * d4psi_p_dxiddel3(i) + dtaua_dtau * d5psi_p_dxidxjddel3(i,j))
            gl%MIXDERIVFNR_dxidxj(12,i,j) = gl%MIXDERIVFNR_dxidxj(12,i,j) * del3 * tau    
        end if

        !d4_alphar_d_delta2_dtau2 * delta^2*tau^2
        if (GETDER(13) == 1) then
            gl%MIXDERIVFNR_dxidxj(13,i,j) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d4taua_dxidxjdtau2(i,j) * d2psi_p_ddel2 + d3taua_dxidtau2(i) * d3psi_p_dxiddel2(j) + &
                                    & d3taua_dxidtau2(j) * d3psi_p_dxiddel2(i) + d2taua_dtau2 * d4psi_p_dxidxjddel2(i,j))
            gl%MIXDERIVFNR_dxidxj(13,i,j) = gl%MIXDERIVFNR_dxidxj(13,i,j) * del2 * tau2
        end if

        !d4_alphar_d_delta_dtau3 * delta*tau^3
        if (GETDER(14) == 1) then
            gl%MIXDERIVFNR_dxidxj(14,i,j) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d5taua_dxidxjdtau3(i,j) * dpsi_p_ddel + d4taua_dxidtau3(i) * d2psi_p_dxiddel(j) + &
                                    & d4taua_dxidtau3(j) * d2psi_p_dxiddel(i) + d3taua_dtau3 * d3psi_p_dxidxjddel(i,j))
            gl%MIXDERIVFNR_dxidxj(14,i,j) = gl%MIXDERIVFNR_dxidxj(14,i,j) * delta * tau3
        end if

        !d4_alphar_d_tau4 * tau^4
        if (GETDER(15) == 1) then
            gl%MIXDERIVFNR_dxidxj(15,i,j) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d6taua_dxidxjdtau4(i,j) * psi_p + d5taua_dxidtau4(i) * dpsi_p_dxi(j) + &
                                    & d5taua_dxidtau4(j) * dpsi_p_dxi(i) + d4taua_dtau4 * d2psi_p_dxidxj(i,j))
            gl%MIXDERIVFNR_dxidxj(15,i,j) = gl%MIXDERIVFNR_dxidxj(15,i,j) * tau4
        end if
    end do
end do

End subroutine MIXDERIVSFNR_dxidxj_CUBIC



!******************************************************************************************
module subroutine MIXDERIVSFNR_dxidxjdxk_CUBIC (gl,TEMPERATURE, DENSITY, GETDER)
!******************************************************************************************
! Andreas, February 2016
!--------------------------------------------------------------------------------------------------
! SUBROUTINE FOR THE CALCULATION OF DERIVATIVES OF THE RESIDUAL PART OF THE 
! HELMHOLTZ FREE ENERGY FOR MIXTURES CALCULATED WITH CUBIC EOS (SRK AND PR)
! ALL THIRD COMPOSITION DERIVATIVES ARE CALCULATED IN THIS ROUTINE
!--------------------------------------------------------------------------------------------------
! INPUT PARAMETERS:     
! TEMPERATURE - T   K
! DENSITY     - D   mol/m^3
! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
!               ALL THE LISTED DERIVATIVES ARE DERIVATED WITH RESPECT TO xi
!                1. NORMALIZED RESIDUAL MIXTURE HELMHOLTZ ENERGY F AS A FUNCTION OF DEL AND TAU AND X 
!                2. 1ST DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
!                X 3. 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
!                4. 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
!                X 5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
!                X 6: 1ST MIXED DERIVATIVE OF F WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
!                X 7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO D, T, AND T, MULTIPLIED BY TAU*TAU*DEL
!                X 8: 3RD DERIVATIVE OF F WITH RESPECT TO D, MULTIPLIED BY DEL^3
!                X 9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
!                X 10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
!                X 11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^4, MULTIPLIED BY DEL^4
!                X 12: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^3 and TAU, MULTIPLIED BY DEL^3*TAU
!                X 13: 4TH DERIVATIVE OF F WITH RESPECT TO DEL^2 and TAU^2, MULTIPLIED BY DEL^2*TAU^2
!                X 14: 4TH DERIVATIVE OF F WITH RESPECT TO DEL and TAU^3, MULTIPLIED BY DEL*TAU^3
!                X 15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU^4, MULTIPLIED BY TAU^4
! 
! OUTPUT PARAMETERS: 
! gl%MIXDERIVFNR_dxidxjdxk   - A FIELD WITH 15 X 30 X 30 X 30 VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
!                           AS INDICATED IN "GETDER"
!                           THE ENTRIES IN THE MATRIX ARE BUILT UP AS FOLLOWS:
!                           E.g.    MIXDERFNR_dxi(1,1,1,1) --> d3alphar_dx1dx1dx1
!                                   MIXDERFNR_dxi(1,2,1,3) --> d3alphar_dx2dx1dx3
!                                   MIXDERFNR_dxi(2,1,1,2) --> d4alphar_ddel_dx1dx1dx2 * del
!                                   MIXDERFNR_dxi(2,3,4,5) --> d4alphar_ddel_dx3dx4dx5 * del
!                                   MIXDERFNR_dxi(4,3,5,7) --> d6alphar_dtau_dx3dx5dx7 * tau
!--------------------------------------------------------------------------------------------------








implicit none

    type(type_gl) :: gl


double precision:: TEMPERATURE, DENSITY  
integer, dimension(15):: GETDER            ! array specifier to indicate, which derivative is needed 
!SH Variable is part of gl
!double precision, dimension(:,:,:,:),allocatable::gl%MIXDERIVFNR_dxidxjdxk  ! array with the computed values for the derivatives

double precision:: tau, delta

double precision:: D1_cubic, D2_cubic

!Supporting variables
double precision, dimension(30):: dpsi_m_dxi, d2psi_m_dxiddel, d3psi_m_dxiddel2, d4psi_m_dxiddel3, d5psi_m_dxiddel4
double precision, dimension(:,:),allocatable :: d2psi_m_dxidxj, d3psi_m_dxidxjddel, d4psi_m_dxidxjddel2, d5psi_m_dxidxjddel3
double precision, dimension(:,:),allocatable :: d6psi_m_dxidxjddel4
double precision, dimension(:,:,:),allocatable :: d3psi_m_dxidxjdxk, d4psi_m_dxidxjdxkddel
double precision:: psi_p 
double precision:: dpsi_p_ddel, d2psi_p_ddel2, d3psi_p_ddel3, d4psi_p_ddel4
double precision, dimension(30):: dpsi_p_dxi, d2psi_p_dxiddel, d3psi_p_dxiddel2, d4psi_p_dxiddel3, d5psi_p_dxiddel4
double precision, dimension(:,:),allocatable :: d2psi_p_dxidxj, d3psi_p_dxidxjddel, d4psi_p_dxidxjddel2, d5psi_p_dxidxjddel3
double precision, dimension(:,:),allocatable :: d6psi_p_dxidxjddel4
double precision, dimension(:,:,:),allocatable :: d3psi_p_dxidxjdxk, d4psi_p_dxidxjdxkddel

Double precision:: Pi_12, dPi12_ddel, d2Pi12_ddel2!, d3Pi12_ddel3
Double precision:: Pi_12_2, Pi_12_3, Pi_12_4

Double precision, dimension(30):: dPI12_dxi, d2PI12_dxiddel, d3PI12_dxiddel2!, d4PI12_dxiddel3, d5PI12_dxiddel4 
Double precision, dimension(:,:),allocatable:: d2PI12_dxidxj, d3PI12_dxidxjddel, d4PI12_dxidxjddel2
Double precision, dimension(:,:,:),allocatable :: d3PI12_dxidxjdxk, d4PI12_dxidxjdxkddel

Double precision:: help_A
Double precision, dimension(30):: dhelpA_dxi
Double precision, dimension(:,:),allocatable:: d2helpA_dxidxj
Double precision, dimension(:,:,:),allocatable:: d3helpA_dxidxjdxk

!Variables for the derivatives of parameter a wrt tau and xi
Double precision:: da_dtau, d2a_dtau2, d3a_dtau3, d4a_dtau4
Double precision, dimension(15):: dnadtaun
double precision, dimension(15,30)::dnadxidtaum       
double precision, dimension(:,:,:),allocatable ::dnadxi2dtaum
integer :: get_a_deriv

!Variables for the derivatives of the product tau*a wrt tau
Double precision:: dtaua_dtau, d2taua_dtau2, d3taua_dtau3, d4taua_dtau4

Double precision:: a_cubic, b_cubic, R_cubic
Double precision:: rhored_cubic, tred_cubic

Double precision:: tau2, tau3, tau4, del2, del3, del4        !Powers of delta and tau

double precision, dimension(30):: dbdxi_cubic
double precision, dimension(:,:),allocatable:: d2b_dxidxj
double precision, dimension(30):: da_dxi, d2a_dxidtau, d3a_dxidtau2, d4a_dxidtau3, d5a_dxidtau4
double precision, dimension(:,:),allocatable:: d2a_dxidxj, d3a_dxidxjdtau, d4a_dxidxjdtau2, d5a_dxidxjdtau3, d4a_dxidxjdtau4
double precision, dimension(:,:),allocatable:: d6a_dxidxjdtau4
double precision, dimension(30):: dtaua_dxi, d2taua_dxidtau, d3taua_dxidtau2, d4taua_dxidtau3, d5taua_dxidtau4
double precision, dimension(:,:),allocatable:: d2taua_dxidxj, d3taua_dxidxjdtau, d4taua_dxidxjdtau2, d5taua_dxidxjdtau3
double precision, dimension(:,:),allocatable:: d6taua_dxidxjdtau4

double precision :: help_c
double precision, dimension(30) :: dhelp_c_dxi
double precision, dimension(:,:),allocatable :: d2help_c_dxidxj
double precision, dimension(:,:,:),allocatable :: d3help_c_dxidxjdxk

integer:: i, j, k

!Supportive variables to speed up calculations
double precision:: denom_psim, denom_psim_2, denom_psim_3, denom_psim_4, denom_psim_5, denom_psim_6

!All CUBIC equations can be written in the same form:
!See Bell et al. (2016), "Helmholtz energy translations for common cubic equations of state for use in one-fluid and multi-fluid mixture models", to be submitted to ??
!$$ p = \frac{R*T}{v-b} - \frac{a}{(v+D1*b)*(v+D2*b)} $$
!For the different cubics, it is:
!   vdW: D1 = 0,        D2 = 0
!   SRK: D1 = 1,        D2 = 0
!   PR:  D1 = 1+2**0.5  D2 = 1 - 2**0.5
!The Helmholtz translation of the cubics can then also be written in a universal way:
!$$ \alpha_r = \psi^(-) - \frac{\tau*a}{R*T_r} * \psi^(+) $$
!Where:
!$$ \psi^(-) = -ln(1-b * \delta * \rho_r) $$
!$$ \psi^(+) = \frac{ln(\frac{D1*b*\rho_r*\delta}{D2*b*\rho_r*\delta})}{b*(D1-D2)} 
!The variable \psi^(-) is called psi_m in the code, this variable is only a function of delta
!The variable \psi^(+) is called psi_p in the code, this variable is only a function of delta

if (.not. allocated(d2psi_m_dxidxj)) then
    allocate(d2psi_m_dxidxj(30,30))
    allocate(d2psi_p_dxidxj, d3psi_p_dxidxjddel, d4psi_p_dxidxjddel2, d5psi_p_dxidxjddel3, d6psi_p_dxidxjddel4, d2PI12_dxidxj, d3PI12_dxidxjddel, d4PI12_dxidxjddel2,d2helpA_dxidxj, &
            d3psi_m_dxidxjddel, d4psi_m_dxidxjddel2, d5psi_m_dxidxjddel3, d6psi_m_dxidxjddel4,&
            d2b_dxidxj,d2a_dxidxj, d3a_dxidxjdtau, d4a_dxidxjdtau2, d5a_dxidxjdtau3, d4a_dxidxjdtau4,&
            d6a_dxidxjdtau4,d2taua_dxidxj, d3taua_dxidxjdtau, d4taua_dxidxjdtau2, d5taua_dxidxjdtau3, d6taua_dxidxjdtau4, d2help_c_dxidxj, mold=d2psi_m_dxidxj)
    allocate(d3psi_m_dxidxjdxk(30,30,30))
    allocate(d4psi_m_dxidxjdxkddel, d3psi_p_dxidxjdxk, d4psi_p_dxidxjdxkddel, d3PI12_dxidxjdxk, d4PI12_dxidxjdxkddel,d3helpA_dxidxjdxk,d3help_c_dxidxjdxk, mold=d3psi_m_dxidxjdxk)
    
    allocate(dnadxi2dtaum(15,30,30))
endif

If ((gl%Mix_type == 2) .or. (gl%Mix_type == 21) .or. (gl%Mix_type == 22)) then
    !SRK
    D1_cubic = 1.D0
    D2_cubic = 0.D0
    rhored_cubic = gl%rhored_SRK
    tred_cubic = gl%tred_SRK
    a_cubic = gl%a_SRK
    b_cubic = gl%b_SRK
    R_cubic = gl%R_SRK
    call db_SRK_dxi(gl,dbdxi_cubic)
    call d2b_SRK_dxidxj(gl,d2b_dxidxj)
elseif ((gl%Mix_type == 3) .or. (gl%Mix_type == 31)) then
    !PR
    D1_cubic = consA_PR
    D2_cubic = -consB_PR    !THE MINUS IS VERY IMPORTANT!! NOTE THE DIFFERENT DEFINITIONS OF THE CONSTANTS IN THE PAPER AND IN THE CODE!!!
    rhored_cubic = gl%rhored_PR
    tred_cubic = gl%tred_PR
    a_cubic = gl%a_PR
    b_cubic = gl%b_PR
    R_cubic = gl%R_PR
    call db_PR_dxi(gl,dbdxi_cubic)
    call d2b_PR_dxidxj(gl,d2b_dxidxj)
end if

delta = Density / rhored_cubic
tau = Tred_cubic / Temperature

!Variables for speeding up computations
tau2 = tau**2
tau3 = tau*tau2
tau4 = tau2*tau2
del2 = delta**2
del3 = del2 * delta
del4 = del2 * del2

!Calculate the derivatives of psi^(-) wrt xi
!psi_m = -dlog(1.D0-b_cubic*delta*rhored_cubic)
!Help variables:
denom_psim = (1.D0 - b_cubic * delta * rhored_cubic)
denom_psim_2 = denom_psim**2
denom_psim_3 = denom_psim_2 * denom_psim
denom_psim_4 = denom_psim_2 * denom_psim_2
denom_psim_5 = denom_psim_2 * denom_psim_3
denom_psim_6 = denom_psim_3 * denom_psim_3
dpsi_m_dxi = delta * dbdxi_cubic * rhored_cubic / denom_psim
d2psi_m_dxiddel = dbdxi_cubic * rhored_cubic / denom_psim_2
d3psi_m_dxiddel2 = 2.D0 * b_cubic* dbdxi_cubic * rhored_cubic**2 / denom_psim_3  
d4psi_m_dxiddel3 = 6.D0 * b_cubic**2 * dbdxi_cubic * rhored_cubic**3 / denom_psim_4
d5psi_m_dxiddel4 = 24.D0 * b_cubic**3 * dbdxi_cubic * rhored_cubic**4 / denom_psim_5

!Calculate the derivatives of psi^(-) wrt xi and xj
do i=1, gl%ncomp-1
    do j=1, gl%ncomp-1        
        if (i <= j) then!Speed up calculations, use symmetry
            d2psi_m_dxidxj(i,j) = delta * rhored_cubic * d2b_dxidxj(i,j) / denom_psim + &
                            & delta**2 * rhored_cubic**2 * dbdxi_cubic(i) * dbdxi_cubic(j) / denom_psim_2 
            d3psi_m_dxidxjddel(i,j) = rhored_cubic * d2b_dxidxj(i,j) / denom_psim_2 + &
                            & 2.D0 * delta * rhored_cubic**2 * dbdxi_cubic(i) * dbdxi_cubic(j) / denom_psim_3       
            d4psi_m_dxidxjddel2(i,j) = 2.D0 * rhored_cubic**2 * d2b_dxidxj(i,j) / denom_psim_3 + &
                            & 2.D0 * rhored_cubic**2 * dbdxi_cubic(i) * dbdxi_cubic(j) / denom_psim_4 * &
                            & (2.D0 * delta * rhored_cubic * b_cubic + 1.D0)
            d5psi_m_dxidxjddel3(i,j) = 6.D0 * rhored_cubic**3 * d2b_dxidxj(i,j) / denom_psim_4 + &
                            & 6.D0 * rhored_cubic**3 * b_cubic * dbdxi_cubic(i) * dbdxi_cubic(j) / denom_psim_5 * &
                            & (2.D0 * delta * rhored_cubic * b_cubic + 2.D0)
            d6psi_m_dxidxjddel4(i,j) = 24.D0 * rhored_cubic**4 * d2b_dxidxj(i,j) / denom_psim_5 + &
                            & 24.D0 * rhored_cubic**4 * b_cubic**2 * dbdxi_cubic(i) * dbdxi_cubic(j) / denom_psim_6 * &
                            & (2.D0 * delta * rhored_cubic * b_cubic + 3.D0)
        else 
            d2psi_m_dxidxj(i,j) = d2psi_m_dxidxj(j,i) 
            d3psi_m_dxidxjddel(i,j) = d3psi_m_dxidxjddel(j,i)
            d4psi_m_dxidxjddel2(i,j) = d4psi_m_dxidxjddel2(j,i) 
            d5psi_m_dxidxjddel3(i,j) = d5psi_m_dxidxjddel3(j,i)
            d6psi_m_dxidxjddel4(i,j) = d6psi_m_dxidxjddel4(j,i)
        end if
       end do
end do


!Calculate the derivatives of psi^(-) wrt xi, xj, and xk
do i=1, gl%ncomp-1
    do j=1, gl%ncomp-1     
        do k=1, gl%ncomp-1
                   
            d3psi_m_dxidxjdxk(i,j,k) = 2.D0 * delta**3 * rhored_cubic**3 * &
                                     & dbdxi_cubic(i) * dbdxi_cubic(j) * dbdxi_cubic(k) / denom_psim_3 + &
                                     & delta**2 * rhored_cubic**2 * (dbdxi_cubic(i) * d2b_dxidxj(j,k) + &
                                     & dbdxi_cubic(j) * d2b_dxidxj(i,k) + dbdxi_cubic(k) * d2b_dxidxj(i,j) )/ denom_psim_2 
            
            
            d4psi_m_dxidxjdxkddel(i,j,k) = 6.D0 * delta**2 * rhored_cubic**3 * &
                                     & dbdxi_cubic(i) * dbdxi_cubic(j) * dbdxi_cubic(k) / denom_psim_4 + &
                                     & 2.D0 * delta * rhored_cubic**2 * (dbdxi_cubic(i) * d2b_dxidxj(j,k) + &
                                     & dbdxi_cubic(j) * d2b_dxidxj(i,k) + dbdxi_cubic(k) * d2b_dxidxj(i,j) )/ denom_psim_3  
            
        end do
    end do
end do


!Calculate the derivatives of psi^(+) wrt delta, xi, and xj
!Calculate supporting variable Pi_12 and its derivatives wrt delta
!---
Pi_12 = (1.D0 + D1_cubic * b_cubic * rhored_cubic * delta) * (1.D0 + D2_cubic * b_cubic * rhored_cubic * delta)
Pi_12_2 = Pi_12 * Pi_12
Pi_12_3 = Pi_12 * Pi_12_2
Pi_12_4 = Pi_12_2 * Pi_12_2
dPi12_ddel = b_cubic * rhored_cubic * (2.D0 * D1_cubic * D2_cubic * b_cubic * delta * rhored_cubic + D1_cubic + D2_cubic)
d2Pi12_ddel2 = 2.D0 * D1_cubic * D2_cubic * b_cubic**2 * rhored_cubic**2
!d3Pi12_ddel3 = 0.D0
!---
!Calculate derivatives of supporting variable Pi_12 wrt xi and delta
!---
dPI12_dxi = delta * rhored_cubic * dbdxi_cubic * (D1_cubic * (D2_cubic * delta * rhored_cubic * b_cubic + 1.D0) & 
          & + D2_cubic * (D1_cubic * delta * rhored_cubic * b_cubic + 1.D0))
d2PI12_dxiddel = rhored_cubic * dbdxi_cubic * (D1_cubic * (D2_cubic * delta * rhored_cubic * b_cubic + 1.D0) & 
          & + D2_cubic * (D1_cubic * delta * rhored_cubic * b_cubic + 1.D0) & 
          & + 2.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * b_cubic)
d3PI12_dxiddel2 = 4.D0 * D1_cubic * D2_cubic * rhored_cubic**2 * b_cubic * dbdxi_cubic 
!---
!Calculate derivatives of supporting variable Pi_12 wrt xi, xj, and delta
!---
do i=1,gl%ncomp-1
    do j=1,gl%ncomp-1
        if (i <= j) then!Speed up calculations, use symmetry
            d2PI12_dxidxj(i,j) = delta * rhored_cubic * & 
                      & (2.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * dbdxi_cubic(i) * dbdxi_cubic(j) + &
                      & (2.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * b_cubic + D1_cubic + D2_cubic) * d2b_dxidxj(i,j))
            d3PI12_dxidxjddel(i,j) = rhored_cubic * & 
                      & (4.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * dbdxi_cubic(i) * dbdxi_cubic(j) + &
                      & (4.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * b_cubic + D1_cubic + D2_cubic) * d2b_dxidxj(i,j))
            d4PI12_dxidxjddel2(i,j) = 4.D0 * D1_cubic * D2_cubic * rhored_cubic**2 * &
                                & (b_cubic * d2b_dxidxj(i,j) + dbdxi_cubic(i) * dbdxi_cubic(j)) 
        else
            d2PI12_dxidxj(i,j) = d2PI12_dxidxj(j,i)
            d3PI12_dxidxjddel(i,j) = d3PI12_dxidxjddel(j,i)
            d4PI12_dxidxjddel2(i,j) = d4PI12_dxidxjddel2(j,i)
        end if
    end do
end do
!---
!Calculate derivatives of supporting variable Pi_12 wrt xi, xj, xk, and delta
!---
do i=1,gl%ncomp-1
    do j=1,gl%ncomp-1
        do k=1,gl%ncomp-1
            d3PI12_dxidxjdxk(i,j,k) = delta * rhored_cubic * & 
                      & (2.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * (dbdxi_cubic(i) * d2b_dxidxj(j,k) + & 
                      & dbdxi_cubic(j) * d2b_dxidxj(i,k) + dbdxi_cubic(k) * d2b_dxidxj(i,j)))
            d4PI12_dxidxjdxkddel(i,j,k) = rhored_cubic * & 
                      & (4.D0 * D1_cubic * D2_cubic * delta * rhored_cubic * (dbdxi_cubic(i) * d2b_dxidxj(j,k) + & 
                      & dbdxi_cubic(j) * d2b_dxidxj(i,k) + dbdxi_cubic(k) * d2b_dxidxj(i,j)))
        end do
    end do
end do
!---
!Calculate supporting variable A
!---
help_A = dlog((delta * rhored_cubic * b_cubic * D1_cubic + 1.D0) / (delta * rhored_cubic * b_cubic * D2_cubic + 1.D0))
dhelpA_dxi = delta * rhored_cubic * dbdxi_cubic * (D1_cubic - D2_cubic) / Pi_12
do i = 1, gl%ncomp-1
    do j=1,gl%ncomp-1
        if (i <= j) then!Speed up calculations, use symmetry
            d2helpA_dxidxj(i,j) = delta * rhored_cubic * (D1_cubic - D2_cubic) / Pi_12_2 * &
                            & (Pi_12 * d2b_dxidxj(i,j) - dPI12_dxi(j) * dbdxi_cubic(i))
        else
            d2helpA_dxidxj(i,j) = d2helpA_dxidxj(j,i)
        end if
    end do
end do
do i = 1, gl%ncomp-1
    do j=1, gl%ncomp-1
        do k=1, gl%ncomp-1
            d3helpA_dxidxjdxk(i,j,k) = delta * rhored_cubic * (D1_cubic - D2_cubic) / Pi_12_3 * &
                            & (-Pi_12 * (dPI12_dxi(j) * d2b_dxidxj(i,k) + dPI12_dxi(k) * d2b_dxidxj(i,j) &
                            & + dbdxi_cubic(i) * d2PI12_dxidxj(j,k)) + 2.D0 * dPI12_dxi(j) * dPI12_dxi(k) * dbdxi_cubic(i))
        end do
    end do
end do
!---
!Calculate supportung variable help_c
!---
help_c = 1.D0 / b_cubic
dhelp_c_dxi = -help_c**2 * dbdxi_cubic
do i = 1, gl%ncomp-1
    do j = 1, gl%ncomp-1
        d2help_c_dxidxj(i,j) = help_c**3 * (2.D0 * dbdxi_cubic(i)  * dbdxi_cubic(j) - b_cubic * d2b_dxidxj(i,j))     
        do k = 1, gl%ncomp-1
            d3help_c_dxidxjdxk(i,j,k) = help_c**4 * (2.D0 * b_cubic * (dbdxi_cubic(i) * d2b_dxidxj(j,k) + &
                                      & dbdxi_cubic(j) * d2b_dxidxj(i,k) + dbdxi_cubic(k) * d2b_dxidxj(i,j)) - &
                                      & 6.D0 * dbdxi_cubic(i) * dbdxi_cubic(j) * dbdxi_cubic(k))
        end do
    end do
end do
!---
!Calculate the delta derivatives of psi plus
psi_p = dlog((D1_cubic*b_cubic*rhored_cubic*delta+1.D0)/(D2_cubic*b_cubic*rhored_cubic*delta+1.D0)) &
      & /(b_cubic*(D1_cubic - D2_cubic))
dpsi_p_ddel = rhored_cubic / Pi_12
d2psi_p_ddel2 = -rhored_cubic * dPi12_ddel / Pi_12_2
d3psi_p_ddel3 = rhored_cubic * (-Pi_12 * d2Pi12_ddel2 + 2.D0 * dPi12_ddel**2) / Pi_12_3
d4psi_p_ddel4 = rhored_cubic * (6.D0 * Pi_12 * dPi12_ddel * d2Pi12_ddel2 - 6.D0 * dPi12_ddel**3) / Pi_12_4

!Calculate the xi and delta derivatives of psi plus
dpsi_p_dxi = (b_cubic * dhelpA_dxi - help_A * dbdxi_cubic) / ((D1_cubic - D2_cubic) * b_cubic**2) 
d2psi_p_dxiddel = -rhored_cubic / Pi_12_2 * dPI12_dxi
d3psi_p_dxiddel2 = -rhored_cubic / Pi_12_2 * (d2PI12_dxiddel + 2.D0 * Pi_12 / rhored_cubic * dPi12_ddel * d2psi_p_dxiddel)
d4psi_p_dxiddel3 = -1.D0 / Pi_12_2 * (rhored_cubic * d3PI12_dxiddel2 &
                 & + 2.D0 * d2psi_p_dxiddel * (dPi12_ddel**2 + Pi_12 * d2Pi12_ddel2) + &
                 & 4.D0 * Pi_12 * dPi12_ddel * d3psi_p_dxiddel2)
d5psi_p_dxiddel4 = -1.D0 / Pi_12_2 * (2.D0 * d2psi_p_dxiddel * 3.D0 * dPi12_ddel * d2Pi12_ddel2 + &
                 & 6.D0 * d3psi_p_dxiddel2 * (dPi12_ddel**2 + Pi_12 * d2Pi12_ddel2) + &
                 & 6.D0 * Pi_12 * dPi12_ddel * d4psi_p_dxiddel3)

!Calculate the xi, xj, and delta derivatives of psi plus
do i = 1, gl%ncomp-1
    do j=1,gl%ncomp-1
        if (i <= j) then!Speed up calculations, use symmetry
            d2psi_p_dxidxj(i,j) = (-b_cubic * (help_A * d2b_dxidxj(i,j) + dhelpA_dxi(i) * dbdxi_cubic(j) + & 
                            & dhelpA_dxi(j) * dbdxi_cubic(i)) + 2.D0 * help_A * dbdxi_cubic(i) * dbdxi_cubic(j) + &
                            & b_cubic**2 * d2helpA_dxidxj(i,j)) / ((D1_cubic - D2_cubic) * b_cubic**3)
            d3psi_p_dxidxjddel(i,j) = -1.D0 / Pi_12_2 * &
                                & (rhored_cubic * d2PI12_dxidxj(i,j) + 2.D0 * Pi_12 * dPI12_dxi(j) * d2psi_p_dxiddel(i))
            d4psi_p_dxidxjddel2(i,j) = -1.D0 / Pi_12_2 * (rhored_cubic * d3PI12_dxidxjddel(i,j) + &
                                 & 2.D0 * d2psi_p_dxiddel(i) * &
                                 & (Pi_12 * d2PI12_dxiddel(j) + dPi12_ddel * dPI12_dxi(j)) + 2.D0 * Pi_12 * dPi12_ddel * &
                                 & d3psi_p_dxidxjddel(i,j) + 2.D0 * Pi_12 * dPI12_dxi(j) * d3psi_p_dxiddel2(i))
            d5psi_p_dxidxjddel3(i,j) = -1.D0 / Pi_12_2 * (rhored_cubic * d4PI12_dxidxjddel2(i,j) + &
                                 & 2.d0 * d3psi_p_dxidxjddel(i,j) * (Pi_12 * d2Pi12_ddel2 + dPi12_ddel**2) + &
                                 & 4.D0 * d3psi_p_dxiddel2(i) * (Pi_12 * d2PI12_dxiddel(j) + dPI12_ddel * dPI12_dxi(j)) + &
                                 & 2.D0 * d2psi_p_dxiddel(i) * (Pi_12 * d3PI12_dxiddel2(j) + &
                                 & 2.D0 * dPi12_ddel * d2PI12_dxiddel(j) + dPI12_dxi(j) * d2Pi12_ddel2) + &
                                 & 4.D0 * Pi_12 * dPi12_ddel * d4psi_p_dxidxjddel2(i,j) + &
                                 & 2.D0 * Pi_12 * dPI12_dxi(j) * d4psi_p_dxiddel3(i))
            d6psi_p_dxidxjddel4(i,j) = -1.D0 / Pi_12_2 * &
                                 & (6.D0 * d4psi_p_dxidxjddel2(i,j) * (Pi_12 * d2Pi12_ddel2 + dPi12_ddel**2) + & 
                                 & 2.D0 * d3psi_p_dxidxjddel(i,j) * 3.D0 * dPi12_ddel * d2Pi12_ddel2 + &
                                 & 6.D0 * d4psi_p_dxiddel3(i) * (Pi_12 * d2PI12_dxiddel(j) + dPi12_ddel * dPi12_dxi(j)) +&
                                 & 6.D0 * d3psi_p_dxiddel2(i) * (Pi_12 * d3PI12_dxiddel2(j) + & 
                                 & 2.D0 * dPi12_ddel * d2PI12_dxiddel(j) + dPi12_dxi(j) * d2Pi12_ddel2) + &
                                 & 2.D0 * d2psi_p_dxiddel(i) * (3.D0 * dPi12_ddel * d3PI12_dxiddel2(j) + & 
                                 & 3.D0 * d2Pi12_ddel2 * d2PI12_dxiddel(j)) + &
                                 & 6.D0 * Pi_12 * dPi12_ddel * d5psi_p_dxidxjddel3(i,j) + &
                                 & 2.D0 * Pi_12 * dPi12_dxi(j) * d5psi_p_dxiddel4(i)) 
        else
            d2psi_p_dxidxj(i,j) = d2psi_p_dxidxj(j,i)
            d3psi_p_dxidxjddel(i,j) = d3psi_p_dxidxjddel(j,i) 
            d4psi_p_dxidxjddel2(i,j) = d4psi_p_dxidxjddel2(j,i)
            d5psi_p_dxidxjddel3(i,j) = d5psi_p_dxidxjddel3(j,i)
            d6psi_p_dxidxjddel4(i,j) = d6psi_p_dxidxjddel4(j,i)
        end if
    end do
end do

!Calculate the xi, xj, xk, and delta derivatives of psi plus
do i = 1, gl%ncomp-1
    do j = 1,gl%ncomp-1
        do k = 1, gl%ncomp-1
            d3psi_p_dxidxjdxk(i,j,k) = (help_A * d3help_c_dxidxjdxk(i,j,k) + help_c * d3helpA_dxidxjdxk(i,j,k) + &
                                     &  dhelpA_dxi(i) * d2help_c_dxidxj(j,k) + dhelpA_dxi(j) * d2help_c_dxidxj(i,k) + &
                                     &  dhelpA_dxi(k) * d2help_c_dxidxj(i,j) + dhelp_c_dxi(i) * d2helpA_dxidxj(j,k) + &
                                     &  dhelp_c_dxi(j) * d2helpA_dxidxj(i,k) + dhelp_c_dxi(k) * d2helpA_dxidxj(i,j)) / &
                                     &  (D1_cubic - D2_cubic)
            
            d4psi_p_dxidxjdxkddel(i,j,k) = -(rhored_cubic * d3PI12_dxidxjdxk(i,j,k) + 2.D0 * (Pi_12 * d2PI12_dxidxj(j,k) + &
                                         &  dPi12_dxi(j) * dPi12_dxi(k)) * d2psi_p_dxiddel(i) + & 
                                         &  2.D0 * Pi_12 * dPi12_dxi(j) * d3psi_p_dxidxjddel(i,k) + & 
                                         &  2.D0 * Pi_12 * dPi12_dxi(k) * d3psi_p_dxidxjddel(i,j)) / Pi_12_2   
    
        end do
    end do
end do

get_a_deriv = 4 !Get all derivatives of a with respect to tau
call da_dxi_dtau_cubic_all(gl,Temperature, density, get_a_deriv, dnadtaun, dnadxidtaum, dnadxi2dtaum)
!Derivatives of a with respect to tau
da_dtau = dnadtaun(1)
d2a_dtau2 = dnadtaun(2)
d3a_dtau3 = dnadtaun(3)
d4a_dtau4 = dnadtaun(4)
!Derivatives of a with respect to tau and xi
da_dxi = dnadxidtaum(1,:)
d2a_dxidtau = dnadxidtaum(2,:)
d3a_dxidtau2 = dnadxidtaum(3,:)
d4a_dxidtau3 = dnadxidtaum(4,:)
d5a_dxidtau4 = dnadxidtaum(5,:)
!Derivatives of a with respect to tau, xi, and xj
do i = 1, gl%ncomp-1
    do j = 1, gl%ncomp-1
        d2a_dxidxj(i,j) = dnadxi2dtaum(1,i,j)
        d3a_dxidxjdtau(i,j) = dnadxi2dtaum(2,i,j)
        d4a_dxidxjdtau2(i,j) = dnadxi2dtaum(3,i,j)
        d5a_dxidxjdtau3(i,j) = dnadxi2dtaum(4,i,j)
        d6a_dxidxjdtau4(i,j) = dnadxi2dtaum(5,i,j)
    end do
end do


!Calculate the derivatives of the product of (a * tau) with respect to tau
dtaua_dtau = tau * da_dtau + a_cubic
d2taua_dtau2 = tau * d2a_dtau2 + 2.D0 * da_dtau
d3taua_dtau3 = tau * d3a_dtau3 + 3.D0 * d2a_dtau2
d4taua_dtau4 = tau * d4a_dtau4 + 4.D0 * d3a_dtau3

!Calculate the derivatives of the product of (a * tau) with respect to tau and xi
dtaua_dxi = tau * da_dxi 
d2taua_dxidtau = tau * d2a_dxidtau + da_dxi
d3taua_dxidtau2 = tau * d3a_dxidtau2 + 2.D0 * d2a_dxidtau
d4taua_dxidtau3 = tau * d4a_dxidtau3 + 3.D0 * d3a_dxidtau2
d5taua_dxidtau4 = tau * d5a_dxidtau4 + 4.D0 * d4a_dxidtau3

!Calculate the derivatives of the product of (a * tau) with respect to tau, xi, and xj
do i = 1, gl%ncomp-1
    do j = 1, gl%ncomp-1
        d2taua_dxidxj(i,j) = tau * d2a_dxidxj(i,j) 
        d3taua_dxidxjdtau(i,j) = tau * d3a_dxidxjdtau(i,j) + d2a_dxidxj(i,j)
        d4taua_dxidxjdtau2(i,j) = tau * d4a_dxidxjdtau2(i,j) + 2.D0 * d3a_dxidxjdtau(i,j)
        d5taua_dxidxjdtau3(i,j) = tau * d5a_dxidxjdtau3(i,j) + 3.D0 * d4a_dxidxjdtau2(i,j)
        d6taua_dxidxjdtau4(i,j) = tau * d6a_dxidxjdtau4(i,j) + 4.D0 * d5a_dxidxjdtau3(i,j)
    end do
end do

!ALL DERIVATIVES ALSO WITH RESPECT TO xi and xj
do i =1, gl%ncomp-1
    do j =1, gl%ncomp-1
        do k = 1, gl%ncomp-1
        
        !alphar
        if (GETDER(1) == 1) then
            gl%MIXDERIVFNR_dxidxjdxk(1,i,j,k) = d3psi_m_dxidxjdxk(i,j,k) - 1.D0 / R_cubic / Tred_cubic * &
                                    & (d2taua_dxidxj(i,j) * dpsi_p_dxi(k) + 0.D0 * psi_p + & 
                                    &  dtaua_dxi(i) * d2psi_p_dxidxj(j,k) + d2taua_dxidxj(i,k) * dpsi_p_dxi(j) + &
                                    &  dtaua_dxi(j) * d2psi_p_dxidxj(i,k) + d2taua_dxidxj(j,k) * dpsi_p_dxi(i) + &
                                       tau * a_cubic * d3psi_p_dxidxjdxk(i,j,k) + dtaua_dxi(k) * d2psi_p_dxidxj(i,j)) 
        end if

        !d_alphar_d_delta * delta
        if (GETDER(2) == 1) then
            gl%MIXDERIVFNR_dxidxjdxk(2,i,j,k) = d4psi_m_dxidxjdxkddel(i,j,k) - 1.D0 / R_cubic / Tred_cubic * &
                                    & (d2taua_dxidxj(i,j) * d2psi_p_dxiddel(k) + 0.D0 * psi_p + & 
                                    &  dtaua_dxi(i) * d3psi_p_dxidxjddel(j,k) + d2taua_dxidxj(i,k) * d2psi_p_dxiddel(j) + &
                                    &  dtaua_dxi(j) * d3psi_p_dxidxjddel(i,k) + d2taua_dxidxj(j,k) * d2psi_p_dxiddel(i) + &
                                       tau * a_cubic * d4psi_p_dxidxjdxkddel(i,j,k) + dtaua_dxi(k) * d3psi_p_dxidxjddel(i,j))
            gl%MIXDERIVFNR_dxidxjdxk(2,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(2,i,j,k) * delta
        end if

        !d2_alphar_d_delta2 * delta^2
        if (GETDER(3) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(3,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(3,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(3,i,j,k) * del2
        end if

        !d_alphar_d_tau * tau
        if (GETDER(4) == 1) then
            gl%MIXDERIVFNR_dxidxjdxk(4,i,j,k) = -1.D0 / R_cubic / Tred_cubic * &
                                    & (d3taua_dxidxjdtau(i,j) * dpsi_p_dxi(k) + 0.D0 * psi_p + & 
                                    &  d2taua_dxidtau(i) * d2psi_p_dxidxj(j,k) + d3taua_dxidxjdtau(i,k) * dpsi_p_dxi(j) + &
                                    &  d2taua_dxidtau(j) * d2psi_p_dxidxj(i,k) + d3taua_dxidxjdtau(j,k) * dpsi_p_dxi(i) + &
                                       dtaua_dtau * d3psi_p_dxidxjdxk(i,j,k) + d2taua_dxidtau(k) * d2psi_p_dxidxj(i,j))
            gl%MIXDERIVFNR_dxidxjdxk(4,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(4,i,j,k) * tau
        end if

        !d2_alphar_d_tau2 * tau^2
        if (GETDER(5) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(5,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(5,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(5,i,j,k) * tau2
        end if

        !d2_alphar_d_delta_d_tau * delta * tau
        if (GETDER(6) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(6,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(6,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(6,i,j,k) * tau * delta
        end if

        !d3_alphar_d_delta_d_tau2 * delta * tau^2
        if (GETDER(7) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(7,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(7,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(7,i,j,k) * tau2 * delta
        end if

        !d3_alphar_d_delta3 * delta^3
        if (GETDER(8) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(8,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(8,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(8,i,j,k) * del3
        end if

        !d3_alphar_d_tau3 * tau^3
        if (GETDER(9) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(9,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(9,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(9,i,j,k) * tau3
        end if

        !d3_alphar_d_delta2_dtau * delta^2 * tau
        if (GETDER(10) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(10,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(10,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(10,i,j,k) * del2 * tau
        end if

        !d4_alphar_d_delta4 * delta^4
        if (GETDER(11) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(11,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(11,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(11,i,j,k) * del4
        end if

        !d4_alphar_d_delta3_dtau * delta^3*tau
        if (GETDER(12) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(12,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(12,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(12,i,j,k) * del3 * tau    
        end if

        !d4_alphar_d_delta2_dtau2 * delta^2*tau^2
        if (GETDER(13) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(13,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(13,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(13,i,j,k) * del2 * tau2
        end if

        !d4_alphar_d_delta_dtau3 * delta*tau^3
        if (GETDER(14) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(14,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(14,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(14,i,j,k) * delta * tau3
        end if

        !d4_alphar_d_tau4 * tau^4
        if (GETDER(15) == 1) then
            !NOT YET IMPLEMENTED!
            gl%MIXDERIVFNR_dxidxjdxk(15,i,j,k) = 0.D0
            gl%MIXDERIVFNR_dxidxjdxk(15,i,j,k) = gl%MIXDERIVFNR_dxidxjdxk(15,i,j,k) * tau4
        end if
    
        end do
    end do
end do

End subroutine MIXDERIVSFNR_dxidxjdxk_CUBIC



    end submodule impl
