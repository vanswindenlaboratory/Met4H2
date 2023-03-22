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

    ! module for file flash_pure.f90
    module flash_pure_module
    !global use inclusion
    use module_all_types
    use calc_functions
    use module_regula_falsi
    use module_regula_falsi_support
    use rhomix_pt_module
    use ancillary_equations_mix_module
    use fniderivs_module
    use fnrderivs_module
    use setup_module

    contains

	




!------------------------------------------------------------------------
!	Subroutine for calculation of VLE for pure fluids 
!	by iterative solution of the phase equilibrium condition 
!	in terms of Helmholtz energy 
!	(Thol, Wiens Bochum 05-2012)
!------------------------------------------------------------------------


subroutine Flash_Pure_PhaseBoundary (gl,press, Temp, rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst)
!DEC$ ATTRIBUTES DLLEXPORT :: Flash_Pure_PhaseBoundary
implicit none

    type(type_gl) :: gl


double precision :: press, Temp              ![MPa, K]
double precision :: rhovap_est, rholiq_est   ![mol/m³]
integer :: iFlash, nrsubst
integer :: errval, iter, n_call

! IN THIS ROUTINE THE FOLLOWING PHASE EQUILIBRIUM CALCULATIONS CAN BE PERFORMED:
!   - SATURATED TEMP: T  --  iFlash = 1
!   - VAPOR PRESSURE: P  --  iFlash = 2
!--------------------------------------------------------------------------

double precision :: dspinliq, dspinvap, pspinliq, pspinvap !P_CALC dspin_calc
double precision :: rhoredmix_orig, tredmix_orig, rhovapest, rholiqest

errval = 0
iter = 0
n_call = 1

rhoredmix_orig = gl%rhoredmix
tredmix_orig = gl%tredmix

gl%rhoredmix = gl%rhored(nrsubst)
gl%tredmix = gl%tred(nrsubst)


! check whether input is between triple and critical point, otherwise set error and return
if ((iFlash == 1)) then
    if(.not. gl%seawater) then
    if((Temp < gl%ttp(nrsubst)) .or. (Temp > gl%tc(nrsubst))) then 
        errval = -9981 
        return
        end if
    elseif (dabs(Temp - gl%tred(nrsubst)) <= 1.d-12) then !SH equality between double precision variables is always risky
        rhovap_est=gl%rhoredmix
        rholiq_est=gl%rhoredmix
        press=P_CALC(gl,Temp,gl%rhoredmix,nrsubst)
        return
    end if
elseif (iFlash == 2) then
    if((press <= gl%ptp(nrsubst)).or.(press > gl%pc(nrsubst))) then
        errval = -9982 
        return
    elseif (press == gl%pc(nrsubst)) then
        rhovap_est=gl%rhoredmix
        rholiq_est=gl%rhoredmix
        temp=gl%tredmix
        return
    end if
end if


! GENERATE STARTING VALUES

if (iFlash == 1) then
    !SATURATED VAPOR PRESSURE - estimate
    if (trim(gl%components(nrsubst)) /= "oil") then    !Special case if SRK for oil is used (Theresa), Andreas March 2013
        press = vp_eq(gl,Temp, nrsubst)
    else
        press = vp_oil(gl,Temp)
    end if
    !Added an error handling, if no (physically reasonable) vapor pressure was found
    !Andreas March 2013
    if (press <= 0.D0) then
        errval = -121212    
    end if
    !end if
else if (iFlash == 2) then
    !SATURATION TEMPERATURE - estimate
    if (Temp < 1.d-12) Temp = Estimate_Tsat(gl,press, nrsubst)
    if (Temp <= 0.d0) then !VORSICHT, FEHLER AUSGESCHALTET!!!
        errval = -121212
        return
    end if
else
    errval = -121213
    return
end if

!SATURATED LIQUID DENSITY - estimate
!Andreas Feb 2014
if (dabs(rholiq_est) < 1.d-12) then
    if (gl%eq_type(nrsubst) == 1) then
        rholiq_est = dl_eq(gl,Temp, nrsubst) 
        if (dabs(rholiq_est) < 1.d-12) then
            rholiq_est = rhomix_calc(gl,TEMP, press, 0.d0, 1, nrsubst)
        end if
    else    !Monika/Andreas Dez 2014. If no equation for saturated liquid densities is given, try using the density solver
        rholiq_est = rhomix_calc(gl,TEMP, press, 0.d0, 1, nrsubst)
        if ((trim(gl%components(nrsubst)) == "oil") .and. (rholiq_est < 300.d0)) rholiq_est = 465.5d0   !460.d0     
    end if
end if

!SATURATED VAPOR DENSITY - estimate
!Andreas Feb 2014
if (dabs(rhovap_est) < 1.d-12) then
    if (gl%eq_type(nrsubst) == 1) then
        rhovap_est = dv_eq(gl,Temp, nrsubst) 
        if (dabs(rhovap_est) < 1.d-12) then
            rhovap_est = rhomix_calc(gl,TEMP, press, 0.d0, 2, nrsubst)
        end if
    else !Monika/Andreas Dez 2014. If no equation for saturated vapor densities is given, try using the density solver
        rhovap_est = rhomix_calc(gl,TEMP, press, 0.d0, 2, nrsubst)
    end if
end if

if ((dabs(rholiq_est-rhovap_est) < 1.d-10) .or. (rholiq_est < rhovap_est)) then !pressure too low or too high so that only one solution was found !Theresa
!if (dabs(rholiq_est-rhovap_est) < 1.d-5) then !pressure too low or too high so that only one solution was found
    !get spinodals to set new pressure between min and max
    gl%startvaluespin=.true.
    dspinliq=DSPIN_CALC (gl,Temp,1,nrsubst)
    if (dspinliq == 0.d0) then
        errval = -2408
        return
    end if
    pspinliq=P_CALC(gl,Temp,dspinliq,nrsubst)
    dspinvap=DSPIN_CALC (gl,Temp,2,nrsubst)
    pspinvap=P_CALC(gl,Temp,dspinvap,nrsubst)
    press=(pspinliq+pspinvap)/2.d0
    if (press < 0.d0) press=pspinvap/2.d0 !liquid side probably far below zero
    !if (press < 0.d0) press = 1.d-10 !vapor side probaly not present or reasonable
    !recalculate new starting values for the densities from new pressure
    rholiq_est = rhomix_calc(gl,TEMP, press, 0.d0, 1, nrsubst)
    rhovap_est = rhomix_calc(gl,TEMP, press, 0.d0, 2, nrsubst)
    if (rhovap_est == 0.d0) rhovap_est = gl%rhoc(nrsubst)/rholiq_est

end if

if (rhovap_est < 1.d-15) then
    dspinvap = DSPIN_CALC (gl,Temp,2,nrsubst)
    rhovap_est = dspinvap*0.9d0
end if

rhovapest = rhovap_est
rholiqest = rholiq_est

call Flash_Pure_PhaseBoundary_calc(gl,press, Temp, rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst, n_call)



if (errval /= 0) then   !retry with same starting values but lower convergence criteria
    errval = 0
    call Flash_Pure_PhaseBoundary_calc(gl,press, Temp, rhovapest, rholiqest, iFlash, errval, iter, nrsubst, n_call)
end if
        

if (iflash == 1 .and. errval /= 0) then
    !generate new start values
    press = vp_eq(gl,Temp, nrsubst)
    rholiq_est = rhomix_calc(gl,TEMP, press, 0.d0, 1, nrsubst)!*1.05d0
    rhovap_est = rhomix_calc(gl,TEMP, press, 0.d0, 2, nrsubst)!*0.95d0
    if ((rholiq_est-rhovap_est) < 1.d-10 .or. gl%accen(nrsubst) == 0) then !pressure too low or too high so that only one solution was found
        !get spinodals to set new pressure between min and max
        gl%startvaluespin=.true.
        dspinliq=DSPIN_CALC (gl,Temp,1,nrsubst)
        pspinliq=P_CALC(gl,Temp,dspinliq,nrsubst)
        dspinvap=DSPIN_CALC (gl,Temp,2,nrsubst)
        pspinvap=P_CALC(gl,Temp,dspinvap,nrsubst)
        press=(pspinliq+pspinvap)/2.d0
        if (press < 0.d0) press = 1.d-10 !vapor side probaly not present or reasonable
        !recalculate new starting values for the densities from new pressure
        rholiq_est = rhomix_calc(gl,TEMP, press, 0.d0, 1, nrsubst)
        rhovap_est = rhomix_calc(gl,TEMP, press, 0.d0, 2, nrsubst)
        if (rholiq_est < 1.d-14) rholiq_est = dspinliq*1.1d0
        dspinvap=DSPIN_CALC (gl,Temp,2,nrsubst)
        if ((rholiq_est-rhovap_est) < 1.d-10) then
            rholiq_est = dspinliq*1.1d0
            rhovap_est = dspinvap*0.9d0
        end if
    end if
        
        
    if ((press > 1.d-14) .and. (rholiq_est > 1.d-14) .and. (rhovap_est > 1.d-14)) then
        errval=0
        call Flash_Pure_PhaseBoundary_calc(gl,press, Temp, rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst, n_call)
    else
        press=0.d0
        rhovap_est=0.d0
        rholiq_est=0.d0
    end if
end if

if (errval /= 0) then
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    return
else
    !calculate vapor pressure from given temperature and iterated density
    !if (iFlash == 1) press = P_CALC(gl,Temp, rholiq_est, nrsubst)
    if (iFlash == 1) press = P_CALC(gl,Temp, rhovap_est, nrsubst)
    
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
end if


end subroutine Flash_Pure_PhaseBoundary

    
subroutine MAXWELL (gl,ts,ps,dl,dv,error,nrsubst)
!



!
!  iterate for saturated liquid and vapor states of a pure fluid given temperature
!
!  inputs:
!    icomp--number of component
!       ts--temperature [K]
!       dl--initial guess for molar density [mol/L] of saturated liquid
!       dv--initial guess molar density [mol/L] of saturated vapor
!
!  outputs:
!       ps--saturation pressure [kPa]
!       dl--molar density [mol/L] of saturated liquid
!       dv--molar density [mol/L] of saturated vapor
!       error--error message
!



implicit none

    type(type_gl) :: gl

!
      double precision :: ts,ps,dl,dv
      integer :: nrsubst, error     
!
      double precision, parameter :: eps_conv=1.d-10
      integer, parameter :: kmax=200
!
      double precision :: tt,qt,vc,dtl,dtg,a1,a2,g1,g2,p1,p2,alpha1,alpha2 
      double precision :: vtl,vtg,v1,v2,v1ini,v2ini,vdif,pm,a,b,c,e,es,delv1,delv2 !p_calc
      integer :: k,is
     
      !get ideal Helmholtz Energy
      integer, dimension(nderivsi)::GETDERIVS
      double precision, dimension(nderivsi) :: SETDERIVS   
      !Get derivatives of the residual part of the Helmholtz energy needed for calculation dpdrho
      integer, dimension(nderivs):: GETDER
      double precision, dimension(nderivs)::FNRDER                   ! array with the computed values for the derivatives
      double precision:: dar_ddel, d2ar_ddel2
     
      !get residual Helmholtz Energy

      GETDERIVS=(/1,0,0,0,0,0,0,0,0,0/)
      GETDER=(/0,1,1,0,0,0,0,0,0,0,0,0,0,0,0/)   ! array specifier to indicate, which derivative is needed 

      
      if (gl%tminfluid(nrsubst)<gl%ttp(nrsubst)) then
        tt=gl%ttp(nrsubst)
      else
        tt=gl%tminfluid(nrsubst)
      endif
      
      
      qt = (ts-tt) / (gl%tc(nrsubst)-tt) !qt=1 for tc=tt
      vc = 1.d0/gl%rhoc(nrsubst)
      dtl = gl%rhotp(nrsubst)
!      dtg=ptp(nrsubst)/(Req(nrsubst)*ttp(nrsubst))*1.d6
      dtg = gl%ptp(nrsubst)/(gl%Req(nrsubst)*gl%ttp(nrsubst))*gl%factortrans
      vtl = 1.d0/dtl * 0.95d0   ! <=  vtr'
      vtg = 1.d0/dtg * 1.2d0    ! >=  vtr"

! --- start values for v1 and v2 -------------------------------------------
      v1 = 1.d0/dl
      v2 = 1.d0/dv
      v1ini = v1
      v2ini = v2
      if (v1ini<vtl) vtl = v1ini*0.95D0
      if (v2ini>vtg) vtg = v2ini*1.05D0
      k=0

! *** iteration algorithm ******************************************************
      do
        k = k+1
        !Moni: überflüssig?
        dl = 1.d0/v1
        dv = 1.d0/v2
        !Andreas, October 2014
        !For mixtures: nrsubst instead of fixed "1" for calculation needed!
        !p1=P_CALC(gl,ts,dl,1)
        p1 = P_CALC(gl,ts,dl,nrsubst)
        
        call FNIDERIVS(gl,ts,dl,GETDERIVS,SETDERIVS,nrsubst)
        !Andreas, October 2014
        !For mixtures: nrsubst instead of fixed "1" for calculation needed!
        !a1=(SETDERIVS(1)+A00_CALC(gl,ts,dl,1))*Req(nrsubst)*ts
        a1 = (SETDERIVS(1) + A00_CALC(gl,ts,dl,nrsubst) ) * gl%Req(nrsubst) * ts   !warum nicht einfach Helmholtz?
        
        !Andreas, October 2014: no option to calculate dp_drho for a specific component in the mixture. 
        !Replaced by other method
        !call dP_drho(ts,dl,alpha1)
        call FNRDERIVS (gl,ts, dl, GETDER, FNRDER, nrsubst)
        dar_ddel = FNRDER(2)
        d2ar_ddel2 = FNRDER(3)
        alpha1 = ts *  gl%Req(nrsubst) * (1.d0 +  2.D0*dar_ddel + d2ar_ddel2)     
      
        !Andreas, October 2014
        !For mixtures: nrsubst instead of fixed "1" for calculation needed!
        !p2=P_CALC(gl,ts,dv,1)
        p2 = P_CALC(gl,ts,dv,nrsubst)
        
        call FNIDERIVS(gl,ts,dv,GETDERIVS,SETDERIVS,nrsubst)
        !Andreas, October 2014
        !For mixtures: nrsubst instead of fixed "1" for calculation needed!
        !a2=(SETDERIVS(1)+A00_CALC(gl,ts,dv,1))*Req(nrsubst)*ts
        a2 = (SETDERIVS(1) + A00_CALC(gl,ts,dv,nrsubst)) * gl%Req(nrsubst) * ts
        
        !Andreas, October 2014: no option to calculate dp_drho for a specific component in the mixture. 
        !Replaced by other method
        !call dP_drho(ts,dv,alpha2)
        call FNRDERIVS (gl,ts, dv, GETDER, FNRDER, nrsubst)
        dar_ddel = FNRDER(2)
        d2ar_ddel2 = FNRDER(3)
        alpha2 = ts *  gl%Req(nrsubst) * (1.d0 +  2.D0*dar_ddel + d2ar_ddel2)
        
        !alpha1=-alpha1*dl*dl*1.D-6
        !alpha2=-alpha2*dv*dv*1.D-6
        
        alpha1 = -alpha1*dl*dl/gl%factortrans
        alpha2 = -alpha2*dv*dv/gl%factortrans

        is = 0
!        ierr1=0
!        if (alpha1>0.d0 .and. dlsp<0.d0) call SPNDL_L (t,x,dlsp,ierr1,herr)
        do while (alpha1>0.d0) !Maxwell loop (liquid)
          is = is+1
          if (is>10) goto 999
!          if (ierr1==0) then
!            v1=1.d0/dlsp*0.995d0
!          else
            v1 = (vtl+v1)*0.5d0
!          endif
          dl = 1.d0/v1
          !Andreas, October 2014: no option to calculate dp_drho for a specific component in the mixture. 
          !Replaced by other method
          !call dP_drho(ts,dl,alpha1)
          call FNRDERIVS (gl,ts, dl, GETDER, FNRDER, nrsubst)
          dar_ddel = FNRDER(2)
          d2ar_ddel2 = FNRDER(3)
          alpha1 = ts *  gl%Req(nrsubst) * (1.d0 +  2.D0*dar_ddel + d2ar_ddel2)           
          !alpha1=-alpha1*dl*dl*1.D-6
          alpha1 = -alpha1*dl*dl/gl%factortrans
          
          if (alpha1<0.d0) then
            !Andreas, October 2014
            !For mixtures: nrsubst instead of fixed "1" for calculation needed!
            !p1=P_CALC(gl,ts,dl,1)
            p1 = P_CALC(gl,ts,dl,nrsubst)
            
            call FNIDERIVS(gl,ts,dl,GETDERIVS,SETDERIVS,nrsubst)
            !Andreas, October 2014
            !For mixtures: nrsubst instead of fixed "1" for calculation needed!
            !a1=(SETDERIVS(1)+A00_CALC(gl,ts,dl,1))*Req(nrsubst)*ts
            a1 = (SETDERIVS(1) + A00_CALC(gl,ts,dl,nrsubst)) * gl%Req(nrsubst) * ts
            
          end if
        end do
        is = 0
!        if (alpha2>0.d0 .and. dvsp<0.d0) call SPNDL_V (t,x,dvsp,ierr1,herr)
        do while (alpha2>0.d0) !Maxwell loop (vapor)
          is = is+1
!          if (ierr1==0) then
!            v2=1.d0/dvsp*1.005d0
!          else
            if (qt>0.9d0) then
              v2 = v2*1.005d0
              if (v2>v2ini) v2 = v2ini
            else
              v2 = v2*1.2d0
            endif
            if (is>10) v2 = v2ini
!          endif
          if (is>20) goto 999
          dv = 1.d0/v2
          !Andreas, October 2014: no option to calculate dp_drho for a specific component in the mixture. 
          !Replaced by other method
          !call dP_drho(ts,dv,alpha2)
          call FNRDERIVS (gl,ts, dv, GETDER, FNRDER, nrsubst)
          dar_ddel = FNRDER(2)
          d2ar_ddel2 = FNRDER(3)
          alpha2 = ts *  gl%Req(nrsubst) * (1.d0 +  2.D0*dar_ddel + d2ar_ddel2)         
          !alpha2=-alpha2*dv*dv*1.D-6
          alpha2 = -alpha2*dv*dv/gl%factortrans
          
          if (alpha2<0.d0) then
            !Andreas, October 2014
            !For mixtures: nrsubst instead of fixed "1" for calculation needed!  
            !p2=P_CALC(gl,ts,dv,1)
            p2 = P_CALC(gl,ts,dv,nrsubst)
            
            call FNIDERIVS(gl,ts,dv,GETDERIVS,SETDERIVS,nrsubst)
            !Andreas, October 2014
            !For mixtures: nrsubst instead of fixed "1" for calculation needed!  
            !a2=(SETDERIVS(1)+A00_CALC(gl,ts,dv,1))*Req(nrsubst)*ts
            a2 = (SETDERIVS(1) + A00_CALC(gl,ts,dv,nrsubst)) * gl%Req(nrsubst) * ts
            
          end if
        end do

        vdif = v2-v1
        !pm=(a1-a2)/vdif*1.d-6
        pm = (a1-a2)/vdif/gl%factortrans
        if (dabs(vdif)<1.d-12 .or. k>kmax) goto 999
        if (dabs(p1-p2)/p2<eps_conv) goto 998

        a = 0.5d0*alpha1*(alpha1-alpha2)
        b = alpha1*(p1-p2-alpha2*(v1-v2))
        c = alpha2*(v1-v2)*(pm-p1)+0.5d0*(p1-p2)*(p1-p2)
        e = 0.25d0*b*b/(a*a)-c/a
        if (e<=1.d-10) then
          es = e*1.d12
          delv1 = -0.5d0*b/a+dsign(1.d-6*dsqrt(dabs(es)),(alpha1-alpha2)/alpha2)
        else
          delv1 = -0.5d0*b/a+dsign(dsqrt(dabs(e)),(alpha1-alpha2)/alpha2)
        end if
        delv2 = (p1-p2+alpha1*delv1)/alpha2

! --- calculate new values for v1 and v2 -----------------------------------
        v1 = v1+delv1
        v2 = v2+delv2
        if (v1>vc)  v1 = 0.9999d0*vc
        if (v1<vtl) v1 = vtl
        if (v2<vc)  v2 = 1.0001d0*vc
        if (v2>=vtg) v2 = vtg
        if (qt<0.9d0 .and. v1>v1ini*1.2d0) v1 = (vtl+v1ini)*0.5d0
      end do

 998  ps = pm
      return

999   continue
      error = -2212
      ps = -1.d0
      dl = -1.d0
      dv = -1.d0

end subroutine MAXWELL


subroutine Flash_Pure_PhaseBoundary_calc (gl,press, Temp, rhovap_est, rholiq_est, iFlash, errval, iter, nrsubst, n_call)


implicit none

    type(type_gl) :: gl


double precision :: press
double precision :: Temp
double precision :: rhovap_est, rholiq_est
integer :: iFlash, n_call
integer :: errval, iter

integer :: nrsubst, nEQN, k, i
double precision, dimension(3) :: var, var_est
double precision, dimension(3) :: EQN
double precision, dimension(3) :: Delta_D
double precision, dimension(3,3) :: JacMatrix
double precision :: rholiq_new, rhovap_new, max_del, Temp_new

!Exit criterion, Andreas July 2013
double precision:: eps_eqn, eps_eqn2, eps_eqn3,eps_eqn4,eps_p 

integer :: irow,icol
double precision :: det_j, mat_j_inv(3,3), mat(3,3), mat_diag(3,3)
double precision :: PLIQ, PVAP, PROZ, fac_del !P_CALC
double precision, dimension(100,3) :: dd

if (n_call == 1) then
    !define the standard exit criterion, Andreas July 2013
    eps_eqn = 1.D-14
    if (gl%mix_type == 6) eps_eqn = 1.d-15

    !define an alternate (softer) exit criterion in case the first criterion cannot be reached out of numeric reasons, Andreas July 2013
    eps_eqn4 = 1.D-10
    if (gl%mix_type == 6) eps_eqn4 = 5.d-8  
    eps_eqn2 = 1.D-8
    
    !define an alternate (even softer) exit criterion in case the alternate criterium is still too strict (In case of Hill's EOS for d2o necessary)
    eps_eqn3 = 1.D-6
    
    !define acceptable deviation of pvap and pliq
    eps_p = 1.D-6
else 
    !define the standard exit criterion, Andreas July 2013
    eps_eqn = 5.D-13
    if (gl%mix_type == 6) eps_eqn = 1.d-9

    !define an alternate (softer) exit criterion in case the first criterion cannot be reached out of numeric reasons, Andreas July 2013
    eps_eqn4 = 1.D-10
    if (gl%mix_type == 6) eps_eqn4 = 5.d-8  
    eps_eqn2 = 1.D-8
    
    !define an alternate (even softer) exit criterion in case the alternate criterium is still too strict (In case of Hill's EOS for d2o necessary)
    eps_eqn3 = 1.D-6
    
    !define acceptable deviation of pvap and pliq
    eps_p = 6.D-2
end if

n_call = n_call + 1

mat_diag= reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))

!If the relative change of the unknowns (T,p) is below eps_del, the calculation is finished

!Stefan Feb 2014
!if (eq_type(nrsubst) == 2) then ! in case of SRK there is no maximum density specified
!    rhomaxfluid(nrsubst) = 1.d+8
!end if

nEQN = 2
EQN = 0.d0
errval = 0
Delta_D = 0.d0
mat_j_inv = 0.d0
fac_del = 1.d0
dd = 0.d0


if (iFlash == 2) nEQN = 3

do iter = 1, 50 ! a maximum of 50 iterations allowed

    var_est = 0.d0
    var_est(1) = rholiq_est
    var_est(2) = rhovap_est
    var_est(3) = Temp


    !SysOfEqs: calculate system of equations (f'=f'', p'=p'', T'=T'')
    if (iFlash == 1) call SysOfEqs_T(gl, press, Temp, rhovap_est, rholiq_est, EQN, errval,nrsubst)
    if (iFlash == 2) call SysOfEqs_p_pure(gl, press, Temp, rhovap_est, rholiq_est, EQN, errval,nrsubst)
    
    if (errval /= 0) then
        return
    end if
    
    ! end calculation if a certain accuracy is reached
    if (iFlash == 1) then
        ! New exit criterion, Andreas July 2013
        if (maxval(dabs(EQN)) < eps_eqn) then
            !safety caution for numerical effects (see oil, Pvap is 1.d-10 --> very small steps)
            PLIQ = P_CALC(gl,Temp, rholiq_est, nrsubst)
            PVAP = P_CALC(gl,Temp, rhoVAP_est, nrsubst)
            PROZ = DABS((PLIQ-PVAP)/PLIQ*100.D0)
            if (PROZ < eps_p) return !saft/oil?
            !if (PROZ < 1.d-6) return
            
            !change to steps in pressure instead of density
            press = (PLIQ + PVAP) * 0.5d0
            if (press < 0.d0) press = 1.d-12
            call flash_step_p(gl,TEMP,press,rholiq_est,rhovap_est,nrsubst,errval)
            return

        end if
    else if (iFlash == 2)then
        ! New exit criterion, Andreas July 2013
        if (maxval(dabs(EQN)) < eps_eqn) return
    end if
    
    
    Delta_D = -EQN

    ! calculation of the derivatives needed for the solution of the system of equations
    if (iFlash == 1) call Jacobi_Pure_T(gl, press, temp, rhovap_est, rholiq_est, iFlash, JacMatrix, errval,nrsubst)
    if (iFlash == 2) call Jacobi_Pure_p(gl, press, temp, rhovap_est, rholiq_est, iFlash, JacMatrix, errval,nrsubst)
    
    if (errval /= 0) then
        return
    end if
    
    if (iFlash == 1) det_j = DET2X2(gl,JacMatrix)
    if (iFlash == 2) det_j = DET3X3(gl,JacMatrix)

    
    do icol=1,iFlash+1
      do irow=1,iFlash+1
        mat(:,:)=JacMatrix(:,:)
        mat(:,irow)=mat_diag(:,icol)
        if (iFlash == 1) mat_j_inv(irow,icol)=DET2X2(gl,mat)/det_j
        if (iFlash == 2) mat_j_inv(irow,icol)=DET3X3(gl,mat)/det_j
      enddo
    enddo

    Delta_D=matmul(Delta_D,mat_j_inv)!*0.5d0
    dd(iter,:) = Delta_D
    
    if (iter == 6 .and. gl%mix_type == 6) then
        if ((dabs((dd(6,1)-dd(4,1))/dd(6,1)) < 30.d-2 .and. dabs(dd(6,1)) > 1.d-2 .and. dd(6,1)*dd(5,1) < 0.d0) &
            & .or. dabs(dd(5,1)) > dabs(dd(1,1)) .or. dabs(dd(6,1)) > dabs(dd(1,1))) then
            fac_del = 0.5d0
        end if
    end if

    ! set the new values
    rholiq_new = rholiq_est + Delta_D(1)*fac_del
    rhovap_new = rhovap_est + Delta_D(2)*fac_del
    if (iFlash == 2) Temp_new = Temp + Delta_D(3)*fac_del
    
    
    var = 0.d0
    var(1) = rholiq_new
    var(2) = rhovap_new
    if (iFlash == 2) var(3) = Temp_new
    
    

    ! safety cautions ------------------------------------------------------
    do while ((var(1) < gl%rhoc(nrsubst).or.(var(1) > gl%rhomaxfluid(nrsubst))) .or. &
            & (var(2) < 0.d0).or.(var(2) > gl%rhoc(nrsubst)))
        do i = 1, nEQN
            Delta_D(i) = Delta_D(i)*0.5d0
            var(i) = var_est(i) + Delta_D(i)
        end do
        if (maxval(dabs(Delta_D)) < 1.D-17) then
        !if (maxval(dabs(Delta_D)) < 1.D-10) then !was not precise enough for low pressure oil
            errval = -898964
            
            PLIQ = P_CALC(gl,Temp, rholiq_est, nrsubst)
            PVAP = P_CALC(gl,Temp, rhoVAP_est, nrsubst)
            PROZ = DABS((PLIQ-PVAP)/PLIQ*100.D0)
            if (PROZ < eps_p) return
            
            !change to steps in pressure instead of density
            if (iFlash == 1) then
                press = (PLIQ + PVAP) * 0.5d0
                if (press < 0.d0) press = 1.d-12
                call flash_step_p(gl,TEMP,press,rholiq_est,rhovap_est,nrsubst,errval)
            end if
            return
        end if
    end do

    !Andreas, Erik, July 2019
    !For other EoS (for example cubic EOS) the calculated Ttp can be smaller than the tabulated values
    if (iFlash == 2) then
        !do while ((var(3) > gl%tc(nrsubst)) .or. (var(3) < gl%ttp(nrsubst)))
        do while ((var(3) > gl%tc(nrsubst)) .or. (var(3) < (gl%ttp(nrsubst) - 2.D0)))
            do i = 1, nEQN
                Delta_D(i) = Delta_D(i)*0.5d0
                var(i) = var_est(i) + Delta_D(i)
            end do
            if (maxval(dabs(Delta_D)) < 1.D-11) then
            !if (maxval(dabs(Delta_D)) < 1.D-10) then
                errval = -898964
                return
            end if
        end do
    end if


    ! END safety cautions -------------------------------------------

    rholiq_est = var(1)
    rhovap_est = var(2)
    if (iFlash == 2) Temp = var(3)

    

    !Second exit criterion: If the maximum relative change of the variables is lower than eps_del, the algorithm converged
    max_del = 0.D0
    Do k = 1, nEQN   
        if(dabs(Delta_D(k) / var(k)) > max_del) then
            max_del = dabs(Delta_D(k) / var(k))
        end if
    end do
    
    if ((max_del) < eps_eqn) then  !14      
        !safety caution for numerical effects (see oil, Pvap is 1.d-10 --> very small steps)
        PLIQ = P_CALC(gl,Temp, rholiq_est, nrsubst)
        PVAP = P_CALC(gl,Temp, rhoVAP_est, nrsubst)
        !if ((PLIQ < 0.d0) .or. (PVAP < 0.d0)) then
        !    errval = -898964
        !endif
        PROZ = DABS((PLIQ-PVAP)/PLIQ*100.D0)
        if (PROZ < eps_p) return 
        
        !change to steps in pressure instead of density
        if (iFlash == 1) then
            press = (PLIQ + PVAP) * 0.5d0
            if (press < 0.d0) press = 1.d-12
            call flash_step_p(gl,TEMP,press,rholiq_est,rhovap_est,nrsubst,errval)
            return
        end if
    end if
    
    if ((iFlash == 2) .and. (press/gl%pc(nrsubst) > 0.99999999d0) .and. (max_del < 6.d-6)) then
        return
    end if


end do



!New Version, Andreas July 2013 (the first criterium remained unchanged, but is in principle not necessary)
if (((max_del) < eps_eqn4) .and. (maxval(dabs(EQN)) < eps_eqn2)) then
    !safety caution for numerical effects (see oil, Pvap is 1.d-10 --> very small steps)
    if (iFlash == 1) then
    
        PLIQ = P_CALC(gl,Temp, rholiq_est, nrsubst)
        PVAP = P_CALC(gl,Temp, rhoVAP_est, nrsubst)
        PROZ = DABS((PLIQ-PVAP)/PLIQ*100.D0)
        if (PROZ < eps_p) return
    
        !change to steps in pressure instead of density
        press = (PLIQ + PVAP) * 0.5d0
        if (press < 0.d0) press = 1.d-12
        call flash_step_p(gl,TEMP,press,rholiq_est,rhovap_est,nrsubst,errval)
    end if
    return

elseif (maxval(dabs(EQN)) < eps_eqn3) then
    !This criterium is dangerously implemented at the moment, since the user does not get any indication which convergence criterium was fulfilled
    !However, in most cases it should be sufficient to ensure that the residua are smaller than 10^-6
        !safety caution for numerical effects (see oil, Pvap is 1.d-10 --> very small steps)
    PLIQ = P_CALC(gl,Temp, rholiq_est, nrsubst)
    PVAP = P_CALC(gl,Temp, rhoVAP_est, nrsubst)
    PROZ = DABS((PLIQ-PVAP)/PLIQ*100.D0)
    if (PROZ < eps_p) return
    
    if (iFlash == 1) then
        !change to steps in pressure instead of density
        press = (PLIQ + PVAP) * 0.5d0
        if (press < 0.d0) press = 1.d-12
        call flash_step_p(gl,TEMP,press,rholiq_est,rhovap_est,nrsubst,errval)
    end if
    return

else
    errval = -898965
    return
end if  

    end subroutine Flash_Pure_PhaseBoundary_calc
    
    
    subroutine flash_step_p(gl,TEMP,press,rholiq_est,rhovap_est,nrsubst,errval)
    






implicit none

    type(type_gl) :: gl

    integer :: errval, nrsubst, i, jj
    double precision :: TEMP,press,rholiq_est,rhovap_est,gliq,gvap,delpp,delp2,eps !G_CALC
    
    if (.not.gl%cpmodel(nrsubst)) then
        errval = -898966
        return
    end if
    
    eps = 1.d-9
    if (gl%mix_type == 6)  eps = 1.d-8
    delp2 = 1.d0
      
    rholiq_est = rhomix_calc(gl,TEMP, press, rholiq_est, 1, nrsubst)
    rhovap_est = rhomix_calc(gl,TEMP, press, 0.d0, 2, nrsubst)
    if (rholiq_est == 0.d0 .or. rhovap_est == 0.d0) then
        errval = -898966
        return
    end if
    if (.NOT.gl%ref_set) call ref_state(gl,errval) !Gibbs Enthalpy will be used 
    
    do i = 1, 40
        gliq = G_CALC(gl,Temp, rholiq_est, nrsubst)
        gvap = G_CALC(gl,Temp, rhovap_est, nrsubst)
        delpp = (gliq-gvap) / (1.d0/rholiq_est-1.d0/rhovap_est) *1.d-6
        if (DABS(delpp)>0.4d0*press) then    ! reducing step size if necessary
            do jj=1,20
                delpp=0.5*delpp
                if (DABS(delpp)<0.4d0*press) exit
            enddo
        end if  
        press=press-delpp
        rholiq_est = rhomix_calc(gl,TEMP, press, rholiq_est, 1, nrsubst)
        rhovap_est = rhomix_calc(gl,TEMP, press, 0.d0, 2, nrsubst)
        if (rholiq_est == 0.d0) rholiq_est = rhomix_calc(gl,TEMP, press, 0.d0, 1, nrsubst)
        if (rholiq_est == 0.d0 .or. rhovap_est == 0.d0) then
            errval = -898966
            return
        end if
        
        !exit criterion
        if (dabs(delpp/press)<eps .or. dabs((delpp-delp2)/delp2)<1.d-11) then
            return
        end if
        if (press < 1.d-16) then
                if (dabs(delpp/press)<0.3d0 .or. dabs((delpp-delp2)/delp2)<1.d-11) then
                    return
                end if
        end if
        delp2 = delpp
    end do
    errval = -898966
        
    end subroutine flash_step_p

    
    
    !*************************************************************************************************************************
    !*************************************************************************************************************************
    !*************************************************************************************************************************
double precision function DET3X3(gl,mat)

implicit none

    type(type_gl) :: gl

    
    double precision mat(3,3)

    DET3X3=mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)) &
          -mat(2,1)*(mat(1,2)*mat(3,3)-mat(3,2)*mat(1,3)) &
          +mat(3,1)*(mat(1,2)*mat(2,3)-mat(2,2)*mat(1,3))

end function DET3X3
    
    
double precision function DET2X2(gl,mat)

implicit none

    type(type_gl) :: gl


    double precision mat(3,3)

    DET2X2=mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2)

end function DET2X2
    
    
    
    !*************************************************************************************************************************
    !*************************************************************************************************************************
    !*************************************************************************************************************************


subroutine SysOfEqs_T (gl,press, Temp, rhovap_est, rholiq_est, EQN, errval,nrsubst)

! SysOfEqs: calculate system of equations (f'=f'', p'=p'')




implicit none

    type(type_gl) :: gl


    double precision :: press, Temp !, FUGCOPURE_CALC
    double precision :: rhovap_est, rholiq_est
    integer :: errval

    integer :: nrsubst
    double precision, dimension(3) :: EQN
    double precision :: press_vap, press_liq, fugcoef_pure_liq,fugcoef_pure_vap !P_CALC
   
    ! get fugacity coefficients, calculate p' and p''
    fugcoef_pure_liq=FUGCOPURE_CALC(gl,Temp,rholiq_est,nrsubst)
    if (errval /= 0) return
    press_liq = P_CALC(gl,Temp, rholiq_est, nrsubst)
    
    fugcoef_pure_vap=FUGCOPURE_CALC(gl,Temp,rhovap_est, nrsubst)
    if (errval /= 0) return
    press_vap = P_CALC(gl,Temp, rhovap_est, nrsubst)
    

    EQN(1) = press_liq*fugcoef_pure_liq - press_vap*fugcoef_pure_vap ! f'-f''
    EQN(2) = press_liq - press_vap                                         ! p'-p''


end subroutine SysOfEqs_T



!subroutine SysOfEqs_p (press, Temp, rhovap_est, rholiq_est, EQN, errval,nrsubst)
!
!    implicit none
!
!    double precision :: press, Temp
!    double precision :: rhovap_est, rholiq_est
!    integer :: errval
!
!    integer :: nrsubst
!    double precision, dimension(3) :: EQN
!    double precision :: fugcoef_pure_liq, fugcoef_pure_vap, FUGCOPURE_CALC
!    double precision :: P_CALC, press_liq, press_vap
!   
!    ! get fugacity coefficients, calculate p' and p''
!    fugcoef_pure_liq=FUGCOPURE_CALC(gl,Temp,rholiq_est, nrsubst)
!    press_liq = P_CALC(gl,Temp, rholiq_est, nrsubst)
!    
!    fugcoef_pure_vap=FUGCOPURE_CALC(gl,Temp,rhovap_est, nrsubst)
!    press_vap = P_CALC(gl,Temp, rhovap_est, nrsubst)
!    
!    if (errval /= 0) return
!    
!    
!    EQN(1) = press_liq*fugcoef_pure_liq - press_vap*fugcoef_pure_vap  ! f'-f''
!    EQN(2) = press - press_vap                                              ! psat-p''calc
!    EQN(3) = press - press_liq                                              ! psat-p'calc
!
!
!end subroutine SysOfEqs_p

subroutine SysOfEqs_p_PURE (gl,press, Temp, rhovap_est, rholiq_est, EQN, errval,nrsubst)




implicit none

    type(type_gl) :: gl


    double precision :: press, Temp
    double precision :: rhovap_est, rholiq_est
    integer :: errval

    integer :: nrsubst
    double precision, dimension(3) :: EQN
    double precision :: fugcoef_pure_liq, fugcoef_pure_vap !, FUGCOPURE_CALC
    double precision :: press_liq, press_vap !P_CALC
   
    ! get fugacity coefficients, calculate p' and p''
    fugcoef_pure_liq=FUGCOPURE_CALC(gl,Temp,rholiq_est, nrsubst)
    press_liq = P_CALC(gl,Temp, rholiq_est, nrsubst)
    
    fugcoef_pure_vap=FUGCOPURE_CALC(gl,Temp,rhovap_est, nrsubst)
    press_vap = P_CALC(gl,Temp, rhovap_est, nrsubst)
    
    if (errval /= 0) return
    
    
    EQN(1) = press_liq*fugcoef_pure_liq - press_vap*fugcoef_pure_vap  ! f'-f''
    EQN(2) = press - press_vap                                              ! psat-p''calc
    EQN(3) = press - press_liq                                              ! psat-p'calc


end subroutine SysOfEqs_p_PURE

subroutine Jacobi_Pure_T (gl,press, Temp, rhovap_est, rholiq_est, iFlash, JacMatrix, errval,nrsubst)








implicit none

    type(type_gl) :: gl


    double precision :: press, Temp
    double precision :: rhovap_est, rholiq_est
    integer :: iFlash
    integer :: errval

    integer :: nrsubst
    double precision, dimension(3,3) :: JacMatrix
    double precision :: press_liq, press_vap ! ,FUGCOPURE_CALC , P_CALC, DPDD_CALC
    double precision :: dpdd_liq, dpdd_vap, fugcoef_pure_liq, fugcoef_pure_vap
!    double precision:: DFugcoefpureDD_CALC
    !double precision,dimension(30) ::  dphidd_liq, dphidd_vap
    double precision::  dphidd_liq, dphidd_vap

    JacMatrix = 0.d0
    errval = 0
    fugcoef_pure_liq = 0.d0
    fugcoef_pure_vap = 0.d0

    ! get fugacity coefficients, calculate p' and p'', calculate dp/dd' and dp/dd''
    fugcoef_pure_liq=FUGCOPURE_CALC(gl,Temp,rholiq_est, nrsubst)
    dphidd_liq=DFugcoefpureDD_CALC(gl,Temp,rholiq_est, nrsubst)
    press_liq = P_CALC(gl,Temp, rholiq_est, nrsubst)
    dpdd_liq = DPDD_CALC(gl,Temp, rholiq_est, nrsubst)
    
    fugcoef_pure_vap=FUGCOPURE_CALC(gl,Temp,rhovap_est, nrsubst)
    dphidd_vap=DFugcoefpureDD_CALC(gl,Temp,rhovap_est, nrsubst)
    press_vap = P_CALC(gl,Temp, rhovap_est, nrsubst)
    dpdd_vap = DPDD_CALC(gl,Temp, rhovap_est, nrsubst)
    
    if (errval /= 0) return

    ! set Jacobi Matrix with the derivatives needed
    JacMatrix(1,1) = dphidd_liq*press_liq + dpdd_liq*fugcoef_pure_liq
    JacMatrix(2,1) = -(dphidd_vap*press_vap + dpdd_vap*fugcoef_pure_vap)
    JacMatrix(1,2) = DPDD_CALC(gl,Temp, rholiq_est, nrsubst)
    JacMatrix(2,2) = -DPDD_CALC(gl,Temp, rhovap_est, nrsubst)

end subroutine Jacobi_Pure_T


subroutine Jacobi_Pure_p (gl,press, Temp, rhovap_est, rholiq_est, iFlash, JacMatrix, errval, nrsubst)








implicit none

    type(type_gl) :: gl


    double precision :: press, Temp
    double precision :: rhovap_est, rholiq_est
    integer :: iFlash
    integer :: errval

    integer :: nrsubst
    double precision, dimension(3,3) :: JacMatrix
!    double precision :: DFugcoef_pureDT ! DPDD_CALC, DPDT_CALC,
    double precision, dimension(30) :: dphidd_liq, dphidd_vap
    double precision :: press_liq, press_vap, dpdd_liq, dpdd_vap, dpdt_liq, dpdt_vap
    double precision :: dphidT_liq, dphidT_vap, fugcoef_pure_liq, fugcoef_pure_vap !, FUGCOPURE_CALC,P_CALC
   ! double precision:: DFugcoefpureDD_CALC
    
    
    JacMatrix = 0.d0
    errval = 0

    ! get fugacity coefficients, dphi/dd' and dphi/dd'', dphi/dT' and dphi/dT'', p' and p'', dp/dd' and dp/dd'', calculate dp/dT' and dp/dT''
    fugcoef_pure_liq=FUGCOPURE_CALC(gl,Temp,rholiq_est, nrsubst)
    dphidd_liq= DFugcoefpureDD_CALC(gl,Temp,rholiq_est, nrsubst)
    dphidT_liq = DFugcoef_pureDT (gl,Temp, rholiq_est, nrsubst)
    press_liq = P_CALC(gl,Temp, rholiq_est, nrsubst)
    dpdd_liq = DPDD_CALC(gl,Temp, rholiq_est, nrsubst)
    dpdt_liq = DPDT_CALC(gl,Temp,rholiq_est, nrsubst)
    
    fugcoef_pure_vap=FUGCOPURE_CALC(gl,Temp,rhovap_est, nrsubst)
    dphidd_vap=DFugcoefpureDD_CALC(gl,Temp,rhovap_est, nrsubst)
    dphidT_vap = DFugcoef_pureDT (gl,Temp, rhovap_est, nrsubst)
    press_vap = P_CALC(gl,Temp, rhovap_est, nrsubst)
    dpdd_vap = DPDD_CALC(gl,Temp, rhovap_est, nrsubst)
    dpdt_vap = DPDT_CALC(gl,Temp,rhovap_est, nrsubst)
   
    if (errval /= 0) return
    
    ! set Jacobi Matrix with the derivatives needed
    JacMatrix(1,1) = dphidd_liq(nrsubst)*press_liq + dpdd_liq*fugcoef_pure_liq
    JacMatrix(2,1) = -(dphidd_vap(nrsubst)*press_vap + dpdd_vap*fugcoef_pure_vap)
    JacMatrix(3,1) = dphidT_liq*press_liq + dpdt_liq*fugcoef_pure_liq &
                    & -(dphidT_vap*press_vap + dpdt_vap*fugcoef_pure_vap)
    JacMatrix(1,2) = 0.d0
    JacMatrix(2,2) = - dpdd_vap
    JacMatrix(3,2) = - dpdt_vap
    JacMatrix(1,3) = - dpdd_liq
    JacMatrix(2,3) = 0.d0
    JacMatrix(3,3) = - dpdt_liq

end subroutine Jacobi_Pure_p









!--------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION Estimate_Tsat(gl,P, nrsubst)
!function calculating the estimate for Tsat
!calls regula falsi on Res_Pdiff_Tsat=P-Psat(Tsat-est)
!(J.P.Jakobsen, Denmark 09-2009)
!--------------------------------------------------------------------------
!USE module_fluid_parameters


implicit none

    type(type_gl) :: gl
  

DOUBLE PRECISION :: Tsat,P
DOUBLE PRECISION :: theta,Tsat_est
DOUBLE PRECISION :: Tstart_min,Tstart_max				!variables for regula falsi
DOUBLE PRECISION :: Tmin_allowed,Tmax_allowed
DOUBLE PRECISION :: Delta_allowed 
integer  ::Max_Iterations,Iterations,errval, i 
integer :: nrsubst


type(type_additional_parameters):: Res_Pdiff_param
!Res_Pdiff_param = 0.d0
!--------------------------------------------------------------------------
!Strating value
Estimate_Tsat=0.d0
i = nrsubst

!Specify parameters for regula falsi
Tmin_allowed = gl%ttp(i)
Tmax_allowed = gl%tc(i)
Delta_allowed = 1.D-6
!Closer and closer to the critical point
!the strating interval for regula falsi is gettin more narrow
!delta_allowed must be decreased in order for the routine to itterate at all
!IF( 0.99d0 < P/pc(i)) Delta_allowed = 1.D-8
!IF( 0.999d0 < P/pc(i)) Delta_allowed = 1.D-10  
!IF( 0.9999d0 < P/pc(i)) Delta_allowed = 1.D-12
IF( 0.9999d0 < P/gl%pc(i)) then 
    Delta_allowed = 1.D-12
ELSEIF( 0.99d0 < P/gl%pc(i)) then 
    Delta_allowed = 1.D-8
ELSEIF( 0.999d0 < P/gl%pc(i)) then 
    Delta_allowed = 1.D-10  
END IF
Max_iterations = 30
Res_Pdiff_param%a_p(1)=P
Res_Pdiff_param%a_p(2)=i
Res_Pdiff_param%a_p(65) = -999.D0 ! Indicate that Regula Falsi is called by the Maxwell iteration
 
!------------------------------------------------------------------------
!Estimate Tsat

Tsat_est=gl%tc(i)/(1.d0-0.428571d0*log10(P/gl%pc(i))/(1.d0+gl%accen(i)))

IF (P > 0.98d0*gl%pc(i)) THEN
!  initial guesses for density;
!  using correlation developed by E.W. Lemmon, NIST
   theta=(1.d0-P/gl%pc(i))*(1.6d0-gl%accen(i))
   Tsat_est=gl%tc(i)*(1.0d0-0.103947d0*theta-4.108265d-2*theta**2)
END IF

Tstart_min = 0.96d0*Tsat_est
Tstart_max = 1.04d0*Tsat_est
IF (Tstart_max > gl%tc(i)-1.d-4) Tstart_max = gl%tc(i)-1.d-4

if (Tstart_min < Tmin_allowed) Tstart_min = Tmin_allowed

!Iterate on Tsat
CALL Regula_Falsi(gl,Res_Pdiff_Tsat,Tsat,Tstart_min,Tstart_max,Delta_allowed,&
Tmin_allowed,Tmax_allowed, Max_iterations,Iterations, errval,Res_Pdiff_param)

Estimate_Tsat=Tsat    


!VORSICHT, FEHLER AUSGESCHALTET!!!
IF ((errval == 1).or.(errval == 1))THEN ! hier nachbessern
  !PRINT*,"Estimate for Tsat was not found by iteration"
  !PRINT*,"Tsat_est used from correlation developed by E.W. Lemmon, NIST"
  Estimate_Tsat=Tsat_est
else
    Estimate_Tsat=gl%tc(i)/(1.d0-0.428571d0*log10(P/gl%pc(i))/(1.d0+gl%accen(i)))

    IF (P > 0.98d0*gl%pc(i)) THEN
    !  initial guesses for density;
    !  using correlation developed by E.W. Lemmon, NIST
       theta=(1.d0-P/gl%pc(i))*(1.6d0-gl%accen(i))
       Estimate_Tsat=gl%tc(i)*(1.0d0-0.103947d0*theta-4.108265d-2*theta**2)
    END IF
END IF

RETURN
END
!--------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION Res_Pdiff_Tsat(gl,Tsat_akt,parameters)
!function calculating the residual for iteration on Tsat_est
!the residual function is RES = P-Psat(Tsat-est)
!(J.P.Jakobsen, Denmark 09-2009)
!--------------------------------------------------------------------------
!USE module_fluid_parameters


implicit none

    type(type_gl) :: gl
  
DOUBLE PRECISION :: Tsat_akt,Press
    type(type_additional_parameters) :: parameters


integer:: i

!Starting value
Res_Pdiff_Tsat=1.d3

Press=parameters%a_p(1)
i=int(parameters%a_p(2))

Res_Pdiff_Tsat=Press-vp_eq(gl,Tsat_akt,i)
    
RETURN
END 



    end module flash_pure_module
