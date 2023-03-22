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

    ! module for file ancillary_equations_mix.f90
    module ancillary_equations_mix_module
    !global use inclusion
    use module_all_types
    use module_regula_falsi
    use module_regula_falsi_support

    contains

	




Double Precision Function vp_eq(gl,T, nrsubst)
!  ----------------------------------------------------------
!  vapor pressure, EQ from E. Lemmon
!  T. Wiens, Denmark, 09-2009
!  ----------------------------------------------------------



    

implicit none

    type(type_gl) :: gl


    ! Declarations:
    !------------------------------------------------------------
    Double Precision T                          ! Input parameter

    Integer j, nrsubst, i                     ! int used in do-loop
    Double Precision Treduc                     ! reduced temperature
    Double Precision sum                        ! reduced pressure
    Double Precision Psat_anci                  ! vapor pressuer calculated only with fluid properties

    !-Integer vptype                            ! part of 'fluid' telling what kind of equation to use
    !-Integer nvpcoeff                          ! amount of coefficients
    !-Double Precision tvpred                     ! reducing temperature
    !-Double Precision vpcoeff                    ! coefficient
    !-Double Precision vpex                       ! exponent
    !-Double Precision vpred                      ! reducing pressure
    !------------------------------------------------------------
    i = nrsubst
    vp_eq = 0

     ! error catching-----------------------------------
    if (T <= 0) then
        vp_eq = 0.d0
        return
    end if

    if (T >= gl%tvpred(i) .and. gl%nvpcoeff(i) /= 0) then    
        vp_eq = gl%pc(i)
        return
        
    else if (gl%nvpcoeff(i) == 0) then      ! Temperature correlation for saturated pressure (from E.Lemmon)
        vp_eq=gl%pc(i)*10.d0**(-2.333333d0*(1.d0+gl%accen(i))*(gl%tc(i)/T-1.d0))
        !vp_eq = Psat_anci
        return
        
    end if
    !--------------------------------------------------- 
    
    if (gl%vptype(i) < 0) gl%vptype(i)=0                        ! error catching
    if (9 < gl%vptype(i)) gl%vptype(i)=0

    if (gl%vptype(i) > 0) Treduc = dabs(1.d0-T/gl%tvpred(i))
    if (mod(gl%vptype(i),2) == 0) Treduc = Treduc**0.5d0 ! even values

    sum = 0.d0

    if (gl%nvpcoeff(i) /= 0) then                                        
        do j=1, gl%nvpcoeff(i)
            sum = sum + gl%vpcoeff(j,i)*Treduc**gl%vpexp(j,i)
        end do
        if (gl%vptype(i) == 1 .or. gl%vptype(i) == 2) sum = 1.d0 + sum
        if (gl%vptype(i) == 3 .or. gl%vptype(i) == 4) sum = exp(sum)
        if (gl%vptype(i) == 5 .or. gl%vptype(i) == 6) sum = exp(gl%tvpred(i)/T*sum)

        vp_eq = gl%vpred(i)*sum
    end if

    ! error catching
    if (vp_eq < 0) vp_eq = 0.d0
    
End Function vp_eq



Double Precision Function dl_eq(gl,T, i)
!  ----------------------------------------------------------
!  saturated liquid density, EQ from E. Lemmon
!  T. Wiens, Denmark, 09-2009
!  ----------------------------------------------------------




implicit none

    type(type_gl) :: gl


    ! Declarations:
    !------------------------------------------------------------
    Double Precision T                          ! Input parameter
    Integer j, i                              ! int used in do-loop
    Double Precision Treduc                     ! reduced temperature
    Double Precision sum                        ! xxxx

    !-Integer dltype                            ! part of 'fluid' telling what kind of equation to use
    !-Integer ndlcoeff                          ! amount of coefficients
    !-Double Precision tdlred                     ! reducing temperature
    !-Double Precision dlcoeff                    ! coefficient
    !-Double Precision dlex                       ! exponent
    !-Double Precision dlred                      ! reducing density
    !------------------------------------------------------------
    
    dl_eq = 0.d0

     ! error catching-----------------------------------
    if (T <= 0) then
        dl_eq = 0.d0
        return
    end if

    
    if (T >= gl%tdlred(i) .and. gl%ndlcoeff(i) /= 0) then    
        dl_eq = gl%rhoc(i)
        return
    end if
    
    if (gl%ndlcoeff(i) == 0) then
    !  if no equation defined using correlation developed by R. Span
        !dl_eq = rhotp(i) - (rhotp(i)-rhoc(i))*(0.6634d0*(T/tc(i))**1.8d0 + 0.2478d0*(T/tc(i))**20)
        return
    end if
    !--------------------------------------------------- 
    
    if (gl%dltype(i) < 0) gl%dltype(i)=0                                                ! error catching
    if (9 < gl%dltype(i)) gl%dltype(i)=0

    if (gl%dltype(i) > 0) Treduc = dabs(1.d0-T/gl%tdlred(i))
    if (mod(gl%dltype(i),2) == 0) Treduc = Treduc**(1.d0/3.d0)                  ! even values

    sum = 0.d0

    if (gl%ndlcoeff(i) /= 0) then                                        
        do j=1, gl%ndlcoeff(i)
            sum = sum + gl%dlcoeff(j,i)*Treduc**gl%dlexp(j,i)
        end do
        if (gl%dltype(i) == 1 .or. gl%dltype(i) == 2) sum = 1.d0 + sum
        if (gl%dltype(i) == 3 .or. gl%dltype(i) == 4) sum = exp(sum)
        if (gl%dltype(i) == 5 .or. gl%dltype(i) == 6) sum = exp(gl%tdlred(i)/T*sum)

        dl_eq = gl%dlred(i)*sum
    end if
    
    ! error catching
    if (dl_eq < 0.d0) dl_eq = 0.d0
    
    
End Function dl_eq



Double Precision Function dv_eq(gl,T, i)
!  ----------------------------------------------------------
!  saturated vapor density, EQ from E. Lemmon
!  T. Wiens, Denmark, 09-2009
!  ----------------------------------------------------------




implicit none

    type(type_gl) :: gl


    ! Declarations:
    !------------------------------------------------------------
    Double Precision T                          ! Input parameter
    Integer j, i                              ! int used in do-loop
    Double Precision Treduc                     ! reduced temperature
    Double Precision sum                        ! xxx

    !-Integer dvtype                            ! part of 'fluid' telling what kind of equation to use
    !-Integer ndvcoeff                          ! amount of coefficients
    !-Double Precision tdvred                     ! reducing temperature
    !-Double Precision dvcoeff                    ! coefficient
    !-Double Precision dvex                       ! exponent
    !-Double Precision dvred                      ! reducing density
    !------------------------------------------------------------

    dv_eq = 0.d0

     ! error catching-----------------------------------
    if (T <= 0) then
        dv_eq = 0.d0
        return
    end if

    if (T >= gl%tdvred(i) .and. gl%ndvcoeff(i) /= 0) then    
        dv_eq = gl%rhoc(i)
        return
    end if
    
    if (gl%ndvcoeff(i) == 0) then
    !  if no equation defined using correlation developed by R. Span
    !  using correlation developed by R. Span
      !  Treduc = (T/tc(i))**(-4.4645d0)
      !  comp = 1.d0-1.691d0*exp(-Treduc)
      !  dv_eq = vp_eq(gl,T, i)/(Req(i)*T*comp)*factortrans!*wm(i)
      !if (dv_eq < 0.d0) dv_eq=1.d-6
      return
    end if
    !--------------------------------------------------- 
    
    
    if (gl%dvtype(i) < 0) gl%dvtype(i)=0                                                ! error catching
    if (9 < gl%dvtype(i)) gl%dvtype(i)=0

    if (gl%dvtype(i) > 0) Treduc = dabs(1.d0-T/gl%tdvred(i))
    if (mod(gl%dvtype(i),2) == 0) Treduc = Treduc**(1.d0/3.d0)                  ! even values

    sum = 0.d0

    if (gl%ndvcoeff(i) /= 0) then                                        
        do j=1, gl%ndvcoeff(i)
            sum = sum + gl%dvcoeff(j,i)*Treduc**gl%dvexp(j,i)
        end do
        if (gl%dvtype(i) == 1 .or. gl%dvtype(i) == 2) sum = 1.d0 + sum
        if (gl%dvtype(i) == 3 .or. gl%dvtype(i) == 4) sum = exp(sum)
        if (gl%dvtype(i) == 5 .or. gl%dvtype(i) == 6) sum = exp(gl%tdvred(i)/T*sum)

        dv_eq = gl%dvred(i)*sum
    end if
    
    ! error catching
    if (dv_eq < 0) dv_eq = 0.d0
    
End Function dv_eq




Double Precision Function pmelt_eq(gl,T, i)
!  ----------------------------------------------------------
!  melting pressure, EQ from E. Lemmon
!  T. Wiens, Denmark, 09-2009
!  ----------------------------------------------------------





implicit none

    type(type_gl) :: gl


    ! Declarations:
    !------------------------------------------------------------
    Double Precision T                                                   ! Input parameter
    Integer i, j, k                                                       ! int used in do-loop
    Double Precision Treduc, treduc4  ! reduced temperature
    Double Precision sum, sum2                                     ! reduced pressure
    Double Precision pi4, lnpi                            ! parameter
    Double Precision pmelt4, pmelt5       ! parameter
    Double Precision p1, p2, p3, p4     ! parameter     
    Double precision:: tnormd2o, pnormd2o, theta_d2o
    
    

    !-Integer*4 pmelttype                         ! part of 'fluid' telling what equation to use
    !-Integer*4 npmeltcoeff(1-3)                  ! 3 amounts of coefficients
    !-Double Precision tpmeltred                  ! reducing temperature
    !-Double Precision pmeltcoeff                 ! coefficient
    !-Double Precision pmeltexp                   ! exponent
    !-Double Precision pmeltred                   ! reducing pressure
    !------------------------------------------------------------

    ! error catching---------------------------------------------
    
    !Andreas Feb 2014
    pmelt_eq = 0.D0
    
    if (T <= 0) then
        pmelt_eq = 0.d0
        return
    end if
    
    if ((gl%pmelttype(i) /= 24) .and. (gl%pmelttype(i) /= 20)) then
    
        if ((gl%npmeltcoeff1(i)+gl%npmeltcoeff2(i)+gl%npmeltcoeff3(i)) == 0) then    !returning zero in case of no equation the iteration would be ruined, so pmaxfluid is returned
            pmelt_eq = gl%pmaxfluid(i)
            return
        end if
        !------------------------------------------------------------
    
        Treduc = T/gl%tpmeltred(i)
        sum = 0.d0

        do j=1, gl%npmeltcoeff1(i)
            sum = sum + gl%pmeltcoeff(j,i)*Treduc**gl%pmeltexp(j,i)
        end do

        do j=1, gl%npmeltcoeff2(i)
            k=j+gl%npmeltcoeff1(i)
            sum = sum + gl%pmeltcoeff(k,i)*(Treduc-1.d0)**gl%pmeltexp(k,i)
        end do

        do j=1, gl%npmeltcoeff3(i)
            k=j+gl%npmeltcoeff1(i)+gl%npmeltcoeff2(i)
            sum = sum + gl%pmeltcoeff(k,i)*log(Treduc)**gl%pmeltexp(k,i)
        end do

        if (gl%pmelttype(i) == 1)  then
            pmelt_eq = gl%pmeltred(i)*sum
        else if (gl%pmelttype(i) == 2)  then
            pmelt_eq = gl%pmeltred(i)*exp(sum)
    
        !---------------------------------------------------------------------------------------------------------------------------------
        !Andreas July 2011: For water ice another if statement is required, since the equationtype is not 1 or 2 but a wagner type equation
        else if (gl%pmelttype(i) == 39) then
            pmelt_eq = gl%pmeltred(i)*sum
        !---------------------------------------------------------------------------------------------------------------------------------
        else
            pmelt_eq = gl%pmeltred(i)*sum
        end if   
            
            
            
        !Melting Pressure equation, IAPWS 2011
     
    
    elseif (gl%pmelttype(i) == 24)then                                           !check if fluid is Water
    
    !Temperature Ranges of interleaving equations Ice Ih, III and V 
               
        !Andreas, Feb 2014
        !if ((256.164d0 >= t) .and. (t >= 251.165d0)) then                         !Ice Ih or Ice III model possible
        if (( (256.164d0 - t) >= -1.D-8) .and. ( (t - 251.165d0) > -1.D-8)) then     
        
            p1 = pmelt1_eq(gl,T,i)
            p2 = pmelt2_eq(gl,T,i)  
            
            if (gl%melt_it) then                                                 !Two pressures can be found due to phase-boundary curve
             
                if (gl%melt_p <= 208.566d0) then
             
                    pmelt_eq=p1
             
                else if ((gl%melt_p > 208.566d0) .and. (gl%melt_p <= 350.1d0)) then
             
                    pmelt_eq=p2
                end if
                
            else if (.not. gl%melt_it) then                                      !Set Warning due to Interleaving equations
            
             if (p1 < p2) then
                pmelt_eq=p1
                gl%pmelt_high=p2
             elseif (p1 > p2) then
                pmelt_eq=p2
                gl%pmelt_high=p1
             end if
            end if  
             
        else if ((273.16d0 >= t) .and. (t >= 256.164d0))  then                    !Ice Ih or Ice V model possible                  !
        
            p3 = pmelt1_eq(gl,T,i)
            p4 = pmelt3_eq(gl,T,i)
            
            if (gl%melt_it) then                                                 !Two pressures can be found due to phase-boundary curve  
            
                if (gl%melt_p <= 208.566d0) then
             
                    pmelt_eq=p3
             
                else if ((gl%melt_P > 350.1d0) .and. (gl%melt_p <= 632.4d0)) then
             
                    pmelt_eq=p4
                end if
                
            else if (.not. gl%melt_it) then                                     !Set Warning due to Interleaving equations
            
              if (p3 < p4) then
                pmelt_eq=p3
                gl%pmelt_high=p4
              elseif (p3 > p4) then
                pmelt_eq=p4
                gl%pmelt_high=p3
              end if 
             
            end if  
        
        else if ((273.31d0 >= t) .and. (t >= 273.16d0)) then                      !Only Ice V model possible                  
            
            pmelt_eq = pmelt3_eq(gl,t,i)                          
  
     
    !Ice VI equation 
     
        elseif  (( 355 >= t) .and.(t >= 273.31)) then
     
            treduc4 = t/273.31d0                                              !different reducing temperature in Ice VI
            pi4 = 0.d0 
     
            pi4= 1+gl%pmeltcoeff(6,1)*(1-treduc4**gl%pmeltexp(6,1))
    
            pmelt4=pi4*632.4d0                                                !melting pressure with use of special reducing pressure
            pmelt_eq = pmelt4
     
      !Ice VII equation  
     
        else if (( 715 >= t) .and. (t >= 355)) then 
    
            treduc = t/355.d0                                                  !different reducing temperature in Ice VI
            lnpi = 0.d0
            sum2 = 0.d0
    
            do j=gl%npmeltcoeff1(i)+gl%npmeltcoeff2(i)+gl%npmeltcoeff3(i)+ gl%npmeltcoeff4(i)+1,&
                gl%npmeltcoeff1(i)+gl%npmeltcoeff2(i)+gl%npmeltcoeff3(i)+ gl%npmeltcoeff4(i)+gl%npmeltcoeff5(i)
                sum2= sum2 + gl%pmeltcoeff(j,i)*(1-Treduc**gl%pmeltexp(j,i))
            end do
    
            lnpi= sum2
    
            pmelt5=dexp(lnpi)*2216.d0                                         !melting pressure with use of special reducing pressure   
            pmelt_eq = pmelt5
    
        end if

    elseif (gl%pmelttype(i) == 20)then                                           !check if fluid is heavy water
        !Ice Ih
        if (((t - 254.415d0) >= 1.D-8) .and. ((t - 276.969d0) < 1.D-8)) then
        !if ((t .ge. 254.415d0) .and. (t .lt. 276.969d0)) then
            tnormd2o = 276.969d0
            pnormd2o = 0.00066159d0
            theta_d2o = T/tnormd2o
            pmelt_eq = (1.d0 - 0.30153d5 * (1.d0 - theta_d2o ** 5.5d0) + 0.692503d6 * (1.d0 - theta_d2o ** 8.2d0)) * pnormd2o
        elseif ((t .ge. 254.415d0) .and. (t .lt. 258.661d0)) then
        !ice III
            tnormd2o = 254.415d0
            pnormd2o = 222.41d0
            theta_d2o = T/tnormd2o
            pmelt_eq = (1.d0 - 0.802871d0 * (1.d0 - theta_d2o ** 33.d0)) * pnormd2o
        elseif ((t .gt. 258.661d0) .and. (t .lt. 275.748d0)) then
        !ice V
            tnormd2o = 258.661d0
            pnormd2o = 352.19d0
            theta_d2o = T/tnormd2o
            pmelt_eq = (1.d0 - 0.1280388d1 * (1.d0 - theta_d2o ** 7.6d0)) * pnormd2o
        elseif ((t .ge. 275.748d0) .and. (t .le. 315.d0)) then
        !ice VI
            tnormd2o = 275.748d0
            pnormd2o = 634.53d0
            theta_d2o = T/tnormd2o
            pmelt_eq = (1.d0 - 0.1276026d1 * (1.d0 - theta_d2o ** 4.d0)) * pnormd2o
        end if
    end if    
    
    ! error catching
    if (pmelt_eq < 0.d0) pmelt_eq = 0.d0    
    
   
End Function pmelt_eq

Double Precision Function pmelt1_eq(gl,T, i) !Function for ice Ih equation 





implicit none

    type(type_gl) :: gl


    ! Declarations:
    !------------------------------------------------------------
    Double Precision :: T                                                   !Input parameter
    Double Precision :: sum1, treduc1, pi1, pmelt1               !Variables for calculation
    Integer :: i, j                                                       !used in do loop
    
    
    
if (gl%pmelttype(i) == 24)then

!Ice Ih equation
            Treduc1=0.d0            
            Treduc1 = T/gl%tpmeltred(i)
            pi1 = 0.d0
            sum1 = 0.d0
    
            do j=1, gl%npmeltcoeff1(i)  
                sum1= sum1 + gl%pmeltcoeff(j,i)*(1-(Treduc1**gl%pmeltexp(j,i)))
            end do
    
            pi1= 1+sum1

            
            !pmelt1=pi1*611.657d-6                                           !melting pressure 
            pmelt1=pi1*gl%pmeltred(i)/gl%factor                                        !melting pressure 
            
            pmelt1_eq = pmelt1
            
            
end if            


End Function pmelt1_eq

Double Precision Function pmelt2_eq(gl,T, i) !Function for ice III equation 





implicit none

    type(type_gl) :: gl


    ! Declarations:
    !------------------------------------------------------------
    Double Precision :: T                                                   !Input parameter
    Double Precision :: treduc2, pi2, pmelt2                          !Variables for calculation
    Integer :: i                                                       !used in do loop
    
    
    
if (gl%pmelttype(i) == 24)then

!Ice III equation
            
    treduc2 = 0.d0
            treduc2=t/251.165d0                                               !different reducing temperature in Ice III
            pi2 = 0.d0   
    
            pi2= 1+gl%pmeltcoeff(4,1)*(1-treduc2**gl%pmeltexp(4,1))
    
            pmelt2=pi2*208.566d0                                              !melting pressure with use of special reducing pressure
            
            pmelt2_eq=pmelt2           
end if            


End Function pmelt2_eq

Double Precision Function pmelt3_eq(gl,T, i) !Function for ice V equation 





implicit none

    type(type_gl) :: gl


    ! Declarations:
    !------------------------------------------------------------
    Double Precision :: T                                                   !Input parameter
    Double Precision :: treduc3, pi3, pmelt3                          !Variables for calculation
    Integer :: i                                                       !used in do loop
    
    
    
if (gl%pmelttype(i) == 24)then

!Ice V equation
            
    treduc3 = t/256.164d0                                                    !different reducing temperature in Ice V
            pi3 = 0.d0 
     
            pi3= 1+gl%pmeltcoeff(5,1)*(1-treduc3**gl%pmeltexp(5,1))
    
            pmelt3=pi3*350.1d0                                               !melting pressure with use of special reducing pressure 
            
            pmelt3_eq=pmelt3 
                    
end if            


End Function pmelt3_eq

Double Precision Function tmelt_eq(gl,pmelt, nrsubst)








implicit none

    type(type_gl) :: gl


INTEGER :: Max_Iterations, Iterations
DOUBLE PRECISION :: Tstart_min,Tstart_max				!variables for regula falsi
DOUBLE PRECISION :: Tmin_allowed,Tmax_allowed
DOUBLE PRECISION :: Delta_allowed,pmelt,tmelt
INTEGER :: errTmelt,nrsubst
type(type_additional_parameters) :: ResTmelt_param

gl%melt_it = .true.
tmelt = 0.d0
!ResTmelt_param = 0.d0
Tstart_min = 0.d0
Tstart_max = 0.d0
Tmin_allowed = 0.d0
Tmax_allowed = 0.d0
Delta_allowed = 0.d0
tmelt = 0.d0
errTmelt=0

!Parameters for iteration:
gl%melt_p = pmelt
Delta_allowed = 1.d-8
Max_Iterations = 50
ResTmelt_param%a_p(1) = pmelt
ResTmelt_param%a_p(2) = nrsubst
!ResPsat_param(65) = -999.D0 ! Indicate that Regula Falsi is called by the Maxwell iteration


if (gl%pmelttype(nrsubst) == 24) then

!Parameter declaration for Water: 

    Tstart_min = 251.165d0
    Tmin_allowed = 251.165d0   
    
    if (pmelt > 208.566d0) then         
        Tstart_max = 715.d0
        Tmax_allowed = 750.d0
        
     else  
        Tstart_max = 273.16d0          !ice Ih only defined in this temperature region
        Tmax_allowed = 273.16d0   
        
     end if         

elseif (gl%pmelttype(nrsubst) == 20) then
!Parameter declaration for heavy water: 
  
    Tstart_min = 254.415d0
    Tmin_allowed = 254.415d0   
    
    if (pmelt > 222.41d0) then         
        Tstart_max = 315.d0
        Tmax_allowed = 350.d0
        
     else  
        Tstart_max = 275.748d0          !ice Ih only defined in this temperature region
        Tmax_allowed = 275.748d0   
        !Tstart_max = 275.74d0          !ice Ih only defined in this temperature region
        !Tmax_allowed = 275.74d0   
        
     end if         
     
else

!Parameter declaration for all fluids:

Tstart_min = gl%ttp(nrsubst)
Tstart_max = gl%tred(nrsubst)

Tmin_allowed = gl%ttp(nrsubst)*0.9d0
Tmax_allowed = gl%tred(nrsubst)*1.1d0

end if

CALL Regula_Falsi(gl,Res_Tmelt,tmelt,Tstart_min,Tstart_max,Delta_allowed,&
Tmin_allowed,Tmax_allowed, Max_iterations,Iterations, errTmelt,ResTmelt_param)

tmelt_eq=Tmelt

gl%melt_it = .false.


End Function tmelt_eq

!--------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION Res_Tmelt(gl,T_akt,parameters)
!function calculating the residual for iteration on Psat
!the residual function is derived from the phase equilibria condition giL=giV
!(J.P.Jakobsen, Denmark 09-2009)
!--------------------------------------------------------------------------
!USE module_fluid_parameters


implicit none

    type(type_gl) :: gl
  

DOUBLE PRECISION :: T_akt, pmelt, pmelt_akt
type(type_additional_parameters) :: parameters
Integer :: i


! Errorflag = 0 = >  Iteration succeeded
! Errorflag = 1 = >  Number of allowed iterations exceeded
! Errorflag = 2 = >  No root found without conflicts with Xmin_allowed or Xmax_allowed
! Errorflag = 3 = >  The initial interval for X had to be enlarged or reduced to find a root

pmelt=parameters%a_p(1)
i=int(parameters%a_p(2))  

pmelt_akt = pmelt_eq(gl,t_akt,i) 

Res_Tmelt=pmelt-pmelt_akt

RETURN
END 

Double Precision Function psub_eq(gl,T, i)
!  ----------------------------------------------------------
!  sublimation pressure, EQ from E. Lemmon
!  T. Wiens, Denmark, 09-2009
!  ----------------------------------------------------------




implicit none

    type(type_gl) :: gl


    ! Declarations:
    !------------------------------------------------------------
    Double Precision T                          ! Input parameter
    Integer k, j, i                              ! int used in do-loop
    Double Precision Treduc                     ! reduced temperature
    Double Precision sum                        ! reduced pressure

    !-Integer*4 psubtype                         ! part of 'fluid' telling what equation to use
    !-Integer*4 npsubcoeff(1-3)                  ! 3 amounts of coefficients
    !-Double Precision tpsubred                  ! reducing temperature
    !-Double Precision psubcoeff                 ! coefficient
    !-Double Precision psubexp                   ! exponent
    !-Double Precision psubred                   ! reducing pressure
    !------------------------------------------------------------

    ! error catching---------------------------------------------
    if (T <= 0) then
        psub_eq = 0.d0
        return
    end if

   if (T >= gl%tpsubred(i) .and. (gl%npsubcoeff1(i)+gl%npsubcoeff2(i)+gl%npsubcoeff3(i)) /= 0) then
        psub_eq = gl%ptp(i)
        return
   end if
   
   if ((gl%npsubcoeff1(i)+gl%npsubcoeff2(i)+gl%npsubcoeff3(i)) ==  0) then !returning zero in case of no equation the iteration would be ruined, so a fictional vapor pressure is returned
        psub_eq = vp_eq(gl,T, i)
        gl%p_fake = .true.  !if a pure fluid is calculated, the sublimation pressure extrapolated from the vapor pressure curve will be returned; this boolean variable enables to check if it is an extrapolated or a real pressure
        return
   end if
!-----------------------------------------------------------------


    Treduc = T/gl%tpsubred(i)
    sum = 0.d0

    do j=1, gl%npsubcoeff1(i)
        sum = sum + gl%psubcoeff(j,i)*Treduc**gl%psubexp(j,i)
    end do

    do j=1, gl%npsubcoeff2(i)
        k=j+gl%npsubcoeff1(i)
        sum = sum + gl%psubcoeff(k,i)*(1.d0-Treduc)**gl%psubexp(k,i)
    end do

    do j=1, gl%npsubcoeff3(i)
        k=j+gl%npsubcoeff1(i)+gl%npsubcoeff2(i)
        sum = sum + gl%psubcoeff(k,i)*log(Treduc)**gl%psubexp(k,i)
    end do

    if (gl%psubtype(i) == 1)  psub_eq = gl%psubred(i)*sum
    if(sum < 200) then
        if (gl%psubtype(i) == 2) psub_eq = gl%psubred(i)*exp(sum)
        if (gl%psubtype(i) == 3 .and. T /= 0.d0)  psub_eq = gl%psubred(i)*exp(gl%tpsubred(i)/T*sum) ! used for CO2
    else
        sum = 0.d0
    end if
    
    ! error catching
    if (psub_eq < 0.d0) psub_eq = 0.d0
    
End Function psub_eq


!Function for the iterative solution of the sublimation temperature at given pressure from the ancillary equation 
!Andreas Feb 2014
Double Precision Function tsub_eq(gl,psub, nrsubst, errval)
!--------------------------------------------------------------------------------------------------------
!INPUT: 
!   psub    -   sublimation pressure in MPa
!   nrsubst -   specifies the component of which the sublimation pressure is to be computed
!
!OUTPUT:
!   tsub_eq -   sublimation temperature in MPa








implicit none

    type(type_gl) :: gl


INTEGER :: Max_Iterations, Iterations
DOUBLE PRECISION :: Tstart_min,Tstart_max				!variables for regula falsi
DOUBLE PRECISION :: Tmin_allowed,Tmax_allowed
DOUBLE PRECISION :: Delta_allowed,psub,tsub
INTEGER:: errval,nrsubst
type(type_additional_parameters) :: parameters

tsub = 0.d0
!parameters = 0.d0
Tstart_min = 0.d0
Tstart_max = 0.d0
Tmin_allowed = 0.d0
Tmax_allowed = 0.d0
Delta_allowed = 0.d0
tsub_eq = 0.D0
errval=0

!Parameters for iteration:
Delta_allowed = 1.d-8
Max_Iterations = 50
parameters%a_p(1) = psub
parameters%a_p(2) = nrsubst

!Parameter declaration for all fluids:

Tstart_min = 100.D0         !Arbitrary temperature
Tstart_max = gl%ttp(nrsubst)

Tmin_allowed = 50.D0        !Arbitrary temperature
Tmax_allowed = gl%ttp(nrsubst)

!errorhandling
if (psub < 1.d-12) then
    errval = -9933
    return
end if

CALL Regula_Falsi(gl,Res_Tsub,tsub,Tstart_min,Tstart_max,Delta_allowed,&
& Tmin_allowed,Tmax_allowed, Max_iterations,Iterations, errval, parameters)

if (tsub < 1.D-12) then
    errval = -2223
end if

tsub_eq=tsub


End Function tsub_eq

!--------------------------------------------------------------------------
DOUBLE PRECISION FUNCTION Res_Tsub(gl,T,parameters)
! Function for the iterative solution of the sublimation temperature
! INPUT: 
!   T           -   Temperature in K
!   parameters  -   Vector with values passed by the regula falsi routine
!--------------------------------------------------------------------------
!USE module_fluid_parameters


implicit none

    type(type_gl) :: gl
  

DOUBLE PRECISION :: T, psub, psub_calc
type(type_additional_parameters) :: parameters
Integer :: nrsubst

psub=parameters%a_p(1)
nrsubst= int(parameters%a_p(2))  

psub_calc = psub_eq(gl,T,nrsubst)

!Do iteration with logarithmic values, unless the calculation of psub fails (negative values or 0)
if (psub_calc > 1.d-12) then
    Res_Tsub=dlog(psub)-dlog(psub_calc) 
else
    Res_Tsub=psub-psub_calc
end if

RETURN
    END 
    
    
    Double Precision Function vp_oil(gl,T)
!  ----------------------------------------------------------
!  vapor pressure, Wagner equation for oil
!  T. Eckermann, Bochum 2014
!  ----------------------------------------------------------



    

implicit none

    type(type_gl) :: gl


    ! Declarations:
    !------------------------------------------------------------
    Double Precision :: T                   ! Input parameter
    Double Precision :: tau                         ! reduced temperature
    Double Precision :: cn(4)

    vp_oil = 0.d0
    
    tau = 1.d0-t/gl%tc(1)
    
        !cn(1) = 4.98621961531955d0
        !cn(2) = -5.1314531433172d0
        !cn(3) = -11.997017515573d0
        !cn(4) = 9.00761772253395d0
        cn(1) = -2.2d0
        cn(2) = -1.92684297681832d0
        cn(3) = -1.23954967794587d0
        cn(4) = -0.311415969687328d0    
    
    vp_oil = gl%pc(1)*dexp(gl%tc(1)/T * ( cn(1)*tau + cn(2)*tau**1.5d0 + cn(3)*tau**3 + cn(4)*tau**6 ))



    ! error catching
    if (vp_oil < 0) vp_oil = 0.d0
    if (T > gl%tc(1)) vp_oil = gl%pc(1) * T/gl%tc(1)
    
End Function vp_oil


    end module ancillary_equations_mix_module
