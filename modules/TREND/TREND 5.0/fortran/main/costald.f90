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

    ! module for file costald.f90
    module costald_module
    !global use inclusion
    use module_all_types
    use flash_pure_module
    use setup_module
    use cubic_eos_module

    contains



double precision function vtaiteq(gl,input, T, p)
!-------------------------------------------------------------------------------
!Taken from
!G. H. Thomson, K. R. Brobst, R. W. Hankinson, AIChE Journal, 28(4):671-676 (1982)
!calculates liquid molar volume or
!in case of input=TLIQ gives out saturated liquid molar volume
!-------------------------------------------------------------------------------





implicit none

    type(type_gl) :: gl


    Character*12 :: input
    Character*255:: fluids, path, EOS_indicator, moles
    double precision:: T, p
    double precision:: C_Tait, B_Tait, tau, psattait, vsatliq, etait, tredtait, pctait
    integer:: error, iter, errval
    double precision:: dvap, dliq
    double precision:: pr, pr0, pr1, alp, bet, zm

    B_Tait=0.d0
    dvap=0.d0
    dliq=0.d0
    error=0
    
    !calculate saturated liquid volume and directly returns it for input='tliq'
    vsatliq=vsattait(gl,T)
    
    !error catching: reducing parameter out of range of validity
    if (vsatliq==-5888.d0) then
        vtaiteq=vsatliq
        return
    end if

    if (trim(input) == 'tliq') then
        vtaiteq = vsatliq
        return !function ends here for saturated liquid volume
    end if
        
    if (gl%ncomp == 1) then
        !Critical pressure calculated using SRK
        pctait = gl%pc(1)
    else
        !Critical pressure for mixture taken from
        !G. H. Thomson, K. R. Brobst, R. W. Hankinson, AIChE Journal, 28(4):671-676 (1982)
        
        !Eq. (18)
        zm = 0.291d0 - 0.08d0 * gl%omegaCOS    
        
        !Eq. (17) critical pressure mixture
        pctait = ((zm * gl%Req(1) * gl%TcCOS) / gl%vstarCOS) * 0.001d0
    end if
    
    !reducing temperature using critical temperature of mixtures
    tredtait=T/gl%TcCOS
    tau = 1.d0 - tredtait
    
    if (gl%ncomp == 1) then
        !Vapor pressure calculation using SRK
        call Flash_Pure_PhaseBoundary(gl,psattait, t, dvap, dliq, 1, errval, iter, 1)
    else
        !Vapor pressure for mixture taken from
        !G. H. Thomson, K. R. Brobst, R. W. Hankinson, AIChE Journal, 28(4):671-676 (1982)
               
        !Eq. (22)
        alp = 35.d0 - 36.d0 / tredtait - 96.736d0 * dlog(tredtait) + tredtait**6
        
        !Eq. (23)
        bet = dlog(tredtait) + 0.03721754d0 * alp
        
        !Eq. (20)
        pr0 = 5.8031817d0 * dlog(tredtait) + 0.07608141d0 * alp
        
        !Eq. (21)
        pr1 = 4.86601d0 * bet
        
        !Eq. (19)
        pr = dexp(pr0 + gl%omegaCOS * pr1)
        
        !Eq. from paper text
        psattait = pr * pctait
    end if
               
    !Eq. (7)
    etait = dexp(ftait+gtait*gl%omegaCOS + htait*gl%omegaCOS**2)  
    
    !Eq. (8)
    C_tait = jtait + ktait * gl%omegaCOS

    !Eq. (6)
    B_Tait = pctait* (-1.d0 + atait*tau**(1.d0/3.d0) + btait*tau**(2.d0/3.d0) + dtait*tau + etait*tau**(4.d0/3.d0))   
        
    !Eq. (5), m³/mol
    vtaiteq = vsatliq * (1.d0 - C_Tait * dlog((B_Tait + p) / (B_Tait + psattait)))
        
    end function
    
    
double precision function vsattait(gl,T)
!-------------------------------------------------------------------------------
!Taken from
!R. W. Hankinson and G. H. Thomson, AIChE Journal, 25(4):653-663 (1979)
!calculates saturated liquid specific volume
!-------------------------------------------------------------------------------





implicit none

    type(type_gl) :: gl

    
    double precision:: T, tau, tredtait, vr0, vrdel
    integer:: nrsubst
    
    !sets the critical parameters for each investigated fluid
    do nrsubst=1,gl%ncomp
        call init_SRK(gl,T, nrsubst)
    end do
    
    !reducing temperature
    tredtait=T/gl%TcCOS
    tau = 1.d0 - tredtait
    
    !Eq. (17)
    if((tredtait>0.25d0) .and. (tredtait<0.95d0)) then
        vr0 = 1.d0 + ataitold*tau**(1.d0/3.d0) + btaitold*tau**(2.d0/3.d0) + ctaitold*tau + dtaitold*tau**(4.d0/3.d0)
    else 
        vsattait=-5888.d0
        return
    end if

    !Eq. (18)
    if((tredtait>0.25d0) .and. (tredtait<1.d0)) then
        vrdel = (etaitold + ftaitold * tredtait + gtaitold * tredtait**2 + htaitold * tredtait**3) / (tredtait - 1.00001d0)
    else
        vsattait=-5888.d0
        return
    end if

    !Eq. (16), m³/mol
    vsattait=gl%vstarCOS * vr0 * (1.d0 - gl%omegaCOS * vrdel) / gl%factor  
    
    end function
    
double precision function D_RKM(gl,T, p)
!-------------------------------------------------------------------------------
!Taken from
!Taken from: R. W. Hankinson and G. H. Thomson, AIChE Journal, 25(4):653-663 (1979)

!calculates liquid molar volume
!-------------------------------------------------------------------------------





implicit none

    type(type_gl) :: gl

    double precision::k1_used
    double precision::k1_T1
    double precision::k1_T2
    double precision::k2_used
    double precision::k2_T1
    double precision::k2_T2
    double precision::V_used(1:8), ps_used
    double precision::M(1:8)
    character(len=100)::filename
    double precision::T
    double precision::xN2, xCH4, xC2H6
    double precision::V_RKM
    double precision::p
    integer :: errorflag
    integer::i
    integer::j
    double precision::SUMxiVi
    double precision::wm_mix
    double precision::RKM_Mmix, RKM_Tc_mix, ERKM_add_term
    
    double precision::RKM_molfractions(1:8)
   
    k1_used = 0.d0
    k2_used = 0.d0

    !##############################################################################################
    !Temperature interpolation for molar volume
    !##############################################################################################
    
    call wm_mix_calc(gl,wm_mix)
    
    if (T > 135.D0) then
        errorflag = -5999
        D_RKM = errorflag
        return
    else if (T < 100.D0) then
        errorflag = -5999
        D_RKM = errorflag
        return
    else if (p > 10.D0) then
        errorflag = -5999
        D_RKM = errorflag
    end if
    
    Do i = 1,100
        if (gl%RKM_V(i,0) == T) then
            V_used(1:8) = gl%RKM_V(i,1:8)
            exit
        else if (gl%RKM_V(i,0) < T) then
            !nothing
        else if (gl%RKM_V(i,0) > T) then
            V_used(1:8) = gl%RKM_V(i-1,1:8) + ((gl%RKM_V(i,1:8) - gl%RKM_V(i-1,1:8)) / (gl%RKM_V(i,0) - gl%RKM_V(i-1,0)) * (T - gl%RKM_V(i-1,0)))
            exit
        end if
    end do

    !##############################################################################################
    !Temperature and molar mass interpolation for correction factor k1 & k2
    !##############################################################################################
   
    RKM_Mmix = wm_mix*gl%factor
    
    !k1
    Do i = 1, 10
        if (gl%k1(i,0) == T) then
                Do j = 1, 10
                    if (gl%k1(0,j) == RKM_Mmix) then
                        k1_used = gl%k1(i,j)
                        exit
                    else if (gl%k1(0,j) < RKM_Mmix) then
                        !nothing
                    else if (gl%k1(0,j) > RKM_Mmix) then
                        k1_used = gl%k1(i,j-1) + (gl%k1(i,j) - gl%k1(i,j-1)) / (gl%k1(0,j) - gl%k1(0,j-1)) * (RKM_Mmix - gl%k1(0,j-1))
                        exit
                    end if
                end do
            exit
        else if (gl%k1(i,0) < T) then
            !nothing
        else if (gl%k1(i,0) > T) then
                Do j = 1, 10
                    if (gl%k1(0,j) == RKM_Mmix) then
                        k1_used = gl%k1(i-1,j) + (gl%k1(i,j) - gl%k1(i-1,j)) / (gl%k1(i,0) - gl%k1(i-1,0)) * (T - gl%k1(i-1,0))
                        exit
                    else if (gl%k1(0,j) < RKM_Mmix) then
                        !nothing
                    else if (gl%k1(0,j) > RKM_Mmix) then
                        k1_T1 = gl%k1(i-1,j-1) + (gl%k1(i,j-1) - gl%k1(i-1,j-1)) / (gl%k1(i,0) - gl%k1(i-1,0)) * (T - gl%k1(i-1,0))
                        k1_T2 = gl%k1(i-1,j) + (gl%k1(i,j) - gl%k1(i-1,j)) / (gl%k1(i,0) - gl%k1(i-1,0)) * (T - gl%k1(i-1,0))
                        
                        k1_used = k1_T1 + (k1_T2 - k1_T1) / (gl%k1(0,j) - gl%k1(0,j-1)) * (RKM_Mmix - gl%k1(0,j-1))
                        exit
                    end if
                end do
            exit
        end if
    end do
    
    !k2
    Do i = 1, 10
        if (gl%k2(i,0) == T) then
            !k1_pre = k1(i,1:8)
                Do j = 1, 10
                    if (gl%k2(0,j) == RKM_Mmix) then
                        k2_used = gl%k2(i,j)
                        exit
                    else if (gl%k2(0,j) < RKM_Mmix) then
                        !nothing
                    else if (gl%k2(0,j) > RKM_Mmix) then
                        k2_used = gl%k2(i,j-1) + (gl%k2(i,j) - gl%k2(i,j-1)) / (gl%k2(0,j) - gl%k2(0,j-1)) * (RKM_Mmix - gl%k2(0,j-1))
                        exit
                    end if
                end do
            exit
        else if (gl%k2(i,0) < T) then
            !nothing
        else if (gl%k2(i,0) > T) then
                Do j = 1, 10
                    if (gl%k2(0,j) == RKM_Mmix) then
                        k2_used = gl%k2(i-1,j) + (gl%k2(i,j) - gl%k2(i-1,j)) / (gl%k2(i,0) - gl%k2(i-1,0)) * (T - gl%k2(i-1,0))
                        exit
                    else if (gl%k2(0,j) < RKM_Mmix) then
                        !nothing
                    else if (gl%k2(0,j) > RKM_Mmix) then
                        k2_T1 = gl%k2(i-1,j-1) + (gl%k2(i,j-1) - gl%k2(i-1,j-1)) / (gl%k2(i,0) - gl%k2(i-1,0)) * (T - gl%k2(i-1,0))
                        k2_T2 = gl%k2(i-1,j) + (gl%k2(i,j) - gl%k2(i-1,j)) / (gl%k2(i,0) - gl%k2(i-1,0)) * (T - gl%k2(i-1,0))
                        
                        k2_used = k2_T1 + (k2_T2 - k2_T1) / (gl%k2(0,j) - gl%k2(0,j-1)) * (RKM_Mmix - gl%k2(0,j-1))
                        exit
                    end if
                end do
            exit
        end if
    end do
    
    !##############################################################################################
    !Sort Mole Fractions
    !##############################################################################################
    
    RKM_molfractions=0.D0
    do i = 1, 8
        do j = 1, 8
            if ((gl%components(j) == gl%RKM_fluids(i,1)) .or. (gl%components(j) == gl%RKM_fluids(i,2)) .or. (gl%components(j) == gl%RKM_fluids(i,3))) then
                RKM_molfractions(i) = gl%molfractions(j)
            end if
        end do
    end do
    
    !##############################################################################################
    !The RKM equation itself
    !##############################################################################################
    SUMxiVi = 0.D0
    
	Do i = 1, 8
		SUMxiVi = (RKM_molfractions(i) * V_used(i)) + SUMxiVi
    end do
    
    xN2 = RKM_molfractions(6)
    xCH4 = RKM_molfractions(1)
    xC2H6 = RKM_molfractions(2)
    
	V_RKM = SUMxiVi - (k1_used + (k2_used - k1_used) * (xN2 / 0.0425)) * xCH4

    D_RKM = (1.D0/V_RKM)*gl%factor

    !##############################################################################################
    !The enhancement of the RKM equation in regard of pressure dependency (ERKM)
    !Tietz, Apr 2016
    !##############################################################################################
    if (gl%eq_type(1) == 81) then
        !Including a pressure dependency by use of an additional term (Enhanced Revised Klosek-McKinley method)
        
        !interpolation of the saturation pressure of pure methane
        Do i = 1,100
            if (gl%RKM_ps_me(i,0) == T) then
                ps_used = gl%RKM_ps_me(i,1)
                exit
            else if (gl%RKM_ps_me(i,0) < T) then
                !nothing
            else if (gl%RKM_ps_me(i,0) > T) then
                ps_used = gl%RKM_ps_me(i-1,1) + ((gl%RKM_ps_me(i,1) - gl%RKM_ps_me(i-1,1)) / (gl%RKM_ps_me(i,0) - gl%RKM_ps_me(i-1,0)) * (T - gl%RKM_ps_me(i-1,0)))
                exit
            end if
        end do
        
        RKM_Tc_mix = 0.D0
        
        !Calculation of a pseudo critical temperature of the mixture
	    Do i = 1, 8
		    RKM_Tc_mix = (RKM_molfractions(i) * gl%RKM_Tc(i)) + RKM_Tc_mix
        end do
    
        if (p == 0.d0) then
            p=ps_used
        end if
        
        if (p >= ps_used*0.9D0) then        
            ERKM_add_term = (p-ps_used-(xN2*0.11d0*(T-90.d0))+(xC2H6*0.05*(T-95.D0)))*0.000406d0*(RKM_Tc_mix/(RKM_Tc_mix-T))**1.77d0
            D_RKM = D_RKM*(1+ERKM_add_term)
        else
            D_RKM = -5999.D0
        end if
    end if
        
    
    end function
    


    end module costald_module
