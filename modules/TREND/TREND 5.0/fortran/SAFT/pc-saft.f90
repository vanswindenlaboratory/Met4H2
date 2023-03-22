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
    !Cite as: Span, R.; Beckm체ller, R.; Eckermann, T.; Herrig, S.; Hielscher, S.; 
	!          J채ger, A.; Mickoleit, E.; Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (2020): 	
    !          TREND. Thermodynamic Reference and Engineering Data 5.0. 
    !          Lehrstuhl f체r Thermodynamik, Ruhr-Universit채t Bochum.

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

    ! module for file pc-saft.f90
    module pc_saft_module
    !global use inclusion
    use module_all_types
    use pc_saft_ancillary_routines_module
    use variables_transformation_module

    contains


subroutine ARDERIVS(gl,T, D_in, getvector, DERAR, nrsubst)

! Henning Markgraf, June 2016
    
    ! a_res: residual part of the Helmholtz free energy
    ! defined by eq. A.3 in Gross, Sadowski 2001:
    ! a_res = a_hc + a_disp
    ! dependent on D and T







implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    ! input
    double precision, intent (inout) :: T, D_in
    integer, dimension (nderivs), intent (in) :: getvector
    integer, intent (inout) :: nrsubst
    ! output
    double precision, dimension (nderivs), intent (out) :: DERAR
    ! working variables
    double precision :: D
    integer, dimension (nderivs) :: get
    integer, dimension (nderivs) :: getprevious
    integer :: i
    logical :: already_calculated

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    ! II. Initializations
    call init_pure_mix(gl,nrsubst)
    call allocate_arrays_PCSAFT(gl,nrsubst) ! TE nrsubst
    
    ! III. Change the units of the density from mol/m^3 to 1/Angstrom
    D = D_in*N_A
    D = D_in*N_A/1.d30
    
    ! IV. Check if any derivatives of this T, rho, molfractions combination is already saved in the module variables
    get = getvector ! get can be changed if some derivatives have already been calculated
    
    ! Adding the necessary additional T derivatives required for the T to tau conversion:
    call add_T_conversion_derivs(gl,get)

    call check_calculated_PCSAFT(gl,T,D,nrsubst,gl%ar_calculated_PCSAFT,get,already_calculated) ! Warning: already_calculated is true, if T,D and molfractions are the same
                                                                                    ! but gives no information if the same derivatives have been calculated
                                                                                    ! Have to look at the output of get, too!
    if (any(get .EQ. 1)) then ! no need to calculate anything if get = 0
        
        ! V. some of the functions require all lower order derivatives, too
        call add_previous_derivs(gl,get,getprevious)
        call check_calculated_PCSAFT(gl,T,D,nrsubst,gl%ar_calculated_PCSAFT,getprevious,already_calculated)
        
        ! VI. Initialize module variables if T, rho or molfractions different
        if (already_calculated .EQV. .FALSE.) then
            call init_derivs_PCSAFT(gl,nrsubst) ! all saved derivatives are deleted
        end if
        
        ! VII. update the module variables
        call update_modulevariables_PCSAFT(gl,T,D,nrsubst,gl%ar_calculated_PCSAFT,get)
    
        ! VIII.   All of the function parts, which A_r = A_disp + A_hc consists of, are called
        ! i.e.: d_i, zeta_n, mmean, ab, g_ii, I_1, I_2, meo1, meo2, C_1
        call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
        ! and a_hs
        call AHSDERIVS(gl,get)
        
        ! IX. calculate a_hc and a_disp, which depend on the previously calculated functions
        call AHCDERIVS(gl,get)
        call ADISPDERIVS(gl,D,get)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !optional supplements
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !dipol
        if (any(gl%MyPCSAFTD(gl%n_start:gl%n_end) > 1.d-12)) then
            call calculate_PCSAFT_functionparts_DD_Trho(gl,T,D,getprevious)
            call ADDDERIVS(gl,T,D,get) 
        end if
        
        !quadrupol        
        if (any(gl%QPCSAFTQ(gl%n_start:gl%n_end) > 1.d-12))  then
            call calculate_PCSAFT_functionparts_QQ_Trho(gl,T,D,getprevious)
            call AQQDERIVS(gl,T,D,get)
        end if
        
        !association
        if (any(gl%kabPCSAFT(gl%n_start:gl%n_end) > 1.d-12)) then 
            call AASSOC(gl,T,D,get)
        end if
        
        ! X. calculate a_res = a_hc + a_disp
        do i = 1, nderivs
            if (get(i) .eq. 1) then
                gl%ar_PCSAFT(i) = gl%ahc_PCSAFT(i) + gl%adisp_PCSAFT(i) + gl%AQQ_PCSAFTQ(i) + gl%ADD_PCSAFTD(i) + gl%AASSOC_PCSAFT(i) 
            end if
        end do
            
    end if
        
    if (nrsubst /= 0) then 
        gl%molfractions(nrsubst) = gl%molfractions_save
        gl%mol_save = .false.
        gl%n_start = 1
        gl%n_end = gl%ncomp
    end if

    ! XI. Change the T derivatives to tau derivatives
    call convert_T_derivs_Trho(gl,getvector,gl%ar_PCSAFT,DERAR)
    

!DEC$ END IF
end subroutine ARDERIVS

subroutine ADISPDERIVS(gl,D,GETDERADISP)

! Henning Markgraf, June 2016

    ! a_disp: dispersion contribution to the Helmholtz free energy
    ! defined by eq. A.10 in Gross, Sadowski 2001:
    ! a_disp = -2*pi*rho*I_1*meo1-pi*rho*mmean*C_1*I_2*meo2
    ! dependent on D and T





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    double precision, intent (in) :: D
    integer, dimension (nderivs), intent (in) :: GETDERADISP
    !output: adisp_PCSAFT (module variable)
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. calculate the derivatives of a_disp
    ! 1: a_disp
    if (GETDERADISP(1) .eq. 1) then
        !calculate a_disp
        gl%adisp_PCSAFT(1) = -2.d0*piPCSAFT*D*gl%i1_PCSAFT(1)*gl%meo1_PCSAFT(1)-piPCSAFT*D*gl%mmean_PCSAFT*gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)
    end if

    ! 2: 1ST DERIVATIVE OF a_disp WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires C_1, i1_PCSAFT, i2_PCSAFT, meo1_PCSAFT, meo2_PCSAFT, c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2)
    if (GETDERADISP(2) .eq. 1) then
        gl%adisp_PCSAFT(2) = -piPCSAFT*gl%mmean_PCSAFT*D*gl%c_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%i2_PCSAFT(2) &
            & - piPCSAFT*gl%mmean_PCSAFT*D*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%c_PCSAFT(2) &
            & - D*piPCSAFT*gl%mmean_PCSAFT*gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1) - 2.d0*piPCSAFT*D*gl%meo1_PCSAFT(1)*gl%i1_PCSAFT(2) &
            & - D*2.d0*piPCSAFT*gl%i1_PCSAFT(1)*gl%meo1_PCSAFT(1)
    end if

    ! 3: 2ND DERIVATIVE OF a_disp WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires C_1, i2_PCSAFT, meo1_PCSAFT, meo2_PCSAFT, c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2), c_PCSAFT(3), i1_PCSAFT(3), i2_PCSAFT(3)
    if (GETDERADISP(3) .eq. 1) then
        gl%adisp_PCSAFT(3) = -piPCSAFT*(gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D + 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D + &
            & gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D + 2.d0*gl%i1_PCSAFT(3)*gl%meo1_PCSAFT(1)*D) &
            & -piPCSAFT*D*(2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1) + 4.d0*gl%i1_PCSAFT(2)*gl%meo1_PCSAFT(1) + 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1))
    end if

    ! 4: 1ST DERIVATIVE OF a_disp WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires C_1, i1_PCSAFT, i2_PCSAFT, meo1_PCSAFT, meo2_PCSAFT, meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4)
    if (GETDERADISP(4) .eq. 1) then
        gl%adisp_PCSAFT(4) = -piPCSAFT*gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(4)*D - piPCSAFT*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D - &
            & piPCSAFT*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D - 2.d0*piPCSAFT*gl%i1_PCSAFT(1)*gl%meo1_PCSAFT(4)*D - 2.d0*piPCSAFT*gl%i1_PCSAFT(4)* &
            & gl%meo1_PCSAFT(1)*D
    end if

    ! 5: 2ND DERIVATIVE OF a_disp WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires C_1, i1_PCSAFT, i2_PCSAFT, meo1_PCSAFT, meo2_PCSAFT, meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4), meo1_PCSAFT(5), meo2_PCSAFT(5), c_PCSAFT(5), i1_PCSAFT(5), i2_PCSAFT(5)
    if (GETDERADISP(5) .eq. 1) then
        gl%adisp_PCSAFT(5) = -piPCSAFT*D*(gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(5) + 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(4) &
            & + gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1) + 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(4) + 2.d0*gl%c_PCSAFT(4)* &
            & gl%i2_PCSAFT(4)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1) + gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1) + 2.d0*gl%i1_PCSAFT(1)*gl%meo1_PCSAFT(5) + 4.d0* &
            & gl%i1_PCSAFT(4)*gl%meo1_PCSAFT(4) + 2.d0*gl%i1_PCSAFT(5)*gl%meo1_PCSAFT(1))
    end if

    ! 6: 2ND MIXED DERIVATIVE OF a_disp WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires C_1, i1_PCSAFT, i2_PCSAFT, meo1_PCSAFT, meo2_PCSAFT, c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2), meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4), c_PCSAFT(6), i1_PCSAFT(6), i2_PCSAFT(6)
    if (GETDERADISP(6) .eq. 1) then
        gl%adisp_PCSAFT(6) = D*(-piPCSAFT*gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(4) - piPCSAFT*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1) &
            & - piPCSAFT*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1) - 2.d0*piPCSAFT*gl%i1_PCSAFT(1)*gl%meo1_PCSAFT(4) - 2.d0*piPCSAFT* &
            & gl%i1_PCSAFT(4)*gl%meo1_PCSAFT(1)) - piPCSAFT*gl%c_PCSAFT(1)*gl%i2_PCSAFT(6)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D - piPCSAFT*gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(4)*D - &
            & piPCSAFT*gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D - piPCSAFT*gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(4)*D - &
            & piPCSAFT*gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D - piPCSAFT*gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D - &
            & 2.d0*piPCSAFT*gl%i1_PCSAFT(6)*gl%meo1_PCSAFT(1)*D - 2.d0*piPCSAFT*gl%i1_PCSAFT(2)*gl%meo1_PCSAFT(4)*D
    end if

    ! 7: 3RD MIXED DERIVATIVE OF a_disp WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires C_1, i1_PCSAFT, i2_PCSAFT, meo1_PCSAFT, meo2_PCSAFT, c_PCSAFT(2), i1_PCSAFT(2), i2_PCSAFT(2), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4), meo1_PCSAFT(4), meo2_PCSAFT(4),
    !          c_PCSAFT(6), i1_PCSAFT(6), i2_PCSAFT(6), c_PCSAFT(5), i1_PCSAFT(5), i2_PCSAFT(5), meo1_PCSAFT(5), meo2_PCSAFT(5),
    !          i1_PCSAFT(7), i2_PCSAFT(7), c_PCSAFT(7)
    ! c_PCSAFT(1,2,4,5,6), I(1,2,4,5,6,7), meo(1,4,5)
    if (GETDERADISP(7) .eq. 1) then
        gl%adisp_PCSAFT(7) = -piPCSAFT*D*(gl%c_PCSAFT(1)*gl%i2_PCSAFT(7)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(6)* &
            & gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(6)* &
            & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(5)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)* &
            & gl%mmean_PCSAFT + gl%c_PCSAFT(7)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT &
            & + 2.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 2.d0* &
            & gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(2)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0* &
            & gl%i1_PCSAFT(6)*gl%meo1_PCSAFT(4) + 2.d0*gl%i1_PCSAFT(7)*gl%meo1_PCSAFT(1) + 2.d0*gl%i1_PCSAFT(2)*gl%meo1_PCSAFT(5) + &
            & gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)* &
            & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT &
            & + gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%i1_PCSAFT(1)*gl%meo1_PCSAFT(5) + 4.d0*gl%i1_PCSAFT(4)*gl%meo1_PCSAFT(4) + 2.d0* &
            & gl%i1_PCSAFT(5)*gl%meo1_PCSAFT(1))
    end if

    ! 8: 3RD DERIVATIVE OF a_disp WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires C_1, i2_PCSAFT, meo1_PCSAFT, meo2_PCSAFT, c_PCSAFT(2), i2_PCSAFT(2), c_PCSAFT(3), i1_PCSAFT(3), i2_PCSAFT(3), i1_PCSAFT(8), i2_PCSAFT(8), c_PCSAFT(8)
    if (GETDERADISP(8) .eq. 1) then
        gl%adisp_PCSAFT(8) = -piPCSAFT*(gl%c_PCSAFT(1)*gl%i2_PCSAFT(8)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D + &
            & 3.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(3)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D + 3.d0*gl%c_PCSAFT(3)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D + &
            & gl%c_PCSAFT(8)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1)*D + 2.d0*gl%i1_PCSAFT(8)*gl%meo1_PCSAFT(1)*D) &
            & -piPCSAFT*D*(3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1) + 6.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1) + &
            & 3.d0*gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%mmean_PCSAFT*gl%meo2_PCSAFT(1) + 6.d0*gl%i1_PCSAFT(3)*gl%meo1_PCSAFT(1))
    end if

    ! 9: 3RD DERIVATIVE OF a_disp WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires C_1, i1_PCSAFT, i2_PCSAFT, meo1_PCSAFT, meo2_PCSAFT, meo1_PCSAFT(4), meo2_PCSAFT(4), c_PCSAFT(4), i1_PCSAFT(4), i2_PCSAFT(4), meo1_PCSAFT(5), meo2_PCSAFT(5), c_PCSAFT(5), i1_PCSAFT(5), i2_PCSAFT(5)
    !          c_PCSAFT(9), i1_PCSAFT(9), i2_PCSAFT(9), meo1_PCSAFT(9), meo2_PCSAFT(9)
    ! c_PCSAFT(1,4,5,9), I(1,4,5,9), meo(1,4,5,9)
    if (GETDERADISP(9) .eq. 1) then
        gl%adisp_PCSAFT(9) = -piPCSAFT*D*(gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(9)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(5)* &
            & gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(9)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0 &
            & *gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 6.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(4)* &
            & gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(5)*gl%i2_PCSAFT(4) &
            & *gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(9)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%i1_PCSAFT(1)*gl%meo1_PCSAFT(9) + 6.d0* &
            & gl%i1_PCSAFT(4)*gl%meo1_PCSAFT(5) + 6.d0*gl%i1_PCSAFT(5)*gl%meo1_PCSAFT(4) + 2.d0*gl%i1_PCSAFT(9)*gl%meo1_PCSAFT(1))
    end if

    ! 10: 3RD MIXED DERIVATIVE OF a_disp WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! c_PCSAFT(1,2,3,4,6,10), I(1,2,3,4,6,10), meo(1,4)
    if (GETDERADISP(10) .eq. 1) then
        gl%adisp_PCSAFT(10) = -piPCSAFT*D*(2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(10)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + &
            & 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 2.d0* &
            & gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(4)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0* &
            & gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + &
            & gl%c_PCSAFT(10)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 2.d0* &
            & gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0* &
            & gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + &
            & gl%c_PCSAFT(3)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0*gl%i1_PCSAFT(6)*gl%meo1_PCSAFT(1) + 2.d0*gl%i1_PCSAFT(10)* &
            & gl%meo1_PCSAFT(1) + 4.d0*gl%i1_PCSAFT(2)*gl%meo1_PCSAFT(4) + 2.d0*gl%i1_PCSAFT(3)*gl%meo1_PCSAFT(4))
    end if

    ! 11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^4
    ! c_PCSAFT(1,2,3,8,11), I(1,2,3,8,11), meo(1)
    if (GETDERADISP(11) .eq. 1) then
        gl%adisp_PCSAFT(11) = -piPCSAFT*D*(4.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(8)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(11)* &
            & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 12.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0*gl%c_PCSAFT(2)* &
            & gl%i2_PCSAFT(8)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 12.d0*gl%c_PCSAFT(3)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 6.d0* &
            & gl%c_PCSAFT(3)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0*gl%c_PCSAFT(8)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT &
            & + 4.d0*gl%c_PCSAFT(8)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(11)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)* &
            & gl%mmean_PCSAFT + 8.d0*gl%i1_PCSAFT(8)*gl%meo1_PCSAFT(1) + 2.d0*gl%i1_PCSAFT(11)*gl%meo1_PCSAFT(1))
    end if

    ! 12: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, DEL AND DEL, MULTIPLIED BY TAU*DEL^3
    ! c_PCSAFT(1,2,3,4,6,8,10,12), I(1,2,3,4,6,8,10,12), meo(1,4)
    if (GETDERADISP(12) .eq. 1) then
        gl%adisp_PCSAFT(12) = -piPCSAFT*D*(3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(10)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(12)*gl%meo2_PCSAFT(1)* &
            & gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(8)*gl%meo2_PCSAFT(4)* &
            & gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(4)*gl%i2_PCSAFT(8)*gl%meo2_PCSAFT(1)* &
            & gl%mmean_PCSAFT + 6.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(3)* &
            & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(10)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(10)*gl%i2_PCSAFT(2) &
            & *gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(12)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 6.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(6)* &
            & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(10)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 6.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(2) &
            & *gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(3)* &
            & gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(3)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(3)* &
            & gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(3)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + &
            & gl%c_PCSAFT(8)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(8)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + &
            & 6.d0*gl%i1_PCSAFT(10)*gl%meo1_PCSAFT(1) + 2.d0*gl%i1_PCSAFT(12)*gl%meo1_PCSAFT(1) + 6.d0*gl%i1_PCSAFT(3)*gl%meo1_PCSAFT(4) &
            & + 2.d0*gl%i1_PCSAFT(8)*gl%meo1_PCSAFT(4))
    end if

    ! 13: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, DEL, AND DEL, MULTIPLIED BY TAU^2*DEL^2
    ! c_PCSAFT(1,2,3,4,5,6,7,10,13), I(1,2,3,4,5,6,7,10,13), meo(1,4,5)
    if (GETDERADISP(13) .eq. 1) then
        gl%adisp_PCSAFT(13) = -piPCSAFT*D*(2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(7)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(13)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT &
            & + 4.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(10)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + &
            & 2.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 4.d0* &
            & gl%c_PCSAFT(4)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(10)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0* &
            & gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 2.d0* &
            & gl%c_PCSAFT(5)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(5)*gl%i2_PCSAFT(3)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0* &
            & gl%c_PCSAFT(7)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(7)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + &
            & gl%c_PCSAFT(13)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 4.d0* &
            & gl%c_PCSAFT(6)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0* &
            & gl%c_PCSAFT(6)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(10)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + &
            & 2.d0*gl%c_PCSAFT(10)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 4.d0* &
            & gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0* &
            & gl%c_PCSAFT(2)*gl%i2_PCSAFT(7)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + &
            & 2.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + gl%c_PCSAFT(3)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + &
            & 2.d0*gl%c_PCSAFT(3)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(3)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + &
            & 4.d0*gl%i1_PCSAFT(7)*gl%meo1_PCSAFT(1) + 2.d0*gl%i1_PCSAFT(13)*gl%meo1_PCSAFT(1) + 8.d0*gl%i1_PCSAFT(6)*gl%meo1_PCSAFT(4) + 4.d0* &
            & gl%i1_PCSAFT(10)*gl%meo1_PCSAFT(4) + 4.d0*gl%i1_PCSAFT(2)*gl%meo1_PCSAFT(5) + 2.d0*gl%i1_PCSAFT(3)*gl%meo1_PCSAFT(5))
    end if

    ! 14: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, TAU AND DEL, MULTIPLIED BY TAU^3*DEL
    ! c_PCSAFT(1,2,4,5,6,7,9,14), I(1,2,4,5,6,7,9,14), meo(1,4,5,9)
    if (GETDERADISP(14) .eq. 1) then
        gl%adisp_PCSAFT(14) = -piPCSAFT*D*(gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(9)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 3.d0* &
            & gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(9)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(14) &
            & *gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(7)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(6)* &
            & gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(9)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)* &
            & gl%mmean_PCSAFT + 6.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + &
            & 3.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(7)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 6.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 3.d0* &
            & gl%c_PCSAFT(4)*gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 3.d0* &
            & gl%c_PCSAFT(5)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(5)*gl%i2_PCSAFT(6)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(5) &
            & *gl%i2_PCSAFT(2)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(9)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(9)*gl%i2_PCSAFT(2)* &
            & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(14)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(7)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)* &
            & gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(7)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)* &
            & gl%mmean_PCSAFT + 6.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(6)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)* &
            & gl%mmean_PCSAFT + gl%c_PCSAFT(2)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(9)*gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(5)* &
            & gl%mmean_PCSAFT + 3.d0*gl%c_PCSAFT(2)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(2)*gl%i2_PCSAFT(9)*gl%meo2_PCSAFT(1)* &
            & gl%mmean_PCSAFT + 2.d0*gl%i1_PCSAFT(1)*gl%meo1_PCSAFT(9) + 6.d0*gl%i1_PCSAFT(4)*gl%meo1_PCSAFT(5) + 6.d0*gl%i1_PCSAFT(5)*gl%meo1_PCSAFT(4) + 2.d0* &
            & gl%i1_PCSAFT(9)*gl%meo1_PCSAFT(1) + 2.d0*gl%i1_PCSAFT(14)*gl%meo1_PCSAFT(1) + 6.d0*gl%i1_PCSAFT(7)*gl%meo1_PCSAFT(4) + 6.d0* &
            & gl%i1_PCSAFT(6)*gl%meo1_PCSAFT(5) + 2.d0*gl%i1_PCSAFT(2)*gl%meo1_PCSAFT(9))
    end if

    ! 15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
    if (GETDERADISP(15) .eq. 1) then
        gl%adisp_PCSAFT(15) = -piPCSAFT*D*(gl%c_PCSAFT(1)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(15)*gl%mmean_PCSAFT + 4.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(9)*gl%mmean_PCSAFT + 6.d0* &
            & gl%c_PCSAFT(1)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 4.d0*gl%c_PCSAFT(1)*gl%i2_PCSAFT(9)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + gl%c_PCSAFT(1)* &
            & gl%i2_PCSAFT(15)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(9)*gl%mmean_PCSAFT + 12.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(4) &
            & *gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 12.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + 4.d0*gl%c_PCSAFT(4)*gl%i2_PCSAFT(9)* &
            & gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 6.d0*gl%c_PCSAFT(5)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(5)*gl%mmean_PCSAFT + 12.d0*gl%c_PCSAFT(5)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(4)* &
            & gl%mmean_PCSAFT + 6.d0*gl%c_PCSAFT(5)*gl%i2_PCSAFT(5)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 4.d0*gl%c_PCSAFT(9)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(4)*gl%mmean_PCSAFT + &
            & 4.d0*gl%c_PCSAFT(9)*gl%i2_PCSAFT(4)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + gl%c_PCSAFT(15)*gl%i2_PCSAFT(1)*gl%meo2_PCSAFT(1)*gl%mmean_PCSAFT + 2.d0*gl%i1_PCSAFT(1)* &
            & gl%meo1_PCSAFT(15) + 8.d0*gl%i1_PCSAFT(4)*gl%meo1_PCSAFT(9) + 12.d0*gl%i1_PCSAFT(5)*gl%meo1_PCSAFT(5) + 8.d0*gl%i1_PCSAFT(9)* &
            & gl%meo1_PCSAFT(4) + 2.d0*gl%i1_PCSAFT(15)*gl%meo1_PCSAFT(1))
    end if
    
    !DEC$ END IF
end subroutine ADISPDERIVS

subroutine CDERIVS(gl,GETDERC)

! Henning Markgraf, June 2016

    ! C_1: abbreviation for the compressibility expression
    ! defined by eq. A.11 and A.31 in Gross, Sadowski 2001:
    ! C_1 = (1 + mmean*(8*eta-2*eta**2)/(1-eta)**4 + (1-mmean) &
    !       & * (20*eta-27*eta**2+12*eta**3-2*eta**4)/((1-eta)*(2-eta))**2)**(-1)
    ! dependent on T and D





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERC
    !output: c_PCSAFT (module variable)
    !working variables
    double precision :: recurring_factor
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! recurring_factor occurs in all derivatives
    recurring_factor = 2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 16.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 51.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 70.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 27.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 12.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - gl%z3_PCSAFT(1)** &
            & 6 + 8.d0*gl%z3_PCSAFT(1)**5 - 27.d0*gl%z3_PCSAFT(1)**4 + 42.d0*gl%z3_PCSAFT(1)**3 - 26.d0*gl%z3_PCSAFT(1)**2 + 4.d0
    
    !calculate the derivatives C_1, C_2, C_3, ...
    ! 1: C_1
    ! requires z3_PCSAFT(1)
    if (GETDERC(1) .eq. 1) then
        gl%c_PCSAFT(1) = 1.d0/(gl%mmean_PCSAFT*(-2.d0*gl%z3_PCSAFT(1)**2 + 8.d0*gl%z3_PCSAFT(1))/(-gl%z3_PCSAFT(1) + 1.d0)**4 + (-gl%mmean_PCSAFT + &
             & 1.d0)*(-2.d0*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z3_PCSAFT(1)**3 - 27.d0*gl%z3_PCSAFT(1)**2 + 20.d0*gl%z3_PCSAFT(1))/(( &
             & -gl%z3_PCSAFT(1) + 1.d0)**2*(-gl%z3_PCSAFT(1) + 2.d0)**2) + 1.d0)
    end if
    
    ! 2: 1ST DERIVATIVE OF C_1 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires z3_PCSAFT(1)
    if (GETDERC(2) .eq. 1) then
        ! relation eta_D*D = eta is used
        gl%c_PCSAFT(2) = -2.d0*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 26.d0*gl%mmean_PCSAFT &
                     & *gl%z3_PCSAFT(1)**4 + 115.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 186.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 96.d0*gl%mmean_PCSAFT &
                     & *gl%z3_PCSAFT(1) + 12.d0*gl%mmean_PCSAFT + gl%z3_PCSAFT(1)**5 + 4.d0*gl%z3_PCSAFT(1)**4 - 35.d0*gl%z3_PCSAFT(1)**3 + 74.d0* &
                     & gl%z3_PCSAFT(1)**2 - 64.d0*gl%z3_PCSAFT(1) + 20.d0)/recurring_factor**2
    end if
    
    ! 3: 2ND DERIVATIVE OF C_1 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires z3_PCSAFT(1)
    if (GETDERC(3) .eq. 1) then
        ! relation eta_D*D = eta is used
        gl%c_PCSAFT(3) = 2.d0*(20.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**4 - 16.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 2.d0)/(gl%z3_PCSAFT(1) - 1.d0)**3 + 2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2/( &
            & gl%z3_PCSAFT(1) - 1.d0)**2 - 3.d0*gl%z3_PCSAFT(1)**3*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0* &
            & gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) - 4.d0 &
            & *gl%z3_PCSAFT(1)**3*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - &
            & 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**3*(gl%z3_PCSAFT(1) - 1.d0)) - 3.d0*gl%z3_PCSAFT(1)**3*(gl%mmean_PCSAFT - 1.d0)*(2.d0* &
            & gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/(gl%z3_PCSAFT(1) - 2.d0)**4 + 4.d0* &
            & gl%z3_PCSAFT(1)**2*(gl%mmean_PCSAFT - 1.d0)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - &
            & 10.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)) + 4.d0*gl%z3_PCSAFT(1)**2*(gl%mmean_PCSAFT - 1.d0)*(4.d0* &
            & gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0)/(gl%z3_PCSAFT(1) - 2.d0)**3 + 4.d0* &
            & gl%z3_PCSAFT(1)**2*(-4.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**3 + 2.d0* &
            & gl%mmean_PCSAFT*(gl%z3_PCSAFT(1) - 2.d0)/(gl%z3_PCSAFT(1) - 1.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0* &
            & gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0)) + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + &
            & 27.d0*gl%z3_PCSAFT(1) - 20.d0)/(gl%z3_PCSAFT(1) - 2.d0)**3 - (gl%mmean_PCSAFT - 1.d0)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0* &
            & gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0)/(gl%z3_PCSAFT(1) - 2.d0)**2)**2/((gl%z3_PCSAFT(1) - 1.d0)**2* &
            & (-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**4 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - &
            & 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2 &
            & *(gl%z3_PCSAFT(1) - 1.d0)**2) + 1.d0)) - (gl%mmean_PCSAFT - 1.d0)*(12.d0*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z3_PCSAFT(1) &
            & **3 + 27.d0*gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 2.d0)**2)/((gl%z3_PCSAFT(1) - 1.d0)**2*(-2.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**4 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0* &
            & gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2) + 1.d0)**2)
    end if
    
    ! 4: 1ST DERIVATIVE OF C_1 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires z3_PCSAFT(1), z3_PCSAFT(4)
    if (GETDERC(4) .eq. 1) then
        gl%c_PCSAFT(4) = (-4.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)*(-2.d0*gl%z3_PCSAFT(1)**2 + 8.d0*gl%z3_PCSAFT(1))/(-gl%z3_PCSAFT(1) + 1.d0)**5 - &
            & gl%mmean_PCSAFT*(-4.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 8.d0*gl%z3_PCSAFT(4))/(-gl%z3_PCSAFT(1) + 1.d0)**4 - 2.d0* &
            & gl%z3_PCSAFT(4)*(-gl%mmean_PCSAFT + 1.d0)*(-2.d0*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z3_PCSAFT(1)**3 - 27.d0*gl%z3_PCSAFT(1)** &
            & 2 + 20.d0*gl%z3_PCSAFT(1))/((-gl%z3_PCSAFT(1) + 1.d0)**2*(-gl%z3_PCSAFT(1) + 2.d0)**3) - 2.d0*gl%z3_PCSAFT(4)*( &
            & -gl%mmean_PCSAFT + 1.d0)*(-2.d0*gl%z3_PCSAFT(1)**4 + 12.d0*gl%z3_PCSAFT(1)**3 - 27.d0*gl%z3_PCSAFT(1)**2 + 20.d0* &
            & gl%z3_PCSAFT(1))/((-gl%z3_PCSAFT(1) + 1.d0)**3*(-gl%z3_PCSAFT(1) + 2.d0)**2) - (-gl%mmean_PCSAFT + 1.d0)*(-8.d0* &
            & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 36.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 54.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + &
            & 20.d0*gl%z3_PCSAFT(4))/((-gl%z3_PCSAFT(1) + 1.d0)**2*(-gl%z3_PCSAFT(1) + 2.d0)**2))/(gl%mmean_PCSAFT*(-2.d0* &
            & gl%z3_PCSAFT(1)**2 + 8.d0*gl%z3_PCSAFT(1))/(-gl%z3_PCSAFT(1) + 1.d0)**4 + (-gl%mmean_PCSAFT + 1.d0)*(-2.d0*gl%z3_PCSAFT(1) &
            & **4 + 12.d0*gl%z3_PCSAFT(1)**3 - 27.d0*gl%z3_PCSAFT(1)**2 + 20.d0*gl%z3_PCSAFT(1))/((-gl%z3_PCSAFT(1) + 1.d0)**2* &
            & (-gl%z3_PCSAFT(1) + 2.d0)**2) + 1.d0)**2
    end if
    
    ! 5: 2ND DERIVATIVE OF C_1 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2 
    ! requires z3_PCSAFT(1), z3_PCSAFT(4), z3_PCSAFT(5)
    if (GETDERC(5) .eq. 1) then
        gl%c_PCSAFT(5) = 2.d0*(20.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**4 - 4.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**3 - 16.d0*gl%mmean_PCSAFT &
            & *gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 2.d0)/(gl%z3_PCSAFT(1) - 1.d0)**3 + 2.d0*gl%mmean_PCSAFT*(gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(5) + gl%z3_PCSAFT(4)**2 - 2.d0*gl%z3_PCSAFT(5))/(gl%z3_PCSAFT(1) - 1.d0)**2 - 3.d0* &
            & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0* &
            & gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) - 4.d0*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)**2*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) &
            & - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**3*(gl%z3_PCSAFT(1) - 1.d0)) - 3.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*( &
            & gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/(gl%z3_PCSAFT(1) &
            & - 2.d0)**4 + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1) &
            & **2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)) + gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(5)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - &
            & 20.d0)/(gl%z3_PCSAFT(1) - 2.d0)**3 + 4.d0*gl%z3_PCSAFT(4)**2*(gl%mmean_PCSAFT - 1.d0)*(4.d0*gl%z3_PCSAFT(1)**3 - &
            & 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)) + 4.d0 &
            & *gl%z3_PCSAFT(4)**2*(gl%mmean_PCSAFT - 1.d0)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) &
            & - 10.d0)/(gl%z3_PCSAFT(1) - 2.d0)**3 + 4.d0*gl%z3_PCSAFT(4)**2*(-4.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - &
            & 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**3 + 2.d0*gl%mmean_PCSAFT*(gl%z3_PCSAFT(1) - 2.d0)/(gl%z3_PCSAFT(1) - 1.d0)**2 + &
            & gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) &
            & /((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)) + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)** &
            & 3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/(gl%z3_PCSAFT(1) - 2.d0)**3 - (gl%mmean_PCSAFT - 1.d0) &
            & *(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0)/(gl%z3_PCSAFT(1) - 2.d0)**2)** &
            & 2/((gl%z3_PCSAFT(1) - 1.d0)**2*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**4 &
            & + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - &
            & 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) + 1.d0)) - (gl%mmean_PCSAFT - 1.d0)*(4.d0* &
            & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 12.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 18.d0*gl%z3_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(5) - 36.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 27.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 27.d0* &
            & gl%z3_PCSAFT(4)**2 - 10.d0*gl%z3_PCSAFT(5))/(gl%z3_PCSAFT(1) - 2.d0)**2)/((gl%z3_PCSAFT(1) - 1.d0)**2*( &
            & -2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**4 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - &
            & 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2 &
            & *(gl%z3_PCSAFT(1) - 1.d0)**2) + 1.d0)**2)
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF C_1 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires z3_PCSAFT(1), z3_PCSAFT(4)
    if (GETDERC(6) .eq. 1) then
        !DERZETA(6) = DERZETA(4) was used
        gl%c_PCSAFT(6) = 2.d0*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**2*(4.d0*gl%z3_PCSAFT(1)*(4.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 2.d0 &
            & )**3*(gl%z3_PCSAFT(1) - 4.d0) + 2.d0*gl%mmean_PCSAFT*(-gl%z3_PCSAFT(1) + 2.d0)**4*(gl%z3_PCSAFT(1) - 1.d0) - &
            & gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(-gl%z3_PCSAFT(1) + 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - &
            & 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0) &
            & **3*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%mmean_PCSAFT - 1.d0)* &
            & (-gl%z3_PCSAFT(1) + 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**3*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0* &
            & gl%z3_PCSAFT(1) - 10.d0))**2 + (-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 2.d0)**2*(gl%z3_PCSAFT(1) - 4.d0 &
            & ) + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1) &
            & **2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (-gl%z3_PCSAFT(1) + 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)*(20.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*(-gl%z3_PCSAFT(1) + 2.d0)**4*(gl%z3_PCSAFT(1) - 4.d0) + 4.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)* &
            & (-5.d0*gl%z3_PCSAFT(1) + 12.d0)*(-gl%z3_PCSAFT(1) + 2.d0)**4*(gl%z3_PCSAFT(1) - 1.d0) + 4.d0*gl%mmean_PCSAFT*( &
            & -gl%z3_PCSAFT(1) + 2.d0)**4*(gl%z3_PCSAFT(1) - 1.d0)**3 - 3.d0*gl%z3_PCSAFT(1)**2*(gl%mmean_PCSAFT - 1.d0)*( &
            & -gl%z3_PCSAFT(1) + 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0* &
            & gl%z3_PCSAFT(1) - 20.d0) + 4.d0*gl%z3_PCSAFT(1)**2*(gl%mmean_PCSAFT - 1.d0)*(-gl%z3_PCSAFT(1) + 2.d0)*(gl%z3_PCSAFT(1) - 1.d0 &
            & )**3*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) - 3.d0*gl%z3_PCSAFT(1)**2* &
            & (gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**4*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0* &
            & gl%z3_PCSAFT(1) - 20.d0) + 3.d0*gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(-gl%z3_PCSAFT(1) + 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0 &
            & )**3*(6.d0*gl%z3_PCSAFT(1)**3 - 28.d0*gl%z3_PCSAFT(1)**2 + 45.d0*gl%z3_PCSAFT(1) - 20.d0) + 3.d0*gl%z3_PCSAFT(1)*( &
            & gl%mmean_PCSAFT - 1.d0)*(-gl%z3_PCSAFT(1) + 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**4*(-6.d0*gl%z3_PCSAFT(1)**3 + 28.d0* &
            & gl%z3_PCSAFT(1)**2 - 45.d0*gl%z3_PCSAFT(1) + 20.d0) - 2.d0*(gl%mmean_PCSAFT - 1.d0)*(-gl%z3_PCSAFT(1) + 2.d0)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0)**4*(8.d0*gl%z3_PCSAFT(1)**3 - 27.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 5.d0)))/(-2.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 2.d0)**2*(gl%z3_PCSAFT(1) - 4.d0) + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0 &
            & )*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + &
            & (-gl%z3_PCSAFT(1) + 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)**3
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF C_1 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires z3_PCSAFT(1), z3_PCSAFT(4), z3_PCSAFT(5)
    if (GETDERC(7) .eq. 1) then
        gl%c_PCSAFT(7) = 2.d0*(8.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - 24.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20* &
            & gl%z3_PCSAFT(4)**2 - 436.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 1488.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + 8032.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - &
            & 26624.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 81156.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 &
            & *gl%z3_PCSAFT(5) + 252424.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 534174.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 1549074.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(4)**2 - 2493493.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 6810320.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + 8666464.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15* &
            & gl%z3_PCSAFT(5) - 22805576.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 23104428.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 60415872.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13* &
            & gl%z3_PCSAFT(4)**2 + 48083876.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 129127636.d0 &
            & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 78714334.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(5) + 223518640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + &
            & 101046784.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 310510056.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 99978008.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) &
            & + 339461080.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 73240270.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 283338378.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - &
            & 36344673.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 172334112.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + 9373392.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - &
            & 70465056.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 716232.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
            & **6*gl%z3_PCSAFT(5) + 16293168.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - &
            & 1123560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 1146888.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 4*gl%z3_PCSAFT(4)**2 + 136296.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 48672.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 31104.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(5) - 62208.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 3456.d0*gl%mmean_PCSAFT**3 &
            & *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 6912.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 360.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) - 1440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
            & gl%z3_PCSAFT(4)**2 - 8640.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 31632.d0*gl%mmean_PCSAFT** &
            & 2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 97110.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - &
            & 322524.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 678354.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 17*gl%z3_PCSAFT(5) + 2047014.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + 3297213.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 9159936.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15* &
            & gl%z3_PCSAFT(4)**2 - 11836992.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 31025256.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + 32524968.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(5) - 83146632.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 69927876.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 180875700.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(4)**2 + 119101458.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - &
            & 322609296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 161260800.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + 470731032.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 &
            & + 172404018.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 555914028.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 142541058.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + &
            & 522618126.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 87204753.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 382590432.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - &
            & 35812368.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 211605408.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 7155504.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - &
            & 84392064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + 1102536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
            & **5*gl%z3_PCSAFT(5) + 22235496.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 - 997128.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 3059808.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)**2 + 174720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 18432.d0*gl%mmean_PCSAFT &
            & **2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 576.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + &
            & 40320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 4608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**2 - 6.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 18.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - &
            & 33.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 324.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 &
            & + 2616.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - 11736.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18* &
            & gl%z3_PCSAFT(4)**2 - 36258.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 134700.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + 278028.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - &
            & 899892.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - 1429902.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(5) + 4142424.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + 5352816.d0*gl%mmean_PCSAFT &
            & *gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 14306808.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - &
            & 15255744.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 39042168.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13* &
            & gl%z3_PCSAFT(4)**2 + 34024818.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 86863422.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 60394209.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(5) + 159997284.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 86049096.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 244648896.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4) &
            & **2 - 98437470.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 308691564.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 89537880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - &
            & 317765040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 63221112.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8* &
            & gl%z3_PCSAFT(5) + 262809792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + 32972832.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 171194880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)** &
            & 2 - 11295912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 85447920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 5*gl%z3_PCSAFT(4)**2 + 1529376.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 31342560.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 629808.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + &
            & 7866048.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 377088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(5) - 1163520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 74304.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 60288.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 3456.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 4224.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**2 - 384.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(5) + 2.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) - 6.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - &
            & 19.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 12.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 152.d0*gl%z3_PCSAFT(1) &
            & **19*gl%z3_PCSAFT(5) + 1288.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 3996.d0*gl%z3_PCSAFT(1)**18* &
            & gl%z3_PCSAFT(5) - 18272.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 36204.d0*gl%z3_PCSAFT(1)**17* &
            & gl%z3_PCSAFT(5) + 131820.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + 202466.d0*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(5) - 632392.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 801680.d0*gl%z3_PCSAFT(1)**15* &
            & gl%z3_PCSAFT(5) + 2246680.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 + 2391612.d0*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(5) - 6280416.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 5568166.d0*gl%z3_PCSAFT(1)**13* &
            & gl%z3_PCSAFT(5) + 14331770.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + 10338117.d0*gl%z3_PCSAFT(1)** &
            & 12*gl%z3_PCSAFT(5) - 27229796.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 15499528.d0*gl%z3_PCSAFT(1) &
            & **11*gl%z3_PCSAFT(5) + 43355472.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 18853488.d0* &
            & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 57691008.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 18539360.d0 &
            & *gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 63628400.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 + 14530884.d0 &
            & *gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 57501888.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 8806880.d0* &
            & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + 41975872.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 3872624.d0* &
            & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 24303520.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 1037184.d0* &
            & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 10884096.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 25840.d0* &
            & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 3632192.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 107776.d0* &
            & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 849920.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 45568.d0* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 124416.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 8576.d0*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(5) + 8576.d0*gl%z3_PCSAFT(4)**2 - 640.d0*gl%z3_PCSAFT(5))/(16.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**24 - 512.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 + 7776.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 - &
            & 74176.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 + 495928.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 - 2456160.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 + 9285352.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 - 27195504.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 + 61924305.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 - 108773576.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 + 144073828.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 - 136765176.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 + 81765502.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 - 16168680.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 - 16180068.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 + 11992104.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9 - 314415.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 - 1968624.d0*gl%mmean_PCSAFT**4 &
            & *gl%z3_PCSAFT(1)**7 + 146016.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 + 186624.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
            & 5 + 20736.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - 32.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 1024.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - 15600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 149728.d0*gl%mmean_PCSAFT**3 &
            & *gl%z3_PCSAFT(1)**21 - 1010520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 5070624.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**19 - 19506196.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 58453880.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**17 - 137174640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 250911504.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**15 - 351689504.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 363681264.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**13 - 253952496.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + 87886064.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**11 + 23665428.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 39522472.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**9 + 11476728.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 4037184.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
            & **7 - 2432592.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 243648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + &
            & 186624.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 27648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 24.d0*gl%mmean_PCSAFT** &
            & 2.d0*gl%z3_PCSAFT(1)**24 - 768.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 11736.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 22 - 113328.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 + 771942.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - &
            & 3923064.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + 15347616.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - &
            & 47009184.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + 113501556.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - &
            & 215552784.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + 317971272.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - &
            & 354036960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 + 279237222.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - &
            & 130019160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 + 3195072.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + &
            & 39411120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 21887832.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + &
            & 465696.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 + 3588576.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 763776.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 271008.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 62208.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**3 + 13824.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 8.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 + 256.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 3924.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 38120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 &
            & - 262016.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 + 1348176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 5360672.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 + 16767120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 41588880.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**16 + 81795408.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 126414924.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 14 + 150209544.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 130931368.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
            & 74326592.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 15325152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 14073184.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 13226976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 3227520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **7 - 1321984.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 890240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 24960.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 77824.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 6912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + &
            & 3072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + gl%z3_PCSAFT(1)**24 - 32.d0*gl%z3_PCSAFT(1)**23 + 492.d0*gl%z3_PCSAFT(1)**22 - &
            & 4808.d0*gl%z3_PCSAFT(1)**21 + 33342.d0*gl%z3_PCSAFT(1)**20 - 173640.d0*gl%z3_PCSAFT(1)**19 + 701372.d0* &
            & gl%z3_PCSAFT(1)**18 - 2238280.d0*gl%z3_PCSAFT(1)**17 + 5695617.d0*gl%z3_PCSAFT(1)**16 - 11575192.d0* &
            & gl%z3_PCSAFT(1)**15 + 18672048.d0*gl%z3_PCSAFT(1)**14 - 23512368.d0*gl%z3_PCSAFT(1)**13 + 22306712.d0 &
            & *gl%z3_PCSAFT(1)**12 - 14677120.d0*gl%z3_PCSAFT(1)**11 + 4997280.d0*gl%z3_PCSAFT(1)**10 + 1273152.d0* &
            & gl%z3_PCSAFT(1)**9 - 2481104.d0*gl%z3_PCSAFT(1)**8 + 1105152.d0*gl%z3_PCSAFT(1)**7 + 22656.d0*gl%z3_PCSAFT(1) &
            & **6 - 207616.d0*gl%z3_PCSAFT(1)**5 + 57984.d0*gl%z3_PCSAFT(1)**4 + 10752.d0*gl%z3_PCSAFT(1)**3 - 6656.d0 &
            & *gl%z3_PCSAFT(1)**2 + 256.d0)
    end if

    ! 8: 3RD DERIVATIVE OF C_1 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires z3_PCSAFT(1)
    if (GETDERC(8) .eq. 1) then
        ! relation eta_D*D = eta is used
        gl%c_PCSAFT(8) = 2.d0*(-120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**5 + 120.d0*gl%mmean_PCSAFT &
            & *gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 2.d0)/(gl%z3_PCSAFT(1) - 1.d0)**4 - 24.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3/( &
            & gl%z3_PCSAFT(1) - 1.d0)**3 + 12.d0*gl%z3_PCSAFT(1)**4*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0* &
            & gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**3) + &
            & 18.d0*gl%z3_PCSAFT(1)**4*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) &
            & - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**3*(gl%z3_PCSAFT(1) - 1.d0)**2) + 18.d0*gl%z3_PCSAFT(1)**4*(gl%mmean_PCSAFT - &
            & 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**4 &
            & *(gl%z3_PCSAFT(1) - 1.d0)) + 12.d0*gl%z3_PCSAFT(1)**4*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0* &
            & gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/(gl%z3_PCSAFT(1) - 2.d0)**5 - 18.d0*gl%z3_PCSAFT(1)**3*( &
            & gl%mmean_PCSAFT - 1.d0)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0)/(( &
            & gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) - 24.d0*gl%z3_PCSAFT(1)**3*(gl%mmean_PCSAFT - 1.d0)*(4.d0* &
            & gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0)/((gl%z3_PCSAFT(1) - 2.d0)**3*( &
            & gl%z3_PCSAFT(1) - 1.d0)) - 18.d0*gl%z3_PCSAFT(1)**3*(gl%mmean_PCSAFT - 1.d0)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1) &
            & **2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0)/(gl%z3_PCSAFT(1) - 2.d0)**4 + 24.d0*gl%z3_PCSAFT(1)**3*(-4.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**3 + 2.d0*gl%mmean_PCSAFT*(gl%z3_PCSAFT(1) - 2.d0)/( &
            & gl%z3_PCSAFT(1) - 1.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 &
            & + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)) + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT &
            & - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/(gl%z3_PCSAFT(1) - 2.d0)** &
            & 3 - (gl%mmean_PCSAFT - 1.d0)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0)/( &
            & gl%z3_PCSAFT(1) - 2.d0)**2)**3/((gl%z3_PCSAFT(1) - 1.d0)**4*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0 &
            & )/(gl%z3_PCSAFT(1) - 1.d0)**4 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1) &
            & **2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) + 1.d0)**2) &
            & + 6.d0*gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(12.d0*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z3_PCSAFT(1)**3 + 27.d0*gl%z3_PCSAFT(1) &
            & **2)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)) + 6.d0*gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(12.d0* &
            & gl%z3_PCSAFT(1)**4 - 36.d0*gl%z3_PCSAFT(1)**3 + 27.d0*gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 2.d0)**3 - 12.d0* &
            & gl%z3_PCSAFT(1)*(-4.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**3 + 2.d0*gl%mmean_PCSAFT* &
            & (gl%z3_PCSAFT(1) - 2.d0)/(gl%z3_PCSAFT(1) - 1.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - &
            & 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)) + &
            & gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) &
            & /(gl%z3_PCSAFT(1) - 2.d0)**3 - (gl%mmean_PCSAFT - 1.d0)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0* &
            & gl%z3_PCSAFT(1) - 10.d0)/(gl%z3_PCSAFT(1) - 2.d0)**2)*(-20.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 4.d0)/ &
            & (gl%z3_PCSAFT(1) - 1.d0)**4 + 16.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 2.d0)/(gl%z3_PCSAFT(1) - 1.d0)** &
            & 3 - 2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1) - 1.d0)**2 + 3.d0*gl%z3_PCSAFT(1)**3*(gl%mmean_PCSAFT - 1.d0) &
            & *(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2) + 4.d0*gl%z3_PCSAFT(1)**3*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0* &
            & gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**3*(gl%z3_PCSAFT(1) - 1.d0)) + 3.d0* &
            & gl%z3_PCSAFT(1)**3*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - &
            & 20.d0)/(gl%z3_PCSAFT(1) - 2.d0)**4 - 4.d0*gl%z3_PCSAFT(1)**2*(gl%mmean_PCSAFT - 1.d0)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0* &
            & gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)) - 4.d0* &
            & gl%z3_PCSAFT(1)**2*(gl%mmean_PCSAFT - 1.d0)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - &
            & 10.d0)/(gl%z3_PCSAFT(1) - 2.d0)**3 + (gl%mmean_PCSAFT - 1.d0)*(12.d0*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z3_PCSAFT(1)**3 + &
            & 27.d0*gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 2.d0)**2)/((gl%z3_PCSAFT(1) - 1.d0)**2*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & *(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0)**4 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 &
            & - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**2 &
            & ) + 1.d0)) - (gl%mmean_PCSAFT - 1.d0)*(24.d0*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z3_PCSAFT(1)**3)/(gl%z3_PCSAFT(1) - 2.d0) &
            & **2)/((gl%z3_PCSAFT(1) - 1.d0)**2*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)/(gl%z3_PCSAFT(1) - 1.d0) &
            & **4 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) &
            & - 20.d0)/((gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) + 1.d0)**2)
    end if
    
    ! 9: 3RD DERIVATIVE OF C_1 WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires z3_PCSAFT(1), z3_PCSAFT(4), z3_PCSAFT(5), z3_PCSAFT(9)
    if (GETDERC(9) .eq. 1) then
        gl%c_PCSAFT(9) = 2.d0*(gl%z3_PCSAFT(1) - 1.d0)*(-120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) &
            & - 2.d0)**5*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*( &
            & gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0* &
            & gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)**2 + 60.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 2.d0)**5*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(gl%z3_PCSAFT(1) &
            & - 4.d0) + gl%z3_PCSAFT(4)**2*(2.d0*gl%z3_PCSAFT(1) - 4.d0))*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)* &
            & (gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 &
            & - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)** &
            & 4)**2 + 2.d0*gl%mmean_PCSAFT*(gl%z3_PCSAFT(1) - 2.d0)**5*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(9) + 3.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 2.d0*gl%z3_PCSAFT(9))*(-2.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) &
            & - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - &
            & 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)**2 + 4.d0*gl%mmean_PCSAFT*(gl%z3_PCSAFT(1) - 2.d0)**5*(gl%z3_PCSAFT(1) - 1.d0) &
            & **2*(-gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9)*(gl%z3_PCSAFT(1) - 4.d0) - 6.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*( &
            & gl%z3_PCSAFT(1) - 2.d0) - 6.d0*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + gl%z3_PCSAFT(4)**2 - 2.d0* &
            & gl%z3_PCSAFT(5)))*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + &
            & gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + &
            & 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)**2 + 12.d0*gl%z3_PCSAFT(1) &
            & *gl%z3_PCSAFT(4)**3*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)**3*(gl%z3_PCSAFT(1) - 1.d0)**2*(-2.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + ( &
            & gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + &
            & 27.d0*gl%z3_PCSAFT(1) - 20.d0) + 18.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0) &
            & **2*(gl%z3_PCSAFT(1) - 1.d0)**3*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 &
            & + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 &
            & + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)**2*(2.d0*gl%z3_PCSAFT(1) &
            & **3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + 18.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3*( &
            & gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**4*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) &
            & - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0* &
            & gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - 2.d0)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0)**4)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) &
            & + 12.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**5*(-2.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) &
            & - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - &
            & 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) &
            & - 20.d0) + 24.d0*gl%z3_PCSAFT(4)**3*(-4.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0) &
            & **3 + 2.d0*gl%mmean_PCSAFT*(gl%z3_PCSAFT(1) - 2.d0)**4*(gl%z3_PCSAFT(1) - 1.d0) + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*( &
            & gl%z3_PCSAFT(1) - 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0* &
            & gl%z3_PCSAFT(1) - 20.d0) + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**3*(2.d0*gl%z3_PCSAFT(1)**3 - &
            & 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) - (gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)*(gl%z3_PCSAFT(1) &
            & - 1.d0)**3*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0))**3 + 9.d0* &
            & gl%z3_PCSAFT(4)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)**3*(gl%z3_PCSAFT(1) - 1.d0)**3*(-gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(5)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + &
            & gl%z3_PCSAFT(4)**2*(-8.d0*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z3_PCSAFT(1)**2 - 54.d0*gl%z3_PCSAFT(1) + 20.d0))*(-2.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + ( &
            & gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)**2 + 12.d0*gl%z3_PCSAFT(4)*(gl%mmean_PCSAFT - 1.d0)*( &
            & gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(-gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(2.d0*gl%z3_PCSAFT(1)**3 - &
            & 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + gl%z3_PCSAFT(4)**2*(-8.d0*gl%z3_PCSAFT(1)**3 + 36.d0* &
            & gl%z3_PCSAFT(1)**2 - 54.d0*gl%z3_PCSAFT(1) + 20.d0))*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*( &
            & gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 &
            & - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)** &
            & 4)**2 + 9.d0*gl%z3_PCSAFT(4)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**5*( &
            & -gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) &
            & + gl%z3_PCSAFT(4)**2*(-8.d0*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z3_PCSAFT(1)**2 - 54.d0*gl%z3_PCSAFT(1) + 20.d0))*(-2.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + ( &
            & gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)**2 + 12.d0*gl%z3_PCSAFT(4)*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & *(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)** &
            & 2.d0*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - 2.d0)**2 &
            & *(gl%z3_PCSAFT(1) - 1.d0)**4)*(-4.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**3 + &
            & 2.d0*gl%mmean_PCSAFT*(gl%z3_PCSAFT(1) - 2.d0)**4*(gl%z3_PCSAFT(1) - 1.d0) + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) &
            & - 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0 &
            & ) + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**3*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)** &
            & 2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) - (gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**3*(4.d0 &
            & *gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0))*(20.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**4 + 2.d0*gl%mmean_PCSAFT*(gl%z3_PCSAFT(1) - 2.d0)** &
            & 4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + gl%z3_PCSAFT(4)**2 - 2.d0*gl%z3_PCSAFT(5) &
            & ) + 4.d0*gl%mmean_PCSAFT*(gl%z3_PCSAFT(1) - 2.d0)**4*(gl%z3_PCSAFT(1) - 1.d0)*(-gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*( &
            & gl%z3_PCSAFT(1) - 4.d0) + gl%z3_PCSAFT(4)**2*(-4.d0*gl%z3_PCSAFT(1) + 8.d0)) - 3.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 &
            & *(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0* &
            & gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) - 4.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(gl%mmean_PCSAFT - 1.d0)*( &
            & gl%z3_PCSAFT(1) - 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**3*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0* &
            & gl%z3_PCSAFT(1) - 20.d0) - 3.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**4*( &
            & 2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) - (gl%mmean_PCSAFT - 1.d0)*( &
            & gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(4.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 12.d0* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 18.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 36.d0*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)**2 + 27.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 27.d0*gl%z3_PCSAFT(4)**2 - 10.d0* &
            & gl%z3_PCSAFT(5)) + (gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z3_PCSAFT(1) &
            & *gl%z3_PCSAFT(5)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + &
            & gl%z3_PCSAFT(4)**2*(16.d0*gl%z3_PCSAFT(1)**3 - 72.d0*gl%z3_PCSAFT(1)**2 + 108.d0*gl%z3_PCSAFT(1) - 40.d0)) + ( &
            & gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(2.d0* &
            & gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + gl%z3_PCSAFT(4)**2*(16.d0* &
            & gl%z3_PCSAFT(1)**3 - 72.d0*gl%z3_PCSAFT(1)**2 + 108.d0*gl%z3_PCSAFT(1) - 40.d0))) - (gl%mmean_PCSAFT - 1.d0)*( &
            & gl%z3_PCSAFT(1) - 2.d0)**3*(gl%z3_PCSAFT(1) - 1.d0)**5*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*( &
            & gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 &
            & - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)** &
            & 4)**2*(4.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + 36.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 18.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 24.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 - &
            & 108.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 27.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) - 36.d0* &
            & gl%z3_PCSAFT(4)**3 + 81.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 10.d0*gl%z3_PCSAFT(9)) + (gl%mmean_PCSAFT &
            & - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)**3*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9)*(2.d0* &
            & gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + 6.d0*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 10.d0) + 6.d0* &
            & gl%z3_PCSAFT(4)*(4.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 12.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 18.d0* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 36.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 27.d0*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(5) + 27.d0*gl%z3_PCSAFT(4)**2 - 10.d0*gl%z3_PCSAFT(5)))*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*( &
            & gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2* &
            & (2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - 2.d0)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0)**4)**2 + (gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**5* &
            & (gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9)*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0 &
            & ) + 6.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(4.d0*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) &
            & - 10.d0) + 6.d0*gl%z3_PCSAFT(4)*(4.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 12.d0*gl%z3_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(4)**2 - 18.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 36.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + &
            & 27.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 27.d0*gl%z3_PCSAFT(4)**2 - 10.d0*gl%z3_PCSAFT(5)))*(-2.d0*gl%mmean_PCSAFT &
            & *gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) &
            & - 1.d0)**2*(2.d0*gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - &
            & 2.d0)**2*(gl%z3_PCSAFT(1) - 1.d0)**4)**2)/((gl%z3_PCSAFT(1) - 2.d0)*(-2.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) &
            & - 4.d0)*(gl%z3_PCSAFT(1) - 2.d0)**2 + gl%z3_PCSAFT(1)*(gl%mmean_PCSAFT - 1.d0)*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0* &
            & gl%z3_PCSAFT(1)**3 - 12.d0*gl%z3_PCSAFT(1)**2 + 27.d0*gl%z3_PCSAFT(1) - 20.d0) + (gl%z3_PCSAFT(1) - 2.d0)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0)**4)**4)
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF C_1 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires z3_PCSAFT(1), z3_PCSAFT(4)
    if (GETDERC(10) .eq. 1) then
        gl%c_PCSAFT(10) = -4.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(12.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - 744.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 19 + 13312.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - 126212.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + &
            & 774537.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - 3405160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + &
            & 11402788.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 - 30207936.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + &
            & 64563818.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 - 111759320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 + &
            & 155255028.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 169730540.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + &
            & 141669189.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 - 86167056.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + &
            & 35232528.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 - 8146584.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 573444.d0 &
            & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 24336.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 31104.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**2 + 3456.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) + 720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - &
            & 15816.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 161262.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - 1023507.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 4579968.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - 15512628.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 41573316.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 - 90437850.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 161304648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11 - 235365516.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 277957014.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9 - 261309063.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 191295216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7 - 105802704.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 42196032.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - 11117748.d0*gl%mmean_PCSAFT &
            & **2*gl%z3_PCSAFT(1)**4 + 1529904.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 9216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
            & **2 - 20160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - 2304.d0*gl%mmean_PCSAFT**2 - 9.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 &
            & - 162.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 + 5868.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 67350.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**17 + 449946.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 2071212.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + &
            & 7153404.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 19521084.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + 43431711.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 79998642.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + 122324448.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**10 - 154345782.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 158882520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 8 - 131404896.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + 85597440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - &
            & 42723960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + 15671280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 3933024.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 581760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 30144.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - &
            & 2112.d0*gl%mmean_PCSAFT + 3.d0*gl%z3_PCSAFT(1)**20 - 6.d0*gl%z3_PCSAFT(1)**19 - 644.d0*gl%z3_PCSAFT(1)**18 + 9136.d0* &
            & gl%z3_PCSAFT(1)**17 - 65910.d0*gl%z3_PCSAFT(1)**16 + 316196.d0*gl%z3_PCSAFT(1)**15 - 1123340.d0* &
            & gl%z3_PCSAFT(1)**14 + 3140208.d0*gl%z3_PCSAFT(1)**13 - 7165885.d0*gl%z3_PCSAFT(1)**12 + 13614898.d0* &
            & gl%z3_PCSAFT(1)**11 - 21677736.d0*gl%z3_PCSAFT(1)**10 + 28845504.d0*gl%z3_PCSAFT(1)**9 - 31814200.d0* &
            & gl%z3_PCSAFT(1)**8 + 28750944.d0*gl%z3_PCSAFT(1)**7 - 20987936.d0*gl%z3_PCSAFT(1)**6 + 12151760.d0* &
            & gl%z3_PCSAFT(1)**5 - 5442048.d0*gl%z3_PCSAFT(1)**4 + 1816096.d0*gl%z3_PCSAFT(1)**3 - 424960.d0*gl%z3_PCSAFT(1) &
            & **2 + 62208.d0*gl%z3_PCSAFT(1) - 4288.d0)/(16.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 - 512.d0*gl%mmean_PCSAFT**4 &
            & *gl%z3_PCSAFT(1)**23 + 7776.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 - 74176.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
            & 21 + 495928.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 - 2456160.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 + &
            & 9285352.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 - 27195504.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 + &
            & 61924305.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 - 108773576.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 + &
            & 144073828.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 - 136765176.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 + &
            & 81765502.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 - 16168680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 - &
            & 16180068.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 + 11992104.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9 - &
            & 314415.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 - 1968624.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 + 146016.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 + 186624.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 + 20736.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**4 - 32.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 1024.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 - &
            & 15600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 149728.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - 1010520.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 5070624.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - 19506196.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 58453880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - 137174640.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 250911504.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 - 351689504.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 363681264.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 - 253952496.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + 87886064.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 + 23665428.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 - 39522472.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 11476728.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 4037184.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - 2432592.d0*gl%mmean_PCSAFT** &
            & 3*gl%z3_PCSAFT(1)**6 - 243648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 186624.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
            & **4 + 27648.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 + 24.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 768.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 + 11736.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 113328.d0*gl%mmean_PCSAFT**2 &
            & *gl%z3_PCSAFT(1)**21 + 771942.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 3923064.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**19 + 15347616.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 47009184.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**17 + 113501556.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 215552784.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**15 + 317971272.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 354036960.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**13 + 279237222.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 130019160.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**11 + 3195072.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 39411120.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**9 - 21887832.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 465696.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
            & **7 + 3588576.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 763776.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5 - &
            & 271008.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 62208.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3 + 13824.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 8.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 + 256.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - &
            & 3924.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 38120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 262016.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**20 + 1348176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 5360672.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 &
            & + 16767120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 41588880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + &
            & 81795408.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 126414924.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
            & 150209544.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 130931368.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
            & 74326592.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 15325152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 14073184.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + 13226976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 3227520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **7 - 1321984.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + 890240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 24960.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 77824.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 6912.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + &
            & 3072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + gl%z3_PCSAFT(1)**24 - 32.d0*gl%z3_PCSAFT(1)**23 + 492.d0*gl%z3_PCSAFT(1)**22 - &
            & 4808.d0*gl%z3_PCSAFT(1)**21 + 33342.d0*gl%z3_PCSAFT(1)**20 - 173640.d0*gl%z3_PCSAFT(1)**19 + 701372.d0* &
            & gl%z3_PCSAFT(1)**18 - 2238280.d0*gl%z3_PCSAFT(1)**17 + 5695617.d0*gl%z3_PCSAFT(1)**16 - 11575192.d0* &
            & gl%z3_PCSAFT(1)**15 + 18672048.d0*gl%z3_PCSAFT(1)**14 - 23512368.d0*gl%z3_PCSAFT(1)**13 + 22306712.d0 &
            & *gl%z3_PCSAFT(1)**12 - 14677120.d0*gl%z3_PCSAFT(1)**11 + 4997280.d0*gl%z3_PCSAFT(1)**10 + 1273152.d0* &
            & gl%z3_PCSAFT(1)**9 - 2481104.d0*gl%z3_PCSAFT(1)**8 + 1105152.d0*gl%z3_PCSAFT(1)**7 + 22656.d0*gl%z3_PCSAFT(1) &
            & **6 - 207616.d0*gl%z3_PCSAFT(1)**5 + 57984.d0*gl%z3_PCSAFT(1)**4 + 10752.d0*gl%z3_PCSAFT(1)**3 - 6656.d0 &
            & *gl%z3_PCSAFT(1)**2 + 256.d0)
    end if
    
    ! 11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^4
    ! requires z3_PCSAFT(1)
    if (GETDERC(11) .eq. 1) then
        gl%c_PCSAFT(11) = 24.d0*gl%z3_PCSAFT(1)**4*(40.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24 - 2800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 + &
            & 64560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22 - 823360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 + 6981434.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 - 43204900.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 + 206213750.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 - 785933560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 + 2444418380.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 - 6282058840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 + &
            & 13412911340.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 - 23787155240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 &
            & + 34858008390.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 - 41774875700.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
            & 11 + 40262377234.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 - 30394270600.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**9 + 17209360320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8 - 6750206160.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**7 + 1525468260.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6 - 80184456.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**5 - 19320300.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 - 10562400.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**3 + 1779840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 + 933120.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) + &
            & 82944.d0*gl%mmean_PCSAFT**4 - 20.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24 + 3800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 23 - 103020.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 + 1402400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 - &
            & 12307903.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 + 77833870.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 - &
            & 377527265.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 + 1459730820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 - &
            & 4609716450.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 + 12063014756.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 &
            & - 26353966730.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 + 48167003760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 13 - 73488180395.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12 + 93036008470.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**11 - 96768354961.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10 + 81478052780.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**9 - 54349602420.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8 + 27808494720.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 - 10366403640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6 + 2564709792.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 - 336047820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 + 2460000.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**3 + 3480480.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 + 63360.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1) - 48384.d0*gl%mmean_PCSAFT**3 - 30.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 - 1500.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**23 + 57780.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 - 870960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 21 + 8020329.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 - 52191390.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 + &
            & 258395535.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 - 1016872980.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 + &
            & 3268903560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 - 8728750344.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 + &
            & 19542317190.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 - 36831700920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 &
            & + 58439761125.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 - 77814647790.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 11 + 86421328395.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 - 79310963580.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**9 + 59346227640.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 - 35524349280.d0*gl%mmean_PCSAFT**2 &
            & *gl%z3_PCSAFT(1)**7 + 16542352020.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 - 5737121112.d0*gl%mmean_PCSAFT**2 &
            & *gl%z3_PCSAFT(1)**5 + 1373794920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 - 192877440.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**3 + 8242080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 + 940800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) + &
            & 5952.d0*gl%mmean_PCSAFT**2 + 25.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 + 50.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 - 12795.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 + 230840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21 - 2280112.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**20 + 15418280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 - 78355210.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 18 + 315029480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 - 1033767955.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 + &
            & 2822334298.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 - 6482953575.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 + &
            & 12600134080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 - 20756290970.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 + &
            & 28941729740.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 - 34026640036.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 + &
            & 33516243600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 - 27402555260.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 + &
            & 18354886960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 - 9889271840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 + &
            & 4173978976.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 - 1326060880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 + &
            & 296960640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 - 41344000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 + 2551040.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1) + 34624.d0*gl%mmean_PCSAFT - 5.d0*gl%z3_PCSAFT(1)**24 + 50.d0*gl%z3_PCSAFT(1)**23 + 795.d0* &
            & gl%z3_PCSAFT(1)**22 - 21520.d0*gl%z3_PCSAFT(1)**21 + 237122.d0*gl%z3_PCSAFT(1)**20 - 1690220.d0* &
            & gl%z3_PCSAFT(1)**19 + 8881270.d0*gl%z3_PCSAFT(1)**18 - 36642000.d0*gl%z3_PCSAFT(1)**17 + 123096495.d0 &
            & *gl%z3_PCSAFT(1)**16 - 344227742.d0*gl%z3_PCSAFT(1)**15 + 811836695.d0*gl%z3_PCSAFT(1)**14 - &
            & 1626323640.d0*gl%z3_PCSAFT(1)**13 + 2775625820.d0*gl%z3_PCSAFT(1)**12 - 4035886360.d0*gl%z3_PCSAFT(1) &
            & **11 + 4988054584.d0*gl%z3_PCSAFT(1)**10 - 5216783960.d0*gl%z3_PCSAFT(1)**9 + 4586032920.d0* &
            & gl%z3_PCSAFT(1)**8 - 3356872320.d0*gl%z3_PCSAFT(1)**7 + 2019468640.d0*gl%z3_PCSAFT(1)**6 - &
            & 980565280.d0*gl%z3_PCSAFT(1)**5 + 374478240.d0*gl%z3_PCSAFT(1)**4 - 108217600.d0*gl%z3_PCSAFT(1)**3 + &
            & 22230080.d0*gl%z3_PCSAFT(1)**2 - 2890880.d0*gl%z3_PCSAFT(1) + 178816.d0)/recurring_factor**5
    end if
    
    ! 12: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, DEL AND DEL, MULTIPLIED BY TAU*DEL^3
    ! requires z3_PCSAFT(1), z3_PCSAFT(4)
    if (GETDERC(12) .eq. 1) then
        gl%c_PCSAFT(12) = 24.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*(16.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25 - 1240.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**24 + 28200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23 - 345400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
            & 22 + 2787416.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21 - 16436122.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20 + &
            & 75472340.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19 - 281266990.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18 + &
            & 872201720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17 - 2277215500.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16 + &
            & 5016529280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15 - 9286935080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14 + &
            & 14337790800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13 - 18267548750.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12 &
            & + 18939928756.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 - 15669224254.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
            & 10 + 10048168200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9 - 4773384420.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
            & **8 + 1557015600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7 - 303207804.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
            & 6 + 27586152.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5 - 2640600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 + &
            & 444960.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3 + 233280.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2 + 20736.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) - 8.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25 + 1820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 24 - 47640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23 + 614300.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22 - &
            & 5072452.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21 + 30219907.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20 - &
            & 139389530.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 + 521042205.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 - &
            & 1624043220.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17 + 4283915594.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16 - &
            & 9607323164.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15 + 18285023730.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14 &
            & - 29381707940.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13 + 39585453835.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 12 - 44338529914.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11 + 40843339517.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**10 - 30494834580.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9 + 18074830020.d0*gl%mmean_PCSAFT** &
            & 3*gl%z3_PCSAFT(1)**8 - 8239931880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7 + 2738412168.d0*gl%mmean_PCSAFT**3 &
            & *gl%z3_PCSAFT(1)**6 - 594260376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5 + 59472120.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**4 + 3510720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3 - 429120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 2 - 245376.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) - 20736.d0*gl%mmean_PCSAFT**3 - 12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
            & **25 - 870.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24 + 29250.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23 - &
            & 405750.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22 + 3458256.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21 - &
            & 20925903.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20 + 97330860.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 - &
            & 366066555.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 + 1149824760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17 - &
            & 3070763076.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16 + 7021470930.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15 - &
            & 13744998660.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14 + 22951046400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13 &
            & - 32518658175.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12 + 38855212800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 11 - 38857687287.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10 + 32204496900.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**9 - 21816704160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8 + 11842176720.d0*gl%mmean_PCSAFT**2 &
            & *gl%z3_PCSAFT(1)**7 - 4997624988.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6 + 1562458608.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**5 - 331952640.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4 + 39444000.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**3 - 968640.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2 - 189312.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) - &
            & 3456.d0*gl%mmean_PCSAFT**2 + 10.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 + 125.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24 - 7620.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23 + 117725.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22 - 1047598.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**21 + 6478370.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20 - 30517090.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19 &
            & + 115860680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18 - 367593130.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17 + &
            & 995281297.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16 - 2320892328.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15 + &
            & 4667327605.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14 - 8071877570.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13 + &
            & 11952019760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12 - 15073911994.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11 + &
            & 16096493958.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10 - 14445188000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 + &
            & 10787655280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8 - 6614548520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7 + &
            & 3267864304.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6 - 1265918752.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5 + &
            & 368931360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4 - 75493120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3 + 9463040.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2 - 483584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) - 14208.d0*gl%mmean_PCSAFT - 2.d0*gl%z3_PCSAFT(1)** &
            & 25 + 5.d0*gl%z3_PCSAFT(1)**24 + 690.d0*gl%z3_PCSAFT(1)**23 - 12625.d0*gl%z3_PCSAFT(1)**22 + 119078.d0* &
            & gl%z3_PCSAFT(1)**21 - 757940.d0*gl%z3_PCSAFT(1)**20 + 3632500.d0*gl%z3_PCSAFT(1)**19 - 13964790.d0* &
            & gl%z3_PCSAFT(1)**18 + 44833170.d0*gl%z3_PCSAFT(1)**17 - 123113363.d0*gl%z3_PCSAFT(1)**16 + &
            & 292460834.d0*gl%z3_PCSAFT(1)**15 - 602609565.d0*gl%z3_PCSAFT(1)**14 + 1074784250.d0*gl%z3_PCSAFT(1)** &
            & 13 - 1652831350.d0*gl%z3_PCSAFT(1)**12 + 2181833656.d0*gl%z3_PCSAFT(1)**11 - 2460311612.d0* &
            & gl%z3_PCSAFT(1)**10 + 2356433760.d0*gl%z3_PCSAFT(1)**9 - 1903119240.d0*gl%z3_PCSAFT(1)**8 + &
            & 1283601040.d0*gl%z3_PCSAFT(1)**7 - 713577760.d0*gl%z3_PCSAFT(1)**6 + 321060192.d0*gl%z3_PCSAFT(1)**5 &
            & - 113915200.d0*gl%z3_PCSAFT(1)**4 + 30659840.d0*gl%z3_PCSAFT(1)**3 - 5879360.d0*gl%z3_PCSAFT(1)**2 + &
            & 715264.d0*gl%z3_PCSAFT(1) - 41472.d0)/recurring_factor**5
    end if
    
    ! 13: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, DEL, AND DEL, MULTIPLIED BY TAU^2*DEL^2
    ! requires z3_PCSAFT(1), z3_PCSAFT(4), z3_PCSAFT(5)
    if (GETDERC(13) .eq. 1) then
        gl%c_PCSAFT(13) = -4.d0*(24.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(5) - 72.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26* &
            & gl%z3_PCSAFT(4)**2 - 1680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 5760.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 + 39140.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - &
            & 130060.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 504200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
            & 24*gl%z3_PCSAFT(5) + 1568200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 4299782.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 12424714.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(4)**2 - 26591508.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) + 72025224.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 125974859.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21* &
            & gl%z3_PCSAFT(5) - 326859181.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 473989210.d0 &
            & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 1213612730.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19* &
            & gl%z3_PCSAFT(4)**2 + 1451755955.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - &
            & 3781454365.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 3677984500.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 9985308500.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 &
            & + 7772982770.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 22326192910.d0*gl%mmean_PCSAFT**4 &
            & *gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - 13721514924.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(5) + 42000095556.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + &
            & 20120913700.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 65905831100.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 24205886460.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(5) + 85399406040.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + &
            & 23406178307.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 90233394229.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 17611040946.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(5) + 76404304578.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + &
            & 9748346703.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 50540662497.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 3517359420.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) &
            & + 25122947100.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + 516452940.d0*gl%mmean_PCSAFT**4 &
            & *gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 8825640660.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 &
            & + 163441872.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 1982688696.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 82448532.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - &
            & 247965444.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + 5537376.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 21380976.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + &
            & 889920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 1779840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4* &
            & gl%z3_PCSAFT(4)**2 + 466560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 933120.d0*gl%mmean_PCSAFT &
            & **4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 41472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - &
            & 82944.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 12.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27* &
            & gl%z3_PCSAFT(5) + 36.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)**2 + 2280.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 8640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - &
            & 62740.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + 223100.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(4)**2 + 865600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5) - 2820200.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 7659439.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23* &
            & gl%z3_PCSAFT(5) + 22775273.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + 48474598.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 132844844.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21* &
            & gl%z3_PCSAFT(4)**2 - 233413868.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 602923312.d0 &
            & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 889915470.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20* &
            & gl%z3_PCSAFT(5) - 2236337760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - &
            & 2760694105.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 6983565215.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 + 7097945568.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) &
            & - 18605547996.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 15292609536.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + 42351329448.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(4)**2 + 27728015284.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - &
            & 81982127096.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - 42226178005.d0*gl%mmean_PCSAFT** &
            & 3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 134064069635.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(4)**2 + 53592767470.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - &
            & 183919955540.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 55931990428.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 210099189056.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(4)**2 + 47037958958.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - &
            & 198022078144.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 - 30902255091.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + 152066752389.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10* &
            & gl%z3_PCSAFT(4)**2 + 15034101220.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - &
            & 93414878900.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - 4813162740.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 44626428540.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 &
            & + 616556016.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - 15813916992.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + 226098396.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) + &
            & 3791660652.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 125748000.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 482580720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 + &
            & 21531600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 467280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4 &
            & *gl%z3_PCSAFT(4)**2 - 264960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 2309760.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 179712.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(5) + 1292544.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 13824.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 110592.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 18.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(5) + 54.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4) &
            & **2 - 900.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 4320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 25*gl%z3_PCSAFT(4)**2 + 35445.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) - 140055.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 543450.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(5) + 1891050.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 5058732.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) - 15690804.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(4)**2 - 32996742.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) + 92558676.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 162235719.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21* &
            & gl%z3_PCSAFT(5) - 421749441.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 628883400.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5) + 1567515930.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
            & gl%z3_PCSAFT(4)**2 + 1981121460.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) - &
            & 4917827100.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 5178993072.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 13245585384.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)** &
            & 2 + 11387345895.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 30741479685.d0*gl%mmean_PCSAFT &
            & **2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 - 21204344346.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(5) + 61265647614.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 + &
            & 33468115440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 104238162960.d0*gl%mmean_PCSAFT**2 &
            & *gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2 - 44577645510.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(5) + 150534303540.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + &
            & 49647118605.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - 183484158195.d0*gl%mmean_PCSAFT**2 &
            & *gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 45588108708.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(5) + 187558015014.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + &
            & 33808748394.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 159418233006.d0*gl%mmean_PCSAFT**2 &
            & *gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - 19606017120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10* &
            & gl%z3_PCSAFT(5) + 111294207840.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 + &
            & 8378553120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 62674507200.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2 - 2273394816.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + &
            & 27712355112.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 + 145595352.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 9229156296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 + &
            & 169556544.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 2161272384.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 - 74981520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - &
            & 311645520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 + 12862080.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 18673920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - &
            & 321984.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 813888.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(4)**2 - 105984.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 85248.d0*gl%mmean_PCSAFT** &
            & 2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 9216.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 9216.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**2 + 15.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(5) - 45.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)**2 + 30.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) - 720.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2 - 7960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5) + &
            & 37760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 + 146350.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(5) - 560000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 - 1464091.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) + 4821497.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 + &
            & 9937740.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5) - 28932480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
            & gl%z3_PCSAFT(4)**2 - 50177152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5) + 132925388.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 + 198551720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20* &
            & gl%z3_PCSAFT(5) - 496612360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 - 637077415.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5) + 1568481365.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) &
            & **2 + 1696961414.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) - 4274726368.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 - 3812048392.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) + &
            & 10113305576.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2 + 7288111966.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) - 20715853664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 - &
            & 11896522965.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) + 36534742455.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **14*gl%z3_PCSAFT(4)**2 + 16546011280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) - &
            & 55166107280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 - 19481516848.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) + 70961955116.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 + &
            & 19213595692.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5) - 77365368056.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **11*gl%z3_PCSAFT(4)**2 - 15630347624.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) + &
            & 71040780376.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 + 10251192160.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) - 54474739520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2 - &
            & 5219318240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) + 34467972880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 8*gl%z3_PCSAFT(4)**2 + 1911003328.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) - &
            & 17696182496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2 - 397114848.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **7*gl%z3_PCSAFT(5) + 7198397664.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2 - 26687872.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) - 2240276032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4) &
            & **2 + 52096320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) + 505055040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **4*gl%z3_PCSAFT(4)**2 - 18156800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) - 74935040.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 3012672.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + &
            & 5914176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 172032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(5) - 86784.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 8448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(5) - 8448.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**2 - 3.d0*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(5) + 9.d0* &
            & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)**2 + 30.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) + 515.d0*gl%z3_PCSAFT(1)**25 &
            & *gl%z3_PCSAFT(5) - 3625.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2 - 14000.d0*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(5) + 61750.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2 + 156056.d0*gl%z3_PCSAFT(1)**23* &
            & gl%z3_PCSAFT(5) - 558412.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2 - 1117040.d0*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(5) + 3430600.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2 + 5832946.d0*gl%z3_PCSAFT(1)**21* &
            & gl%z3_PCSAFT(5) - 15962054.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2 - 23670000.d0*gl%z3_PCSAFT(1)**20 &
            & *gl%z3_PCSAFT(5) + 60118740.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2 + 77609045.d0*gl%z3_PCSAFT(1)** &
            & 19*gl%z3_PCSAFT(5) - 191389975.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2 - 211092426.d0* &
            & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5) + 527587752.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2 + &
            & 484907751.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5) - 1269857253.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)** &
            & 2 - 951217432.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5) + 2664439958.d0*gl%z3_PCSAFT(1)**15* &
            & gl%z3_PCSAFT(4)**2 + 1601522470.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) - 4847183030.d0*gl%z3_PCSAFT(1) &
            & **14*gl%z3_PCSAFT(4)**2 - 2313984580.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5) + 7603003520.d0* &
            & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2 + 2856447652.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5) - &
            & 10234554284.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2 - 2988050648.d0*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(5) + 11773819024.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2 + 2617328304.d0* &
            & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5) - 11521274256.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2 - &
            & 1887085840.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5) + 9531629600.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)** &
            & 2 + 1090692480.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5) - 6610913760.d0*gl%z3_PCSAFT(1)**8* &
            & gl%z3_PCSAFT(4)**2 - 482004480.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5) + 3799462080.d0*gl%z3_PCSAFT(1)** &
            & 7*gl%z3_PCSAFT(4)**2 + 145793408.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5) - 1780567744.d0*gl%z3_PCSAFT(1) &
            & **6*gl%z3_PCSAFT(4)**2 - 18173696.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) + 665317504.d0*gl%z3_PCSAFT(1) &
            & **5*gl%z3_PCSAFT(4)**2 - 7990720.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5) - 191949760.d0*gl%z3_PCSAFT(1) &
            & **4*gl%z3_PCSAFT(4)**2 + 5466880.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5) + 40743040.d0*gl%z3_PCSAFT(1)** &
            & 3*gl%z3_PCSAFT(4)**2 - 1588352.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 5879936.d0*gl%z3_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(4)**2 + 248832.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 497664.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
            & **2 - 17152.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 17152.d0*gl%z3_PCSAFT(4)**2)/recurring_factor**5
    end if
    
    ! 14: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, TAU AND DEL, MULTIPLIED BY TAU^3*DEL
    ! requires z3_PCSAFT(1), z3_PCSAFT(4), z3_PCSAFT(5), z3_PCSAFT(9)
    if (GETDERC(14) .eq. 1) then
        gl%c_PCSAFT(14) = 2.d0*(16.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(9) - 144.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) + 192.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**3 + 10080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 23448.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(9) - 14880.d0 &
            & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 - 234840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 313620.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(9) + &
            & 338400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 + 3025200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
            & **23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 2807212.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(9) &
            & - 4144800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 - 25798692.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 18246642.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(9) + 33448992.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 + 159549048.d0 &
            & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 90364242.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(9) - 197233464.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 &
            & - 755849154.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 351527431.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) + 905668080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19* &
            & gl%z3_PCSAFT(4)**3 + 2843935260.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
            & 1095821600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) - 3375203880.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 - 8710535730.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 2772663215.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(9) + &
            & 10466420640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 + 22067907000.d0*gl%mmean_PCSAFT** &
            & 4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 5735183160.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17* &
            & gl%z3_PCSAFT(9) - 27326586000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**3 - &
            & 46637896620.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 9716828902.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(9) + 60198351360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15* &
            & gl%z3_PCSAFT(4)**3 + 82329089544.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
            & + 13430529548.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) - 111443220960.d0*gl%mmean_PCSAFT &
            & **4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**3 - 120725482200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 14964967460.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(9) + &
            & 172053489600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 + 145235318760.d0*gl%mmean_PCSAFT &
            & **4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 13117667042.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
            & 13*gl%z3_PCSAFT(9) - 219210585000.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**3 - &
            & 140437069842.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 8615783839.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) + 227279145072.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11 &
            & *gl%z3_PCSAFT(4)**3 + 105666245676.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) + 3786214464.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) - &
            & 188030691048.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 - 58490080218.d0*gl%mmean_PCSAFT &
            & **4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 703782987.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10 &
            & *gl%z3_PCSAFT(9) + 120578018400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**3 + &
            & 21104156520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 292610820.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) - 57280613040.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8* &
            & gl%z3_PCSAFT(4)**3 - 3098717640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
            & 216928512.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) + 18684187200.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 - 980651232.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 29751048.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(9) - 3638493648.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 + 494691192.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 11803752.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) + &
            & 331033824.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3 - 33224256.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 2233440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5* &
            & gl%z3_PCSAFT(9) - 31687200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 - 5339520.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 466560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4 &
            & *gl%z3_PCSAFT(9) + 5339520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 - 2799360.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 41472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(9) + 2799360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 - 248832.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 248832.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)**3 - 8.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(9) + 72.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1220.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26* &
            & gl%z3_PCSAFT(9) - 96.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**3 - 13680.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 34776.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25* &
            & gl%z3_PCSAFT(9) + 21840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 + 376440.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 508340.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(9) - 571680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 - 5193600.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 4795114.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
            & **23*gl%z3_PCSAFT(9) + 7371600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 + &
            & 45956634.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 32321997.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(9) - 60869424.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) &
            & **3 - 290847588.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
            & 164706564.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(9) + 362638884.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 + 1400483208.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) + 656906812.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) - 1672674360.d0 &
            & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 - 5339492820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19 &
            & *gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 2097393910.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) + &
            & 6252506460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 + 16564164630.d0*gl%mmean_PCSAFT**3 &
            & *gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 5440998823.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18* &
            & gl%z3_PCSAFT(9) - 19488518640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 - &
            & 42587673408.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 11573850044.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(9) + 51406987128.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(4)**3 + 91755657216.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
            & + 20274139992.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(9) - 115287877968.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**3 - 166368091704.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 29225546750.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) + &
            & 219420284760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**3 + 253357068030.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 34432325155.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 14*gl%z3_PCSAFT(9) - 352580495280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 - &
            & 321556604820.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 32641175388.d0 &
            & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(9) + 475025446020.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 12*gl%z3_PCSAFT(4)**3 + 335591942568.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) + 24148643532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) - &
            & 532062358968.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**3 - 282227753748.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 13095571138.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 11*gl%z3_PCSAFT(9) + 490120074204.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 + &
            & 185413530546.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 4406466157.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) - 365938014960.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9* &
            & gl%z3_PCSAFT(4)**3 - 90204607320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
            & 238099340.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) + 216897960240.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 + 28878976440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 582000084.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) - 98879182560.d0 &
            & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 - 3699336096.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 265953384.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(9) + &
            & 32860946016.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 - 1356590376.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 25358760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6* &
            & gl%z3_PCSAFT(9) - 7131124512.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3 + &
            & 754488000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 12446208.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) + 713665440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4) &
            & **3 - 129189600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 2567520.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(9) + 42128640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)**3 + 1589760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
            & 131328.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) - 5149440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2 &
            & *gl%z3_PCSAFT(4)**3 + 1078272.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
            & 13824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) - 2944512.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)**3 + 82944.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 248832.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4)**3 - 12.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(9) + 108.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 330.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26* &
            & gl%z3_PCSAFT(9) - 144.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**3 + 5400.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 16974.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25* &
            & gl%z3_PCSAFT(9) - 10440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 - 212670.d0*gl%mmean_PCSAFT &
            & **2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 291585.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(9) + 351000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 + 3260700.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 2975382.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
            & **23*gl%z3_PCSAFT(9) - 4869000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 - &
            & 30352392.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 21059748.d0*gl%mmean_PCSAFT &
            & **2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(9) + 41499072.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4) &
            & **3 + 197980452.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
            & 111205302.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(9) - 251110836.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 - 973414314.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 456751023.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(9) + 1167970320.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 + 3773300400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19 &
            & *gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1498111950.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) - &
            & 4392798660.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 - 11886728760.d0*gl%mmean_PCSAFT**2 &
            & *gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 3992387952.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18* &
            & gl%z3_PCSAFT(9) + 13797897120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 + &
            & 31073958432.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 8742162690.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(9) - 36849156912.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(4)**3 - 68324075370.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
            & - 15830834583.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(9) + 84257651160.d0*gl%mmean_PCSAFT &
            & **2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**3 + 127226066076.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 23754627810.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) - &
            & 164939983920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**3 - 200808692640.d0*gl%mmean_PCSAFT &
            & **2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 29447075160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 14*gl%z3_PCSAFT(9) + 275412556800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 + &
            & 267465873060.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 29869230330.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(9) - 390223898100.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 12*gl%z3_PCSAFT(4)**3 - 297882711630.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 24308313837.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) + &
            & 466262553600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**3 + 273528652248.d0*gl%mmean_PCSAFT &
            & **2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 15273858102.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 11*gl%z3_PCSAFT(9) - 466292247444.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 - &
            & 202852490364.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 6805123770.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) + 386453962800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
            & gl%z3_PCSAFT(4)**3 + 117636102720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
            & + 1613062320.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) - 261800449920.d0*gl%mmean_PCSAFT** &
            & 2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 - 50271318720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) &
            & *gl%z3_PCSAFT(5) + 259332324.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) + &
            & 142106120640.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 + 13640368896.d0*gl%mmean_PCSAFT** &
            & 2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 377279232.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
            & gl%z3_PCSAFT(9) - 59971499856.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 - &
            & 873572112.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 127467552.d0*gl%mmean_PCSAFT &
            & **2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) + 18749503296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5* &
            & gl%z3_PCSAFT(4)**3 - 1017339264.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
            & 8103456.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) - 3983431680.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 + 449889120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 6300000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(9) + 473328000.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 - 77172480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1524096.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) - &
            & 11623680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 + 1931904.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
            & **2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 49536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) - &
            & 2271744.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 + 635904.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 4608.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) - 41472.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**3 + 55296.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 10.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(9) - 90.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 85.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) + 120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25 &
            & *gl%z3_PCSAFT(4)**3 - 180.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 2616.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(9) + 1500.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 + &
            & 47760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 67140.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 24*gl%z3_PCSAFT(9) - 91440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**3 - 878100.d0*gl%mmean_PCSAFT &
            & *gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 782666.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23* &
            & gl%z3_PCSAFT(9) + 1412700.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 + 8784546.d0*gl%mmean_PCSAFT &
            & *gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 5941965.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(9) - 12571176.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 - 59626440.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 32862936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
            & gl%z3_PCSAFT(9) + 77740440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 + 301062912.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 139881596.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20* &
            & gl%z3_PCSAFT(9) - 366205080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 - 1191310320.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 473260930.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19* &
            & gl%z3_PCSAFT(9) + 1390328160.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 + 3822464490.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1299237369.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 18*gl%z3_PCSAFT(9) - 4411117560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**3 - &
            & 10181768484.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 2933689208.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(9) + 11943375564.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(4)**3 + 22872290352.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
            & 5495114348.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(9) - 27850707936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **15*gl%z3_PCSAFT(4)**3 - 43728671796.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 8574673798.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) + 56007931260.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**3 + 71379137790.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 11145121875.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(9) - &
            & 96862530840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 - 99276067680.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 12001604728.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13* &
            & gl%z3_PCSAFT(9) + 143424237120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**3 + &
            & 116889101088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 10570972388.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) - 180886943928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11* &
            & gl%z3_PCSAFT(4)**3 - 115281574152.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - &
            & 7425908312.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) + 193157927496.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 + 93782085744.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) + 3951946396.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) - 173342256000.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**3 - 61507152960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1397041120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) + &
            & 129451863360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 + 31315909440.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 155930688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8* &
            & gl%z3_PCSAFT(9) - 79374582240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**3 - &
            & 11466019968.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 152260288.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(9) + 39214371648.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 + &
            & 2382689088.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 101185136.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) - 15191025024.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3 + &
            & 160127232.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 25990208.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) + 4427176320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**3 - &
            & 312577920.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 117440.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(9) - 905917440.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 + &
            & 108940800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1705088.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + 113556480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 - &
            & 18076032.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 392832.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **2*gl%z3_PCSAFT(9) - 5803008.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 + 1032192.d0*gl%mmean_PCSAFT &
            & *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 21504.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) - &
            & 170496.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**3 + 50688.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1536.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(9) - 2.d0*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(9) + 18.d0*gl%z3_PCSAFT(1)**26* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 35.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(9) - 24.d0*gl%z3_PCSAFT(1)**25* &
            & gl%z3_PCSAFT(4)**3 - 180.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 54.d0*gl%z3_PCSAFT(1)**25* &
            & gl%z3_PCSAFT(9) + 60.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**3 - 3090.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4) &
            & *gl%z3_PCSAFT(5) - 4615.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(9) + 8280.d0*gl%z3_PCSAFT(1)**23* &
            & gl%z3_PCSAFT(4)**3 + 84000.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 71426.d0*gl%z3_PCSAFT(1) &
            & **23*gl%z3_PCSAFT(9) - 151500.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**3 - 936336.d0*gl%z3_PCSAFT(1)** &
            & 22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 605880.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(9) + 1428936.d0* &
            & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**3 + 6702240.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
            & 3570708.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(9) - 9095280.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**3 - &
            & 34997676.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 15896174.d0*gl%z3_PCSAFT(1)**20* &
            & gl%z3_PCSAFT(9) + 43590000.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**3 + 142020000.d0*gl%z3_PCSAFT(1)** &
            & 19*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 55790690.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(9) - &
            & 167577480.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**3 - 465654270.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 158375661.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(9) + 537998040.d0*gl%z3_PCSAFT(1)** &
            & 17*gl%z3_PCSAFT(4)**3 + 1266554556.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
            & 369691514.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(9) - 1477360356.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4) &
            & **3 - 2909446506.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 717213891.d0*gl%z3_PCSAFT(1) &
            & **16*gl%z3_PCSAFT(9) + 3509530008.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**3 + 5707304592.d0* &
            & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1163621030.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(9) &
            & - 7231314780.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**3 - 9609134820.d0*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1582094710.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(9) + &
            & 12897411000.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**3 + 13883907480.d0*gl%z3_PCSAFT(1)**13* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1798178232.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(9) - &
            & 19833976200.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**3 - 17138685912.d0*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 1694152872.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(9) + &
            & 26182003872.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**3 + 17928303888.d0*gl%z3_PCSAFT(1)**11* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1300126312.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(9) - &
            & 29523739344.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**3 - 15703969824.d0*gl%z3_PCSAFT(1)**10* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 785162152.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(9) + &
            & 28277205120.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**3 + 11322515040.d0*gl%z3_PCSAFT(1)**9* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 345574560.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(9) - 22837430880.d0 &
            & *gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**3 - 6544154880.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
            & - 85916320.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(9) + 15403212480.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4) &
            & **3 + 2892026880.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 10458528.d0*gl%z3_PCSAFT(1)** &
            & 7*gl%z3_PCSAFT(9) - 8562933120.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**3 - 874760448.d0* &
            & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 20644832.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(9) + &
            & 3852722304.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**3 + 109042176.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 9101440.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(9) - 1366982400.d0*gl%z3_PCSAFT(1)**4* &
            & gl%z3_PCSAFT(4)**3 + 47944320.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 1665600.d0* &
            & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(9) + 367918080.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3 - 32801280.d0 &
            & *gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 181248.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) - &
            & 70552320.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 + 9530112.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 165632.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 8583168.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
            & **3 - 1492992.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 34304.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) &
            & - 497664.d0*gl%z3_PCSAFT(4)**3 + 102912.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 2560.d0* &
            & gl%z3_PCSAFT(9))/recurring_factor**5
    end if
    
    ! 15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
    ! requires z3_PCSAFT(1), z3_PCSAFT(4), z3_PCSAFT(5), z3_PCSAFT(9), z3_PCSAFT(15)
    if (GETDERC(15) .eq. 1) then
        gl%c_PCSAFT(15) = -2.d0*(8.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(15) - 96.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 72.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5)**2 - 440.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(15) + 576.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 5760.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 4320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5)**2 + 10132.d0*gl%mmean_PCSAFT &
            & **4*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(15) - 480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**4 &
            & - 37440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 134320.d0*gl%mmean_PCSAFT &
            & **4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 100740.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(5)**2 - 138100.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(15) + 33600.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**4 + 872640.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 1806880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 1355160.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5)**2 + 1281042.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(15) - 774720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(4)**4 - 11471040.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) &
            & - 16353016.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 12264762.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5)**2 - 8699406.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(15) + 9880320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**4 + 100656432.d0 &
            & *gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 107784192.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 80838144.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21* &
            & gl%z3_PCSAFT(5)**2 + 45143719.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(15) - &
            & 83777208.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**4 - 642450672.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 542031844.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 406523883.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5)**2 &
            & - 183816499.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(15) + 518458800.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**4 + 3137793840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4) &
            & **2*gl%z3_PCSAFT(5) + 2141375720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 1606031790.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5)**2 + &
            & 596855765.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(15) - 2474565000.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**4 - 12111997680.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4) &
            & **2*gl%z3_PCSAFT(5) - 6770709460.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 5078032095.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5)**2 - &
            & 1559038965.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(15) + 9431202720.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**4 + 37733199840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4) &
            & **2*gl%z3_PCSAFT(5) + 17326808720.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 12995106540.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5)**2 + &
            & 3284598190.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(15) - 29333020560.d0*gl%mmean_PCSAFT**4 &
            & *gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**4 - 96116240160.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 36079125400.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 27059344050.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5)**2 - &
            & 5562096494.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(15) + 75384706080.d0*gl%mmean_PCSAFT**4 &
            & *gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**4 + 201513169440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 61115701584.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 45836776188.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5)**2 + &
            & 7488991932.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(15) - 160954936080.d0*gl%mmean_PCSAFT** &
            & 4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**4 - 348005283840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 83678085920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 62758564440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5)**2 - &
            & 7843799980.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(15) + 285445862880.d0*gl%mmean_PCSAFT** &
            & 4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**4 + 492485222160.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 91235069760.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 68426302320.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5)**2 + &
            & 6121051607.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(15) - 418296100680.d0*gl%mmean_PCSAFT** &
            & 4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**4 - 564175846800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 76954874596.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 57716155947.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5)**2 - &
            & 3223267907.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(15) + 501298508400.d0*gl%mmean_PCSAFT** &
            & 4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**4 + 511738763472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 47356206984.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 35517155238.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5)**2 + &
            & 787823073.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(15) - 483148526808.d0*gl%mmean_PCSAFT**4 &
            & *gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**4 - 353401112304.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 18296150148.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 13722112611.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5)**2 + &
            & 267498927.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(15) + 364731247200.d0*gl%mmean_PCSAFT**4 &
            & *gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**4 + 171868610880.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
            & **2*gl%z3_PCSAFT(5) + 1745136240.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) &
            & + 1308852180.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5)**2 - 286936800.d0*gl%mmean_PCSAFT**4 &
            & *gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(15) - 206512323840.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) &
            & **4 - 47443721760.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 2318190480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 1738642860.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5)**2 + 65676312.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**8* &
            & gl%z3_PCSAFT(15) + 81002473920.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**4 - &
            & 757136160.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 1130419296.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 847814472.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5)**2 + 18386676.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(15) &
            & - 18305619120.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**4 + 5352560352.d0*gl%mmean_PCSAFT** &
            & 4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 45457488.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 34093116.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5)**2 - &
            & 9577224.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(15) + 962213472.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**4 - 1125754848.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2 &
            & *gl%z3_PCSAFT(5) + 85523904.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + &
            & 64142928.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5)**2 - 453600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
            & **5*gl%z3_PCSAFT(15) + 231843600.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**4 - &
            & 190123200.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 7119360.d0* &
            & gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 5339520.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1) &
            & **4*gl%z3_PCSAFT(5)**2 + 466560.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(15) + &
            & 126748800.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**4 + 32037120.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 3732480.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 2799360.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5)**2 + &
            & 41472.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(15) - 21358080.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)** &
            & 2*gl%z3_PCSAFT(4)**4 + 16796160.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) &
            & - 331776.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 248832.d0*gl%mmean_PCSAFT**4 &
            & *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)**2 - 11197440.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**4 + &
            & 1492992.d0*gl%mmean_PCSAFT**4*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 995328.d0*gl%mmean_PCSAFT**4* &
            & gl%z3_PCSAFT(4)**4 - 4.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(15) + 48.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 36.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5) &
            & **2 + 460.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(15) - 288.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 25*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 6720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 5040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5)**2 - 13244.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(15) + 240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**4 &
            & + 47520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 192080.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 144060.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(5)**2 + 202260.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(15) - 45600.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**4 - 1329120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 2842400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 2131800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5)**2 - 2016349.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(15) + 1236240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(4)**4 + 18914400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) &
            & + 27245852.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 20434389.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5)**2 + 14431331.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(15) - 16828800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**4 - &
            & 173650824.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 187013312.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 140259984.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5)**2 - 78155038.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**21* &
            & gl%z3_PCSAFT(15) + 147694836.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**4 + &
            & 1142735112.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 971446408.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 728584806.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5)**2 + 330554948.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**20* &
            & gl%z3_PCSAFT(15) - 934006440.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**4 - &
            & 5715305640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 3949847040.d0 &
            & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 2962385280.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5)**2 - 1113258445.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**19* &
            & gl%z3_PCSAFT(15) + 4530327180.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**4 + &
            & 22528526760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 12842609420.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 9631957065.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5)**2 + 3019253933.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**18 &
            & *gl%z3_PCSAFT(15) - 17516769840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**4 - &
            & 71656157520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 33841011024.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 25380758268.d0 &
            & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5)**2 - 6626827552.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 17*gl%z3_PCSAFT(15) + 55316597400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**4 + &
            & 186698379888.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 72802710384.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 54602032788.d0 &
            & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5)**2 + 11761804784.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 16*gl%z3_PCSAFT(15) - 144756177072.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**4 - &
            & 401919445584.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 128143779104.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - &
            & 96107834328.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5)**2 - 16761220155.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(15) + 316247600760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(4)**4 + 717167520720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) + 183947067620.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) &
            & + 137960300715.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5)**2 + 18876016645.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(15) - 578004045120.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 13*gl%z3_PCSAFT(4)**4 - 1058555338920.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) - 213233367200.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) &
            & - 159925025400.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5)**2 - 16286484278.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(15) + 881858164740.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 12*gl%z3_PCSAFT(4)**4 + 1282813311240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) + 195710638664.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) &
            & + 146782978998.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5)**2 + 10082824036.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(15) - 1116432101640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 11*gl%z3_PCSAFT(4)**4 - 1258315801128.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) - 136925870272.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) &
            & - 102694402704.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5)**2 - 3711709451.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(15) + 1161220259532.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)** &
            & 10*gl%z3_PCSAFT(4)**4 + 975233118312.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) + 67229122356.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + &
            & 50421841767.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5)**2 + 26762023.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(15) - 977736633360.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4) &
            & **4 - 572514428160.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 17732912720.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 13299684540.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5)**2 + 798318560.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**9* &
            & gl%z3_PCSAFT(15) + 652195229040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**4 + &
            & 233607952800.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 2240876880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 1680657660.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5)**2 - 381960684.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**8* &
            & gl%z3_PCSAFT(15) - 333701936640.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**4 - &
            & 51035322240.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 3855843072.d0 &
            & *gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 2891882304.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5)**2 + 24373332.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(15) &
            & + 124396843680.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**4 - 4168857024.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 1161306864.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1) &
            & **6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 870980148.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5) &
            & **2 + 37082520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(15) - 30776517504.d0*gl%mmean_PCSAFT &
            & **3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**4 + 6197101344.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 46895040.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 35171280.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5)**2 - 9145872.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(15) + 4032573840.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4* &
            & gl%z3_PCSAFT(4)**4 - 1368290880.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) &
            & + 86368320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 64776240.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5)**2 - 1317600.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**4* &
            & gl%z3_PCSAFT(15) - 29520000.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**4 - 725760.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 4999680.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 3749760.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(5)**2 + 442368.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(15) - 41765760.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**4 + 11819520.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 2294784.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 1721088.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)**2 + 41472.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(15) - 760320.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
            & **4 + 4727808.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 221184.d0* &
            & gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 165888.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(5)**2 + 580608.d0*gl%mmean_PCSAFT**3*gl%z3_PCSAFT(4)**4 + 497664.d0*gl%mmean_PCSAFT**3* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 6.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(15) + 72.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 54.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26* &
            & gl%z3_PCSAFT(5)**2 - 30.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(15) - 432.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 1440.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 1080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5)**2 + &
            & 4641.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(15) + 360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(4)**4 - 15120.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 86460.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 64845.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5)**2 - 95385.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(15) + &
            & 18000.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**4 + 684720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 23*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 1547880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 1160910.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5)**2 + 1090032.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(15) - 693360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(4)**4 - 11165040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) &
            & - 16261656.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 12196242.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5)**2 - 8474484.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(15) + 10451520.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**4 + &
            & 109489752.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 118136928.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 88602696.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5)**2 + 48688329.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**21* &
            & gl%z3_PCSAFT(15) - 96243948.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**4 - &
            & 750371688.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 639574524.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 479680893.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5)**2 - 215978427.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**20* &
            & gl%z3_PCSAFT(15) + 626296680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**4 + &
            & 3865552200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 2690917800.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 2018188350.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5)**2 + 758924190.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**19* &
            & gl%z3_PCSAFT(15) - 3100746420.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**4 - &
            & 15619354200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 9028144560.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 6771108420.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5)**2 - 2144593512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**18 &
            & *gl%z3_PCSAFT(15) + 12202475760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**4 + &
            & 50857891200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 24547925856.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 18410944392.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5)**2 + 4912183935.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 17*gl%z3_PCSAFT(15) - 39226842720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**4 - &
            & 135791694432.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 54617386500.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 40963039875.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5)**2 - 9138225591.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 16*gl%z3_PCSAFT(15) + 104745004128.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**4 + &
            & 300500310240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 99876240696.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 74907180522.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5)**2 + 13754686200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 15*gl%z3_PCSAFT(15) - 234507806280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**4 - &
            & 554080854240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 150037256040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - &
            & 112527942030.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5)**2 - 16567258560.d0*gl%mmean_PCSAFT &
            & **2*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(15) + 441980411040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13* &
            & gl%z3_PCSAFT(4)**4 + 851729153400.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) + 184057334880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) &
            & + 138043001160.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5)**2 + 15615882255.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(15) - 701277133500.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 12*gl%z3_PCSAFT(4)**4 - 1087103750760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) - 181940450340.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) &
            & - 136455337755.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5)**2 - 11010130041.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(15) + 933775773480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 11*gl%z3_PCSAFT(4)**4 + 1141586774280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) + 141273775512.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) &
            & + 105955331634.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5)**2 + 5200854984.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(15) - 1037055940740.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)** &
            & 10*gl%z3_PCSAFT(4)**4 - 970878631032.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) - 81898852344.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - &
            & 61424139258.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5)**2 - 1004149350.d0*gl%mmean_PCSAFT** &
            & 2*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(15) + 951731562960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9* &
            & gl%z3_PCSAFT(4)**4 + 651401537760.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) + 31237092480.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + &
            & 23427819360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5)**2 - 591264840.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(15) - 712154731680.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4) &
            & **4 - 328983482880.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 4087189920.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 3065392440.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5)**2 + 535070244.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**8* &
            & gl%z3_PCSAFT(15) + 426292191360.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**4 + &
            & 112804207200.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 3177610272.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 2383207704.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5)**2 - 145643496.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**7* &
            & gl%z3_PCSAFT(15) - 198508224240.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**4 - &
            & 17747906976.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 2091690912.d0 &
            & *gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 1568768184.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5)**2 - 18873696.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(15) &
            & + 68845453344.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**4 - 4527928512.d0*gl%mmean_PCSAFT** &
            & 2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 434375424.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 325781568.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5)**2 + &
            & 20681136.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(15) - 16485539040.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**4 + 3337804800.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2 &
            & *gl%z3_PCSAFT(5) - 50310720.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - &
            & 37733040.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5)**2 - 2743200.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(15) + 2314529280.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**4 &
            & - 748846080.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 36172800.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 27129600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1) &
            & **3*gl%z3_PCSAFT(5)**2 - 811584.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(15) - &
            & 98904960.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**4 + 45826560.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 2850048.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 2137536.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)**2 + &
            & 162432.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(15) - 11289600.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)**4 + 4686336.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 451584.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 338688.d0*gl%mmean_PCSAFT**2* &
            & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)**2 + 13824.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(15) - 71424.d0* &
            & gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**4 + 82944.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 36864.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 27648.d0*gl%mmean_PCSAFT**2*gl%z3_PCSAFT(5)**2 &
            & + 5.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(15) - 60.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 45.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5)**2 - 95.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 26*gl%z3_PCSAFT(15) + 360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 720.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25* &
            & gl%z3_PCSAFT(5)**2 + 226.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(15) - 300.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**4 - 1800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) + 9560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 7170.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5)**2 + 12740.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(15) - 600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**4 - 124200.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 319520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 239640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5)**2 - &
            & 218271.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(15) + 153540.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(4)**4 + 2714760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 4003748.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 3002811.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5)**2 + 2000115.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(15) - &
            & 2770080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**4 - 29580336.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 31768320.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 23826240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5)**2 - 12636152.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(15) + 27361344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4) &
            & **4 + 214557840.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 181996352.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 136497264.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5)**2 + 59904284.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(15) - &
            & 185019360.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**4 - 1148114880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **19*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 799143520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 599357640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5)**2 - 222011305.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(15) + 940262520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18* &
            & gl%z3_PCSAFT(4)**4 + 4780051200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 2781088940.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 2085816705.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5)**2 + 657907219.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**18* &
            & gl%z3_PCSAFT(15) - 3780353760.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**4 - &
            & 15988195800.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 7828578352.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 5871433764.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 17*gl%z3_PCSAFT(5)**2 - 1578384894.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(15) + &
            & 12405215460.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**4 + 43849272024.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 18048296408.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 13536222306.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5)**2 + &
            & 3082386636.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(15) - 33868011576.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**4 - 99889469928.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2 &
            & *gl%z3_PCSAFT(5) - 34310003936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - &
            & 25732502952.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5)**2 - 4897630037.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(15) + 77795442900.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**4 &
            & + 190387355400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 53889215340.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 40416911505.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(5)**2 + 6288286805.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(15) - 151201608960.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**4 - &
            & 304425921600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 69733634720.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 52300226040.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **13*gl%z3_PCSAFT(5)**2 - 6422334508.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(15) + &
            & 249075491640.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**4 + 407753039520.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 73695756944.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 55271817708.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5)**2 + &
            & 5054931204.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(15) - 347300756880.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**4 - 454865473008.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)** &
            & 2*gl%z3_PCSAFT(5) - 62503614368.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - &
            & 46877710776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(5)**2 - 2858478584.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(15) + 408319680432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)** &
            & 4 + 418073991408.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 41137547584.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 30853160688.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5)**2 + 931324804.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**10* &
            & gl%z3_PCSAFT(15) - 402194923200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**4 - &
            & 310976814240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 19533084800.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 14649813600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9 &
            & *gl%z3_PCSAFT(5)**2 + 66912720.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(15) + &
            & 328830663120.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**4 + 181613560320.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 5320513600.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 3990385200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5)**2 - &
            & 254704672.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(15) - 220258643520.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1) &
            & **7*gl%z3_PCSAFT(4)**4 - 78593359680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) + 395095936.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + &
            & 296321952.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5)**2 + 124570304.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 7*gl%z3_PCSAFT(15) + 118671262080.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**4 + &
            & 21746752128.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 1107322368.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 830491776.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6* &
            & gl%z3_PCSAFT(5)**2 - 16068432.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(15) - 50087747712.d0 &
            & *gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**4 - 1443411072.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 469014272.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 351760704.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5)**2 - 9788288.d0* &
            & gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(15) + 15912730560.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4* &
            & gl%z3_PCSAFT(4)**4 - 1727297280.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 64807680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 48605760.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5)**2 + 4450240.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(15) - &
            & 3563527680.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**4 + 819578880.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)** &
            & 3*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 17331200.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 12998400.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5)**2 - 271168.d0*gl%mmean_PCSAFT &
            & *gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(15) + 496128000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**4 - &
            & 165888000.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 7905024.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 5928768.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) &
            & **2 - 192384.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(15) - 30612480.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)**4 + 12436992.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 801792.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 601344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(5)**2 + 29952.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(15) - 415488.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(4)**4 + 340992.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 33792.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 25344.d0*gl%mmean_PCSAFT*gl%z3_PCSAFT(5)**2 + 1536.d0*gl%mmean_PCSAFT* &
            & gl%z3_PCSAFT(15) - gl%z3_PCSAFT(1)**27*gl%z3_PCSAFT(15) + 12.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 9.d0*gl%z3_PCSAFT(1)**26*gl%z3_PCSAFT(5)**2 + 25.d0*gl%z3_PCSAFT(1)**26* &
            & gl%z3_PCSAFT(15) - 72.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 240.d0*gl%z3_PCSAFT(1)** &
            & 25*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 180.d0*gl%z3_PCSAFT(1)**25*gl%z3_PCSAFT(5)**2 - 251.d0*gl%z3_PCSAFT(1) &
            & **25*gl%z3_PCSAFT(15) + 60.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)**4 + 1080.d0*gl%z3_PCSAFT(1)**24* &
            & gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 1220.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 915.d0 &
            & *gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(5)**2 + 825.d0*gl%z3_PCSAFT(1)**24*gl%z3_PCSAFT(15) - 600.d0* &
            & gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**4 + 2520.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 15160.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 11370.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(5) &
            & **2 + 8606.d0*gl%z3_PCSAFT(1)**23*gl%z3_PCSAFT(15) - 9540.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**4 - &
            & 213480.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 320128.d0*gl%z3_PCSAFT(1)**22* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 240096.d0*gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(5)**2 - 141640.d0* &
            & gl%z3_PCSAFT(1)**22*gl%z3_PCSAFT(15) + 258240.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**4 + 2833056.d0* &
            & gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 2990080.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 2242560.d0*gl%z3_PCSAFT(1)**21*gl%z3_PCSAFT(5)**2 + 1093886.d0*gl%z3_PCSAFT(1)**21 &
            & *gl%z3_PCSAFT(15) - 2845464.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)**4 - 22374720.d0*gl%z3_PCSAFT(1)** &
            & 20*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 18658376.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 13993782.d0*gl%z3_PCSAFT(1)**20*gl%z3_PCSAFT(5)**2 - 5788246.d0*gl%z3_PCSAFT(1)** &
            & 20*gl%z3_PCSAFT(15) + 20282640.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**4 + 125970480.d0* &
            & gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 86737680.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 65053260.d0*gl%z3_PCSAFT(1)**19*gl%z3_PCSAFT(5)**2 + 23126915.d0*gl%z3_PCSAFT(1)** &
            & 19*gl%z3_PCSAFT(15) - 106575240.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**4 - 544253040.d0* &
            & gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 315670420.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 236752815.d0*gl%z3_PCSAFT(1)**18*gl%z3_PCSAFT(5)**2 - 72760671.d0*gl%z3_PCSAFT(1) &
            & **18*gl%z3_PCSAFT(15) + 439704000.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**4 + 1878319800.d0* &
            & gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 924545328.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 693408996.d0*gl%z3_PCSAFT(1)**17*gl%z3_PCSAFT(5)**2 + 184152457.d0*gl%z3_PCSAFT(1) &
            & **17*gl%z3_PCSAFT(15) - 1477157940.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**4 - 5306745096.d0 &
            & *gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 2215375884.d0*gl%z3_PCSAFT(1)**16* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 1661531913.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(5)**2 - &
            & 378996827.d0*gl%z3_PCSAFT(1)**16*gl%z3_PCSAFT(15) + 4130732904.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4) &
            & **4 + 12465020664.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 4384842872.d0* &
            & gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 3288632154.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(5) &
            & **2 + 636539640.d0*gl%z3_PCSAFT(1)**15*gl%z3_PCSAFT(15) - 9742040340.d0*gl%z3_PCSAFT(1)**14* &
            & gl%z3_PCSAFT(4)**4 - 24569137800.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - &
            & 7200642680.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 5400482010.d0*gl%z3_PCSAFT(1)**14 &
            & *gl%z3_PCSAFT(5)**2 - 870030730.d0*gl%z3_PCSAFT(1)**14*gl%z3_PCSAFT(15) + 19515883680.d0* &
            & gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**4 + 40820197680.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) + 9808501760.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + &
            & 7356376320.d0*gl%z3_PCSAFT(1)**13*gl%z3_PCSAFT(5)**2 + 957815972.d0*gl%z3_PCSAFT(1)**13* &
            & gl%z3_PCSAFT(15) - 33307509840.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**4 - 57193320240.d0* &
            & gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 11023976816.d0*gl%z3_PCSAFT(1)**12* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 8267982612.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(5)**2 - &
            & 830741176.d0*gl%z3_PCSAFT(1)**12*gl%z3_PCSAFT(15) + 48430636320.d0*gl%z3_PCSAFT(1)**11* &
            & gl%z3_PCSAFT(4)**4 + 67349302272.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 10099576192.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 7574682144.d0*gl%z3_PCSAFT(1)** &
            & 11*gl%z3_PCSAFT(5)**2 + 541743704.d0*gl%z3_PCSAFT(1)**11*gl%z3_PCSAFT(15) - 59856655008.d0* &
            & gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**4 - 66155336352.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) - 7367480064.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - &
            & 5525610048.d0*gl%z3_PCSAFT(1)**10*gl%z3_PCSAFT(5)**2 - 235258168.d0*gl%z3_PCSAFT(1)**10* &
            & gl%z3_PCSAFT(15) + 62601407520.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**4 + 53510379840.d0* &
            & gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 4081681280.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 3061260960.d0*gl%z3_PCSAFT(1)**9*gl%z3_PCSAFT(5)**2 + 35468160.d0*gl%z3_PCSAFT(1) &
            & **9*gl%z3_PCSAFT(15) - 55032395040.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**4 - 34890073920.d0 &
            & *gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 1524170880.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 1143128160.d0*gl%z3_PCSAFT(1)**8*gl%z3_PCSAFT(5)**2 + 33945760.d0*gl%z3_PCSAFT(1) &
            & **8*gl%z3_PCSAFT(15) + 40282467840.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**4 + 17660822400.d0 &
            & *gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 207882240.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 155911680.d0*gl%z3_PCSAFT(1)**7*gl%z3_PCSAFT(5)**2 - 28498784.d0*gl%z3_PCSAFT(1)** &
            & 7*gl%z3_PCSAFT(15) - 24233623680.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**4 - 6407700480.d0* &
            & gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 155829248.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 116871936.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(5)**2 + 8924064.d0*gl%z3_PCSAFT(1)**6 &
            & *gl%z3_PCSAFT(15) + 11766783360.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**4 + 1282033152.d0* &
            & gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 118275584.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) - 88706688.d0*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(5)**2 + 232640.d0*gl%z3_PCSAFT(1)**5* &
            & gl%z3_PCSAFT(15) - 4493738880.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**4 + 136742400.d0*gl%z3_PCSAFT(1) &
            & **4*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 35475200.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(9) + 26606400.d0*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(5)**2 - 1155520.d0*gl%z3_PCSAFT(1)**4* &
            & gl%z3_PCSAFT(15) + 1298611200.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**4 - 202314240.d0*gl%z3_PCSAFT(1) &
            & **3*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) - 2040320.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) &
            & - 1530240.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5)**2 + 334208.d0*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(15) - &
            & 266760960.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**4 + 71723520.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) - 2061824.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 1546368.d0* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)**2 - 256.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(15) + 34690560.d0* &
            & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**4 - 12874752.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5) + &
            & 663552.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 497664.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)**2 - &
            & 17152.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(15) - 2145792.d0*gl%z3_PCSAFT(4)**4 + 995328.d0*gl%z3_PCSAFT(4) &
            & **2*gl%z3_PCSAFT(5) - 68608.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) - 51456.d0*gl%z3_PCSAFT(5)**2 + &
            & 2560.d0*gl%z3_PCSAFT(15))/recurring_factor**5
    end if

!DEC$ END IF
end subroutine CDERIVS
    
subroutine MEODERIVS(gl,T, GETDERMEO, index)
    
    ! meo1, meo2: abbreviations
    ! defined by eq. A.12 and A.13 in Gross, Sadowski 2001:
    ! meo1 = Sum_i(Sum_j(x_i*x_j*m_i*m_j*(epsilon_ij/(kB*T))*sigma_ij**3))     if index==1
    ! meo2 = Sum_i(Sum_j(x_i*x_j*m_i*m_j*(epsilon_ij/(kB*T))**2*sigma_ij**3))  if index==2
    ! combining rules for sigma_ij and epsilon_ij
    ! sigma_ij = 1/2*(sigma_i+sigma_j)
    ! epsilon_ij = sqrt(epsilon_i*epsilon_j)*(1-kij)
    ! dependent only on T






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    double precision, intent (in) :: T
    integer, dimension (nderivs), intent (in) :: GETDERMEO
    integer, intent (in) :: index
    !output: meo1_PCSAFT, meo2_PCSAFT (module variables)
    !working variables
    double precision :: sigma_ij, epsk_ij, basiceqn, derivative
    integer :: i, j

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! Basic Equation, all T derivatives are just this equation times a factor
    basiceqn = 0.d0
    do i = gl%n_start, gl%n_end
        do j = gl%n_start, gl%n_end
            sigma_ij = 1.d0/2.d0*(gl%sigPCSAFT(i)+gl%sigPCSAFT(j))
            epsk_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j)) !TE: kij = 0 for PCP-SAFT Gross 2005
            basiceqn = basiceqn + gl%molfractions(i)*gl%molfractions(j)*gl%mPCSAFT(i)*gl%mPCSAFT(j) &
                & * (epsk_ij/T)**index * sigma_ij**3
        end do
    end do

    !calculate the derivatives of meo1, meo2
    ! 1: meo1 or meo2
    if (GETDERMEO(1) .eq. 1) then
        derivative = basiceqn
        select case (index)
            case (1)
                gl%meo1_PCSAFT(1) = derivative
            case (2)
                gl%meo2_PCSAFT(1) = derivative
        end select
    end if
    
    ! 2, 3: derivatives with respect to D equal zero
    
    ! 4: 1ST DERIVATIVE OF meo1, meo2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! T*d(meo_index)/dT = -index*meo_index
    if (GETDERMEO(4) .eq. 1) then
        derivative = -index*basiceqn
        select case (index)
            case (1)
                gl%meo1_PCSAFT(4) = derivative
            case (2)
                gl%meo2_PCSAFT(4) = derivative
        end select
    end if
    
    ! 5: 2ND DERIVATIVE OF meo1, meo2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2 
    if (GETDERMEO(5) .eq. 1) then
        derivative = index*(index+1)*basiceqn
        select case (index)
            case (1)
                gl%meo1_PCSAFT(5) = derivative
            case (2)
                gl%meo2_PCSAFT(5) = derivative
        end select
    end if
    
    ! 6, 7, 8: derivatives with respect to D equal zero
    
    ! 9: 3RD DERIVATIVE OF meo1, meo2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^3 
    if (GETDERMEO(9) .eq. 1) then
        derivative = -index*(index+1)*(index+2)*basiceqn
        select case (index)
            case (1)
                gl%meo1_PCSAFT(9) = derivative
            case (2)
                gl%meo2_PCSAFT(9) = derivative
        end select
    end if
    
    ! 10, 11, 12, 13, 14: derivatives with respect to D equal zero
    
    ! 15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
    if (GETDERMEO(15) .eq. 1) then
        derivative = index*(index+1)*(index+2)*(index+3)*basiceqn
        select case (index)
            case (1)
                gl%meo1_PCSAFT(15) = derivative
            case (2)
                gl%meo2_PCSAFT(15) = derivative
        end select
    end if
    
!DEC$ END IF
end subroutine MEODERIVS

subroutine IDERIVS(gl,GETDERI, index)

! Henning Markgraf, June 2016
    
    ! I_1, I_2: integrals of the perturbation theory
    ! defined by eq. A.14 and A.15 in Gross, Sadowski 2001:
    ! I_1 = Sum_0_6(a_i*eta**i)     if index==1
    ! I_2 = Sum_0_6(b_i*eta**i)     if index==2
    ! dependent on D and T
    ! derivative stored in module variables





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    integer, intent (in) :: index
    !output: i1_PCSAFT, i2_PCSAFT (module variables)
    integer :: i
    double precision :: sum
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !calculate the derivatives of I_1, I_2
    ! 1: I_index
    if (GETDERI(1) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**i
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(1) = sum
            case (2)
                gl%i2_PCSAFT(1) = sum
        end select
    end if
    
    !  2: 1ST DERIVATIVE OF I_1, I_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + i*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**i
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(2) = sum
            case (2)
                gl%i2_PCSAFT(2) = sum
        end select
    end if
    
    ! 3: 2ND DERIVATIVE OF I_1, I_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + i*(i-1)*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**i
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(3) = sum
            case (2)
                gl%i2_PCSAFT(3) = sum
        end select
    end if
    
    ! 4: 1ST DERIVATIVE OF I_1, I_2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(1)**(i-1)*gl%z3_PCSAFT(4)
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(4) = sum
            case (2)
                gl%i2_PCSAFT(4) = sum
        end select
    end if
    
    ! 5: 2ND DERIVATIVE OF I_1, I_2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + gl%ab_PCSAFT(index,i)*i*gl%z3_PCSAFT(1)**(i-2)*(i*gl%z3_PCSAFT(4)**2+gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)-gl%z3_PCSAFT(4)**2)
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(5) = sum
            case (2)
                gl%i2_PCSAFT(5) = sum
        end select
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF I_1, I_2 WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    if (GETDERI(6) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + gl%ab_PCSAFT(index,i)*i**2*gl%z3_PCSAFT(1)**(i-1)*gl%z3_PCSAFT(4)
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(6) = sum
            case (2)
                gl%i2_PCSAFT(6) = sum
        end select
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF I_1, I_2 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! 7 equals 5, but with i**2 instead of i in the sum
    ! requires zeta(5,4,1)  (number derivatives)
    if (GETDERI(7) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + gl%ab_PCSAFT(index,i)*i**2*gl%z3_PCSAFT(1)**(i-2)*(i*gl%z3_PCSAFT(4)**2+gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)-gl%z3_PCSAFT(4)**2)
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(7) = sum
            case (2)
                gl%i2_PCSAFT(7) = sum
        end select
    end if
    
    ! 8: 3RD DERIVATIVE OF I_1, I_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    if (GETDERI(8) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + i*(i-1)*(i-2)*gl%ab_PCSAFT(index,i)*gl%z3_PCSAFT(1)**i
        end do 
        select case (index)
            case (1)
                gl%i1_PCSAFT(8) = sum
            case (2)
                gl%i2_PCSAFT(8) = sum
        end select
    end if
    
    ! 9: 3RD DERIVATIVE OF I_1, I_2 WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires zeta(9,5,4,1)  (number derivatives)
    if (GETDERI(9) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + i*gl%z3_PCSAFT(1)**(i-3)*(i**2*gl%z3_PCSAFT(4)**3 + 3.d0*i*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
                & - 3.d0*i*gl%z3_PCSAFT(4)**3 + gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) - 3.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
                & gl%z3_PCSAFT(5) + 2.d0*gl%z3_PCSAFT(4)**3)*gl%ab_PCSAFT(index,i)
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(9) = sum
            case (2)
                gl%i2_PCSAFT(9) = sum
        end select
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF I_1, I_2 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires zeta(4,1)  (number derivatives)
    if (GETDERI(10) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + i**2*gl%z3_PCSAFT(1)**(i - 1)*gl%z3_PCSAFT(4)*(i - 1)*gl%ab_PCSAFT(index,i)
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(10) = sum
            case (2)
                gl%i2_PCSAFT(10) = sum
        end select
    end if
    
    ! 11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^4
    ! requires zeta(1)  (number derivatives)
    if (GETDERI(11) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + gl%z3_PCSAFT(1)**i*(i**4 - 6*i**3 + 11*i**2 - 6*i)*gl%ab_PCSAFT(index,i)
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(11) = sum
            case (2)
                gl%i2_PCSAFT(11) = sum
        end select
    end if
    
    ! 12: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, DEL AND DEL, MULTIPLIED BY TAU*DEL^3
    ! requires zeta(4,1)  (number derivatives)
    if (GETDERI(12) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + gl%z3_PCSAFT(1)**(i - 1)*gl%z3_PCSAFT(4)*(i**4 - 3*i**3 + 2*i**2)*gl%ab_PCSAFT(index,i)
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(12) = sum
            case (2)
                gl%i2_PCSAFT(12) = sum
        end select
    end if
    
    ! 13: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, DEL, AND DEL, MULTIPLIED BY TAU^2*DEL^2
    ! requires zeta(5,4,1)  (number derivatives)
    if (GETDERI(13) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + i**2*(i-1)*gl%ab_PCSAFT(index,i)*(gl%z3_PCSAFT(1)**(i - 2)*gl%z3_PCSAFT(4)**2*(i-1) + gl%z3_PCSAFT(1)**(i-1)* &
                & gl%z3_PCSAFT(5))
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(13) = sum
            case (2)
                gl%i2_PCSAFT(13) = sum
        end select
    end if
    
    ! 14: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, TAU AND DEL, MULTIPLIED BY TAU^3*DEL
    ! requires zeta(9,5,4,1)  (number derivatives)
    if (GETDERI(14) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + i**2*gl%ab_PCSAFT(index,i)*(gl%z3_PCSAFT(1)**(i-3)*gl%z3_PCSAFT(4)**3*(i-1)*(i-2) + &
                & 3.d0*gl%z3_PCSAFT(1)**(i-2)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(i-1) + gl%z3_PCSAFT(1)**(i-1) &
                & *gl%z3_PCSAFT(9))
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(14) = sum
            case (2)
                gl%i2_PCSAFT(14) = sum
        end select
    end if
    
    ! 15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
    ! requires zeta(15,9,5,4,1)  (number derivatives)
    if (GETDERI(15) .eq. 1) then
        sum = 0.d0
        do i = 0, 6
            sum = sum + i*gl%ab_PCSAFT(index,i)*(gl%z3_PCSAFT(1)**(i-4)*gl%z3_PCSAFT(4)**4*(i-1)*(i-2)*(i-3) &
                & + 6.d0*gl%z3_PCSAFT(1)**(i-3)*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5)*(i-1)*(i-2) + 4.d0*gl%z3_PCSAFT(1)**(i-2)*gl%z3_PCSAFT(4)* &
                & gl%z3_PCSAFT(9)*(i-1) + 3.d0*gl%z3_PCSAFT(1)**(i-2)*gl%z3_PCSAFT(5)**2*(i-1) + gl%z3_PCSAFT(1)**(i-1)* &
                & gl%z3_PCSAFT(15))
        end do
        select case (index)
            case (1)
                gl%i1_PCSAFT(15) = sum
            case (2)
                gl%i2_PCSAFT(15) = sum
        end select
    end if
    
    !DEC$ END IF
end subroutine IDERIVS
 
subroutine AHCDERIVS(gl,GETDERAHC)
    
    ! a_hc: Hard chain reference contribution to the Helmholtz free energy
    ! defined by eq. A.4 in Gross, Sadowski 2001:
    ! a_hc = mmean*a_hs - Sum( x_i*(m_i-1)*ln(g_ii_hs) )
    ! dependent on D and T






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAHC
    !output: ahc_PCSAFT (module variable)
    !working variables
    double precision :: sum
    integer :: i
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !Initializations
    sum = 0.d0
    
    !calculate the derivatives of a_hc
    ! 1: a_hc
    ! requires g_ii_hs, ahs_PCSAFT
    if (GETDERAHC(1) .eq. 1) then
        !iterate over all components of the mixture to calculate sum x_i*(m_i-1)*ln(g_ii_hs)
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + gl%molfractions(i)*(gl%mPCSAFT(i)-1.d0)*dlog(gl%gii_PCSAFT(i,1))
        end do
        !calculate a_hc
        gl%ahc_PCSAFT(1) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(1)-sum
    end if
    
    !  2: 1ST DERIVATIVE OF a_hc WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires g_ii_hs, g_ii_hs_D, a_hs_D
    if (GETDERAHC(2) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + gl%molfractions(i)*(gl%mPCSAFT(i)-1.d0)/gl%gii_PCSAFT(i,1)*gl%gii_PCSAFT(i,2)
        end do
        ! factor D is already included in ahs_PCSAFT(2) and gii_PCSAFT(i,2)
        gl%ahc_PCSAFT(2) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(2)-sum
    end if
    
    ! 3: 2ND DERIVATIVE OF a_hc WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires g_ii_hs, g_ii_hs_D, g_ii_hs_DD, a_hs_DD
    if (GETDERAHC(3) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + (gl%gii_PCSAFT(i,3) - gl%gii_PCSAFT(i,2)**2/gl%gii_PCSAFT(i,1))*(gl%mPCSAFT(i)-1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)
        end do
        ! factor D^2 is already included in ahs_PCSAFT(3) and gii_PCSAFT(i,3) and gii_PCSAFT(i,2)**2
        gl%ahc_PCSAFT(3) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(3)-sum
    end if
    
    ! 4: 1ST DERIVATIVE OF a_hc WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires g_ii_hs, g_ii_hs_T, a_hs_T
    if (GETDERAHC(4) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + gl%gii_PCSAFT(i,4)*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)
        end do
        ! factor T is already included in ahs_PCSAFT(4) and gii_PCSAFT(i,4)
        gl%ahc_PCSAFT(4) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(4)-sum
    end if
    
    ! 5: 2ND DERIVATIVE OF a_hc WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires g_ii_hs, g_ii_hs_T, g_ii_hs_TT, a_hs_TT
    if (GETDERAHC(5) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + (gl%gii_PCSAFT(i,5) - gl%gii_PCSAFT(i,4)**2/gl%gii_PCSAFT(i,1))*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)
        end do
        ! factor T^2 is already included in ahs_PCSAFT(5) and gii_PCSAFT(i,5)
        gl%ahc_PCSAFT(5) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(5)-sum
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_hc WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires g_ii_hs, g_ii_hs_T, g_ii_hs_D, g_ii_hs_TD, a_hs_TD
    if (GETDERAHC(6) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + (gl%gii_PCSAFT(i,6)*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1) - &
                & gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,2)*(gl%mPCSAFT(i) - 1.d0)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**2)
        end do
        ! factor *T*D is already included in ahs_PCSAFT(6) and gii_PCSAFT(i,6)
        gl%ahc_PCSAFT(6) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(6)-sum
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_hc WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires ahs_PCSAFT(7) and g_ii_hs(7,6,5,4,2,1)  (derivative numbers)
    if (GETDERAHC(7) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + gl%molfractions(i)/gl%gii_PCSAFT(i,1)**3*(gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,1)**2*gl%gii_PCSAFT(i,7) - &
                    & 2.d0*gl%gii_PCSAFT(i,1)*gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,6) + gl%gii_PCSAFT(i,4)**2*gl%gii_PCSAFT(i,2) - &
                    & gl%gii_PCSAFT(i,2)*(gl%gii_PCSAFT(i,1)*gl%gii_PCSAFT(i,5) - gl%gii_PCSAFT(i,4)**2))
        end do
        ! factor *T*T*D is already included in the GII and AHS derivatives
        gl%ahc_PCSAFT(7) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(7)-sum
    end if
    
    ! 8: 3RD DERIVATIVE OF a_hc WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires g_ii_hs, g_ii_hs_D, g_ii_hs_DD, g_ii_hs_DDD, a_hs_DDD
    if (GETDERAHC(8) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + (gl%mPCSAFT(i)-1.d0)*(gl%gii_PCSAFT(i,8) - 3.d0*gl%gii_PCSAFT(i,2)*gl%gii_PCSAFT(i,3)/gl%gii_PCSAFT(i,1) + &
                & 2.d0*gl%gii_PCSAFT(i,2)**3/gl%gii_PCSAFT(i,1)**2)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)
        end do
        ! factor D^3 is already included in ahs_PCSAFT(8) and gii_PCSAFT(i,8) and gii_PCSAFT(i,2)*gii_PCSAFT(i,3)
        gl%ahc_PCSAFT(8) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(8)-sum
    end if
    
    ! 9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires ahs_PCSAFT(9) and g_ii_hs(9,5,4,1)  (derivative numbers)
    if (GETDERAHC(9) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + gl%molfractions(i)/gl%gii_PCSAFT(i,1)*(gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,9) - 3.d0*gl%gii_PCSAFT(i,4)* &
                & gl%gii_PCSAFT(i,5)/gl%gii_PCSAFT(i,1) + 2.d0*gl%gii_PCSAFT(i,4)**3/gl%gii_PCSAFT(i,1)**2)
        end do
        ! factor T^3 is already included in ahs_PCSAFT and DERGII
        gl%ahc_PCSAFT(9) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(9)-sum
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires ahs_PCSAFT(10) and g_ii_hs(10,6,4,3,2,1)  (derivative numbers)
    if (GETDERAHC(10) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + gl%molfractions(i)/gl%gii_PCSAFT(i,1)**3*(gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,1)**2* &
                & gl%gii_PCSAFT(i,10) - gl%gii_PCSAFT(i,1)*(gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,3) + 2.d0* &
                & gl%gii_PCSAFT(i,6)*gl%gii_PCSAFT(i,2)) + 2.d0*gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,2)**2)
        end do
        ! factor T^3 is already included in ahs_PCSAFT and DERGII
        gl%ahc_PCSAFT(10) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(10)-sum
    end if
    
    ! 11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^4
    ! requires ahs_PCSAFT(11) and g_ii_hs(11,8,3,2,1)  (derivative numbers)
    if (GETDERAHC(11) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,11) - 4.d0 &
                & *gl%gii_PCSAFT(i,2)*gl%gii_PCSAFT(i,8)/gl%gii_PCSAFT(i,1) - 3.d0*gl%gii_PCSAFT(i,3)**2/ &
                & gl%gii_PCSAFT(i,1) + 12.d0*gl%gii_PCSAFT(i,2)**2*gl%gii_PCSAFT(i,3)/gl%gii_PCSAFT(i,1)**2 - 6.d0* &
                & gl%gii_PCSAFT(i,2)**4/gl%gii_PCSAFT(i,1)**3)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)
        end do
        gl%ahc_PCSAFT(11) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(11)-sum
    end if
    
    ! 12: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, DEL AND DEL, MULTIPLIED BY TAU*DEL^3
    ! requires ahs_PCSAFT(12) and g_ii_hs(12,10,8,6,4,3,2,1)  (derivative numbers)
    if (GETDERAHC(12) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,12) - &
                & gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,8)/gl%gii_PCSAFT(i,1) - 3.d0*gl%gii_PCSAFT(i,6)* &
                & gl%gii_PCSAFT(i,3)/gl%gii_PCSAFT(i,1) - 3.d0*gl%gii_PCSAFT(i,10)*gl%gii_PCSAFT(i,2)/gl%gii_PCSAFT(i,1) + &
                & 6.d0*gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,2)*gl%gii_PCSAFT(i,3)/gl%gii_PCSAFT(i,1)**2 + 6.d0* &
                & gl%gii_PCSAFT(i,6)*gl%gii_PCSAFT(i,2)**2/gl%gii_PCSAFT(i,1)**2 - 6.d0*gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,2) &
                & **3/gl%gii_PCSAFT(i,1)**3)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)
        end do
        gl%ahc_PCSAFT(12) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(12)-sum
    end if
    
    ! 13: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, DEL, AND DEL, MULTIPLIED BY TAU^2*DEL^2
    ! requires ahs_PCSAFT(13) and g_ii_hs(13,10,7,6,5,4,3,2,1)  (derivative numbers)
    if (GETDERAHC(13) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,13) - 2.d0* &
                & gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,10)/gl%gii_PCSAFT(i,1) - 2.d0*gl%gii_PCSAFT(i,6)**2/gl%gii_PCSAFT(i,1) - 2.d0 &
                & *gl%gii_PCSAFT(i,2)*(gl%gii_PCSAFT(i,7) - 2.d0*gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,6)/gl%gii_PCSAFT(i,1) + &
                & gl%gii_PCSAFT(i,4)**2*gl%gii_PCSAFT(i,2)/gl%gii_PCSAFT(i,1)**2)/gl%gii_PCSAFT(i,1) - gl%gii_PCSAFT(i,3)*( &
                & gl%gii_PCSAFT(i,5) - gl%gii_PCSAFT(i,4)**2/gl%gii_PCSAFT(i,1))/gl%gii_PCSAFT(i,1) + gl%gii_PCSAFT(i,4)**2* &
                & gl%gii_PCSAFT(i,3)/gl%gii_PCSAFT(i,1)**2 + 4.d0*gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,6)*gl%gii_PCSAFT(i,2)/ &
                & gl%gii_PCSAFT(i,1)**2 + 2.d0*gl%gii_PCSAFT(i,2)**2*(gl%gii_PCSAFT(i,5) - gl%gii_PCSAFT(i,4)**2/gl%gii_PCSAFT(i,1)) &
                & /gl%gii_PCSAFT(i,1)**2 - 2.d0*gl%gii_PCSAFT(i,4)**2*gl%gii_PCSAFT(i,2)**2/gl%gii_PCSAFT(i,1)**3)*gl%molfractions(i)/ &
                & gl%gii_PCSAFT(i,1)
        end do
        gl%ahc_PCSAFT(13) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(13)-sum
    end if
    
    ! 14: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, TAU AND DEL, MULTIPLIED BY TAU^3*DEL
    ! requires ahs_PCSAFT(14) and g_ii_hs(14,9,7,6,5,4,2,1)  (derivative numbers)
    if (GETDERAHC(14) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,14) - 3.d0*gl%gii_PCSAFT(i,4)* &
                & gl%gii_PCSAFT(i,7)/gl%gii_PCSAFT(i,1) - 3.d0*gl%gii_PCSAFT(i,5)*gl%gii_PCSAFT(i,6)/gl%gii_PCSAFT(i,1) + 6.d0* &
                & gl%gii_PCSAFT(i,4)**2*gl%gii_PCSAFT(i,6)/gl%gii_PCSAFT(i,1)**2 + 3.d0*gl%gii_PCSAFT(i,4)*gl%gii_PCSAFT(i,5)* &
                & gl%gii_PCSAFT(i,2)/gl%gii_PCSAFT(i,1)**2 - 4.d0*gl%gii_PCSAFT(i,4)**3*gl%gii_PCSAFT(i,2)/gl%gii_PCSAFT(i,1)**3)*gl%molfractions &
                & (i)/gl%gii_PCSAFT(i,1) - gl%gii_PCSAFT(i,2)*(gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,9) - 3.d0*gl%gii_PCSAFT(i,4)* &
                & gl%gii_PCSAFT(i,5)/gl%gii_PCSAFT(i,1) + 2.d0*gl%gii_PCSAFT(i,4)**3/gl%gii_PCSAFT(i,1)**2)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)**2
        end do
        gl%ahc_PCSAFT(14) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(14)-sum
    end if
    
    ! 15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
    ! requires ahs_PCSAFT(15) and g_ii_hs(15,9,5,4,1)  (derivative numbers)
    if (GETDERAHC(15) .eq. 1) then
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + (gl%mPCSAFT(i) - 1.d0)*(gl%gii_PCSAFT(i,15) - 4.d0*gl%gii_PCSAFT(i,4)* &
                & gl%gii_PCSAFT(i,9)/gl%gii_PCSAFT(i,1) - 3.d0*gl%gii_PCSAFT(i,5)**2/gl%gii_PCSAFT(i,1) + 12.d0*gl%gii_PCSAFT(i,4)**2* &
                & gl%gii_PCSAFT(i,5)/gl%gii_PCSAFT(i,1)**2 - 6.d0*gl%gii_PCSAFT(i,4)**4/gl%gii_PCSAFT(i,1)**3)*gl%molfractions(i)/gl%gii_PCSAFT(i,1)
        end do
        gl%ahc_PCSAFT(15) = gl%mmean_PCSAFT*gl%ahs_PCSAFT(15)-sum
    end if
    
!DEC$ END IF
end subroutine AHCDERIVS
    
subroutine AHSDERIVS(gl,GETDERAHS)

! Henning Markgraf, June 2016
    
    ! a_hs: Helmholtz free energy of the hard sphere fluid
    ! defined by eq. A.6 in Gross, Sadowski 2001:
    ! a_hs = 1/zeta_0 * ( 3*zeta_1*zeta_2/(1-zeta_3) + zeta_2**3/(zeta_3*(1-zeta_3)**2)
    !           + (zeta_2**3/zeta_3**2 - zeta_0)*ln(1-zeta_3) )
    ! dependent on D and T





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAHS
    !output: ahs_PCSAFT (module variable)
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !calculate the derivatives of a_hs
    ! 1: a_hs
    if (GETDERAHS(1) .eq. 1) then
        !calculate a_hs
        gl%ahs_PCSAFT(1) = 1.d0/gl%z0_PCSAFT(1) *(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(1.d0-gl%z3_PCSAFT(1)) + gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*(1.d0-gl%z3_PCSAFT(1))**2) &
              & + (gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2 - gl%z0_PCSAFT(1))*dlog(1.d0-gl%z3_PCSAFT(1)))
    end if
    
    ! 2: 1ST DERIVATIVE OF a_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, d_i, d_j
    if (GETDERAHS(2) .eq. 1) then
        ! relation zeta_n_D*D = zeta_n was used to simplify the derivative
        gl%ahs_PCSAFT(2) = (3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0)**2 + 3.d0*gl%z1_PCSAFT(1)* &
     & gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0) + 2.d0* &
     & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0)**3) - gl%z2_PCSAFT(1)**3* &
     & gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1)**2*(-gl%z3_PCSAFT(1) + 1.d0)**2) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)/ &
     & (gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0)**2) - gl%z3_PCSAFT(1)*(-gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/ &
     & gl%z3_PCSAFT(1)**2)/(-gl%z3_PCSAFT(1) + 1.d0) + (-gl%z0_PCSAFT(1) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/ &
     & gl%z3_PCSAFT(1)**3 + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0 &
     & ))/gl%z0_PCSAFT(1) - gl%z0_PCSAFT(1)*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0) + gl%z2_PCSAFT(1)**3 &
     & /(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0)**2) + (-gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog( &
     & -gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1)**2
    end if
    
    ! 3: 2ND DERIVATIVE OF a_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT
    if (GETDERAHS(3) .eq. 1) then
        ! relation zeta_n_D*D = zeta_n was used to simplify the derivative
        gl%ahs_PCSAFT(3) = (-6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1) - 1.d0)**3 + 6.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**2 - 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0 &
            & ) + 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4) + 4.d0*gl%z2_PCSAFT(1) &
            & **3*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3) + 2.d0*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2) - 12.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3) - 6.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) + 6.d0*gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)**2/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) + gl%z3_PCSAFT(1)**2*(gl%z0_PCSAFT(1) - &
            & gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0)**2 - 2.d0*gl%z3_PCSAFT(1)*(gl%z0_PCSAFT(1) + &
            & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)/gl%z3_PCSAFT(1) &
            & **2)/(gl%z3_PCSAFT(1) - 1.d0) - (-6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2/gl%z3_PCSAFT(1)**4 + 12.d0* &
            & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/gl%z3_PCSAFT(1)**3 - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(1)** &
            & 2/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0) + 2.d0*gl%z0_PCSAFT(1)*(-3.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - &
            & 1.d0) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/( &
            & gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) &
            & - 1.d0)**2) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) + &
            & gl%z3_PCSAFT(1)*(gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0) + ( &
            & gl%z0_PCSAFT(1) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(1)/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1) - 2.d0*gl%z0_PCSAFT(1) &
            & **2*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) - gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0 &
            & )**2) + (gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/ &
            & gl%z0_PCSAFT(1)**2)/gl%z0_PCSAFT(1)
    end if
    
    ! 4: 1ST DERIVATIVE OF a_hs WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T
    if (GETDERAHS(4) .eq. 1) then
        gl%ahs_PCSAFT(4) = (3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)/(-gl%z3_PCSAFT(1) + 1.d0)**2 + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)/( &
            & -gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0)**3) - gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)** &
            & 2*(-gl%z3_PCSAFT(1) + 1.d0)**2) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0) &
            & **2) - gl%z3_PCSAFT(4)*(-gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(-gl%z3_PCSAFT(1) + 1.d0) + ( &
            & -gl%z0_PCSAFT(4) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4) &
            & /gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1) - gl%z0_PCSAFT(4)*(3.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0) + gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0)**2) + ( &
            & -gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1)**2
    end if
    
    ! 5: 2ND DERIVATIVE OF a_hs WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T, z2_PCSAFT_TT, z3_PCSAFT_TT
    if (GETDERAHS(5) .eq. 1) then
        gl%ahs_PCSAFT(5) = (-6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0) + 6.d0*gl%z2_PCSAFT(1)**3 &
            & *gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 + 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)**2 - 2.d0* &
            & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)* &
            & gl%z3_PCSAFT(4)) - 2.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**3*(2.d0*gl%z2_PCSAFT(1) &
            & *gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + 3.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0) &
            & **2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2) + gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4* &
            & (6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0 &
            & *gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4) &
            & **2))*dlog(-gl%z3_PCSAFT(1) + 1.d0) - 3.d0*gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4) + gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)) + 3.d0*gl%z3_PCSAFT(1) &
            & **4*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4) &
            & *gl%z3_PCSAFT(4) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) - gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**3*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) &
            & - 1.d0)**2*(-gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + &
            & gl%z3_PCSAFT(4)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)))/(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**4*( &
            & gl%z3_PCSAFT(1) - 1.d0)**4)
    end if
    
    
    ! 6: 2ND MIXED DERIVATIVE OF a_hs WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T
    if (GETDERAHS(6) .eq. 1) then
        gl%ahs_PCSAFT(6) = -(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)/(-gl%z3_PCSAFT(1) + 1.d0)**2 + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)/( &
            & -gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0)**3) - gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)** &
            & 2*(-gl%z3_PCSAFT(1) + 1.d0)**2) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0) &
            & **2) - gl%z3_PCSAFT(4)*(-gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(-gl%z3_PCSAFT(1) + 1.d0) + ( &
            & -gl%z0_PCSAFT(4) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4) &
            & /gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1) + (6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)/(-gl%z3_PCSAFT(1) + 1.d0)**3 + 9.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)/( &
            & -gl%z3_PCSAFT(1) + 1.d0)**2 + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0)**2 + 6.d0* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)/(-gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(-gl%z3_PCSAFT(1) &
            & + 1.d0)**2 + 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0) + 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/ &
            & (-gl%z3_PCSAFT(1) + 1.d0)**4 + 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0)**3 &
            & ) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)**2*(-gl%z3_PCSAFT(1) + 1.d0)**2) + 6.d0*gl%z2_PCSAFT(1) &
            & **2*gl%z2_PCSAFT(4)/(-gl%z3_PCSAFT(1) + 1.d0)**3 + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/(gl%z3_PCSAFT(1)*( &
            & -gl%z3_PCSAFT(1) + 1.d0)**2) - gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(-gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2 &
            & )/(-gl%z3_PCSAFT(1) + 1.d0)**2 - gl%z3_PCSAFT(1)*(-gl%z0_PCSAFT(4) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/ &
            & gl%z3_PCSAFT(1)**3 + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/gl%z3_PCSAFT(1)**2)/(-gl%z3_PCSAFT(1) + 1.d0) - 2.d0* &
            & gl%z3_PCSAFT(4)*(-gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(-gl%z3_PCSAFT(1) + 1.d0) + ( &
            & -gl%z0_PCSAFT(4) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4) &
            & /gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1) + gl%z0_PCSAFT(4)*(3.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0) + gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0)**2) + ( &
            & -gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1)**2 - &
            & gl%z0_PCSAFT(4)*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0)**2 + 6.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3/(-gl%z3_PCSAFT(1) + 1.d0)**3 + 2.d0*gl%z2_PCSAFT(1)**3 &
            & /(gl%z3_PCSAFT(1)*(-gl%z3_PCSAFT(1) + 1.d0)**2) - gl%z3_PCSAFT(1)*(-gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)** &
            & 2)/(-gl%z3_PCSAFT(1) + 1.d0) + (-gl%z0_PCSAFT(1) + gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + &
            & 1.d0))/gl%z0_PCSAFT(1)**2
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_hs WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T, z0_PCSAFT_TT, z1_PCSAFT_TT, z2_PCSAFT_TT, z3_PCSAFT_TT
    if (GETDERAHS(7) .eq. 1) then
        gl%ahs_PCSAFT(7) = (gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 2.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)** &
            & 2 - 3.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 4.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)**2 + 3.d0*gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 2.d0*gl%z0_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)**2 - gl%z0_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
            & **2 + 12.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 18.d0*gl%z0_PCSAFT(1)**2 &
            & *gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 6.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) &
            & - 12.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 24.d0*gl%z0_PCSAFT(1)**2* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 12.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)* &
            & gl%z3_PCSAFT(4) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**3 - 9.d0*gl%z0_PCSAFT(1)**2* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)* &
            & gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5) - 12.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 24.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4) - 12.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%z0_PCSAFT(1)**2* &
            & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)* &
            & gl%z3_PCSAFT(1)**2 + 18.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) - 6.d0*gl%z0_PCSAFT(1)**2* &
            & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - 9.d0* &
            & gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 9.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(5)* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 3.d0*gl%z0_PCSAFT(1)**2*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1) - 2.d0*gl%z0_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)**2 + 10.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 30.d0* &
            & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 8.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(5) - 12.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + &
            & 60.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) - 48.d0*gl%z0_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 3.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)* &
            & gl%z3_PCSAFT(1)**3 - 15.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2 + 21.d0* &
            & gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) - 9.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(5) + 6.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)**3 - 30.d0*gl%z0_PCSAFT(1) &
            & **2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)**2 + 42.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)** &
            & 2*gl%z3_PCSAFT(1) - 18.d0*gl%z0_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2 + 12.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 24.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 12.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(4) - 6.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 18.d0* &
            & gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 18.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4) - 6.d0* &
            & gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)* &
            & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 18.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1) + 6.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1) + 4.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)* &
            & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 20.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 16.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4) - 6.d0* &
            & gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 30.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4) &
            & *gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 42.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 18.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(4)*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4) - 3.d0* &
            & gl%z0_PCSAFT(1)*gl%z0_PCSAFT(5)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 9.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(5)* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 9.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(5)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(5)*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) - gl%z0_PCSAFT(1)*gl%z0_PCSAFT(5)* &
            & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3 + 5.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(5)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2 - 7.d0 &
            & *gl%z0_PCSAFT(1)*gl%z0_PCSAFT(5)*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(1)*gl%z0_PCSAFT(5)*gl%z2_PCSAFT(1) &
            & **3 + 6.d0*gl%z0_PCSAFT(4)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 - 18.d0*gl%z0_PCSAFT(4)**2* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 18.d0*gl%z0_PCSAFT(4)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 6.d0 &
            & *gl%z0_PCSAFT(4)**2*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1) + 2.d0*gl%z0_PCSAFT(4)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3 - &
            & 10.d0*gl%z0_PCSAFT(4)**2*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2 + 14.d0*gl%z0_PCSAFT(4)**2*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(1) - 6.d0*gl%z0_PCSAFT(4)**2*gl%z2_PCSAFT(1)**3)/(gl%z0_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**5)
    end if
    
    ! 8: 3RD DERIVATIVE OF a_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT
    if (GETDERAHS(8) .eq. 1) then
        ! relation zeta_n_D*D = zeta_n was used to simplify the derivative
        gl%ahs_PCSAFT(8) = (18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3/(gl%z3_PCSAFT(1) - 1.d0)**4 - 18.d0*gl%z1_PCSAFT(1) &
            & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1) - 1.d0)**3 - 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1) - 1.d0)**3 + 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**2 - 24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**5) - 18.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - &
            & 1.d0)**4) - 12.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3/(gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**3) - &
            & 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3/(gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**2) + 54.d0*gl%z2_PCSAFT(1) &
            & **2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4) + 36.d0*gl%z2_PCSAFT(1) &
            & **2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3) + 18.d0* &
            & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2) - &
            & 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3) - 18.d0* &
            & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) + 6.d0* &
            & gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) - 2.d0*gl%z3_PCSAFT(1)**3*(gl%z0_PCSAFT(1) &
            & - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0)**3 + 3.d0*gl%z3_PCSAFT(1)**2*( &
            & gl%z0_PCSAFT(1) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(1)/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0)**2 - 3.d0*gl%z3_PCSAFT(1)*(-6.d0*gl%z2_PCSAFT(1)** &
            & 3*gl%z3_PCSAFT(1)**2/gl%z3_PCSAFT(1)**4 + 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/ &
            & gl%z3_PCSAFT(1)**3 - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(1)**2/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0) - (24.d0* &
            & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3/gl%z3_PCSAFT(1)**5 - 54.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)**2/gl%z3_PCSAFT(1)**4 + 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)/ &
            & gl%z3_PCSAFT(1)**3 - 6.d0*gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0) + 3.d0* &
            & gl%z0_PCSAFT(1)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1) - 1.d0)**3 - 6.d0* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**2 - 6.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/( &
            & gl%z3_PCSAFT(1) - 1.d0) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4) &
            & - 4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3) - 2.d0* &
            & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2) + 12.d0*gl%z2_PCSAFT(1)** &
            & 2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3) + 6.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) - 6.d0*gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)**2/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) - gl%z3_PCSAFT(1)**2*(gl%z0_PCSAFT(1) - &
            & gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0*gl%z3_PCSAFT(1)*(gl%z0_PCSAFT(1) + &
            & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)/gl%z3_PCSAFT(1) &
            & **2)/(gl%z3_PCSAFT(1) - 1.d0) + (-6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2/gl%z3_PCSAFT(1)**4 + 12.d0* &
            & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/gl%z3_PCSAFT(1)**3 - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(1)** &
            & 2/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1) - 6.d0*gl%z0_PCSAFT(1)**2*(-3.d0* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/( &
            & gl%z3_PCSAFT(1) - 1.d0) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/( &
            & gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2) + gl%z3_PCSAFT(1)*(gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/( &
            & gl%z3_PCSAFT(1) - 1.d0) + (gl%z0_PCSAFT(1) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/gl%z3_PCSAFT(1)**3 - 3.d0* &
            & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(1)/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1)**2 + &
            & 6.d0*gl%z0_PCSAFT(1)**3*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) - gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1) &
            & *(gl%z3_PCSAFT(1) - 1.d0)**2) + (gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + &
            & 1.d0))/gl%z0_PCSAFT(1)**3)/gl%z0_PCSAFT(1)
    end if
    
    ! 9: 3RD DERIVATIVE OF g_ij WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T, z0_PCSAFT_TT, z1_PCSAFT_TT, z2_PCSAFT_TT, z3_PCSAFT_TT
    !           z0_PCSAFT_TTT, z1_PCSAFT_TTT, z2_PCSAFT_TTT, z3_PCSAFT_TTT
    if (GETDERAHS(9) .eq. 1) then
        gl%ahs_PCSAFT(9) = (18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**3/(gl%z3_PCSAFT(1) - 1.d0)**4 - 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)/(gl%z3_PCSAFT(1) - 1.d0)**3 + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(9)/(gl%z3_PCSAFT(1) - 1.d0)**2 - 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2/( &
            & gl%z3_PCSAFT(1) - 1.d0)**3 + 9.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 9.d0 &
            & *gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0)**2 - 3.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(9)/(gl%z3_PCSAFT(1) - 1.d0) - 18.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2/( &
            & gl%z3_PCSAFT(1) - 1.d0)**3 + 9.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)/(gl%z3_PCSAFT(1) - 1.d0)**2 + &
            & 18.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0)**2 - 9.d0*gl%z1_PCSAFT(4)* &
            & gl%z2_PCSAFT(5)/(gl%z3_PCSAFT(1) - 1.d0) + 9.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - &
            & 1.d0)**2 - 9.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0) - 3.d0*gl%z1_PCSAFT(9)* &
            & gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) - 24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - &
            & 1.d0)**5) + 18.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4 &
            & ) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(9)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3) - 18.d0*gl%z2_PCSAFT(1) &
            & **3*gl%z3_PCSAFT(4)**3/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**4) + 12.d0*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3) - gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(9)/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) - 12.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4) &
            & **3/(gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**3) + 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
            & /(gl%z3_PCSAFT(1)**3*(gl%z3_PCSAFT(1) - 1.d0)**2) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3/(gl%z3_PCSAFT(1)**4 &
            & *(gl%z3_PCSAFT(1) - 1.d0)**2) + 54.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2/(gl%z3_PCSAFT(1)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**4) - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5)/(gl%z3_PCSAFT(1)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**3) + 36.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2/(gl%z3_PCSAFT(1)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0)**3) - 9.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5)/(gl%z3_PCSAFT(1)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2) + 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2/(gl%z3_PCSAFT(1)**3*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2) - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**3) - 9.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)**2*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)** &
            & 2) - 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3) - 18.d0 &
            & *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) + 18.d0* &
            & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) + 6.d0*gl%z2_PCSAFT(4) &
            & **3/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) - 2.d0*gl%z3_PCSAFT(4)**3*(gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/ &
            & gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0)**3 + 3.d0*gl%z3_PCSAFT(4)**2*(gl%z0_PCSAFT(4) + 2.d0*gl%z2_PCSAFT(1)** &
            & 3*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - &
            & 1.d0)**2 + 3.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/( &
            & gl%z3_PCSAFT(1) - 1.d0)**2 - 3.d0*gl%z3_PCSAFT(4)*(gl%z0_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5)/ &
            & gl%z3_PCSAFT(1)**3 - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2/gl%z3_PCSAFT(1)**4 + 12.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)/gl%z3_PCSAFT(1)**2 - &
            & 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0) - 3.d0*gl%z3_PCSAFT(5)*( &
            & gl%z0_PCSAFT(4) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/ &
            & gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0) - gl%z3_PCSAFT(9)*(gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1) &
            & **2)/(gl%z3_PCSAFT(1) - 1.d0) - (gl%z0_PCSAFT(9) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(9)/ &
            & gl%z3_PCSAFT(1)**3 - 18.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)/gl%z3_PCSAFT(1)**4 + 24.d0* &
            & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**3/gl%z3_PCSAFT(1)**5 + 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)* &
            & gl%z3_PCSAFT(5)/gl%z3_PCSAFT(1)**3 - 54.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2/gl%z3_PCSAFT(1)** &
            & 4 + 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(9)/gl%z3_PCSAFT(1)**2 + 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 &
            & - 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)/gl%z3_PCSAFT(1)**2 - 6.d0*gl%z2_PCSAFT(4)**3/gl%z3_PCSAFT(1) &
            & **2)*dlog(-gl%z3_PCSAFT(1) + 1.d0) + 3.d0*gl%z0_PCSAFT(4)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
            & **2/(gl%z3_PCSAFT(1) - 1.d0)**3 - 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)/(gl%z3_PCSAFT(1) - 1.d0)**2 &
            & - 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 3.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(5)/(gl%z3_PCSAFT(1) - 1.d0) - 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0) &
            & **2 + 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0) + 3.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)/( &
            & gl%z3_PCSAFT(1) - 1.d0) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4) + &
            & 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3) - 4.d0*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)**2/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5)/( &
            & gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2/(gl%z3_PCSAFT(1)**3*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2) + 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) &
            & - 1.d0)**3) + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0) &
            & **2) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) - 6.d0*gl%z2_PCSAFT(1) &
            & *gl%z2_PCSAFT(4)**2/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) - gl%z3_PCSAFT(4)**2*(gl%z0_PCSAFT(1) - &
            & gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0*gl%z3_PCSAFT(4)*(gl%z0_PCSAFT(4) + 2.d0* &
            & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/gl%z3_PCSAFT(1)**2)/( &
            & gl%z3_PCSAFT(1) - 1.d0) + gl%z3_PCSAFT(5)*(gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - &
            & 1.d0) + (gl%z0_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5)/gl%z3_PCSAFT(1)**3 - 6.d0*gl%z2_PCSAFT(1)**3 &
            & *gl%z3_PCSAFT(4)**2/gl%z3_PCSAFT(1)**4 + 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 &
            & - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)/gl%z3_PCSAFT(1)**2 - 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2/gl%z3_PCSAFT(1) &
            & **2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1) + 3.d0*gl%z0_PCSAFT(5)*(-3.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0) &
            & + 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**3) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2) &
            & - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) + gl%z3_PCSAFT(4)*( &
            & gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0) + (gl%z0_PCSAFT(4) + 2.d0*gl%z2_PCSAFT(1) &
            & **3*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/gl%z3_PCSAFT(1)**2)*dlog( &
            & -gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1) + gl%z0_PCSAFT(9)*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) &
            & - 1.d0) - gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) + (gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/ &
            & gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1) - 6.d0*gl%z0_PCSAFT(4)**2*(-3.d0* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)/( &
            & gl%z3_PCSAFT(1) - 1.d0) + 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**3) + gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1)**2 &
            & *(gl%z3_PCSAFT(1) - 1.d0)**2) - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2 &
            & ) + gl%z3_PCSAFT(4)*(gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0) + ( &
            & gl%z0_PCSAFT(4) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/ &
            & gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1)**2 - 6.d0*gl%z0_PCSAFT(4)* &
            & gl%z0_PCSAFT(5)*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) - gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**2) + (gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + &
            & 1.d0))/gl%z0_PCSAFT(1)**2 + 6.d0*gl%z0_PCSAFT(4)**3*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) - &
            & gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) + (gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)** &
            & 2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1)**3)/gl%z0_PCSAFT(1)
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T
    if (GETDERAHS(10) .eq. 1) then
        gl%ahs_PCSAFT(10) = (18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0)**4 - 24.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0)**3 + 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
            & /(gl%z3_PCSAFT(1) - 1.d0)**2 - 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1) - 1.d0)**3 + &
            & 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**2 - 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1) - 1.d0)**3 + 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) &
            & **2 - 24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0)**5 + 18.d0*gl%z2_PCSAFT(1)**3 &
            & *gl%z3_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0)**4 + 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) &
            & - 1.d0)**4 - 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)/(gl%z3_PCSAFT(1) - 1.d0)**3 - 2.d0*gl%z3_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(4)*(gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0)**3 + gl%z3_PCSAFT(1) &
            & **2*(gl%z0_PCSAFT(4) + 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)/gl%z3_PCSAFT(1)**3 - 3.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(4)/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z0_PCSAFT(1) - &
            & gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0*gl%z0_PCSAFT(4)*(3.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) - gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) + ( &
            & gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0))/gl%z0_PCSAFT(1) - 2.d0* &
            & gl%z0_PCSAFT(4)*(-3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 6.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0) + 2.d0*gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1) - 1.d0)**3 - 2.d0*gl%z2_PCSAFT(1)**3/( &
            & gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) + gl%z3_PCSAFT(1)*(gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/( &
            & gl%z3_PCSAFT(1) - 1.d0) + (gl%z0_PCSAFT(1) - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)*dlog(-gl%z3_PCSAFT(1) + 1.d0 &
            & ))/gl%z0_PCSAFT(1) + gl%z0_PCSAFT(4)*(6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2/(gl%z3_PCSAFT(1) - 1.d0)**3 - &
            & 12.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)/(gl%z3_PCSAFT(1) &
            & - 1.d0) - 6.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)/(gl%z3_PCSAFT(1) - 1.d0)**4 + 8.d0*gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1) - &
            & 1.d0)**3 - 2.d0*gl%z2_PCSAFT(1)**3/(gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**2) - gl%z3_PCSAFT(1)**2*(gl%z0_PCSAFT(1) &
            & - gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0)**2 + 2.d0*gl%z3_PCSAFT(1)*(gl%z0_PCSAFT(1) - &
            & gl%z2_PCSAFT(1)**3/gl%z3_PCSAFT(1)**2)/(gl%z3_PCSAFT(1) - 1.d0))/gl%z0_PCSAFT(1))/gl%z0_PCSAFT(1)
    end if
    
    ! 11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^4
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT
    if (GETDERAHS(11) .eq. 1) then
        gl%ahs_PCSAFT(11) = 6.d0*gl%z3_PCSAFT(1)**2*(20.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1) - 4.d0*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(3.d0*gl%z1_PCSAFT(1) &
            & *gl%z3_PCSAFT(1)**2 + 4.d0*gl%z2_PCSAFT(1)**2) + (gl%z3_PCSAFT(1) - 1.d0)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + &
            & 12.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - gl%z2_PCSAFT(1)**3))/(gl%z0_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**6)
    end if

    ! 12: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, DEL AND DEL, MULTIPLIED BY TAU*DEL^3
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T
    if (GETDERAHS(12) .eq. 1) then
        gl%ahs_PCSAFT(12) = 2.d0*(60.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 12.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0) &
            & *(3.d0*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 5.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + (gl%z3_PCSAFT(1) - 1.d0)**3*(-18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4) - 9.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 9.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
            & **2 + gl%z2_PCSAFT(1)**2*(-2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) - 3.d0* &
            & gl%z3_PCSAFT(4)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(18.d0* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + &
            & 3.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 3.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 9.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3 &
            & )))/(gl%z0_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**6)
    end if
    
    ! 13: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, DEL, AND DEL, MULTIPLIED BY TAU^2*DEL^2
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T, z0_PCSAFT_TT, z1_PCSAFT_TT, z2_PCSAFT_TT, z3_PCSAFT_TT
    if (GETDERAHS(13) .eq. 1) then
        gl%ahs_PCSAFT(13) = (120.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2 - 24.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) &
            & - 1.d0)*(3.d0*gl%z1_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) &
            & + 4.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)) - &
            & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)**4*(6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
            & *(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2)) + 6.d0*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z1_PCSAFT(1) &
            & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5) &
            & *gl%z3_PCSAFT(1) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
            & + gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)) - 2.d0*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**3*(12.d0* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + &
            & 24.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2 &
            & + 24.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
            & **2 + 3.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 6.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5) + 12.d0* &
            & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2 + gl%z3_PCSAFT(5)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + &
            & 6.d0*gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) &
            & + 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
            & **2*gl%z3_PCSAFT(4) + 6.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)**3 &
            & *gl%z3_PCSAFT(5) + 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1) + gl%z3_PCSAFT(4)**2*( &
            & gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + 4.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)** &
            & 3*(-gl%z2_PCSAFT(1)**2*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + 2.d0* &
            & gl%z3_PCSAFT(4)*(-gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + gl%z2_PCSAFT(1)**3)) + 2.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3) + 2.d0*gl%z3_PCSAFT(4)*( &
            & gl%z3_PCSAFT(1) - 1.d0)**4*(2.d0*gl%z2_PCSAFT(1)**2*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)* &
            & gl%z3_PCSAFT(1)) + gl%z3_PCSAFT(4)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)))/(gl%z0_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)**2*(gl%z3_PCSAFT(1) - 1.d0)**6)
    end if
    
    ! 14: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, TAU AND DEL, MULTIPLIED BY TAU^3*DEL
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z0_PCSAFT_T, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T, z0_PCSAFT_TT, z1_PCSAFT_TT, z2_PCSAFT_TT, z3_PCSAFT_TT,
    !          z0_PCSAFT_TTT, z1_PCSAFT_TTT, z2_PCSAFT_TTT, z3_PCSAFT_TTT
    if (GETDERAHS(14) .eq. 1) then
        gl%ahs_PCSAFT(14) = (gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(9) - 6.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 4.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + 6.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(4)**3 + 18.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 6.d0*gl%z0_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) - 12.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 - 18.d0*gl%z0_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 4.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) + 6.d0*gl%z0_PCSAFT(1) &
            & *gl%z3_PCSAFT(4)**3 + 6.d0*gl%z0_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + gl%z0_PCSAFT(1)*gl%z3_PCSAFT(9) - &
            & 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + 54.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 &
            & *gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) - 72.d0* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 - 108.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) - 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) + 72.d0* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**3 + 54.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
            & 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9) - 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(5) + 54.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 54.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 108.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)**2 - 54.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 54.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2 + 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) - 18.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 54.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(4) - 54.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 18.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(1)**4 - 12.d0*gl%z1_PCSAFT(1)* &
            & gl%z2_PCSAFT(9)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(1)**2 - 12.d0*gl%z1_PCSAFT(1) &
            & *gl%z2_PCSAFT(9)*gl%z3_PCSAFT(1) + 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(9) - 18.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 54.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + &
            & 54.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 108.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 54.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 54.d0* &
            & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 18.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) - 36.d0* &
            & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 108.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 108.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 36.d0* &
            & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 9.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**4 - 36.d0* &
            & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**3 + 54.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2 - &
            & 36.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) + 9.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(5) - 18.d0* &
            & gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 54.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
            & **2*gl%z3_PCSAFT(4) - 54.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 18.d0*gl%z1_PCSAFT(5) &
            & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 9.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z1_PCSAFT(5)* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3 + 54.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 - 36.d0* &
            & gl%z1_PCSAFT(5)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 9.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(4) + 3.d0*gl%z1_PCSAFT(9)* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**4 - 12.d0*gl%z1_PCSAFT(9)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3 + 18.d0*gl%z1_PCSAFT(9) &
            & *gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - 12.d0*gl%z1_PCSAFT(9)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 3.d0*gl%z1_PCSAFT(9)* &
            & gl%z2_PCSAFT(1) - 2.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(9) + 18.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)** &
            & 2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 12.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) - 24.d0* &
            & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 - 108.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) - 18.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) + 144.d0*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(4)**3 + 90.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 8.d0*gl%z2_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(9) - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(5) + 54.d0* &
            & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 108.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4) &
            & *gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) - 324.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 &
            & - 162.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 270.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2 + 72.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) - 18.d0* &
            & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 108.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)* &
            & gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 162.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 72.d0 &
            & *gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(1)**4 &
            & - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(1)**3 + 36.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)* &
            & gl%z3_PCSAFT(1)**2 - 30.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(1) + 9.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z2_PCSAFT(9) - 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 216.d0* &
            & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 324.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 144.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4) + 18.d0*gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**4 - 108.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)* &
            & gl%z3_PCSAFT(1)**3 + 216.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2 - 180.d0*gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) + 54.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5) + 6.d0* &
            & gl%z2_PCSAFT(4)**3*gl%z3_PCSAFT(1)**4 - 36.d0*gl%z2_PCSAFT(4)**3*gl%z3_PCSAFT(1)**3 + 72.d0*gl%z2_PCSAFT(4)**3 &
            & *gl%z3_PCSAFT(1)**2 - 60.d0*gl%z2_PCSAFT(4)**3*gl%z3_PCSAFT(1) + 18.d0*gl%z2_PCSAFT(4)**3)/(gl%z0_PCSAFT(1)*( &
            & gl%z3_PCSAFT(1)**6 - 6.d0*gl%z3_PCSAFT(1)**5 + 15.d0*gl%z3_PCSAFT(1)**4 - 20.d0*gl%z3_PCSAFT(1)**3 + 15.d0*gl%z3_PCSAFT(1) &
            & **2 - 6.d0*gl%z3_PCSAFT(1) + 1.d0))
    end if
    
    ! 15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
    ! requires z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT, z1_PCSAFT_T, z2_PCSAFT_T, z3_PCSAFT_T, z1_PCSAFT_TT, z2_PCSAFT_TT, z3_PCSAFT_TT,
    !          z1_PCSAFT_TTT, z2_PCSAFT_TTT, z3_PCSAFT_TTT, z1_PCSAFT_TTTT, z2_PCSAFT_TTTT, z3_PCSAFT_TTTT
    if (GETDERAHS(15) .eq. 1) then
        gl%ahs_PCSAFT(15) = (-72.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**4*(gl%z3_PCSAFT(1) - 1.d0) + 120.d0*gl%z2_PCSAFT(1) &
            & **3*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**4 + 96.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(4)**4*( &
            & gl%z3_PCSAFT(1) - 1.d0) + 72.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**4*(gl%z3_PCSAFT(1) - 1.d0)**2 &
            & + 48.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**4*(gl%z3_PCSAFT(1) - 1.d0)**3 - 144.d0*gl%z2_PCSAFT(1) &
            & **2*gl%z3_PCSAFT(1)**5*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 8.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - &
            & 1.d0)**3*(-9.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) - 18.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%z3_PCSAFT(4)*(-2.d0* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1))) - 4.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3* &
            & gl%z3_PCSAFT(9)*(gl%z3_PCSAFT(1) - 1.d0)**5*(2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
            & ) - 36.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 6.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0 &
            & )**2*(4.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5)**2 &
            & + 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)* &
            & gl%z3_PCSAFT(4)**2 + 36.d0*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)**2) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**3*( &
            & gl%z3_PCSAFT(1) - 1.d0)**4*(4.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 3.d0*gl%z2_PCSAFT(1)**2* &
            & gl%z3_PCSAFT(5)**2 + 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 18.d0*gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)**2 + 6.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(2.d0*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(4) - 3.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)) + 36.d0*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)**2) + 6.d0* &
            & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5)*(gl%z3_PCSAFT(1) - 1.d0)**5*(6.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) &
            & **2 - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + &
            & 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**2)) + 6.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
            & *gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)**4*(4.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - gl%z3_PCSAFT(1)*(6.d0 &
            & *gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 3.d0*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)** &
            & 2))) + 36.d0*gl%z3_PCSAFT(1)**6*gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(5) + 2.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 2.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**5*(gl%z1_PCSAFT(1)*gl%z2_PCSAFT(15) + 4.d0* &
            & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(9) + 6.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(5) + 4.d0*gl%z1_PCSAFT(9)* &
            & gl%z2_PCSAFT(4) + gl%z1_PCSAFT(15)*gl%z2_PCSAFT(1)) + 3.d0*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**4*( &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(15) + 4.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(9) + 6.d0* &
            & gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(5) + 4.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(4) + 4.d0* &
            & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9) + 12.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 12.d0* &
            & gl%z1_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 6.d0*gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 12.d0* &
            & gl%z1_PCSAFT(5)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 4.d0*gl%z1_PCSAFT(9)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) - 6.d0* &
            & gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0)**3*(4.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + &
            & 3.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)**2 + 12.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)* &
            & gl%z3_PCSAFT(5) + 6.d0*gl%z1_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)**2 + 12.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 12.d0*gl%z1_PCSAFT(4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2 + 6.d0* &
            & gl%z1_PCSAFT(5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2) + 3.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - 1.d0)**4*( &
            & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(15) + 8.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(9) + 6.d0*gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(5)**2 + 12.d0*gl%z2_PCSAFT(4)**2*gl%z2_PCSAFT(5)) - 2.d0*gl%z3_PCSAFT(1)**5*(gl%z3_PCSAFT(1) - &
            & 1.d0)**3*(gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(15) + 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(9) + &
            & 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(5) + 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)* &
            & gl%z3_PCSAFT(4) + 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 72.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
            & gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 24.d0*gl%z2_PCSAFT(4)**3*gl%z3_PCSAFT(4)) + 6.d0*gl%z3_PCSAFT(1)**4* &
            & gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)**2*(-18.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5) - 36.d0*gl%z2_PCSAFT(1) &
            & **2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%z3_PCSAFT(4)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3 &
            & )) - gl%z3_PCSAFT(1)**4*gl%z3_PCSAFT(15)*(gl%z3_PCSAFT(1) - 1.d0)**5*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
            & gl%z2_PCSAFT(1)**3) + gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**4*(-gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(15) - &
            & 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(9) - 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)* &
            & gl%z3_PCSAFT(5) - 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(4) - 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) &
            & **2*gl%z3_PCSAFT(5) - 72.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) - 24.d0* &
            & gl%z2_PCSAFT(4)**3*gl%z3_PCSAFT(4) + 4.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - &
            & gl%z2_PCSAFT(1)**3) + 3.d0*gl%z3_PCSAFT(5)**2*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1)**3)) + 4.d0* &
            & gl%z3_PCSAFT(1)**4*(gl%z3_PCSAFT(1) - 1.d0)**3*(4.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 3.d0* &
            & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(5)**2 + 36.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
            & + 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)**2 + 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(4)**2 - 3.d0*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(5)*(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**2 - gl%z2_PCSAFT(1) &
            & **3)) - 4.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)**5*(24.d0*gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4) &
            & **3 - 18.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 3.d0*gl%z2_PCSAFT(4) &
            & *gl%z3_PCSAFT(4)) + 2.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 9.d0*gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 9.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 18.d0*gl%z2_PCSAFT(4)**2 &
            & *gl%z3_PCSAFT(4)) - 3.d0*gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9) + 6.d0*gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)**3)) + (gl%z3_PCSAFT(1) - 1.d0)**6*(120.d0* &
            & gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(4)**4 - 144.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(gl%z2_PCSAFT(1)* &
            & gl%z3_PCSAFT(5) + 2.d0*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)) + 6.d0*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*(4.d0*gl%z2_PCSAFT(1)** &
            & 2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + 3.d0*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5)**2 + 36.d0*gl%z2_PCSAFT(1)* &
            & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 18.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)**2 + &
            & 36.d0*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)**2) + 3.d0*gl%z3_PCSAFT(1)**4*(gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(15) &
            & + 8.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(9) + 6.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)**2 + 12.d0* &
            & gl%z2_PCSAFT(4)**2*gl%z2_PCSAFT(5)) - 2.d0*gl%z3_PCSAFT(1)**3*(gl%z2_PCSAFT(1)**3*gl%z3_PCSAFT(15) + 12.d0* &
            & gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(9) + 18.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(5) &
            & + 12.d0*gl%z2_PCSAFT(1)**2*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(4) + 36.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)**2* &
            & gl%z3_PCSAFT(5) + 72.d0*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 24.d0*gl%z2_PCSAFT(4)** &
            & 3*gl%z3_PCSAFT(4)))*dlog(-gl%z3_PCSAFT(1) + 1.d0))/(gl%z0_PCSAFT(1)*gl%z3_PCSAFT(1)**6*(gl%z3_PCSAFT(1) - 1.d0) &
            & **6)
    end if
    
!DEC$ END IF
end subroutine AHSDERIVS
    
subroutine GIIDERIVS(gl,GETDERGII)

! Henning Markgraf, June 2016
    
    ! g_ii_hs: site-site radial distribution function of hard sphere fluid
    ! defined by eq. A.7 in Gross, Sadowski 2001:
    ! g_ij_hs = 1/(1-zeta_3) + (d_i*d_j/(d_i+d_j)) * 3*zeta_2/(1-zeta_3)**2
    !           + (d_i*d_j/(d_i+d_j))**2 * 2*zeta_2**2/(1-zeta_3)**3
    ! dependent on D and T






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERGII
    !output: gii_PCSAFT (module variable)
    ! working variable
    integer :: i
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! III. calculate the derivatives of g_ii_hs for every component i
    do i = gl%n_start, gl%n_end
        !1: g_ii_hs
        if (GETDERGII(1) .eq. 1) then
            gl%gii_PCSAFT(i,1) = (1.d0/2.d0)*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2/(-gl%z3_PCSAFT(1) + 1.d0)**3 + (3.d0/2.d0)* &
                            & gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)/(-gl%z3_PCSAFT(1) + 1.d0)**2 + 1.d0/(-gl%z3_PCSAFT(1) + 1.d0)
        end if
    
        ! 2: 1ST DERIVATIVE OF g_ij_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
        ! requires d_nrsubst, z2_PCSAFT, z3_PCSAFT
        if (GETDERGII(2) .eq. 1) then
            gl%gii_PCSAFT(i,2) = 0.5d0*(3.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) - 2.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - &
                             & 1.d0)*(gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1) + 3.d0*gl%z3_PCSAFT(1)) + (gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1) + 2.d0* &
                             & gl%z3_PCSAFT(1)))/(gl%z3_PCSAFT(1) - 1.d0)**4
        end if
    
        ! 3: 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGII(3) .eq. 1) then
            gl%gii_PCSAFT(i,3) = (-6.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)* &
                             & (2.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1) + 3.d0*gl%z3_PCSAFT(1)) - (gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2 + 6.d0 &
                             & *gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)**2))/(gl%z3_PCSAFT(1) - 1.d0)**5
        end if
    
        ! 4: 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst, z3_PCSAFT(4), z2_PCSAFT(4), d_nrsubst(i,4)
        if (GETDERGII(4) .eq. 1) then
            gl%gii_PCSAFT(i,4) = (1.d0/2.d0)*(3.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 2.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) &
                             & - 1.d0)*(gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4) + gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1) + 3.d0*gl%z3_PCSAFT(4)) + (gl%z3_PCSAFT(1) - 1.d0)**2 &
                             & *(3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(4)))/(gl%z3_PCSAFT(1) - 1.d0)**4
        end if
    
        ! 5: 2ND DERIVATIVE OF g_ij WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst, z3_PCSAFT(4), z2_PCSAFT(4), d_nrsubst(i,4), z3_PCSAFT(5), z2_PCSAFT(5), d_nrsubst(i,5)
        if (GETDERGII(5) .eq. 1) then
            gl%gii_PCSAFT(i,5) = (1.d0/2.d0)*(-12.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*( &
                             & gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 4.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 4.d0* &
                             & gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%z3_PCSAFT(4)**2) + (gl%z3_PCSAFT(1) - 1.d0)**3*(3.d0*gl%di_PCSAFT(i,1)* &
                             & gl%z2_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(5)) - 2.d0 &
                             & *(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2 + &
                             & 4.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1) &
                             & *gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)**2 + 6.d0* &
                             & gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%z3_PCSAFT(4)**2))/(gl%z3_PCSAFT(1) - 1.d0)**5
        end if
    
        ! 6: 2ND MIXED DERIVATIVE OF g_ij WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst, z3_PCSAFT(4), z2_PCSAFT(4), d_nrsubst(i,4)
        if (GETDERGII(6) .eq. 1) then
            gl%gii_PCSAFT(i,6) = (1.d0/2.d0)*(-12.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1) &
                             & *(gl%z3_PCSAFT(1) - 1.d0)*(3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 2.d0* &
                             & gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 6.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)) + (gl%z3_PCSAFT(1) - 1.d0)**3*(3.d0*gl%di_PCSAFT(i,1) &
                             & *gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(4)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(2.d0* &
                             & gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + 2.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2 + 6.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 3.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 2.d0* &
                             & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)))/(gl%z3_PCSAFT(1) - 1.d0)**5
        end if
    
        ! 7: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO rho, T, AND T, MULTIPLIED BY T*T*rho
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst, z3_PCSAFT(4), z2_PCSAFT(4), z3_PCSAFT(5), z2_PCSAFT(5), d_nrsubst(i,4), d_nrsubst(i,5)
        if (GETDERGII(7) .eq. 1) then
            gl%gii_PCSAFT(i,7) = (1.d0/2.d0)*(36.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - 6.d0*gl%di_PCSAFT(i,1)* &
                             & gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 8.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(4)**2 + 4.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2) + 2.d0*gl%z3_PCSAFT(1)*(12.d0*gl%di_PCSAFT(i,1)**2* &
                             & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(5) + 4.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0 &
                             & *gl%z3_PCSAFT(4)**2) - (gl%z3_PCSAFT(1) - 1.d0)**3*(3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,4)* &
                             & gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(5)) + 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*( &
                             & gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2 + 4.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1) &
                             & *gl%z2_PCSAFT(4) + gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1) &
                             & *gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)**2 + 6.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
                             & + 2.d0*gl%z3_PCSAFT(4)**2)) + (gl%z3_PCSAFT(1) - 1.d0)**4*(3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,4)* &
                             & gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(5)) - 4.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*( &
                             & gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2 + 4.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1) &
                             & *gl%z2_PCSAFT(4) + gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1) &
                             & *gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)**2 + 6.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
                             & + 2.d0*gl%z3_PCSAFT(4)**2) + (gl%z3_PCSAFT(1) - 1.d0)**2*(9.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5) &
                             & + 36.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)* &
                             & gl%z3_PCSAFT(1) + 2.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1) + 36.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(4) + 8.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 2.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)* &
                             & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 54.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(4)**2 + 12.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1) &
                             & **2*gl%z3_PCSAFT(1) + 12.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
                             & **2))/(gl%z3_PCSAFT(1) - 1.d0)**6
        end if
    
        ! 8: 3RD DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGII(8) .eq. 1) then
            gl%gii_PCSAFT(i,8) = 3.d0*gl%z3_PCSAFT(1)*(10.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 - 12.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
                             & gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1) + gl%z3_PCSAFT(1)) + (gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(i,1)**2* &
                             & gl%z2_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)**2))/(gl%z3_PCSAFT(1) - 1.d0)**6
        end if
    
        ! 9: 3RD DERIVATIVE OF g_ij WITH RESPECT TO T, MULTIPLIED BY T^3
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst, z3_PCSAFT(4), z2_PCSAFT(4), z3_PCSAFT(5), z2_PCSAFT(5),
        ! z3_PCSAFT(9), z2_PCSAFT(9), d_nrsubst(i,4), d_nrsubst(i,5), d_nrsubst(i,9)
        if (GETDERGII(9) .eq. 1) then
            gl%gii_PCSAFT(i,9) = (1.d0/2.d0)*(60.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 - 36.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)* &
                             & gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%z3_PCSAFT(4)**2) + (gl%z3_PCSAFT(1) - 1.d0 &
                             & )**4*(3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(9) + 9.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(5) + 9.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(4) + &
                             & 3.d0*gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(9)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(i,1)**2* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(9) + 3.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)**2 + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)* &
                             & gl%z2_PCSAFT(4) + gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9) + 9.d0* &
                             & gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 9.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,4)**2* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(5) + 18.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 9.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
                             & + 6.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(9) + 6.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1) &
                             & *gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4) &
                             & *gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 24.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 6.d0* &
                             & gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
                             & + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2 + 6.d0*gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + &
                             & 18.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 4.d0*gl%z3_PCSAFT(4)**3))/(gl%z3_PCSAFT(1) - 1.d0)**6

        end if
    
        ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst, z3_PCSAFT(4), z2_PCSAFT(4), d_nrsubst(i,4)
        if (GETDERGII(10) .eq. 1) then
            gl%gii_PCSAFT(i,10) = (30.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 12.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*( &
                             & gl%z3_PCSAFT(1) - 1.d0)*(3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(i,4)* &
                             & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 3.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(i,1)**2* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0 &
                             & *gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 3.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
                             & ) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(i,1)**2* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 4.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) + 12.d0*gl%di_PCSAFT(i,1)* &
                             & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,4)* &
                             & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 2.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)))/(gl%z3_PCSAFT(1) - 1.d0)**6
        end if
    
        ! 11: 4TH DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^4
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGII(11) .eq. 1) then
            gl%gii_PCSAFT(i,11) = 12.d0*gl%z3_PCSAFT(1)**2*(-15.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 + 5.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
                             & *(gl%z3_PCSAFT(1) - 1.d0)*(4.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1) + 3.d0*gl%z3_PCSAFT(1)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0* &
                             & gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2 + 6.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + gl%z3_PCSAFT(1)**2))/(gl%z3_PCSAFT(1) - 1.d0) &
                             & **7
        end if
    
        ! 12: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, DEL AND DEL, MULTIPLIED BY TAU*DEL^3
        ! requires z2_PCSAFT, z3_PCSAFT, di_PCSAFT, z2_PCSAFT(4), z3_PCSAFT(4), di_PCSAFT(i,4)
        if (GETDERGII(12) .eq. 1) then
            gl%gii_PCSAFT(i,12) = 3.d0*(-60.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**3*gl%z3_PCSAFT(4) + 10.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) &
                             & **2*(gl%z3_PCSAFT(1) - 1.d0)*(9.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + &
                             & 2.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 6.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)) - 4.d0*gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0 &
                             & )**2*(9.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
                             & gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)* &
                             & gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + &
                             & 2.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 2.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)** &
                             & 2*gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)** &
                             & 2 + 3.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 2.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)))/(gl%z3_PCSAFT(1) - 1.d0 &
                             & )**7
        end if
    
        ! 13: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, DEL, AND DEL, MULTIPLIED BY TAU^2*DEL^2
        ! requires z2_PCSAFT, z3_PCSAFT, di_PCSAFT, z2_PCSAFT(4), z3_PCSAFT(4), di_PCSAFT(i,4), z2_PCSAFT(5), z3_PCSAFT(5), di_PCSAFT(i,5)
        if (GETDERGII(13) .eq. 1) then
            gl%gii_PCSAFT(i,13) = (30.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 2.d0*gl%di_PCSAFT(i,1)* &
                             & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 - gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - 8.d0*gl%di_PCSAFT(i,1)* &
                             & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 4.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 4.d0*gl%di_PCSAFT(i,1)* &
                             & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) - 4.d0* &
                             & gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 - 6.d0*gl%z3_PCSAFT(1) &
                             & *gl%z3_PCSAFT(4)**2) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + gl%di_PCSAFT(i,1) &
                             & **2*gl%z2_PCSAFT(4)**2 + 4.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1) &
                             & **2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,1)* &
                             & gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)**2 + 6.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + &
                             & 6.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 3.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1)* &
                             & gl%z3_PCSAFT(5) + 2.d0*gl%z3_PCSAFT(4)**2) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(3.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)** &
                             & 2*gl%z3_PCSAFT(5) + 12.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(i,1)**2* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) + 4.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1) + 12.d0*gl%di_PCSAFT(i,1)* &
                             & gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 16.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + &
                             & 4.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) + 12.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + &
                             & 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 24.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0* &
                             & gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2 + 4.d0*gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) + 24.d0*gl%di_PCSAFT(i,4)* &
                             & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,5)* &
                             & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 + 2.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 8.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2) &
                             & - 12.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(3.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1) &
                             & **2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 12.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)* &
                             & gl%z3_PCSAFT(4) + gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)**2 + gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2 &
                             & *gl%z3_PCSAFT(1)**2 + 12.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4) &
                             & *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2 + gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)**2 + 3.d0* &
                             & gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 &
                             & + 6.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)** &
                             & 2 + 6.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 2.d0*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2)) &
                             & /(gl%z3_PCSAFT(1) - 1.d0)**7
        end if
    
        ! 14: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, TAU AND DEL, MULTIPLIED BY TAU^3*DEL
        ! requires z2_PCSAFT, z3_PCSAFT, di_PCSAFT, z2_PCSAFT(4), z3_PCSAFT(4), di_PCSAFT(i,4), z2_PCSAFT(5), z3_PCSAFT(5), di_PCSAFT(i,5),
        !          z2_PCSAFT(9), z3_PCSAFT(9), di_PCSAFT(i,9)
        if (GETDERGII(14) .eq. 1) then
            gl%gii_PCSAFT(i,14) = (1.d0/2.d0)*(-240.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3 + 12.d0*gl%di_PCSAFT(i,1)* &
                             & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%z3_PCSAFT(1) - 1.d0)*(9.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 25.d0* &
                             & gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 18.d0* &
                             & gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 18.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2) + 2.d0*gl%z3_PCSAFT(1)* &
                             & (-60.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**3 + 36.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*( &
                             & gl%z3_PCSAFT(1) - 1.d0)*(gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 2.d0* &
                             & gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%z3_PCSAFT(4)**2) - (gl%z3_PCSAFT(1) - 1.d0)**4*(3.d0*gl%di_PCSAFT(i,1)* &
                             & gl%z2_PCSAFT(9) + 9.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(5) + 9.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,9)* &
                             & gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(9)) + 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)* &
                             & gl%z2_PCSAFT(9) + 3.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)* &
                             & gl%z2_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)**2 + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)* &
                             & gl%z2_PCSAFT(4) + gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9) + 9.d0* &
                             & gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 9.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,4)**2* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(5) + 18.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 9.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
                             & + 6.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)) - 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(9) + 6.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1) &
                             & *gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4) &
                             & *gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 24.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 6.d0* &
                             & gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) &
                             & + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2 + 6.d0*gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + &
                             & 18.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 4.d0*gl%z3_PCSAFT(4)**3)) + (gl%z3_PCSAFT(1) - 1.d0)**5*( &
                             & 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(9) + 9.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(5) + 9.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(4) + 3.d0* &
                             & gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(9)) - 4.d0*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1) &
                             & *gl%z2_PCSAFT(9) + 3.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)* &
                             & gl%z2_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)**2 + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)* &
                             & gl%z2_PCSAFT(4) + gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9) + 9.d0* &
                             & gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 9.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,4)**2* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2 + 9.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(5) + 18.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 9.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
                             & + 6.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)) + (gl%z3_PCSAFT(1) - 1.d0)**3*(9.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(9) + 54.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 54.d0*gl%di_PCSAFT(i,1)**2* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(1) + &
                             & 54.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)* &
                             & gl%z3_PCSAFT(1) + 54.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 216.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)* &
                             & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 12.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) + 12.d0*gl%di_PCSAFT(i,1) &
                             & *gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1) + 54.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 12.d0 &
                             & *gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 2.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(9) + 162.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
                             & *gl%z3_PCSAFT(5) + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 162.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)* &
                             & gl%z3_PCSAFT(4)**2 + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 54.d0*gl%di_PCSAFT(i,4)**2* &
                             & gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 12.d0*gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 6.d0*gl%di_PCSAFT(i,4) &
                             & *gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1) + 18.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 162.d0 &
                             & *gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 36.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + &
                             & 18.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 12.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
                             & 36.d0*gl%z3_PCSAFT(4)**3) - 6.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)* &
                             & gl%z3_PCSAFT(9) + 24.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(i,1)**2* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 48.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)* &
                             & gl%z3_PCSAFT(4)**2 + 6.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,1)** &
                             & 2*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)* &
                             & gl%z3_PCSAFT(5) + 48.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 24.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)* &
                             & gl%z3_PCSAFT(4) + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 48.d0*gl%di_PCSAFT(i,1)* &
                             & gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**3 + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 6.d0*gl%di_PCSAFT(i,4) &
                             & **2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 18.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)** &
                             & 2 + 4.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**3))/(gl%z3_PCSAFT(1) - 1.d0)**7
        end if
    
        ! 15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
        ! requires z2_PCSAFT, z3_PCSAFT, di_PCSAFT, z2_PCSAFT(4), z3_PCSAFT(4), di_PCSAFT(i,4), z2_PCSAFT(5), z3_PCSAFT(5), di_PCSAFT(i,5)
        !          z2_PCSAFT(9), z3_PCSAFT(9), di_PCSAFT(i,9), z2_PCSAFT(15), z3_PCSAFT(15), di_PCSAFT(i,15)
        if (GETDERGII(15) .eq. 1) then
            gl%gii_PCSAFT(i,15) = (1.d0/2.d0)*(-360.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**4 + 120.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(4)**2*(gl%z3_PCSAFT(1) - 1.d0)*(3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 4.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4) &
                             & *gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0*gl%z3_PCSAFT(4)**2) + (gl%z3_PCSAFT(1) - &
                             & 1.d0)**5*(3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(15) + 12.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(9) + 18.d0*gl%di_PCSAFT(i,5)* &
                             & gl%z2_PCSAFT(5) + 12.d0*gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(i,15)*gl%z2_PCSAFT(1) + 2.d0* &
                             & gl%z3_PCSAFT(15)) - 2.d0*(gl%z3_PCSAFT(1) - 1.d0)**4*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(15) + 4.d0* &
                             & gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(9) + 3.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(5)**2 + 8.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(9) + 24.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5) + 12.d0*gl%di_PCSAFT(i,1)* &
                             & gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + 12.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(4)**2 + 8.d0*gl%di_PCSAFT(i,1)* &
                             & gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,15)*gl%z2_PCSAFT(1)**2 + 3.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)* &
                             & gl%z3_PCSAFT(15) + 12.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(9) + 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(5)* &
                             & gl%z3_PCSAFT(5) + 12.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(4) + 12.d0*gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)* &
                             & gl%z2_PCSAFT(5) + 12.d0*gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(4)**2 + 24.d0*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)* &
                             & gl%z2_PCSAFT(4) + 4.d0*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(1)**2 + 12.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9) &
                             & + 36.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 36.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 3.d0* &
                             & gl%di_PCSAFT(i,5)**2*gl%z2_PCSAFT(1)**2 + 18.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 36.d0*gl%di_PCSAFT(i,5)* &
                             & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 12.d0*gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 8.d0*gl%z3_PCSAFT(4)* &
                             & gl%z3_PCSAFT(9) + 6.d0*gl%z3_PCSAFT(5)**2) + 3.d0*(gl%z3_PCSAFT(1) - 1.d0)**3*(gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)** &
                             & 2*gl%z3_PCSAFT(15) + 8.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(9) + 12.d0*gl%di_PCSAFT(i,1)**2* &
                             & gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(5) + 8.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(9)*gl%z3_PCSAFT(4) &
                             & + 12.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(5) + 24.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5)* &
                             & gl%z3_PCSAFT(4) + 8.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 48.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1) &
                             & *gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) + 48.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 48.d0 &
                             & *gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4) + 12.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(5) + 48.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 8.d0*gl%di_PCSAFT(i,1)* &
                             & gl%di_PCSAFT(i,9)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4) + 24.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(9) + &
                             & 18.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)**2 + 72.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + &
                             & 36.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4)**2 + 12.d0*gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5) + &
                             & 48.d0*gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 24.d0*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(4) + 72.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 72.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(4) &
                             & *gl%z3_PCSAFT(4)**2 + 36.d0*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 24.d0*gl%z3_PCSAFT(4)**2* &
                             & gl%z3_PCSAFT(5)) - 12.d0*(gl%z3_PCSAFT(1) - 1.d0)**2*(4.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)* &
                             & gl%z3_PCSAFT(9) + 3.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(5)**2 + 24.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)* &
                             & gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 12.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) &
                             & **2 + 12.d0*gl%di_PCSAFT(i,1)**2*gl%z2_PCSAFT(4)**2*gl%z3_PCSAFT(4)**2 + 24.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5) + 48.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**2 + &
                             & 12.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,5)*gl%z2_PCSAFT(1)**2*gl%z3_PCSAFT(4)**2 + 36.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2* &
                             & gl%z3_PCSAFT(5) + 24.d0*gl%di_PCSAFT(i,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)**3 + 12.d0*gl%di_PCSAFT(i,4)**2*gl%z2_PCSAFT(1)**2* &
                             & gl%z3_PCSAFT(4)**2 + 24.d0*gl%di_PCSAFT(i,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**3 + 4.d0*gl%z3_PCSAFT(4)**4))/( &
                             & gl%z3_PCSAFT(1) - 1.d0)**7
        end if
    end do
    
!DEC$ END IF
end subroutine GIIDERIVS


subroutine GIJDERIVS(gl,D, GETDERGIJ)

! Henning Markgraf, June 2016
    
    ! g_ii_hs: site-site radial distribution function of hard sphere fluid
    ! defined by eq. A.7 in Gross, Sadowski 2001:
    ! g_ij_hs = 1/(1-z3_PCSAFT(1)) + (d_i*d_j/(d_i+d_j)) * 3*z2_PCSAFT(1)/(1-z3_PCSAFT(1))**2
    !           + (d_i*d_j/(d_i+d_j))**2 * 2*z2_PCSAFT(1)**2/(1-z3_PCSAFT(1))**3
    ! dependent on D and T






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    double precision, intent (in) :: D
    integer, dimension (nderivs), intent (in) :: GETDERGIJ
    !output: gij_PCSAFT (module variable)
    ! working variable
    integer :: i
    double precision :: z311,z312,z313,z314,z315,z316, di2,dj2,dipj,dipj2,dipj3,dipj4,dipj5,dimj, ditpdjt,ditpdjt2,ditpdjt3
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    z311 = gl%z3_PCSAFT(1) - 1.d0
    z312 = z311*z311
    z313 = z311*z312
    z314 = z311*z313
    z315 = z311*z314
    z316 = z311*z315
    
    di2 = gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,1)
    dj2 = gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,1)
    dipj = gl%di_PCSAFT(1,1)+gl%di_PCSAFT(2,1)
    dipj2 = dipj*dipj
    dipj3 = dipj*dipj2
    dipj4 = dipj*dipj3
    dipj5 = dipj*dipj4
    dimj = gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,1)
    ditpdjt = gl%di_PCSAFT(1,4)+gl%di_PCSAFT(2,4)
    ditpdjt2 = ditpdjt*ditpdjt
    ditpdjt3 = ditpdjt*ditpdjt2
   
    ! III. calculate the derivatives of g_ij_hs for every component i
    do i = 1, 1 ! n_start, n_end nur bin�r
        
        !1: g_ii_hs
        if (GETDERGIJ(1) .eq. 1) then
            gl%gij_PCSAFT(i,1) = 2.d0*(gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,1)/dipj)**2 *gl%z2_PCSAFT(1)**2/(-z313) + 3.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,1)/dipj*gl%z2_PCSAFT(1)/z312 + 1.d0/(-z311)
        end if
    
        ! 2: 1ST DERIVATIVE OF g_ij_hs WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
        ! requires d_nrsubst, z2_PCSAFT, z3_PCSAFT
        if (GETDERGIJ(2) .eq. 1) then
            gl%gij_PCSAFT(i,2) = D*(6.d0*di2 * dj2 * gl%z2_PCSAFT(1)**2 * gl%z3_PCSAFT(2) / (dipj2 *z314) + 4.d0*di2*dj2 * gl%z2_PCSAFT(1)*gl%z2_PCSAFT(2)/(dipj2 * (-z313)) &
            & + 6.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(2)/(dipj*(-z313)) + 3.d0*dimj*gl%z2_PCSAFT(2)/(dipj*z312) + gl%z3_PCSAFT(2)/z312)
        end if
    
        ! 3: 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGIJ(3) .eq. 1) then      
            gl%gij_PCSAFT(i,3) = 2.d0*(-12.d0*di2 * dj2 * gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)**2 + 12.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0) - 2.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *z312 &
            & + 9.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 *dipj*(gl%z3_PCSAFT(1) - 1.d0) - 6.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*dipj*z312 &
            & - gl%z3_PCSAFT(1)**2 *dipj2 *z312)/(dipj2*z315)
        end if
        
        ! 4: 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGIJ(4) .eq. 1) then      
            gl%gij_PCSAFT(i,4) = (6.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4)*dipj + 4.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *ditpdjt*(gl%z3_PCSAFT(1) - 1.d0) &
            & - 6.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*dipj2 *(gl%z3_PCSAFT(1) - 1.d0) - 3.d0*dimj*gl%z2_PCSAFT(1)*dipj*ditpdjt*z312 &
            & - 4.d0*dimj*gl%z2_PCSAFT(1)*dipj*(gl%z3_PCSAFT(1) - 1.d0)*(dimj*gl%z2_PCSAFT(4) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1) + gl%di_PCSAFT(1,4)* gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)) + gl%z3_PCSAFT(4)*dipj3 *z312 &
            & + 3.d0*dipj2 *z312 *(dimj*gl%z2_PCSAFT(4) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1) + gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)))/(dipj3 *z314)
        end if
        
        ! 5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGIJ(5) .eq. 1) then      
            gl%gij_PCSAFT(i,5) = (-24.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4)**2 *dipj2 - 24.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4)*dipj*ditpdjt*(gl%z3_PCSAFT(1) - 1.d0)&
            & - 12.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *ditpdjt2 *z312 + 18.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 *dipj3 *(gl%z3_PCSAFT(1) - 1.d0)&
            & + 6.d0*dimj*gl%z2_PCSAFT(1)*dipj2 *(gl%z3_PCSAFT(1) - 1.d0)*(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 4.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4))&
            & + 6.d0*dimj*gl%z2_PCSAFT(1)*dipj*ditpdjt2 *z313 + 4.d0*dimj*gl%z2_PCSAFT(1)*dipj*z312 *(dimj*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) &
            & + 4.d0*dimj*gl%z2_PCSAFT(4)*ditpdjt + 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*ditpdjt + 4.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*ditpdjt) &
            & - 2.d0*gl%z3_PCSAFT(4)**2 *dipj4 *z312 + gl%z3_PCSAFT(5)*dipj4 *z313 &
            & + 3.d0*dipj3 *z313 *(dimj*gl%z2_PCSAFT(5) + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(4) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1) + gl%di_PCSAFT(1,5)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)) &
            & - 6.d0*dipj3 *z312 *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) &
            & - 3.d0*dipj2 *z313 *(dimj*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) + 2.d0*dimj*gl%z2_PCSAFT(4)*ditpdjt + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*ditpdjt + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*ditpdjt) &
            & + 4.d0*dipj2 *z312 *(-di2 *dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) - di2 *dj2 *gl%z2_PCSAFT(4)**2 - 4.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) - di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)**2 &
            & - di2 *gl%di_PCSAFT(2,4)**2 *gl%z2_PCSAFT(1)**2 - 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) - 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 - gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,5)*dj2 *gl%z2_PCSAFT(1)**2 + 3.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*ditpdjt - gl%di_PCSAFT(1,4)**2*dj2*gl%z2_PCSAFT(1)**2))&
            & /(dipj4 *z315)
        end if
        
        ! 6: 2ND MIXED DERIVATIVE OF F WITH RESPECT TO D AND T, MULTIPLIED BY T*D
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGIJ(5) .eq. 1) then      
            gl%gij_PCSAFT(i,6) = (-24.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*dipj - 12.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)*ditpdjt*(gl%z3_PCSAFT(1) - 1)&
            & + 8.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *ditpdjt*z312 &
            & + 18.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*dipj2 *(gl%z3_PCSAFT(1) - 1.d0)&
            & - 3.d0*dimj*gl%z2_PCSAFT(1)*dipj*ditpdjt*z313 &
            & + 2.d0*dimj*gl%z2_PCSAFT(1)*dipj*z312 *(-4.d0*dimj*gl%z2_PCSAFT(4) - 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1) - 4.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1) + 3.d0*gl%z3_PCSAFT(1)*ditpdjt)&
            & + 6.d0*dimj*gl%z2_PCSAFT(1)*dipj*(gl%z3_PCSAFT(1) - 1.d0)*(3.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1))&
            & - 2.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*dipj3 *z312 + gl%z3_PCSAFT(4)*dipj3 *z313 &
            & + 3.d0*dipj2 *z313 *(dimj*gl%z2_PCSAFT(4) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1) + gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)) &
            & - 6.d0*dipj2 *z312 *(2.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)))&
            & /(dipj3 *z315)
        end if
        
        ! 7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGIJ(7) .eq. 1) then      
            gl%gij_PCSAFT(i,7) = (72.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 *dipj2 &
            & + 48.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*dipj*ditpdjt*(gl%z3_PCSAFT(1) - 1.d0) &
            & + 12.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)*ditpdjt2 *z312 &
            & - 24.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *ditpdjt2 *z313 &
            & - 36.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 *dipj3 *(gl%z3_PCSAFT(1) - 1.d0) &
            & - 12.d0*dimj*gl%z2_PCSAFT(1)*dipj2 *(gl%z3_PCSAFT(1) - 1.d0)*(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 8.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 &
            & + 4.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4))&
            & + 6.d0*dimj*gl%z2_PCSAFT(1)*dipj*ditpdjt2 *z314 &
            & + 8.d0*dimj*gl%z2_PCSAFT(1)*dipj*(gl%z3_PCSAFT(1) - 1)**3 *(dimj*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) &
            & + 4.d0*dimj*gl%z2_PCSAFT(4)*ditpdjt + 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*ditpdjt &
            & + 4.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*ditpdjt) &
            & - 4.d0*dimj*gl%z2_PCSAFT(1)*dipj*(gl%z3_PCSAFT(1) - 1)**2 *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) & 
            & + 18.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*ditpdjt + 4.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*ditpdjt &
            & + 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*ditpdjt + 4.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*ditpdjt) &
            & + 2.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2 *dipj4 *z312 &
            & + 2.d0*gl%z3_PCSAFT(1)*(24.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4)**2 *dipj2 &
            & + 24.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4)*dipj*ditpdjt*(gl%z3_PCSAFT(1) - 1.d0) &
            & + 12.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *ditpdjt2 *z312 &
            & - 18.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 *dipj3 *(gl%z3_PCSAFT(1) - 1.d0) &
            & - 6.d0*dimj*gl%z2_PCSAFT(1)*dipj2 *(gl%z3_PCSAFT(1) - 1.d0) &
            & *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 4.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 4.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) &
            & - 6.d0*dimj*gl%z2_PCSAFT(1)*dipj*ditpdjt2 *z313 &
            & - 4.d0*dimj*gl%z2_PCSAFT(1)*dipj*z312 &
            & * (dimj*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) + 4.d0*dimj*gl%z2_PCSAFT(4)*ditpdjt + 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*ditpdjt + 4.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*ditpdjt) &
            & + 2.d0*gl%z3_PCSAFT(4)**2 *dipj4 *z312 - gl%z3_PCSAFT(5)*dipj4 *z313 &
            & - 3.d0*dipj3 *z313 &
            & *(dimj*gl%z2_PCSAFT(5) + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(4) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1) + gl%di_PCSAFT(1,5)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)) &
            & + 6.d0*dipj3 *z312 &
            & *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) &
            & + 3.d0*dipj2 *z313 *(dimj*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) + 2.d0*dimj*gl%z2_PCSAFT(4)*ditpdjt &
            & + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*ditpdjt + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*ditpdjt) &
            & + 4.d0*dipj2 *z312 *(di2 *dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) + di2 *dj2 *gl%z2_PCSAFT(4)**2 &
            & + 4.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)**2 + di2 *gl%di_PCSAFT(2,4)**2 *gl%z2_PCSAFT(1)**2 &
            & + 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) + 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,5)*dj2 *gl%z2_PCSAFT(1)**2 &
            & - 3.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*ditpdjt + gl%di_PCSAFT(1,4)**2 *dj2 *gl%z2_PCSAFT(1)**2)) &
            & - 4.d0*gl%z3_PCSAFT(4)**2 *dipj4 *z313 + gl%z3_PCSAFT(5)*dipj4 *z314 &
            & + 3.d0*dipj3 *z314 &
            & *(dimj*gl%z2_PCSAFT(5) + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(4) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1) + gl%di_PCSAFT(1,5)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)) &
            & - 12.d0*dipj3 *z313 &
            & *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) &
            & + 6.d0*dipj3 *z312 &
            & *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) + 9.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 + 2.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
            & + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)) &
            & - 3.d0*dipj2 *z314 *(dimj*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) + 2.d0*dimj*gl%z2_PCSAFT(4)*ditpdjt &
            & + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*ditpdjt + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*ditpdjt) &
            & + 8.d0*dipj2 *z313 *(-di2 *dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) - di2 *dj2 *gl%z2_PCSAFT(4)**2 &
            & - 4.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) - di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)**2 - di2 *gl%di_PCSAFT(2,4)**2 *gl%z2_PCSAFT(1)**2 &
            & - 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) - 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 - gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,5)*dj2 *gl%z2_PCSAFT(1)**2 &
            & + 3.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*ditpdjt - gl%di_PCSAFT(1,4)**2 *dj2 *gl%z2_PCSAFT(1)**2) &
            & + 2.d0*dipj2 *z312 *(9.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(5) &
            & + 36.d0*di2 *dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 2*di2 *dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(1) + 2.d0*di2 *dj2 *gl%z2_PCSAFT(4)**2 *gl%z3_PCSAFT(1) &
            & + 36.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4) + 8.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 2.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1) &
            & + 2.d0*di2 *gl%di_PCSAFT(2,4)**2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1) + 36.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4) + 8.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) &
            & + 8.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1) + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,5)*dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1) &
            & - 6.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*ditpdjt + 2.d0*gl%di_PCSAFT(1,4)**2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)))&
            & /(dipj4 *z316)
        end if
        
        ! 8: 3RD DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGIJ(8) .eq. 1) then      
            gl%gij_PCSAFT(i,8) = 6.d0*gl%z3_PCSAFT(1)*(20.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)**2 - 24.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0) &
            & + 6.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *z312 - 12.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 *dipj*(gl%z3_PCSAFT(1) - 1.d0) &
            & + 9.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*dipj*z312 + gl%z3_PCSAFT(1)**2 *dipj2 *z312) &
            & /(dipj2 *z316)
        end if
        
        ! 9: 3RD DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^3
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGIJ(9) .eq. 1) then      
            gl%gij_PCSAFT(i,9) = (120.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4)**3 *dipj3 &
            & + 144.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4)**2 *dipj2 *ditpdjt*(gl%z3_PCSAFT(1) - 1.d0) &
            & + 108.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4)*dipj*ditpdjt2 *z312 &
            & + 48.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *ditpdjt3 *z313 &
            & - 72.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**3 *dipj4 *(gl%z3_PCSAFT(1) - 1.d0) &
            & - 72.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*dipj3 *(gl%z3_PCSAFT(1) - 1.d0) &
            & *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 2.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) &
            & + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) &
            & - 36.d0*dimj*gl%z2_PCSAFT(1)*dipj2 *z312 *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) &
            & + dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*ditpdjt + 4.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*ditpdjt &
            & + 4.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*ditpdjt + 4.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*ditpdjt) &
            & - 18.d0*dimj*gl%z2_PCSAFT(1)*dipj*ditpdjt3 *z314 &
            & - 36.d0*dimj*gl%z2_PCSAFT(1)*dipj*ditpdjt*z313 &
            & *(dimj*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) + 2.d0*dimj*gl%z2_PCSAFT(4)*ditpdjt &
            & + 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*ditpdjt + 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*ditpdjt) &
            & + 6.d0*gl%z3_PCSAFT(4)**3 *dipj5 *z312 - 6.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*dipj5 *z313 &
            & + 54.d0*gl%z3_PCSAFT(4)*dipj4 *z312 &
            & *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) &
            & + gl%z3_PCSAFT(9)*dipj5 *z314 + 3.d0*dipj4 *z314 &
            & *(dimj*gl%z2_PCSAFT(9) + 3.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(5) + 3.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(4) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,9)*gl%z2_PCSAFT(1) &
            & + 3.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(5) + 6.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(4) + 3.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1) + 3.d0*gl%di_PCSAFT(1,5)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(4) &
            & + 3.d0*gl%di_PCSAFT(1,5)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1) + gl%di_PCSAFT(1,9)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)) &
            & - 6.d0*dipj4 *z313 *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(9) + 3.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) &
            & + 3.d0*dimj*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) &
            & + 3.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5) + 6.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) &
            & + 6.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 3.d0*gl%di_PCSAFT(1,5)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)) &
            & - 3.d0*dipj3 *z314 *(dimj*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,9) + gl%di_PCSAFT(2,9)) &
            & + 3.d0*dimj*gl%z2_PCSAFT(4)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) + 3.d0*dimj*gl%z2_PCSAFT(5)*ditpdjt &
            & + 3.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) + 6.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(4)*ditpdjt &
            & + 3.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)*ditpdjt + 3.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) &
            & + 6.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(4)*ditpdjt + 6.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*ditpdjt &
            & + 3.d0*gl%di_PCSAFT(1,5)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*ditpdjt) + 2.d0*dipj3 *z313 &
            & *(-2.d0*di2 *dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(9) - 6.d0*di2* dj2 *gl%z2_PCSAFT(4)*gl%z2_PCSAFT(5) &
            & - 12.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) - 12.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(4)**2 &
            & - 12.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) - 2.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,9)*gl%z2_PCSAFT(1)**2 &
            & - 12.d0*di2 *gl%di_PCSAFT(2,4)**2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) - 6.d0*di2 *gl%di_PCSAFT(2,4)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)**2 &
            & - 12.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5) - 12.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(4)**2 &
            & - 48.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) - 12.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)**2 &
            & - 12.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,4)**2 *gl%z2_PCSAFT(1)**2 - 12.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,5)*dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) &
            & - 12.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,5)*gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 - 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,9)*dj2 *gl%z2_PCSAFT(1)**2 &
            & + 9.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) + 9.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(5)*ditpdjt &
            & + 18.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4)*ditpdjt + 18.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*ditpdjt &
            & - 12.d0*gl%di_PCSAFT(1,4)**2 *dj2*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4) - 12.d0*gl%di_PCSAFT(1,4)**2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 - 6.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(1,5)*dj2 *gl%z2_PCSAFT(1)**2 &
            & + 18.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*ditpdjt) + 6.d0*dipj3 *z312 &
            & *(di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(9) + 6.d0*di2 *dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(5) &
            & + 6.d0*di2 *dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*gl%z3_PCSAFT(4) + 6.d0*di2 *dj2 *gl%z2_PCSAFT(4)**2 *gl%z3_PCSAFT(4) &
            & + 6.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(5) + 24.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) &
            & + 6.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4) + 6.d0*di2 *gl%di_PCSAFT(2,4)**2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4) &
            & + 6.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(5) + 24.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(4) &
            & + 24.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4) + 6.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,5)*dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4) &
            & - 9.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)**2 *ditpdjt + 6.d0*gl%di_PCSAFT(1,4)**2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(4)) &
            & + 18.d0*dipj2 *ditpdjt*z314 &
            & *(dimj*gl%z2_PCSAFT(1)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) + dimj*gl%z2_PCSAFT(4)*ditpdjt + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*ditpdjt + gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*ditpdjt) &
            & + 4.d0*dipj2 *z313 *(di2 *dj2 *gl%z2_PCSAFT(1)**2 *(gl%di_PCSAFT(1,9) + gl%di_PCSAFT(2,9)) &
            & + 6.d0*di2 *dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) + 6.d0*di2 *dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(5)*ditpdjt &
            & + 6.d0*di2 *dj2 *gl%z2_PCSAFT(4)**2 *ditpdjt + 6.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 *(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) &
            & + 24.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*ditpdjt + 6.d0*di2 *gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,5)*gl%z2_PCSAFT(1)**2 *ditpdjt &
            & + 6.d0*di2 *gl%di_PCSAFT(2,4)**2 *gl%z2_PCSAFT(1)**2 *ditpdjt + 6.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(1)**2 *(gl%di_PCSAFT(1,5) + gl%di_PCSAFT(2,5)) &
            & + 24.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*dj2 *gl%z2_PCSAFT(1)*gl%z2_PCSAFT(4)*ditpdjt + 24.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)**2 *ditpdjt &
            & + 6.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(1,5)*dj2 *gl%z2_PCSAFT(1)**2 *ditpdjt - 9.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4)*ditpdjt2 &
            & + 6.d0*gl%di_PCSAFT(1,4)**2 *dj2 *gl%z2_PCSAFT(1)**2 *ditpdjt)) & 
            & /(dipj5 *z316)
        end if
        
        ! 10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGIJ(10) .eq. 1) then      
            gl%gij_PCSAFT(i,10) = 2.d0*(60.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)**2 *gl%z3_PCSAFT(4)*dipj &
            & + 24.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)**2 *ditpdjt*(gl%z3_PCSAFT(1) - 1.d0) &
            & - 24.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)*ditpdjt*z312 &
            & + 4.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *ditpdjt*z313 &
            & - 36.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 *gl%z3_PCSAFT(4)*dipj2 *(gl%z3_PCSAFT(1) - 1.d0) &
            & - 24.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*dipj*(gl%z3_PCSAFT(1) - 1.d0) &
            & *(3.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)) &
            & + 2.d0*dimj*gl%z2_PCSAFT(1)*dipj*z313 &
            & *(-2.d0*dimj*gl%z2_PCSAFT(4) - 2.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1) - 2.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1) + 3.d0*gl%z3_PCSAFT(1)*ditpdjt) &
            & + 3.d0*dimj*gl%z2_PCSAFT(1)*dipj*(gl%z3_PCSAFT(1) - 1)**2 & 
            & *(6.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + 8.d0*dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + 8.d0*gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + 8.d0*gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) - 3.d0*gl%z3_PCSAFT(1)**2 *ditpdjt) &
            & + 3*gl%z3_PCSAFT(1)**2 *gl%z3_PCSAFT(4)*dipj3 *z312 - 2.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*dipj3 *z313 &
            & + 9.d0*gl%z3_PCSAFT(1)*dipj2 *z312 &
            & *(4.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)) &
            & - 6.d0*dipj2 *z313 &
            & *(dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(4) + dimj*gl%z2_PCSAFT(4)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(1,1)*gl%di_PCSAFT(2,4)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1) + gl%di_PCSAFT(1,4)*gl%di_PCSAFT(2,1)*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1))) &
            & /(dipj3 *z316)
        end if
        
        ! 11: 4TH DERIVATIVE OF F WITH RESPECT TO D, MULTIPLIED BY D^4
        ! requires z2_PCSAFT, z3_PCSAFT, d_nrsubst
        if (GETDERGIJ(11) .eq. 1) then      
            gl%gij_PCSAFT(i,11) = 24.d0*gl%z3_PCSAFT(1)**2 *(-30.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)**2 + 40.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *gl%z3_PCSAFT(1)*(gl%z3_PCSAFT(1) - 1.d0) &
            & - 12.d0*di2 *dj2 *gl%z2_PCSAFT(1)**2 *z312 + 15.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)**2 *dipj*(gl%z3_PCSAFT(1) - 1.d0) &
            & - 12.d0*dimj*gl%z2_PCSAFT(1)*gl%z3_PCSAFT(1)*dipj*z312 - gl%z3_PCSAFT(1)**2 *dipj2 *z312)&
            & /(dipj2 *(gl%z3_PCSAFT(1) - 1.d0)**7)
        end if
        
    end do
    
!DEC$ END IF
end subroutine GIJDERIVS

subroutine ZETADERIVS(gl,D, GETDERZETA, n_zeta)

! Henning Markgraf, June 2016   

    ! zeta: abbreviation defined by eq. 9 in Gross, Sadowski 2001:
    ! zeta_n = pi/6*rho*Sum(x_i*m_i*d_i**n) , n element of {0, 1, 2, 3}
    ! dependent on D and T






implicit none

    type(type_gl) :: gl


    !Declarations
    !input
    integer, intent (in) :: n_zeta
    double precision, intent (in) :: D
    integer, dimension (nderivs), intent (in) :: GETDERZETA
    !output: z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT (module variables)
    !working variables
    double precision :: derivative
    double precision :: sum
    integer :: i
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !calculate the derivatives of zeta_0, zeta_1, zeta_2 and zeta_3
    ! 1: zeta_n
    ! 2: 1ST DERIVATIVE OF zeta_n WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D == zeta_n
    if (GETDERZETA(1) .eq. 1 .OR. GETDERZETA(2) .eq. 1) then
        !iterate over all components of the mixture to find the sum of x_i*m_i*d_i
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + gl%molfractions(i)*gl%mPCSAFT(i)*gl%di_PCSAFT(i,1)**n_zeta
        end do
        derivative = piPCSAFT/6.d0*D*sum
        ! d(zeta_n)/d(D) * D = zeta_n
        select case (n_zeta)
            case (0)
                gl%z0_PCSAFT(1) = derivative
                gl%z0_PCSAFT(2) = derivative
            case (1)
                gl%z1_PCSAFT(1) = derivative
                gl%z1_PCSAFT(2) = derivative
            case (2)
                gl%z2_PCSAFT(1) = derivative
                gl%z2_PCSAFT(2) = derivative
            case (3)
                gl%z3_PCSAFT(1) = derivative
                gl%z3_PCSAFT(2) = derivative
        end select
    end if    
    
    ! 3: 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! DERZETA(3) = 0.d0 (do nothing, already initialized as 0.d0)
    
    
    ! 4: 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERZETA(4) .eq. 1 .OR. GETDERZETA(6) .EQ. 1) then
        !iterate over all components of the mixture to find the sum n*d_i**(n-1)+m_i*x_i*d_i_T
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + n_zeta*gl%molfractions(i)*gl%mPCSAFT(i)*gl%di_PCSAFT(i,1)**(n_zeta-1)*gl%di_PCSAFT(i,4)
        end do
        
        derivative = piPCSAFT*D/6.d0*sum
        ! Derivative number 6 the same as 4
        select case (n_zeta)
            case (0)
                gl%z0_PCSAFT(4) = derivative
                gl%z0_PCSAFT(6) = derivative
            case (1)
                gl%z1_PCSAFT(4) = derivative
                gl%z1_PCSAFT(6) = derivative
            case (2)
                gl%z2_PCSAFT(4) = derivative
                gl%z2_PCSAFT(6) = derivative
            case (3)
                gl%z3_PCSAFT(4) = derivative
                gl%z3_PCSAFT(6) = derivative
        end select
    end if
    
    ! 5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERZETA(5) .eq. 1 .OR. GETDERZETA(7) .EQ. 1) then
        !iterate over all components of the mixture to find the sum n*d_i**(n-1)*m_i*x_i* &&
        !                       && (n*d_i_T**2/d_i + d_i_TT - d_i_T**2/d_i)
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + n_zeta*gl%molfractions(i)*gl%mPCSAFT(i)*gl%di_PCSAFT(i,1)**(n_zeta-1) * &
                    & (n_zeta*gl%di_PCSAFT(i,4)**2/gl%di_PCSAFT(i,1) + gl%di_PCSAFT(i,5) - gl%di_PCSAFT(i,4)**2/gl%di_PCSAFT(i,1))
        end do
        ! factor T**2 not needed because included in DERDI(4)**2 and DERDI(5)
        derivative = piPCSAFT*D/6.d0*sum
        ! same as DERZETA(7)
        select case (n_zeta)
            case (0)
                gl%z0_PCSAFT(5) = derivative
                gl%z0_PCSAFT(7) = derivative
            case (1)
                gl%z1_PCSAFT(5) = derivative
                gl%z1_PCSAFT(7) = derivative
            case (2)
                gl%z2_PCSAFT(5) = derivative
                gl%z2_PCSAFT(7) = derivative
            case (3)
                gl%z3_PCSAFT(5) = derivative
                gl%z3_PCSAFT(7) = derivative
        end select
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF F WITH RESPECT TO rho AND T, MULTIPLIED BY T*rho
    ! DERZETA(6) = DERZETA(4)
    
    ! 7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO rho, T, AND T, MULTIPLIED BY T*T*rho
    ! DERZETA(7) = DERZETA(5)
    
    ! 8: 3RD DERIVATIVE OF F WITH RESPECT TO rho, MULTIPLIED BY rho^3
    ! DERZETA(8) = 0.d0 (do nothing, already initialized as 0.d0)
    
    ! 9: 3RD DERIVATIVE OF F WITH RESPECT TO T, MULTIPLIED BY T^3
    if (GETDERZETA(9) .eq. 1 .OR. GETDERZETA(14) .EQ. 1) then
        !iterate over all components of the mixture to find the sum
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + n_zeta*gl%di_PCSAFT(i,1)**(n_zeta-1)*gl%molfractions(i)*gl%mPCSAFT(i)* &
                & (n_zeta**2*gl%di_PCSAFT(i,4)**3/gl%di_PCSAFT(i,1)**2 + 3.d0*n_zeta*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,5)/gl%di_PCSAFT(i,1) - &
                & 3.d0*n_zeta*gl%di_PCSAFT(i,4)**3/gl%di_PCSAFT(i,1)**2 + gl%di_PCSAFT(i,9) - 3.d0*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,5)/gl%di_PCSAFT(i,1) + &
                & 2.d0*gl%di_PCSAFT(i,4)**3/gl%di_PCSAFT(i,1)**2)
        end do
        derivative = piPCSAFT*D/6.d0*sum
        ! Derivative number 14 the same as 9
        select case (n_zeta)
            case (0)
                gl%z0_PCSAFT(9) = derivative
                gl%z0_PCSAFT(14) = derivative
            case (1)
                gl%z1_PCSAFT(9) = derivative
                gl%z1_PCSAFT(14) = derivative
            case (2)
                gl%z2_PCSAFT(9) = derivative
                gl%z2_PCSAFT(14) = derivative
            case (3)
                gl%z3_PCSAFT(9) = derivative
                gl%z3_PCSAFT(14) = derivative
        end select
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO T, rho, AND rho, MULTIPLIED BY T*rho*rho
    ! DERZETA(10) = 0.d0 (do nothing, already initialized as 0.d0)

    ! 11: all derivatives with respect to rho higher than second order of zeta_n equal zero
    ! 12: all derivatives with respect to rho higher than second order of zeta_n equal zero 
    ! 13:  all derivatives with respect to rho higher than second order of zeta_n equal zero
    
    ! 14: 4TH MIXED DERIVATIVE OF F WITH RESPECT TO TAU, TAU, TAU AND DEL, MULTIPLIED BY TAU^3*DEL
    ! DERZETA(14) = DERZETA(9)
    
    ! 15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
    if (GETDERZETA(15) .eq. 1) then
        !iterate over all components of the mixture to find the sum
        sum = 0.d0
        do i = gl%n_start, gl%n_end
            sum = sum + gl%di_PCSAFT(i,1)**(n_zeta - 4)*n_zeta*(gl%di_PCSAFT(i,1)**3* &
                & gl%di_PCSAFT(i,15) + gl%di_PCSAFT(i,1)**2*(4.d0*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,9)*n_zeta - 4.d0*gl%di_PCSAFT(i,4)*gl%di_PCSAFT(i,9) + 3.d0* &
                & gl%di_PCSAFT(i,5)**2*n_zeta - 3.d0*gl%di_PCSAFT(i,5)**2) + 6.d0*gl%di_PCSAFT(i,1)*gl%di_PCSAFT(i,4)**2*gl%di_PCSAFT(i,5)*(n_zeta** &
                & 2 - 3.d0*n_zeta + 2.d0) + gl%di_PCSAFT(i,4)**4*(n_zeta**3 - 6.d0*n_zeta**2 + 11.d0*n_zeta &
                & - 6.d0))*gl%mPCSAFT(i)*gl%molfractions(i)
        end do
        derivative = piPCSAFT*D/6.d0*sum
        select case (n_zeta)
            case (0)
                gl%z0_PCSAFT(15) = derivative
            case (1)
                gl%z1_PCSAFT(15) = derivative
            case (2)
                gl%z2_PCSAFT(15) = derivative
            case (3)
                gl%z3_PCSAFT(15) = derivative
        end select
    end if
    
!DEC$ END IF
end subroutine

subroutine DIDERIVS(gl,T, GETDERDI)

! Henning Markgraf, June 2016

    ! d_i: segment diameter as in eq. A.9 in Gross, Sadowski 2001
    ! d_i = sigma_i*(1-0.12*exp(-3*epsilon_i/(k*T)))
    ! dependent only on T






implicit none

    type(type_gl) :: gl


    ! Declarations
    ! input
    double precision, intent (in) :: T
    integer, dimension (nderivs), intent (in) :: GETDERDI
    ! output: di_PCSAFT (module variable)
    ! working variables
    integer :: i
    integer:: errorfld 
    
!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! calculate d_i for every component and all derivatives with respect to T
    do i = gl%n_start, gl%n_end
        ! 1: d_i
        if (GETDERDI(1) .eq. 1) then
            gl%di_PCSAFT(i,1) = gl%sigPCSAFT(i)*(1.d0-0.12d0 * dexp(-3.d0*gl%epskPCSAFT(i)/T))
        end if
    
        ! 2: all rho derivatives of d_i equal zero
        ! 3: all rho derivatives of d_i equal zero
    
        ! 4: 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T 
        if (GETDERDI(4) .eq. 1) then
            gl%di_PCSAFT(i,4) = -0.36d0*gl%sigPCSAFT(i)*gl%epskPCSAFT(i)/T*dexp(-3.d0*gl%epskPCSAFT(i)/T)
        end if
    
        ! 5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2 
        if (GETDERDI(5) .eq. 1) then
            gl%di_PCSAFT(i,5) = (0.72d0 - 1.08d0*gl%epskPCSAFT(i)/T)*dexp(-3.d0*gl%epskPCSAFT(i)/T)* &
         & gl%epskPCSAFT(i)*gl%sigPCSAFT(i)/T
        end if
    
        ! 6: all rho derivatives of d_i equal zero
        ! 7: all rho derivatives of d_i equal zero
        ! 8: all rho derivatives of d_i equal zero
    
        ! 9: 3RD DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^3 
        if (GETDERDI(9) .eq. 1) then
            gl%di_PCSAFT(i,9) = (-2.16d0 + 6.48d0*gl%epskPCSAFT(i)/T - 3.24d0*gl%epskPCSAFT(i)**2/T**2 &
         & )*dexp(-3.d0*gl%epskPCSAFT(i)/T)*gl%epskPCSAFT(i)*gl%sigPCSAFT(i)/T
        end if
    
        ! 10: all rho derivatives of d_i equal zero    
        ! 11: all rho derivatives of d_i equal zero
        ! 12: all rho derivatives of d_i equal zero
        ! 13: all rho derivatives of d_i equal zero
        ! 14: all rho derivatives of d_i equal zero
    
        ! 15: 4TH DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^4
        if (GETDERDI(15) .eq. 1) then
            gl%di_PCSAFT(i,15) = (8.64d0 - 38.88d0*gl%epskPCSAFT(i)/T + 38.88d0*gl%epskPCSAFT( &
                & i)**2/T**2 - 9.72d0*gl%epskPCSAFT(i)**3/T**3)*dexp(-3.0d0 &
                & *gl%epskPCSAFT(i)/T)*gl%epskPCSAFT(i)*gl%sigPCSAFT(i)/T
        end if
    end do
    
!DEC$ END IF
end subroutine DIDERIVS


subroutine AQQDERIVS(gl,T,DENS,GETDERAQQ)
    
    ! a_qq: contribution to the residual helmholz engergy from quadrupol-quadrupol interactions 
    ! defined by eq. A.8 in Gross, Sadowski 2005:
    ! a_qq = a_2qq*a_2qq/(a_2qq-a_3qq)
    ! dependent on D and T





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAQQ
    double precision, intent (in) :: T, DENS
    !output: aqq_PCSAFT (module variable)
    !working variables

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !calculate the derivatives of a_qq
    if (GETDERAQQ(1) .eq. 1) then
        gl%AQQ_PCSAFTQ(1) = gl%A2_PCSAFTQ(1)/(1.d0-(gl%A3_PCSAFTQ(1)/gl%A2_PCSAFTQ(1)))
    end if
    
    !  2: 1ST DERIVATIVE OF a_qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(2) .eq. 1) then
        gl%AQQ_PCSAFTQ(2) = (2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(2)*(gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1))-gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*(gl%A2_PCSAFTQ(2)-gl%A3_PCSAFTQ(2))) &
            & /((gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1))*(gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)))
    end if
    ! I. Schuelling 05/17
    ! 3: 2ND DERIVATIVE OF a_qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(3) .eq. 1) then
    
        gl%AQQ_PCSAFTQ(3) = ( gl%A2_PCSAFTQ(1) * gl%A2_PCSAFTQ(1) * (-(gl%A2_PCSAFTQ(3)-gl%A3_PCSAFTQ(3))) + 2.d0 * gl%A2_PCSAFTQ(1) * gl%A2_PCSAFTQ(3) * (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) &
            & + 2.d0 * gl%A2_PCSAFTQ(2) * gl%A2_PCSAFTQ(2) * (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1))) &
            & / ((gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) * (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1))) &
            & - ( 2.d0 * (2.d0 * gl%A2_PCSAFTQ(1) * gl%A2_PCSAFTQ(2)* (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) - gl%A2_PCSAFTQ(1) * gl%A2_PCSAFTQ(1)* (gl%A2_PCSAFTQ(2)-gl%A3_PCSAFTQ(2)))* (gl%A2_PCSAFTQ(2)-gl%A3_PCSAFTQ(2)) ) &
            & / ((gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) * (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) *(gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)))
     end if
     
    ! 4: 1ST DERIVATIVE OF a_qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(4) .eq. 1) then
         gl%AQQ_PCSAFTQ(4) = (2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*(gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1))-gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*(gl%A2_PCSAFTQ(4)-gl%A3_PCSAFTQ(4))) &
            & /((gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1))*(gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)))
    end if
   
    ! 5: 2ND DERIVATIVE OF a_qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires a_2qq, a_3qq
     if (GETDERAQQ(5) .eq. 1) then
        gl%AQQ_PCSAFTQ(5) = ( gl%A2_PCSAFTQ(1) * gl%A2_PCSAFTQ(1) * (-(gl%A2_PCSAFTQ(5)-gl%A3_PCSAFTQ(5))) + 2.d0 * gl%A2_PCSAFTQ(1) * gl%A2_PCSAFTQ(5) * (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) &
            & + 2.d0 * gl%A2_PCSAFTQ(4) * gl%A2_PCSAFTQ(4) * (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1))) / ((gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) * (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1))) &
            & - ( 2.d0 * (2.d0 * gl%A2_PCSAFTQ(1) * gl%A2_PCSAFTQ(4)* (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) - gl%A2_PCSAFTQ(1) * gl%A2_PCSAFTQ(1)* (gl%A2_PCSAFTQ(4)-gl%A3_PCSAFTQ(4)))* (gl%A2_PCSAFTQ(4)-gl%A3_PCSAFTQ(4)) ) &
            & / ((gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) * (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) *(gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)))
     end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(6) .eq. 1) then
        gl%AQQ_PCSAFTQ(6) = (gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(6) + gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A3_PCSAFTQ(6) &
            & - 3.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(6)*gl%A3_PCSAFTQ(1) - gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) &
            & * gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(6) + 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(2) - 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(2) &
            & + 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(6)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) - 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4) + 2.d0*gl%A2_PCSAFTQ(4) &
            & * gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1))/(gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) - 3.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) &
            & + 3.d0*gl%A2_PCSAFTQ(1)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) - gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1))
    end if

    ! 7: 3RD MIXED DERIVATIVE OF a_hc WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires a_2qq, a_3qq
    
    if (GETDERAQQ(7) .eq. 1) then
        gl%AQQ_PCSAFTQ(7) = (gl%A2_PCSAFTQ(1)**4*gl%A2_PCSAFTQ(7) + gl%A2_PCSAFTQ(1)**4*gl%A3_PCSAFTQ(7) - 4.d0*gl%A2_PCSAFTQ(1)**3*gl%A2_PCSAFTQ(7)*gl%A3_PCSAFTQ(1) &
            & - 2.d0*gl%A2_PCSAFTQ(1)**3 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(7) + 4.d0*gl%A2_PCSAFTQ(1)**3*gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(6) + 2.d0*gl%A2_PCSAFTQ(1)**3*gl%A3_PCSAFTQ(5)*gl%A3_PCSAFTQ(2) &
            & - 4.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(6) - 4.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(2) &
            & - 2.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(5)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(2) + 5.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(7)*gl%A3_PCSAFTQ(1)**2 &
            & - 4.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(6)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4) - 2.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(5) &
            & - 2.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(4)**2 + gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(7) - 4.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(6)&
            & - 2.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(5)*gl%A3_PCSAFTQ(2) + 6.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(4)**2 *gl%A3_PCSAFTQ(2) &
            & + 4.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)**2 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(2) + 4.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(6)*gl%A3_PCSAFTQ(1)**2 &
            & + 8.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4) + 4.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(6) &
            & - 8.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(2) + 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(5)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)**2 &
            & + 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(5)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(2) - 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(7)*gl%A3_PCSAFTQ(1)**3 &
            & + 4.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(6)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(4) + 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(5) &
            & - 4.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4)**2 - 6.d0*gl%A2_PCSAFTQ(4)**2 *gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)**2 + 2.d0*gl%A2_PCSAFTQ(4)**2*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(2) &
            & - 4.d0*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(6)*gl%A3_PCSAFTQ(1)**3 + 4.d0*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(4) - 2.d0*gl%A2_PCSAFTQ(5)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)**3) &
            & /(gl%A2_PCSAFTQ(1)**4 - 4.d0*gl%A2_PCSAFTQ(1)**3 *gl%A3_PCSAFTQ(1) + 6.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)**2 - 4.d0*gl%A2_PCSAFTQ(1)*gl%A3_PCSAFTQ(1)**3 + gl%A3_PCSAFTQ(1)**4)
 
    end if
    
    ! I. Schuelling 05/17
    ! 8: 3ND DERIVATIVE OF a_qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(8) .eq. 1) then
        gl%AQQ_PCSAFTQ(8) = 1.d0/((gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) * (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) *(gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)) * (gl%A2_PCSAFTQ(1)-gl%A3_PCSAFTQ(1)))&
            & * ( - 4.d0 * gl%A2_PCSAFTQ(8) * gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1) &
            & + 5.d0 *gl%A2_PCSAFTQ(8) * gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) &
            & - 2.d0 * gl%A2_PCSAFTQ(8) * gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) &
            & + gl%A2_PCSAFTQ(8) * gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) &
            & - 6.d0 * gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) *gl%A3_PCSAFTQ(1) * gl%A2_PCSAFTQ(3) * gl%A3_PCSAFTQ(2) &
            & + 6.d0 * gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) * gl%A2_PCSAFTQ(3) * gl%A3_PCSAFTQ(2) &
            & - 6.d0 * gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1) * gl%A2_PCSAFTQ(2) * gl%A3_PCSAFTQ(3) &
            & + 6.d0 * gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) * gl%A2_PCSAFTQ(2) * gl%A3_PCSAFTQ(3) &
            & - 6.d0 * gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) * gl%A2_PCSAFTQ(2) * gl%A3_PCSAFTQ(2)*gl%A3_PCSAFTQ(2) &
            & - 12.d0 * gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1) * gl%A2_PCSAFTQ(2) * gl%A3_PCSAFTQ(2)*gl%A3_PCSAFTQ(2) &
            & + 12.d0 * gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1) * gl%A2_PCSAFTQ(2)*gl%A2_PCSAFTQ(2) * gl%A3_PCSAFTQ(2) &
            & + 6.d0 * gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) * gl%A2_PCSAFTQ(2)*gl%A2_PCSAFTQ(2) * gl%A3_PCSAFTQ(2) &
            & - 6.d0 * gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) * gl%A2_PCSAFTQ(2)*gl%A2_PCSAFTQ(2)*gl%A2_PCSAFTQ(2) &
            & + 6.d0 * gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1)* gl%A3_PCSAFTQ(1) *gl%A2_PCSAFTQ(2) * gl%A2_PCSAFTQ(3) &
            & - 6.d0 * gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) * gl%A2_PCSAFTQ(2) * gl%A2_PCSAFTQ(3) &
            & + gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(8) &
            & - 2.d0 * gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1) * gl%A3_PCSAFTQ(8) &
            & + gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(1) * gl%A3_PCSAFTQ(8) &
            & + 6.d0 * gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(2)*gl%A3_PCSAFTQ(2)*gl%A3_PCSAFTQ(2) &
            & + 6.d0 * gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(2) * gl%A3_PCSAFTQ(3) &
            & - 6.d0 * gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(1) * gl%A3_PCSAFTQ(1) * gl%A3_PCSAFTQ(2) * gl%A3_PCSAFTQ(3) ) 
    end if   
    
    !9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    if (GETDERAQQ(10) .eq. 1) then
        gl%AQQ_PCSAFTQ(9) = (gl%A2_PCSAFTQ(1)**4 *gl%A2_PCSAFTQ(9) + gl%A2_PCSAFTQ(1)**4 *gl%A3_PCSAFTQ(9) - 4.d0*gl%A2_PCSAFTQ(1)**3 *gl%A2_PCSAFTQ(9)*gl%A3_PCSAFTQ(1) &
            & - 2.d0*gl%A2_PCSAFTQ(1)**3 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(9) + 6.d0*gl%A2_PCSAFTQ(1)**3 *gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(5) - 6.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(5) &
            & - 6.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(4)**2 - 6.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(5)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4) &
            & + 5.d0*gl%A2_PCSAFTQ(1)**2* gl%A2_PCSAFTQ(9)*gl%A3_PCSAFTQ(1)**2 + gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(9) - 6.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(5) &
            & + 6.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(4)**3 + 12.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)**2 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4) + 6.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(5)*gl%A3_PCSAFTQ(1)**2 &
            & + 6.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(5) - 12.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4)**2 &
            & + 6.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(5)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(4) - 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(9)*gl%A3_PCSAFTQ(1)**3 - 6.d0*gl%A2_PCSAFTQ(4)**3 *gl%A3_PCSAFTQ(1)**2 &
            & + 6.d0*gl%A2_PCSAFTQ(4)**2 *gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(4) - 6.d0*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(5)*gl%A3_PCSAFTQ(1)**3)/(gl%A2_PCSAFTQ(1)**4 &
            & - 4.d0*gl%A2_PCSAFTQ(1)**3 *gl%A3_PCSAFTQ(1) + 6.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)**2 - 4.d0*gl%A2_PCSAFTQ(1)*gl%A3_PCSAFTQ(1)**3 + gl%A3_PCSAFTQ(1)**4)
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF AQQ WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires a_2qq, a_3qq
    if (GETDERAQQ(10) .eq. 1) then
        gl%AQQ_PCSAFTQ(10) = (gl%A2_PCSAFTQ(1)**4 *gl%A2_PCSAFTQ(10) + gl%A2_PCSAFTQ(1)**4 *gl%A3_PCSAFTQ(10) - 4.d0*gl%A2_PCSAFTQ(1)**3*gl%A2_PCSAFTQ(10)*gl%A3_PCSAFTQ(1) &
            & - 2.d0*gl%A2_PCSAFTQ(1)**3 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(10) + 2.d0*gl%A2_PCSAFTQ(1)**3 *gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(3) + 4.d0*gl%A2_PCSAFTQ(1)**3 *gl%A3_PCSAFTQ(6)*gl%A3_PCSAFTQ(2) &
            & - 2.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(3) - 2.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(2)**2 &
            & - 4.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(6)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(2) + 5.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(10)*gl%A3_PCSAFTQ(1)**2 &
            & - 4.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(6) - 4.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(2) &
            & - 2.d0*gl%A2_PCSAFTQ(1)**2 *gl%A2_PCSAFTQ(3)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4) + gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(10) &
            & - 2.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(3) - 4.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(6)*gl%A3_PCSAFTQ(2) &
            & + 6.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(2)**2 + 8.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(2) &
            & + 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(3)*gl%A3_PCSAFTQ(1)**2 + 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(3) &
            & - 4.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(4)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(2)**2 + 4.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(6)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)**2 &
            & + 4.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(6)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(2) - 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(10)*gl%A3_PCSAFTQ(1)**3 &
            & + 4.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(2)**2 *gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4) + 4.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(6) &
            & - 8.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)*gl%A3_PCSAFTQ(4)*gl%A3_PCSAFTQ(2) + 2.d0*gl%A2_PCSAFTQ(1)*gl%A2_PCSAFTQ(3)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(4) &
            & - 6.d0*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(2)**2 *gl%A3_PCSAFTQ(1)**2 + 4.d0*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(2) &
            & - 2.d0*gl%A2_PCSAFTQ(4)*gl%A2_PCSAFTQ(3)*gl%A3_PCSAFTQ(1)**3 - 4.d0*gl%A2_PCSAFTQ(6)*gl%A2_PCSAFTQ(2)*gl%A3_PCSAFTQ(1)**3 + 2.d0*gl%A2_PCSAFTQ(2)**2 *gl%A3_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(4)) &
            & /(gl%A2_PCSAFTQ(1)**4 - 4.d0*gl%A2_PCSAFTQ(1)**3 *gl%A3_PCSAFTQ(1) + 6.d0*gl%A2_PCSAFTQ(1)**2 *gl%A3_PCSAFTQ(1)**2 - 4.d0*gl%A2_PCSAFTQ(1)*gl%A3_PCSAFTQ(1)**3 + gl%A3_PCSAFTQ(1)**4)
    
    end if

!DEC$ END IF
end subroutine AQQDERIVS

subroutine J2DERIVS(gl,T, DENS, GETDERI)

! T.Eckermann, 02/2017
    
    ! J2: integral of the perturbation theory
    ! defined by eq. eq 12 in Gross 2005:
    ! J2 = Sum_0_4(a_nij+b_nij*eps/T)*eta**n
    ! dependent on D and T
    ! derivative stored in module variable





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (inout) :: GETDERI
    double precision, intent (in) :: T, DENS
    !output: j2_PCSAFT (module variables)
    integer :: i,j,n
    double precision :: eps_ij, sum1, sum2, sum3, sum4, sum5
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    eps_ij = 0.d0
    gl%J2_PCSAFTQ = 0.d0
    call add_T_conversion_derivs(gl,GETDERI)
    
    !calculate the derivatives of J_2
    ! 1: J_index
    if (GETDERI(1) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTQ(1,i,j) = gl%J2_PCSAFTQ(1,i,j) + (gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T)*gl%z3_PCSAFT(1)**n  
                end do
            end do
        end do
    end if
    
    
    !  2: 1ST DERIVATIVE OF J_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                gl%J2_PCSAFTQ(2,i,j) = gl%J2_PCSAFTQ(2,i,j) + n*gl%z3_PCSAFT(1)**n*(gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T)
                end do
            end do
        end do
    end if
    
    
    ! 3: 2ND DERIVATIVE OF J_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTQ(3,i,j) = gl%J2_PCSAFTQ(3,i,j) + n*(n-1)*gl%z3_PCSAFT(1)**n*(gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T)
                 end do
            end do
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
          do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTQ(4,i,j) = gl%J2_PCSAFTQ(4,i,j) - gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T*gl%z3_PCSAFT(1)**n + n*gl%z3_PCSAFT(1)**(n - 1)*gl%z3_PCSAFT(4)*(gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T)
                end do
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
          do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTQ(5,i,j) = gl%J2_PCSAFTQ(5,i,j) + (gl%z3_PCSAFT(1)**(n - 2)*(T**2 *n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                                & + T**2 *n*gl%z3_PCSAFT(4)**2 *(n - 1)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                                & - 2.d0*T**2 *gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + T**2 *2.d0*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2)/T**3)
                end do
                
            end do
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF J2 WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    if (GETDERI(6) .eq. 1) then
     do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTQ(6,i,j) = gl%J2_PCSAFTQ(6,i,j) + (n**2 *gl%z3_PCSAFT(1)**(n - 2) *gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*(gl%ab_PCSAFTQ(1,n,i,j)+ gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T) &
                    & - n*gl%z3_PCSAFT(1)**(n - 2) *gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*(gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T) &
                    & + n*gl%z3_PCSAFT(1)**(n - 1) *gl%z3_PCSAFT(6)*(gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T) &
                    & - gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T*n*gl%z3_PCSAFT(1)**(n - 1) *gl%z3_PCSAFT(2))
                 end do
            end do
        end do   
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF J2 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires
    if (GETDERI(7) .eq. 1) then 
    do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTQ(7,i,j) = gl%J2_PCSAFTQ(7,i,j) + (n*gl%z3_PCSAFT(1)**(n - 3) *(T*(T*gl%z3_PCSAFT(1)**2 *gl%z3_PCSAFT(7)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                     & + T*gl%z3_PCSAFT(1)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(2.d0*n*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(6) - 2.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(6) - gl%z3_PCSAFT(5)*gl%z3_PCSAFT(2)) &
                     & + 2.d0*T*gl%z3_PCSAFT(4)**2 *gl%z3_PCSAFT(2)*(-n + 1)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) & 
                     & - 2.d0*T*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2 *gl%z3_PCSAFT(6) + 2.d0*T*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)) & 
                     & + gl%z3_PCSAFT(2)*(T**2 *n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                     & + T**2 *n*gl%z3_PCSAFT(4)**2 *(n - 1)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                     & - 2.d0*T**2 *gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + 2.d0*T**2 *gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2))/T**3)
                end do
            end do
        end do
    end if
    
    ! 8: 3RD DERIVATIVE OF J_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    if (GETDERI(8) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTQ(8,i,j) = gl%J2_PCSAFTQ(8,i,j) + n*gl%z3_PCSAFT(1)**n*(gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij/T)&
                        & *(n**2 *gl%z3_PCSAFT(1) - 3.d0*n*gl%z3_PCSAFT(1) + 2.d0*gl%z3_PCSAFT(1))/gl%z3_PCSAFT(1)
                        
                 end do
            end do
        end do
    end if
    
    !! 9: 3RD DERIVATIVE OF J_2 WITH RESPECT TO D AT CONSTANT D, MULTIPLIED BY T^3
    if (GETDERI(9) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4                   
                    gl%J2_PCSAFTQ(9,i,j) = gl%J2_PCSAFTQ(9,i,j) + (gl%z3_PCSAFT(1)**(n - 3) *(T**3 *n*gl%z3_PCSAFT(1)**2 *gl%z3_PCSAFT(9)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                        & + 3.d0*T**3 *n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(n - 1)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                        & + T**3 *n*gl%z3_PCSAFT(4)**3 *(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(n**2 - 3*n + 2) &
                        & - 3.d0*T**3 *gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)**2 *gl%z3_PCSAFT(5) + 3.d0*T**3 *gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(-n + 1) &
                        & + 6.d0*T**3 *gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)**2 *gl%z3_PCSAFT(4) - 6.d0*T**3*gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**3)/T**4)
                 end do
            end do
        end do
    end if
    
    !10: 3RD MIXED DERIVATIVE OF J2 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    if (GETDERI(10) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTQ(10,i,j) = gl%J2_PCSAFTQ(10,i,j) + (n*gl%z3_PCSAFT(1)**(n - 3.d0) *gl%z3_PCSAFT(2)*(2.d0*T*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(6)*(n - 1)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij) &
                        & + T*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*(T*gl%ab_PCSAFTQ(1,n,i,j) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij)*(n**2 - 3*n + 2) + gl%ab_PCSAFTQ(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(2)*T*(-n + 1))/T**2)
                 end do
            end do
        end do
    end if
 

!DEC$ END IF
end subroutine J2DERIVS    

subroutine J3DERIVS(gl,T,DENS,GETDERI)

! T.Eckermann, 02/2017
    
    ! J3: integral of the perturbation theory
    ! defined by eq. eq 13 in Gross 2005:
    ! J3 = Sum_0_4c_nijk*eta**n
    ! dependent on D and T
    




implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    double precision, intent (in) :: DENS, T
    !output: j3_PCSAFT (module variables)
    integer :: i,j,k,n
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    gl%J3_PCSAFTQ = 0.d0
        
    !calculate the derivatives of J_3
    ! 1: J_index 
    if (GETDERI(1) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTQ(1,i,j,k) = gl%J3_PCSAFTQ(1,i,j,k) + gl%c_PCSAFTQ(n,i,j,k)*gl%z3_PCSAFT(1)**n
                    end do
                end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF J_3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTQ(2,i,j,k) = gl%J3_PCSAFTQ(2,i,j,k) + gl%c_PCSAFTQ(n,i,j,k)*n*gl%z3_PCSAFT(1)**n
                    end do
                end do
            end do
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTQ(3,i,j,k) = gl%J3_PCSAFTQ(3,i,j,k) + n*(n-1) * gl%c_PCSAFTQ(n,i,j,k)*gl%z3_PCSAFT(1)**(n-2)*gl%z3_PCSAFT(2)*gl%z3_PCSAFT(2)
                    end do
                end do
            end do
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTQ(4,i,j,k) = gl%J3_PCSAFTQ(4,i,j,k) + gl%c_PCSAFTQ(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 1)*gl%z3_PCSAFT(4)
                    end do
                end do
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTQ(5,i,j,k) = gl%J3_PCSAFTQ(5,i,j,k) + gl%c_PCSAFTQ(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 2)*(n*gl%z3_PCSAFT(4)**2 &
                            & + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - gl%z3_PCSAFT(4)**2)
                    end do
                end do
            end do
        end do    
    end if
    
     ! 6: 2ND MIXED DERIVATIVE OF J3 WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    if (GETDERI(6) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTQ(6,i,j,k) = gl%J3_PCSAFTQ(6,i,j,k) + gl%c_PCSAFTQ(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 2)*(gl%z3_PCSAFT(1)*gl%z3_PCSAFT(6) &
                            & + gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*(n - 1))
                    end do
                end do
            end do
        end do   
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF J3 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires 
    if (GETDERI(7) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTQ(7,i,j,k) = gl%J3_PCSAFTQ(7,i,j,k) + gl%c_PCSAFTQ(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 3)*(gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(7) + 2* &
                         & gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(6)*(n - 1) + gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(2)*(-n &
                         & + 1) + gl%z3_PCSAFT(2)*(n - 1)*(n*gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - &
                         & gl%z3_PCSAFT(4)**2))
                    end do
                end do
            end do
        end do   
    end if
    
    ! 8: 3RD DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    if (GETDERI(8) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTQ(8,i,j,k) = gl%J3_PCSAFTQ(8,i,j,k) + n*(n-1)*(n-2) * gl%c_PCSAFTQ(n,i,j,k)*gl%z3_PCSAFT(1)**(n-3)*gl%z3_PCSAFT(2)*gl%z3_PCSAFT(2)*gl%z3_PCSAFT(2)
                    end do
                end do
            end do
        end do
    end if
    
    ! 9: 3RD DERIVATIVE OF J3 WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires
    if (GETDERI(9) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTQ(9,i,j,k) = gl%J3_PCSAFTQ(9,i,j,k) + gl%c_PCSAFTQ(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 3)*(gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9) + 3*gl%z3_PCSAFT(1)* &
                            & gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(n - 1) + gl%z3_PCSAFT(4)**3*(n**2 - 3*n + 2))
                    end do
                end do
            end do
        end do   
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF J3 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires 
    if (GETDERI(10) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTQ(10,i,j,k) = gl%J3_PCSAFTQ(10,i,j,k) + gl%c_PCSAFTQ(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 3)*gl%z3_PCSAFT(2)*(2*gl%z3_PCSAFT(1)* &
                            & gl%z3_PCSAFT(6)*(n - 1) + gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*(n**2 - 3*n + 2))
                    end do
                end do
            end do
        end do   
    end if
        
!DEC$ END IF
end subroutine J3DERIVS    

subroutine A2DERIVS(gl,T,DENS,GETDERAQQ)
    ! a_2qq: second-order perturbation term of the quadrupolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.9 in Gross, Sadowski 2005:
    ! a_2qq =-pi*(3/4)**2*rho*sum_j*sum_i(x_i*x_j*epsk_ii/T*espk_jj/T*sig_ii**5*sig_jj**5/sig_ij**7*n_Qi*n_Qj*Q_i**2*Q_j**2*J_2ij
    ! dependent on D and T







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAQQ
    double precision, intent (in) :: T, DENS
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sum, sum3, sum1,sum1a, sum1b, sum2, sum2a,sum2b, sum4a, sum4b, sum5a, sum5b, sum5c, sum6, sum6b, sum7c,sum8, sum9,sum10, sigma_ij, factor_lc
    integer :: i, j
    double precision, dimension(2) :: qq
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !Initializations
    sum3 = 0.d0
    sum1 = 0.d0
    sum2 = 0.d0
    
    sigma_ij = 0.d0
    factor_lc = -piPCSAFT*0.5625d0*dens
    

    !calculate the derivatives of a_qq
    ! 1: a_2qq
    if (GETDERAQQ(1) .eq. 1) then
            sum1 = 0.d0
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                    sum1 = sum1 + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7  * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                        & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(1,i,j)
                end do
            end do
        
            gl%A2_PCSAFTQ(1) = factor_lc * sum1 *1.d-38 ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19
    end if
    
    !  2: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires j_2qq
    if (GETDERAQQ(2) .eq. 1) then
        sum2 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                sum2 = sum2 + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(2,i,j)
            end do
        end do
    
        gl%A2_PCSAFTQ(2) = gl%A2_PCSAFTQ(1) + sum2*factor_lc*1.d-38 ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19

    end if
    
   ! I. Schuelling 05/17:
   ! 3: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
   ! requires j_2qq
    if (GETDERAQQ(3) .eq. 1) then
        sum3 = 0.d0
        sum2 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                
                sum2 = sum2 + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(2,i,j)
                
                sum3 = sum3 + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(3,i,j)

            end do
        end do
        gl%A2_PCSAFTQ(3) = 2.d0*sum2*factor_lc*1.d-38 + sum3*factor_lc*1.d-38 
        !A2_PCSAFTQ(3) = 2.d0*sum2(incl. *DENS)*factor/DENS *DENS*1.d-38 + sum3(incl. DEL**2 from J2)*factor*1.d-38 
    end if
 
    ! 4: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires j_2qq
     if (GETDERAQQ(4) .eq. 1) then
        sum4a = 0.d0
        sum4b = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                sum4a = sum4a - 2.d0 * gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(1,i,j)
                
                sum4b = sum4b + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(4,i,j)
            end do
        end do
        
        gl%A2_PCSAFTQ(4) = factor_lc * sum4a *1.d-38 + factor_lc * sum4b *1.d-38 ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19
     end if
    
    ! 5: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
     !requires j_2qq
     if (GETDERAQQ(5) .eq. 1) then
        sum5a = 0.d0
        sum5b = 0.d0
        sum5c = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                
                sum5a = sum5a + 6.d0 * gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(1,i,j)
                
                sum5b = sum5b - 2.d0 * gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(4,i,j)
               
                sum5c = sum5c + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(5,i,j)
            end do
        end do
        !A2_PCSAFTQ(5) = factor * sum5a *1.d-38 + 2.d0 * factor * sum5b *1.d-38 + factor * sum5c *1.d-38! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19
        gl%A2_PCSAFTQ(5) = factor_lc * (sum5a + 2.d0*sum5b + sum5c) *1.d-38 ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19 !shortend the term
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_2qq
    if (GETDERAQQ(6) .eq. 1) then
     sum1a = 0.d0
     sum1b = 0.d0
     sum2a = 0.d0
     sum2b = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                sum1a = sum1a + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(4,i,j)
                sum1b = sum1b - 2.d0 * gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(1,i,j)
                sum2a = sum2a + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) & 
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(6,i,j)
                sum2b = sum2b - 2.d0 * gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(2,i,j)
           
            end do
        end do
     
        gl%A2_PCSAFTQ(6) = (sum1a+sum1b+sum2a+sum2b) *factor_lc*1.d-38 ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
    end if
    
    ! Schuelling 06/17:
    ! 7: 3RD MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires j_2qq
    if (GETDERAQQ(7) .eq. 1) then
        sum5a = 0.d0
        sum5b = 0.d0
        sum5c = 0.d0
        sum2a = 0.d0
        sum6b = 0.d0
        sum7c = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                
                sum5a = sum5a + 6.d0 * gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(1,i,j)
                
                sum5b = sum5b - 2.d0 * gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(4,i,j)
               
                sum5c = sum5c + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(5,i,j)
           
                 sum2a = sum2a + 6.d0 * gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(2,i,j)
                
                sum6b = sum6b - 2.d0 * gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(6,i,j)
               
                sum7c = sum7c + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(7,i,j)
    
            end do
        end do
        !A2_PCSAFTQ(5) = factor * sum5a *1.d-38 + 2.d0 * factor * sum5b *1.d-38 + factor * sum5c *1.d-38! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19
        gl%A2_PCSAFTQ(7) = (sum5a + sum2a + 2.d0*sum5b + 2.d0*sum6b + sum5c + sum7c)*factor_lc*1.d-38  ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19 !shortend the term 
    end if 
    
    ! I. Schuelling 05/17:
    ! 8: 3RD DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires j2_qq
    if (GETDERAQQ(8) .eq. 1) then
        sum3 = 0.d0
        sum2 = 0.d0
        sum8 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                
                sum2 = sum2 + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(2,i,j)
                
                sum3 = sum3 + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(3,i,j)
               
                sum8 = sum8 + gl%molfractions(i)*gl%molfractions(j) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%sigPCSAFT(i)**5 * gl%sigPCSAFT(j)**5 / sigma_ij**7 * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) &
                    & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%J2_PCSAFTQ(8,i,j)

            end do
        end do
        gl%A2_PCSAFTQ(8) = 3.d0*sum3*factor_lc*1.d-38 + sum8*factor_lc*1.d-38  
    end if
   
    ! 9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires 
    if (GETDERAQQ(9) .eq. 1) then
    sum1 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0                     
                     sum1 = sum1 + gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) *gl%molfractions(i)*gl%molfractions(j) / (kbol**2 *T**2 *sigma_ij**7 *gl%mPCSAFT(i)*gl%mPCSAFT(j)) &
                     & *(-24.d0*gl%J2_PCSAFTQ(1,i,j) + 18.d0*gl%J2_PCSAFTQ(4,i,j) - 6.d0*gl%J2_PCSAFTQ(5,i,j) + gl%J2_PCSAFTQ(9,i,j))
            end do
        end do
        gl%A2_PCSAFTQ(9) = factor_lc*sum1*1.d-38 
    end if
    
    !Schuelling: 06/17
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires j_2qq
    if (GETDERAQQ(10) .eq. 1) then
            sum1 = 0.d0
            sum2 = 0.d0
            sum = 0.d0
            
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                    
                    sum = gl%QPCSAFTQ(i)**2*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%molfractions(i)*gl%molfractions(j)/(T**2*kbol**2*sigma_ij**7*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                    
                    sum1 = sum1 + (gl%J2_PCSAFTQ(10,i,j) - 2.d0 * gl%J2_PCSAFTQ(3,i,j)) * sum
                   
                    sum2 = sum2 + (gl%J2_PCSAFTQ(6,i,j) - 2.d0 * gl%J2_PCSAFTQ(2,i,j)) * sum
                end do
            end do
            gl%A2_PCSAFTQ(10) = factor_lc*1.d-38 * (Sum1 + 2.d0 * Sum2)
    end if 
    
!DEC$ END IF
end subroutine A2DERIVS
   
subroutine A3DERIVS (gl,T, DENS, GETDERAQQ)

    ! a_3qq: thrid-order perturbation term of the quadrupolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.10 in Gross, Sadowski 2005:
    ! a_3qq = pi/3*(3/4)**3*rho*sum_i*sum_j*x_i*x_j*(espk_ii/T)**1.5*(espk_jj/T)**1.5*sig_ii**(15/2)*sig_jj**(15/2)/sig_ij**12*n_Qi*n_Qj*Q_i**3*Q_j**3*J_3ij + &       
    ! & 4 *pi**2/3*(3/4)**3*rho**2*sum_j*sum_i*sum_k(x_i*x_j*x_k*epsk_ii/T*espk_jj/T*espk_kk/T*sig_ii**5*sig_jj**5*sig_kk**5/sig_ij**3/sig_ik**3/sig_jk**3*n_Qi*n_Qj*n_Qk*Q_i**2*Q_j**2*Q_k**2*J_3ijk
    ! J3,ij = 0, which means that the first summation equals zero and is therefore left out
    ! dependent on D and T 







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAQQ
    double precision, intent (in) :: T, DENS
    !working variables
    double precision :: sum1, sum1a, sum1b, sum2, sum2a, sum2b, sum3, sum3a, sum3b, sum4, sum4a, sum4b, sum5, sum5a, sum5b, sum5c, sum6, sum6b, sum7, sum7c,sum8, sum9,sum10, factor2
    double precision, dimension(gl%ncomp,gl%ncomp) :: sigma_ij
    integer :: i, j, k
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !Initializations
    sigma_ij = 0.d0
    factor2 = piPCSAFT*piPCSAFT*0.5625d0*dens*dens * 1.d-57  ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19
    do j = gl%n_start, gl%n_end
        do i = gl%n_start, gl%n_end
            sigma_ij(i,j) = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
        end do
    end do
    
    !eq. 10 in Gross (2005). J3,ij = 0, which means that the first summation equals zero and is therefore left out
    !calculate the derivatives of a_3qq
    if (GETDERAQQ(1) .eq. 1) then
        sum1= 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                        & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                        & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                        & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(1,i,j,k)
                end do
            end do
        end do
        gl%A3_PCSAFTQ(1) = factor2 * sum1
    end if
    
    !  2: 1ST DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires J_3ijk
    if (GETDERAQQ(2) .eq. 1) then
        sum2 = 0.d0
        sum1 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                        & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                        & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                        & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(1,i,j,k)
                    sum2 = sum2 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                        & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                        & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                        & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(2,i,j,k)
                    
            end do
        end do
        end do
        
        gl%A3_PCSAFTQ(2) = sum1*factor2*2.d0 + factor2 * sum2
        !A3_PCSAFTQ(2) = A3_PCSAFTQ(1)*2.d0/DENS + factor2 * sum2 ! /DENS raus, da die Gleichung noch mit DENS multipliziert werden muss.In sum2 steckt zeta3 und in ZETADERIVS wurde die Gleichung schon mit dens multipliziert

    end if
        
    ! 3: 2ND DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires j_3ijk
    if (GETDERAQQ(3) .eq. 1) then
       sum2 = 0.d0
       sum1 = 0.d0
       sum3 = 0.d0
       do j = gl%n_start, gl%n_end
           do i = gl%n_start, gl%n_end
               do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                        & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5/sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                        & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                        & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(1,i,j,k)
    
                    sum2 = sum2 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                       & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                       & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                       & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(2,i,j,k)
                   
                   sum3 = sum3 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                       & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                       & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                       & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(3,i,j,k)
                end do
           end do
       end do
        gl%A3_PCSAFTQ(3) = 2.d0*sum1*factor2 + 4.d0*sum2*factor2 + factor2*sum3
        end if
    
    ! 4: 1ST DERIVATIVE OF a_3qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires j_3ijk
    if (GETDERAQQ(4) .eq. 1) then
        sum4= 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum4 = sum4 + (-3.d0*gl%J3_PCSAFTQ(1,i,j,k)*gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) &
                        & *gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)/(kbol**3 *T**3 *sigma_ij(i,j)**3 *sigma_ij(i,k)**3 *sigma_ij(j,k)** 3*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k)) &
                        & + gl%J3_PCSAFTQ(4,i,j,k)*gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) &
                        & /(kbol**3 *T**3 *sigma_ij(i,j)**3 *sigma_ij(i,k)**3 *sigma_ij(j,k)**3 *gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k)))
                end do
            end do
        end do
        gl%A3_PCSAFTQ(4) = factor2*sum4
    end if
    
    ! 5: 2ND DERIVATIVE OF a_3qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires j_3ijk
      if (GETDERAQQ(5) .eq. 1) then
        sum5= 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                   sum5 = sum5 + gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) &
                        & /(T**3 *kbol**3 *sigma_ij(i,j)**3 *sigma_ij(i,k)**3 *sigma_ij(j,k)**3 *gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k)) &
                        & *(12.d0*gl%J3_PCSAFTQ(1,i,j,k) - 6.d0*gl%J3_PCSAFTQ(4,i,j,k) + gl%J3_PCSAFTQ(5,i,j,k))
                end do
            end do
        end do
        gl%A3_PCSAFTQ(5) = factor2 * sum5
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_3 WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_3ijk
    if (GETDERAQQ(6) .eq. 1) then
      sum1 = 0.d0
      sum2 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                sum1 = sum1 + (gl%J3_PCSAFTQ(6,i,j,k)*gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) &
                    & /(T**3 *kbol**3 *sigma_ij(i,j)**3 *sigma_ij(i,k)**3 *sigma_ij(j,k)**3 *gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k)) &
                    & - 3.d0*gl%J3_PCSAFTQ(2,i,j,k)*gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k)* &
                    & gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)/(T**3 *kbol**3 *sigma_ij(i,j)**3 *sigma_ij(i,k)**3 *sigma_ij(j,k)**3 *gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k)))
                sum2 = sum2 + (-3.d0*gl%J3_PCSAFTQ(1,i,j,k)*gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) &
                    & /(T**3 *kbol**3 *sigma_ij(i,j)**3 *sigma_ij(i,k)**3 *sigma_ij(j,k)**3 *gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k)) &
                    & + gl%J3_PCSAFTQ(4,i,j,k)*gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) &
                    & /(T**3 *kbol**3 *sigma_ij(i,j)**3 *sigma_ij(i,k)**3 *sigma_ij(j,k)**3 *gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k)))
                   
                end do
            end do
        end do
        gl%A3_PCSAFTQ(6) = factor2*sum1 + 2.d0*factor2*sum2

    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_3qq WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires j_3ijk
    if (GETDERAQQ(7) .eq. 1) then
        sum1 = 0.d0
        sum2 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k)*(gl%J3_PCSAFTQ(7,i,j,k) &
                        & - 6.d0*gl%J3_PCSAFTQ(6,i,j,k) + 12.d0*gl%J3_PCSAFTQ(2,i,j,k))*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) &
                        & /(T**3 *kbol**3 *sigma_ij(i,j)**3 *sigma_ij(i,k)**3 *sigma_ij(j,k)**3 *gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                    sum2 = sum2 + gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) &
                        & *(12.d0*gl%J3_PCSAFTQ(1,i,j,k) - 6.d0*gl%J3_PCSAFTQ(4,i,j,k) + gl%J3_PCSAFTQ(5,i,j,k))*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) &
                        & /(T**3 *kbol**3 *sigma_ij(i,j)**3 *sigma_ij(i,k)**3 *sigma_ij(j,k)**3 *gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                end do
            end do
        end do
        gl%A3_PCSAFTQ(7) = factor2*sum1 + 2.d0*factor2*sum2
    end if
    
    ! 8: 3RD DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires j_3ijk
    if (GETDERAQQ(8) .eq. 1) then
       sum2 = 0.d0
       sum8 = 0.d0
       sum3 = 0.d0
       do j = gl%n_start, gl%n_end
           do i = gl%n_start, gl%n_end
               do k = gl%n_start, gl%n_end
                    sum8 = sum8 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                        & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3* gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) &
                        & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 &
                        & * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(8,i,j,k)
                   
                    sum2 = sum2 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                       & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                       & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                       & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(2,i,j,k)
                   
                   sum3 = sum3 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                       & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3* gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) &
                       & * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 &
                       & * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(3,i,j,k)
                end do
           end do
       end do
        gl%A3_PCSAFTQ(8) = 2.d0*sum2*factor2 + 4.d0*sum2*factor2 +6.d0*factor2*sum3 + factor2*sum8

    end if
   ! 9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires 
    if (GETDERAQQ(9) .eq. 1) then
        sum1 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) &
                         & *(-60.d0*gl%J3_PCSAFTQ(1,i,j,k) + 36.d0*gl%J3_PCSAFTQ(4,i,j,k) - 9.d0*gl%J3_PCSAFTQ(5,i,j,k) + gl%J3_PCSAFTQ(9,i,j,k))*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) &
                         & /(T**3 *kbol**3 *sigma_ij(i,j)**3 *sigma_ij(i,k)**3 *sigma_ij(j,k)**3*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                end do
            end do
        end do
       gl%A3_PCSAFTQ(9) = factor2*sum1
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires j_3ijk
    if (GETDERAQQ(10) .eq. 1) then
       sum2a = 0.d0
       sum1a = 0.d0
       sum3a = 0.d0
       sum2b = 0.d0
       sum1b = 0.d0
       sum3b = 0.d0
       do j = gl%n_start, gl%n_end
           do i = gl%n_start, gl%n_end
               do k = gl%n_start, gl%n_end
                    sum1a = sum1a + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                        & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                        & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                        & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(4,i,j,k)
                    sum1b = sum1b - 3.d0 * gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                        & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                        & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                        & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(1,i,j,k)
                   
                    sum2a = sum2a -3.d0 * gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                       & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                       & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                       & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(2,i,j,k)
                    sum2b = sum2b + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                       & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                       & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                       & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(6,i,j,k)
                   
                   sum3a = sum3a - 3.d0 * gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                       & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                       & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                       & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(3,i,j,k)
                   sum3b = sum3b + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%epskPCSAFT(i)/T * gl%epskPCSAFT(j)/T * gl%epskPCSAFT(k)/T &
                       & * gl%sigPCSAFT(i)**5 *gl%sigPCSAFT(j)**5 *gl%sigPCSAFT(k)**5 /sigma_ij(i,j)**3 /sigma_ij(i,k)**3 /sigma_ij(j,k)**3 &
                       & * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)*gl%QPCSAFTQ(i)/gl%mPCSAFT(i)/gl%epskPCSAFT(i)/kbol/gl%sigPCSAFT(i)**5 &
                       & * gl%QPCSAFTQ(j)*gl%QPCSAFTQ(j)/gl%mPCSAFT(j)/gl%epskPCSAFT(j)/kbol/gl%sigPCSAFT(j)**5 * gl%QPCSAFTQ(k)*gl%QPCSAFTQ(k)/gl%mPCSAFT(k)/gl%epskPCSAFT(k)/kbol/gl%sigPCSAFT(k)**5 * gl%J3_PCSAFTQ(10,i,j,k)
               
                end do
           end do
       end do
        gl%A3_PCSAFTQ(10) = 2.d0*(sum1b+sum1a)*factor2+ 4.d0*(sum2a+sum2b)*factor2 + factor2*(sum3a+sum3b)    
    end if

!DEC$ END IF
end subroutine A3DERIVS


subroutine ADDDERIVS(gl,T,DENS,GETDERADD)
    
    ! a_dd: contribution to the residual helmholz engergy from dipolar interactions 
    ! defined by eq. A.7 in Gross, Sadowski 2006:
    ! a_dd = a_2dd*a_2dd/(a_2dd-a_3dd)
    ! dependent on D and T





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERADD
    double precision, intent (inout) :: T, DENS
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sum1, sum2, factor1, factor2
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !calculate the derivatives of a_qq
    if (GETDERADD(1) .eq. 1) then
        gl%ADD_PCSAFTD(1) = gl%A2_PCSAFTD(1) / (1.d0 - (gl%A3_PCSAFTD(1) / gl%A2_PCSAFTD(1)))
    end if
    
    !  2: 1ST DERIVATIVE OF a_qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires a_2qq, a_3qq
    if (GETDERADD(2) .eq. 1) then
          gl%ADD_PCSAFTD(2) = (2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(2) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) - gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * (gl%A2_PCSAFTD(2) - gl%A3_PCSAFTD(2))) &
            & / ((gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)))
    end if
    
    ! I. Schuelling 05/17
    ! 3: 2ND DERIVATIVE OF a_qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires a_2qq, a_3qq
    if (GETDERADD(3) .eq. 1) then
        gl%ADD_PCSAFTD(3) = (gl%A2_PCSAFTD(1)**3 * gl%A2_PCSAFTD(3) + gl%A2_PCSAFTD(1)**3 * gl%A3_PCSAFTD(3) - 3.d0 * gl%A2_PCSAFTD(1)**2 * gl%A2_PCSAFTD(3) * gl%A3_PCSAFTD(1)  &
         & - gl%A2_PCSAFTD(1)**2 * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(3) + 2.d0 * gl%A2_PCSAFTD(1)**2 * gl%A3_PCSAFTD(2)**2 &
         & - 4.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(2) + 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(3) * gl%A3_PCSAFTD(1)**2 &
         & + 2.d0 * gl%A2_PCSAFTD(2)**2 * gl%A3_PCSAFTD(1)**2) &
         & / (gl%A2_PCSAFTD(1)**3 - 3.d0 * gl%A2_PCSAFTD(1)**2 * gl%A3_PCSAFTD(1) + 3.d0 * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1)**2 - gl%A3_PCSAFTD(1)**3)
     end if
     
    ! 4: 1ST DERIVATIVE OF a_qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires a_2qq, a_3qq
    if (GETDERADD(4) .eq. 1) then
        gl%ADD_PCSAFTD(4) = (2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(4)* ( gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) - gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))) &
            & / ((gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)))
    end if
   
    ! 5: 2ND DERIVATIVE OF a_qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires a_2qq, a_3qq
     if (GETDERADD(5) .eq. 1) then
        gl%ADD_PCSAFTD(5) = ( gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * (-(gl%A2_PCSAFTD(5) - gl%A3_PCSAFTD(5))) + 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(5) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) &
            & + 2.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))) &
            &/ ((gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1))) &
            & - ( 2.d0 * (2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(4) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) - gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4))) * (gl%A2_PCSAFTD(4) - gl%A3_PCSAFTD(4)) ) &
            &/ ((gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) * (gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)) *(gl%A2_PCSAFTD(1) - gl%A3_PCSAFTD(1)))
     end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires a_2qq, a_3qq
    if (GETDERADD(6) .eq. 1) then
        gl%ADD_PCSAFTD(6) = (gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(6) + gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(6) &
            & - 3.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(6) * gl%A3_PCSAFTD(1)  &
            & - gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(6) + 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(4) * gl%A3_PCSAFTD(2) &
            & - 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(4) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(2) &
            & + 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(6) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) - 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(4) &
            & + 2.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1)) &
            & / (gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) - 3.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) + 3.d0 * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) &
            & - gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1))
    end if

    ! 7: 3RD MIXED DERIVATIVE OF a_hc WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires a_2qq, a_3qq
    
    if (GETDERADD(7) .eq. 1) then
        gl%ADD_PCSAFTD(7) = (gl%A2_PCSAFTD(1)**4 * gl%A2_PCSAFTD(7) &
            & + gl%A2_PCSAFTD(1)**4 * gl%A3_PCSAFTD(7) &
            & - 4.d0 * gl%A2_PCSAFTD(1)**3 * gl%A2_PCSAFTD(7) * gl%A3_PCSAFTD(1)  &
            & - 2.d0 * gl%A2_PCSAFTD(1)**3 * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(7) &
            & + 4.d0 * gl%A2_PCSAFTD(1)**3 * gl%A3_PCSAFTD(4) * gl%A3_PCSAFTD(6) &
            & + 2.d0 * gl%A2_PCSAFTD(1)**3 * gl%A3_PCSAFTD(5) * gl%A3_PCSAFTD(2)  &
            & - 4.d0 * gl%A2_PCSAFTD(1)**2 * gl%A2_PCSAFTD(4) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(6) &
            & - 4.d0 * gl%A2_PCSAFTD(1)**2 * gl%A2_PCSAFTD(4) * gl%A3_PCSAFTD(4) * gl%A3_PCSAFTD(2) &
            & - 2.d0 * gl%A2_PCSAFTD(1)**2 * gl%A2_PCSAFTD(5) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(2) &
            & + 5.d0 * gl%A2_PCSAFTD(1)**2 * gl%A2_PCSAFTD(7) * gl%A3_PCSAFTD(1)**2 &
            & - 4.d0 * gl%A2_PCSAFTD(1)**2 * gl%A2_PCSAFTD(6) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(4) &
            & - 2.d0 * gl%A2_PCSAFTD(1)**2 * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(5) &
            & - 2.d0 * gl%A2_PCSAFTD(1)**2 * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(4)**2 &
            & +        gl%A2_PCSAFTD(1)**2 * gl%A3_PCSAFTD(1)**2 * gl%A3_PCSAFTD(7)  &
            & - 4.d0 * gl%A2_PCSAFTD(1)**2 * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(4) * gl%A3_PCSAFTD(6) &
            & - 2.d0 * gl%A2_PCSAFTD(1)**2 * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(5) * gl%A3_PCSAFTD(2)  &
            & + 6.d0 * gl%A2_PCSAFTD(1)**2 * gl%A3_PCSAFTD(4)**2 * gl%A3_PCSAFTD(2) &
            & + 4.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(4)**2 * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(2)  &
            & + 4.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(6) * gl%A3_PCSAFTD(1)**2 &
            & + 8.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(4)  &
            & + 4.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(4) * gl%A3_PCSAFTD(1)**2 * gl%A3_PCSAFTD(6) &
            & - 8.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(4) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(4) * gl%A3_PCSAFTD(2)  &
            & + 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(5) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1)**2 &
            & + 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(5) * gl%A3_PCSAFTD(1)**2 * gl%A3_PCSAFTD(2) &
            & - 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(7) * gl%A3_PCSAFTD(1)**3 &
            & + 4.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(6) * gl%A3_PCSAFTD(1)**2 * gl%A3_PCSAFTD(4)   &
            & + 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1)**2 * gl%A3_PCSAFTD(5) &
            & - 4.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(4)**2   &
            & - 6.d0 * gl%A2_PCSAFTD(4)**2 * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1)**2 &
            & + 2.d0 * gl%A2_PCSAFTD(4)**2 * gl%A3_PCSAFTD(1)**2 * gl%A3_PCSAFTD(2) &
            & - 4.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(6) * gl%A3_PCSAFTD(1)**3 &
            & + 4.d0 * gl%A2_PCSAFTD(4) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1)**2 * gl%A3_PCSAFTD(4) &
            & - 2.d0 * gl%A2_PCSAFTD(5) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(1)**3) &
            & /(gl%A2_PCSAFTD(1)**4 &
            & - 4.d0 * gl%A2_PCSAFTD(1)**3 * gl%A3_PCSAFTD(1) &
            & + 6.d0 * gl%A2_PCSAFTD(1)**2 * gl%A3_PCSAFTD(1)**2 &
            & - 4.d0 * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1)**3 &
            & + gl%A3_PCSAFTD(1)**4)
 
    end if
    
    ! I. Schuelling 05/17
    ! 8: 3ND DERIVATIVE OF a_qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires a_2qq, a_3qq
    if (GETDERADD(8) .eq. 1) then
        gl%ADD_PCSAFTD(8) = 1/((gl%A2_PCSAFTD(1)-gl%A3_PCSAFTD(1)) * (gl%A2_PCSAFTD(1)-gl%A3_PCSAFTD(1)) * (gl%A2_PCSAFTD(1)-gl%A3_PCSAFTD(1)) * (gl%A2_PCSAFTD(1)-gl%A3_PCSAFTD(1)))&
            & * ( - 4.d0 * gl%A2_PCSAFTD(8) * gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) &
            & + 5.d0 *gl%A2_PCSAFTD(8) * gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(1) &
            & - 2.d0 * gl%A2_PCSAFTD(8) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(1) &
            & + gl%A2_PCSAFTD(8) * gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(1) &
            & - 6.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2_PCSAFTD(3) * gl%A3_PCSAFTD(2) &
            & + 6.d0 * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2_PCSAFTD(3) * gl%A3_PCSAFTD(2) &
            & - 6.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(3) &
            & + 6.d0 * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(3) &
            & - 6.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(2) * gl%A3_PCSAFTD(2) &
            & - 12.d0 * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(2) * gl%A3_PCSAFTD(2) &
            & + 12.d0 * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(2) &
            & + 6.d0 * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A2_PCSAFTD(2) * gl%A3_PCSAFTD(2) &
            & - 6.d0 * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A2_PCSAFTD(2) * gl%A2_PCSAFTD(2) &
            & + 6.d0 * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A2_PCSAFTD(3) &
            & - 6.d0 * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A2_PCSAFTD(2) * gl%A2_PCSAFTD(3) &
            & + gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(8) &
            & - 2.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(8) &
            & + gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(8) &
            & + 6.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(2) * gl%A3_PCSAFTD(2) * gl%A3_PCSAFTD(2) &
            & + 6.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(2) * gl%A3_PCSAFTD(3) &
            & - 6.d0 * gl%A2_PCSAFTD(1) * gl%A2_PCSAFTD(1) * gl%A3_PCSAFTD(1) * gl%A3_PCSAFTD(2) * gl%A3_PCSAFTD(3) ) 
    end if   
    
    !9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    if (GETDERADD(10) .eq. 1) then
        gl%ADD_PCSAFTD(9) = (gl%A2_PCSAFTD(1)**4 *gl%A2_PCSAFTD(9) + gl%A2_PCSAFTD(1)**4 *gl%A3_PCSAFTD(9) - 4.d0*gl%A2_PCSAFTD(1)**3 *gl%A2_PCSAFTD(9)*gl%A3_PCSAFTD(1)  &
            & - 2.d0*gl%A2_PCSAFTD(1)**3 *gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(9) + 6.d0*gl%A2_PCSAFTD(1)**3 *gl%A3_PCSAFTD(4)*gl%A3_PCSAFTD(5) &
            & - 6.d0*gl%A2_PCSAFTD(1)**2 *gl%A2_PCSAFTD(4)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(5) - 6.d0*gl%A2_PCSAFTD(1)**2 *gl%A2_PCSAFTD(4)*gl%A3_PCSAFTD(4)**2  &
            & - 6.d0*gl%A2_PCSAFTD(1)**2 *gl%A2_PCSAFTD(5)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(4) + 5.d0*gl%A2_PCSAFTD(1)**2 *gl%A2_PCSAFTD(9)*gl%A3_PCSAFTD(1)**2  &
            & + gl%A2_PCSAFTD(1)**2 *gl%A3_PCSAFTD(1)**2 *gl%A3_PCSAFTD(9) - 6.d0*gl%A2_PCSAFTD(1)**2 *gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(4)*gl%A3_PCSAFTD(5) &
            & + 6.d0*gl%A2_PCSAFTD(1)**2 *gl%A3_PCSAFTD(4)**3 + 12.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(4)**2 *gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(4) &
            & + 6.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(4)*gl%A2_PCSAFTD(5)*gl%A3_PCSAFTD(1)**2 + 6.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(4)*gl%A3_PCSAFTD(1)**2 *gl%A3_PCSAFTD(5) &
            & - 12*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(4)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(4)**2 + 6*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(5)*gl%A3_PCSAFTD(1)**2*gl%A3_PCSAFTD(4) &
            & - 2.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(9)*gl%A3_PCSAFTD(1)**3 - 6.d0*gl%A2_PCSAFTD(4)**3 *gl%A3_PCSAFTD(1)**2 + 6.d0*gl%A2_PCSAFTD(4)**2 *gl%A3_PCSAFTD(1)**2 *gl%A3_PCSAFTD(4) &
            & - 6.d0*gl%A2_PCSAFTD(4)*gl%A2_PCSAFTD(5)*gl%A3_PCSAFTD(1)**3) &
            & /(gl%A2_PCSAFTD(1)**4 - 4.d0*gl%A2_PCSAFTD(1)**3 *gl%A3_PCSAFTD(1) + 6.d0*gl%A2_PCSAFTD(1)**2 *gl%A3_PCSAFTD(1)**2 - 4.d0*gl%A2_PCSAFTD(1)*gl%A3_PCSAFTD(1)**3 + gl%A3_PCSAFTD(1)**4)
            &
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF ADD WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires a_2qq, a_3qq
    if (GETDERADD(10) .eq. 1) then
        gl%ADD_PCSAFTD(10) = (gl%A2_PCSAFTD(1)**4 *gl%A2_PCSAFTD(10) + gl%A2_PCSAFTD(1)**4 *gl%A3_PCSAFTD(10) - 4.d0*gl%A2_PCSAFTD(1)**3 *gl%A2_PCSAFTD(10)*gl%A3_PCSAFTD(1)  &
            & - 2.d0*gl%A2_PCSAFTD(1)**3 *gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(10) + 2.d0*gl%A2_PCSAFTD(1)**3 *gl%A3_PCSAFTD(4)*gl%A3_PCSAFTD(3) &
            & + 4.d0*gl%A2_PCSAFTD(1)**3 *gl%A3_PCSAFTD(6)*gl%A3_PCSAFTD(2) - 2.d0*gl%A2_PCSAFTD(1)**2 *gl%A2_PCSAFTD(4)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(3) &
            & - 2.d0*gl%A2_PCSAFTD(1)**2 *gl%A2_PCSAFTD(4)*gl%A3_PCSAFTD(2)**2 - 4.d0*gl%A2_PCSAFTD(1)**2 *gl%A2_PCSAFTD(6)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(2)  &
            & + 5.d0*gl%A2_PCSAFTD(1)**2 *gl%A2_PCSAFTD(10)*gl%A3_PCSAFTD(1)**2 - 4.d0*gl%A2_PCSAFTD(1)**2*gl%A2_PCSAFTD(2)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(6)  &
            & - 4.d0*gl%A2_PCSAFTD(1)**2 *gl%A2_PCSAFTD(2)*gl%A3_PCSAFTD(4)*gl%A3_PCSAFTD(2) - 2.d0*gl%A2_PCSAFTD(1)**2*gl%A2_PCSAFTD(3)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(4)  &
            & + gl%A2_PCSAFTD(1)**2 *gl%A3_PCSAFTD(1)**2 *gl%A3_PCSAFTD(10) - 2.d0*gl%A2_PCSAFTD(1)**2*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(4)*gl%A3_PCSAFTD(3)  &
            & - 4.d0*gl%A2_PCSAFTD(1)**2 *gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(6)*gl%A3_PCSAFTD(2) + 6.d0*gl%A2_PCSAFTD(1)**2 *gl%A3_PCSAFTD(4)*gl%A3_PCSAFTD(2)**2 &
            & + 8.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(4)* gl%A2_PCSAFTD(2)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(2) + 2.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(4)*gl%A2_PCSAFTD(3)*gl%A3_PCSAFTD(1)**2&
            & + 2.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(4) *gl%A3_PCSAFTD(1)**2 *gl%A3_PCSAFTD(3) - 4.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(4)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(2)**2  &
            & + 4.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(6)*gl%A2_PCSAFTD(2)*gl%A3_PCSAFTD(1)**2 + 4.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(6)*gl%A3_PCSAFTD(1)**2 *gl%A3_PCSAFTD(2)  &
            & - 2.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(10)*gl%A3_PCSAFTD(1)**3 + 4.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(2)**2 *gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(4)  &
            & + 4.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(2)*gl%A3_PCSAFTD(1)**2 *gl%A3_PCSAFTD(6) - 8.d0*gl%A2_PCSAFTD(1)*gl%A2_PCSAFTD(2)*gl%A3_PCSAFTD(1)*gl%A3_PCSAFTD(4)*gl%A3_PCSAFTD(2) &
            & + 2.d0*gl%A2_PCSAFTD(1)* gl%A2_PCSAFTD(3)*gl%A3_PCSAFTD(1)**2 *gl%A3_PCSAFTD(4) - 6.d0*gl%A2_PCSAFTD(4)*gl%A2_PCSAFTD(2)**2*gl%A3_PCSAFTD(1)**2  &
            & + 4.d0*gl%A2_PCSAFTD(4)*gl%A2_PCSAFTD(2)*gl%A3_PCSAFTD(1)**2*gl%A3_PCSAFTD(2) - 2.d0*gl%A2_PCSAFTD(4)*gl%A2_PCSAFTD(3)*gl%A3_PCSAFTD(1)**3 &
            & - 4.d0*gl%A2_PCSAFTD(6)*gl%A2_PCSAFTD(2)*gl%A3_PCSAFTD(1)**3 + 2.d0*gl%A2_PCSAFTD(2)**2*gl%A3_PCSAFTD(1)**2*gl%A3_PCSAFTD(4))/(gl%A2_PCSAFTD(1)**4  &
            & - 4.d0*gl%A2_PCSAFTD(1)**3 *gl%A3_PCSAFTD(1) + 6.d0*gl%A2_PCSAFTD(1)**2*gl%A3_PCSAFTD(1)**2 - 4.d0*gl%A2_PCSAFTD(1)*gl%A3_PCSAFTD(1)**3 + gl%A3_PCSAFTD(1)**4)
    
    end if 
    
!DEC$ END IF
end subroutine ADDDERIVS

subroutine A2DERIVS_D(gl,T,DENS,GETDERADD) !Muss noch korriegiert werden! Siehe Sympy
   
    ! a_2dd: second-order perturbation term of the dipolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.8 in Gross, Sadowski 2006:
    ! a_2dd =-pi**rho*sum_i*sum_j(x_i*x_j*epsk_ii/T*espk_jj/T*sig_ii**3*sig_jj**3/sig_ij**3*n_Di*n_Dj*My_i**2*My_j**2*J2_D_ij
    ! dependent on D and T







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERADD
    double precision, intent (in) :: T, DENS
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sum3, sum1,sum1a, sum1b, sum2, sum2a,sum2b, sum4, sum4b, sum5, sum5b, sum5c, sum6, sum6b, sum7c,sum8, sum9,sum10, sigma_ij, factor_lc
    integer :: i, j
    double precision, dimension(2) :: dd
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !Initializations
    sum3 = 0.d0
    sum1 = 0.d0
    sum2 = 0.d0
    
    sigma_ij = 0.d0
    factor_lc = -piPCSAFT*dens*1.d-38 ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19
    
   !calculate the derivatives of a_qq
    ! 1: a_2qq
    if (GETDERADD(1) .eq. 1) then
            sum1 = 0.d0
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                    sum1 = sum1 + gl%J2_PCSAFTD(1,i,j)*gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 *gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                        & /(T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                end do
            end do
        
            gl%A2_PCSAFTD(1) = factor_lc * sum1 
    end if
    
    !  2: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires j_2qq
    if (GETDERADD(2) .eq. 1) then
        sum1 = 0.d0
        sum2 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                sum2 = sum2 + gl%J2_PCSAFTD(2,i,j)*gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 *gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                    & /(T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                sum1 = sum1 + gl%J2_PCSAFTD(1,i,j)*gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 *gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                    & /(T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))    
            end do
        end do
    
        gl%A2_PCSAFTD(2) = factor_lc*sum2 + factor_lc*sum1 

    end if
    
   ! I. Schuelling 05/17:
   ! 3: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
   ! requires j_2qq
    if (GETDERADD(3) .eq. 1) then
        sum3 = 0.d0
        sum2 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                
                sum2 = sum2 + gl%J2_PCSAFTD(3,i,j)*gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j)/(T**2 &
                    & *kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                sum3 = sum3 + gl%J2_PCSAFTD(2,i,j)*gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                    & /(T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                
            end do
        end do
        gl%A2_PCSAFTD(3) = factor_lc*(sum2 + 2.d0*sum3)
     
    end if
 
    ! 4: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires j_2qq
     if (GETDERADD(4) .eq. 1) then
        sum4 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                sum4 = sum4 + (-2.d0*gl%J2_PCSAFTD(1,i,j)+ gl%J2_PCSAFTD(4,i,j))*gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                    & /(T**2 *kbol**2 *sigma_ij**3 *gl%mPCSAFT(i)*gl%mPCSAFT(j)) 
            end do
        end do
        
        gl%A2_PCSAFTD(4) = factor_lc*sum4
     end if
    
    ! 5: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
     !requires j_2qq
     if (GETDERADD(5) .eq. 1) then
        sum5 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                sum5 = sum5 + gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 *gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*(6.d0*gl%J2_PCSAFTD(1,i,j) - 4.d0*gl%J2_PCSAFTD(4,i,j) + gl%J2_PCSAFTD(5,i,j))*gl%molfractions(i)*gl%molfractions(j)&
                    & /(T**2 *kbol**2 *sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
            end do
        end do
        gl%A2_PCSAFTD(5) = factor_lc*sum5
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_2qq
    if (GETDERADD(6) .eq. 1) then
     sum1 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                sum1 = sum1 + gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 *gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                    & *(gl%J2_PCSAFTD(6,i,j) - 2.d0*gl%J2_PCSAFTD(2,i,j) - 2.d0*gl%J2_PCSAFTD(1,i,j) + gl%J2_PCSAFTD(4,i,j)) &
                    & /(T**2 *kbol**2 *sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j)) 
            end do
        end do
     
        gl%A2_PCSAFTD(6) =  factor_lc*sum1 
    end if
    
    ! Schuelling 06/17:
    ! 7: 3RD MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires j_2qq
    if (GETDERADD(7) .eq. 1) then
        sum1 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                
                sum1 = sum1 + gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 *gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                    & *(gl%J2_PCSAFTD(7,i,j) - 4.d0*gl%J2_PCSAFTD(6,i,j) + 6.d0*gl%J2_PCSAFTD(2,i,j) + 6.d0*gl%J2_PCSAFTD(1,i,j) - 4.d0*gl%J2_PCSAFTD(4,i,j) + gl%J2_PCSAFTD(5,i,j)) &
                    & /(T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
            end do
        end do
        gl%A2_PCSAFTD(7) = factor_lc*sum1
    end if 
    
    ! I. Schuelling 05/17:
    ! 8: 3RD DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires j2_qq
    if (GETDERADD(8) .eq. 1) then
        sum2 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                
                sum2 = sum2 + (gl%J2_PCSAFTD(8,i,j) + 3.d0*gl%J2_PCSAFTD(3,i,j))*gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 *gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                    & /(T**2 *kbol**2 *sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
            end do
        end do
        gl%A2_PCSAFTD(8) = sum2*factor_lc
    end if
   
    ! 9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires 
    if (GETDERADD(9) .eq. 1) then
    sum1 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0                     
                    sum1 = sum1 + gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 *gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                        & *(-24.d0*gl%J2_PCSAFTD(1,i,j) + 18.d0*gl%J2_PCSAFTD(4,i,j) - 6.d0*gl%J2_PCSAFTD(5,i,j) + gl%J2_PCSAFTD(9,i,j)) &
                        & /(T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j)) 
                          
            end do
        end do
        gl%A2_PCSAFTD(9) = factor_lc*sum1
    end if
    
    !Schuelling: 06/17
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires j_2qq
    if (GETDERADD(10) .eq. 1) then
            sum3 = 0.d0
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                
                    sum3 = sum3 + gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 *gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                        & *(gl%J2_PCSAFTD(10,i,j) - 2.d0*gl%J2_PCSAFTD(3,i,j) + 2.d0*gl%J2_PCSAFTD(6,i,j) - 4.d0*gl%J2_PCSAFTD(2,i,j)) &
                        & /(T**2 *kbol**2 *sigma_ij**3 *gl%mPCSAFT(i)*gl%mPCSAFT(j))
    
                end do
            end do
            gl%A2_PCSAFTD(10) = factor_lc*sum3
    end if   
!DEC$ END IF
end subroutine A2DERIVS_D

subroutine A3DERIVS_D (gl,T, DENS, GETDERADD) !Muss noch korriegiert werden! Siehe Sympy
    
    ! a_3dd: second-order perturbation term of the dipolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.9 in Gross, Sadowski 2006:
    ! a_3dd =-4*pi**2/3*rho**2*sum_i*sum_j*sum_k(x_i*x_j*x_k*epsk_ii/T*espk_jj/T*espk_kk/T*sig_ii**3*sig_jj**3*sig_kk**3/sig_ij/sig_ik/sig_jk*n_Di*n_Dj*n_Dk*My_i**2*My_j**2*My_k**2*J3_D_ijk
    ! dependent on D and T







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERADD
    double precision, intent (in) :: T, DENS
    !working variables
    double precision :: sum1, sum1a, sum1b, sum2, sum2a, sum2b, sum3, sum3a, sum3b, sum4, sum4a, sum4b, sum5, sum5a, sum5b, sum5c, sum6, sum6b, sum7c,sum8, sum9,sum10, factor2
    double precision, dimension(gl%ncomp,gl%ncomp) :: sigma_ij
    integer :: i, j, k
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !Initializations
    sigma_ij = 0.d0
    factor2 = -piPCSAFT*piPCSAFT*4.d0/3.d0*dens*dens * 1.d-57  ! 1.d-57: comes from units in dimensionless Myi*^2 = 1.d-19
    do j = gl%n_start, gl%n_end
        do i = gl%n_start, gl%n_end
            sigma_ij(i,j) = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
        end do
    end do
   
    !eq. 10 in Gross (2005). J3,ij = 0, which means that the first summation equals zero and is therefore left out
    !calculate the derivatives of a_3qq
    if (GETDERADD(1) .eq. 1) then
        sum1= 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)*gl%J3_PCSAFTD(1,i,j,k) &
                         & /(T**3*kbol**3*sigma_ij(i,j)*sigma_ij(i,k)*sigma_ij(j,k)*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                end do
            end do
        end do
        gl%A3_PCSAFTD(1) = sum1*factor2
    end if
    
    !  2: 1ST DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires J_3ijk
    if (GETDERADD(2) .eq. 1) then
        sum1 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)*(2.d0*gl%J3_PCSAFTD(1,i,j,k) + gl%J3_PCSAFTD(2,i,j,k)) &
                         & /(T**3*kbol**3*sigma_ij(i,j)*sigma_ij(i,k)*sigma_ij(j,k)*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                end do
            end do
        end do
        
        gl%A3_PCSAFTD(2) = sum1*factor2
    end if
        
    ! 3: 2ND DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires j_3ijk
    if (GETDERADD(3) .eq. 1) then
       sum1 = 0.d0
       do j = gl%n_start, gl%n_end
           do i = gl%n_start, gl%n_end
               do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)*(2.d0*gl%J3_PCSAFTD(1,i,j,k) + 4.d0*gl%J3_PCSAFTD(2,i,j,k) + gl%J3_PCSAFTD(3,i,j,k)) &
                         & /(T**3*kbol**3*sigma_ij(i,j)*sigma_ij(i,k)*sigma_ij(j,k)*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                
                end do
           end do
       end do
        gl%A3_PCSAFTD(3) = sum1*factor2
        end if
    
    ! 4: 1ST DERIVATIVE OF a_3qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires j_3ijk
    if (GETDERADD(4) .eq. 1) then
        sum4= 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum4 = sum4 + gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)*(-3.d0*gl%J3_PCSAFTD(1,i,j,k) + gl%J3_PCSAFTD(4,i,j,k)) &
                         & /(T**3*kbol**3*sigma_ij(i,j)*sigma_ij(i,k)*sigma_ij(j,k)*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                
                end do
            end do
        end do
        gl%A3_PCSAFTD(4) = factor2*sum4
    end if
    
    ! 5: 2ND DERIVATIVE OF a_3qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires j_3ijk
      if (GETDERADD(5) .eq. 1) then
        sum5= 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                   sum5 = sum5 + gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)*(12.d0*gl%J3_PCSAFTD(1,i,j,k) - 6.d0*gl%J3_PCSAFTD(4,i,j,k) + gl%J3_PCSAFTD(5,i,j,k)) &
                         & /(T**3*kbol**3*sigma_ij(i,j)*sigma_ij(i,k)*sigma_ij(j,k)*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                end do
            end do
        end do
        gl%A3_PCSAFTD(5) = factor2 * sum5
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF a_3 WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_3ijk
    if (GETDERADD(6) .eq. 1) then
      sum1 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                sum1 = sum1 + gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)*(gl%J3_PCSAFTD(6,i,j,k) - 3.d0*gl%J3_PCSAFTD(2,i,j,k) - 6.d0*gl%J3_PCSAFTD(1,i,j,k) + 2.d0*gl%J3_PCSAFTD(4,i,j,k)) &
                         & /(T**3*kbol**3*sigma_ij(i,j)*sigma_ij(i,k)*sigma_ij(j,k)*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                end do
            end do
        end do
        gl%A3_PCSAFTD(6) = factor2*sum1 

    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_3qq WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires j_3ijk
    if (GETDERADD(7) .eq. 1) then
        sum1 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)*(gl%J3_PCSAFTD(7,i,j,k) - 6.d0*gl%J3_PCSAFTD(6,i,j,k) + 12.d0*gl%J3_PCSAFTD(2,i,j,k) + 24.d0*gl%J3_PCSAFTD(1,i,j,k) - 12.d0*gl%J3_PCSAFTD(4,i,j,k) + 2.d0*gl%J3_PCSAFTD(5,i,j,k)) &
                         & /(T**3*kbol**3*sigma_ij(i,j)*sigma_ij(i,k)*sigma_ij(j,k)*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                end do
            end do
        end do
        gl%A3_PCSAFTD(7) = factor2*sum1 
    end if
    
    ! 8: 3RD DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires j_3ijk
    if (GETDERADD(8) .eq. 1) then
       sum8 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum8 = sum8 + gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)*(gl%J3_PCSAFTD(8,i,j,k) + 6.d0*gl%J3_PCSAFTD(3,i,j,k) + 6.d0*gl%J3_PCSAFTD(2,i,j,k)) &
                         & /(T**3*kbol**3*sigma_ij(i,j)*sigma_ij(i,k)*sigma_ij(j,k)*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                end do
            end do
        end do
        gl%A3_PCSAFTD(8) = sum8*factor2

    end if
   ! 9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires 
    if (GETDERADD(9) .eq. 1) then
        sum1 = 0.d0
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)*(-60.d0*gl%J3_PCSAFTD(1,i,j,k) + 36.d0*gl%J3_PCSAFTD(4,i,j,k) - 9.d0*gl%J3_PCSAFTD(5,i,j,k) + gl%J3_PCSAFTD(9,i,j,k)) &
                         & /(T**3*kbol**3*sigma_ij(i,j)*sigma_ij(i,k)*sigma_ij(j,k)*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                end do
            end do
        end do
       gl%A3_PCSAFTD(9) = factor2*sum1
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires j_3ijk
    if (GETDERADD(10) .eq. 1) then
       sum1 = 0.d0
       
       do j = gl%n_start, gl%n_end
           do i = gl%n_start, gl%n_end
               do k = gl%n_start, gl%n_end
                    sum1 = sum1 + gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)*(gl%J3_PCSAFTD(10,i,j,k) - 3.d0*gl%J3_PCSAFTD(3,i,j,k) + 4.d0*gl%J3_PCSAFTD(6,i,j,k) - 12.d0*gl%J3_PCSAFTD(2,i,j,k) - 6.d0*gl%J3_PCSAFTD(1,i,j,k) + 2.d0*gl%J3_PCSAFTD(4,i,j,k)) &
                         & /(T**3*kbol**3*sigma_ij(i,j)*sigma_ij(i,k)*sigma_ij(j,k)*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
               end do
           end do
       end do
        gl%A3_PCSAFTD(10) = factor2*sum1
    end if 

 !DEC$ END IF
end subroutine A3DERIVS_D

subroutine J2DERIVS_D(gl,T, DENS, GETDERI)

! I.Schuelling 07/17
    
    ! J2: integral of the perturbation theory
    ! defined by eq. eq 10 in Gross 2006:
    ! J2_D = Sum_0_4(a_nij+b_nij*eps/T)*eta**n
    ! dependent on D and T
    ! derivative stored in module variable





implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (inout) :: GETDERI
    double precision, intent (in) :: T, DENS
    !output: j2_PCSAFT (module variables)
    integer :: i,j,n
    double precision :: eps_ij, sum1, sum2, sum3, sum4, sum5
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    gl%J2_PCSAFTD = 0.d0
    eps_ij = 0.d0
    call add_T_conversion_derivs(gl,GETDERI)
    
  !calculate the derivatives of J_2
    ! 1: J_index
    if (GETDERI(1) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTD(1,i,j) = gl%J2_PCSAFTD(1,i,j) + (gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T)*gl%z3_PCSAFT(1)**n  
                end do
            end do
        end do
    end if
    
    
    !  2: 1ST DERIVATIVE OF J_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                gl%J2_PCSAFTD(2,i,j) = gl%J2_PCSAFTD(2,i,j) + n*gl%z3_PCSAFT(1)**n*(gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T)
                end do
            end do
        end do
    end if
    
    
    ! 3: 2ND DERIVATIVE OF J_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTD(3,i,j) = gl%J2_PCSAFTD(3,i,j) + n*(n-1)*gl%z3_PCSAFT(1)**n*(gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T)
                 end do
            end do
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
          do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTD(4,i,j) = gl%J2_PCSAFTD(4,i,j) - gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T*gl%z3_PCSAFT(1)**n + n*gl%z3_PCSAFT(1)**(n - 1)*gl%z3_PCSAFT(4)*(gl%ab_PCSAFTD(1,n,i,j) &
                        & + gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T)
                end do
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF J2 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
          do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                
                gl%J2_PCSAFTD(5,i,j) = gl%J2_PCSAFTD(5,i,j) + (gl%z3_PCSAFT(1)**(n - 2)*(T**2*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j) &
                    & * eps_ij) + T**2*n*gl%z3_PCSAFT(4)**2*(n - 1)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) - 2.d0*T**2*gl%ab_PCSAFTD(2,n,i,j) &
                    & * eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) + T**2*2.d0*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2)/T**3)
                
                end do
                
            end do
        end do
    end if
    
    ! 6: 2ND MIXED DERIVATIVE OF J2 WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    if (GETDERI(6) .eq. 1) then
     do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                   
                gl%J2_PCSAFTD(6,i,j) = gl%J2_PCSAFTD(6,i,j) + (n**2*gl%z3_PCSAFT(1)**(n - 2)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*(gl%ab_PCSAFTD(1,n,i,j)+ gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T) &
                    & - n*gl%z3_PCSAFT(1)**(n - 2)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*(gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T)&
                    & + n*gl%z3_PCSAFT(1)**(n - 1)*gl%z3_PCSAFT(6)*(gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T) &
                    & - gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T*n*gl%z3_PCSAFT(1)**(n - 1)*gl%z3_PCSAFT(2))
                 end do
            end do
        end do   
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF J2 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires
    if (GETDERI(7) .eq. 1) then 
    do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTD(7,i,j) = gl%J2_PCSAFTD(7,i,j) + (n*gl%z3_PCSAFT(1)**(n - 3)*(T*(T*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(7)*(T* &
                        & gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) &
                        & + T*gl%z3_PCSAFT(1)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(2.d0*n*gl%z3_PCSAFT(4) &
                        & * gl%z3_PCSAFT(6) - 2.d0*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(6) - gl%z3_PCSAFT(5)*gl%z3_PCSAFT(2)) &
                        & + 2.d0*T*gl%z3_PCSAFT(4)**2*gl%z3_PCSAFT(2)*(-n + 1)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) & 
                        & - 2.d0*T*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(6) &
                        & + 2.d0*T*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)) & 
                        & + gl%z3_PCSAFT(2)*(T**2*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) &
                        & + T**2*n*gl%z3_PCSAFT(4)**2*(n - 1)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) &
                        & - 2.d0*T**2*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4) &
                        & + 2.d0*T**2*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**2))/T**3)
                end do
            end do
        end do
    end if
    
    ! 8: 3RD DERIVATIVE OF J_2 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    if (GETDERI(8) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                    gl%J2_PCSAFTD(8,i,j) = gl%J2_PCSAFTD(8,i,j) + n*gl%z3_PCSAFT(1)**n*(gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij/T)&
                        &*(n**2*gl%z3_PCSAFT(1) - 3*n*gl%z3_PCSAFT(1) + 2*gl%z3_PCSAFT(1))/gl%z3_PCSAFT(1)
                        
                 end do
            end do
        end do
    end if
    
    !! 9: 3RD DERIVATIVE OF J_2 WITH RESPECT TO D AT CONSTANT D, MULTIPLIED BY T^3
    if (GETDERI(9) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4                   
                    gl%J2_PCSAFTD(9,i,j) = gl%J2_PCSAFTD(9,i,j) + (gl%z3_PCSAFT(1)**(n - 3)*(T**3*n*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9)*(T*gl%ab_PCSAFTD(1,n,i,j) &
                        & + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) + 3.d0*T**3*n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(n - 1)*(T*gl%ab_PCSAFTD(1,n,i,j)  &
                        & + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) + T**3*n*gl%z3_PCSAFT(4)**3*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(n**2 - 3*n + 2) &
                        & - 3.d0*T**3*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(5) + 3.d0*T**3*gl%ab_PCSAFTD(2,n,i,j)*eps_ij &
                        & *n*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)**2*(-n + 1) + 6.d0*T**3*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*n*gl%z3_PCSAFT(1)**2* &
                        & gl%z3_PCSAFT(4) - 6.d0*T**3*gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)**3)/T**4)
                 end do
            end do
        end do
    end if
    
    !10: 3RD MIXED DERIVATIVE OF J2 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    if (GETDERI(10) .eq. 1) then
        do j = gl%n_start, gl%n_end
            do i = gl%n_start, gl%n_end
                eps_ij = dsqrt(gl%epskPCSAFT(i)*gl%epskPCSAFT(j))*(1.d0-gl%kij_PCSAFT(i,j))
                do n = 0, 4
                   gl%J2_PCSAFTD(10,i,j) = gl%J2_PCSAFTD(10,i,j) + (n*gl%z3_PCSAFT(1)**(n - 3.d0)*gl%z3_PCSAFT(2)*(2.d0*T*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(6) &
                        & * (n - 1)*(T*gl%ab_PCSAFTD(1,n,i,j) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij) + T*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*(T*gl%ab_PCSAFTD(1,n,i,j) &
                        & + gl%ab_PCSAFTD(2,n,i,j)*eps_ij)*(n**2 - 3*n + 2) + gl%ab_PCSAFTD(2,n,i,j)*eps_ij*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(2)*T*(-n + 1))/T**2) 
                end do
            end do
        end do
    end if

!DEC$ END IF
end subroutine J2DERIVS_D

subroutine J3DERIVS_D(gl,T,DENS,GETDERI)

! I. Schuelling 07/2017
    
    ! J3_D: integral over the reference-fluid pair-correlation function and over three-body correlation functions
    ! defined by eq. eq 11 in Gross 2006:
    ! J3 = Sum_0_4c_nijk*eta**n
    ! dependent on D and T
    




implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    integer, dimension (nderivs), intent (in) :: GETDERI
    double precision, intent (in) :: DENS,T
    !output: j1_PCSAFT, j2_PCSAFT (module variables)
    integer :: i,j,k,n
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    gl%J3_PCSAFTD = 0.d0
    
    !calculate the derivatives of J_3
    ! 1: J_index 
    if (GETDERI(1) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTD(1,i,j,k) = gl%J3_PCSAFTD(1,i,j,k) + gl%c_PCSAFTD(n,i,j,k)*gl%z3_PCSAFT(1)**n
                    end do
                end do
            end do
        end do
    end if
    
    !  2: 1ST DERIVATIVE OF J_3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    if (GETDERI(2) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTD(2,i,j,k) = gl%J3_PCSAFTD(2,i,j,k) + gl%c_PCSAFTD(n,i,j,k)*n*gl%z3_PCSAFT(1)**n
                    end do
                end do
            end do
        end do
    end if
    
    ! 3: 2ND DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    if (GETDERI(3) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTD(3,i,j,k) = gl%J3_PCSAFTD(3,i,j,k) + n*(n-1)*gl%c_PCSAFTD(n,i,j,k)*gl%z3_PCSAFT(1)**(n-2) *gl%z3_PCSAFT(2)*gl%z3_PCSAFT(2)
                    end do
                end do
            end do
        end do
    end if
    
    ! 4: 1ST DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    if (GETDERI(4) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTD(4,i,j,k) = gl%J3_PCSAFTD(4,i,j,k) + gl%c_PCSAFTD(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 1) *gl%z3_PCSAFT(4)
                    end do
                end do
            end do
        end do
    end if
    
    ! 5: 2ND DERIVATIVE OF J3 WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    if (GETDERI(5) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTD(5,i,j,k) = gl%J3_PCSAFTD(5,i,j,k) + gl%c_PCSAFTD(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 2)*(n*gl%z3_PCSAFT(4)**2 &
                            & + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - gl%z3_PCSAFT(4)**2)
                    end do
                end do
            end do
        end do    
    end if
    
     ! 6: 2ND MIXED DERIVATIVE OF J3 WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    if (GETDERI(6) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTD(6,i,j,k) = gl%J3_PCSAFTD(6,i,j,k) + gl%c_PCSAFTD(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 2)*(gl%z3_PCSAFT(1)*gl%z3_PCSAFT(6) &
                            & + gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*(n - 1))
                    end do
                end do
            end do
        end do   
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF J3 WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires 
    if (GETDERI(7) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTD(7,i,j,k) = gl%J3_PCSAFTD(7,i,j,k) + gl%c_PCSAFTD(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 3) *(gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(7) &
                         & + 2.d0* gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(6)*(n - 1) + gl%z3_PCSAFT(4)**2 *gl%z3_PCSAFT(2)*(-n+ 1)  &
                         & + gl%z3_PCSAFT(2)*(n - 1)*(n*gl%z3_PCSAFT(4)**2 + gl%z3_PCSAFT(1)*gl%z3_PCSAFT(5) - gl%z3_PCSAFT(4)**2)) 
                         
                    end do
                end do
            end do
        end do   
    end if
    
    ! 8: 3RD DERIVATIVE OF J3 WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    if (GETDERI(8) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTD(8,i,j,k) = gl%J3_PCSAFTD(8,i,j,k) + n*(n-1)*(n-2) * gl%c_PCSAFTD(n,i,j,k)*gl%z3_PCSAFT(1)**(n-3)*gl%z3_PCSAFT(2)*gl%z3_PCSAFT(2)*gl%z3_PCSAFT(2)
                    end do
                end do
            end do
        end do
    end if
    
    ! 9: 3RD DERIVATIVE OF J3 WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires
    if (GETDERI(9) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTD(9,i,j,k) = gl%J3_PCSAFTD(9,i,j,k) + gl%c_PCSAFTD(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 3)*(gl%z3_PCSAFT(1)**2*gl%z3_PCSAFT(9)  &
                            & + 3.d0*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(4)*gl%z3_PCSAFT(5)*(n - 1) + gl%z3_PCSAFT(4)**3*(n**2 - 3*n + 2))
                    end do
                end do
            end do
        end do   
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF J3 WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires 
    if (GETDERI(10) .eq. 1) then
        do k = gl%n_start, gl%n_end
            do j = gl%n_start, gl%n_end
                do i = gl%n_start, gl%n_end
                    do n = 0, 4
                        gl%J3_PCSAFTD(10,i,j,k) = gl%J3_PCSAFTD(10,i,j,k) + gl%c_PCSAFTD(n,i,j,k)*n*gl%z3_PCSAFT(1)**(n - 3) *gl%z3_PCSAFT(2) &
                            & *(2*gl%z3_PCSAFT(1)*gl%z3_PCSAFT(6)*(n - 1) + gl%z3_PCSAFT(4)*gl%z3_PCSAFT(2)*(n**2 - 3*n + 2))
                    end do
                end do
            end do
        end do   
    end if

    
!DEC$ END IF
end subroutine J3DERIVS_D



subroutine AASSOC_base(gl,TEMP, DENS, AASSOC)
! will be called to build up numerical differentiation






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    double precision, intent (inout) :: TEMP, DENS
    integer, dimension (nderivs) :: GETDERFUNC
    !output: AASSOC_PCSAFT (module variable)
    integer :: i,A
    double precision :: AASSOC, sum, sum2,DP,DM,TP,TM,DP2,DM2,TP2,TM2,DP3,DM3,TP3,TM3
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    GETDERFUNC = 0
    GETDERFUNC(1) = 1
    
    
        call calculate_PCSAFT_functionparts_Trho(gl,TEMP,DENS,GETDERFUNC)
        
        !calculate association strenghts of all molecules and all sites  
        call get_delta_AB(gl,TEMP, DENS)  
        !set starting values for all components and sites of molfraction of sites, not bonded
        gl%xA_PCSAFT(:,:) = 1.d0
        !calculate all molfractions of sites, not bonded
        call get_xA(gl,TEMP,DENS) 
        
        sum2 = 0.d0
        do i = gl%n_start, gl%n_end
            sum = 0.d0
            do A = 1, gl%n_sites(i)
                sum = sum + dlog(gl%xA_PCSAFT(A,i)) - gl%xA_PCSAFT(A,i)*0.5d0 !0.130192440750397
            end do
            sum2 = sum2 + gl%molfractions(i) * (sum + 0.5d0*gl%n_sites(i))
        end do
        
        AASSOC = sum2
        
!DEC$ END IF
end subroutine AASSOC_base


subroutine AASSOC(gl,TEMP, DENS, GETDERAASSOC)

    ! a_assoc: association contribution to the Helmholtz free energy
    ! defined by eq. 21 in Chapman et al. 1990:
    ! dependent on T and D






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    double precision, intent (inout) :: TEMP, DENS
    integer, dimension (nderivs), intent (in) :: GETDERAASSOC
    integer, dimension (nderivs) :: GETDERFUNC
    !output: AASSOC_PCSAFT (module variable)
    integer :: i,A
    double precision :: sum, sum2
    double precision :: DELTA,DP,DM,TP,TM,DP2,DM2,TP2,TM2,DP3,DM3,TP3,TM3
    double precision :: DELTAx,DPx,DMx,TPx,TMx,DP2x,DM2x,TP2x,TM2x,DP3x,DM3x,TP3x,TM3x
    double precision :: DELTAxx,DPxx,DMxx,TPxx,TMxx,DP2xx,DM2xx,TP2xx,TM2xx,DP3xx,DM3xx,TP3xx,TM3xx
    double precision :: AASSOC_P,AASSOC_M, AASSOC_0, AASSOC_P2,AASSOC_M2,AASSOC_PP,AASSOC_MM,AASSOC_PM,AASSOC_MP, AASSOC_P3,AASSOC_M3
    double precision :: AASSOC_21,AASSOC_01,AASSOC_M21,AASSOC_2M1,AASSOC_0M1,AASSOC_M2M1
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    GETDERFUNC = 0
    GETDERFUNC(1) = 1
    
    DELTA = 1.d-5
    DP = DENS * (1.d0 + DELTA)
    DP2 = DENS * (1.d0 + 2.d0*DELTA)
    DM = DENS * (1.d0 - DELTA)
    DM2 = DENS * (1.d0 - 2.d0*DELTA)
    TP = TEMP * (1.d0 + DELTA)
    TP2 = TEMP * (1.d0 + 2.d0*DELTA)
    TM = TEMP * (1.d0 - DELTA)
    TM2 = TEMP * (1.d0 - 2.d0*DELTA)  
    

    DELTAx = 5.d-4
    DPx = DENS * (1.d0 + DELTAx)
    DP2x = DENS * (1.d0 + 2.d0*DELTAx)
    DMx = DENS * (1.d0 - DELTAx)
    DM2x = DENS * (1.d0 - 2.d0*DELTAx)

    DELTAxx = 5.d-3
    DPxx = DENS * (1.d0 + DELTAxx)
    DP2xx = DENS * (1.d0 + 2.d0*DELTAxx)
    DP3xx = DENS * (1.d0 + 3.d0*DELTAxx)
    DMxx = DENS * (1.d0 - DELTAxx)
    DM2xx = DENS * (1.d0 - 2.d0*DELTAxx)
    DM3xx = DENS * (1.d0 - 3.d0*DELTAxx)
    TPxx = TEMP * (1.d0 + DELTAxx)
    TP2xx = TEMP * (1.d0 + 2.d0*DELTAxx)
    TP3xx = TEMP * (1.d0 + 3.d0*DELTAxx)
    TMxx = TEMP * (1.d0 - DELTAxx)
    TM2xx = TEMP * (1.d0 - 2.d0*DELTAxx)
    TM3xx = TEMP * (1.d0 - 3.d0*DELTAxx)
    
    !_____________________________________________________________________________________________________
    !---------basic eq--------------
    if (GETDERAASSOC(1) .eq. 1) then
    !_____________________________________________________________________________________________________
        
        call AASSOC_base(gl,TEMP, DENS, AASSOC_0)
        
        gl%AASSOC_PCSAFT(1) = AASSOC_0
    end if
    
    !_____________________________________________________________________________________________________
    !---------D-Deriv---------------
    if (GETDERAASSOC(2) .eq. 1) then
    !_____________________________________________________________________________________________________

        call AASSOC_base(gl,TEMP, DP, AASSOC_P)
        call AASSOC_base(gl,TEMP, DP2, AASSOC_P2)
        call AASSOC_base(gl,TEMP, DM, AASSOC_M)
        call AASSOC_base(gl,TEMP, DM2, AASSOC_M2)
        
        gl%AASSOC_PCSAFT(2) = (8.d0*AASSOC_P - 8.d0*AASSOC_M - AASSOC_P2 + AASSOC_M2) / (12.d0 * DELTA) 
     end if
        
    !_____________________________________________________________________________________________________
    !---------DD-Deriv--------------
    if (GETDERAASSOC(3) .eq. 1) then
    !_____________________________________________________________________________________________________
    
        call AASSOC_base(gl,TEMP, DPX, AASSOC_P)
        call AASSOC_base(gl,TEMP, DP2x, AASSOC_P2)
        call AASSOC_base(gl,TEMP, DMx, AASSOC_M)
        call AASSOC_base(gl,TEMP, DM2x, AASSOC_M2)
        call AASSOC_base(gl,TEMP, DENS, AASSOC_0)
        
        gl%AASSOC_PCSAFT(3) = (-AASSOC_M2 + 16.d0*AASSOC_M - 30.d0*AASSOC_0 + 16.d0*AASSOC_P - AASSOC_P2) / (12.d0*DELTAx**2) 
    end if        
    
    !_____________________________________________________________________________________________________
    !---------T-Deriv--------------
    if (GETDERAASSOC(4) .eq. 1) then
    !_____________________________________________________________________________________________________
     
        call AASSOC_base(gl,TP, DENS, AASSOC_P)
        call AASSOC_base(gl,TP2, DENS, AASSOC_P2)
        call AASSOC_base(gl,TM, DENS, AASSOC_M)
        call AASSOC_base(gl,TM2, DENS, AASSOC_M2)         

        gl%AASSOC_PCSAFT(4) = (8.d0*AASSOC_P - 8.d0*AASSOC_M - AASSOC_P2 + AASSOC_M2) / (12.d0 * DELTA) 
     end if
    
    
    !_____________________________________________________________________________________________________
    !---------TT-Deriv--------------
    if (GETDERAASSOC(5) .eq. 1) then
    !_____________________________________________________________________________________________________
    
        call AASSOC_base(gl,TP, DENS, AASSOC_P)
        call AASSOC_base(gl,TP2, DENS, AASSOC_P2)
        call AASSOC_base(gl,TM, DENS, AASSOC_M)
        call AASSOC_base(gl,TM2, DENS, AASSOC_M2)
        call AASSOC_base(gl,TEMP, DENS, AASSOC_0)
        
        gl%AASSOC_PCSAFT(5) = (-AASSOC_M2 + 16.d0*AASSOC_M - 30.d0*AASSOC_0 + 16.d0*AASSOC_P - AASSOC_P2) / (12.d0*DELTA**2) 
    end if
    
    !_____________________________________________________________________________________________________
    !---------TD-Deriv--------------
    if (GETDERAASSOC(6) .eq. 1) then
    !_____________________________________________________________________________________________________
    
        call AASSOC_base(gl,TP, DP, AASSOC_PP)
        call AASSOC_base(gl,TP, DM, AASSOC_PM)
        call AASSOC_base(gl,TM, DP, AASSOC_MP)
        call AASSOC_base(gl,TM, DM, AASSOC_MM)
        
        gl%AASSOC_PCSAFT(6) = (AASSOC_PP - AASSOC_PM - AASSOC_MP + AASSOC_MM) / (4.d0*DELTA**2) 
    end if    

    !---------TTD-Deriv-------------
    if (GETDERAASSOC(7) .eq. 1) then
    !_____________________________________________________________________________________________________
   
        call AASSOC_base(gl,TP2xx, DPxx, AASSOC_21)
        call AASSOC_base(gl,TEMP, DPxx, AASSOC_01)
        call AASSOC_base(gl,TM2xx, DPxx, AASSOC_M21)
        call AASSOC_base(gl,TP2xx, DMxx, AASSOC_2M1)
        call AASSOC_base(gl,TEMP, DMxx, AASSOC_0M1)
        call AASSOC_base(gl,TM2xx, DMxx, AASSOC_M2M1)
        
        gl%AASSOC_PCSAFT(7) = (AASSOC_21 - 2.d0*AASSOC_01 + AASSOC_M21 - AASSOC_2M1 + 2.d0*AASSOC_0M1 - AASSOC_M2M1) / (8.d0*DELTAxx**3) 
    end if 
    
    !_____________________________________________________________________________________________________
    !---------DDD-Deriv-------------
    if (GETDERAASSOC(8) .eq. 1) then
    !_____________________________________________________________________________________________________
        
        call AASSOC_base(gl,TEMP, DPxx, AASSOC_P)
        call AASSOC_base(gl,TEMP, DP2xx, AASSOC_P2)
        call AASSOC_base(gl,TEMP, DP3xx, AASSOC_P3)
        call AASSOC_base(gl,TEMP, DMxx, AASSOC_M)
        call AASSOC_base(gl,TEMP, DM2xx, AASSOC_M2)
        call AASSOC_base(gl,TEMP, DM3xx, AASSOC_M3)
        
        gl%AASSOC_PCSAFT(8) = (-AASSOC_P3 + 8.d0*AASSOC_P2 - 13.d0*AASSOC_P + 13.d0*AASSOC_M - 8.d0*AASSOC_M2 + AASSOC_M3) / (8.d0*DELTAxx**3) 
    end if      

        !_____________________________________________________________________________________________________
    !---------TTT-Deriv-------------
    if (GETDERAASSOC(9) .eq. 1) then
    !_____________________________________________________________________________________________________
    
        call AASSOC_base(gl,TPxx, DENS, AASSOC_P)
        call AASSOC_base(gl,TP2xx, DENS, AASSOC_P2)
        call AASSOC_base(gl,TP3xx, DENS, AASSOC_P3)
        call AASSOC_base(gl,TMxx, DENS, AASSOC_M)
        call AASSOC_base(gl,TM2xx, DENS, AASSOC_M2)
        call AASSOC_base(gl,TM3xx, DENS, AASSOC_M3)

        gl%AASSOC_PCSAFT(9) = (-AASSOC_P3 + 8.d0*AASSOC_P2 - 13.d0*AASSOC_P + 13.d0*AASSOC_M - 8.d0*AASSOC_M2 + AASSOC_M3) / (8.d0*DELTAxx**3) 
    end if 

    !---------TDD-Deriv-------------
    if (GETDERAASSOC(10) .eq. 1) then
    !_____________________________________________________________________________________________________
    
        call AASSOC_base(gl,TPxx, DP2xx, AASSOC_21)
        call AASSOC_base(gl,TPxx, DENS, AASSOC_01)
        call AASSOC_base(gl,TPxx, DM2xx, AASSOC_M21)
        call AASSOC_base(gl,TMxx, DP2xx, AASSOC_2M1)
        call AASSOC_base(gl,TMxx, DENS, AASSOC_0M1)
        call AASSOC_base(gl,TMxx, DM2xx, AASSOC_M2M1)
        
        gl%AASSOC_PCSAFT(10) = (AASSOC_21 - 2.d0*AASSOC_01 + AASSOC_M21 - AASSOC_2M1 + 2.d0*AASSOC_0M1 - AASSOC_M2M1) / (8.d0*DELTAxx**3) 
    end if     
!DEC$ END IF
end subroutine AASSOC


subroutine get_xA(gl,TEMP,DENS)
    ! get_xA: molfraction of molecules not bonded to site A
    ! defined by eq. 22 in Chapman et al. 1990:
    ! dependent on T and D







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    double precision, intent (in) :: TEMP,DENS
    !output: xA_PCSAFT (module variable)
    integer :: A,i,B,j,n
    double precision :: sum,P_CALC,P,RHO_J,rhomix_calc,D
    double precision, dimension(4,2) :: xA_PCSAFT_old    
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    if (gl%n_start == gl%n_end) then ! pure fluid
        if (gl%n_sites(gl%n_start) == 2) then
            gl%xA_PCSAFT(1:2,gl%n_start) = (-1.d0 + dsqrt(1.d0 + 4.d0*DENS*gl%delta_AB(gl%n_start,gl%n_start,1,2))) / (2.d0*DENS*gl%delta_AB(gl%n_start,gl%n_start,1,2)) 
        elseif (gl%n_sites(gl%n_start) == 4) then
            gl%xA_PCSAFT(1:4,gl%n_start) = (-1.d0 + dsqrt(1.d0 + 8.d0*DENS*gl%delta_AB(gl%n_start,gl%n_start,1,3))) / (4.d0*DENS*gl%delta_AB(gl%n_start,gl%n_start,1,3)) 
        end if        
    else
        sum = 0.d0
        n = 1
        xA_PCSAFT_old = gl%xA_PCSAFT ! already inititalized on the outside
    
        do while (n < 100) ! maximum of iterations
            do i = gl%n_start, gl%n_end
                do A = 1, gl%n_sites(i)
        
                    do j = gl%n_start, gl%n_end
                        do B = 1, gl%n_sites(j)
                            !sum = sum + xA_PCSAFT(B,j,1)*delta_AB(A,i,B,j,1)
                            sum = sum + gl%molfractions(j)*DENS*gl%xA_PCSAFT(B,j)*gl%delta_AB(i,j,A,B)
                        end do
                    end do
            
                    gl%xA_PCSAFT(A,i) = 1.d0 / (1.d0 + sum)   
                    sum = 0.d0
                end do
            end do

        
            if (maxval(dabs(gl%xA_PCSAFT-xA_PCSAFT_old)) < 1.d-13) then
                exit
            else
                xA_PCSAFT_old = gl%xA_PCSAFT
                n = n+1
            end if
        end do
    end if
    
!DEC$ END IF
end subroutine get_xA


subroutine get_delta_AB(gl,TEMP, DENS)
    ! get_delta_AB: association energy
    ! defined by eq. 24 in Chapman et al. 1990:
    ! dependent on T and D






implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !input
    double precision, intent (inout) :: TEMP, DENS
    !output: delta_AB (module variable)
    integer :: i, j, A, B
    double precision :: kaibjPCSAFT,epsaibj,dij_c,G_IJ_c,dij,G_IJ,K_AB,EPS_AB
    integer, dimension(15) :: GETDERAASSOC
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    GETDERAASSOC = 1
    
    if (gl%n_end > gl%n_start) then !two associating components
    
        kaibjPCSAFT = dsqrt(gl%kabPCSAFT(1)*gl%kabPCSAFT(2)) * (dsqrt(gl%sigPCSAFT(1)*gl%sigPCSAFT(2)) * 2.d0/(gl%sigPCSAFT(1)+gl%sigPCSAFT(2)))**3     !cross parameter KAB: Gross, Sadowski (2002)
        epsaibj = 0.5d0 * (gl%epsabkPCSAFT(1)+gl%epsabkPCSAFT(2))
        
        if (any(gl%components == "water") .and. any(gl%components == "ammonia")) then ! andere Mischungsregeln K. Mejbri, A. Bellagi, Int. J. Refrig. (2006) 
            kaibjPCSAFT = dsqrt(gl%kabPCSAFT(1)*gl%kabPCSAFT(2)) * (dsqrt(gl%sigPCSAFT(1)*gl%sigPCSAFT(2)) * 2.d0/(gl%sigPCSAFT(1)+gl%sigPCSAFT(2)))**3     !to do!!!
            epsaibj = 0.5d0 * (gl%epsabkPCSAFT(1)+gl%epsabkPCSAFT(2))
        end if
        
        !dij_c = 0.5d0*(di_PCSAFT(1,1)+di_PCSAFT(2,1)) !Abwandlung
        
        dij_c = 0.5d0*(gl%sigPCSAFT(1)+gl%sigPCSAFT(2)) ! orginalversion wie in Huang Rados
        call GIJDERIVS(gl,DENS, GETDERAASSOC) ! von i und j erg�nzen!
        G_IJ_c = gl%gij_PCSAFT(1,1)
        
    !else !vermutlich �berfl�ssig
    !    kaibjPCSAFT = kabPCSAFT(n_start)
    !    epsaibj = epsabkPCSAFT(n_start)
    !    dij_c = di_PCSAFT(n_start,1)
    !    G_IJ_c = gii_PCSAFT(n_start,1)
    end if
   
    do i = gl%n_start, gl%n_end
        do A = 1, gl%n_sites(i)
            do j = gl%n_start, gl%n_end            
                do B = 1, gl%n_sites(j)
                    if (i == j) then    ! pure fluid 
                        dij = gl%di_PCSAFT(i,1) ! Abwandlung
                        dij = gl%sigPCSAFT(i) !orginal ver�ffentlichung Hunag and Rados
                        G_IJ = gl%gii_PCSAFT(i,1)
                        K_AB = gl%kabPCSAFT(i)
                        EPS_AB = gl%epsabkPCSAFT(i)
                        
                        if ((gl%n_sites(i) == 4) .and. &  ! 4C site, e.g. water 
                            & (((A==1).and.(B==1)) .or. ((A==1).and.(B==2)) .or. ((A==2).and.(B==1)) .or. ((A==2).and.(B==2)) &
                            & .or. ((A==3).and.(B==3)) .or. ((A==3).and.(B==4)) .or. ((A==4).and.(B==3)) .or. ((A==4).and.(B==4)))) then
                            K_AB = 0.d0
                        else if ((gl%n_sites(i) == 3) .and. &  ! 3B site, e.g. ethanol
                            & (((A==1).and.(B==1)) .or. ((A==2).and.(B==2)) .or. ((A==2).and.(B==3)) .or. ((A==3).and.(B==2)) .or. ((A==3).and.(B==3)))) then
                            K_AB = 0.d0  
                        else if ((gl%n_sites(i) == 2) .and. &  ! 2B site, e.g. ethanol
                            & (((A==1).and.(B==1)) .or. ((A==2).and.(B==2)))) then
                            K_AB = 0.d0  
                        end if
                            
                    else ! cross assoc
                        dij = dij_c
                        G_IJ = G_IJ_c
                        G_IJ_c = gl%gij_PCSAFT(gl%n_start,1)
                        K_AB = kaibjPCSAFT
                        EPS_AB = epsaibj
                        
                        !4C and 3B
                        if ((gl%n_sites(i) == 4 .and. gl%n_sites(j) == 3) .and. &
                            &(((A==1).and.(B==1)) .or. ((A==2).and.(B==1)) .or. ((A==3).and.(B==2)) .or. ((A==3).and.(B==3)) .or. ((A==4).and.(B==2)) .or. ((A==4).and.(B==3)))) then
                            K_AB = 0.d0
                        else if ((gl%n_sites(i) == 3 .and. gl%n_sites(j) == 4) .and. & 
                            &(((A==1).and.(B==1)) .or. ((A==1).and.(B==2)) .or. ((A==2).and.(B==3)) .or. ((A==3).and.(B==3)) .or. ((A==2).and.(B==4)) .or. ((A==3).and.(B==4)))) then
                            K_AB = 0.d0
                        end if
                        
                        !4C and 2B
                        if ((gl%n_sites(i) == 2 .and. gl%n_sites(j) == 4) .and. &
                            &(((A==1).and.(B==1)) .or. ((A==1).and.(B==2)) .or. ((A==2).and.(B==3)) .or. ((A==2).and.(B==4)))) then
                            K_AB = 0.d0
                        else if ((gl%n_sites(j) == 2 .and. gl%n_sites(i) == 4) .and. & 
                            &(((A==1).and.(B==1)) .or. ((A==2).and.(B==1)) .or. ((A==3).and.(B==2)) .or. ((A==4).and.(B==2)))) then
                            K_AB = 0.d0
                        end if
                        
                        !2B and 2B
                        if ((gl%n_sites(i) == 2 .and. gl%n_sites(j) == 2) .and. &
                            &(((A==1).and.(B==1)) .or. ((A==2).and.(B==2)))) then
                            K_AB = 0.d0
                        end if
                            
                    end if
                                
                    gl%delta_AB(i,j,A,B) = dij**3 * G_IJ * K_AB * (dexp(EPS_AB/TEMP) - 1.d0)
                end do
            end do
        end do
    end do
    
    


!DEC$ END IF
end subroutine get_delta_AB
 



    end module pc_saft_module
