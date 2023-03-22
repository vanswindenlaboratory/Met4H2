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

    ! module for file pc-saft_ARX_derivs.f90
    module pc_saft_ARX_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module
    use pc_saft_ancillary_routines_module
    use variables_transformation_module

    use pc_saft_ADISPX_derivs_module
    use pc_saft_AHCX_derivs_module
    use pc_saft_AHSX_derivs_module
    use pc_saft_ANCX_derivs_module

    contains



subroutine ARX1DERIVS(gl,T, D_in, getvector, x1)

! Henning Markgraf, June 2016
    
    ! a_res: residual part of the Helmholtz free energy
    ! defined by eq. A.3 in Gross, Sadowski 2001:
    ! a_res = a_hc + a_disp
    ! dependent on x_i and D [1/Angstrom] and T [K]
    
!--------------------------------------------------------------------------------------------------
! All first, second and third composition derivatives of the residual part of the Helmholtz free
! energy are calculated in this subroutine
!--------------------------------------------------------------------------------------------------
! INPUT PARAMETERS:     
! TEMPERATURE - T       K
! DENSITY     - D_in   mol/m^3
! GETDER      - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
!		only numbers 1 to 10 are calculated in this subroutine, the 15 getvector entries are retained
!		only to grant compliance with the other PC-SAFT routines
!
!               All the listed derivatives up to number 10 are derived by x_i and saved in output matrix x1(10,ncomp)
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
! 
! OUTPUT PARAMETERS: 
! x1	     - a matrix with dimensions (10,ncomp), WITH VALUES
!                   EITHER "0.d0" OR THE RESULTS OF THE DERIVATIVES AS INDICATED IN "GETDER"
!--------------------------------------------------------------------------------------------------







implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    double precision, intent (inout) :: T, D_in
    integer, dimension (nderivs), intent (inout) :: getvector
    integer :: nrsubst
    !output
    double precision, dimension (nx1derivs,gl%ncomp), intent (out) :: x1
    !working variables
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
    nrsubst = 0
    call allocate_arrays_PCSAFT(gl,nrsubst)
    
    if (gl%factor == 1.d3) then   !real fluids are calculated => unit conversion is necessary
    ! III. Change the units of the density from mol/m^3 to 1/Angstrom
        D = D_in*N_A/1.d30
    end if
    
    ! IV. Check if any derivatives of this T, rho, molfractions combination is already saved in the module variables
    get = getvector ! get can be changed if some derivatives have already been calculated
    
    ! Adding the necessary additional T derivatives required for the T to tau conversion:
    call add_T_conversion_derivs(gl,get)
    
    call check_calculated_PCSAFT(gl,T,D,nrsubst,gl%arx1_calculated_PCSAFT,get,already_calculated)
    if (any(get .EQ. 1)) then ! no need to calculate anything if get = 0
        
        ! V. some of the functions require all lower order derivatives, too
        call add_previous_derivs(gl,get,getprevious)
        call check_calculated_PCSAFT(gl,T,D,nrsubst,gl%arx1_calculated_PCSAFT,getprevious,already_calculated)
        
        ! VI. Initialize module variables if T, rho or molfractions different
        if (already_calculated .EQV. .FALSE.) then
            call init_derivs_PCSAFT(gl,nrsubst) ! all saved derivatives are deleted
        end if
        
        ! VII. update the module variables
        call update_modulevariables_PCSAFT(gl,T,D,nrsubst,gl%arx1_calculated_PCSAFT,get)
        
        ! VIII.  The T, rho derivatives of all the function parts, which A_disp and A_hc consists of
        !       are needed for the x1 derivatives
        ! i.e.: d_i, zeta_n, mmean, ab, g_ii, I_1, I_2, meo1, meo2, C_1
        call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
        ! and a_hs
        call AHSDERIVS(gl,get)
        
        ! IX. X1 DERIVATIVES:
        ! The composition derivatives of all the function parts, which A_disp and a_hc consist of
        ! i.e.: zeta_n_x, ab_x, gii_x, I_x, meo_x, C_x
        call calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)
        ! and ahs_x
        call AHSX1DERIVS(gl,get)
    
        ! X. calculate ahc_x and adisp_x, which depend on the previously calculated functions
        call AHCX1DERIVS(gl,get)
        call ADISPX1DERIVS(gl,D,get)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !optional supplements
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !dipol
        if (any(gl%nPCSAFTD(gl%n_start:gl%n_end) /= 0)) then
            call calculate_PCSAFT_functionparts_DD_Trho(gl,T,D,getprevious)
            call calculate_PCSAFT_functionparts_x1_DD(gl,T,D,getprevious)
            call ADDX1DERIVS(gl,T,D,get)
        end if
        
        !quadrupol        
        if (any(gl%nPCSAFTQ(gl%n_start:gl%n_end) /= 0))  then
            call calculate_PCSAFT_functionparts_QQ_Trho(gl,T,D,getprevious)
            call calculate_PCSAFT_functionparts_x1_QQ(gl,T,D,getprevious)
            call AQQX1DERIVS(gl,T,D,get)
        end if
        
        !association
        if (any(gl%kabPCSAFT(gl%n_start:gl%n_end) > 1.d-12)) then 
            call AASSOCX1(gl,T, D, get)
        end if
    
        ! XI. calculate the first composition derivatives of a_res
        do i = 1, nx1derivs
            if (get(i) .eq. 1) then
                gl%arx1_PCSAFT(i,1:gl%ncomp) = gl%ahcx1_PCSAFT(i,1:gl%ncomp) + gl%adispx1_PCSAFT(i,1:gl%ncomp) + gl%ADDx1_PCSAFTD(i,1:gl%ncomp) + gl%AQQx1_PCSAFTQ(i,1:gl%ncomp) + gl%AASSOCx1_PCSAFT(i,1:gl%ncomp)
            end if
        end do
    end if
    
    ! XII. Change the T derivatives to tau derivatives
    call convert_T_derivs_x1(gl,getvector,gl%arx1_PCSAFT,x1)

!DEC$ END IF
end subroutine ARX1DERIVS



subroutine ARX2DERIVS(gl,T, D_in, getvector, x2)

! Henning, March 2016
    
    ! a_res: residual part of the Helmholtz free energy
    ! defined by eq. A.3 in Gross, Sadowski 2001:
    ! a_res = a_hc + a_disp
    ! dependent on x_i and D [1/Angstrom] and T [K]
    
!--------------------------------------------------------------------------------------------------
! All first, second and third composition derivatives of the residual part of the Helmholtz free
! energy are calculated in this subroutine
!--------------------------------------------------------------------------------------------------
! INPUT PARAMETERS:     
! TEMPERATURE - T       K
! DENSITY     - D_in   mol/m^3
! get         - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
!		only numbers 1 to 6 are calculated in this subroutine, the 15 getvector entries are retained
!		only to grant compliance with the other PC-SAFT routines
!
!               All the listed derivatives up to number 10 are derived by x_i and x_j and saved in output matrix x2(6,ncomp,ncomp)
!                1. NORMALIZED RESIDUAL MIXTURE HELMHOLTZ ENERGY F AS A FUNCTION OF DEL AND TAU AND X 
!                2. 1ST DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
!                3. 2ND DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL^2
!                4. 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
!                5: 2ND DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU^2
!                6: 1ST MIXED DERIVATIVE OF F WITH RESPECT TO D AND T, MULTIPLIED BY TAU*DEL
! 
! OUTPUT PARAMETERS: 
! x2	     - a matrix with dimensions (6,ncomp,ncomp), WITH VALUES
!                   EITHER "0.d0" OR THE RESULTS OF THE DERIVATIVES AS INDICATED IN "GETDER"
!--------------------------------------------------------------------------------------------------







implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    double precision, intent (inout) :: T, D_in
    integer, dimension (nderivs), intent (inout) :: getvector
    integer :: nrsubst
    !output
    double precision, dimension (nx2derivs,gl%ncomp,gl%ncomp), intent (out) :: x2
    !working variables
    double precision :: D
    integer, dimension (nderivs) :: get
    integer, dimension (nderivs) :: getprevious
    integer :: i, xi, xj
    logical :: already_calculated
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. Initializations
    nrsubst = 0
    call allocate_arrays_PCSAFT(gl,nrsubst)

    
    if (gl%factor == 1.d3) then   !real fluids are calculated => unit conversion is necessary
    ! III. Change the units of the density from mol/m^3 to 1/Angstrom
        D = D_in*N_A/1.d30
    end if
    
    ! IV. Check if any derivatives of this T, rho, molfractions combination is already saved in the module variables
    get = getvector ! get can be changed if some derivatives have already been calculated
    
    ! Adding the necessary additional T derivatives required for the T to tau conversion:
    call add_T_conversion_derivs(gl,get)
    
    call check_calculated_PCSAFT(gl,T,D,nrsubst,gl%arx2_calculated_PCSAFT,get,already_calculated)
    if (any(get .EQ. 1)) then ! no need to calculate anything if get = 0
        
        ! V. some of the functions require all lower order derivatives, too
        call add_previous_derivs(gl,get,getprevious)
        call check_calculated_PCSAFT(gl,T,D,nrsubst,gl%arx2_calculated_PCSAFT,getprevious,already_calculated)
        
        ! VI. Initialize module variables if T, rho or molfractions different
        if (already_calculated .EQV. .FALSE.) then
            call init_derivs_PCSAFT(gl,nrsubst) ! all saved derivatives are deleted
        end if
        
        ! VII. update the module variables
        call update_modulevariables_PCSAFT(gl,T,D,nrsubst,gl%arx2_calculated_PCSAFT,get)
        
        ! VIII.  The T, rho derivatives of all the function parts, which A_disp and A_hc consists of
        !       are needed for the x2 derivatives
        ! i.e.: d_i, zeta_n, mmean, ab, g_ii, I_1, I_2, meo1, meo2, C_1
        call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
        ! and a_hs
        call AHSDERIVS(gl,get)
        
        ! IX. The x1 derivatives are also needed for the x2 derivatives
        ! i.e.: zeta_n_x, ab_x, gii_x, I_x, meo_x, C_x
        call calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)
        ! and ahs_x
        call AHSX1DERIVS(gl,get)
        
        ! X. X2 DERIVATIVES:
        ! ab_xx, gii_xx, I_xx, meo_xx, C_xx
        call calculate_PCSAFT_functionparts_x2(gl,T,D,getprevious)
        ! and ahs_xx
        call AHSX2DERIVS(gl,get)
        
        ! XI. calculate ahc_xx and adisp_xx, which depend on the previously calculated functions
        call AHCX2DERIVS(gl,get)
        call ADISPX2DERIVS(gl,D,get)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !optional supplements
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !dipol
        if (any(gl%nPCSAFTD(gl%n_start:gl%n_end) /= 0)) then
            call calculate_PCSAFT_functionparts_DD_Trho(gl,T,D,getprevious)
            call calculate_PCSAFT_functionparts_x1_DD(gl,T,D,getprevious)
            call calculate_PCSAFT_functionparts_x2_DD(gl,T,D,getprevious)
            call ADDX2DERIVS(gl,T,D,get)
        end if
        
        !quadrupol        
        if (any(gl%nPCSAFTQ(gl%n_start:gl%n_end) /= 0))  then
            call calculate_PCSAFT_functionparts_QQ_Trho(gl,T,D,getprevious)
            call calculate_PCSAFT_functionparts_x1_QQ(gl,T,D,getprevious)
            call calculate_PCSAFT_functionparts_x2_QQ(gl,T,D,getprevious)
            call AQQX2DERIVS(gl,T,D,get)
        end if
        
        !association
        if (any(gl%kabPCSAFT(gl%n_start:gl%n_end) > 1.d-12)) then 
            call AASSOCX2(gl,T, D, get)
        end if
        
        ! XII. calculate the second composition derivatives of a_res
        do i = 1, nx2derivs
            if (get(i) .eq. 1) then
                gl%arx2_PCSAFT(i,1:gl%ncomp,1:gl%ncomp) = gl%ahcx2_PCSAFT(i,1:gl%ncomp,1:gl%ncomp) + gl%adispx2_PCSAFT(i,1:gl%ncomp,1:gl%ncomp) + gl%ADDx2_PCSAFTD(i,1:gl%ncomp,1:gl%ncomp) + gl%AQQx2_PCSAFTQ(i,1:gl%ncomp,1:gl%ncomp) + gl%AASSOCx2_PCSAFT(i,1:gl%ncomp,1:gl%ncomp)
            end if
        end do

        ! XIII. Fill the double spots of the symmetric matrix
        ! up to here, the symmetric part has not been calculated to save computing time
        ! note: is there a more elegant way to do this, as I iterate over all elements again; maybe transpose the matrix?
        do i = 1, nx2derivs
            if (get(i) .eq. 1) then ! the symmetric part is only copied for the derivatives, which have been calculated
                do xj = 1, gl%ncomp
                    do xi = 1, gl%ncomp
	                    if (xj .GT. xi) then
                            gl%arx2_PCSAFT(1:nx2derivs,xi,xj) = gl%arx2_PCSAFT(1:nx2derivs,xj,xi)
	                    end if
	                end do
                end do
            end if
        end do
        
    end if
    
    ! XIV. Change the T derivatives to tau derivatives
    call convert_T_derivs_x2(gl,getvector,gl%arx2_PCSAFT,x2)

!DEC$ END IF
end subroutine ARX2DERIVS

    
subroutine ARX3DERIVS(gl,T, D_in, getvector, x3)

! Henning Markgraf, June 2016
    
    ! a_res: residual part of the Helmholtz free energy
    ! defined by eq. A.3 in Gross, Sadowski 2001:
    ! a_res = a_hc + a_disp
    ! dependent on x_i and D [1/Angstrom] and T [K]
    
!--------------------------------------------------------------------------------------------------
! All third composition derivatives of the residual part of the Helmholtz free
! energy are calculated in this subroutine
!--------------------------------------------------------------------------------------------------
! INPUT PARAMETERS:     
! TEMPERATURE - T       K
! DENSITY     - D_in   mol/m^3
! get         - AN ARRAY WITH 15 ENTRIES WITH VALUES EITHER "1" OR "0", INDICATING WHICH DERIVATIVES ARE NEEDED:
!		only Nr. 1, 2 and 4 are calculated in this subroutine, the 15 entries are there only to grant
!		compliance with the other PC-SAFT routines
!               Derivatives number 1, 2 and 4 are derived by x_i, x_j and x_k and saved in output matrix x3(3,ncomp,ncomp,ncomp)
!                1. NORMALIZED RESIDUAL MIXTURE HELMHOLTZ ENERGY F AS A FUNCTION OF DEL AND TAU AND X 
!                2. 1ST DERIVATIVE OF F WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY DEL
!                4. 1ST DERIVATIVE OF F WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY TAU
! 
! OUTPUT PARAMETERS: 
! x3   - a matrix with dimensions (3,ncomp,ncomp,ncomp), WITH VALUES
!                   EITHER "0.d0" OR THE RESULTS OF THE DERIVATIVES AS INDICATED IN "GETDER"
!--------------------------------------------------------------------------------------------------







implicit none

    type(type_gl) :: gl

    
    ! I. Declarations
    !input
    double precision, intent (inout) :: T, D_in
    integer, dimension (nderivs), intent (inout) :: getvector
    integer :: nrsubst
    !output
    double precision, dimension (nx3derivs,gl%ncomp,gl%ncomp,gl%ncomp), intent (out) :: x3
    !working variables
    double precision :: D
    integer, dimension (nderivs) :: get
    integer, dimension (nderivs) :: getprevious
    integer :: i, xi, xj, xk, lowest, second, highest
    logical :: already_calculated
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    ! II. Initializations
    nrsubst = 0
    call allocate_arrays_PCSAFT(gl,nrsubst)

    if (gl%factor == 1.d3) then   !real fluids are calculated => unit conversion is necessary
    ! III. Change the units of the density from mol/m^3 to 1/Angstrom
        D = D_in*N_A/1.d30
    end if
    
    ! IV. Check if any derivatives of this T, rho, molfractions combination is already saved in the module variables
    get = getvector ! get can be changed if some derivatives have already been calculated
    
    ! Adding the necessary additional T derivatives required for the T to tau conversion:
    call add_T_conversion_derivs(gl,get)
    
    call check_calculated_PCSAFT(gl,T,D,nrsubst,gl%arx3_calculated_PCSAFT,get,already_calculated)
    if (any(get .EQ. 1)) then ! no need to calculate anything if get = 0
        
        ! V. some of the functions require all lower order derivatives, too
        call add_previous_derivs(gl,get,getprevious)
        call check_calculated_PCSAFT(gl,T,D,nrsubst,gl%arx3_calculated_PCSAFT,getprevious,already_calculated)
        
        ! VI. Initialize module variables if T, rho or molfractions different
        if (already_calculated .EQV. .FALSE.) then
            call init_derivs_PCSAFT(gl,nrsubst) ! all saved derivatives are deleted
        end if
        
        ! VII. update the module variables
        call update_modulevariables_PCSAFT(gl,T,D,nrsubst,gl%arx3_calculated_PCSAFT,get)
        
        ! VIII.  The T, rho derivatives of all the function parts, which A_disp and A_hc consists of
        !       are needed for the x3 derivatives
        ! i.e.: d_i, zeta_n, mmean, ab, g_ii, I_1, I_2, meo1, meo2, C_1
        call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
        ! and a_hs
        call AHSDERIVS(gl,get)
        
        ! IX. The x1 derivatives are also needed for the x3 derivatives
        ! i.e.: zeta_n_x, ab_x, gii_x, I_x, meo_x, C_x
        call calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)
        ! and ahs_x
        call AHSX1DERIVS(gl,get)
        
        ! X. The x2 derivatives are also needed for the x3 derivatives
        ! ab_xx, gii_xx, I_xx, meo_xx, C_xx
        call calculate_PCSAFT_functionparts_x2(gl,T,D,getprevious)
        ! and ahs_xx
        call AHSX2DERIVS(gl,get)
        
        ! XI. X3 DERIVATIVES:
        ! ab_xxx, gii_xxx, I_xxx, C_xxx
        call calculate_PCSAFT_functionparts_x3(gl,T,D,getprevious)
        ! and ahs_xx
        call AHSX3DERIVS(gl,get)
        
        ! XII. calculate ahc_xxx and adisp_xxx, which depend on the previously calculated functions
        call AHCX3DERIVS(gl,get)
        call ADISPX3DERIVS(gl,D,get)
    
        ! XIII. calculate the third composition derivatives of a_res
        do i = 1, 2
            if (get(i) .eq. 1) then
                gl%arx3_PCSAFT(i,1:gl%ncomp,1:gl%ncomp,1:gl%ncomp) = gl%ahcx3_PCSAFT(i,1:gl%ncomp,1:gl%ncomp,1:gl%ncomp) + gl%adispx3_PCSAFT(i,1:gl%ncomp,1:gl%ncomp,1:gl%ncomp)
            end if
        end do
        if (get(4) .EQ. 1) then ! problem is that in my x3-routines the arrays only have 3 spots, and deriv nr. 4 (d/dT) is on number 3
            gl%arx3_PCSAFT(3,1:gl%ncomp,1:gl%ncomp,1:gl%ncomp) = gl%ahcx3_PCSAFT(3,1:gl%ncomp,1:gl%ncomp,1:gl%ncomp) + gl%adispx3_PCSAFT(3,1:gl%ncomp,1:gl%ncomp,1:gl%ncomp)
        end if

        ! IX. Fill the double spots of the symmetric matrix
        ! up to here, the symmetric part has not been calculated to save computing time
        ! all double spots are just 0.d0 up to the calculation of a_res
        ! note: is there a more elegant way to do this, as I iterate over all elements again; maybe transpose the matrix?
        
        do xk = 1, gl%ncomp
            do xj = 1, gl%ncomp
	            do xi = 1, gl%ncomp
	                if ( xj .GT. xi .OR. xk .GT. xj ) then
		                ! order the integers to descending order (because e.g. derivative wrt components x3,x1,x4 is the same as derivative wrt x4,x3,x1)
		                highest = max(xi,xj,xk)
		                lowest = min(xi,xj,xk)
		                second = xi+xj+xk-highest-lowest
		                gl%arx3_PCSAFT(1:nx3derivs,xi,xj,xk) = gl%arx3_PCSAFT(1:nx3derivs,highest,second,lowest)
	                end if
	            end do
	        end do
        end do
        
    end if
    
    ! XIV. Change the T derivatives to tau derivatives
    call convert_T_derivs_x3(gl,getvector,gl%arx3_PCSAFT,x3)

!DEC$ END IF
end subroutine ARX3DERIVS



    end module pc_saft_ARX_derivs_module
