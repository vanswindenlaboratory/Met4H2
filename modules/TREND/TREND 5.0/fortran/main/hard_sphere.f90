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

    ! module for file hard_sphere.f90
    module hard_sphere_module
    !global use inclusion
    use module_all_types


    contains

	




    subroutine HSDERIVS (gl,TEMP, DENS, GETHSDER, HSDER, NRSUBST)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF DERIVATIVES OF THE RESIDUAL PART OF THE
    ! HARD SPHERE TERM BY M. Thol 2013/08
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    ! GETHSDER    - AN ARRAY WITH 10 ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. NORMALIZED RESIDUAL HELMHOLTZ ENERGY F AS A FUNCTION OF D AND T
    !                2. 1ST DERIVATIVE OF F WITH RESPECT TO DEL AT CONSTANT TAU, MULTIPLIED BY del
    !                3. 2ND DERIVATIVE OF F WITH RESPECT TO DEL AT CONSTANT TAU, MULTIPLIED BY del^2
    !                4. 1ST DERIVATIVE OF F WITH RESPECT TO TAU AT CONSTANT DEL, MULTIPLIED BY TAU
    !                5: 2ND DERIVATIVE OF F WITH RESPECT TO TAU AT CONSTANT DEL, MULTIPLIED BY TAU^2
    !                6: 2ND MIXED DERIVATIVE OF F WITH RESPECT TO DEL AND TAU, MULTIPLIED BY TAU*DEL
    !                7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO DEL, TAU, AND TAU, MULTIPLIED BY TAU*TAU*DEL
    !                8: 3RD DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^3
    !                9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !               10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    ! FLUID:      - INTEGER: THE NUMBER OF THE FLUID THE HELMHOLTZ ENERGY IS CALCULATED FOR IN THE ORDER GIVEN
    !                BY THE USER INPUT
    !
    ! OUTPUT PARAMETERS:
    ! HSDER       - AN ARRAY WITH 10 ENTRIES WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
    !               AS INDECATED IN "GETDER"
    !--------------------------------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION:: TAU, DEL, TEMP, DENS
    INTEGER:: NRSUBST
    DOUBLE PRECISION,DIMENSION(nderivshs):: HSDER, PACKFRACDER,PACKFRACDERNUM
    INTEGER, DIMENSION(nderivshs):: GETHSDER, GETFRAC,GETFRACPF
    DOUBLE PRECISION:: PF
    DOUBLE PRECISION:: HS,HSD,HSDD,HST,HSTT,HSTD,HSDTT,HSDDT,HSDDD,HSTTT,FH,DFHDK,D2FHDK,D3FHDK
    DOUBLE PRECISION:: DPFDD,DPF2D,DPFDT,DPFDTD, DPFDTT,DPFDTTD,DPF3D,DPFDTTT,DPFDDDT

    HSDER = (/0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0/)

    DEL=DENS/gl%RHOC(NRSUBST)
    TAU=gl%TC(NRSUBST)/TEMP

    !calculate packing fraction
    !GETFRAC = (/1,0,0,0,0,0,0,0,0,0/)   !Vector tells the subroutine which derivative is neccessary to calculate the property (0 or 1)

    PACKFRACDER=(/0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0/)
    GETFRACPF=(/1,0,0,0,0,0,0,0,0,0/)
    CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRACPF,PACKFRACDER,NRSUBST)
    PF=PACKFRACDER(1)


    !TN, Dezember 2018
    !For the calculation of virial coefficients, the density needs to be set to 0. Elsewise, the virials could be inaccurate due to rounding errors
    !IMPORTANT!!!:  Furthermore, as the zero density limit of the delta-derivatives of alpha_r are needed for the calculation of the virial coefficients,
    !   the derivatives must not be multiplied by delta, delta², delta3!!!
    !ATTENTION!!    --> If the density is 0, then the derivatives are not longer multiplied by delta, delta² and so on!!!
    !Unfortunately, the system of calculating the derivatives in a clever manner cannot be applied then anymore. Thus, all derivatives are programmed in a different way for
    !the evaluation of the virial coefficients (see the following else-statement)

    if (del > 1.D-14) then


        !Derivatives of the hard sphere term with respect to the packing fraction:
        FH=(4.D0*PF-3.D0*PF**2)/(1.D0-PF)**2
        DFHDK=(4.d0-2.d0*PF)/(1.d0-PF)**3
        D2FHDK=(10.d0-4.d0*PF)/(1.d0-PF)**4
        D3FHDK=(36.d0-12.d0*PF)/(1.d0-PF)**5


        !*****************************************************
        !   HS         HARD SPHERE TERM [-]
        !*****************************************************
        HS = 0.D0

        HS=FH

        HSDER(1)=HS


        !******************************************************************
        !   HSD         1ST DERIVATIVE OF THE HARD SPHERE TERM WRT DEL [-]
        !******************************************************************
        if (GETHSDER(2) == 1) then

            HSD = 0.D0
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETHSDER,PACKFRACDER,NRSUBST)
            DPFDD=PACKFRACDER(2)

            HSD=DFHDK*DPFDD

            HSDER(2)=HSD*DEL
        END IF


        !******************************************************************
        !   HSDD         2ND DERIVATIVE OF THE HARD SPHERE TERM WRT DEL [-]
        !******************************************************************
        if (GETHSDER(3) == 1) then
            HSDD = 0.D0

            GETFRAC = (/1,1,1,0,0,0,0,0,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)
            DPFDD=PACKFRACDER(2)
            DPF2D=PACKFRACDER(3)

            HSDD=D2FHDK*DPFDD**2+DFHDK*DPF2D

            HSDER(3)=HSDD*DEL**2
        END IF


        !******************************************************************
        !   HST         1ST DERIVATIVE OF THE HARD SPHERE TERM WRT TAU [-]
        !******************************************************************
        if (GETHSDER(4) == 1) then
            HST = 0.D0
            GETFRAC = (/1,0,0,1,0,0,0,0,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            DPFDT=PACKFRACDER(4)

            HST=DFHDK*DPFDT

            HSDER(4)=HST*TAU
        END IF

        !******************************************************************
        !   HSTT         2ND DERIVATIVE OF THE HARD SPHERE TERM WRT TAU [-]
        !******************************************************************
        if (GETHSDER(5) == 1) then
            HSTT = 0.D0

            GETFRAC = (/1,0,0,1,1,0,0,0,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            DPFDT=PACKFRACDER(4)
            DPFDTT=PACKFRACDER(5)

            HSTT = D2FHDK*DPFDT**2+DFHDK*DPFDTT

            HSDER(5)=HSTT*TAU**2
        END IF

        !********************************************************************************
        !   HSDT         2ND MIXED DERIVATIVE OF THE HARD SPHERE TERM WRT TAU AND DEL [-]
        !********************************************************************************
        if (GETHSDER(6) == 1) then
            HSTD = 0.D0
            GETFRAC = (/1,1,0,1,0,1,0,0,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            DPFDD=PACKFRACDER(2)
            DPFDT=PACKFRACDER(4)
            DPFDTD=PACKFRACDER(6)

            HSTD = D2FHDK*DPFDT*DPFDD+DFHDK*DPFDTD

            HSDER(6)=HSTD*DEL*TAU
        END IF

        !***********************************************************************************
        !   HSDTT         3RD MIXED DERIVATIVE OF THE HARD SPHERE TERM WRT TAU^2 AND DEL [-]
        !***********************************************************************************
        if (GETHSDER(7) == 1) then
            HSDTT = 0.D0
            GETFRAC = (/1,1,0,1,1,1,1,0,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            DPFDD=PACKFRACDER(2)
            DPFDT=PACKFRACDER(4)
            DPFDTT=PACKFRACDER(5)
            DPFDTD=PACKFRACDER(6)
            DPFDTTD=PACKFRACDER(7)

            HSDTT = D3FHDK*DPFDD*DPFDT**2+D2FHDK*DPFDTT*DPFDD+D2FHDK*DPFDT*DPFDT+DFHDK*DPFDTTD
            !HSDTT = D3FHDK*DPFDD*DPFDT**2+2.D0*D2FHDK*DPFDTD+D2FHDK*DPFDD*DPFDTT+DFHDK*DPFDTTD

            HSDER(7)=HSDTT*DEL*TAU**2
        END IF

        !******************************************************************
        !   HSDDD         3RD DERIVATIVE OF THE HARD SPHERE TERM WRT DEL [-]
        !******************************************************************
        if (GETHSDER(8) == 1) then
            HSDDD = 0.D0
            GETFRAC = (/1,1,1,0,0,0,0,1,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            DPFDD=PACKFRACDER(2)
            DPF2D=PACKFRACDER(3)
            DPF3D=PACKFRACDER(8)

            HSDDD = D3FHDK*DPFDD**3+2.D0*D2FHDK*DPF2D+D2FHDK*DPFDD*DPF2D+DFHDK*DPF3D

            HSDER(8)=HSDDD*DEL**3
        END IF

        !******************************************************************
        !   HSTTT         3RD DERIVATIVE OF THE HARD SPHERE TERM WRT TAU [-]
        !******************************************************************
        if (GETHSDER(9) == 1) then
            HSTTT = 0.D0
            GETFRAC = (/1,0,0,1,1,0,0,0,1,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            DPFDT=PACKFRACDER(4)
            DPFDTT=PACKFRACDER(5)
            DPFDTTT=PACKFRACDER(9)

            HSTTT = D3FHDK*DPFDT**3+2.D0*D2FHDK*DPFDTT*DPFDT+D2FHDK*DPFDT*DPFDTT+DFHDK*DPFDTTT

            HSDER(9)=HSTTT*TAU**3
        END IF

        !***********************************************************************************
        !   HSDDT         3RD MIXED DERIVATIVE OF THE HARD SPHERE TERM WRT TAU AND DEL^2 [-]
        !***********************************************************************************
        if (GETHSDER(10) == 1) then
            HSDDT = 0.D0
            GETFRAC = (/1,1,1,1,0,1,0,0,0,1/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            DPFDD=PACKFRACDER(2)
            DPF2D=PACKFRACDER(3)
            DPFDT=PACKFRACDER(4)
            DPFDTD=PACKFRACDER(6)
            DPFDDDT=PACKFRACDER(10)

            HSDDT = D3FHDK*DPFDT*DPFDD**2+2.D0*D2FHDK*DPFDTD+D2FHDK*DPFDT*DPF2D+DFHDK*DPFDDDT

            HSDER(10)=HSDDT*DEL**2*TAU

        END IF



    else ! Ab hier Ableitungen für delta = 0

        !Derivatives of the hard sphere term with respect to the packing fraction:
        FH= 0.d0
        DFHDK= 4.d0
        D2FHDK= 10.d0
        D3FHDK= 36.d0


        !*****************************************************
        !   HS         HARD SPHERE TERM [-]
        !*****************************************************
        HS = 0.D0

        HS=FH

        HSDER(1)=HS


        !******************************************************************
        !   HSD         1ST DERIVATIVE OF THE HARD SPHERE TERM WRT DEL [-]
        !******************************************************************
        if (GETHSDER(2) == 1) then

            HSD = 0.D0

            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETHSDER,PACKFRACDER,NRSUBST)
            DPFDD=PACKFRACDER(2)

            HSD=DFHDK*DPFDD

            HSDER(2)=HSD
        END IF


        !******************************************************************
        !   HSDD         2ND DERIVATIVE OF THE HARD SPHERE TERM WRT DEL [-]
        !******************************************************************
        if (GETHSDER(3) == 1) then
            HSDD = 0.D0

            GETFRAC = (/0,1,0,0,0,0,0,0,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)
            DPFDD=PACKFRACDER(2)

            HSDD=D2FHDK*DPFDD**2

            HSDER(3)=HSDD
        END IF


        !******************************************************************
        !   HST         1ST DERIVATIVE OF THE HARD SPHERE TERM WRT TAU [-]
        !******************************************************************
        if (GETHSDER(4) == 1) then

            HST = 0.D0

            HSDER(4)=HST
        END IF

        !******************************************************************
        !   HSTT         2ND DERIVATIVE OF THE HARD SPHERE TERM WRT TAU [-]
        !******************************************************************
        if (GETHSDER(5) == 1) then

            HSTT = 0.D0

            HSDER(5)=HSTT
        END IF

        !********************************************************************************
        !   HSDT         2ND MIXED DERIVATIVE OF THE HARD SPHERE TERM WRT TAU AND DEL [-]
        !********************************************************************************
        if (GETHSDER(6) == 1) then
            HSTD = 0.D0
            GETFRAC = (/0,0,0,0,0,1,0,0,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            DPFDTD=PACKFRACDER(6)

            HSTD = DFHDK*DPFDTD

            HSDER(6)=HSTD
        END IF

        !***********************************************************************************
        !   HSDTT         3RD MIXED DERIVATIVE OF THE HARD SPHERE TERM WRT TAU^2 AND DEL [-]
        !***********************************************************************************
        if (GETHSDER(7) == 1) then
            HSDTT = 0.D0
            GETFRAC = (/0,0,0,0,0,0,1,0,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            
            DPFDTTD=PACKFRACDER(7)

            HSDTT = DFHDK*DPFDTTD
            !HSDTT = D3FHDK*DPFDD*DPFDT**2+2.D0*D2FHDK*DPFDTD+D2FHDK*DPFDD*DPFDTT+DFHDK*DPFDTTD

            HSDER(7)=HSDTT
        END IF

        !******************************************************************
        !   HSDDD         3RD DERIVATIVE OF THE HARD SPHERE TERM WRT DEL [-]
        !******************************************************************
        if (GETHSDER(8) == 1) then
            HSDDD = 0.D0
            GETFRAC = (/0,1,0,0,0,0,0,0,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            DPFDD=PACKFRACDER(2)

            HSDDD = D3FHDK*DPFDD**3

            HSDER(8)=HSDDD
        END IF

        !******************************************************************
        !   HSTTT         3RD DERIVATIVE OF THE HARD SPHERE TERM WRT TAU [-]
        !******************************************************************
        if (GETHSDER(9) == 1) then

            HSTTT = 0.D0

            HSDER(9)=HSTTT
        END IF

        !***********************************************************************************
        !   HSDDT         3RD MIXED DERIVATIVE OF THE HARD SPHERE TERM WRT TAU AND DEL^2 [-]
        !***********************************************************************************
        if (GETHSDER(10) == 1) then
            HSDDT = 0.D0
            GETFRAC = (/0,1,0,0,0,1,0,0,0,0/)
            CALL PACKFRACKDERIVS(gl,TEMP,DENS,GETFRAC,PACKFRACDER,NRSUBST)

            
            DPFDTD=PACKFRACDER(6)
            

            HSDDT = 2.D0*D2FHDK*DPFDTD

            HSDER(10)=HSDDT

        END IF


    END IF

    END SUBROUTINE





    subroutine PACKFRACKDERIVS (gl,T, D, GETPFDER, PFDER, NRSUBST)
    !**************************************************************************
    !--------------------------------------------------------------------------------------------------
    ! SUBROUTINE FOR THE CALCULATION OF DERIVATIVES OF THE RESIDUAL PART OF THE
    ! HARD SPHERE TERM BY M. Thol 2013/08
    !--------------------------------------------------------------------------------------------------
    ! INPUT PARAMETERS:
    ! TEMPERATURE - T   K
    ! DENSITY     - D   MOL/M^3
    ! GETPFDER    - AN ARRAY WITH 10 ENTRIES WITH VALUES EITHER "1" OR "0",
    !                INDICATING WHICH DERIVATIVES ARE NEEDED:
    !                1. NORMALIZED RESIDUAL HELMHOLTZ ENERGY F AS A FUNCTION OF D AND T
    !                2. 1ST DERIVATIVE OF F WITH RESPECT TO DEL AT CONSTANT TAU, MULTIPLIED BY del
    !                3. 2ND DERIVATIVE OF F WITH RESPECT TO DEL AT CONSTANT TAU, MULTIPLIED BY del^2
    !                4. 1ST DERIVATIVE OF F WITH RESPECT TO TAU AT CONSTANT DEL, MULTIPLIED BY TAU
    !                5: 2ND DERIVATIVE OF F WITH RESPECT TO TAU AT CONSTANT DEL, MULTIPLIED BY TAU^2
    !                6: 2ND MIXED DERIVATIVE OF F WITH RESPECT TO DEL AND TAU, MULTIPLIED BY TAU*DEL
    !                7: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO DEL, TAU, AND TAU, MULTIPLIED BY TAU*TAU*DEL
    !                8: 3RD DERIVATIVE OF F WITH RESPECT TO DEL, MULTIPLIED BY DEL^3
    !                9: 3RD DERIVATIVE OF F WITH RESPECT TO TAU, MULTIPLIED BY TAU^3
    !               10: 3RD MIXED DERIVATIVE OF F WITH RESPECT TO TAU, DEL, AND DEL, MULTIPLIED BY TAU*DEL*DEL
    ! FLUID:      - INTEGER: THE NUMBER OF THE FLUID THE HELMHOLTZ ENERGY IS CALCULATED FOR IN THE ORDER GIVEN
    !                BY THE USER INPUT
    !
    ! OUTPUT PARAMETERS:
    ! PFDER       - AN ARRAY WITH 10 ENTRIES WITH VALUES EITHER "0" OR THE RESULTS OF THE DERIVATIVES
    !               AS INDECATED IN "GETPFDER"
    !--------------------------------------------------------------------------------------------------






    implicit none

    type(type_gl) :: gl


    DOUBLE PRECISION:: TAU, DEL, T, D, XDENOM, H134, ENUM1,ENUM2
    DOUBLE PRECISION:: NUMER, DENOM, DNUMERDT, DDENOMDT, NUMER1, NUMER2, NUMER3
    INTEGER:: NRSUBST
    DOUBLE PRECISION:: PF,PFDD,PFT,PFTT,PFTD,PFTTD,PFDDT,PFDDD,PFTTT
    DOUBLE PRECISION,DIMENSION(nderivshs):: PFDER
    INTEGER, DIMENSION(nderivshs):: GETPFDER


    DEL=D/gl%RHOC(NRSUBST)
    TAU=gl%TC(NRSUBST)/T

    H134 = gl%eos_coeff%h1(NRSUBST)*gl%eos_coeff%h3(NRSUBST)*gl%eos_coeff%h4(NRSUBST)
    XDENOM=gl%eos_coeff%h2(NRSUBST)+gl%eos_coeff%h3(NRSUBST)*TAU**gl%eos_coeff%h4(NRSUBST)

    if (del > 1.d-14) then
        !*****************************************************
        !   PF         PACKING FRACTION [-]
        !*****************************************************
        PF = 0.D0

        PF=gl%eos_coeff%h1(NRSUBST)*DEL/XDENOM

        PFDER(1)=PF


        !***********************************************************************
        !   PFD         1ST DERIVATIVE OF THE PACKING FRACTION WRT DEL * DEL [-]
        !***********************************************************************
        IF (GETPFDER(2) == 1) THEN

            PFDER(2)=PF/DEL

        END IF


        !***************************************************************************
        !   PFDD         2ND DERIVATIVE OF THE PACKING FRACTION WRT DEL * DEL**2 [-]
        !***************************************************************************
        IF (GETPFDER(3) == 1) THEN
            PFDD = 0.D0
        END IF

        !***********************************************************************
        !   PFT         1ST DERIVATIVE OF THE PACKING FRACTION WRT TAU * TAU [-]
        !***********************************************************************
        IF (GETPFDER(4) == 1) THEN
            PFT = 0.D0
            PFT = -H134*TAU**(gl%eos_coeff%h4(NRSUBST)-1.D0)*DEL/XDENOM**2

            PFDER(4)=PFT

        END IF

        !***************************************************************************
        !   PFTT         2ND DERIVATIVE OF THE PACKING FRACTION WRT TAU * TAU^2 [-]
        !***************************************************************************
        IF (GETPFDER(5) == 1) THEN
            PFTT = 0.D0
            ENUM1 = gl%eos_coeff%h2(NRSUBST)*(1.D0-gl%eos_coeff%h4(NRSUBST))*TAU**(gl%eos_coeff%h4(NRSUBST)-2.D0)
            ENUM2 = (1.D0+gl%eos_coeff%h4(NRSUBST))*gl%eos_coeff%h3(NRSUBST)*TAU**(2.D0*gl%eos_coeff%h4(NRSUBST)-2.D0)

            PFTT = DEL*H134*(ENUM1+ENUM2)/XDENOM**3

            PFDER(5)=PFTT
        END IF

        !********************************************************************************************
        !   PFTD         2ND MIXED DERIVATIVE OF THE PACKING FRACTION WRT TAU AND DEL * TAU * DEL [-]
        !********************************************************************************************
        IF (GETPFDER(6) == 1) THEN
            PFTD = 0.D0
            PFTD = -H134*TAU**(gl%eos_coeff%h4(NRSUBST)-1.D0)/XDENOM**2

            PFDER(6)=PFTD

        END IF

        !*************************************************************************************************
        !   PFTTD         3RD MIXED DERIVATIVE OF THE PACKING FRACTION WRT TAU^2 AND DEL * TAU^2 * DEL [-]
        !*************************************************************************************************
        IF (GETPFDER(7) == 1) THEN
            PFTTD = 0.D0
            ENUM1 = gl%eos_coeff%h2(NRSUBST)*(1.D0-gl%eos_coeff%h4(NRSUBST))*TAU**(gl%eos_coeff%h4(NRSUBST)-2.D0)
            ENUM2 = (1.D0+gl%eos_coeff%h4(NRSUBST))*gl%eos_coeff%h3(NRSUBST)*TAU**(2.D0*gl%eos_coeff%h4(NRSUBST)-2.D0)

            PFTTD = H134*(ENUM1+ENUM2)/XDENOM**3

            PFDER(7)=PFTTD
        END IF

        !***************************************************************************
        !   PFDDD         3RD DERIVATIVE OF THE PACKING FRACTION WRT DEL * DEL^3 [-]
        !***************************************************************************
        IF (GETPFDER(8) == 1) THEN
            PFDDD = 0.D0

            PFDER(8)=PFDDD
        END IF

        !**************************************************************************
        !   PFTTT         3RD DERIVATIVE OF THE PACKING FRACTION WRT TAU * TAU^3[-]
        !**************************************************************************
        IF (GETPFDER(9) == 1) THEN
            PFTTT = 0.D0
            NUMER = 0.D0
            NUMER1 = 0.D0
            NUMER2 = 0.D0
            NUMER3 = 0.D0
            DENOM = 0.D0
            DNUMERDT = 0.D0
            DDENOMDT = 0.D0

            !------------------------------------------------------------------------------------------------
            !DERIVE
            !NUMER1 = H3(NRSUBST)**2*TAU**(2.d0*H4(NRSUBST))*(H4(NRSUBST)+1.D0)*(H4(NRSUBST)+2.D0)
            !NUMER2 = 4.D0*H2(NRSUBST)*H3(NRSUBST)*TAU**H4(NRSUBST)*(H4(NRSUBST)+1.d0)*(1.D0-H4(NRSUBST))
            !NUMER3 = H2(NRSUBST)**2*(H4(NRSUBST)-1.d0)*(H4(NRSUBST)-2.d0)
            !
            !PFTTT = H134*DEL*TAU**(H4(nrsubst)-3.D0)*(NUMER1+NUMER2+NUMER3)/XDENOM**4

            !------------------------------------------------------------------------------------------------
            !Mirco
            NUMER1 = gl%eos_coeff%h3(NRSUBST)**2*TAU**(2.d0*gl%eos_coeff%h4(NRSUBST))*(-gl%eos_coeff%h4(NRSUBST)-1.D0)*(gl%eos_coeff%h4(NRSUBST)+2.D0)
            NUMER2 = 4.D0*gl%eos_coeff%h2(NRSUBST)*gl%eos_coeff%h3(NRSUBST)*TAU**gl%eos_coeff%h4(NRSUBST)*(gl%eos_coeff%h4(NRSUBST)**2-1.d0)
            NUMER3 = gl%eos_coeff%h2(NRSUBST)**2*(1.D0-gl%eos_coeff%h4(NRSUBST))*(gl%eos_coeff%h4(NRSUBST)-2.d0)

            PFTTT = H134*DEL*TAU**(gl%eos_coeff%h4(nrsubst)-3.D0)*(NUMER1+NUMER2+NUMER3)/XDENOM**4

            PFDER(9)=PFTTT
        END IF

        !*************************************************************************************************
        !   PFDDT         3RD MIXED DERIVATIVE OF THE PACKING FRACTION WRT TAU AND DEL^2 * TAU * DEL^2 [-]
        !*************************************************************************************************
        IF (GETPFDER(10) == 1) THEN
            PFDDT = 0.D0

            PFDER(10)=PFDDT
        END IF
        
        

    else ! Ab hier Packing Fraction für delta = 0

        !TN: Nur PFD,PFTD,PFTTD sind nicht null. Daher wurden alle anderen Ableitungen weggelassen

        !***********************************************************************
        !   PFD         1ST DERIVATIVE OF THE PACKING FRACTION WRT DEL * DEL [-]
        !***********************************************************************
        IF (GETPFDER(2) == 1) THEN

            PFDER(2)=gl%eos_coeff%h1(NRSUBST)/XDENOM

        END IF


        
        !********************************************************************************************
        !   PFTD         2ND MIXED DERIVATIVE OF THE PACKING FRACTION WRT TAU AND DEL * TAU * DEL [-]
        !********************************************************************************************
        IF (GETPFDER(6) == 1) THEN
            PFTD = 0.D0
            PFTD = -H134*TAU**(gl%eos_coeff%h4(NRSUBST)-1.D0)/XDENOM**2

            PFDER(6)=PFTD

        END IF

        !*************************************************************************************************
        !   PFTTD         3RD MIXED DERIVATIVE OF THE PACKING FRACTION WRT TAU^2 AND DEL * TAU^2 * DEL [-]
        !*************************************************************************************************
        IF (GETPFDER(7) == 1) THEN
            PFTTD = 0.D0
            ENUM1 = gl%eos_coeff%h2(NRSUBST)*(1.D0-gl%eos_coeff%h4(NRSUBST))*TAU**(gl%eos_coeff%h4(NRSUBST)-2.D0)
            ENUM2 = (1.D0+gl%eos_coeff%h4(NRSUBST))*gl%eos_coeff%h3(NRSUBST)*TAU**(2.D0*gl%eos_coeff%h4(NRSUBST)-2.D0)

            PFTTD = H134*(ENUM1+ENUM2)/XDENOM**3

            PFDER(7)=PFTTD
        END IF
        

    end if

    END SUBROUTINE


    end module hard_sphere_module
