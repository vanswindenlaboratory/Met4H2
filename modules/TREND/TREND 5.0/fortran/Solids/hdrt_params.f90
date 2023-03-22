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

    ! module for file hdrt_params.f90
    module hdrt_params_module
    !global use inclusion
    use module_all_types


    contains



    
! Vaclav - October + November 2016 (new parameters for the mixed hydrates model)
    
    ! universal reference state conditions T0 = T0a = 273.15 K, p0 = p0a = 1 Pa
    ! new lattice parameter correlation - conference article EFM 2016
    ! REOS - Reference EoS (not GERG)
    ! BS - Langmuir from - multilayered water shell by Ballard & Sloan
    ! M4 - compressibility from Murnaghan EoS
    ! CONS - constant shell radii for all hydrate formers
    ! constant gw0 and hw0 for both sI and sII hydrates
    
!******************************************************************************
subroutine hdrt_params_argon(gl,J,hdrt_structure_flag)
!******************************************************************************
!  Definition of the argon-hydrate properties




implicit none

    type(type_gl) :: gl


! Input arguments
integer, intent(in) :: J, hdrt_structure_flag

    !Pure Hydrate? -> structure is known
    !if (N_guests == 1) hdrt_structure = 's2'

    ! Kihara potential parameters
    ! ---------------------------
    gl%a_hc(J) = 0.1909d0     ! [A] hard core radius ... Our 2nd virial coef.
    !! gw_B0 and hw_B0 CONSTANT for sII from 'z_Kiharas_3_phase_HdrNr_2016.xlsx'
    !    gw_B0 = 904.887d0      ! Gibbs energy of water in empty lattice [J/mol] 
    !    hw_B0 = -5044.369d0    ! Enthalpy of empty lattice [J.mol^-1]
    ! Kihara parameters fitted to CPO and hdrt_composition data with uncer  
    gl%sigma(J) =3.01185d0!!3.01537d0!(fitted for double occ)     3.02306d0!VHIw!3.01190d0!VLwH !3.01659d0!5VHIw  !      3.08714 params single occ Andy      !3.1353d0!3.24495d0!3.1164d0!3.12611d0
    gl%eps_k(J) =141.176d0!!141.175d0!(fitted for double occ)     141.172d0!VHIw!141.172d0!VLwH !141.169d0!5VHIw  !       136.701 params single occ Andy      141.074d0!123.1565d0 !124.305d0!
    
    gl%sigmad(J) =2.9953d0!!2.99704d0!!(fitted for double occ)   params sec vir Andy       3.00047d0!VHIw!2.99842d0!VLwH !2.99770d0!5VHIw  !     2.92552d0!2.9805d0!3.529d0     ! taken from Abolala_j2015
    gl%epsd_k(J) = 148.520d0!!148.519d0!(fitted for double occ)  params sec vir Andy      148.516d0!VHIw!148.518d0 !VLwH!148.518d0!5VHIw !     !124.29d0
    
    gl%r_vdw(J) =1.7144d0!!1.7173d0!1.71733d0!(fitted for double occ)   !1.88d0 param lit       1.72076d0!VHIw!1.71857d0 !VLwH !1.71682d0!5VHIw  !                 !1.5798d0!1.88d0
    !//////////////////////////////////////////////////////////////  

if (hdrt_structure_flag==1) then ! sI hydrate   
    
    ! a0_hdrt need to be fitted
    gl%a0_hdrt(J) = 12.d0         ! [A] lattice param @ p0_a = 1 Pa, T0_a = 273.15 K
    gl%B0_hdrt(1,J) = 1.d10       ! [Pa] Bulk modulus
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/0.d0, 0.d0, 0.d0/)
        gl%KS_l(:,J) = (/0.d0, 0.d0, 0.d0/)
!DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) '! Error no T-fit for CiJ for argon'
!DEC$ END IF ! WO_WRITE_TO_CONSOLE
    end if
        
elseif (hdrt_structure_flag==2) then ! sII hydrate  
    
    gl%a0_hdrt(J) = 17.232725d0   ! [A] lattice param @ p0_a = 1 Pa, T0_a = 273.15 K - October 2016
    gl%B0_hdrt(1,J) = 1.d10        ! [Pa] Bulk modulus    
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/0.d0, 0.d0, 0.d0/)
        gl%KS_l(:,J) = (/0.d0, 0.d0, 0.d0/)
!DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) '! Error no T-fit for CiJ for argon'
!DEC$ END IF ! WO_WRITE_TO_CONSOLE
    end if

end if
    
end subroutine hdrt_params_argon
!******************************************************************************
!******************************************************************************



!******************************************************************************
subroutine hdrt_params_co(gl,J,hdrt_structure_flag)
!******************************************************************************
!  Definition of the co-hydrate properties




implicit none

    type(type_gl) :: gl


! Input arguments
integer, intent(in) :: J, hdrt_structure_flag

    !Pure Hydrate? -> structure is known
    !if (N_guests == 1) hdrt_structure = 's1'

    ! Kihara potential parameters
    ! ---------------------------
    gl%a_hc(J) = 0.4385d0       ! [A] hard core radius ... Our 2nd virial coef.
    !! gw_B0 and hw_B0 CONSTANT for sI from 'z_Kiharas_3_phase_HdrNr_2016.xlsx'
    !    gw_B0 = 987.030d0     ! Gibbs energy of water in empty lattice [J/mol] 
    !    hw_B0 = -4790.157d0    ! Enthalpy of empty lattice [J.mol^-1]
    !! Kihara parameters fitted to CPO and hdrt_composition data with uncer  
    !    sigma = 3.07688d0     ! sigma taken from 2 virial. coef. - paper I 2016
    !    eps_k = 132.864d0     ! only eps_k correlated!
    ! Kihara parameters - in second step refitted only sigma  
    gl%sigma(J) = 3.07376d0     
    gl%eps_k(J) = 132.864d0 
    
    gl%sigmad(J) = 2.58937d0
    gl%epsd_k(J) = 163.842d0
    
    gl%r_vdw(J) = 1.5697d0
    
if (hdrt_structure_flag==1) then ! sI hydrate

    ! a0_hdrt need to be fitted
    gl%a0_hdrt(J) = 12.017111d0   ! [A] lattice param @ p0_a = 1 Pa, T0_a = 273.15 K - October 2016
    gl%B0_hdrt(1,J) = 1.d10         ! [Pa] Bulk modulus           
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/0.d0, 0.d0, 0.d0/)
        gl%KS_l(:,J) = (/0.d0, 0.d0, 0.d0/)
!DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) '! Error no T-fit for CiJ for CO'
!DEC$ END IF ! WO_WRITE_TO_CONSOLE
    end if
    
elseif (hdrt_structure_flag==2) then ! sII hydrate  
    
    ! a0_hdrt need to be fitted
    gl%a0_hdrt(J) = 17.3d0      ! [A] lattice param @ p0_a = 1 Pa, T0_a = 273.15 K 
    gl%B0_hdrt(1,J) = 1.d10     ! [Pa] Bulk modulus    
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/0.d0, 0.d0, 0.d0/)
        gl%KS_l(:,J) = (/0.d0, 0.d0, 0.d0/)
!DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) '! Error no T-fit for CiJ for CO'
!DEC$ END IF ! WO_WRITE_TO_CONSOLE
    end if
    
end if         

end subroutine hdrt_params_co
!******************************************************************************
!******************************************************************************




!******************************************************************************
subroutine hdrt_params_co2(gl,J,hdrt_structure_flag)
!******************************************************************************
!  Definition of the co2-hydrate properties




implicit none

    type(type_gl) :: gl


! Input arguments
integer, intent(in) :: J, hdrt_structure_flag

    !Pure Hydrate? -> structure is known
    !if (N_guests == 1) hdrt_structure = 's1'
    
    ! Kihara potential parameters
    ! ---------------------------
    gl%a_hc(J) = 0.6805d0      ! [A] hard core radius (considered as constant) ... Tee et al. 1966
    !! gw_B0 and hw_B0 CONSTANT for sI from 'z_Kiharas_3_phase_HdrNr_2016.xlsx'
    !    gw_B0 = 987.030d0     ! Gibbs energy of water in empty lattice [J/mol] 
    !    hw_B0 = -4790.157d0    ! Enthalpy of empty lattice [J.mol^-1]
    ! Kihara parameters fitted to CPO and hdrt_composition data with uncer  
    gl%sigma(J) = 2.99447d0!0.298925967289d01!
    gl%eps_k(J) = 175.739d0
    
    gl%sigmad(J) = 2.24269d0
    gl%epsd_k(J) = 443.132d0
    
    gl%r_vdw(J) = 1.675d0
    
if (hdrt_structure_flag==1) then ! sI hydrate
      
    gl%a0_hdrt(J) = 11.979675d0   ! [A] lattice param @ p0_a = 1 Pa, T0_a = 273.15 K - October 2016  
    gl%B0_hdrt(1,J) = 10.d9         ! [Pa] Bulk modulus = 10 GPa 
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/-24.9824d0, 2743.7375d0, 31948.6496d0/)
        gl%KS_l(:,J) = (/-22.4037d0, 3171.7604d0, 0.d0/)
    end if    
    
elseif (hdrt_structure_flag==2) then ! sII hydrate  
    
    ! a0_hdrt need to be fitted
    gl%a0_hdrt(J) = 17.3d0      ! [A] lattice param @ p0_a = 0 Pa, T0_a = 298.15 K
    gl%B0_hdrt(1,J) = 1.d10     ! [Pa] Bulk modulus    

    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/-25.1752d0, 3089.4741d0, 48259.6778d0/)
        gl%KS_l(:,J) = (/-21.0917d0, 2405.3662d0, 28783.0000d0/)
    end if
        
end if

end subroutine hdrt_params_co2
!******************************************************************************
!******************************************************************************



!******************************************************************************
subroutine hdrt_params_ethane(gl,J, hdrt_structure_flag)
!******************************************************************************
!  Definition of the ethane-hydrate properties
!  Note: VLwLe equilibrium is not predicted well by fluid EoSs




implicit none

    type(type_gl) :: gl


! Input arguments
integer, intent(in) :: J, hdrt_structure_flag
    
    !Pure Hydrate? -> structure is known
    !if (N_guests == 1) hdrt_structure = 's1'

    ! Kihara potential parameters
    ! ---------------------------
    gl%a_hc(J) = 0.5651d0      ! [A] hard core radius ... Tee et al. 1966
    !! gw_B0 and hw_B0 CONSTANT for sI from 'z_Kiharas_3_phase_HdrNr_2016.xlsx'
    !    gw_B0 = 987.030d0     ! Gibbs energy of water in empty lattice [J/mol] 
    !    hw_B0 = -4790.157d0    ! Enthalpy of empty lattice [J.mol^-1]
    !! Kihara parameters fitted to CPO and hdrt_composition data with uncer  
    gl%sigma(J) = 3.20797d0
    gl%eps_k(J) = 186.227d0
    
    gl%sigmad(J) = 3.404131d0
    gl%epsd_k(J) = 331.697d0
    
    gl%r_vdw(J) = 1.8679d0

if (hdrt_structure_flag==1) then ! sI hydrate
     
    gl%a0_hdrt(J) = 12.067015d0   ! [A] lattice param @ p0_a = 1 Pa, T0_a = 273.15 K - October 2016
    gl%B0_hdrt(1,J) = 10.d9         ! [Pa] Bulk modulus = 10 GPa
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/-23.1806d0,   13.7469d0, 65052.4158d0/)
        gl%KS_l(:,J) = (/-23.7290d0, 3843.2773d0,  8882.4254d0/)
    end if
    
elseif (hdrt_structure_flag==2) then ! sII hydrate  
    
    ! a0_hdrt need to be fitted
    gl%a0_hdrt(J) = 17.3d0      ! [A] lattice param @ p0_a = 0 Pa, T0_a = 298.15 K
    gl%B0_hdrt(1,J) = 1.d10     ! [Pa] Bulk modulus    
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/-24.8378d0,  926.9897d0, 43614.8526d0/)
        gl%KS_l(:,J) = (/-22.2291d0, 3534.8896d0, -3371.3000d0/)
    end if
    
end if

end subroutine hdrt_params_ethane
!******************************************************************************
!******************************************************************************




!******************************************************************************
subroutine hdrt_params_methane(gl,J,hdrt_structure_flag)
!******************************************************************************
!  Definition of the methane-hydrate properties




implicit none

    type(type_gl) :: gl


! Input arguments
integer, intent(in) :: J, hdrt_structure_flag

    !Pure Hydrate? -> structure is known
    !if (N_guests == 1) hdrt_structure = 's1'

    ! Kihara potential parameters
    ! ---------------------------
    gl%a_hc(J) = 0.3834d0      ! [A] hard core radius ... Tee et al. 1966
    !! gw_B0 and hw_B0 CONSTANT for sI from 'z_Kiharas_3_phase_HdrNr_2016.xlsx'
    !    gw_B0 = 987.030d0     ! Gibbs energy of water in empty lattice [J/mol] 
    !    hw_B0 = -4790.157d0    ! Enthalpy of empty lattice [J.mol^-1]
    !! Kihara parameters fitted to CPO and hdrt_composition data with uncer  
        gl%sigma(J) = 3.22061d0
        gl%eps_k(J) = 154.293d0
        
        gl%sigmad(J) = 3.134344d0
        gl%epsd_k(J) = 194.468d0
        
        gl%r_vdw(J) = 1.5982d0
    
if (hdrt_structure_flag==1) then ! sI hydrate
     
    gl%a0_hdrt(J) = 11.970005d0   ! [A] lattice param @ p0_a = 1 Pa, T0_a = 273.15 K - October 2016
    gl%B0_hdrt(1,J) = 10.d9         ! [Pa] Bulk modulus = 10 GPa 
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/-23.6453d0, 2714.5643d0, 0.d0/)
        gl%KS_l(:,J) = (/-22.0651d0, 2760.1604d0, 0.d0/)
    end if
    
elseif (hdrt_structure_flag==2) then ! sII hydrate  
    
    ! a0_hdrt need to be fitted
    gl%a0_hdrt(J) = 17.3d0      ! [A] lattice param @ p0_a = 0 Pa, T0_a = 298.15 K
    gl%B0_hdrt(1,J) = 1.d10     ! [Pa] Bulk modulus  
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/-23.5746d0, 2708.8070d0, 0.d0/)
        gl%KS_l(:,J) = (/-20.6991d0, 2147.6899d0, -12013.6211d0/)
    end if
    
end if

end subroutine hdrt_params_methane
!******************************************************************************
!******************************************************************************




!******************************************************************************
subroutine hdrt_params_nitrogen(gl,J,hdrt_structure_flag)
!******************************************************************************
!  Definition of the nitrogen-hydrate properties




implicit none

    type(type_gl) :: gl


! Input arguments
integer, intent(in) :: J, hdrt_structure_flag

    !Pure Hydrate? -> structure is known
    !if (N_guests == 1) hdrt_structure = 's2'

    ! Kihara potential parameters
    ! ---------------------------
    gl%a_hc(J) = 0.3526d0      ! [A] hard core radius ... Sherwood and Prausnitz 1964 
    !! gw_B0 and hw_B0 CONSTANT for sII from 'z_Kiharas_3_phase_HdrNr_2016.xlsx'
    !    gw_B0 = 904.887d0      ! Gibbs energy of water in empty lattice [J/mol] 
    !    hw_B0 = -5044.369d0    ! Enthalpy of empty lattice [J.mol^-1]
    !! Kihara parameters fitted to CPO and hdrt_composition data with uncer  
        gl%sigma(J) = 2.9978d0!!3.04964d0!(fitted for double occ)  3.04610d0!nur VLwH3.12517d0 !nur VHIw3.04610d0!nur VLwH   3.12337 params single occ Andy      !3.08340d0!3.12337d0! Param Jäger        !3.0701d0!3.107066d0!3.1408d0!3.15185d0   2.99163d0!3.12337d0!3.056304d0!
        gl%eps_k(J) = 128.308d0!!128.37d0!127.028d0!(fitted for double occ)    128.368  params single occ Andy  127.037d0!nur VLwH 127.130d0!nur VHIw127.037d0!nur VLwH    !127.030d0!128.368d0! Param Jäger        !126.4523d0!127.0415d0!128.203d0        128.5788d0!128.368d0!126.981d0!
        
        gl%sigmad(J) = 2.4621d0!!2.41249d0!(fitted for double occ)     params sec vir Andy   2.41250d0! nur VLWH 2.31068d0!nur VHIw2.41250d0! nur VLWH   !2.13352d0!3.035315d0!2.35278d0!2.3178d0!  2.7d0
        gl%epsd_k(J) = 141.737d0!!141.787d0!(fitted for double occ)     params sec vir Andy   141.179d0!nur VLwH 141.785d0!nur VHIw141.179d0!nur VLwH   !141.791d0!157.995d0!0.419382d3!255.9403d0!141.8d0!141.791d0     !Sherwood & Prausnitz  !taken from Rasoolzadeh j2017, will be fitted! 126.981d0 !
        
        
    gl%r_vdw(J) =1.6181d0!!1.56844d0!1.56844d0!(fitted for double occ)     1.55 param lit  1.57196d0!nur VLwH1.46663d0!nur VHIw1.57196d0!nur VLwH    !1.50191!1.5083d0!1.48d0 !1.55d0           ![A] distance between two molecules (double occupancy)


if (hdrt_structure_flag==1) then ! sI hydrate   
    
    ! need to be fitted
    gl%a0_hdrt(J) = 11.7d0         ! [A] lattice param @ p0_a = 0 Pa, T0_a = 275.15 K
    gl%B0_hdrt(1,J) = 1.d10       ! [Pa] Bulk modulus
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/-23.0646d0, 2475.8673d0, 0.d0/)
        gl%KS_l(:,J) = (/-21.8424d0, 2337.1765d0, 0.d0/)
    end if
       
elseif (hdrt_structure_flag==2) then ! sI hydrate
     
    gl%a0_hdrt(J) = 17.275606d0   ! [A] lattice param @ p0_a = 1 Pa, T0_a = 273.15 K - October 2016
    gl%B0_hdrt(1,J) = 10.d9         ! [Pa] Bulk modulus = 10 GPa

    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/-22.9726d0, 2499.2232d0, 0.d0/)
        gl%KS_l(:,J) = (/-20.6160d0, 2033.6043d0, -12672.192d0/)
    end if
    
end if

end subroutine hdrt_params_nitrogen
!******************************************************************************
!******************************************************************************




!******************************************************************************
subroutine hdrt_params_oxygen(gl,J, hdrt_structure_flag)
!******************************************************************************
!  Definition of the oxygen-hydrate properties




implicit none

    type(type_gl) :: gl


! Input arguments
integer, intent(in) :: J, hdrt_structure_flag

    !Pure Hydrate? -> structure is known
    !if (N_guests == 1) hdrt_structure = 's2'

    ! Kihara potential parameters
    ! ---------------------------
    gl%a_hc(J) = 0.22980d0     ! [A] hard core radius ... Our 2nd virial coef.
    !! gw_B0 and hw_B0 CONSTANT for sII from 'z_Kiharas_3_phase_HdrNr_2016.xlsx'
    !    gw_B0 = 904.887d0      ! Gibbs energy of water in empty lattice [J/mol] 
    !    hw_B0 = -5044.369d0    ! Enthalpy of empty lattice [J.mol^-1]
    !! Kihara parameters fitted to CPO and hdrt_composition data with uncer  
    gl%sigma(J) = 3.1647d0!!3.17740d0!(fitted for double occ)    3.28373 params single occ Andy      !3.214567d0!3.2954304d0!3.295027d0!3.244906d0!3.28373d0    !3.3292d0!3.29938d0!3.29210d0
    gl%eps_k(J) =128.926d0!!128.940d0!(fitted for double occ)    133.338 params single occ Andy      !128.077d0!132.226d0 !130.26278d0 !133.338d0    !132.175d0
    
    gl%sigmad(J) = 2.4028d0!!2.39580d0!(fitted for double occ)  2.90318d0   params sec vir Andy      !2.29887d0!2.5002745d0!2.58963d0!2.90318d0!(sec vir)    !2.16009d0!2.1563d0!2.68601d0
    gl%epsd_k(J) = 152.961d0!!152.968d0!(fitted for double occ)  152.963d0   params sec vir Andy       !(sec vir)      !289.528d0!167.619d0       !taken from rasoolzadeh j2017, will be fitted
    
    gl%r_vdw(J) =1.52204d0!1.51574d0!(fitted for double occ)    1.52 param lit      !1.413625d0!1.5d0!1.3712d0 !1.3736d0!1,44 !1.52d0


if (hdrt_structure_flag==1) then ! sI hydrate   
    
    ! need to be fitted
    gl%a0_hdrt(J) = 12.d0         ! [A] lattice param @ p0_a = 0 Pa, T0_a = 275.15 K
    gl%B0_hdrt(1,J) = 1.d10       ! [Pa] Bulk modulus
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/0.d0, 0.d0, 0.d0/)
        gl%KS_l(:,J) = (/0.d0, 0.d0, 0.d0/)
!DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) '! Error no T-fit for CiJ for oxygen'
!DEC$ END IF ! WO_WRITE_TO_CONSOLE
    end if
    
elseif (hdrt_structure_flag==2) then ! sI hydrate
    
    gl%a0_hdrt(J) = 17.262769d0   ! [A] lattice param @ p0_a = 1 Pa, T0_a = 273.15 K - October 2016
    gl%B0_hdrt(1,J) = 10.d9         ! [Pa] Bulk modulus = 10 GPa

    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/0.d0, 0.d0, 0.d0/)
        gl%KS_l(:,J) = (/0.d0, 0.d0, 0.d0/)
!DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) '! Error no T-fit for CiJ for oxygen'
!DEC$ END IF ! WO_WRITE_TO_CONSOLE
    end if
    
end if

end subroutine hdrt_params_oxygen
!******************************************************************************
!******************************************************************************




!******************************************************************************
subroutine hdrt_params_propane(gl,J, hdrt_structure_flag)
!******************************************************************************
!  Definition of the propane-hydrate properties




implicit none

    type(type_gl) :: gl


! Input arguments
integer, intent(in) :: J, hdrt_structure_flag

    !Pure Hydrate? -> structure is known
    !if (N_guests == 1) hdrt_structure = 's2'

    ! Kihara potential parameters
    ! ---------------------------
    gl%a_hc(J) = 0.6502d0          ! [A] hard core radius ... Tee et al. 1966
    !! gw_B0 and hw_B0 CONSTANT for sII from 'z_Kiharas_3_phase_HdrNr_2016.xlsx'
    !    gw_B0 = 904.887d0      ! Gibbs energy of water in empty lattice [J/mol] 
    !    hw_B0 = -5044.369d0    ! Enthalpy of empty lattice [J.mol^-1]
    !! Kihara parameters fitted to CPO and hdrt_composition data with uncer  
    gl%sigma(J) = 3.47023d0     ! for a_ref = 17.31d0
    gl%eps_k(J) = 187.718d0     ! for a_ref = 17.31d0 - taken from correlation of CSMGem occup_L
    
    gl%sigmad(J) = 3.19624d0
    gl%epsd_k(J) = 484.879d0
    
    gl%r_vdw(J) = 2.0778d0

if (hdrt_structure_flag==1) then ! sI hydrate   
    
    ! need to be fitted
    gl%a0_hdrt(J) = 12.d0         ! [A] lattice param @ p0_a = 0 Pa, T0_a = 275.15 K
    gl%B0_hdrt(1,J) = 1.d10       ! [Pa] Bulk modulus

        
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/0.d0, 0.d0, 0.d0/)
        gl%KS_l(:,J) = (/0.d0, 0.d0, 0.d0/)
!DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) '! Warning: propane does not form sI - see T-fit for CiJ'
!DEC$ END IF ! WO_WRITE_TO_CONSOLE
    end if
    
elseif (hdrt_structure_flag==2) then ! sII hydrate

    gl%a0_hdrt(J) = 17.288930d0   ! [A] lattice param @ p0_a = 1 Pa, T0_a = 273.15 K - October 2016
    gl%B0_hdrt(1,J) = 1.d10             ! [Pa] Bulk modulus    
    
    if (gl%Langmuir=='ks') then    ! Klauda & Sandler (2003) temperature fit for CiJ
        gl%KS_s(:,J) = (/0.d0, 0.d0, 0.d0/)
        gl%KS_l(:,J) = (/-23.1307d0, 4176.1979d0, 45939.4593d0/)
    end if
    
end if

    end subroutine hdrt_params_propane
    
    
!!>Refstate for gibbs energy and enthalpy. Currently only implemented for pure hydrate formation.
!subroutine hdrt_params_refstate_gh(J)
!

!
!implicit none
!
!!>Guest Nr
!integer, intent(in) :: J
!
!    if (hdrt_structure == 's1') then
!        gw_B0 = -0.16240D0*a0_hdrt(J)**3+1272.4D0   !991.6437618717270d0     ! Gibbs energy of water in empty lattice [J/mol] 
!        hw_B0 = 5.5880D0*a0_hdrt(J)**3-14476.0D0    !-4815.495944206963d0    ! Enthalpy of empty lattice [J.mol^-1]
!        ! Kihara parameters fitted to CPO and hdrt_composition data with uncer
!        sigma = 2.97175d0
!        eps_k = 176.242d0
!    elseif (hdrt_structure == 's2') then
!        ! gw_B0 and hw_B0 from linear function for sI based on CSMGem data ... 23.8.14
!        gw_B0 = 917.816d0     ! Gibbs energy of water in empty lattice [J/mol] 
!        hw_B0 = -5094.306d0    ! Enthalpy of empty lattice [J.mol^-1]
!        ! Kihara parameters fitted to CPO and hdrt_composition data with uncer
!        sigma = 3.52825d0     ! for a_ref = 17.31d0
!        eps_k = 188.238d0     ! for a_ref = 17.31d0
!    endif
!
!end subroutine hdrt_params_refstate_gh
!!******************************************************************************
!!******************************************************************************


    end module hdrt_params_module
