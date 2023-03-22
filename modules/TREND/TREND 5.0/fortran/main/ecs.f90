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

    ! module for file ecs.f90
    module ecs_module
    !global use inclusion
    use module_all_types
    use calc_functions


    contains

	



    
    
!subroutine initialize_ecs
!

!
!!reference fluid variables
!input_0 = 'td'
!t_0 = 300.d0
!d_0 = 10.d0
!fluids_0 = ""
!moles_0 = "1"
!molev_0 = 0.d0
!eos_indicator_0 = ""
!eos_indicator_0_call = "1"
!tc_0 = 0.d0
!dc_0 = 0.d0
!m_0 = 0.d0
!sigma_0 = 0.d0
!epsilon_0 = 0.d0
!
!!fluid variables
!tc_fluid = 0.d0
!dc_fluid = 0.d0
!w_fluid = 0.d0
!m_fluid = 0.d0
!ncomp_fluid = 0
!alpha_temp = 0.d0
!dalpha_ddelta_temp = 0.d0
!d2alpha_ddelta2_temp = 0.d0
!dalpha_dtau_temp = 0.d0
!d2alpha_ddelta_dtau_temp = 0.d0
!z_temp = 0.d0
!residual_ecs = 0.d0
!    
!end subroutine initialize_ecs
    

    
!Calculation of the ideal visosity of fluid "1" [mycroPas]. 
!See McLinden et al. (gl,1999): "An Extended Corresponding States Model for the Thermal Conductivity of Refrigerants and Refrigerant Mixtures".
double precision function viscosity_ideal(gl,t, comp)



implicit none

    type(type_gl) :: gl




double precision :: t
integer :: comp

double precision :: sigma_ecs    !Lennard-Jones-Parameter. See Huber et al (1992): "Prediction of the Viscosity of Refrigerants and Refrigerant Mixtures"
double precision :: epsilon_ecs    !Lennard-Jones-Parameter. See Huber et al (1992): "Prediction of the Viscosity of Refrigerants and Refrigerant Mixtures"
double precision :: t_red    !reduced temperature used in the "collision integral" of the Lennard-Jones_Fluids = (k * T / epsilon_ecs)
double precision :: collisionintegral    !"collision integral" of the LF fluid. See Neufeld et al. (1972): "Empirical Equations to Calculate 16 of the Transport Collision Integrals..."
double precision :: a_0_ecs
double precision :: b_0_ecs
double precision :: c_0_ecs
double precision :: d_0_ecs 
double precision :: e_0_ecs 
double precision :: f_0_ecs 
double precision :: g_0_ecs 
double precision :: h_0_ecs 
double precision :: p_0_ecs 
double precision :: r_0_ecs 
double precision :: s_0_ecs 
double precision :: v_0_ecs 


!initializing
sigma_ecs = 0.d0  
epsilon_ecs = 0.d0  
t_red = 0.d0   
collisionintegral = 0.d0
a_0_ecs = 1.16145d0
b_0_ecs = 0.14874d0
c_0_ecs = 0.52487d0
d_0_ecs = 0.7732d0
e_0_ecs = 2.16178d0
f_0_ecs = 2.43787d0
g_0_ecs = 0.d0
h_0_ecs = 0.d0
p_0_ecs = 7.27371d0
r_0_ecs = -6.435d-4
s_0_ecs = 18.0323d0
v_0_ecs = -0.7683d0


!calculation
sigma_ecs = gl%sigma_0 * 1.d9 * (gl%dc_0 / gl%dc_fluid(comp)) ** (1.d0 / 3.d0)    !Sigma in [nm]
epsilon_ecs = gl%epsilon_0 * (gl%tc_fluid(comp) / gl%tc_0)
t_red = t / epsilon_ecs
collisionintegral = (a_0_ecs / t_red ** b_0_ecs) + (c_0_ecs / dexp(d_0_ecs * t_red)) + (e_0_ecs / dexp(f_0_ecs * t_red)) + &
                  & (g_0_ecs / dexp(h_0_ecs * t_red))    !without additional sinus term, in order to improve the results of the (2, 2)-CollisionIntegral

viscosity_ideal = 26.69d-3 * ((gl%m_fluid(comp) * 1.d3) * t) ** (0.5d0) / (sigma_ecs ** 2 * collisionintegral)    !M in [g/mol]

!check result
if (viscosity_ideal <= 0) then
    viscosity_ideal = -5301
end if

end function 
    
    
    
!Calculation of the ideal visosity of the fluid mixture [mycroPas]. 
!See McLinden et al. (1999): "An Extended Corresponding States Model for the Thermal Conductivity of Refrigerants and Refrigerant Mixtures" und Hirschfelder et al. (1954): "Molecular Theory of Gases and Liquids" bzw. vergleiche mit REFPROP Quellcode.
double precision function viscosity_ideal_mix(gl,t, x)



implicit none

    type(type_gl) :: gl



!Declaration
!Input
double precision :: t
double precision, dimension(30) :: x
double precision, dimension(30) :: sigma_i
double precision, dimension(30) :: epsilon_i
double precision, dimension(30) :: tred_i
double precision, dimension(30, 30) :: m_ij
double precision, dimension(30, 30) :: sigma_ij
double precision, dimension(30, 30) :: tred_ij
double precision :: collisionintegral11
double precision :: collisionintegral22
double precision :: term
double precision :: rho_ik
double precision :: eta_ik
double precision, dimension(30) :: eta_i
double precision :: d_ij
double precision :: d_jk
double precision :: sum_k
double precision, dimension(31, 31) :: h
double precision :: eta_mix
!constants
double precision :: a_11
double precision :: b_11
double precision :: c_11
double precision :: d_11
double precision :: e_11
double precision :: f_11
double precision :: g_11
double precision :: h_11
double precision :: a_22
double precision :: b_22
double precision :: c_22
double precision :: d_22
double precision :: e_22
double precision :: f_22
double precision :: g_22
double precision :: h_22

integer :: i
integer :: j
integer :: k


!initializing
sigma_i = 0.d0
epsilon_i = 0.d0
tred_i = 0.d0
m_ij = 0.d0
sigma_ij = 0.d0
tred_ij = 0.d0
collisionintegral11 = 0.d0
collisionintegral22 = 0.d0
term = 0.d0
rho_ik = 0.d0
eta_ik = 0.d0
eta_i = 0.d0
d_ij = 0.d0
d_jk = 0.d0
sum_k = 0.d0
h = 1.d0
eta_mix = 0.d0
viscosity_ideal_mix = 0.d0

a_11 = 1.06036d0
b_11 = 0.1561d0
c_11 = 0.193d0
d_11 = 0.47635d0
e_11 = 1.03587d0
f_11 = 1.52996d0
g_11 = 1.76474d0
h_11 = 3.89411d0
a_22 = 1.16145d0
b_22 = 0.14874d0
c_22 = 0.52487d0
d_22 = 0.7732d0
e_22 = 2.16178d0
f_22 = 2.43787d0
g_22 = 0.d0
h_22 = 0.d0

i = 0
j = 0
k = 0


!calculation
!get pure fluid results
do i = 1, gl%ncomp_fluid
    sigma_i(i) = gl%sigma_0 * 1.d9 * (gl%dc_0 / gl%dc_fluid(i)) ** (1.d0 / 3.d0)
    epsilon_i(i) = gl%epsilon_0 * (gl%tc_fluid(i) / gl%tc_0)
    tred_i(i) = t / epsilon_i(i)
end do

!mixture contribution
do i = 1, gl%ncomp_fluid
    do j = 1, gl%ncomp_fluid
        m_ij(i, j) = 2.d0 * gl%m_fluid(i) * gl%m_fluid(j) / (gl%m_fluid(i) + gl%m_fluid(j)) * 1.d3
        sigma_ij(i, j) = 0.5d0 * (sigma_i(i) + sigma_i(j))
        tred_ij(i, j) = (tred_i(i) * tred_i(j)) ** 0.5d0
    end do
end do

!calculate viscosity
do i = 1, gl%ncomp_fluid
    do j = 1, gl%ncomp_fluid
        sum_k = 0.d0
        if (i == j) then
            d_ij = 1.d0
        else
            d_ij = 0.d0
        end if
        do k = 1, gl%ncomp_fluid
            collisionintegral11 = 0.d0
            collisionintegral22 = 0.d0
            term = 0.d0
            rho_ik = 0.d0
            eta_ik = 0.d0
            if (j == k) then
                d_jk = 1.d0
            else
                d_jk = 0.d0
            end if
            collisionintegral11 = (a_11 / tred_ij(i, k) ** b_11) + (c_11 / dexp(d_11 * tred_ij(i, k))) + (e_11 / dexp(f_11 * tred_ij(i, k))) + &
                                & (g_11 / dexp(h_11 * tred_ij(i, k)))
            collisionintegral22 = (a_22 / tred_ij(i, k) ** b_22) + (c_22 / dexp(d_22 * tred_ij(i, k))) + (e_22 / dexp(f_22 * tred_ij(i, k))) + &
                                & (g_22 / dexp(h_22 * tred_ij(i, k)))
            term = (m_ij(i, k) * t) ** 0.5d0 / sigma_ij(i, k) ** 2
            rho_ik = 0.3203d0 * term / collisionintegral11
            eta_ik = 26.69d-3 * term / collisionintegral22
            term = (m_ij(j, j) / (rho_ik * m_ij(k, k))) * (d_ij - d_jk) + 0.5d0 * (d_ij + d_jk) / eta_ik
            sum_k = sum_k + x(k) * m_ij(j, k) ** 2 * term
        end do
        h(i, j) = sum_k / (m_ij(i, i) * m_ij(j, j))
    end do
end do

do i = 1, gl%ncomp_fluid
    do j = i + 1, gl%ncomp_fluid + 1
        h(i, j) = h(i, j) / h(i, i)
    end do
    h(i, i) = 1.d0
    do j = 1, gl%ncomp_fluid
        if (i /= j) then
            do k = i + 1, gl%ncomp_fluid + 1
                h(j, k) = h(j, k) - h(j, i) * h(i, k)
            end do
            h(j, i) = 0.d0
        end if
    end do
end do

do i = 1, gl%ncomp_fluid
    eta_mix = eta_mix + x(i) * h(i, gl%ncomp_fluid + 1)
end do

viscosity_ideal_mix = eta_mix


!check result
if (viscosity_ideal_mix <= 0) then
    viscosity_ideal_mix = -5301
end if
    
end function 

    
!Calculation of viscosity according to the Extended Corresponding States Model. 
!See Klein et al. (1997): "An Improved Extended Corresponding States Method for Estimation of Viscosity of pure Refrigerants and Mixtures".
!double precision function viscosity_ecs(input, prop1, prop2, fluids, fluid_ref, moles, eos_indicator, eos_indicator_ref, path)    ![mycroPa s]
!!DEC$ ATTRIBUTES DLLEXPORT :: viscosity_ecs
!



!implicit none
!
!    type(type_gl) :: gl
!
!
!
!!Input variables
!character(12) :: input
!double precision :: prop1
!double precision :: prop2
!character(255) :: fluids
!character(255) :: fluid_ref
!character(255) :: moles
!character(255) :: eos_indicator
!character(255) :: eos_indicator_ref
!character(255) :: path
!!ancillary variables
!double precision, dimension(30) :: molev
!integer :: errorflag
!double precision :: t,d,p,h,s,dvap,dliq
!double precision :: A00_CALC
!double precision :: A01_CALC
!double precision :: A02_CALC
!double precision :: A10_CALC
!double precision :: A11_CALC
!double precision :: eta_ideal
!double precision :: viscosity_ideal
!double precision :: viscosity_ideal_mix
!double precision :: t0_loc
!double precision :: d0_loc
!double precision :: fx
!double precision :: hx
!double precision :: gx
!double precision :: eta_multiplikator_f
!double precision :: multiplikator_f_viscosity
!double precision :: multiplikator_f_mix
!double precision :: dummy_eta
!double precision :: visdyn_calc
!double precision :: eta_residual
!integer :: phase, nrsubst
!double precision :: vapfrac
!integer :: numberofphases
!double precision, dimension(5) :: rho                
!double precision, dimension(30, 5) :: x_phase(30,5)   
!integer, dimension(5) :: phasetype
!integer :: calling_prop
!double precision, dimension(5) :: phasefrac    ! Phase amounts for all phases
!integer, dimension(5) :: phase_type
!
!
!
! 
!
!!initializing
!!call initialize_ecs
!viscosity_ecs = 0.d0
!molev = 0.d0
!errorflag = 0
!t = 0.d0
!p = 0.d0
!d = 0.d0
!eta_ideal = 0.d0
!t0_loc = 0.d0
!d0_loc = 0.d0
!fx = 0.d0
!hx = 0.d0
!gx = 0.d0
!eta_multiplikator_f = 0.d0
!dummy_eta = 0.d0
!eta_residual = 0.d0
!p = 0.d0
!dvap = 0.d0
!dliq = 0.d0
!phase = 0
!vapfrac = 0.d0
!numberofphases = 0
!rho = 0.d0
!x_phase = 0.d0
!phasetype = 0
!!set module variables needed to call reference fluid
!gl%fluids_0 = fluid_ref
!gl%eos_indicator_0 = eos_indicator_ref
!
!
!
!!start of calculation
!call setup(gl,input, prop1, prop2, fluids, moles, molev, path, eos_indicator, errorflag)
!if (errorflag /= 0) then      
!    viscosity_ecs = errorflag
!    return
!end if 
!
!!set fluid module variables
!gl%tc_fluid = gl%tc
!gl%dc_fluid = gl%rhoc
!gl%w_fluid = gl%accen
!gl%m_fluid = gl%wm
!gl%ncomp_fluid = gl%ncomp
!
!call inpt_handle(gl,input, prop1, prop2, p, t, d, dvap, dliq, molev, rho, x_Phase, phasetype, phasefrac, numberofphases, nrsubst, errorflag)
!
!!    
!!if (input == 'td') then
!!    t = prop1
!!    d = prop2
!!    if (ncomp == 1) then
!!        nrsubst = 1
!!        call PhaseDet_pure(p, t, d, dvap, dliq, phase, vapfrac, numberofphases, errorflag)
!!    else
!!        nrsubst = 0
!!        call PhaseDet_td(p, t, molev, rho, x_Phase, phasetype, vapfrac, d, numberofphases, errorflag)
!!        d = rho(phasetype(1)) ! if only one phase present, it can be found in phasetype(1), otherwise numberofphases > 1 -- >  exit function
!!    end if
!! 
!!else if (input == 'tp') then
!!    t = prop1
!!    p = prop2
!!    if (ncomp == 1) then
!!        nrsubst = 1
!!        call PhaseDet_pure_tp(p, t, d, dvap, dliq, phase, vapfrac, numberofphases, errorflag)
!!    else
!!        nrsubst = 0
!!        call PhaseDet(p, t, molev, rho, x_Phase, phasetype, vapfrac, numberofphases, errorflag)
!!        d = rho(phasetype(1)) ! if only one phase present, it can be found in phasetype(1), otherwise numberofphases > 1 -- >  exit function
!!    end if 
!!
!!else if (input == 'ph') then
!!    p = prop1
!!    h = prop2
!!    if (ncomp == 1) then
!!        nrsubst = 1
!!        call PhaseDet_ph_pure(p, t, d, dvap, dliq, phase, vapfrac, h, numberofphases, errorflag)
!!    else
!!        nrsubst = 0
!!        call PhaseDet_ph(p, t, molev, rho, x_phase, phasetype, vapfrac, h, numberofphases, errorflag)
!!        d = rho(phasetype(1)) ! if only one phase present, it can be found in phasetype(1), otherwise numberofphases > 1 -- >  exit function
!!    end if
!!
!!else if (input == 'ps') then
!!    p = prop1
!!    s = prop2
!!    if (ncomp == 1) then
!!        nrsubst = 1
!!        call PhaseDet_ps_pure(p, t, d, dvap, dliq, phase, vapfrac, s, numberofphases, errorflag)
!!    else
!!        nrsubst = 0
!!        call PhaseDet_ps(p, t, molev, rho, x_phase, phasetype, vapfrac, s, numberofphases, errorflag)
!!        d = rho(phasetype(1)) ! if only one phase present, it can be found in phasetype(1), otherwise numberofphases > 1 -- >  exit function
!!    end if
!!
!!else
!!    errorflag = -290885
!!    return
!!end if
!!    
!    
!
!    
!if (errorflag /= 0) then 
!    viscosity_ecs = errorflag
!    return
!elseif (numberofphases > 1) then
!    viscosity_ecs = -6666.d0
!    return
!end if
!
!
!
!!get derivatives for Newton method
!gl%alpha_temp = A00_CALC(gl,t, d, nrsubst)
!gl%dalpha_ddelta_temp = A01_CALC(gl,t, d, nrsubst)
!gl%d2alpha_ddelta2_temp = A02_CALC(gl,t, d, nrsubst)
!gl%dalpha_dtau_temp = A10_CALC(gl,t, d, nrsubst)
!gl%d2alpha_ddelta_dtau_temp = A11_CALC(gl,t, d, nrsubst)
!gl%z_temp = 1.d0 + gl%dalpha_ddelta_temp
!
!! save phase of the given fluid
!phase_type = gl%phase_id
!
!!set module variables of reference fluid
!call setup(gl,gl%input_0, gl%t_0, gl%d_0, gl%fluids_0, gl%moles_0, gl%molev_0, path, gl%eos_indicator_0, errorflag)
!if (errorflag /= 0) then      
!    viscosity_ecs = errorflag
!    return
!end if
!gl%tc_0 = gl%tc(1)
!gl%dc_0 = gl%rhoc(1)
!gl%m_0 = gl%wm(1)
!if (trim(gl%fluids_0) == "propane") then
!    gl%sigma_0 = 0.47d-9             
!    gl%epsilon_0 = 358.9d0
!elseif (trim(gl%fluids_0) == "r134a") then
!    gl%sigma_0 = 0.468932d-9           
!    gl%epsilon_0 = 299.363d0
!elseif (trim(gl%fluids_0) == "nitrogen") then
!    gl%sigma_0 = 0.3656d-9           
!    gl%epsilon_0 = 98.94d0
!elseif (trim(gl%fluids_0) == "octane") then
!    gl%sigma_0 = 0.63617d-9           
!    gl%epsilon_0 = 452.090d0
!elseif (trim(gl%fluids_0) == "hydrogen") then
!    gl%sigma_0 = 0.2827d-9           
!    gl%epsilon_0 = 59.7d0
!else !invalid reference fluid
!    viscosity_ecs = -5300
!    return
!end if
!    
!
!    
!if (gl%ncomp_fluid == 1) then        
!        
!    !viscosity
!    !ideal part
!    eta_ideal = viscosity_ideal(gl,t, 1)
!    if (eta_ideal <= 0.d0) then
!        viscosity_ecs = -5301
!        return
!    end if
!    
!    !find t and d reference parameters
!    call reference_parameters(gl,t, d, t0_loc, d0_loc, fx, hx, errorflag)
!    if (errorflag /= 0) then
!        viscosity_ecs = errorflag
!        return
!    end if
!    
!    !multiplicator
!    eta_multiplikator_f = multiplikator_f_viscosity(gl,fx, hx)
!    if (eta_multiplikator_f <= 0.d0) then
!        viscosity_ecs = -5301
!        return
!    end if
!
!
!!mixtures       
!else   
!    !viscosity
!    !ideal part
!    eta_ideal = viscosity_ideal_mix(gl,t, molev)
!    if (eta_ideal <= 0.d0) then
!        viscosity_ecs = -5301
!        return
!    end if
!    
!    !find t and d reference parameters
!    calling_prop = 1
!    call reference_parameters_mix(gl,t, d, molev, t0_loc, d0_loc, fx, hx, gx, calling_prop, errorflag)
!    if (errorflag /= 0) then
!        viscosity_ecs = errorflag
!        return
!    end if
!    
!    !Multiplicator
!    eta_multiplikator_f = multiplikator_f_mix(gl,fx, hx, gx)
!    if (eta_multiplikator_f <= 0.d0) then
!        viscosity_ecs = -5301
!        return
!    end if
!end if
!
!
!!Residual part (using viscosity model of the reference fluid)
!if (trim(gl%eos_indicator_0) /= "1") then    !the standard fluid file needs to be read in
!    call setup(gl,gl%input_0, gl%t_0, gl%d_0, gl%fluids_0, gl%moles_0, gl%molev_0, path, gl%eos_indicator_0_call, errorflag)
!    if (errorflag /= 0) then      
!        viscosity_ecs = errorflag
!        return
!    end if
!end if
!dummy_eta = visdyn_calc(gl,t0_loc, d0_loc, 1)
!eta_residual = gl%residual_ecs
!gl%residual_ecs = 0.d0
!    
!!summation of individual parts
!viscosity_ecs = eta_ideal + eta_multiplikator_f * eta_residual
!if (viscosity_ecs <= 0.d0) then
!    viscosity_ecs = -5301
!end if
!
!! save phase of the given fluid
!gl%phase_id = phase_type
!
!end function
!    
    
    
!Calculation of viscosity according to the Extended Corresponding States Model. 
!See McLinden et al. (1999): "An Extended Corresponding States Model for the Thermal Conductivity of Refrigerants and Refrigerant Mixtures".
!double precision function thermal_conductivity_ecs(input, prop1, prop2, fluids, fluid_ref, moles, eos_indicator, eos_indicator_ref, path)    ![W / (m K)]
!!DEC$ ATTRIBUTES DLLEXPORT :: thermal_conductivity_ecs
!



!implicit none
!
!    type(type_gl) :: gl
!
!
!
!!Input variables
!character(12) :: input
!double precision :: prop1
!double precision :: prop2
!character(255) :: fluids
!character(255) :: fluid_ref
!character(255) :: moles
!character(255) :: eos_indicator
!character(255) :: eos_indicator_ref
!character(255) :: path
!!ancillary variables
!double precision, dimension(30) :: molev
!integer :: errorflag
!double precision :: t,d,p,h,s,dvap,dliq
!double precision :: A00_CALC
!double precision :: A01_CALC
!double precision :: A02_CALC
!double precision :: A10_CALC
!double precision :: A11_CALC
!double precision :: lambda_int
!double precision :: lambda_ideal
!double precision :: thermal_conductivity_int
!double precision :: thermal_conductivity_ideal
!double precision :: thermal_conductivity_ideal_mix
!double precision :: t0_loc
!double precision :: d0_loc
!double precision :: fx
!double precision :: hx
!double precision :: gx
!double precision :: lambda_multiplikator_f
!double precision :: multiplikator_f_thermal_conductivity
!double precision :: multiplikator_f_mix
!double precision :: dlamtc1
!double precision :: lambda_residual
!integer :: phase
!double precision :: vapfrac
!integer :: numberofphases
!double precision, dimension(5) :: rho                
!double precision, dimension(30, 5) :: x_phase(30,5)   
!integer, dimension(5) :: phasetype
!integer :: calling_prop, nrsubst
!
! 
!
!!initializing
!!call initialize_ecs
!thermal_conductivity_ecs = 0.d0
!molev = 0.d0
!errorflag = 0
!t = 0.d0
!p = 0.d0
!d = 0.d0
!lambda_int = 0.d0
!lambda_ideal = 0.d0
!t0_loc = 0.d0
!d0_loc = 0.d0
!fx = 0.d0
!hx = 0.d0
!gx = 0.d0
!lambda_multiplikator_f = 0.d0
!lambda_residual = 0.d0
!p = 0.d0
!dvap = 0.d0
!dliq = 0.d0
!phase = 0
!vapfrac = 0.d0
!numberofphases = 0
!rho = 0.d0
!x_phase = 0.d0
!phasetype = 0
!!set module variables needed to call reference fluid
!gl%fluids_0 = fluid_ref
!gl%eos_indicator_0 = eos_indicator_ref
!
!
!
!!set module variables of reference fluid
!call setup(gl,gl%input_0, gl%t_0, gl%d_0, gl%fluids_0, gl%moles_0, gl%molev_0, path, gl%eos_indicator_0, errorflag)
!if (errorflag /= 0) then      
!    thermal_conductivity_ecs = errorflag
!    return
!end if
!    
!gl%tc_0 = gl%tc(1)
!gl%dc_0 = gl%rhoc(1)
!gl%m_0 = gl%wm(1)
!
!if (trim(gl%fluids_0) == "propane") then
!    gl%sigma_0 = 0.47d-9             
!    gl%epsilon_0 = 358.9d0
!elseif (trim(gl%fluids_0) == "r134a") then
!    gl%sigma_0 = 0.468932d-9           
!    gl%epsilon_0 = 299.363d0
!elseif (trim(gl%fluids_0) == "nitrogen") then
!    gl%sigma_0 = 0.3656d-9           
!    gl%epsilon_0 = 98.94d0
!elseif (trim(gl%fluids_0) == "octane") then
!    gl%sigma_0 = 0.63617d-9           
!    gl%epsilon_0 = 452.090d0
!elseif (trim(gl%fluids_0) == "hydrogen") then
!    gl%sigma_0 = 0.2827d-9           
!    gl%epsilon_0 = 59.7d0
!else !invalid reference fluid
!    thermal_conductivity_ecs = -5300
!    return
!end if
!
!
!    
!
!!start of calculation
!!read in given fluid or fluid mixture
!call setup(gl,input, prop1, prop2, fluids, moles, molev, path, eos_indicator, errorflag)
!if (errorflag /= 0) then      
!    thermal_conductivity_ecs = errorflag
!    return
!end if 
!
!!set fluid module variables
!gl%tc_fluid = gl%tc
!gl%dc_fluid = gl%rhoc
!gl%w_fluid = gl%accen
!gl%m_fluid = gl%wm
!gl%ncomp_fluid = gl%ncomp
!
!
!if (input == 'td') then
!    t = prop1
!    d = prop2
!    if (gl%ncomp == 1) then
!        nrsubst = 1
!        call PhaseDet_pure(gl,p, t, d, dvap, dliq, phase, vapfrac, numberofphases, errorflag)
!    else
!        nrsubst = 0
!        call PhaseDet_td(gl,p, t, molev, rho, x_Phase, phasetype, vapfrac, d, numberofphases, errorflag)
!        d = rho(phasetype(1)) ! if only one phase present, it can be found in phasetype(1), otherwise numberofphases > 1 -- >  exit function
!    end if
! 
!else if (input == 'tp') then
!    t = prop1
!    p = prop2
!    if (gl%ncomp == 1) then
!        nrsubst = 1
!        call PhaseDet_pure_tp(gl,p, t, d, dvap, dliq, phase, vapfrac, numberofphases, errorflag)
!    else
!        nrsubst = 0
!        call PhaseDet(gl,p, t, molev, rho, x_Phase, phasetype, vapfrac, numberofphases, errorflag)
!        d = rho(phasetype(1)) ! if only one phase present, it can be found in phasetype(1), otherwise numberofphases > 1 -- >  exit function
!    end if 
!
!else if (input == 'ph') then
!    p = prop1
!    h = prop2
!    if (gl%ncomp == 1) then
!        nrsubst = 1
!        call PhaseDet_ph_pure(gl,p, t, d, dvap, dliq, phase, vapfrac, h, numberofphases, errorflag)
!    else
!        nrsubst = 0
!        call PhaseDet_ph(gl,p, t, molev, rho, x_phase, phasetype, vapfrac, h, numberofphases, errorflag)
!        d = rho(phasetype(1)) ! if only one phase present, it can be found in phasetype(1), otherwise numberofphases > 1 -- >  exit function
!    end if
!
!else if (input == 'ps') then
!    p = prop1
!    s = prop2
!    if (gl%ncomp == 1) then
!        nrsubst = 1
!        call PhaseDet_ps_pure(gl,p, t, d, dvap, dliq, phase, vapfrac, s, numberofphases, errorflag)
!    else
!        nrsubst = 0
!        call PhaseDet_ps(gl,p, t, molev, rho, x_phase, phasetype, vapfrac, s, numberofphases, errorflag)
!        d = rho(phasetype(1)) ! if only one phase present, it can be found in phasetype(1), otherwise numberofphases > 1 -- >  exit function
!    end if
!
!else
!    errorflag = -290885
!    return
!end if
!
!
!    
!if (errorflag /= 0) then 
!    thermal_conductivity_ecs = errorflag
!    return
!elseif (numberofphases > 1) then
!    thermal_conductivity_ecs = -6666.d0
!    return
!end if
!        
!        
!!get derivatives for Newton method
!gl%alpha_temp = A00_CALC(gl,t, d, nrsubst)
!gl%dalpha_ddelta_temp = A01_CALC(gl,t, d, nrsubst)
!gl%d2alpha_ddelta2_temp = A02_CALC(gl,t, d, nrsubst)
!gl%dalpha_dtau_temp = A10_CALC(gl,t, d, nrsubst)
!gl%d2alpha_ddelta_dtau_temp = A11_CALC(gl,t, d, nrsubst)
!gl%z_temp = 1.d0 + gl%dalpha_ddelta_temp
!
!
!                 
!if (gl%ncomp_fluid == 1) then  
!    
!    !reference fluid must have been set in advance (tc_0,dc_0)
!    !actual fluid needs to be current in order to calculate cp0 correctly
!    lambda_int = thermal_conductivity_int(gl,t, 1)  
!    if (lambda_int <= 0.d0) then
!        thermal_conductivity_ecs = -5301
!        return
!    end if    
!            
!    
!    !set module variables of reference fluid
!    call setup(gl,gl%input_0, gl%t_0, gl%d_0, gl%fluids_0, gl%moles_0, gl%molev_0, path, gl%eos_indicator_0, errorflag)
!    if (errorflag /= 0) then      
!        thermal_conductivity_ecs = errorflag
!        return
!    end if
!
!
!    !ideal part
!    lambda_ideal = thermal_conductivity_ideal(gl,t, 1)
!        
!    !find reference parameters t and d
!    call reference_parameters(gl,t, d, t0_loc, d0_loc, fx, hx, errorflag)
!    if (errorflag /= 0) then
!        thermal_conductivity_ecs = errorflag
!        return
!    end if
!    
!    !Multiplicator
!    lambda_multiplikator_f = multiplikator_f_thermal_conductivity(gl,fx, hx)
!    if (lambda_multiplikator_f <= 0.d0) then
!        thermal_conductivity_ecs = -5301
!        return
!    end if
!    
!    !Residual part (using viscosity model of the reference fluid)
!    if (trim(gl%eos_indicator_0) /= "1") then    !the standard fluid file needs to be read in
!        call setup(gl,gl%input_0, gl%t_0, gl%d_0, gl%fluids_0, gl%moles_0, gl%molev_0, path, gl%eos_indicator_0_call, errorflag)
!        if (errorflag /= 0) then      
!            thermal_conductivity_ecs = errorflag
!            return
!        end if
!    end if
!    lambda_residual = dlamtc1(gl,t0_loc, d0_loc, 1)    !not yet implemented in TREND
!     
!    !summation of individual parts
!    thermal_conductivity_ecs = lambda_int + lambda_ideal + lambda_multiplikator_f * lambda_residual
!    if (thermal_conductivity_ecs <= 0.d0) then
!        thermal_conductivity_ecs = -5301
!        return
!    end if
!    
!    
!    
!
!!mixtures       
!else   
!    !ideal part
!    lambda_ideal = thermal_conductivity_ideal_mix(gl,t, molev)
!    if (lambda_ideal <= 0.d0) then
!        thermal_conductivity_ecs = -5301
!        return
!    end if
!    
!    !set module variables of reference fluid
!    call setup(gl,gl%input_0, gl%t_0, gl%d_0, gl%fluids_0, gl%moles_0, gl%molev_0, path, gl%eos_indicator_0, errorflag)
!    if (errorflag /= 0) then      
!        thermal_conductivity_ecs = errorflag
!        return
!    end if
!    
!    !find reference parameters t and d
!    calling_prop = 2
!    call reference_parameters_mix(gl,t, d, molev, t0_loc, d0_loc, fx, hx, gx, calling_prop, errorflag)
!    if (errorflag /= 0) then
!        thermal_conductivity_ecs = errorflag
!        return
!    end if
!    
!    !Multiplicator
!    lambda_multiplikator_f = multiplikator_f_mix(gl,fx, hx, gx)
!    if (lambda_multiplikator_f <= 0.d0) then
!        thermal_conductivity_ecs = -5301
!        return
!    end if
!    
!    !Residual part (using viscosity model of the reference fluid)
!    if (trim(gl%eos_indicator_0) /= "1") then    !the standard fluid file needs to be read in
!        call setup(gl,gl%input_0, gl%t_0, gl%d_0, gl%fluids_0, gl%moles_0, gl%molev_0, path, gl%eos_indicator_0_call, errorflag)
!        if (errorflag /= 0) then      
!            thermal_conductivity_ecs = errorflag
!            return
!        end if
!    end if
!    lambda_residual = dlamtc1(gl,t0_loc, d0_loc, 1)    !not yet implemented in TREND
!     
!    !summation of individual parts
!    thermal_conductivity_ecs = lambda_ideal + lambda_multiplikator_f * lambda_residual
!    if (thermal_conductivity_ecs <= 0.d0) then
!        thermal_conductivity_ecs = -5301
!    end if
!end if
!
!end function
    


!calculation of the contribution of internal molecule movements to the thermal conductivity
!See McLinden et al. (1999): "An Extended Corresponding States Model for the Thermal Conductivity of Refrigerants and Refrigerant Mixtures".
double precision function thermal_conductivity_int(gl,t, comp)






implicit none

    type(type_gl) :: gl



!Declaration
double precision :: t
integer :: comp

double precision :: f_int
double precision :: eta_ideal
double precision :: cp0
!double precision :: CP0_CALC


!Initializing
f_int = 1.32d-3
eta_ideal = 0.d0
cp0 = 0.d0


!calculation
eta_ideal = viscosity_ideal(gl,t, comp)
if (eta_ideal <= 0.d0) then
    thermal_conductivity_int = -5301
    return
end if

cp0 = CP0_CALC(gl,t, comp)
if (cp0 <= 0.d0) then
    thermal_conductivity_int = -5301
    return
end if

thermal_conductivity_int = f_int * eta_ideal / (gl%m_fluid(comp) * 1.d3) * (cp0 - 5.d0 / 2.d0 * 8.314472d0)


!check result
if (thermal_conductivity_int < 0.d0) then    !=0 possible for argon, since cp0 of argon as an "ideal rotund gas" always = 5 / 2 * R
    thermal_conductivity_int = -5301
end if

end function
    
    
    
!Calculation of the ideal thermal conductivity.
!See: McLinden et al. (1999): "An Extended Corresponding States Model for the Thermal Conductivity of Refrigerants and Refrigerant Mixtures".
double precision function thermal_conductivity_ideal(gl,t, comp)




implicit none

    type(type_gl) :: gl



!Declaration
double precision :: t
integer :: comp

double precision :: eta_ideal


!Initializing
eta_ideal = 0.d0


!Calculation
eta_ideal = viscosity_ideal(gl,t, comp)
if (eta_ideal <= 0.d0) then
    thermal_conductivity_ideal = -5301
    return
end if

thermal_conductivity_ideal = 15.d0 * 8.314472d0 * eta_ideal * 1.d-6 / (4.d0 * gl%m_fluid(comp))


!check results
if (thermal_conductivity_ideal <= 0.d0) then
    thermal_conductivity_ideal = -5301
end if

end function
    
    
    
!Calculation of mixture ideal thermal conductivity - combining int and ideal part
!See: McLinden et al. (1999): "An Extended Corresponding States Model for the Thermal Conductivity of Refrigerants and Refrigerant Mixtures".
double precision function thermal_conductivity_ideal_mix(gl,t, x)



implicit none

    type(type_gl) :: gl



!Declaration
!Input variables
double precision :: t
double precision, dimension(30) :: x
!Variables
double precision, dimension(30) :: eta_ideal
double precision, dimension(30) :: lambda_int
double precision, dimension(30) :: lambda_ideal
double precision, dimension(30, 30) :: phi_ji
double precision :: temp_1
double precision :: temp_2

integer :: i
integer :: j


!Initializing
eta_ideal = 0.d0
lambda_int = 0.d0
lambda_ideal = 0.d0
phi_ji = 0.d0
temp_1 = 0.d0
temp_2 = 0.d0

i = 0
j = 0


!calculation
!calculate eta_ideal for all pure fluids in mixture
do i = 1, gl%ncomp_fluid
    eta_ideal(i) = viscosity_ideal(gl,t, i)
end do
!calculate lambda_int for all pure fluids in mixture
do i = 1, gl%ncomp_fluid
    lambda_int(i) = thermal_conductivity_int(gl,t, i)
end do
!calculate lambda_ideal for all pure fluids in mixture
do i = 1, gl%ncomp_fluid
    lambda_ideal(i) = thermal_conductivity_ideal(gl,t, i)
end do
do j = 1, gl%ncomp_fluid
    do i = 1, gl%ncomp_fluid
        phi_ji = (1.d0 + (eta_ideal(j) / eta_ideal(i)) ** 0.5d0 * (gl%m_fluid(j) / gl%m_fluid(i)) ** 0.25d0) ** 2 / (8.d0 * (1.d0 + gl%m_fluid(j) / gl%m_fluid(i))) ** 0.5d0
        temp_1 = temp_1 + x(i) * phi_ji(j, i)
    end do
    temp_2 = temp_2 + x(j) * (lambda_int(j) + lambda_ideal(j)) / (temp_1)
    temp_1 = 0.d0
end do
thermal_conductivity_ideal_mix = temp_2


!check result
if (thermal_conductivity_ideal_mix <= 0.d0) then
    thermal_conductivity_ideal_mix = -5301
end if

end function
    
    

!Calculation of reference parameters t0 and d0 as well as the combined parameters f and h for pure fluids.
!See: Klein et al. (1997): "An Improved Extended Corresponding States Method for Estimation of Viscosity of pure Refrigerants and Mixtures" und McLinden et al. (1999): "An Extended Corresponding States Model for the Thermal Conductivity of Refrigerants and Refrigerant Mixtures".
subroutine reference_parameters(gl,t, d, t0_loc, d0, fx, hx, errorflag)





implicit none

    type(type_gl) :: gl



!Declaration
double precision :: t
double precision :: d
double precision :: t0_loc
double precision :: d0
double precision :: fx
double precision :: hx
integer :: errorflag

double precision :: t0_est
double precision :: d0_est
double precision :: d0_est_i
double precision :: f_est
double precision :: h_est
double precision :: theta
double precision :: phi

double precision :: multiplicator
integer :: i


!Initializing
t0_loc = 0.d0
d0 = 0.d0
fx = 0.d0
hx = 0.d0
errorflag = 0

t0_est = 0.d0
d0_est = 0.d0
d0_est_i = 0.d0
f_est = 0.d0
h_est = 0.d0
theta = 0.d0
phi = 0.d0 

multiplicator = 0.d0
i = 0


!Calculation
!estimate start values
theta = theta_calc(gl,t, gl%tc_fluid(1), d, gl%dc_fluid(1), gl%w_fluid(1))
phi = phi_calc(gl,t, gl%tc_fluid(1), d, gl%dc_fluid(1), gl%w_fluid(1))
f_est = gl%tc_fluid(1) / gl%tc_0 * theta
h_est = gl%dc_0 / gl%dc_fluid(1) * phi
t0_est = t / f_est
d0_est = d * h_est

!get true values of t0 and d0
call Get_T0Rho0(gl,t, d, t0_est, d0_est, t0_loc, d0, fx, hx, errorflag)

!check whether calculation was successful
if (errorflag == 0) then
    return

!try different start values for d0 in Newton method
else
    multiplicator = 0.1d0
    do i = 1, 19
        d0_est_i = d0_est * multiplicator
        multiplicator = multiplicator + 0.1d0
        if (i /= 10) then
            errorflag = 0
            call Get_T0Rho0(gl,t, d, t0_est, d0_est_i, t0_loc, d0, fx, hx, errorflag)
            if (errorflag == 0) then
                return
            end if
        end if
    end do
end if
    
!in case of failure
errorflag = -5302

end subroutine
    
    
    
    
!Calculation of reference parameters t0 and d0 as well as the combined parameters f and h for mixtures.
!See: McLinden et al. (1999): "An Extended Corresponding States Model for the Thermal Conductivity of Refrigerants and Refrigerant Mixtures".
!See: Klein et al. (1997): "An Improved Extended Corresponding States Method for Estimation of Viscosity of pure Refrigerants and Mixtures".
subroutine reference_parameters_mix(gl,t, d, x, t0_loc, d0, fx, hx, gx, calling_prop, errorflag)





implicit none

    type(type_gl) :: gl



!Declaration
double precision :: t
double precision :: d
double precision, dimension(30) :: x
double precision :: t0_loc
double precision :: d0
double precision :: fx
double precision :: hx
double precision :: gx
integer :: errorflag

double precision :: t0_est
double precision :: d0_est
double precision :: d0_est_i
double precision, dimension(30) :: f_i_est
double precision, dimension(30) :: h_i_est
double precision, dimension(30, 30) :: f_ij_est
double precision, dimension(30, 30) :: h_ij_est
double precision :: f_x_h_x_est
double precision :: f_x_est
double precision :: h_x_est
double precision :: theta
double precision :: phi
double precision, dimension(30, 30) :: m_ij
double precision :: temp

double precision :: multiplicator
integer :: i
integer :: j
integer :: calling_prop ! = 1 if viscosity, = 2 if thermal conductivity



!Initializing
t0_loc = 0.d0
d0 = 0.d0
fx = 0.d0
hx = 0.d0
gx = 0.d0
errorflag = 0

t0_est = 0.d0
d0_est = 0.d0
d0_est_i = 0.d0
f_i_est = 0.d0
h_i_est = 0.d0
f_ij_est = 0.d0
h_ij_est = 0.d0
f_x_h_x_est = 0.d0
f_x_est = 0.d0
h_x_est = 0.d0
theta = 0.d0
phi = 0.d0
m_ij = 0.d0
temp = 0.d0

multiplicator = 0.d0
i = 0
j = 0



!Calculation
!estimate start values
!get t0_est and d0_est for each component
do i = 1, gl%ncomp_fluid
    theta = theta_calc(gl,t, gl%tc_fluid(i), d, gl%dc_fluid(i), gl%w_fluid(i))
    phi = phi_calc(gl,t, gl%tc_fluid(i), d, gl%dc_fluid(i), gl%w_fluid(i))
    f_i_est(i) = gl%tc_fluid(i) / gl%tc_0 * theta
    h_i_est(i) = gl%dc_0 / gl%dc_fluid(i) * phi
end do
do i = 1, gl%ncomp_fluid
    do j = 1, gl%ncomp_fluid
        f_ij_est(i, j) = (f_i_est(i) * f_i_est(j)) ** 0.5d0
        h_ij_est(i, j) = ((h_i_est(i) ** (1.d0 / 3.d0) + h_i_est(j) ** (1.d0 / 3.d0)) ** 3) / 8.d0
        h_x_est = h_x_est + x(i) * x(j) * h_ij_est(i, j)
        f_x_h_x_est = f_x_h_x_est + x(i) * x(j) *  f_ij_est(i, j) * h_ij_est(i, j)
    end do
end do
f_x_est = f_x_h_x_est / h_x_est
t0_est = t / f_x_est
d0_est = d * h_x_est

!get true values of t0 and d0
call Get_T0Rho0(gl,t, d, t0_est, d0_est, t0_loc, d0, fx, hx, errorflag)

!check whether calculation was successful
if (errorflag /= 0) then
    multiplicator = 0.1d0
    do i = 1, 19
        d0_est_i = d0_est * multiplicator
        multiplicator = multiplicator + 0.1d0
        if (i /= 10) then
            errorflag = 0
            call Get_T0Rho0(gl,t, d, t0_est, d0_est_i, t0_loc, d0, fx, hx, errorflag)
            if (errorflag == 0) then
                exit
            end if
        end if
    end do
end if
  
!if successful get gx
if (errorflag == 0) then
    do i = 1, gl%ncomp_fluid
        do j = 1, gl%ncomp_fluid
            m_ij(i, j) = 2.d0 / ((1.d0 / (gl%m_fluid(i) * 1000.d0)) + (1.d0 / (gl%m_fluid(j) * 1000.d0)))
            if (calling_prop == 1) then
                 temp = temp + x(i) * x(j) * f_ij_est(i, j) ** (0.5d0) * h_ij_est(i, j) ** (4.d0 / 3.d0) * m_ij(i, j) ** (0.5d0)    !REFPROP calculates f_ij and h_ij exactly. Here estimated values are used.
           else if (calling_prop == 2) then
                temp = temp + x(i) * x(j) * f_ij_est(i, j) ** (0.5d0) * h_ij_est(i, j) ** (4.d0 / 3.d0) * m_ij(i, j) ** (-0.5d0)    !REFPROP calculates f_ij and h_ij exactly. Here estimated values are used.
            end if
        end do
    end do
    
    if (calling_prop == 1) then
        gx = (temp / (fx ** (0.5d0) * hx ** (4.d0 / 3.d0))) ** 2 / (gl%m_0 * 1000.d0)    !Factor "1.d0 / (m_0 * 1000.d0)" is missing in Klein-Paper
    else if (calling_prop == 2) then
        gx = (temp / (fx ** (0.5d0) * hx ** (4.d0 / 3.d0))) ** 2 * (gl%m_0 * 1000.d0)    
    else
        errorflag = -5301
    end if
else
    errorflag = -5302
    return
end if
    
end subroutine
    
       

    
!Calculation of viscosity residual part multiplicator
double precision function multiplikator_f_viscosity(gl,fx, hx)



implicit none

    type(type_gl) :: gl



!Declaration
double precision :: fx
double precision :: hx


multiplikator_f_viscosity = fx ** (0.5d0) * hx ** (-2.d0 / 3.d0) * (gl%m_fluid(1) / gl%m_0) ** (0.5d0)

!check result
if (multiplikator_f_viscosity <= 0) then
    multiplikator_f_viscosity = -5301
end if

end function
    
    
    
!Calculation of thermal conductivity residual part multiplicator
double precision function multiplikator_f_thermal_conductivity(gl,fx, hx)



implicit none

    type(type_gl) :: gl



!Declaration
double precision :: fx
double precision :: hx


multiplikator_f_thermal_conductivity = fx ** (0.5d0) * hx ** (-2.d0 / 3.d0) * (gl%m_0 / gl%m_fluid(1)) ** (0.5d0)

!check result
if (multiplikator_f_thermal_conductivity <= 0) then
    multiplikator_f_thermal_conductivity = -5301
end if

end function
    
 
    
!Calculation of thermal conductivity residual part multiplicator for mixtures
double precision function multiplikator_f_mix(gl,fx, hx, gx)



implicit none

    type(type_gl) :: gl



!Declaration
double precision :: fx
double precision :: hx
double precision :: gx

multiplikator_f_mix =0.d0

multiplikator_f_mix = fx ** 0.5d0 * hx ** (-2.d0 / 3.d0) * gx ** 0.5d0     

!check result
if (multiplikator_f_mix <= 0) then
    multiplikator_f_mix = -1
end if

end function
    
    
    
!Determination of exact reference values t0 and d0
subroutine Get_T0Rho0(gl,t, d, t0_est, d0_est, t0_loc, d0, fx, hx, errorflag)


implicit none

    type(type_gl) :: gl



!Declaration
double precision :: t  
double precision :: d   
double precision :: t0_est   
double precision :: d0_est    
double precision :: t0_loc    
double precision :: d0   
double precision :: fx
double precision :: hx
integer :: errorflag

double precision, dimension(101) :: x
double precision, dimension(101) :: y
double precision, dimension(2) :: f
double precision :: residuum
double precision, dimension(2) :: f_final
double precision, dimension(2, 2) :: A
double precision, dimension(2, 2) :: A_final
double precision :: dx
double precision :: dy

integer :: i
integer :: j


!Initializing
t0_loc = 0.d0
d0 = 0.d0   
fx = 0.d0
hx = 0.d0
errorflag = 0

x = 0.d0
y = 0.d0
i = 0
j = 0


!set first estimation
x(1) = t0_est    !x: vector containing T0 values
y(1) = d0_est    !y: vector containing Rho0 values

!Newton method
do i = 1, 100
    !Initializing
    f = 0.d0
    residuum = 0.d0
    f_final = 0.d0
    A = 0.d0
    A_final = 0.d0
    dx = 0.d0
    dy = 0.d0
    
    !Calculate vector f
    call Get_f (gl,x(i), y(i), f)
    
    !check whether iteration has already been succesfull
    residuum = 0.5d0 * (f(1) ** 2 + f(2) ** 2)
    if (Residuum < 1.d-8) then
        t0_loc = x(i)
        d0 = y(i)
        fx = t / t0_loc
        hx = d0 / d
        if ((t0_loc <= 0.d0) .or. (d0 <= 0.d0)) then
            errorflag = -1
            return
        end if
        return
    end if
    
    !Calculation of final vector f_final
    f_final = -f

    !Calculation of final matrix A_final
    call Get_A (gl,x(i), y(i), A)
    call Get_Inverse(gl,A, A_final)

    !Calculation of dx and dy 
    dx = A_final(1, 1) * f_final(1) + A_final(1, 2) * f_final(2)
    dy = A_final(2, 1) * f_final(1) + A_final(2, 2) * f_final(2)
    !write (*,*) dx, dy
    
    !Calculation der neuen Werte für x und y
    !Calculate new x and y
    x(i + 1) = x(i) + dx
    y(i + 1) = y(i) + dy
end do

errorflag = -5302

end subroutine
    
    
    
!Calculation of vector f
subroutine Get_f(gl,t0_loc, d0, f)





implicit none

    type(type_gl) :: gl



!Declaration
double precision :: t0_loc
double precision :: d0
double precision, dimension(2) :: f

!double precision :: A00_CALC
!double precision :: A01_CALC
double precision :: a0
double precision :: dalpha_ddelta_0
double precision :: z0



!Initializing
f = 0.d0

a0 = 0.d0
dalpha_ddelta_0 = 0.d0
z0 = 0.d0



!Calculation of a0 and z0 of the reference fluid at the equivalent temperature and pressure
a0 = A00_CALC(gl,t0_loc, d0, 0)
dalpha_ddelta_0 = A01_CALC(gl,t0_loc, d0, 0)
z0 = 1.d0 + dalpha_ddelta_0

!Calculation of f1 and f2
f(1) = a0 - gl%alpha_temp
f(2) = z0 - gl%z_temp


end subroutine
    
    
    
!Calculation of matrix A containing analytical derivatives
subroutine Get_A(gl,t0_loc, d0, A)




implicit none

    type(type_gl) :: gl



!Declaration
double precision :: t0_loc
double precision :: d0
double precision, dimension(2, 2) :: A

!double precision :: A01_CALC    !dAlpha_dDelta * Delta
!double precision :: A02_CALC    !d2Alpha_dDelta2 * Delta2
!double precision :: A10_CALC    !dAlpha_dTau * Tau
!double precision :: A11_CALC    !d2Alpha_dDelta_dTau * Delta * Tau
double precision :: dalpha_ddelta
double precision :: d2alpha_ddelta2
double precision :: dalpha_dtau
double precision :: d2alpha_ddelta_dtau



!Initializing
A = 0.d0

dalpha_ddelta = 0.d0
d2alpha_ddelta2 = 0.d0
dalpha_dtau = 0.d0
d2alpha_ddelta_dtau = 0.d0



!Calculation of requested derivatives of the reference fluid
dalpha_ddelta = A01_CALC(gl,t0_loc, d0, 0)
d2alpha_ddelta2 = A02_CALC(gl,t0_loc, d0, 0)
dalpha_dtau = A10_CALC(gl,t0_loc, d0, 0)
d2alpha_ddelta_dtau = A11_CALC(gl,t0_loc, d0, 0)

!set matrix A
A(1, 1) = dalpha_dtau * (- 1.d0 / t0_loc)
A(1, 2) = dalpha_ddelta * (1.d0 / d0)
A(2, 1) = d2alpha_ddelta_dtau * (- 1.d0 / t0_loc)
A(2, 2) = (dalpha_ddelta + d2alpha_ddelta2) * (1.d0 / d0)


end subroutine
    
    
    
!Building the inverse of a square matrix
subroutine Get_Inverse(gl,A, A_Invertiert)


implicit none

    type(type_gl) :: gl



!Declaration
double precision, dimension(2, 2) :: A
double precision, dimension(2, 2) :: A_Invertiert
double precision :: Det


!Initializing
A_Invertiert = 0.d0
Det = 0.d0


!Inverting
Det = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
A_Invertiert(1, 1) = A(2, 2) / Det
A_Invertiert(1, 2) = -A(1, 2) / Det
A_Invertiert(2, 1) = -A(2, 1) / Det
A_Invertiert(2, 2) = A(1, 1) / Det


end subroutine
    
    
    
!Function calculating first shape parameter theta.
!See: Ely et al. (1981): "Prediction of Transport Properties. 1. Viscotity of Fluids and Mixtures"
double precision function Theta_calc(gl,T, Tc, Rho, Rhoc, w)



implicit none

    type(type_gl) :: gl



!Declaration
double precision :: T
double precision :: Tc
double precision :: Rho
double precision :: Rhoc
double precision :: w

double precision :: T_red
double precision :: v_red
double precision :: Tplus
double precision :: vplus
double precision :: F
double precision :: w_ref
double precision :: a_ref
double precision :: b_ref
double precision :: c_ref
double precision :: d_ref


!Initializing
T_red = 0.d0
v_red = 0.d0
Tplus = 0.d0
vplus = 0.d0
F = 0.d0
w_ref = 0.01142d0
a_ref = 0.090569d0
b_ref = -0.862762d0
c_ref = 0.316636d0
d_ref = -0.465684d0


!Calculation
T_red = T / Tc
v_red = Rhoc / Rho
if (T_red <= 0.5d0) then
    Tplus = 0.5d0
elseif (T_red >= 2.d0) then
    Tplus = 2.d0
else
    Tplus = T_red
end if
if (v_red <= 0.5d0) then
    vplus = 0.5d0
elseif (v_red >= 2.d0) then
    vplus = 2.d0
else
    vplus = v_red
end if
F = a_ref + b_ref * dlog(Tplus) + (c_ref + d_ref / Tplus) * (vplus - 0.5d0)
Theta_calc = 1.d0 + (w - w_ref) * F

end function
    
    
    
!Function calculating second shape parameter phi.
!See: Ely et al. (1981): "Prediction of Transport Properties. 1. Viscotity of Fluids and Mixtures"
double precision function Phi_calc(gl,T, Tc, Rho, Rhoc, w)


implicit none

    type(type_gl) :: gl



!Declaration
double precision :: T
double precision :: Tc
double precision :: Rho
double precision :: Rhoc
double precision :: w

double precision :: T_red
double precision :: v_red
double precision :: Tplus
double precision :: vplus
double precision :: G
double precision :: Z
double precision :: w_ref
double precision :: Z_ref
double precision :: a_ref
double precision :: b_ref
double precision :: c_ref
double precision :: d_ref


!Initializing
T_red = 0.d0
v_red = 0.d0
Tplus = 0.d0
vplus = 0.d0
G = 0.d0
Z = 0.d0
w_ref = 0.01142d0
Z_ref = 0.28953d0    !See Lee et al. (1975): "A Generalized Thermodynamic Correlation Based on Three-Parameter Corresponding States"
a_ref = 0.394301d0
b_ref = -1.023545d0
c_ref = -0.932813d0
d_ref = -0.754639d0


!Calculation
T_red = T / Tc
v_red = Rhoc / Rho
if (T_red <= 0.5d0) then
    Tplus = 0.5d0
elseif (T_red >= 2.d0) then
    Tplus = 2.d0
else
    Tplus = T_red
end if
if (v_red <= 0.5d0) then
    vplus = 0.5d0
elseif (v_red >= 2.d0) then
    vplus = 2.d0
else
    vplus = v_red
end if
G = a_ref * (vplus + b_ref) + c_ref * (vplus + d_ref) * dlog(Tplus)
Z = 0.2905d0 - 0.085d0 * w    !See Lee et al. (1975): "A Generalized Thermodynamic Correlation Based on Three-Parameter Corresponding States"
Phi_calc = (1.d0 + (w - w_ref) * G) * Z_ref / Z

end function

    



    end module ecs_module
