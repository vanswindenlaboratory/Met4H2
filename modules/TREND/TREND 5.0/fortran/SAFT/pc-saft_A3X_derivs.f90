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

    ! module for file pc-saft_A3X_derivs.f90
    module pc_saft_A3X_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module
    use pc_saft_ancillary_routines_module
    use variables_transformation_module

    use pc_saft_JX_derivs_module


    contains


   
subroutine A3X1DERIVS (gl,T, DENS, GETDERAQQ)

    ! a_3qq: thrid-order perturbation term of the quadrupolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.10 in Gross, Sadowski 2005:
    ! a_3qq = pi/3*(3/4)**3*rho*sum_i*sum_j*x_i*x_j*(espk_ii/T)**1.5*(espk_jj/T)**1.5*sig_ii**(15/2)*sig_jj**(15/2)/sig_ij**12*n_Qi*n_Qj*Q_i**3*Q_j**3*J_3ij + &       
    ! & 4 *pi**2/3*(3/4)**3*rho**2*sum_j*sum_i*sum_k(x_i*x_j*x_k*epsk_ii/T*espk_jj/T*espk_kk/T*sig_ii**5*sig_jj**5*sig_kk**5/sig_ij**3/sig_ik**3/sig_jk**3*n_Qi*n_Qj*n_Qk*Q_i**2*Q_j**2*Q_k**2*J_3ijk
    ! J3,ij = 0, which means that the first summation equals zero and is therefore left out
    ! dependent on D and T and x_i







implicit none

    type(type_gl) :: gl


     ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAQQ
    double precision, intent (in) :: T, DENS
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sigma_ij, sigma_jk, sigma_oj, sigma_ok, sigma_ik
    double precision :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16, sum17, sum18, sum19, sum20, sum21
    double precision :: part13, part3, part3a, part3b, part3c, sum, help, help2, help3, part0, part00, part1, part2, part4, part5, suma , sumb, sumc
    integer :: i, j, k, o
    double precision, dimension(2) :: qq

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
   ! !eq. 10 in Gross (2005). J3,ij = 0, which means that the first summation equals zero and is therefore left out
   ! !calculate the derivatives of a_3qq
    if (GETDERAQQ(1) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0   
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_oj**6 * gl%mPCSAFT(j))
                sum1 = sum1 - suma * gl%molfractions(o)**2 * gl%J3X1_PCSAFTQ(1,o,j,o,o)
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_jk**3 * sigma_ok**3 * sigma_oj**3 * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum3 = sum3 + sumb * gl%J3_PCSAFTQ(1,j,k,o)
                    sum4 = sum4 + sumb * gl%J3X1_PCSAFTQ(1,j,k,o,o)
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sum5 = sum5 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k) * gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 *gl%QPCSAFTQ(k)**2  &
                                & / (kbol**3 * T**3 * sigma_ij**3 * sigma_ik**3 * sigma_jk**3 * gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k)) * gl%J3X1_PCSAFTQ(1,i,j,k,o)
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTQ(o)**3 * gl%QPCSAFTQ(o)**6 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**9 * gl%mPCSAFT(o)**3) 
            part1 = 1.d0 * gl%molfractions(o)**3 * help * gl%J3X1_PCSAFTQ(1,o,o,o,o)
            
            gl%A3X1_PCSAFTQ(1,o) = piPCSAFT**2 * 0.5625d0 * DENS * DENS * 1.d-57 * (3.d0 * gl%QPCSAFTQ(o)**4 *gl%nPCSAFTQ(o)**2 / (T**3*kbol**3*gl%mPCSAFT(o)**2 * gl%sigPCSAFT(o)**3) * sum1  &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 *gl%nPCSAFTQ(o) / (T**3*kbol**3*gl%mPCSAFT(o)) * (sum3 + gl%molfractions(o)*sum4)  &
                & + part1 + sum5)
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
  
    !  2: 1ST DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires J_3ijk
    if (GETDERAQQ(2) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_oj**6 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTQ(2,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTQ(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_jk**3 * sigma_ok**3 * sigma_oj**3 * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum3 = sum3 + gl%J3X1_PCSAFTQ(2,j,k,o,o) * sumb
                    sum4 = sum4 + gl%J3_PCSAFTQ(2,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTQ(1,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTQ(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%QPCSAFTQ(i)**2*gl%QPCSAFTQ(j)**2*gl%QPCSAFTQ(k)**2*gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%nPCSAFTQ(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)/(T**3*kbol**3*sigma_ij**3*sigma_ik**3*sigma_jk**3*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                            sum5 = sum5 + gl%J3X1_PCSAFTQ(2,i,j,k,o) * sumc
                            sum8 = sum8 + gl%J3X1_PCSAFTQ(1,i,j,k,o) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTQ(o)**3 * gl%QPCSAFTQ(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**9 * gl%mPCSAFT(o)**3) 
            part1 = gl%J3X1_PCSAFTQ(2,o,o,o,o) * help
            part2 = gl%J3X1_PCSAFTQ(1,o,o,o,o) * help
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTQ(2,o) = 0.5625d0 * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
            & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum3 + Sum4)/ help3  &
            & + Sum5 + 2.d0 * (part2 - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
            & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 + Sum8))
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
    
    ! 3: 2ND DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires j_3ijk
    if (GETDERAQQ(3) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_oj**6 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTQ(1,o,j,o,o)*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(j)*gl%molfractions(j)/(sigma_oj**6*gl%mPCSAFT(j))
                sum2 = sum2 + gl%J3X1_PCSAFTQ(3,o,j,o,o)*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(j)*gl%molfractions(j)/(sigma_oj**6*gl%mPCSAFT(j))
                sum9 = sum9 + gl%J3X1_PCSAFTQ(2,o,j,o,o)*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(j)*gl%molfractions(j)/(sigma_oj**6*gl%mPCSAFT(j))
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_jk**3 * sigma_ok**3 * sigma_oj**3 * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum3 = sum3 + gl%J3X1_PCSAFTQ(1,j,k,o,o) * sumb
                    sum4 = sum4 + gl%J3_PCSAFTQ(1,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTQ(3,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTQ(3,j,k,o) * sumb
                    sum10 = sum10 + gl%J3X1_PCSAFTQ(2,j,k,o,o) * sumb
                    sum11 = sum11 + gl%J3_PCSAFTQ(2,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%nPCSAFTQ(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij**3 * sigma_ik**3 * sigma_jk**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum5 = sum5 + gl%J3X1_PCSAFTQ(3,i,j,k,o) * sumc
                            sum8 = sum8 + gl%J3X1_PCSAFTQ(2,i,j,k,o) * sumc
                            sum12 = sum12 + gl%J3X1_PCSAFTQ(1,i,j,k,o) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTQ(o)**3 * gl%QPCSAFTQ(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**9 * gl%mPCSAFT(o)**3) 
            part1 = gl%J3X1_PCSAFTQ(1,o,o,o,o) * help
            part2 = gl%J3X1_PCSAFTQ(3,o,o,o,o) * help
            part3 = gl%J3X1_PCSAFTQ(2,o,o,o,o) * help
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTQ(3,o) = piPCSAFT**2 * 0.5625d0 * DENS * DENS * 1.d-57 * (2.d0 * part1 &
             & - 6.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
             & + 6.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum3 + Sum4) / help3 &
             & + part2 - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
             & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 &
             & + Sum5 + 4.d0 *(part3 - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum9 / help2 &
             & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum10 + Sum11) / help3 &
             & + Sum8) + 2.d0 * Sum12)
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if

    
    
   ! ! 4: 1ST DERIVATIVE OF a_3qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
   ! ! requires j_3ijk
    if (GETDERAQQ(4) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_oj**6 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTQ(4,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTQ(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_jk**3 * sigma_ok**3 * sigma_oj**3 * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum3 = sum3 + gl%J3X1_PCSAFTQ(4,j,k,o,o) * sumb
                    sum4 = sum4 + gl%J3_PCSAFTQ(4,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTQ(1,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTQ(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%nPCSAFTQ(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij**3 * sigma_ik**3 * sigma_jk**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum5 = sum5 + (-3.d0 * gl%J3X1_PCSAFTQ(1,i,j,k,o) + gl%J3X1_PCSAFTQ(4,i,j,k,o)) * sumc ! T für Jq_3_ijko_xo
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTQ(o)**3 * gl%QPCSAFTQ(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**9 * gl%mPCSAFT(o)**3) 
            part1 = - 3.d0 * gl%J3X1_PCSAFTQ(1,o,o,o,o) * help ! T
            part2 = gl%J3X1_PCSAFTQ(4,o,o,o,o) * help
            
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTQ(4,o) =  0.5625d0 * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 + part2 + Sum5 &
                & - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
                & + 9.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 & !T
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum3 + Sum4) / help3 &
                & - 9.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 ) !T
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if

   
    
    
    ! 5: 2ND DERIVATIVE OF a_3qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires j_3ijk
    if (GETDERAQQ(5) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_oj**6 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTQ(5,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTQ(4,o,j,o,o) * suma
                sum3 = sum3 + gl%J3X1_PCSAFTQ(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_jk**3 * sigma_ok**3 * sigma_oj**3 * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum4 = sum4 + gl%J3X1_PCSAFTQ(5,j,k,o,o) * sumb
                    sum5 = sum5 + gl%J3_PCSAFTQ(5,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTQ(4,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTQ(4,j,k,o) * sumb
                    sum8 = sum8 + gl%J3X1_PCSAFTQ(1,j,k,o,o) * sumb
                    sum9 = sum9 + gl%J3_PCSAFTQ(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%nPCSAFTQ(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij**3 * sigma_ik**3 * sigma_jk**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum10 = sum10 + (12.d0 * gl%J3X1_PCSAFTQ(1,i,j,k,o) - 6.d0 * gl%J3X1_PCSAFTQ(4,i,j,k,o) + gl%J3X1_PCSAFTQ(5,i,j,k,o)) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTQ(o)**3 * gl%QPCSAFTQ(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**9 * gl%mPCSAFT(o)**3) 
            part1 = (12.d0 * gl%J3X1_PCSAFTQ(1,o,o,o,o) - 6.d0 * gl%J3X1_PCSAFTQ(4,o,o,o,o) + gl%J3X1_PCSAFTQ(5,o,o,o,o)) * help
                        
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTQ(5,o) =  0.5625d0 * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 + Sum10 &
                & - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
                & + 18.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & - 36.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum3 / help2 &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum4 + Sum5) / help3 &
                & - 18.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 &
                & + 36.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum8 + Sum9) / help3 )
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
  
    
    ! 6: 2ND MIXED DERIVATIVE OF a_3 WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_3ijk
    if (GETDERAQQ(6) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        sum13 = 0.d0
        sum14 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_oj**6 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTQ(6,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTQ(2,o,j,o,o) * suma
                sum3 = sum3 + gl%J3X1_PCSAFTQ(4,o,j,o,o) * suma
                sum8 = sum8 + gl%J3X1_PCSAFTQ(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_jk**3 * sigma_ok**3 * sigma_oj**3 * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum4 = sum4 + gl%J3X1_PCSAFTQ(6,j,k,o,o) * sumb
                    sum5 = sum5 + gl%J3_PCSAFTQ(6,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTQ(2,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTQ(2,j,k,o) * sumb
                    sum9 = sum9 + gl%J3X1_PCSAFTQ(4,j,k,o,o) * sumb
                    sum10 = sum10 + gl%J3_PCSAFTQ(4,j,k,o) * sumb
                    sum11 = sum11 + gl%J3X1_PCSAFTQ(1,j,k,o,o) * sumb
                    sum12 = sum12 + gl%J3_PCSAFTQ(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%nPCSAFTQ(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij**3 * sigma_ik**3 * sigma_jk**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum13 = sum13 + (gl%J3X1_PCSAFTQ(6,i,j,k,o) - 3.d0 * gl%J3X1_PCSAFTQ(2,i,j,k,o)) * sumc
                            sum14 = sum14 + (-3.d0 * gl%J3X1_PCSAFTQ(1,i,j,k,o) + gl%J3X1_PCSAFTQ(4,i,j,k,o)) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTQ(o)**3 * gl%QPCSAFTQ(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**9 * gl%mPCSAFT(o)**3) 
            part1 = gl%J3X1_PCSAFTQ(6,o,o,o,o) * help
            part2 = 3.d0 * gl%J3X1_PCSAFTQ(2,o,o,o,o) * help
            part3 = - 3.d0 * gl%J3X1_PCSAFTQ(1,o,o,o,o) * help
            part4 = gl%J3X1_PCSAFTQ(4,o,o,o,o) * help
                        
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTQ(6,o) =  0.5625d0 * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 - part2  &
                & - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
                & + 9.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum4 + Sum5) / help3 &
                & - 9.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 &
                & + Sum13 + 2.d0 * (part3 + part4 &
                & - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum3 / help2 &
                & + 9.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum8 / help2 &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum9 + Sum10) / help3 &
                & - 9.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum11 + Sum12) / help3 &
                & + Sum14))
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_3qq WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires j_3ijk
    if (GETDERAQQ(7) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        sum13 = 0.d0
        sum14 = 0.d0
        sum15 = 0.d0
        sum16 = 0.d0
        sum17 = 0.d0
        sum18 = 0.d0
        sum19 = 0.d0
        sum20 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_oj**6 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTQ(7,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTQ(6,o,j,o,o) * suma
                sum3 = sum3 + gl%J3X1_PCSAFTQ(2,o,j,o,o) * suma
                sum8 = sum8 + gl%J3X1_PCSAFTQ(5,o,j,o,o) * suma
                sum15 = sum15 + gl%J3X1_PCSAFTQ(4,o,j,o,o) * suma
                sum16 = sum16 + gl%J3X1_PCSAFTQ(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_jk**3 * sigma_ok**3 * sigma_oj**3 * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum4 = sum4 + gl%J3X1_PCSAFTQ(7,j,k,o,o) * sumb
                    sum5 = sum5 + gl%J3_PCSAFTQ(7,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTQ(6,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTQ(6,j,k,o) * sumb
                    sum9 = sum9 + gl%J3X1_PCSAFTQ(2,j,k,o,o) * sumb
                    sum10 = sum10 + gl%J3_PCSAFTQ(2,j,k,o) * sumb
                    sum11 = sum11 + gl%J3X1_PCSAFTQ(5,j,k,o,o) * sumb
                    sum12 = sum12 + gl%J3_PCSAFTQ(5,j,k,o) * sumb
                    sum17 = sum17 + gl%J3X1_PCSAFTQ(4,j,k,o,o) * sumb
                    sum18 = sum18 + gl%J3_PCSAFTQ(4,j,k,o) * sumb
                    sum19 = sum19 + gl%J3X1_PCSAFTQ(1,j,k,o,o) * sumb
                    sum20 = sum20 + gl%J3_PCSAFTQ(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%nPCSAFTQ(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij**3 * sigma_ik**3 * sigma_jk**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum13 = sum13 + (gl%J3X1_PCSAFTQ(7,i,j,k,o) - 6.d0 * gl%J3X1_PCSAFTQ(6,i,j,k,o) + 12.d0 * gl%J3X1_PCSAFTQ(2,i,j,k,o)) * sumc
                            sum14 = sum14 + (12.d0 * gl%J3X1_PCSAFTQ(1,i,j,k,o) - 6.d0 * gl%J3X1_PCSAFTQ(4,i,j,k,o) + gl%J3X1_PCSAFTQ(5,i,j,k,o)) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTQ(o)**3 * gl%QPCSAFTQ(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**9 * gl%mPCSAFT(o)**3) 
            part1 = (gl%J3X1_PCSAFTQ(7,o,o,o,o) - 6.d0 * gl%J3X1_PCSAFTQ(6,o,o,o,o) + 12.d0 * gl%J3X1_PCSAFTQ(2,o,o,o,o)) * help
            part2 = (12.d0 * gl%J3X1_PCSAFTQ(1,o,o,o,o) - 6.d0 * gl%J3X1_PCSAFTQ(4,o,o,o,o) + gl%J3X1_PCSAFTQ(5,o,o,o,o)) * help
                        
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTQ(7,o) =  0.5625d0 * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 &
                & - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * sum1 / help2 &
                & + 18.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & - 36.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum3 / help2 &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum4 + Sum5) / help3 &
                & - 18.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 &
                & + 36.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum9 + Sum10) / help3 &
                & + Sum13 + 2.d0 * (part2 &
                & - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum8 / help2 &
                & + 18.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum15 / help2 &
                & - 36.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum16 / help2 &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum11 + Sum12) / help3 &
                & - 18.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum17 + Sum18) / help3 &
                & + 36.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum19 + Sum20) / help3 &
                & + Sum14))
                ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if   
    
    
    ! 8: 3RD DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires j_3ijk
    if (GETDERAQQ(8) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_oj**6 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTQ(2,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTQ(8,o,j,o,o) * suma
                sum9 = sum9 + gl%J3X1_PCSAFTQ(3,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) &
                        & / (sigma_jk**3 * sigma_ok**3 * sigma_oj**3 * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum3 = sum3 + gl%J3X1_PCSAFTQ(2,j,k,o,o) * sumb
                    sum4 = sum4 + gl%J3_PCSAFTQ(2,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTQ(8,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTQ(8,j,k,o) * sumb
                    sum10 = sum10 + gl%J3X1_PCSAFTQ(3,j,k,o,o) * sumb
                    sum11 = sum11 + gl%J3_PCSAFTQ(3,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%nPCSAFTQ(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij**3 * sigma_ik**3 * sigma_jk**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum5 = sum5 + gl%J3X1_PCSAFTQ(8,i,j,k,o) * sumc
                            sum8 = sum8 + gl%J3X1_PCSAFTQ(3,i,j,k,o) * sumc
                            sum12 = sum12 + gl%J3X1_PCSAFTQ(2,i,j,k,o) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTQ(o)**3 * gl%QPCSAFTQ(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**9 * gl%mPCSAFT(o)**3) 
            part1 = 6.d0 * gl%J3X1_PCSAFTQ(2,o,o,o,o) * help
            part2 = gl%J3X1_PCSAFTQ(8,o,o,o,o) * help
            part3 = gl%J3X1_PCSAFTQ(3,o,o,o,o) * help
            
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTQ(8,o) = piPCSAFT**2 * 0.5625d0 * DENS * DENS * 1.d-57 * (part1 &
                & - 18.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum1 / help2  &
                & + 18.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum3 + Sum4) / help3 &
                & + part2 - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3  &
                & + Sum5 + 6.d0 * (part3 - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum9 / help2 &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum10 + Sum11) / help3 &
                & + Sum8) + 6.d0 * Sum12)
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
  
     
   ! 9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires 
    if (GETDERAQQ(9) .eq. 1) then
      do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        sum13 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_oj**6 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTQ(9,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTQ(5,o,j,o,o) * suma
                sum3 = sum3 + gl%J3X1_PCSAFTQ(4,o,j,o,o) * suma
                sum4 = sum4 + gl%J3X1_PCSAFTQ(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_jk**3 * sigma_ok**3 * sigma_oj**3 * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum5 = sum5 + gl%J3X1_PCSAFTQ(9,j,k,o,o) * sumb
                    sum6 = sum6 + gl%J3_PCSAFTQ(9,j,k,o) * sumb
                    sum7 = sum7 + gl%J3X1_PCSAFTQ(5,j,k,o,o) * sumb
                    sum8 = sum8 + gl%J3_PCSAFTQ(5,j,k,o) * sumb
                    sum9 = sum9 + gl%J3X1_PCSAFTQ(4,j,k,o,o) * sumb
                    sum10 = sum10 + gl%J3_PCSAFTQ(4,j,k,o) * sumb
                    sum11 = sum11 + gl%J3X1_PCSAFTQ(1,j,k,o,o) * sumb
                    sum12 = sum12 + gl%J3_PCSAFTQ(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%nPCSAFTQ(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij**3 * sigma_ik**3 * sigma_jk**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum13 = sum13 + (-60.d0 * gl%J3X1_PCSAFTQ(1,i,j,k,o) + 36.d0 * gl%J3X1_PCSAFTQ(4,i,j,k,o) - 9.d0 * gl%J3X1_PCSAFTQ(5,i,j,k,o) + gl%J3X1_PCSAFTQ(9,i,j,k,o)) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTQ(o)**3 * gl%QPCSAFTQ(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**9 * gl%mPCSAFT(o)**3) 
            part1 = (-60.d0 * gl%J3X1_PCSAFTQ(1,o,o,o,o) + 36.d0 * gl%J3X1_PCSAFTQ(4,o,o,o,o) - 9.d0 * gl%J3X1_PCSAFTQ(5,o,o,o,o) + gl%J3X1_PCSAFTQ(9,o,o,o,o)) * help
            
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTQ(9,o) =  0.5625d0 * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 &
                & - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
                & + 27.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & - 108.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum3 / help2 &
                & + 180.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum4 / help2 &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum5 + Sum6) / help3 &
                & - 27.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum7 + Sum8) / help3 &
                & + 108.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum9 + Sum10) / help3 &
                & - 180.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum11 + Sum12) / help3 &
                & + Sum13)
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if

    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires j_3ijk
    if (GETDERAQQ(10) .eq. 1) then
    do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        sum13 = 0.d0
        sum14 = 0.d0
        sum15 = 0.d0
        sum16 = 0.d0
        sum17 = 0.d0
        sum18 = 0.d0
        sum19 = 0.d0
        sum20 = 0.d0
        sum21 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_oj**6 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTQ(4,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTQ(1,o,j,o,o) * suma
                sum3 = sum3 + gl%J3X1_PCSAFTQ(10,o,j,o,o) * suma
                sum15 = sum15 + gl%J3X1_PCSAFTQ(3,o,j,o,o) * suma
                sum16 = sum16 + gl%J3X1_PCSAFTQ(6,o,j,o,o) * suma
                sum17 = sum17 + gl%J3X1_PCSAFTQ(2,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) &
                        & / (sigma_jk**3 * sigma_ok**3 * sigma_oj**3 * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum4 = sum4 + gl%J3X1_PCSAFTQ(4,j,k,o,o) * sumb
                    sum5 = sum5 + gl%J3_PCSAFTQ(4,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTQ(1,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTQ(1,j,k,o) * sumb
                    sum8 = sum8 + gl%J3X1_PCSAFTQ(10,j,k,o,o) * sumb
                    sum9 = sum9 + gl%J3_PCSAFTQ(10,j,k,o) * sumb
                    sum10 = sum10 + gl%J3X1_PCSAFTQ(3,j,k,o,o) * sumb
                    sum11 = sum11 + gl%J3_PCSAFTQ(3,j,k,o) * sumb
                    sum12 = sum12 + gl%J3X1_PCSAFTQ(6,j,k,o,o) * sumb
                    sum13 = sum13 + gl%J3_PCSAFTQ(6,j,k,o) * sumb
                    sum18 = sum18 + gl%J3X1_PCSAFTQ(2,j,k,o,o) * sumb
                    sum19 = sum19 + gl%J3_PCSAFTQ(2,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%nPCSAFTQ(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij**3 * sigma_ik**3 * sigma_jk**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                           
                            sum14 = sum14 + (gl%J3X1_PCSAFTQ(10,i,j,k,o) - 3.d0 * gl%J3X1_PCSAFTQ(3,i,j,k,o)) * sumc
                            sum20 = sum20 + (gl%J3X1_PCSAFTQ(6,i,j,k,o) - 3.d0 * gl%J3X1_PCSAFTQ(2,i,j,k,o)) * sumc
                            sum21 = sum21 + (-3.d0 * gl%J3X1_PCSAFTQ(1,i,j,k,o) + gl%J3X1_PCSAFTQ(4,i,j,k,o)) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTQ(o)**3 * gl%QPCSAFTQ(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**9 * gl%mPCSAFT(o)**3) 
            part1 = (-6.d0 * gl%J3X1_PCSAFTQ(1,o,o,o,o) + 2.d0 * gl%J3X1_PCSAFTQ(4,o,o,o,o)) * help
            part2 = (gl%J3X1_PCSAFTQ(10,o,o,o,o) - 3.d0 * gl%J3X1_PCSAFTQ(3,o,o,o,o)) * help
            part3 = (gl%J3X1_PCSAFTQ(6,o,o,o,o) - 3.d0 * gl%J3X1_PCSAFTQ(2,o,o,o,o)) * help
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTQ(10,o) = piPCSAFT**2 * 0.5625d0 * DENS * DENS * 1.d-57 * (part1  &
                & - 6.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 & 
                & + 18.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & + 6.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum4 + Sum5) / help3 &
                & - 18.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 &
                & + part2 &
                & - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum3 / help2 &
                & + 9.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum15 / help2 &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum8 + Sum9) / help3 &
                & - 9.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum10 + Sum11) / help3 &
                & + Sum14 &
                & + 4.d0 * (part3 - 3.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum16 / help2 &
                & + 9.d0 * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 * gl%molfractions(o)**2 * Sum17 / help2 &
                & + 3.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum12 + Sum13) / help3 &
                & - 9.d0 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * (gl%molfractions(o) * Sum18 + Sum19) / help3 &
                & + Sum20) + 2.d0 * Sum21)
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
        
!DEC$ END IF
end subroutine A3X1DERIVS


subroutine A3X2DERIVS(gl,T,DENS,GETDERAQQ)
    ! a_3dd: second-order perturbation term of the dipolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.9 in Gross, Sadowski 2006:
    ! a_3dd =-4*pi**2/3*rho**2*sum_i*sum_j*sum_k(x_i*x_j*x_k*epsk_ii/T*espk_jj/T*espk_kk/T*sig_ii**3*sig_jj**3*sig_kk**3/sig_ij/sig_ik/sig_jk*n_Di*n_Dj*n_Dk*My_i**2*My_j**2*My_k**2*J3_D_ijk
    ! dependent on D and T and xi, xj and xk
   




   

implicit none

    type(type_gl) :: gl

   
    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAQQ
    double precision, intent (in) :: T, DENS
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sigma_ij, sigma_jk, sigma_ik, sigma_ou, sigma_io, sigma_iu, sigma_jo, sigma_ju, sigma_ku
    double precision :: part1, sum1a, sum1b, sum1c, sum1d, sum2a, sum2b, sum3a, sum3b, sum3c, sum3d, sum4a, sum4b, sum4c, sum4d, sum4e, sum4f, sum4g, sum4h, sum, suma, sumb
    integer :: i, j, xi, k, o, u
    double precision :: DP, DP2, DM, DM2, TP, TP2, TM, TM2
    double precision :: DPx, DP2x, DMx, DM2x, TPx, TP2x, TMx, TM2x
    double precision :: DPxx, DP2xx, DMxx, DM2xx, TPxx, TP2xx, TMxx, TM2xx
    double precision :: DELTA, DEL, DELT
    double precision ::  A_u, A_o, A_ou, A_oij, A_iou, A_i, A_j, A_k, A_io, A_ju, A_ku, A_iu 
    double precision, dimension(gl%ncomp,gl%ncomp) :: A3X2, A3X2DP, A3X2DP2, A3X2DM, A3X2DM2, A3X2TP, A3X2TP2, A3X2TM, A3X2TM2, A3X2PP, A3X2PM, A3X2MP, A3X2MM
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    DELTA = 1.d-7
    
    DP = DENS * (1.d0 + DELTA)
    DP2 = DENS * (1.d0 + 2.d0*DELTA)
    DM = DENS * (1.d0 - DELTA)
    DM2 = DENS * (1.d0 - 2.d0*DELTA)
    TP = T * (1.d0 + DELTA)
    TP2 = T * (1.d0 + 2.d0*DELTA)
    TM = T * (1.d0 - DELTA)
    TM2 = T * (1.d0 - 2.d0*DELTA)
    
    DEL = 5.d-4 !gut für 2. Temperauturableitung 
    
    DPx = DENS * (1.d0 + DEL)
    DP2x = DENS * (1.d0 + 2.d0*DEL)
    DMx = DENS * (1.d0 - DEL)
    DM2x = DENS * (1.d0 - 2.d0*DEL)
    TPx = T * (1.d0 + DEL)
    TP2x = T * (1.d0 + 2.d0*DEL)
    TMx = T * (1.d0 - DEL)
    TM2x = T * (1.d0 - 2.d0*DEL)
    
    DELT = 5.d-4 ! Gemischte Ableitung & 2. Dichteableitung
    
    DPxx = DENS * (1.d0 + DELT)
    DP2xx = DENS * (1.d0 + 2.d0*DELT)
    DMxx = DENS * (1.d0 - DELT)
    DM2xx = DENS * (1.d0 - 2.d0*DELT)
    TPxx = T * (1.d0 + DELT)
    TP2xx = T * (1.d0 + 2.d0*DELT)
    TMxx = T * (1.d0 - DELT)
    TM2xx = T * (1.d0 - 2.d0*DELT)
   
    !calculate the derivatives of a_qq
    ! 1: a_2qq
    if (GETDERAQQ(1) .eq. 1) then

        call base_A3X2(gl,T,DENS,A3X2)
        
        gl%A3X2_PCSAFTQ(1,:,:) = A3X2
    end if
    
    !  2: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires j_2qq
    if (GETDERAQQ(2) .eq. 1) then  
    
        call base_A3X2(gl,T,DP,A3X2DP)
        call base_A3X2(gl,T,DP2,A3X2DP2)
        call base_A3X2(gl,T,DM,A3X2DM)
        call base_A3X2(gl,T,DM2,A3X2DM2)
        
        gl%A3X2_PCSAFTQ(2,:,:) = (8.d0*A3X2DP - 8.d0*A3X2DM - A3X2DP2 + A3X2DM2) / (12.d0 * DELTA) 
    end if
   
   ! 3: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
   ! requires j_2qq
    if (GETDERAQQ(3) .eq. 1) then
    
        call base_A3X2(gl,T,DENS,A3X2)    
        call base_A3X2(gl,T,DPxx,A3X2DP)
        call base_A3X2(gl,T,DP2xx,A3X2DP2)
        call base_A3X2(gl,T,DMxx,A3X2DM)
        call base_A3X2(gl,T,DM2xx,A3X2DM2)
        
        gl%A3X2_PCSAFTQ(3,:,:) = (-A3X2DM2 + 16.d0*A3X2DM - 30.d0*A3X2 + 16.d0*A3X2DP - A3X2DP2) / (12.d0*DELT**2) 
    
    end if
     
    ! 4: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires j_2qq
    if (GETDERAQQ(4) .eq. 1) then
    
        call base_A3X2(gl,TP,DENS,A3X2TP)
        call base_A3X2(gl,TP2,DENS,A3X2TP2)
        call base_A3X2(gl,TM,DENS,A3X2TM)
        call base_A3X2(gl,TM2,DENS,A3X2TM2)
        
        gl%A3X2_PCSAFTQ(4,:,:) = (8.d0*A3X2TP - 8.d0*A3X2TM - A3X2TP2 + A3X2TM2) / (12.d0 * DELTA) 
    end if
   
    ! 5: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    !requires j_2qq
    if (GETDERAQQ(5) .eq. 1) then
    
        call base_A3X2(gl,T,DENS,A3X2)    
        call base_A3X2(gl,TPx,DENS,A3X2TP)
        call base_A3X2(gl,TP2x,DENS,A3X2TP2)
        call base_A3X2(gl,TMx,DENS,A3X2TM)
        call base_A3X2(gl,TM2x,DENS,A3X2TM2)
        
        gl%A3X2_PCSAFTQ(5,:,:) = (-A3X2TM2 + 16.d0*A3X2TM - 30.d0*A3X2 + 16.d0*A3X2TP - A3X2TP2) / (12.d0*DEL**2) 
    end if
       
    ! 6: 2ND MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_2qq
    if (GETDERAQQ(6) .eq. 1) then
    
        call base_A3X2(gl,TP,DPxx,A3X2PP)
        call base_A3X2(gl,TM,DPxx,A3X2MP)
        call base_A3X2(gl,TM,DMxx,A3X2MM)
        call base_A3X2(gl,TP,DMxx,A3X2PM)
        
        gl%A3X2_PCSAFTQ(6,:,:) = (A3X2PP - A3X2PM - A3X2MP + A3X2MM) / (4.d0*DELT*DELTA) 
    end if
    
!DEC$ END IF
end subroutine A3X2DERIVS   


subroutine A3X2DERIVS_D(gl,T,DENS,GETDERADD)
    ! a_3dd: second-order perturbation term of the dipolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.9 in Gross, Sadowski 2006:
    ! a_3dd =-4*pi**2/3*rho**2*sum_i*sum_j*sum_k(x_i*x_j*x_k*epsk_ii/T*espk_jj/T*espk_kk/T*sig_ii**3*sig_jj**3*sig_kk**3/sig_ij/sig_ik/sig_jk*n_Di*n_Dj*n_Dk*My_i**2*My_j**2*My_k**2*J3_D_ijk
    ! dependent on D and T and xi, xj and xk
   




   

implicit none

    type(type_gl) :: gl

   
  ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERADD
    double precision, intent (in) :: T, DENS
    !output: add_PCSAFT (module variable)
    !working variables
    double precision :: sigma_ij, sigma_jk, sigma_ik, sigma_ou, sigma_io, sigma_iu, sigma_jo, sigma_ju, sigma_ku
    double precision :: part1, sum1a, sum1b, sum1c, sum1d, sum2a, sum2b, sum3a, sum3b, sum3c, sum3d, sum4a, sum4b, sum4c, sum4d, sum4e, sum4f, sum4g, sum4h, sum, suma, sumb
    integer :: i, j, xi, k, o, u
    double precision :: DP, DP2, DM, DM2, TP, TP2, TM, TM2
    double precision :: DPx, DP2x, DMx, DM2x, TPx, TP2x, TMx, TM2x
    double precision :: DPxx, DP2xx, DMxx, DM2xx, TPxx, TP2xx, TMxx, TM2xx
    double precision :: DELTA, DEL, DELT
    double precision :: A_u, A_o, A_ou, A_oij, A_iou, A_i, A_j, A_k, A_io, A_ju, A_ku, A_iu 
    double precision, dimension(gl%ncomp,gl%ncomp) :: A3X2, A3X2DP, A3X2DP2, A3X2DM, A3X2DM2, A3X2TP, A3X2TP2, A3X2TM, A3X2TM2, A3X2PP, A3X2PM, A3X2MP, A3X2MM
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    DELTA = 1.d-7
    
    DP = DENS * (1.d0 + DELTA)
    DP2 = DENS * (1.d0 + 2.d0*DELTA)
    DM = DENS * (1.d0 - DELTA)
    DM2 = DENS * (1.d0 - 2.d0*DELTA)
    TP = T * (1.d0 + DELTA)
    TP2 = T * (1.d0 + 2.d0*DELTA)
    TM = T * (1.d0 - DELTA)
    TM2 = T * (1.d0 - 2.d0*DELTA)
    
     
    DEL = 5.d-4 !gut für 2. Temperauturableitung & 2. Dichteableitung
    
    DPx = DENS * (1.d0 + DEL)
    DP2x = DENS * (1.d0 + 2.d0*DEL)
    DMx = DENS * (1.d0 - DEL)
    DM2x = DENS * (1.d0 - 2.d0*DEL)
    TPx = T * (1.d0 + DEL)
    TP2x = T * (1.d0 + 2.d0*DEL)
    TMx = T * (1.d0 - DEL)
    TM2x = T * (1.d0 - 2.d0*DEL)
    
    DELT = 5.d-4 ! Gemischte Ableitung
    
    DPxx = DENS * (1.d0 + DELT)
    DP2xx = DENS * (1.d0 + 2.d0*DELT)
    DMxx = DENS * (1.d0 - DELT)
    DM2xx = DENS * (1.d0 - 2.d0*DELT)
    TPxx = T * (1.d0 + DELT)
    TP2xx = T * (1.d0 + 2.d0*DELT)
    TMxx = T * (1.d0 - DELT)
    TM2xx = T * (1.d0 - 2.d0*DELT)
   
    !calculate the derivatives of a_dd
    ! 1: a_3dd
    if (GETDERADD(1) .eq. 1) then

        call base_A3X2_D(gl,T,DENS,A3X2)   
        
        gl%A3X2_PCSAFTD(1,:,:) = A3X2
        
   end if
    
    
    !  2: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires j_3qq
    if (GETDERADD(2) .eq. 1) then   
               
        call base_A3X2_D(gl,T,DP,A3X2DP)
        call base_A3X2_D(gl,T,DP2,A3X2DP2)
        call base_A3X2_D(gl,T,DM,A3X2DM)
        call base_A3X2_D(gl,T,DM2,A3X2DM2)
        
        gl%A3X2_PCSAFTD(2,:,:) = (8.d0*A3X2DP - 8.d0*A3X2DM - A3X2DP2 + A3X2DM2) / (12.d0 * DELTA) 
    
    end if
    
   
   ! 3: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
   ! requires j_3qq
    if (GETDERADD(3) .eq. 1) then
        
        call base_A3X2_D(gl,T,DENS,A3X2)    
        call base_A3X2_D(gl,T,DPx,A3X2DP)
        call base_A3X2_D(gl,T,DP2x,A3X2DP2)
        call base_A3X2_D(gl,T,DMx,A3X2DM)
        call base_A3X2_D(gl,T,DM2x,A3X2DM2)
        
        gl%A3X2_PCSAFTD(3,:,:) = (-A3X2DM2 + 16.d0*A3X2DM - 30.d0*A3X2 + 16.d0*A3X2DP - A3X2DP2) / (12.d0*DEL**2) 
    
    end if
    
   
     
    ! 4: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires j_3qq
    if (GETDERADD(4) .eq. 1) then
        call base_A3X2_D(gl,TP,DENS,A3X2TP)
        call base_A3X2_D(gl,TP2,DENS,A3X2TP2)
        call base_A3X2_D(gl,TM,DENS,A3X2TM)
        call base_A3X2_D(gl,TM2,DENS,A3X2TM2)
        
        gl%A3X2_PCSAFTD(4,:,:) = (8.d0*A3X2TP - 8.d0*A3X2TM - A3X2TP2 + A3X2TM2) / (12.d0 * DELTA) 
    end if
    
   
    ! 5: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    !requires j_3qq
    if (GETDERADD(5) .eq. 1) then
    
        call base_A3X2_D(gl,T,DENS,A3X2)    
        call base_A3X2_D(gl,TPx,DENS,A3X2TP)
        call base_A3X2_D(gl,TP2x,DENS,A3X2TP2)
        call base_A3X2_D(gl,TMx,DENS,A3X2TM)
        call base_A3X2_D(gl,TM2x,DENS,A3X2TM2)
        
        gl%A3X2_PCSAFTD(5,:,:) = (-A3X2TM2 + 16.d0*A3X2TM - 30.d0*A3X2 + 16.d0*A3X2TP - A3X2TP2) / (12.d0*DEL**2) 
    end if
       
     
    ! 6: 2ND MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_3qq
    if (GETDERADD(6) .eq. 1) then
    
        call base_A3X2_D(gl,TPxx,DPxx,A3X2PP)
        call base_A3X2_D(gl,TMxx,DPxx,A3X2MP)
        call base_A3X2_D(gl,TMxx,DMxx,A3X2MM)
        call base_A3X2_D(gl,TPxx,DMxx,A3X2PM)
        
        gl%A3X2_PCSAFTD(6,:,:) = (A3X2PP - A3X2PM - A3X2MP + A3X2MM) / (4.d0*DELT**2) 
    end if
!DEC$ END IF
end subroutine A3X2DERIVS_D


subroutine A3X1DERIVS_D (gl,T, DENS, GETDERADD)
    
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
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sigma_ij, sigma_jk, sigma_oj, sigma_ok, sigma_ik
    double precision :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16, sum17, sum18, sum19, sum20, sum21
    double precision :: part13, part3, part3a, part3b, part3c, sum, help, help2, help3, part0, part00, part1, part2, part4, part5, suma , sumb, sumc
    integer :: i, j, k, o
    double precision, dimension(2) :: dd
   
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
   ! !eq. 10 in Gross (2005). J3,ij = 0, which means that the first summation equals zero and is therefore left out
   ! !calculate the derivatives of a_3qq
    if (GETDERADD(1) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0   
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_oj**2 * gl%mPCSAFT(j))
                sum1 = sum1 - suma * gl%molfractions(o)**2 * gl%J3X1_PCSAFTD(1,o,j,o,o)
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_jk * sigma_ok * sigma_oj * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum3 = sum3 + sumb * gl%J3_PCSAFTD(1,j,k,o)
                    sum4 = sum4 + sumb * gl%J3X1_PCSAFTD(1,j,k,o,o)
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sum5 = sum5 + gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k) * gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k) * gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 *gl%MyPCSAFTD(k)**2  &
                                & / (kbol**3 * T**3 * sigma_ij * sigma_ik * sigma_jk * gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k)) * gl%J3X1_PCSAFTD(1,i,j,k,o)
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTD(o)**3 * gl%MyPCSAFTD(o)**6 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**3) 
            part1 = 1.d0 * gl%molfractions(o)**3 * help * gl%J3X1_PCSAFTD(1,o,o,o,o)
            
            gl%A3X1_PCSAFTD(1,o) = -piPCSAFT**2 * (4.d0/3.d0) * DENS * DENS * 1.d-57 * (3.d0 * gl%MyPCSAFTD(o)**4 *gl%nPCSAFTD(o)**2 / (T**3*kbol**3*gl%mPCSAFT(o)**2 * gl%sigPCSAFT(o)) * sum1  &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 *gl%nPCSAFTD(o) / (T**3*kbol**3*gl%mPCSAFT(o)) * (sum3 + gl%molfractions(o)*sum4)  &
                & + part1 + sum5)
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
  
    !  2: 1ST DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires J_3ijk
    if (GETDERADD(2) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_oj**2 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTD(2,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTD(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_jk * sigma_ok * sigma_oj * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum3 = sum3 + gl%J3X1_PCSAFTD(2,j,k,o,o) * sumb
                    sum4 = sum4 + gl%J3_PCSAFTD(2,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTD(1,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTD(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%nPCSAFTD(k)*gl%molfractions(i)*gl%molfractions(j)*gl%molfractions(k)/(T**3*kbol**3*sigma_ij*sigma_ik*sigma_jk*gl%mPCSAFT(i)*gl%mPCSAFT(j)*gl%mPCSAFT(k))
                            sum5 = sum5 + gl%J3X1_PCSAFTD(2,i,j,k,o) * sumc
                            sum8 = sum8 + gl%J3X1_PCSAFTD(1,i,j,k,o) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTD(o)**3 * gl%MyPCSAFTD(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**3) 
            part1 = gl%J3X1_PCSAFTD(2,o,o,o,o) * help
            part2 = gl%J3X1_PCSAFTD(1,o,o,o,o) * help
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o) * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTD(2,o) = -(4.d0/3.d0)  * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
            & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum3 + Sum4)/ help3  &
            & + Sum5 + 2.d0 * (part2 - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
            & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 + Sum8))
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
    
    ! 3: 2ND DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
    ! requires j_3ijk
    if (GETDERADD(3) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_oj**2 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTD(1,o,j,o,o)*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(j)*gl%molfractions(j)/(sigma_oj**2*gl%mPCSAFT(j))
                sum2 = sum2 + gl%J3X1_PCSAFTD(3,o,j,o,o)*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(j)*gl%molfractions(j)/(sigma_oj**2*gl%mPCSAFT(j))
                sum9 = sum9 + gl%J3X1_PCSAFTD(2,o,j,o,o)*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(j)*gl%molfractions(j)/(sigma_oj**2*gl%mPCSAFT(j))
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_jk * sigma_ok * sigma_oj * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum3 = sum3 + gl%J3X1_PCSAFTD(1,j,k,o,o) * sumb
                    sum4 = sum4 + gl%J3_PCSAFTD(1,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTD(3,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTD(3,j,k,o) * sumb
                    sum10 = sum10 + gl%J3X1_PCSAFTD(2,j,k,o,o) * sumb
                    sum11 = sum11 + gl%J3_PCSAFTD(2,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%nPCSAFTD(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij * sigma_ik * sigma_jk * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum5 = sum5 + gl%J3X1_PCSAFTD(3,i,j,k,o) * sumc
                            sum8 = sum8 + gl%J3X1_PCSAFTD(2,i,j,k,o) * sumc
                            sum12 = sum12 + gl%J3X1_PCSAFTD(1,i,j,k,o) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTD(o)**3 * gl%MyPCSAFTD(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**3) 
            part1 = gl%J3X1_PCSAFTD(1,o,o,o,o) * help
            part2 = gl%J3X1_PCSAFTD(3,o,o,o,o) * help
            part3 = gl%J3X1_PCSAFTD(2,o,o,o,o) * help
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o) * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTD(3,o) = - piPCSAFT**2 *(4.d0/3.d0)  * DENS * DENS * 1.d-57 * (2.d0 * part1 &
             & - 6.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
             & + 6.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum3 + Sum4) / help3 &
             & + part2 - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
             & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 &
             & + Sum5 + 4.d0 *(part3 - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum9 / help2 &
             & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum10 + Sum11) / help3 &
             & + Sum8) + 2.d0 * Sum12)
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if

    
    
   ! ! 4: 1ST DERIVATIVE OF a_3qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
   ! ! requires j_3ijk
    if (GETDERADD(4) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_oj**2 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTD(4,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTD(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_jk * sigma_ok * sigma_oj * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum3 = sum3 + gl%J3X1_PCSAFTD(4,j,k,o,o) * sumb
                    sum4 = sum4 + gl%J3_PCSAFTD(4,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTD(1,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTD(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%nPCSAFTD(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij * sigma_ik * sigma_jk * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum5 = sum5 + (-3.d0 * gl%J3X1_PCSAFTD(1,i,j,k,o) + gl%J3X1_PCSAFTD(4,i,j,k,o)) * sumc ! T für Jq_3_ijko_xo
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTD(o)**3 * gl%MyPCSAFTD(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**3) 
            part1 = - 3.d0 * gl%J3X1_PCSAFTD(1,o,o,o,o) * help ! T
            part2 = gl%J3X1_PCSAFTD(4,o,o,o,o) * help
            
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o) * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTD(4,o) =  -(4.d0/3.d0)  * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 + part2 + Sum5 &
                & - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
                & + 9.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 & !T
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum3 + Sum4) / help3 &
                & - 9.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 ) !T
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if

   
    
    
    ! 5: 2ND DERIVATIVE OF a_3qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    ! requires j_3ijk
    if (GETDERADD(5) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_oj**2 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTD(5,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTD(4,o,j,o,o) * suma
                sum3 = sum3 + gl%J3X1_PCSAFTD(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_jk * sigma_ok * sigma_oj * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum4 = sum4 + gl%J3X1_PCSAFTD(5,j,k,o,o) * sumb
                    sum5 = sum5 + gl%J3_PCSAFTD(5,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTD(4,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTD(4,j,k,o) * sumb
                    sum8 = sum8 + gl%J3X1_PCSAFTD(1,j,k,o,o) * sumb
                    sum9 = sum9 + gl%J3_PCSAFTD(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%nPCSAFTD(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij * sigma_ik * sigma_jk * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum10 = sum10 + (12.d0 * gl%J3X1_PCSAFTD(1,i,j,k,o) - 6.d0 * gl%J3X1_PCSAFTD(4,i,j,k,o) + gl%J3X1_PCSAFTD(5,i,j,k,o)) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTD(o)**3 * gl%MyPCSAFTD(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**3) 
            part1 = (12.d0 * gl%J3X1_PCSAFTD(1,o,o,o,o) - 6.d0 * gl%J3X1_PCSAFTD(4,o,o,o,o) + gl%J3X1_PCSAFTD(5,o,o,o,o)) * help
                        
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o) * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTD(5,o) =  - (4.d0/3.d0)  * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 + Sum10 &
                & - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
                & + 18.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & - 36.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum3 / help2 &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum4 + Sum5) / help3 &
                & - 18.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 &
                & + 36.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum8 + Sum9) / help3 )
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
  
    
    ! 6: 2ND MIXED DERIVATIVE OF a_3 WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_3ijk
    if (GETDERADD(6) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        sum13 = 0.d0
        sum14 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_oj**2 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTD(6,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTD(2,o,j,o,o) * suma
                sum3 = sum3 + gl%J3X1_PCSAFTD(4,o,j,o,o) * suma
                sum8 = sum8 + gl%J3X1_PCSAFTD(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_jk * sigma_ok * sigma_oj * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum4 = sum4 + gl%J3X1_PCSAFTD(6,j,k,o,o) * sumb
                    sum5 = sum5 + gl%J3_PCSAFTD(6,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTD(2,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTD(2,j,k,o) * sumb
                    sum9 = sum9 + gl%J3X1_PCSAFTD(4,j,k,o,o) * sumb
                    sum10 = sum10 + gl%J3_PCSAFTD(4,j,k,o) * sumb
                    sum11 = sum11 + gl%J3X1_PCSAFTD(1,j,k,o,o) * sumb
                    sum12 = sum12 + gl%J3_PCSAFTD(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%nPCSAFTD(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij * sigma_ik * sigma_jk * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum13 = sum13 + (gl%J3X1_PCSAFTD(6,i,j,k,o) - 3.d0 * gl%J3X1_PCSAFTD(2,i,j,k,o)) * sumc
                            sum14 = sum14 + (-3.d0 * gl%J3X1_PCSAFTD(1,i,j,k,o) + gl%J3X1_PCSAFTD(4,i,j,k,o)) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTD(o)**3 * gl%MyPCSAFTD(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**3) 
            part1 = gl%J3X1_PCSAFTD(6,o,o,o,o) * help
            part2 = 3.d0 * gl%J3X1_PCSAFTD(2,o,o,o,o) * help
            part3 = - 3.d0 * gl%J3X1_PCSAFTD(1,o,o,o,o) * help
            part4 = gl%J3X1_PCSAFTD(4,o,o,o,o) * help
                        
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o) * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTD(6,o) =  - (4.d0/3.d0)  * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 - part2  &
                & - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
                & + 9.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum4 + Sum5) / help3 &
                & - 9.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 &
                & + Sum13 + 2.d0 * (part3 + part4 &
                & - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum3 / help2 &
                & + 9.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum8 / help2 &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum9 + Sum10) / help3 &
                & - 9.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum11 + Sum12) / help3 &
                & + Sum14))
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
    
    ! 7: 3RD MIXED DERIVATIVE OF a_3qq WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires j_3ijk
    if (GETDERADD(7) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        sum13 = 0.d0
        sum14 = 0.d0
        sum15 = 0.d0
        sum16 = 0.d0
        sum17 = 0.d0
        sum18 = 0.d0
        sum19 = 0.d0
        sum20 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_oj**2 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTD(7,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTD(6,o,j,o,o) * suma
                sum3 = sum3 + gl%J3X1_PCSAFTD(2,o,j,o,o) * suma
                sum8 = sum8 + gl%J3X1_PCSAFTD(5,o,j,o,o) * suma
                sum15 = sum15 + gl%J3X1_PCSAFTD(4,o,j,o,o) * suma
                sum16 = sum16 + gl%J3X1_PCSAFTD(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_jk * sigma_ok * sigma_oj * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum4 = sum4 + gl%J3X1_PCSAFTD(7,j,k,o,o) * sumb
                    sum5 = sum5 + gl%J3_PCSAFTD(7,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTD(6,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTD(6,j,k,o) * sumb
                    sum9 = sum9 + gl%J3X1_PCSAFTD(2,j,k,o,o) * sumb
                    sum10 = sum10 + gl%J3_PCSAFTD(2,j,k,o) * sumb
                    sum11 = sum11 + gl%J3X1_PCSAFTD(5,j,k,o,o) * sumb
                    sum12 = sum12 + gl%J3_PCSAFTD(5,j,k,o) * sumb
                    sum17 = sum17 + gl%J3X1_PCSAFTD(4,j,k,o,o) * sumb
                    sum18 = sum18 + gl%J3_PCSAFTD(4,j,k,o) * sumb
                    sum19 = sum19 + gl%J3X1_PCSAFTD(1,j,k,o,o) * sumb
                    sum20 = sum20 + gl%J3_PCSAFTD(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%nPCSAFTD(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij * sigma_ik * sigma_jk * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum13 = sum13 + (gl%J3X1_PCSAFTD(7,i,j,k,o) - 6.d0 * gl%J3X1_PCSAFTD(6,i,j,k,o) + 12.d0 * gl%J3X1_PCSAFTD(2,i,j,k,o)) * sumc
                            sum14 = sum14 + (12.d0 * gl%J3X1_PCSAFTD(1,i,j,k,o) - 6.d0 * gl%J3X1_PCSAFTD(4,i,j,k,o) + gl%J3X1_PCSAFTD(5,i,j,k,o)) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTD(o)**3 * gl%MyPCSAFTD(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**3) 
            part1 = (gl%J3X1_PCSAFTD(7,o,o,o,o) - 6.d0 * gl%J3X1_PCSAFTD(6,o,o,o,o) + 12.d0 * gl%J3X1_PCSAFTD(2,o,o,o,o)) * help
            part2 = (12.d0 * gl%J3X1_PCSAFTD(1,o,o,o,o) - 6.d0 * gl%J3X1_PCSAFTD(4,o,o,o,o) + gl%J3X1_PCSAFTD(5,o,o,o,o)) * help
                        
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o) * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTD(7,o) =  -(4.d0/3.d0)  * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 &
                & - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * sum1 / help2 &
                & + 18.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & - 36.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum3 / help2 &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum4 + Sum5) / help3 &
                & - 18.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 &
                & + 36.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum9 + Sum10) / help3 &
                & + Sum13 + 2.d0 * (part2 &
                & - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum8 / help2 &
                & + 18.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum15 / help2 &
                & - 36.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum16 / help2 &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum11 + Sum12) / help3 &
                & - 18.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum17 + Sum18) / help3 &
                & + 36.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum19 + Sum20) / help3 &
                & + Sum14))
                ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if   
    
    
    ! 8: 3RD DERIVATIVE OF a_3qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires j_3ijk
    if (GETDERADD(8) .eq. 1) then
        do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_oj**2 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTD(2,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTD(8,o,j,o,o) * suma
                sum9 = sum9 + gl%J3X1_PCSAFTD(3,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) &
                        & / (sigma_jk * sigma_ok * sigma_oj * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum3 = sum3 + gl%J3X1_PCSAFTD(2,j,k,o,o) * sumb
                    sum4 = sum4 + gl%J3_PCSAFTD(2,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTD(8,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTD(8,j,k,o) * sumb
                    sum10 = sum10 + gl%J3X1_PCSAFTD(3,j,k,o,o) * sumb
                    sum11 = sum11 + gl%J3_PCSAFTD(3,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%nPCSAFTD(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij * sigma_ik * sigma_jk * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum5 = sum5 + gl%J3X1_PCSAFTD(8,i,j,k,o) * sumc
                            sum8 = sum8 + gl%J3X1_PCSAFTD(3,i,j,k,o) * sumc
                            sum12 = sum12 + gl%J3X1_PCSAFTD(2,i,j,k,o) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTD(o)**3 * gl%MyPCSAFTD(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**3) 
            part1 = 6.d0 * gl%J3X1_PCSAFTD(2,o,o,o,o) * help
            part2 = gl%J3X1_PCSAFTD(8,o,o,o,o) * help
            part3 = gl%J3X1_PCSAFTD(3,o,o,o,o) * help
            
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o) * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTD(8,o) = - piPCSAFT**2 * (4.d0/3.d0)  * DENS * DENS * 1.d-57 * (part1 &
                & - 18.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum1 / help2  &
                & + 18.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum3 + Sum4) / help3 &
                & + part2 - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3  &
                & + Sum5 + 6.d0 * (part3 - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum9 / help2 &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum10 + Sum11) / help3 &
                & + Sum8) + 6.d0 * Sum12)
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
  
     
   ! 9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires 
    if (GETDERADD(9) .eq. 1) then
      do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        sum13 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_oj**2 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTD(9,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTD(5,o,j,o,o) * suma
                sum3 = sum3 + gl%J3X1_PCSAFTD(4,o,j,o,o) * suma
                sum4 = sum4 + gl%J3X1_PCSAFTD(1,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_jk * sigma_ok * sigma_oj * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum5 = sum5 + gl%J3X1_PCSAFTD(9,j,k,o,o) * sumb
                    sum6 = sum6 + gl%J3_PCSAFTD(9,j,k,o) * sumb
                    sum7 = sum7 + gl%J3X1_PCSAFTD(5,j,k,o,o) * sumb
                    sum8 = sum8 + gl%J3_PCSAFTD(5,j,k,o) * sumb
                    sum9 = sum9 + gl%J3X1_PCSAFTD(4,j,k,o,o) * sumb
                    sum10 = sum10 + gl%J3_PCSAFTD(4,j,k,o) * sumb
                    sum11 = sum11 + gl%J3X1_PCSAFTD(1,j,k,o,o) * sumb
                    sum12 = sum12 + gl%J3_PCSAFTD(1,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%nPCSAFTD(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij * sigma_ik * sigma_jk * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                            sum13 = sum13 + (-60.d0 * gl%J3X1_PCSAFTD(1,i,j,k,o) + 36.d0 * gl%J3X1_PCSAFTD(4,i,j,k,o) - 9.d0 * gl%J3X1_PCSAFTD(5,i,j,k,o) + gl%J3X1_PCSAFTD(9,i,j,k,o)) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTD(o)**3 * gl%MyPCSAFTD(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**3) 
            part1 = (-60.d0 * gl%J3X1_PCSAFTD(1,o,o,o,o) + 36.d0 * gl%J3X1_PCSAFTD(4,o,o,o,o) - 9.d0 * gl%J3X1_PCSAFTD(5,o,o,o,o) + gl%J3X1_PCSAFTD(9,o,o,o,o)) * help
            
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o) * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTD(9,o) = -(4.d0/3.d0)  * piPCSAFT**2 * DENS**2 * 1.d-57 * (part1 &
                & - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 &
                & + 27.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & - 108.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum3 / help2 &
                & + 180.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum4 / help2 &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum5 + Sum6) / help3 &
                & - 27.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum7 + Sum8) / help3 &
                & + 108.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum9 + Sum10) / help3 &
                & - 180.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum11 + Sum12) / help3 &
                & + Sum13)
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if

    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires j_3ijk
    if (GETDERADD(10) .eq. 1) then
    do o = 1, gl%ncomp
        sum = 0.d0
        sum1 = 0.d0
        sum2 = 0.d0
        sum3 = 0.d0 
        sum4 = 0.d0
        sum5 = 0.d0
        sum6 = 0.d0
        sum7 = 0.d0
        sum8 = 0.d0
        sum9 = 0.d0
        sum10 = 0.d0
        sum11 = 0.d0
        sum12 = 0.d0
        sum13 = 0.d0
        sum14 = 0.d0
        sum15 = 0.d0
        sum16 = 0.d0
        sum17 = 0.d0
        sum18 = 0.d0
        sum19 = 0.d0
        sum20 = 0.d0
        sum21 = 0.d0
        do j = 1, gl%ncomp
                sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                suma = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_oj**2 * gl%mPCSAFT(j))
                sum1 = sum1 + gl%J3X1_PCSAFTD(4,o,j,o,o) * suma
                sum2 = sum2 + gl%J3X1_PCSAFTD(1,o,j,o,o) * suma
                sum3 = sum3 + gl%J3X1_PCSAFTD(10,o,j,o,o) * suma
                sum15 = sum15 + gl%J3X1_PCSAFTD(3,o,j,o,o) * suma
                sum16 = sum16 + gl%J3X1_PCSAFTD(6,o,j,o,o) * suma
                sum17 = sum17 + gl%J3X1_PCSAFTD(2,o,j,o,o) * suma
                do k = 1, gl%ncomp
                    sigma_oj = (gl%sigPCSAFT(o)+gl%sigPCSAFT(j))*0.5d0
                    sigma_ok = (gl%sigPCSAFT(o)+gl%sigPCSAFT(k))*0.5d0
                    sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                   
                    sumb = gl%molfractions(j) * gl%molfractions(k) * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) &
                        & / (sigma_jk * sigma_ok * sigma_oj * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                    sum4 = sum4 + gl%J3X1_PCSAFTD(4,j,k,o,o) * sumb
                    sum5 = sum5 + gl%J3_PCSAFTD(4,j,k,o) * sumb
                    sum6 = sum6 + gl%J3X1_PCSAFTD(1,j,k,o,o) * sumb
                    sum7 = sum7 + gl%J3_PCSAFTD(1,j,k,o) * sumb
                    sum8 = sum8 + gl%J3X1_PCSAFTD(10,j,k,o,o) * sumb
                    sum9 = sum9 + gl%J3_PCSAFTD(10,j,k,o) * sumb
                    sum10 = sum10 + gl%J3X1_PCSAFTD(3,j,k,o,o) * sumb
                    sum11 = sum11 + gl%J3_PCSAFTD(3,j,k,o) * sumb
                    sum12 = sum12 + gl%J3X1_PCSAFTD(6,j,k,o,o) * sumb
                    sum13 = sum13 + gl%J3_PCSAFTD(6,j,k,o) * sumb
                    sum18 = sum18 + gl%J3X1_PCSAFTD(2,j,k,o,o) * sumb
                    sum19 = sum19 + gl%J3_PCSAFTD(2,j,k,o) * sumb
                    do i = 1, gl%ncomp
                        if (i /= o .and. j /= o .and. k /= o)  then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_ik = (gl%sigPCSAFT(i)+gl%sigPCSAFT(k))*0.5d0
                            sigma_jk = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                            sumc = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%nPCSAFTD(k) * gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) &
                                & / (T**3 * kbol**3 * sigma_ij * sigma_ik * sigma_jk * gl%mPCSAFT(i) * gl%mPCSAFT(j) * gl%mPCSAFT(k))
                           
                            sum14 = sum14 + (gl%J3X1_PCSAFTD(10,i,j,k,o) - 3.d0 * gl%J3X1_PCSAFTD(3,i,j,k,o)) * sumc
                            sum20 = sum20 + (gl%J3X1_PCSAFTD(6,i,j,k,o) - 3.d0 * gl%J3X1_PCSAFTD(2,i,j,k,o)) * sumc
                            sum21 = sum21 + (-3.d0 * gl%J3X1_PCSAFTD(1,i,j,k,o) + gl%J3X1_PCSAFTD(4,i,j,k,o)) * sumc
                        end if
                    end do
                end do
            end do
            help =  gl%nPCSAFTD(o)**3 * gl%MyPCSAFTD(o)**6 * gl%molfractions(o)**3 / (kbol**3 * T**3 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**3) 
            part1 = (-6.d0 * gl%J3X1_PCSAFTD(1,o,o,o,o) + 2.d0 * gl%J3X1_PCSAFTD(4,o,o,o,o)) * help
            part2 = (gl%J3X1_PCSAFTD(10,o,o,o,o) - 3.d0 * gl%J3X1_PCSAFTD(3,o,o,o,o)) * help
            part3 = (gl%J3X1_PCSAFTD(6,o,o,o,o) - 3.d0 * gl%J3X1_PCSAFTD(2,o,o,o,o)) * help
            help2 = T**3 * kbol**3 * gl%sigPCSAFT(o) * gl%mPCSAFT(o)**2
            help3 = T**3 * kbol**3 * gl%mPCSAFT(o)
            
            gl%A3X1_PCSAFTD(10,o) = - piPCSAFT**2 * (4.d0/3.d0)  * DENS * DENS * 1.d-57 * (part1  &
                & - 6.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum1 / help2 & 
                & + 18.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum2 / help2 &
                & + 6.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum4 + Sum5) / help3 &
                & - 18.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum6 + Sum7) / help3 &
                & + part2 &
                & - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum3 / help2 &
                & + 9.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum15 / help2 &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum8 + Sum9) / help3 &
                & - 9.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum10 + Sum11) / help3 &
                & + Sum14 &
                & + 4.d0 * (part3 - 3.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum16 / help2 &
                & + 9.d0 * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 * gl%molfractions(o)**2 * Sum17 / help2 &
                & + 3.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum12 + Sum13) / help3 &
                & - 9.d0 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * (gl%molfractions(o) * Sum18 + Sum19) / help3 &
                & + Sum20) + 2.d0 * Sum21)
            ! 1.d-57: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do
    end if
!DEC$ END IF
end subroutine A3X1DERIVS_D
 
 
subroutine base_A3X2(gl,T,D,A3X2)
    ! a_3dd: second-order perturbation term of the dipolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.9 in Gross, Sadowski 2006:
    ! a_3dd =-4*pi**2/3*rho**2*sum_i*sum_j*sum_k(x_i*x_j*x_k*epsk_ii/T*espk_jj/T*espk_kk/T*sig_ii**3*sig_jj**3*sig_kk**3/sig_ij/sig_ik/sig_jk*n_Di*n_Dj*n_Dk*My_i**2*My_j**2*My_k**2*J3_D_ijk
    ! dependent on D and T and xi, xj and xk
   




   

implicit none

    type(type_gl) :: gl

   
    ! I. Declarations
    !inputs
    double precision :: T, D
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sigma_ij, sigma_jk, sigma_ik, sigma_ou, sigma_io, sigma_iu, sigma_jo, sigma_ju, sigma_ku
    double precision :: part1, sum1a, sum1b, sum1c, sum1d, sum2a, sum2b, sum3a, sum3b, sum3c, sum3d, sum4a, sum4b, sum4c, sum4d, sum4e, sum4f, sum4g, sum4h, sum, suma, sumb
    integer :: i, j, xi, k, o, u
    double precision :: DELTA, A_u, A_o, A_ou, A_oij, A_iou, A_i, A_j, A_k, A_io, A_ju, A_ku, A_iu 
    double precision, dimension(gl%ncomp,gl%ncomp) :: A3X2
    integer, dimension(nderivs) :: getprevious
    integer :: nrsubst
 
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    getprevious = 1
    
    !call init_derivs_PCSAFT(nrsubst)
    call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_QQ_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1_QQ(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x2(gl,T,D,getprevious)
    call J3X1DERIVS(gl,T,D,getprevious)
    call J3X2DERIVS(gl,T,D,getprevious)
    
   do u = 1, gl%ncomp
            A_u = gl%QPCSAFTQ(u)**2 * gl%nPCSAFTQ(u)  / (T * kbol * gl%sigPCSAFT(u)**3 * gl%mPCSAFT(u))
            do o = 1, gl%ncomp 
                part1 = 0.d0
                sum1a = 0.d0
                sum1b = 0.d0
                sum1c = 0.d0
                sum1d = 0.d0
                sum2a = 0.d0
                sum2b = 0.d0
                sum3a = 0.d0
                sum3b = 0.d0
                sum3c = 0.d0
                sum3d = 0.d0
                sum4a = 0.d0
                sum4b = 0.d0
                sum4c = 0.d0
                sum4d = 0.d0
                sum4e = 0.d0
                sum4f = 0.d0
                sum4g = 0.d0
                sum4h = 0.d0
                
                A_o = gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o)  / (T * kbol * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o))
                
                sigma_ou = (gl%sigPCSAFT(o)+gl%sigPCSAFT(u))*0.5d0
                A_ou = gl%QPCSAFTQ(u)**2 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(u) * gl%nPCSAFTQ(o)/ (T**2 * kbol**2 * sigma_ou**6 * gl%mPCSAFT(u) * gl%mPCSAFT(o))
                
                sum1a = 3.d0 * A_u*A_ou * gl%molfractions(u) * (2.d0 * gl%J3_PCSAFTQ(1,o,u,u) + gl%molfractions(u) * gl%J3X1_PCSAFTQ(1,o,u,u,u))
                
                if (u == o) then
                    part1 = gl%molfractions(o)**2 * (3.d0 * gl%J3X1_PCSAFTQ(1,o,o,o,o) + gl%molfractions(o) * gl%J3X2_PCSAFTQ(1,o,o,o,o,u)) * A_o**3  
                    
                else
                    part1 = gl%molfractions(o)**3 * gl%J3X2_PCSAFTQ(1,o,o,o,o,u) * A_o**3
                    
                    sum2a = sum2a - 3.d0 * A_o*A_ou * gl%molfractions(o)**2 * (gl%J3X1_PCSAFTQ(1,o,u,o,o) +  gl%molfractions(u) * gl%J3X2_PCSAFTQ(1,o,u,o,o,u))
                    sum3a = sum3a + 3.d0 * A_u*A_ou * gl%molfractions(o) * gl%molfractions(u) * (2.d0 *gl%J3X1_PCSAFTQ(1,o,u,u,o) + gl%molfractions(u) * gl%J3X2_PCSAFTQ(1,o,u,u,o,u))
                end if
                
                do i = 1, gl%ncomp
                    sigma_io = (gl%sigPCSAFT(i)+gl%sigPCSAFT(o))*0.5d0
                    sigma_iu = (gl%sigPCSAFT(i)+gl%sigPCSAFT(u))*0.5d0
                    sigma_ou = (gl%sigPCSAFT(u)+gl%sigPCSAFT(o))*0.5d0
                    A_i = gl%QPCSAFTQ(i)**2 * gl%nPCSAFTQ(i) / (T * kbol * gl%sigPCSAFT(i)**3 * gl%mPCSAFT(i))
                    A_io = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(o)/ (T**2 * kbol**2 * sigma_io**6 * gl%mPCSAFT(i) * gl%mPCSAFT(o))
                    A_iou = gl%QPCSAFTQ(o)**2 * gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(u)**2 * gl%nPCSAFTQ(o) * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(u) / (T**3 * kbol**3 * sigma_io**3 * sigma_iu**3 * sigma_ou**3 * gl%mPCSAFT(o) * gl%mPCSAFT(u) * gl%mPCSAFT(i))
                    
                    if (i /= u) then 
                    
                        sum1b = sum1b + 3.d0 * gl%molfractions(i) * (gl%J3_PCSAFTQ(1,i,u,o) + gl%molfractions(u) * gl%J3X1_PCSAFTQ(1,i,u,o,u)) * A_iou
                        
                        if (u == o) then
                            sum2b = sum2b - 3.d0 * gl%molfractions(i) * gl%molfractions(o) * (2.d0 * gl%J3X1_PCSAFTQ(1,i,o,o,o) + gl%molfractions(o) * gl%J3X2_PCSAFTQ(1,i,o,o,o,u)) * A_io*A_o
                            sum3b = sum3b + 3.d0 * gl%molfractions(i) * gl%molfractions(o) * (2.d0 * gl%J3X1_PCSAFTQ(1,i,u,o,o) + gl%molfractions(u) * gl%J3X2_PCSAFTQ(1,i,u,o,o,u)) * A_iou
                        else
                            sum2b = sum2b - 3.d0 * gl%molfractions(i) * gl%molfractions(o)**2 * gl%J3X2_PCSAFTQ(1,i,o,o,o,u) * A_io*A_o
                            sum3b = sum3b + 3.d0 * gl%molfractions(i) * gl%molfractions(o) * (gl%J3X1_PCSAFTQ(1,i,u,o,o) + gl%molfractions(u) * gl%J3X2_PCSAFTQ(1,i,u,o,o,u)) * A_iou      !entspricht auch sum3c, unten mal 2 
                        end if
                        
                    end if
                    
                    do j = 1, gl%ncomp
                        sigma_ju = (gl%sigPCSAFT(j)+gl%sigPCSAFT(u))*0.5d0
                        A_ju =  gl%QPCSAFTQ(u)**2 * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(u) * gl%nPCSAFTQ(j) / (T**2 * kbol**2 * sigma_ju**6 * gl%mPCSAFT(j) * gl%mPCSAFT(u))
                        A_j = gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (T * kbol * gl%sigPCSAFT(j)**3 * gl%mPCSAFT(j))
                        if (j /= u .and. i /= u) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_io = (gl%sigPCSAFT(i)+gl%sigPCSAFT(o))*0.5d0
                            sigma_jo = (gl%sigPCSAFT(j)+gl%sigPCSAFT(o))*0.5d0
                            
                            A_oij = gl%QPCSAFTQ(o)**2 * gl%QPCSAFTQ(j)**2 * gl%QPCSAFTQ(i)**2 * gl%nPCSAFTQ(o) * gl%nPCSAFTQ(j) * gl%nPCSAFTQ(i) / (T**3 * kbol**3 * sigma_ij**3 * sigma_io**3 * sigma_jo**3 * gl%mPCSAFT(o) * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            sum1d = sum1d + 3.d0 * gl%molfractions(j)* gl%molfractions(i) * gl%J3X1_PCSAFTQ(1,i,j,o,u) * A_oij   
                                
                            if (u == o) then
                                sum3d = sum3d + 3.d0 * gl%molfractions(j)* gl%molfractions(i)* (gl%J3X1_PCSAFTQ(1,i,j,o,o) + gl%molfractions(o) * gl%J3X2_PCSAFTQ(1,i,j,o,o,u)) * A_oij                            
                            else
                                sum3d = sum3d + 3.d0 * gl%molfractions(j)* gl%molfractions(i)* gl%molfractions(o) * gl%J3X2_PCSAFTQ(1,i,j,o,o,u) * A_oij                                 
                            end if
                        end if
                                   !Bis hier programmiert
                        do k = 1, gl%ncomp
                            sigma_ku = (gl%sigPCSAFT(k)+gl%sigPCSAFT(u))*0.5d0
                            A_k = gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) / (T * kbol * gl%sigPCSAFT(k)**3 * gl%mPCSAFT(k))
                            A_ku = gl%QPCSAFTQ(u)**2 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(u) * gl%nPCSAFTQ(k) / (T**2 * kbol**2 * sigma_ku**6 * gl%mPCSAFT(k) * gl%mPCSAFT(u))
                            
                            if (i /= o .and. j /= o .and. k /= o) then
                                if (i == u .and. j == u .and. k == u) then 
                                    sum4a = sum4a + gl%molfractions(u)**2 * (3.d0 * gl%J3X1_PCSAFTQ(1,u,u,u,o) + gl%molfractions(u) * gl%J3X2_PCSAFTQ(1,u,u,u,o,u)) * A_u**3
                                    
                                elseif (i /= u .and. j == u .and. k == u) then 
                                    sum4b = sum4b + gl%molfractions(u) * (2.d0 * gl%molfractions(i) * gl%J3X1_PCSAFTQ(1,i,u,u,o) + gl%molfractions(u) * gl%molfractions(i) * gl%J3X2_PCSAFTQ(1,i,u,u,o,u)) * A_iu*A_u
                                    
                                elseif (i == u .and. j /= u .and. k == u) then 
                                    sum4c = sum4c + gl%molfractions(u) * (2.d0 * gl%molfractions(j) * gl%J3X1_PCSAFTQ(1,u,j,u,o) + gl%molfractions(u) * gl%molfractions(j) * gl%J3X2_PCSAFTQ(1,u,j,u,o,u)) * A_ju*A_u
                                    
                                elseif (i == u .and. j == u .and. k /= u) then 
                                    sum4d = sum4d + gl%molfractions(u) * (2.d0 * gl%molfractions(k) * gl%J3X1_PCSAFTQ(1,u,u,k,o) + gl%molfractions(u) * gl%molfractions(k) * gl%J3X2_PCSAFTQ(1,u,u,k,o,u)) * A_ku*A_u
                                    
                                elseif (i == u .and. j /= u .and. k /= u) then 
                                    sum4e = sum4e + gl%molfractions(j) * gl%molfractions(k) * (gl%J3X1_PCSAFTQ(1,u,j,k,o) + gl%molfractions(u) * gl%J3X2_PCSAFTQ(1,u,j,k,o,u)) * A_ku*A_j
                                    
                                elseif (i /= u .and. j == u .and. k /= u) then 
                                    sum4f = sum4f + gl%molfractions(i) * gl%molfractions(k) * (gl%J3X1_PCSAFTQ(1,i,u,k,o) + gl%molfractions(u) * gl%J3X2_PCSAFTQ(1,i,u,k,o,u)) * A_ku*A_i
                                    
                                elseif (i /= u .and. j /= u .and. k == u) then 
                                    sum4g = sum4g + gl%molfractions(j) * gl%molfractions(i) * (gl%J3X1_PCSAFTQ(1,i,j,u,o) + gl%molfractions(u) * gl%J3X2_PCSAFTQ(1,i,j,u,o,u)) * A_ju*A_i
                                    
                                elseif (i /= u .and. j /= u .and. k /= u) then 
                                    sum4h = sum4h + gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) * gl%J3X2_PCSAFTQ(1,i,j,k,o,u) * A_i*A_j*A_k
                                    
                                end if
                            end if
                        end do
                    end do
                end do              
                A3X2(o,u) =  piPCSAFT**2 * 0.5625d0 * D * D * 1.d-57 * (part1 + sum1a + 2.d0*sum1b + sum1c + sum1d + sum2a + sum2b + sum3a + 2.d0*sum3b + sum3c + sum3d + sum4a + sum4b +sum4c + sum4d + sum4e + sum4f + sum4g + sum4h)
            end do 
        end do 

!DEC$ END IF
end subroutine

 
subroutine base_A3X2_D(gl,T,D,A3X2_D)
    ! a_3dd: second-order perturbation term of the dipolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.9 in Gross, Sadowski 2006:
    ! a_3dd =-4*pi**2/3*rho**2*sum_i*sum_j*sum_k(x_i*x_j*x_k*epsk_ii/T*espk_jj/T*espk_kk/T*sig_ii**3*sig_jj**3*sig_kk**3/sig_ij/sig_ik/sig_jk*n_Di*n_Dj*n_Dk*My_i**2*My_j**2*My_k**2*J3_D_ijk
    ! dependent on D and T and xi, xj and xk
   




   

implicit none

    type(type_gl) :: gl

   
      ! I. Declarations
    !inputs
    double precision :: T, D
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sigma_ij, sigma_jk, sigma_ik, sigma_ou, sigma_io, sigma_iu, sigma_jo, sigma_ju, sigma_ku
    double precision :: part1, sum1a, sum1b, sum1c, sum1d, sum2a, sum2b, sum3a, sum3b, sum3c, sum3d, sum4a, sum4b, sum4c, sum4d, sum4e, sum4f, sum4g, sum4h, sum, suma, sumb
    integer :: i, j, xi, k, o, u
    double precision :: DELTA, A_u, A_o, A_ou, A_oij, A_iou, A_i, A_j, A_k, A_io, A_ju, A_ku, A_iu 
    double precision, dimension(gl%ncomp,gl%ncomp) :: A3X2_D
    integer, dimension(nderivs) :: getprevious
    integer :: nrsubst

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    getprevious = 1
    
    !call init_derivs_PCSAFT(nrsubst)
    call calculate_PCSAFT_functionparts_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_DD_Trho(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x1_DD(gl,T,D,getprevious)
    call calculate_PCSAFT_functionparts_x2(gl,T,D,getprevious)
    call J3X1DERIVS_D(gl,T,D,getprevious)
    call J3X2DERIVS_D(gl,T,D,getprevious)
    
   do u = 1, gl%ncomp
            A_u = gl%MyPCSAFTD(u)**2 * gl%nPCSAFTD(u)  / (T * kbol * gl%sigPCSAFT(u) * gl%mPCSAFT(u))
            do o = 1, gl%ncomp 
                part1 = 0.d0
                sum1a = 0.d0
                sum1b = 0.d0
                sum1c = 0.d0
                sum1d = 0.d0
                sum2a = 0.d0
                sum2b = 0.d0
                sum3a = 0.d0
                sum3b = 0.d0
                sum3c = 0.d0
                sum3d = 0.d0
                sum4a = 0.d0
                sum4b = 0.d0
                sum4c = 0.d0
                sum4d = 0.d0
                sum4e = 0.d0
                sum4f = 0.d0
                sum4g = 0.d0
                sum4h = 0.d0
                
                A_o = gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o)  / (T * kbol * gl%sigPCSAFT(o) * gl%mPCSAFT(o))
                
                sigma_ou = (gl%sigPCSAFT(o)+gl%sigPCSAFT(u))*0.5d0
                A_ou = gl%MyPCSAFTD(u)**2 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(u) * gl%nPCSAFTD(o)/ (T**2 * kbol**2 * sigma_ou**2 * gl%mPCSAFT(u) * gl%mPCSAFT(o))
                
                sum1a = 3.d0 * A_u*A_ou * gl%molfractions(u) * (2.d0 * gl%J3_PCSAFTD(1,o,u,u) + gl%molfractions(u) * gl%J3X1_PCSAFTD(1,o,u,u,u))
                
                if (u == o) then
                    part1 = gl%molfractions(o)**2 * (3.d0 * gl%J3X1_PCSAFTD(1,o,o,o,o) + gl%molfractions(o) * gl%J3X2_PCSAFTD(1,o,o,o,o,u)) * A_o**3  
                    
                else
                    part1 = gl%molfractions(o)**3 * gl%J3X2_PCSAFTD(1,o,o,o,o,u) * A_o**3
                    
                    sum2a = sum2a - 3.d0 * A_o*A_ou * gl%molfractions(o)**2 * (gl%J3X1_PCSAFTD(1,o,u,o,o) +  gl%molfractions(u) * gl%J3X2_PCSAFTD(1,o,u,o,o,u))
                    sum3a = sum3a + 3.d0 * A_u*A_ou * gl%molfractions(o) * gl%molfractions(u) * (2.d0 *gl%J3X1_PCSAFTD(1,o,u,u,o) + gl%molfractions(u) * gl%J3X2_PCSAFTD(1,o,u,u,o,u))
                end if
                
                do i = 1, gl%ncomp
                    sigma_io = (gl%sigPCSAFT(i)+gl%sigPCSAFT(o))*0.5d0
                    sigma_iu = (gl%sigPCSAFT(i)+gl%sigPCSAFT(u))*0.5d0
                    sigma_ou = (gl%sigPCSAFT(u)+gl%sigPCSAFT(o))*0.5d0
                    A_i = gl%MyPCSAFTD(i)**2 * gl%nPCSAFTD(i) / (T * kbol * gl%sigPCSAFT(i) * gl%mPCSAFT(i))
                    A_io = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(o)/ (T**2 * kbol**2 * sigma_io**2 * gl%mPCSAFT(i) * gl%mPCSAFT(o))
                    A_iou = gl%MyPCSAFTD(o)**2 * gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(u)**2 * gl%nPCSAFTD(o) * gl%nPCSAFTD(i) * gl%nPCSAFTD(u) / (T**3 * kbol**3 * sigma_io * sigma_iu * sigma_ou * gl%mPCSAFT(o) * gl%mPCSAFT(u) * gl%mPCSAFT(i))
                    
                    if (i /= u) then ! i == u schon darüber erledigt
                    
                        sum1b = sum1b + 3.d0 * gl%molfractions(i) * (gl%J3_PCSAFTD(1,i,u,o) + gl%molfractions(u) * gl%J3X1_PCSAFTD(1,i,u,o,u)) * A_iou
                        
                        if (u == o) then
                            sum2b = sum2b - 3.d0 * gl%molfractions(i) * gl%molfractions(o) * (2.d0 * gl%J3X1_PCSAFTD(1,i,o,o,o) + gl%molfractions(o) * gl%J3X2_PCSAFTD(1,i,o,o,o,u)) * A_io*A_o
                            sum3b = sum3b + 3.d0 * gl%molfractions(i) * gl%molfractions(o) * (2.d0 * gl%J3X1_PCSAFTD(1,i,u,o,o) + gl%molfractions(u) * gl%J3X2_PCSAFTD(1,i,u,o,o,u)) * A_iou
                        else
                            sum2b = sum2b - 3.d0 * gl%molfractions(i) * gl%molfractions(o)**2 * gl%J3X2_PCSAFTD(1,i,o,o,o,u) * A_io*A_o
                            sum3b = sum3b + 3.d0 * gl%molfractions(i) * gl%molfractions(o) * (gl%J3X1_PCSAFTD(1,i,u,o,o) + gl%molfractions(u) * gl%J3X2_PCSAFTD(1,i,u,o,o,u)) * A_iou      !entspricht auch sum3c, unten mal 2 
                        end if
                        
                    end if
                    
                    do j = 1, gl%ncomp
                        sigma_ju = (gl%sigPCSAFT(j)+gl%sigPCSAFT(u))*0.5d0
                        A_ju =  gl%MyPCSAFTD(u)**2 * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(u) * gl%nPCSAFTD(j) / (T**2 * kbol**2 * sigma_ju**2 * gl%mPCSAFT(j) * gl%mPCSAFT(u))
                        A_j = gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (T * kbol * gl%sigPCSAFT(j) * gl%mPCSAFT(j))
                        if (j /= u .and. i /= u) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sigma_io = (gl%sigPCSAFT(i)+gl%sigPCSAFT(o))*0.5d0
                            sigma_jo = (gl%sigPCSAFT(j)+gl%sigPCSAFT(o))*0.5d0
                            
                            A_oij = gl%MyPCSAFTD(o)**2 * gl%MyPCSAFTD(j)**2 * gl%MyPCSAFTD(i)**2 * gl%nPCSAFTD(o) * gl%nPCSAFTD(j) * gl%nPCSAFTD(i) / (T**3 * kbol**3 * sigma_ij * sigma_io * sigma_jo * gl%mPCSAFT(o) * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            sum1d = sum1d + 3.d0 * gl%molfractions(j)* gl%molfractions(i) * gl%J3X1_PCSAFTD(1,i,j,o,u) * A_oij   
                                
                            if (u == o) then
                                sum3d = sum3d + 3.d0 * gl%molfractions(j)* gl%molfractions(i)* (gl%J3X1_PCSAFTD(1,i,j,o,o) + gl%molfractions(o) * gl%J3X2_PCSAFTD(1,i,j,o,o,u)) * A_oij                            
                            else
                                sum3d = sum3d + 3.d0 * gl%molfractions(j)* gl%molfractions(i)* gl%molfractions(o) * gl%J3X2_PCSAFTD(1,i,j,o,o,u) * A_oij                                 
                            end if
                        end if
                                   !Bis hier programmiert
                        do k = 1, gl%ncomp
                            sigma_ku = (gl%sigPCSAFT(k)+gl%sigPCSAFT(u))*0.5d0
                            A_k = gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) / (T * kbol * gl%sigPCSAFT(k) * gl%mPCSAFT(k))
                            A_ku = gl%MyPCSAFTD(u)**2 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(u) * gl%nPCSAFTD(k) / (T**2 * kbol**2 * sigma_ku**2 * gl%mPCSAFT(k) * gl%mPCSAFT(u))
                            
                            if (i /= o .and. j /= o .and. k /= o) then
                                if (i == u .and. j == u .and. k == u) then 
                                    sum4a = sum4a + gl%molfractions(u)**2 * (3.d0 * gl%J3X1_PCSAFTD(1,u,u,u,o) + gl%molfractions(u) * gl%J3X2_PCSAFTD(1,u,u,u,o,u)) * A_u**3
                                    
                                elseif (i /= u .and. j == u .and. k == u) then 
                                    sum4b = sum4b + gl%molfractions(u) * (2.d0 * gl%molfractions(i) * gl%J3X1_PCSAFTD(1,i,u,u,o) + gl%molfractions(u) * gl%molfractions(i) * gl%J3X2_PCSAFTD(1,i,u,u,o,u)) * A_iu*A_u
                                    
                                elseif (i == u .and. j /= u .and. k == u) then 
                                    sum4c = sum4c + gl%molfractions(u) * (2.d0 * gl%molfractions(j) * gl%J3X1_PCSAFTD(1,u,j,u,o) + gl%molfractions(u) * gl%molfractions(j) * gl%J3X2_PCSAFTD(1,u,j,u,o,u)) * A_ju*A_u
                                    
                                elseif (i == u .and. j == u .and. k /= u) then 
                                    sum4d = sum4d + gl%molfractions(u) * (2.d0 * gl%molfractions(k) * gl%J3X1_PCSAFTD(1,u,u,k,o) + gl%molfractions(u) * gl%molfractions(k) * gl%J3X2_PCSAFTD(1,u,u,k,o,u)) * A_ku*A_u
                                    
                                elseif (i == u .and. j /= u .and. k /= u) then 
                                    sum4e = sum4e + gl%molfractions(j) * gl%molfractions(k) * (gl%J3X1_PCSAFTD(1,u,j,k,o) + gl%molfractions(u) * gl%J3X2_PCSAFTD(1,u,j,k,o,u)) * A_ku*A_j
                                    
                                elseif (i /= u .and. j == u .and. k /= u) then 
                                    sum4f = sum4f + gl%molfractions(i) * gl%molfractions(k) * (gl%J3X1_PCSAFTD(1,i,u,k,o) + gl%molfractions(u) * gl%J3X2_PCSAFTD(1,i,u,k,o,u)) * A_ku*A_i
                                    
                                elseif (i /= u .and. j /= u .and. k == u) then 
                                    sum4g = sum4g + gl%molfractions(j) * gl%molfractions(i) * (gl%J3X1_PCSAFTD(1,i,j,u,o) + gl%molfractions(u) * gl%J3X2_PCSAFTD(1,i,j,u,o,u)) * A_ju*A_i
                                    
                                elseif (i /= u .and. j /= u .and. k /= u) then 
                                    sum4h = sum4h + gl%molfractions(i) * gl%molfractions(j) * gl%molfractions(k) * gl%J3X2_PCSAFTD(1,i,j,k,o,u) * A_i*A_j*A_k
                                    
                                end if
                            end if
                        end do
                    end do
                end do              
                A3X2_D(o,u) =  - piPCSAFT**2 * (4.d0/3.d0) * D * D * 1.d-57 * (part1 + sum1a + 2.d0*sum1b + sum1c + sum1d + sum2a + sum2b + sum3a + 2.d0*sum3b + sum3c + sum3d + sum4a + sum4b +sum4c + sum4d + sum4e + sum4f + sum4g + sum4h)
            end do 
        end do 

!DEC$ END IF
end subroutine


    end module pc_saft_A3X_derivs_module
