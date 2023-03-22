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

    ! module for file pc-saft_A2X_derivs.f90
    module pc_saft_A2X_derivs_module
    !global use inclusion
    use module_all_types
    use pc_saft_module
    use pc_saft_ancillary_routines_module
    use variables_transformation_module
    
    use pc_saft_JX_derivs_module

    contains



subroutine A2X1DERIVS(gl,T,DENS,GETDERAQQ)
    ! a_2qq: second-order perturbation term of the quadrupolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.9 in Gross, Sadowski 2005:
    ! a_2qq =-pi*(3/4)**2*rho*sum_j*sum_i(x_i*x_j*epsk_ii/T*espk_jj/T*sig_ii**5*sig_jj**5/sig_ij**7*n_Qi*n_Qj*Q_i**2*Q_j**2*J_2ij
    ! dependent on D and T and x_i







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAQQ
    double precision, intent (in) :: T, DENS
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sigma_ij, sigma_kj
    double precision :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, sum14, part13, part3, part3a, part3b, part3c, sum, help,help2, help3, part1, part2
    integer :: i, j, xi,  k
    double precision, dimension(2) :: qq
    
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !calculate the derivatives of a_qq
    ! 1: a_2qq
    if (GETDERAQQ(1) .eq. 1) then
        do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
                do j = 1, gl%ncomp
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (sigma_kj**7 * gl%mPCSAFT(j))
                    sum1 = sum1 + sum*gl%J2X1_PCSAFTQ(1,j,k,k) !*J2X1_PCSAFTQ(1,j,k,k) sollte egal sein
                    sum2 = sum2 + sum*gl%J2_PCSAFTQ(1,j,k)
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sum3 = sum3 + gl%molfractions(i)*gl%molfractions(j) * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) * gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(i)*gl%mPCSAFT(j)) * gl%J2X1_PCSAFTQ(1,i,j,k)
                        end if
                    end do
                end do
            part3 = gl%molfractions(k)**2 * gl%nPCSAFTQ(k)**2 * gl%QPCSAFTQ(k)**4 / (kbol**2 * T**2 * gl%sigPCSAFT(k)**7 * gl%mPCSAFT(k)**2) * gl%J2X1_PCSAFTQ(1,k,k,k)
            gl%A2X1_PCSAFTQ(1,k) = - piPCSAFT*0.5625d0*DENS*1.d-38  * (2.d0 * gl%QPCSAFTQ(k)**2 *gl%nPCSAFTQ(k) / (T**2*kbol**2*gl%mPCSAFT(k)) * (gl%molfractions(k)*sum1 + sum2) - part3 + sum3)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
    
    !  2: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires j_2qq
    if (GETDERAQQ(2) .eq. 1) then   
        do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
            sum4 = 0.d0
            sum5 = 0.d0
            sum6 = 0.d0
                do j = 1, gl%ncomp
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (T*kbol*sigma_kj**7 * gl%mPCSAFT(j))
                    sum1 = sum1 + gl%J2X1_PCSAFTQ(2,j,k,k) * sum
                    sum2 = sum2 + gl%J2_PCSAFTQ(2,j,k) * sum
                    sum3 = sum3 + gl%J2X1_PCSAFTQ(1,j,k,k) * sum
                    sum4 = sum4 + gl%J2_PCSAFTQ(1,j,k) * sum
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            part1 = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            sum6 = sum6 + gl%J2X1_PCSAFTQ(2,i,j,k) * part1
                            sum5 = sum5 + gl%J2X1_PCSAFTQ(1,i,j,k) * part1
                        end if
                    end do
                end do
            part3a = - gl%J2X1_PCSAFTQ(2,k,k,k) * gl%QPCSAFTQ(k)**4 * gl%nPCSAFTQ(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**7 * gl%mPCSAFT(k)**2) 
            part3b = - gl%J2X1_PCSAFTQ(1,k,k,k) * gl%QPCSAFTQ(k)**4 * gl%nPCSAFTQ(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**7 * gl%mPCSAFT(k)**2)
            part2 = 2.d0 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTQ(2,k) = - 0.5625d0 * piPCSAFT * DENS * 1.d-38 * (part3a + part2 * (gl%molfractions(k) * Sum1 + Sum2) + Sum6 + part3b + part2 * (gl%molfractions(k) * Sum3 + Sum4)  + Sum5)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
   


     
   ! 3: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
   ! requires j_2qq
    if (GETDERAQQ(3) .eq. 1) then
        do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
            sum4 = 0.d0
            sum5 = 0.d0
            sum6 = 0.d0
                do j = 1, gl%ncomp
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (T*kbol*sigma_kj**7 * gl%mPCSAFT(j))
                    sum1 = sum1 + gl%J2X1_PCSAFTQ(2,j,k,k)*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(j)*gl%molfractions(j)/(T*kbol*sigma_kj**7*gl%mPCSAFT(j))
                    sum2 = sum2 + gl%J2_PCSAFTQ(2,j,k)*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(j)*gl%molfractions(j)/(T*kbol*sigma_kj**7*gl%mPCSAFT(j))
                    sum3 = sum3 + gl%J2X1_PCSAFTQ(3,j,k,k)*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(j)*gl%molfractions(j)/(T*kbol*sigma_kj**7*gl%mPCSAFT(j))
                    sum4 = sum4 + gl%J2_PCSAFTQ(3,j,k)*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(j)*gl%molfractions(j)/(T*kbol*sigma_kj**7*gl%mPCSAFT(j))
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            part1 = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            sum6 = sum6 + gl%J2X1_PCSAFTQ(2,i,j,k)*gl%QPCSAFTQ(i)**2*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%molfractions(i)*gl%molfractions(j)/(T**2*kbol**2*sigma_ij**7*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                            sum5 = sum5 + gl%J2X1_PCSAFTQ(3,i,j,k)*gl%QPCSAFTQ(i)**2*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%molfractions(i)*gl%molfractions(j)/(T**2*kbol**2*sigma_ij**7*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                        end if
                    end do
                end do
            part3a = - gl%J2X1_PCSAFTQ(2,k,k,k)*gl%QPCSAFTQ(k)**4*gl%nPCSAFTQ(k)**2*gl%molfractions(k)**2/(T**2*kbol**2*gl%sigPCSAFT(k)**7*gl%mPCSAFT(k)**2)
            part3b = - gl%J2X1_PCSAFTQ(3,k,k,k) * gl%QPCSAFTQ(k)**4 * gl%nPCSAFTQ(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**7 * gl%mPCSAFT(k)**2)
            part2 = 2.d0 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTQ(3,k) = - piPCSAFT * 0.5625d0 * DENS * 1.d-38 *(2.d0 * part3a + 2.d0 * part2 * (gl%molfractions(k) * Sum1 + Sum2) + (part3b + part2 * (gl%molfractions(k) * Sum3 + Sum4) + Sum5) + 2.d0 * Sum6)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
   
     
    ! 4: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires j_2qq
    if (GETDERAQQ(4) .eq. 1) then
      do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
            sum4 = 0.d0
            sum5 = 0.d0
            
                do j = 1, gl%ncomp
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (T*kbol*sigma_kj**7 * gl%mPCSAFT(j))
                    sum1 = sum1 + (gl%J2X1_PCSAFTQ(4,j,k,k)- gl%J2X1_PCSAFTQ(1,j,k,k)) * sum
                    sum2 = sum2 + (gl%J2_PCSAFTQ(4,j,k) - gl%J2_PCSAFTQ(1,j,k)) * sum 
                    sum3 = sum3 + gl%J2X1_PCSAFTQ(1,j,k,k) * sum
                    sum4 = sum4 + gl%J2_PCSAFTQ(1,j,k) * sum
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sum5 = sum5 + (gl%J2X1_PCSAFTQ(4,i,j,k) - 2.d0 * gl%J2X1_PCSAFTQ(1,i,j,k)) * gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%molfractions(i) * gl%molfractions(j) &
                                &  / (T**2 * kbol**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                        end if
                    end do
                end do
            part3 = (2.d0 * gl%J2X1_PCSAFTQ(1,k,k,k) - gl%J2X1_PCSAFTQ(4,k,k,k)) * gl%QPCSAFTQ(k)**4 * gl%nPCSAFTQ(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**7 * gl%mPCSAFT(k)**2)
            part2 = 2.d0 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTQ(4,k) = -0.5625d0 * piPCSAFT * DENS * 1.d-38 * (part3 + part2 * (gl%molfractions(k) * Sum1 + Sum2)   &
                & - part2 * (gl%molfractions(k) * Sum3 + Sum4) + Sum5)
             !1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    

    ! 5: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    !requires j_2qq
    if (GETDERAQQ(5) .eq. 1) then
        do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
            sum4 = 0.d0
            sum5 = 0.d0
            sum6 = 0.d0
            sum7 = 0.d0
            
            do j = 1, gl%ncomp
                sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                sum = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (T*kbol*sigma_kj**7 * gl%mPCSAFT(j))
                sum1 = sum1 + (2.d0 * gl%J2X1_PCSAFTQ(1,j,k,k) - 2.d0 * gl%J2X1_PCSAFTQ(4,j,k,k) + gl%J2X1_PCSAFTQ(5,j,k,k)) * sum
                sum2 = sum2 + (2.d0 * gl%J2_PCSAFTQ(1,j,k) - 2.d0 * gl%J2_PCSAFTQ(4,j,k) + gl%J2_PCSAFTQ(5,j,k)) * sum
                sum3 = sum3 + (- gl%J2X1_PCSAFTQ(1,j,k,k) + gl%J2X1_PCSAFTQ(4,j,k,k)) * sum
                sum4 = sum4 + (- gl%J2_PCSAFTQ(1,j,k) + gl%J2_PCSAFTQ(4,j,k)) * sum 
                sum5 = sum5 + gl%J2X1_PCSAFTQ(1,j,k,k) * sum
                sum6 = sum6 + gl%J2_PCSAFTQ(1,j,k) * sum
                do i = 1, gl%ncomp
                    if (i /= k .and. j /= k) then 
                        sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                        sum7 = sum7 + gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * (6.d0 * gl%J2X1_PCSAFTQ(1,i,j,k) - 4.d0 * gl%J2X1_PCSAFTQ(4,i,j,k) + gl%J2X1_PCSAFTQ(5,i,j,k)) &
                            & * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                    end if
                end do
            end do
            part3 = (- 6.d0 * gl%J2X1_PCSAFTQ(1,k,k,k) + 4.d0 * gl%J2X1_PCSAFTQ(4,k,k,k) - gl%J2X1_PCSAFTQ(5,k,k,k)) / (T**2 * kbol**2 * gl%sigPCSAFT(k)**7 * gl%mPCSAFT(k)**2) &
            & * gl%QPCSAFTQ(k)**4 * gl%nPCSAFTQ(k)**2 * gl%molfractions(k)**2
            part2 = 2.d0 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTQ(5,k) =  -0.5625d0 * piPCSAFT * DENS * 1.d-38 * (part3 + part2 * (gl%molfractions(k) * Sum1 + Sum2) &
                & - 2.d0 * part2 * (gl%molfractions(k) * Sum3 + Sum4) &
                & + 2.d0 * part2 * (gl%molfractions(k) * Sum5 + Sum6) + Sum7)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
       
     
    ! 6: 2ND MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_2qq
    if (GETDERAQQ(6) .eq. 1) then
         do k = 1, gl%ncomp
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
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (T*kbol*sigma_kj**7 * gl%mPCSAFT(j))
                    sum1 = sum1 + (gl%J2X1_PCSAFTQ(6,j,k,k) - gl%J2X1_PCSAFTQ(2,j,k,k)) * sum
                        !Jq_2_kjk_xk_D bereits mit T multipliziert

                    sum2 = sum2 + (gl%J2_PCSAFTQ(6,j,k) - gl%J2_PCSAFTQ(2,j,k)) * sum
                        !Jq_2_kj_D bereits mit T multipliziert
                    
                    sum3 = sum3 + gl%J2X1_PCSAFTQ(2,j,k,k) * sum
                    sum4 = sum4 + gl%J2_PCSAFTQ(2,j,k) * sum
                    sum5 = sum5 + (-gl%J2X1_PCSAFTQ(1,j,k,k) + gl%J2X1_PCSAFTQ(4,j,k,k)) * sum 
                        !Jq_2_kjk_xk bereits mit T multipliziert
                    
                    sum6 = sum6 + (-gl%J2_PCSAFTQ(1,j,k) +  gl%J2_PCSAFTQ(4,j,k)) * sum
                        !Jq_2_kj bereits mit T multipliziert
                    
                    sum7 = sum7 + gl%J2X1_PCSAFTQ(1,j,k,k) * sum
                    sum8 = sum8 + gl%J2_PCSAFTQ(1,j,k) * sum
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sum9 = sum9 + (gl%J2X1_PCSAFTQ(6,i,j,k) - 2.d0 * gl%J2X1_PCSAFTQ(2,i,j,k)) * gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            !Jq_2_ijk_xk_D bereits mit T multipliziert
                            sum10 = sum10 + (-2.d0 * gl%J2X1_PCSAFTQ(1,i,j,k) + gl%J2X1_PCSAFTQ(4,i,j,k)) * gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(j)) 
                            !Jq_2_ijk_xk bereits mit T multipliziert
                            
                        end if
                    end do
                end do
            part1 = gl%QPCSAFTQ(k)**4 * gl%nPCSAFTQ(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**7 * gl%mPCSAFT(k)**2)
            part3a = (- gl%J2X1_PCSAFTQ(6,k,k,k) + 2.d0 * gl%J2X1_PCSAFTQ(2,k,k,k)) * part1
                    !Jq_2_kkk_xk_D bereits mit T multipliziert
            part3b = (2.d0 * gl%J2X1_PCSAFTQ(1,k,k,k) - gl%J2X1_PCSAFTQ(4,k,k,k)) * part1
                    !Jq_2_kkk_xk bereits mit T multipliziert
            part2 = 2.d0 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTQ(6,k) =  -0.5625d0 * piPCSAFT * DENS * 1.d-38 * (part3a + part2 * (gl%molfractions(k) * Sum1 + Sum2) &
                & - part2 * (gl%molfractions(k) * Sum3 + Sum4) + Sum9 + part3b &
                & + part2 * (gl%molfractions(k) * Sum5 + Sum6) - part2 * (gl%molfractions(k) * Sum7 + Sum8) + Sum10)
             !1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
    
    ! 7: 3RD MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires j_2qq
    if (GETDERAQQ(7) .eq. 1) then
      do k = 1, gl%ncomp
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
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (T*kbol*sigma_kj**7 * gl%mPCSAFT(j))
                    sum1 = sum1 + (gl%J2X1_PCSAFTQ(7,j,k,k) - 2.d0 * gl%J2X1_PCSAFTQ(6,j,k,k) + 2.d0 * gl%J2X1_PCSAFTQ(2,j,k,k)) * sum
                    sum2 = sum2 + (gl%J2_PCSAFTQ(7,j,k) - 2.d0 * gl%J2_PCSAFTQ(6,j,k) + 2.d0 * gl%J2_PCSAFTQ(2,j,k)) * sum
                    sum3 = sum3 + (gl%J2X1_PCSAFTQ(6,j,k,k) - gl%J2X1_PCSAFTQ(2,j,k,k)) * sum
                    sum4 = sum4 + (gl%J2_PCSAFTQ(6,j,k) - gl%J2_PCSAFTQ(2,j,k)) * sum
                    sum5 = sum5 + gl%J2X1_PCSAFTQ(2,j,k,k) * sum
                    sum6 = sum6 + gl%J2_PCSAFTQ(2,j,k) * sum
                    sum7 = sum7 + (2.d0 * gl%J2X1_PCSAFTQ(1,j,k,k) - 2.d0 * gl%J2X1_PCSAFTQ(4,j,k,k) + gl%J2X1_PCSAFTQ(5,j,k,k)) * sum
                    sum8 = sum8 + (2.d0 * gl%J2_PCSAFTQ(1,j,k) - 2.d0 * gl%J2_PCSAFTQ(4,j,k) + gl%J2_PCSAFTQ(5,j,k)) * sum
                    sum9 = sum9 + (-gl%J2X1_PCSAFTQ(1,j,k,k) + gl%J2X1_PCSAFTQ(4,j,k,k)) * sum
                    sum10 = sum10 + (-gl%J2_PCSAFTQ(1,j,k) + gl%J2_PCSAFTQ(4,j,k)) * sum
                    sum11 = sum11 + gl%J2X1_PCSAFTQ(1,j,k,k) * sum
                    sum12 = sum12 + gl%J2_PCSAFTQ(1,j,k) * sum
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            part13 = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            sum13 = sum13 + (gl%J2X1_PCSAFTQ(7,i,j,k) - 4.d0 * gl%J2X1_PCSAFTQ(6,i,j,k) + 6.d0 * gl%J2X1_PCSAFTQ(2,i,j,k)) * part13
                            sum14 = sum14 + (6.d0 * gl%J2X1_PCSAFTQ(1,i,j,k) - 4.d0 * gl%J2X1_PCSAFTQ(4,i,j,k) + gl%J2X1_PCSAFTQ(5,i,j,k)) * part13
                            
                        end if
                    end do
                end do
            part1 = gl%QPCSAFTQ(k)**4 * gl%nPCSAFTQ(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**7 * gl%mPCSAFT(k)**2)
            part3a = (-gl%J2X1_PCSAFTQ(7,k,k,k) + 4.d0 * gl%J2X1_PCSAFTQ(6,k,k,k) - 6.d0 * gl%J2X1_PCSAFTQ(2,k,k,k)) * part1
            part3b = (-6.d0 * gl%J2X1_PCSAFTQ(1,k,k,k) + 4.d0 * gl%J2X1_PCSAFTQ(4,k,k,k) - gl%J2X1_PCSAFTQ(5,k,k,k)) * part1
            part2 = 2.d0 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTQ(7,k) =   -0.5625d0 * piPCSAFT * DENS * 1.d-38 * (part3a + part2 * (gl%molfractions(k) * Sum1 + Sum2) - 2.d0 * part2 * (gl%molfractions(k) * Sum3 + Sum4) &
                & + 2.d0 * part2 * (gl%molfractions(k) * Sum5 + Sum6) + Sum13 + part3b + part2 * (gl%molfractions(k) * Sum7 + Sum8) - 2.d0 * part2 * (gl%molfractions(k) * Sum9 + Sum10) &
                & + 2.d0 * part2 * (gl%molfractions(k)*Sum11 + Sum12) + Sum14)
             !1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
    
    ! 8: 3RD DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires j2_qq
    if (GETDERAQQ(8) .eq. 1) then
       do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
            sum4 = 0.d0
            sum5 = 0.d0
            sum6 = 0.d0
                do j = 1, gl%ncomp
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (T*kbol*sigma_kj**7 * gl%mPCSAFT(j))
                    sum1 = sum1 + gl%J2X1_PCSAFTQ(3,j,k,k) * sum
                    sum2 = sum2 + gl%J2_PCSAFTQ(3,j,k) * sum
                    sum3 = sum3 + gl%J2X1_PCSAFTQ(8,j,k,k) * sum
                    sum4 = sum4 + gl%J2_PCSAFTQ(8,j,k) * sum
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            part1 = gl%QPCSAFTQ(i)**2 * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(i) * gl%nPCSAFTQ(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            sum6 = sum6 + gl%J2X1_PCSAFTQ(3,i,j,k) * part1
                            sum5 = sum5 + gl%J2X1_PCSAFTQ(8,i,j,k) * part1
                        end if
                    end do
                end do
            part1 = gl%QPCSAFTQ(k)**4*gl%nPCSAFTQ(k)**2*gl%molfractions(k)**2/(T**2*kbol**2*gl%sigPCSAFT(k)**7*gl%mPCSAFT(k)**2)    
            part3a = - 3.d0*gl%J2X1_PCSAFTQ(3,k,k,k) * part1
            part3b = - gl%J2X1_PCSAFTQ(8,k,k,k) * part1
            part2 = 2.d0 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTQ(8,k) = - piPCSAFT * 0.5625d0 * 1.d-38 * DENS * (part3a + 3.d0 * part2 * (gl%molfractions(k) * Sum1 + Sum2) &
                & + part3b + part2 * (gl%molfractions(k) * Sum3 + Sum4) + Sum5 + 3.d0 * Sum6)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
   
    
    
    ! 9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires 
    if (GETDERAQQ(9) .eq. 1) then
          do k = 1, gl%ncomp
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
                sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                sum = gl%molfractions(j) * gl%QPCSAFTQ(j)**2 * gl%nPCSAFTQ(j) / (T*kbol*sigma_kj**7 * gl%mPCSAFT(j))
                sum1 = sum1 + (- 6.d0 * gl%J2X1_PCSAFTQ(1,j,k,k) + 6.d0 * gl%J2X1_PCSAFTQ(4,j,k,k) - 3.d0 * gl%J2X1_PCSAFTQ(5,j,k,k) + gl%J2X1_PCSAFTQ(9,j,k,k)) * sum
                sum2 = sum2 + (- 6.d0 * gl%J2_PCSAFTQ(1,j,k) + 6.d0 * gl%J2_PCSAFTQ(4,j,k) - 3.d0 * gl%J2_PCSAFTQ(5,j,k) + gl%J2_PCSAFTQ(9,j,k)) * sum
                sum3 = sum3 + (2.d0 * gl%J2X1_PCSAFTQ(1,j,k,k) - 2.d0 * gl%J2X1_PCSAFTQ(4,j,k,k) + gl%J2X1_PCSAFTQ(5,j,k,k)) * sum
                sum4 = sum4 + (2.d0 * gl%J2_PCSAFTQ(1,j,k) - 2.d0 * gl%J2_PCSAFTQ(4,j,k) + gl%J2_PCSAFTQ(5,j,k)) * sum
                sum5 = sum5 + (- gl%J2X1_PCSAFTQ(1,j,k,k) + gl%J2X1_PCSAFTQ(4,j,k,k)) * sum
                sum6 = sum6 + (- gl%J2_PCSAFTQ(1,j,k) + gl%J2_PCSAFTQ(4,j,k)) * sum
                sum7 = sum7 + gl%J2X1_PCSAFTQ(1,j,k,k) * sum
                sum8 = sum8 + gl%J2_PCSAFTQ(1,j,k) * sum
                do i = 1, gl%ncomp
                    if (i /= k .and. j /= k) then 
                        sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                        sum9 = sum9 + (- 24.d0 * gl%J2X1_PCSAFTQ(1,i,j,k) + 18.d0 * gl%J2X1_PCSAFTQ(4,i,j,k) - 6.d0 * gl%J2X1_PCSAFTQ(5,i,j,k) + gl%J2X1_PCSAFTQ(9,i,j,k)) &
                            & *gl%QPCSAFTQ(i)**2*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%molfractions(i)*gl%molfractions(j) &
                            & /(T**2*kbol**2*sigma_ij**7*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                    end if
                end do
            end do
            part3 = (24.d0 * gl%J2X1_PCSAFTQ(1,k,k,k) - 18.d0 * gl%J2X1_PCSAFTQ(4,k,k,k) + 6.d0 * gl%J2X1_PCSAFTQ(5,k,k,k) - gl%J2X1_PCSAFTQ(9,k,k,k)) * gl%QPCSAFTQ(k)**4 * gl%nPCSAFTQ(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**7 * gl%mPCSAFT(k)**2) 
            part2 = 2.d0 * gl%QPCSAFTQ(k)**2 * gl%nPCSAFTQ(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTQ(9,k) =  -0.5625d0 * piPCSAFT * DENS * 1.d-38 * (part3 + part2 * (gl%molfractions(k) * Sum1 + Sum2) &
                & - 3.d0 * part2 * (gl%molfractions(k)*Sum3 + Sum4) + 6.d0 * part2 *(gl%molfractions(k)*Sum5 + Sum6) &
                & - 6.d0 * part2 * (gl%molfractions(k)*Sum7 + Sum8) + Sum9)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires j_2qq
    if (GETDERAQQ(10) .eq. 1) then
        do k = 1, gl%ncomp
            help = 0.d0
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
            sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
            help = gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(j)*gl%molfractions(j)/(T*kbol*sigma_kj**7*gl%mPCSAFT(j))
            sum1 = sum1 + (gl%J2X1_PCSAFTQ(6,j,k,k) - gl%J2X1_PCSAFTQ(2,j,k,k)) * help
            sum2 = sum2 + (gl%J2_PCSAFTQ(6,j,k) - gl%J2_PCSAFTQ(2,j,k)) * help
            sum3 = sum3 + gl%J2X1_PCSAFTQ(2,j,k,k) * help
            sum4 = sum4 + gl%J2_PCSAFTQ(2,j,k) * help
            sum5 = sum5 + (gl%J2X1_PCSAFTQ(10,j,k,k) - gl%J2X1_PCSAFTQ(3,j,k,k)) * help
            sum6 = sum6 + (gl%J2_PCSAFTQ(10,j,k) - gl%J2_PCSAFTQ(3,j,k)) * help
            sum7 = sum7 + gl%J2X1_PCSAFTQ(3,j,k,k) * help
            sum8 = sum8 + gl%J2_PCSAFTQ(3,j,k) * help
                do i = 1, gl%ncomp
                
                    if (i /= k .and. j /= k) then 
                    sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                    sum9 = sum9 + (gl%J2X1_PCSAFTQ(10,i,j,k) - 2.d0 * gl%J2X1_PCSAFTQ(3,i,j,k))*gl%QPCSAFTQ(i)**2*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%molfractions(i)*gl%molfractions(j) &
                        & / (T**2*kbol**2*sigma_ij**7*gl%mPCSAFT(i)*gl%mPCSAFT(j))  
                    sum10 = sum10 + (gl%J2X1_PCSAFTQ(6,i,j,k) - 2.d0 * gl%J2X1_PCSAFTQ(2,i,j,k))*gl%QPCSAFTQ(i)**2*gl%QPCSAFTQ(j)**2*gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j)*gl%molfractions(i)*gl%molfractions(j) &
                        & /(T**2*kbol**2*sigma_ij**7*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                    end if
                end do
            end do
            part3a = (- 2.d0 * gl%J2X1_PCSAFTQ(6,k,k,k) + 4.d0 * gl%J2X1_PCSAFTQ(2,k,k,k)) *gl%QPCSAFTQ(k)**4*gl%nPCSAFTQ(k)**2*gl%molfractions(k)**2 &
                        & /(T**2*kbol**2*gl%sigPCSAFT(k)**7*gl%mPCSAFT(k)**2) 
            part3b = (- gl%J2X1_PCSAFTQ(10,k,k,k) + 2.d0 * gl%J2X1_PCSAFTQ(3,k,k,k)) *gl%QPCSAFTQ(k)**4*gl%nPCSAFTQ(k)**2*gl%molfractions(k)**2 &
                        & /(T**2*kbol**2*gl%sigPCSAFT(k)**7*gl%mPCSAFT(k)**2)
                        
            help2 = 2.d0 * gl%QPCSAFTQ(k)**2*gl%nPCSAFTQ(k) / (T*kbol*gl%mPCSAFT(k))
    
            gl%A2X1_PCSAFTQ(10,k) =  -piPCSAFT * 0.5625d0 * DENS * 1.d-38 * (part3a + 2.d0 * help2 * (gl%molfractions(k) * Sum1 + Sum2) &
                                & - 2.d0 * help2 *(gl%molfractions(k)*Sum3 + Sum4) &
                                & + part3b + help2 * (gl%molfractions(k)*Sum5 + Sum6) &
                                & - help2 * (gl%molfractions(k)*Sum7 + Sum8) &
                                & + Sum9 + 2.d0*Sum10)
            
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
    !! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    !! requires j_2qq
    !if (GETDERAQQ(10) .eq. 1) then
    !    do k = 1, ncomp
    !   
    !        sum1 = 0.d0
    !        sum2 = 0.d0
    !        sum3 = 0.d0
    !        sum4 = 0.d0
    !        sum5 = 0.d0
    !        sum6 = 0.d0
    !        sum7 = 0.d0
    !        sum8 = 0.d0
    !        sum9 = 0.d0
    !        sum10 = 0.d0
    !            do j = 1, ncomp
    !                sigma_kj = (sigPCSAFT(k)+sigPCSAFT(j))*0.5d0
    !                sum = molfractions(j) * QPCSAFTQ(j)**2 * nPCSAFTQ(j) / (T*kbol*sigma_kj**7 * mPCSAFT(j))
    !                sum1 = sum1 + (J2X1_PCSAFTQ(6,j,k,k) - J2X1_PCSAFTQ(2,j,k,k)) * sum
    !                sum2 = sum2 + (J2_PCSAFTQ(6,j,k) - J2_PCSAFTQ(2,j,k)) * sum
    !                sum3 = sum3 + J2X1_PCSAFTQ(2,j,k,k) * sum
    !                sum4 = sum4 + J2_PCSAFTQ(2,j,k) * sum
    !                sum5 = sum5 + (J2X1_PCSAFTQ(10,j,k,k) - J2X1_PCSAFTQ(3,j,k,k)) * sum
    !                sum6 = sum6 + (J2_PCSAFTQ(10,j,k) - J2_PCSAFTQ(3,j,k)) * sum
    !                sum7 = sum7 + J2X1_PCSAFTQ(3,j,k,k) * sum
    !                sum8 = sum8 + J2_PCSAFTQ(3,j,k) * sum
    !                do i = 1, ncomp
    !                    if (i /= k .and. j /= k) then 
    !                        sigma_ij = (sigPCSAFT(i)+sigPCSAFT(j))*0.5d0
    !                        part1 = QPCSAFTQ(i)**2 * QPCSAFTQ(j)**2 * nPCSAFTQ(i) * nPCSAFTQ(j) * molfractions(i) * molfractions(j) / (T**2 * kbol**2 * sigma_ij**7 * mPCSAFT(i) * mPCSAFT(j))
    !                        sum9 = sum9 + (J2X1_PCSAFTQ(10,i,j,k) - 2.d0 * J2X1_PCSAFTQ(3,i,j,k)) * part1
    !                        sum10 = sum10 + (J2X1_PCSAFTQ(6,i,j,k) - 2.d0 * J2X1_PCSAFTQ(2,i,j,k)) * part1
    !                    end if
    !                end do
    !            end do
    !        part3 = QPCSAFTQ(k)**4 * nPCSAFTQ(k)**2 * molfractions(k)**2 / (T**2 * kbol**2 * sigPCSAFT(k)**7 * mPCSAFT(k)**2)        
    !        part3a = (- 2.d0 * J2X1_PCSAFTQ(6,k,k,k) + 4.d0 * J2X1_PCSAFTQ(2,k,k,k)) * part3
    !        part3b = (- J2X1_PCSAFTQ(10,k,k,k) + 2.d0 * J2X1_PCSAFTQ(3,k,k,k)) * part3
    !        part2 = 2.d0 * QPCSAFTQ(k)**2 * nPCSAFTQ(k) / (T * kbol * mPCSAFT(k))
    !        
    !        A2X1_PCSAFTQ(10,k) = -0.5625d0 * piPCSAFT * 1.d-38 * DENS * (&
    !        & part3a &
    !        & + 2.d0 * part2 * (molfractions(k) * Sum1 + Sum2) &
    !        & - 2.d0 * part2 * (molfractions(k) * Sum3 + Sum4) &
    !        & + part3b &
    !        & + part2 * (molfractions(k)*Sum5 + Sum6) &
    !        & - part2 * (molfractions(k) * Sum7 + Sum8) &
    !        & + Sum9 + 2.d0 * Sum10)
    !        
    !        ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
    !    end do
    !end if
!DEC$ END IF
end subroutine A2X1DERIVS


subroutine A2X2DERIVS(gl,T,DENS,GETDERAQQ)
    ! a_2qq: second-order perturbation term of the quadrupolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.9 in Gross, Sadowski 2005:
    ! a_2qq =-pi*(3/4)**2*rho*sum_j*sum_i(x_i*x_j*epsk_ii/T*espk_jj/T*sig_ii**5*sig_jj**5/sig_ij**7*n_Qi*n_Qj*Q_i**2*Q_j**2*J_2ij
    ! dependent on D and T and x_i







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERAQQ
    double precision, intent (in) :: T, DENS
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sigma_ij, sigma_kj
    double precision :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, sum14, part13, part3, part3a, part3b, part3c, sum, help, help3, part1, part2
    integer :: i, j, xi, k, o, u
    double precision :: DP, DP2, DM, DM2, TP, TP2, TM, TM2
    double precision :: DPx, DP2x, DMx, DM2x, TPx, TP2x, TMx, TM2x
    double precision :: DPxx, DP2xx, DMxx, DM2xx, TPxx, TP2xx, TMxx, TM2xx
    double precision :: DELTA, DEL, DELT
    double precision, dimension(gl%ncomp,gl%ncomp) :: A2X2, A2X2DP, A2X2DP2, A2X2DM, A2X2DM2, A2X2TP, A2X2TP2, A2X2TM, A2X2TM2, A2X2PP, A2X2PM, A2X2MP, A2X2MM

    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    DELTA = 1.d-7
    DELTA = 1.d-5
    
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
    
    DELT = 1.d-4 ! Gemischte Ableitung
    
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
        
        call base_A2X2(gl,T,DENS,A2X2)
        
        gl%A2X2_PCSAFTQ(1,:,:) = A2X2

    end if
    
    
    !  2: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires j_2qq
    if (GETDERAQQ(2) .eq. 1) then   
               
        call base_A2X2(gl,T,DP,A2X2DP)
        call base_A2X2(gl,T,DP2,A2X2DP2)
        call base_A2X2(gl,T,DM,A2X2DM)
        call base_A2X2(gl,T,DM2,A2X2DM2)
        
        gl%A2X2_PCSAFTQ(2,:,:) = (8.d0*A2X2DP - 8.d0*A2X2DM - A2X2DP2 + A2X2DM2) / (12.d0 * DELTA) 
    

    end if
    

   ! 3: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
   ! requires j_2qq
    if (GETDERAQQ(3) .eq. 1) then
        
        call base_A2X2(gl,T,DENS,A2X2)    
        call base_A2X2(gl,T,DPx,A2X2DP)
        call base_A2X2(gl,T,DP2x,A2X2DP2)
        call base_A2X2(gl,T,DMx,A2X2DM)
        call base_A2X2(gl,T,DM2x,A2X2DM2)
        
        gl%A2X2_PCSAFTQ(3,:,:) = (-A2X2DM2 + 16.d0*A2X2DM - 30.d0*A2X2 + 16.d0*A2X2DP - A2X2DP2) / (12.d0*DEL**2) 
    
    end if
    
   
     
    ! 4: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires j_2qq
    if (GETDERAQQ(4) .eq. 1) then
        call base_A2X2(gl,TP,DENS,A2X2TP)
        call base_A2X2(gl,TP2,DENS,A2X2TP2)
        call base_A2X2(gl,TM,DENS,A2X2TM)
        call base_A2X2(gl,TM2,DENS,A2X2TM2)
        
        gl%A2X2_PCSAFTQ(4,:,:) = (8.d0*A2X2TP - 8.d0*A2X2TM - A2X2TP2 + A2X2TM2) / (12.d0 * DELTA) 
    end if
    

    ! 5: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    !requires j_2qq
    if (GETDERAQQ(5) .eq. 1) then
    
        call base_A2X2(gl,T,DENS,A2X2)    
        call base_A2X2(gl,TPx,DENS,A2X2TP)
        call base_A2X2(gl,TP2x,DENS,A2X2TP2)
        call base_A2X2(gl,TMx,DENS,A2X2TM)
        call base_A2X2(gl,TM2x,DENS,A2X2TM2)
        
        gl%A2X2_PCSAFTQ(5,:,:) = (-A2X2TM2 + 16.d0*A2X2TM - 30.d0*A2X2 + 16.d0*A2X2TP - A2X2TP2) / (12.d0*DEL**2) !https://en.wikipedia.org/wiki/Five-point_stencil
    
    end if
       
     
    ! 6: 2ND MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_2qq
    if (GETDERAQQ(6) .eq. 1) then
    
        call base_A2X2(gl,TP,DPxx,A2X2PP)
        call base_A2X2(gl,TM,DPxx,A2X2MP)
        call base_A2X2(gl,TM,DMxx,A2X2MM)
        call base_A2X2(gl,TP,DMxx,A2X2PM)
        
        gl%A2X2_PCSAFTQ(6,:,:) = (A2X2PP - A2X2PM - A2X2MP + A2X2MM) / (4.d0*DELT*DELTA) 
    end if
    

!DEC$ END IF
end subroutine A2X2DERIVS


subroutine A2X1DERIVS_D(gl,T,DENS,GETDERADD)
   
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
    !output: ADD_PCSAFT (module variable)
    !working variables
    double precision :: sigma_ij, sigma_kj
    double precision :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, sum14, part13, part3, part3a, part3b, part3c, sum, help,help2, help3, part1, part2
    integer :: i, j, xi,  k
    double precision, dimension(2) :: dd
    
 
    integer:: errorfld

!DEC$ IF DEFINED(WO_PCSAFT)
    errorfld = -7878
!DEC$ ELSE
    
    !calculate the derivatives of a_dd
    ! 1: a_2dd
    if (GETDERADD(1) .eq. 1) then
        do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
                do j = 1, gl%ncomp
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (sigma_kj**3 * gl%mPCSAFT(j))
                    sum1 = sum1 + sum*gl%J2X1_PCSAFTD(1,j,k,k) !*J2X1_PCSAFTD(1,j,k,k) sollte egal sein
                    sum2 = sum2 + sum*gl%J2_PCSAFTD(1,j,k)
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sum3 = sum3 + gl%molfractions(i)*gl%molfractions(j) * gl%nPCSAFTD(i)*gl%nPCSAFTD(j) * gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(i)*gl%mPCSAFT(j)) * gl%J2X1_PCSAFTD(1,i,j,k)
                        end if
                    end do
                end do
            part3 = gl%molfractions(k)**2 * gl%nPCSAFTD(k)**2 * gl%MyPCSAFTD(k)**4 / (kbol**2 * T**2 * gl%sigPCSAFT(k)**3 * gl%mPCSAFT(k)**2) * gl%J2X1_PCSAFTD(1,k,k,k)
            gl%A2X1_PCSAFTD(1,k) = - piPCSAFT*DENS*1.d-38  * (2.d0 * gl%MyPCSAFTD(k)**2 *gl%nPCSAFTD(k) / (T**2*kbol**2*gl%mPCSAFT(k)) * (gl%molfractions(k)*sum1 + sum2) - part3 + sum3)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
    
    !  2: 1ST DERIVATIVE OF a_2dd WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires j_2dd
    if (GETDERADD(2) .eq. 1) then   
        do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
            sum4 = 0.d0
            sum5 = 0.d0
            sum6 = 0.d0
                do j = 1, gl%ncomp
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (T*kbol*sigma_kj**3 * gl%mPCSAFT(j))
                    sum1 = sum1 + gl%J2X1_PCSAFTD(2,j,k,k) * sum
                    sum2 = sum2 + gl%J2_PCSAFTD(2,j,k) * sum
                    sum3 = sum3 + gl%J2X1_PCSAFTD(1,j,k,k) * sum
                    sum4 = sum4 + gl%J2_PCSAFTD(1,j,k) * sum
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            part1 = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            sum6 = sum6 + gl%J2X1_PCSAFTD(2,i,j,k) * part1
                            sum5 = sum5 + gl%J2X1_PCSAFTD(1,i,j,k) * part1
                        end if
                    end do
                end do
            part3a = - gl%J2X1_PCSAFTD(2,k,k,k) * gl%MyPCSAFTD(k)**4 * gl%nPCSAFTD(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**3 * gl%mPCSAFT(k)**2) 
            part3b = - gl%J2X1_PCSAFTD(1,k,k,k) * gl%MyPCSAFTD(k)**4 * gl%nPCSAFTD(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**3 * gl%mPCSAFT(k)**2)
            part2 = 2.d0 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTD(2,k) = - piPCSAFT * DENS * 1.d-38 * (part3a + part2 * (gl%molfractions(k) * Sum1 + Sum2) + Sum6 + part3b + part2 * (gl%molfractions(k) * Sum3 + Sum4)  + Sum5)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
   


     
   ! 3: 2ND DERIVATIVE OF a_2dd WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
   ! requires j_2dd
    if (GETDERADD(3) .eq. 1) then
        do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
            sum4 = 0.d0
            sum5 = 0.d0
            sum6 = 0.d0
                do j = 1, gl%ncomp
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (T*kbol*sigma_kj**3 * gl%mPCSAFT(j))
                    sum1 = sum1 + gl%J2X1_PCSAFTD(2,j,k,k)*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(j)*gl%molfractions(j)/(T*kbol*sigma_kj**3*gl%mPCSAFT(j))
                    sum2 = sum2 + gl%J2_PCSAFTD(2,j,k)*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(j)*gl%molfractions(j)/(T*kbol*sigma_kj**3*gl%mPCSAFT(j))
                    sum3 = sum3 + gl%J2X1_PCSAFTD(3,j,k,k)*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(j)*gl%molfractions(j)/(T*kbol*sigma_kj**3*gl%mPCSAFT(j))
                    sum4 = sum4 + gl%J2_PCSAFTD(3,j,k)*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(j)*gl%molfractions(j)/(T*kbol*sigma_kj**3*gl%mPCSAFT(j))
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            part1 = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            sum6 = sum6 + gl%J2X1_PCSAFTD(2,i,j,k)*gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j)/(T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                            sum5 = sum5 + gl%J2X1_PCSAFTD(3,i,j,k)*gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j)/(T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                        end if
                    end do
                end do
            part3a = - gl%J2X1_PCSAFTD(2,k,k,k)*gl%MyPCSAFTD(k)**4*gl%nPCSAFTD(k)**2*gl%molfractions(k)**2/(T**2*kbol**2*gl%sigPCSAFT(k)**3*gl%mPCSAFT(k)**2)
            part3b = - gl%J2X1_PCSAFTD(3,k,k,k) * gl%MyPCSAFTD(k)**4 * gl%nPCSAFTD(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**3 * gl%mPCSAFT(k)**2)
            part2 = 2.d0 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTD(3,k) = - piPCSAFT * DENS * 1.d-38 *(2.d0 * part3a + 2.d0 * part2 * (gl%molfractions(k) * Sum1 + Sum2) + (part3b + part2 * (gl%molfractions(k) * Sum3 + Sum4) + Sum5) + 2.d0 * Sum6)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
   
     
    ! 4: 1ST DERIVATIVE OF a_2dd WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires j_2dd
    if (GETDERADD(4) .eq. 1) then
      do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
            sum4 = 0.d0
            sum5 = 0.d0
            
                do j = 1, gl%ncomp
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (T*kbol*sigma_kj**3 * gl%mPCSAFT(j))
                    sum1 = sum1 + (gl%J2X1_PCSAFTD(4,j,k,k)- gl%J2X1_PCSAFTD(1,j,k,k)) * sum
                    sum2 = sum2 + (gl%J2_PCSAFTD(4,j,k) - gl%J2_PCSAFTD(1,j,k)) * sum 
                    sum3 = sum3 + gl%J2X1_PCSAFTD(1,j,k,k) * sum
                    sum4 = sum4 + gl%J2_PCSAFTD(1,j,k) * sum
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sum5 = sum5 + (gl%J2X1_PCSAFTD(4,i,j,k) - 2.d0 * gl%J2X1_PCSAFTD(1,i,j,k)) * gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%molfractions(i) * gl%molfractions(j) &
                                &  / (T**2 * kbol**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                        end if
                    end do
                end do
            part3 = (2.d0 * gl%J2X1_PCSAFTD(1,k,k,k) - gl%J2X1_PCSAFTD(4,k,k,k)) * gl%MyPCSAFTD(k)**4 * gl%nPCSAFTD(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**3 * gl%mPCSAFT(k)**2)
            part2 = 2.d0 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTD(4,k) = -piPCSAFT * DENS * 1.d-38 * (part3 + part2 * (gl%molfractions(k) * Sum1 + Sum2)   &
                & - part2 * (gl%molfractions(k) * Sum3 + Sum4) + Sum5)
             !1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    

    ! 5: 2ND DERIVATIVE OF a_2dd WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    !requires j_2dd
    if (GETDERADD(5) .eq. 1) then
        do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
            sum4 = 0.d0
            sum5 = 0.d0
            sum6 = 0.d0
            sum7 = 0.d0
            
            do j = 1, gl%ncomp
                sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                sum = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (T*kbol*sigma_kj**3 * gl%mPCSAFT(j))
                sum1 = sum1 + (2.d0 * gl%J2X1_PCSAFTD(1,j,k,k) - 2.d0 * gl%J2X1_PCSAFTD(4,j,k,k) + gl%J2X1_PCSAFTD(5,j,k,k)) * sum
                sum2 = sum2 + (2.d0 * gl%J2_PCSAFTD(1,j,k) - 2.d0 * gl%J2_PCSAFTD(4,j,k) + gl%J2_PCSAFTD(5,j,k)) * sum
                sum3 = sum3 + (- gl%J2X1_PCSAFTD(1,j,k,k) + gl%J2X1_PCSAFTD(4,j,k,k)) * sum
                sum4 = sum4 + (- gl%J2_PCSAFTD(1,j,k) + gl%J2_PCSAFTD(4,j,k)) * sum 
                sum5 = sum5 + gl%J2X1_PCSAFTD(1,j,k,k) * sum
                sum6 = sum6 + gl%J2_PCSAFTD(1,j,k) * sum
                do i = 1, gl%ncomp
                    if (i /= k .and. j /= k) then 
                        sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                        sum7 = sum7 + gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * (6.d0 * gl%J2X1_PCSAFTD(1,i,j,k) - 4.d0 * gl%J2X1_PCSAFTD(4,i,j,k) + gl%J2X1_PCSAFTD(5,i,j,k)) &
                            & * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                    end if
                end do
            end do
            part3 = (- 6.d0 * gl%J2X1_PCSAFTD(1,k,k,k) + 4.d0 * gl%J2X1_PCSAFTD(4,k,k,k) - gl%J2X1_PCSAFTD(5,k,k,k)) / (T**2 * kbol**2 * gl%sigPCSAFT(k)**3 * gl%mPCSAFT(k)**2) &
            & * gl%MyPCSAFTD(k)**4 * gl%nPCSAFTD(k)**2 * gl%molfractions(k)**2
            part2 = 2.d0 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTD(5,k) =  -piPCSAFT * DENS * 1.d-38 * (part3 + part2 * (gl%molfractions(k) * Sum1 + Sum2) &
                & - 2.d0 * part2 * (gl%molfractions(k) * Sum3 + Sum4) &
                & + 2.d0 * part2 * (gl%molfractions(k) * Sum5 + Sum6) + Sum7)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
       
     
    ! 6: 2ND MIXED DERIVATIVE OF a_2dd WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_2dd
    if (GETDERADD(6) .eq. 1) then
         do k = 1, gl%ncomp
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
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (T*kbol*sigma_kj**3 * gl%mPCSAFT(j))
                    sum1 = sum1 + (gl%J2X1_PCSAFTD(6,j,k,k) - gl%J2X1_PCSAFTD(2,j,k,k)) * sum
                        !Jq_2_kjk_xk_D bereits mit T multipliziert

                    sum2 = sum2 + (gl%J2_PCSAFTD(6,j,k) - gl%J2_PCSAFTD(2,j,k)) * sum
                        !Jq_2_kj_D bereits mit T multipliziert
                    
                    sum3 = sum3 + gl%J2X1_PCSAFTD(2,j,k,k) * sum
                    sum4 = sum4 + gl%J2_PCSAFTD(2,j,k) * sum
                    sum5 = sum5 + (-gl%J2X1_PCSAFTD(1,j,k,k) + gl%J2X1_PCSAFTD(4,j,k,k)) * sum 
                        !Jq_2_kjk_xk bereits mit T multipliziert
                    
                    sum6 = sum6 + (-gl%J2_PCSAFTD(1,j,k) +  gl%J2_PCSAFTD(4,j,k)) * sum
                        !Jq_2_kj bereits mit T multipliziert
                    
                    sum7 = sum7 + gl%J2X1_PCSAFTD(1,j,k,k) * sum
                    sum8 = sum8 + gl%J2_PCSAFTD(1,j,k) * sum
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            sum9 = sum9 + (gl%J2X1_PCSAFTD(6,i,j,k) - 2.d0 * gl%J2X1_PCSAFTD(2,i,j,k)) * gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            !Jq_2_ijk_xk_D bereits mit T multipliziert
                            sum10 = sum10 + (-2.d0 * gl%J2X1_PCSAFTD(1,i,j,k) + gl%J2X1_PCSAFTD(4,i,j,k)) * gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j)) 
                            !Jq_2_ijk_xk bereits mit T multipliziert
                            
                        end if
                    end do
                end do
            part1 = gl%MyPCSAFTD(k)**4 * gl%nPCSAFTD(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**3 * gl%mPCSAFT(k)**2)
            part3a = (- gl%J2X1_PCSAFTD(6,k,k,k) + 2.d0 * gl%J2X1_PCSAFTD(2,k,k,k)) * part1
                    !Jq_2_kkk_xk_D bereits mit T multipliziert
            part3b = (2.d0 * gl%J2X1_PCSAFTD(1,k,k,k) - gl%J2X1_PCSAFTD(4,k,k,k)) * part1
                    !Jq_2_kkk_xk bereits mit T multipliziert
            part2 = 2.d0 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTD(6,k) =  -piPCSAFT * DENS * 1.d-38 * (part3a + part2 * (gl%molfractions(k) * Sum1 + Sum2) &
                & - part2 * (gl%molfractions(k) * Sum3 + Sum4) + Sum9 + part3b &
                & + part2 * (gl%molfractions(k) * Sum5 + Sum6) - part2 * (gl%molfractions(k) * Sum7 + Sum8) + Sum10)
             !1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
    
    ! 7: 3RD MIXED DERIVATIVE OF a_2dd WITH RESPECT TO D, T, AND T, MULTIPLIED BY T*T*D
    ! requires j_2dd
    if (GETDERADD(7) .eq. 1) then
      do k = 1, gl%ncomp
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
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (T*kbol*sigma_kj**3 * gl%mPCSAFT(j))
                    sum1 = sum1 + (gl%J2X1_PCSAFTD(7,j,k,k) - 2.d0 * gl%J2X1_PCSAFTD(6,j,k,k) + 2.d0 * gl%J2X1_PCSAFTD(2,j,k,k)) * sum
                    sum2 = sum2 + (gl%J2_PCSAFTD(7,j,k) - 2.d0 * gl%J2_PCSAFTD(6,j,k) + 2.d0 * gl%J2_PCSAFTD(2,j,k)) * sum
                    sum3 = sum3 + (gl%J2X1_PCSAFTD(6,j,k,k) - gl%J2X1_PCSAFTD(2,j,k,k)) * sum
                    sum4 = sum4 + (gl%J2_PCSAFTD(6,j,k) - gl%J2_PCSAFTD(2,j,k)) * sum
                    sum5 = sum5 + gl%J2X1_PCSAFTD(2,j,k,k) * sum
                    sum6 = sum6 + gl%J2_PCSAFTD(2,j,k) * sum
                    sum7 = sum7 + (2.d0 * gl%J2X1_PCSAFTD(1,j,k,k) - 2.d0 * gl%J2X1_PCSAFTD(4,j,k,k) + gl%J2X1_PCSAFTD(5,j,k,k)) * sum
                    sum8 = sum8 + (2.d0 * gl%J2_PCSAFTD(1,j,k) - 2.d0 * gl%J2_PCSAFTD(4,j,k) + gl%J2_PCSAFTD(5,j,k)) * sum
                    sum9 = sum9 + (-gl%J2X1_PCSAFTD(1,j,k,k) + gl%J2X1_PCSAFTD(4,j,k,k)) * sum
                    sum10 = sum10 + (-gl%J2_PCSAFTD(1,j,k) + gl%J2_PCSAFTD(4,j,k)) * sum
                    sum11 = sum11 + gl%J2X1_PCSAFTD(1,j,k,k) * sum
                    sum12 = sum12 + gl%J2_PCSAFTD(1,j,k) * sum
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            part13 = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            sum13 = sum13 + (gl%J2X1_PCSAFTD(7,i,j,k) - 4.d0 * gl%J2X1_PCSAFTD(6,i,j,k) + 6.d0 * gl%J2X1_PCSAFTD(2,i,j,k)) * part13
                            sum14 = sum14 + (6.d0 * gl%J2X1_PCSAFTD(1,i,j,k) - 4.d0 * gl%J2X1_PCSAFTD(4,i,j,k) + gl%J2X1_PCSAFTD(5,i,j,k)) * part13
                            
                        end if
                    end do
                end do
            part1 = gl%MyPCSAFTD(k)**4 * gl%nPCSAFTD(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**3 * gl%mPCSAFT(k)**2)
            part3a = (-gl%J2X1_PCSAFTD(7,k,k,k) + 4.d0 * gl%J2X1_PCSAFTD(6,k,k,k) - 6.d0 * gl%J2X1_PCSAFTD(2,k,k,k)) * part1
            part3b = (-6.d0 * gl%J2X1_PCSAFTD(1,k,k,k) + 4.d0 * gl%J2X1_PCSAFTD(4,k,k,k) - gl%J2X1_PCSAFTD(5,k,k,k)) * part1
            part2 = 2.d0 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTD(7,k) =   -piPCSAFT * DENS * 1.d-38 * (part3a + part2 * (gl%molfractions(k) * Sum1 + Sum2) - 2.d0 * part2 * (gl%molfractions(k) * Sum3 + Sum4) &
                & + 2.d0 * part2 * (gl%molfractions(k) * Sum5 + Sum6) + Sum13 + part3b + part2 * (gl%molfractions(k) * Sum7 + Sum8) - 2.d0 * part2 * (gl%molfractions(k) * Sum9 + Sum10) &
                & + 2.d0 * part2 * (gl%molfractions(k)*Sum11 + Sum12) + Sum14)
             !1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
    
    ! 8: 3RD DERIVATIVE OF a_2dd WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^3
    ! requires j2_dd
    if (GETDERADD(8) .eq. 1) then
       do k = 1, gl%ncomp
            sum = 0.d0
            sum1 = 0.d0
            sum2 = 0.d0
            sum3 = 0.d0
            sum4 = 0.d0
            sum5 = 0.d0
            sum6 = 0.d0
                do j = 1, gl%ncomp
                    sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                    sum = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (T*kbol*sigma_kj**3 * gl%mPCSAFT(j))
                    sum1 = sum1 + gl%J2X1_PCSAFTD(3,j,k,k) * sum
                    sum2 = sum2 + gl%J2_PCSAFTD(3,j,k) * sum
                    sum3 = sum3 + gl%J2X1_PCSAFTD(8,j,k,k) * sum
                    sum4 = sum4 + gl%J2_PCSAFTD(8,j,k) * sum
                    do i = 1, gl%ncomp
                        if (i /= k .and. j /= k) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                            part1 = gl%MyPCSAFTD(i)**2 * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(i) * gl%nPCSAFTD(j) * gl%molfractions(i) * gl%molfractions(j) / (T**2 * kbol**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(j))
                            sum6 = sum6 + gl%J2X1_PCSAFTD(3,i,j,k) * part1
                            sum5 = sum5 + gl%J2X1_PCSAFTD(8,i,j,k) * part1
                        end if
                    end do
                end do
            part1 = gl%MyPCSAFTD(k)**4*gl%nPCSAFTD(k)**2*gl%molfractions(k)**2/(T**2*kbol**2*gl%sigPCSAFT(k)**3*gl%mPCSAFT(k)**2)    
            part3a = - 3.d0*gl%J2X1_PCSAFTD(3,k,k,k) * part1
            part3b = - gl%J2X1_PCSAFTD(8,k,k,k) * part1
            part2 = 2.d0 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTD(8,k) = - piPCSAFT * 1.d-38 * DENS * (part3a + 3.d0 * part2 * (gl%molfractions(k) * Sum1 + Sum2) &
                & + part3b + part2 * (gl%molfractions(k) * Sum3 + Sum4) + Sum5 + 3.d0 * Sum6)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
   
    
    
    ! 9: 3RD DERIVATIVE OF a_hc WITH RESPECT TO T, MULTIPLIED BY T^3
    ! requires 
    if (GETDERADD(9) .eq. 1) then
          do k = 1, gl%ncomp
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
                sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
                sum = gl%molfractions(j) * gl%MyPCSAFTD(j)**2 * gl%nPCSAFTD(j) / (T*kbol*sigma_kj**3 * gl%mPCSAFT(j))
                sum1 = sum1 + (- 6.d0 * gl%J2X1_PCSAFTD(1,j,k,k) + 6.d0 * gl%J2X1_PCSAFTD(4,j,k,k) - 3.d0 * gl%J2X1_PCSAFTD(5,j,k,k) + gl%J2X1_PCSAFTD(9,j,k,k)) * sum
                sum2 = sum2 + (- 6.d0 * gl%J2_PCSAFTD(1,j,k) + 6.d0 * gl%J2_PCSAFTD(4,j,k) - 3.d0 * gl%J2_PCSAFTD(5,j,k) + gl%J2_PCSAFTD(9,j,k)) * sum
                sum3 = sum3 + (2.d0 * gl%J2X1_PCSAFTD(1,j,k,k) - 2.d0 * gl%J2X1_PCSAFTD(4,j,k,k) + gl%J2X1_PCSAFTD(5,j,k,k)) * sum
                sum4 = sum4 + (2.d0 * gl%J2_PCSAFTD(1,j,k) - 2.d0 * gl%J2_PCSAFTD(4,j,k) + gl%J2_PCSAFTD(5,j,k)) * sum
                sum5 = sum5 + (- gl%J2X1_PCSAFTD(1,j,k,k) + gl%J2X1_PCSAFTD(4,j,k,k)) * sum
                sum6 = sum6 + (- gl%J2_PCSAFTD(1,j,k) + gl%J2_PCSAFTD(4,j,k)) * sum
                sum7 = sum7 + gl%J2X1_PCSAFTD(1,j,k,k) * sum
                sum8 = sum8 + gl%J2_PCSAFTD(1,j,k) * sum
                do i = 1, gl%ncomp
                    if (i /= k .and. j /= k) then 
                        sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                        sum9 = sum9 + (- 24.d0 * gl%J2X1_PCSAFTD(1,i,j,k) + 18.d0 * gl%J2X1_PCSAFTD(4,i,j,k) - 6.d0 * gl%J2X1_PCSAFTD(5,i,j,k) + gl%J2X1_PCSAFTD(9,i,j,k)) &
                            & *gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                            & /(T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                    end if
                end do
            end do
            part3 = (24.d0 * gl%J2X1_PCSAFTD(1,k,k,k) - 18.d0 * gl%J2X1_PCSAFTD(4,k,k,k) + 6.d0 * gl%J2X1_PCSAFTD(5,k,k,k) - gl%J2X1_PCSAFTD(9,k,k,k)) * gl%MyPCSAFTD(k)**4 * gl%nPCSAFTD(k)**2 * gl%molfractions(k)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(k)**3 * gl%mPCSAFT(k)**2) 
            part2 = 2.d0 * gl%MyPCSAFTD(k)**2 * gl%nPCSAFTD(k) / (T * kbol * gl%mPCSAFT(k))
            
            gl%A2X1_PCSAFTD(9,k) =  -piPCSAFT * DENS * 1.d-38 * (part3 + part2 * (gl%molfractions(k) * Sum1 + Sum2) &
                & - 3.d0 * part2 * (gl%molfractions(k)*Sum3 + Sum4) + 6.d0 * part2 *(gl%molfractions(k)*Sum5 + Sum6) &
                & - 6.d0 * part2 * (gl%molfractions(k)*Sum7 + Sum8) + Sum9)
            ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
        end do
    end if
    
    ! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    ! requires j_2qq
    if (GETDERADD(10) .eq. 1) then
        do k = 1, gl%ncomp
            help = 0.d0
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
            sigma_kj = (gl%sigPCSAFT(k)+gl%sigPCSAFT(j))*0.5d0
            help = gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(j)*gl%molfractions(j)/(T*kbol*sigma_kj**3*gl%mPCSAFT(j))
            sum1 = sum1 + (gl%J2X1_PCSAFTD(6,j,k,k) - gl%J2X1_PCSAFTD(2,j,k,k)) * help
            sum2 = sum2 + (gl%J2_PCSAFTD(6,j,k) - gl%J2_PCSAFTD(2,j,k)) * help
            sum3 = sum3 + gl%J2X1_PCSAFTD(2,j,k,k) * help
            sum4 = sum4 + gl%J2_PCSAFTD(2,j,k) * help
            sum5 = sum5 + (gl%J2X1_PCSAFTD(10,j,k,k) - gl%J2X1_PCSAFTD(3,j,k,k)) * help
            sum6 = sum6 + (gl%J2_PCSAFTD(10,j,k) - gl%J2_PCSAFTD(3,j,k)) * help
            sum7 = sum7 + gl%J2X1_PCSAFTD(3,j,k,k) * help
            sum8 = sum8 + gl%J2_PCSAFTD(3,j,k) * help
                do i = 1, gl%ncomp
                
                    if (i /= k .and. j /= k) then 
                    sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                    sum9 = sum9 + (gl%J2X1_PCSAFTD(10,i,j,k) - 2.d0 * gl%J2X1_PCSAFTD(3,i,j,k))*gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                        & / (T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))  
                    sum10 = sum10 + (gl%J2X1_PCSAFTD(6,i,j,k) - 2.d0 * gl%J2X1_PCSAFTD(2,i,j,k))*gl%MyPCSAFTD(i)**2*gl%MyPCSAFTD(j)**2*gl%nPCSAFTD(i)*gl%nPCSAFTD(j)*gl%molfractions(i)*gl%molfractions(j) &
                        & /(T**2*kbol**2*sigma_ij**3*gl%mPCSAFT(i)*gl%mPCSAFT(j))
                    end if
                end do
            end do
            part3a = (- 2.d0 * gl%J2X1_PCSAFTD(6,k,k,k) + 4.d0 * gl%J2X1_PCSAFTD(2,k,k,k)) *gl%MyPCSAFTD(k)**4*gl%nPCSAFTD(k)**2*gl%molfractions(k)**2 &
                        & /(T**2*kbol**2*gl%sigPCSAFT(k)**3*gl%mPCSAFT(k)**2) 
            part3b = (- gl%J2X1_PCSAFTD(10,k,k,k) + 2.d0 * gl%J2X1_PCSAFTD(3,k,k,k)) *gl%MyPCSAFTD(k)**4*gl%nPCSAFTD(k)**2*gl%molfractions(k)**2 &
                        & /(T**2*kbol**2*gl%sigPCSAFT(k)**3*gl%mPCSAFT(k)**2)
                        
            help2 = 2.d0 * gl%MyPCSAFTD(k)**2*gl%nPCSAFTD(k) / (T*kbol*gl%mPCSAFT(k))
    
            gl%A2X1_PCSAFTD(10,k) =  -piPCSAFT * DENS * 1.d-38 * (part3a + 2.d0 * help2 * (gl%molfractions(k) * Sum1 + Sum2) &
                                & - 2.d0 * help2 *(gl%molfractions(k)*Sum3 + Sum4) &
                                & + part3b + help2 * (gl%molfractions(k)*Sum5 + Sum6) &
                                & - help2 * (gl%molfractions(k)*Sum7 + Sum8) &
                                & + Sum9 + 2.d0*Sum10)
            
            ! 1.d-38: comes from units in dimensionless Myi*^2 = 1.d-19  
        end do
    end if
    !! 10: 3RD MIXED DERIVATIVE OF g_ij WITH RESPECT TO T, D, AND D, MULTIPLIED BY T*D*D
    !! requires j_2dd
    !if (GETDERADD(10) .eq. 1) then
    !    do k = 1, ncomp
    !   
    !        sum1 = 0.d0
    !        
    !        sum2 = 0.d0
    !        sum3 = 0.d0
    !        sum4 = 0.d0
    !        sum5 = 0.d0
    !        sum6 = 0.d0
    !        sum7 = 0.d0
    !        sum8 = 0.d0
    !        sum9 = 0.d0
    !        sum10 = 0.d0
    !            do j = 1, ncomp
    !                sigma_kj = (sigPCSAFT(k)+sigPCSAFT(j))*0.5d0
    !                sum = molfractions(j) * MyPCSAFTD(j)**2 * nPCSAFTD(j) / (T*kbol*sigma_kj**3 * mPCSAFT(j))
    !                sum1 = sum1 + (J2X1_PCSAFTD(6,j,k,k) - J2X1_PCSAFTD(2,j,k,k)) * sum
    !                sum2 = sum2 + (J2_PCSAFTD(6,j,k) - J2_PCSAFTD(2,j,k)) * sum
    !                sum3 = sum3 + J2X1_PCSAFTD(2,j,k,k) * sum
    !                sum4 = sum4 + J2_PCSAFTD(2,j,k) * sum
    !                sum5 = sum5 + (J2X1_PCSAFTD(10,j,k,k) - J2X1_PCSAFTD(3,j,k,k)) * sum
    !                sum6 = sum6 + (J2_PCSAFTD(10,j,k) - J2_PCSAFTD(3,j,k)) * sum
    !                sum7 = sum7 + J2X1_PCSAFTD(3,j,k,k) * sum
    !                sum8 = sum8 + J2_PCSAFTD(3,j,k) * sum
    !                do i = 1, ncomp
    !                    if (i /= k .and. j /= k) then 
    !                        sigma_ij = (sigPCSAFT(i)+sigPCSAFT(j))*0.5d0
    !                        part1 = MyPCSAFTD(i)**2 * MyPCSAFTD(j)**2 * nPCSAFTD(i) * nPCSAFTD(j) * molfractions(i) * molfractions(j) / (T**2 * kbol**2 * sigma_ij**3 * mPCSAFT(i) * mPCSAFT(j))
    !                        sum9 = sum9 + (J2X1_PCSAFTD(10,i,j,k) - 2.d0 * J2X1_PCSAFTD(3,i,j,k)) * part1
    !                        sum10 = sum10 + (J2X1_PCSAFTD(6,i,j,k) - 2.d0 * J2X1_PCSAFTD(2,i,j,k)) * part1
    !                    end if
    !                end do
    !            end do
    !        part3 = MyPCSAFTD(k)**4 * nPCSAFTD(k)**2 * molfractions(k)**2 / (T**2 * kbol**2 * sigPCSAFT(k)**3 * mPCSAFT(k)**2)        
    !        part3a = (- 2.d0 * J2X1_PCSAFTD(6,k,k,k) + 4.d0 * J2X1_PCSAFTD(2,k,k,k)) * part3
    !        part3b = (- J2X1_PCSAFTD(10,k,k,k) + 2.d0 * J2X1_PCSAFTD(3,k,k,k)) * part3
    !        part2 = 2.d0 * MyPCSAFTD(k)**2 * nPCSAFTD(k) / (T * kbol * mPCSAFT(k))
    !        
    !        A2X1_PCSAFTD(10,k) = - piPCSAFT * 1.d-38 * DENS * (&
    !        & part3a &
    !        & + 2.d0 * part2 * (molfractions(k) * Sum1 + Sum2) &
    !        & - 2.d0 * part2 * (molfractions(k) * Sum3 + Sum4) &
    !        & + part3b &
    !        & + part2 * (molfractions(k)*Sum5 + Sum6) &
    !        & - part2 * (molfractions(k) * Sum7 + Sum8) &
    !        & + Sum9 + 2.d0 * Sum10)
    !        ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
    !    end do
    !end if
       
!DEC$ END IF
end subroutine A2X1DERIVS_D


subroutine A2X2DERIVS_D(gl,T,DENS,GETDERADD)
    ! a_2qq: second-order perturbation term of the quadrupolar contribution to the residual Helmholtz free energy
    ! defined by eq. A.9 in Gross, Sadowski 2005:
    ! a_2qq =-pi*(3/4)**2*rho*sum_j*sum_i(x_i*x_j*epsk_ii/T*espk_jj/T*sig_ii**5*sig_jj**5/sig_ij**7*n_Qi*n_Qj*Q_i**2*Q_j**2*J_2ij
    ! dependent on D and T and x_i







implicit none

    type(type_gl) :: gl


    ! I. Declarations
    !inputs
    integer, dimension (nderivs), intent (in) :: GETDERADD
    double precision, intent (in) :: T, DENS
    !output: aqq_PCSAFT (module variable)
    !working variables
    double precision :: sigma_ij, sigma_kj
    double precision :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, sum14, part13, part3, part3a, part3b, part3c, sum, help, help3, part1, part2
    integer :: i, j, xi, k, o, u
    double precision :: DP, DP2, DM, DM2, TP, TP2, TM, TM2
    double precision :: DPx, DP2x, DMx, DM2x, TPx, TP2x, TMx, TM2x
    double precision :: DPxx, DP2xx, DMxx, DM2xx, TPxx, TP2xx, TMxx, TM2xx
    double precision :: DELTA, DEL, DELT
    double precision, dimension(gl%ncomp,gl%ncomp) :: A2X2, A2X2DP, A2X2DP2, A2X2DM, A2X2DM2, A2X2TP, A2X2TP2, A2X2TM, A2X2TM2, A2X2PP, A2X2PM, A2X2MP, A2X2MM

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
    
    DELT = 5.d-4 
    
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
    if (GETDERADD(1) .eq. 1) then
        
        call base_A2X2_D(gl,T,DENS,A2X2)
        
        gl%A2X2_PCSAFTD(1,:,:) = A2X2

    end if
    
    
    !  2: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D
    ! requires j_2qq
    if (GETDERADD(2) .eq. 1) then   
               
        call base_A2X2_D(gl,T,DP,A2X2DP)
        call base_A2X2_D(gl,T,DP2,A2X2DP2)
        call base_A2X2_D(gl,T,DM,A2X2DM)
        call base_A2X2_D(gl,T,DM2,A2X2DM2)
        
        gl%A2X2_PCSAFTD(2,:,:) = (8.d0*A2X2DP - 8.d0*A2X2DM - A2X2DP2 + A2X2DM2) / (12.d0 * DELTA) 
    

    end if
    

   ! 3: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO D AT CONSTANT T, MULTIPLIED BY D^2
   ! requires j_2qq
    if (GETDERADD(3) .eq. 1) then
        
        call base_A2X2_D(gl,T,DENS,A2X2)    
        call base_A2X2_D(gl,T,DPxx,A2X2DP)
        call base_A2X2_D(gl,T,DP2xx,A2X2DP2)
        call base_A2X2_D(gl,T,DMxx,A2X2DM)
        call base_A2X2_D(gl,T,DM2xx,A2X2DM2)
        
        gl%A2X2_PCSAFTD(3,:,:) = (-A2X2DM2 + 16.d0*A2X2DM - 30.d0*A2X2 + 16.d0*A2X2DP - A2X2DP2) / (12.d0*DELT**2) 
    
    end if
    
   
     
    ! 4: 1ST DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T
    ! requires j_2qq
    if (GETDERADD(4) .eq. 1) then
        call base_A2X2_D(gl,TP,DENS,A2X2TP)
        call base_A2X2_D(gl,TP2,DENS,A2X2TP2)
        call base_A2X2_D(gl,TM,DENS,A2X2TM)
        call base_A2X2_D(gl,TM2,DENS,A2X2TM2)
        
        gl%A2X2_PCSAFTD(4,:,:) = (8.d0*A2X2TP - 8.d0*A2X2TM - A2X2TP2 + A2X2TM2) / (12.d0 * DELTA) 
    end if
    

    ! 5: 2ND DERIVATIVE OF a_2qq WITH RESPECT TO T AT CONSTANT D, MULTIPLIED BY T^2
    !requires j_2qq
    if (GETDERADD(5) .eq. 1) then
    
        call base_A2X2_D(gl,T,DENS,A2X2)    
        call base_A2X2_D(gl,TPx,DENS,A2X2TP)
        call base_A2X2_D(gl,TP2x,DENS,A2X2TP2)
        call base_A2X2_D(gl,TMx,DENS,A2X2TM)
        call base_A2X2_D(gl,TM2x,DENS,A2X2TM2)
        
        gl%A2X2_PCSAFTD(5,:,:) = (-A2X2TM2 + 16.d0*A2X2TM - 30.d0*A2X2 + 16.d0*A2X2TP - A2X2TP2) / (12.d0*DEL**2) !https://en.wikipedia.org/wiki/Five-point_stencil
    
    end if
       
     
    ! 6: 2ND MIXED DERIVATIVE OF a_2qq WITH RESPECT TO D AND T, MULTIPLIED BY T*D
    ! requires j_2qq
    if (GETDERADD(6) .eq. 1) then
    
        call base_A2X2_D(gl,TPxx,DPxx,A2X2PP)
        call base_A2X2_D(gl,TMxx,DPxx,A2X2MP)
        call base_A2X2_D(gl,TMxx,DMxx,A2X2MM)
        call base_A2X2_D(gl,TPxx,DMxx,A2X2PM)
        
        gl%A2X2_PCSAFTD(6,:,:) = (A2X2PP - A2X2PM - A2X2MP + A2X2MM) / (4.d0*DELT**2) 
    end if
    
!DEC$ END IF
end subroutine A2X2DERIVS_D

 
subroutine base_A2X2(gl,T,D,A2X2)





    

implicit none

    type(type_gl) :: gl

    double precision :: sigma_ij, sigma_kj
    double precision :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum, part1
    integer :: i, j, xi, k, o, u
    double precision :: T,D
    double precision, dimension(gl%ncomp,gl%ncomp) :: A2X2
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
    call calculate_PCSAFT_functionparts_x2(gl,T,D,getprevious) !Danger: Hier darf A2X2 nicht selber aufgerufen werden --> Endlosschleife... 
    call J2X2DERIVS(gl,T,D,getprevious)
    
        do u = 1, gl%ncomp
            do o = 1, gl%ncomp 
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
            
                sigma_ij = (gl%sigPCSAFT(o)+gl%sigPCSAFT(u))*0.5d0
                sum = gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * gl%QPCSAFTQ(u)**2 * gl%nPCSAFTQ(u)  / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(o) * gl%mPCSAFT(u))
                sum10 = sum * (gl%molfractions(u) * gl%J2X1_PCSAFTQ(1,u,o,u) + gl%J2_PCSAFTQ(1,u,o))
                
                if (u == o) then
                    part1 = (- 2.d0 * gl%molfractions(o) * gl%J2X1_PCSAFTQ(1,o,o,o) - gl%molfractions(o)**2 * gl%J2X2_PCSAFTQ(1,o,o,o,o)) * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2  / (T**2 * kbol**2 * gl%sigPCSAFT(o)**7 * gl%mPCSAFT(o)**2)  
                   
                    sigma_ij = (gl%sigPCSAFT(o)+gl%sigPCSAFT(u))*0.5d0
                    sum = gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * gl%QPCSAFTQ(u)**2 * gl%nPCSAFTQ(u)  / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(o) * gl%mPCSAFT(u))
                    sum11 = sum * (gl%molfractions(u)**2 * gl%J2X2_PCSAFTQ(1,u,o,o,u) + 2.d0*gl%molfractions(u)*gl%J2X1_PCSAFTQ(1,u,o,o))
                else
                    part1 = - gl%molfractions(o)**2 * gl%J2X2_PCSAFTQ(1,o,o,o,u) * gl%QPCSAFTQ(o)**4 * gl%nPCSAFTQ(o)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(o)**7 * gl%mPCSAFT(o)**2)
                    
                    sigma_ij = (gl%sigPCSAFT(o)+gl%sigPCSAFT(u))*0.5d0
                    sum = gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * gl%QPCSAFTQ(u)**2 * gl%nPCSAFTQ(u)  / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(o) * gl%mPCSAFT(u))
                    sum11 = sum * gl%molfractions(o) * (gl%molfractions(u) * gl%J2X2_PCSAFTQ(1,u,o,o,u) + gl%J2X1_PCSAFTQ(1,u,o,o))
                end if
                
                do i = 1, gl%ncomp
                    if (i /= u) then 
                        sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(o))*0.5d0
                        sum1 = sum1 + gl%molfractions(i) * gl%J2X1_PCSAFTQ(1,i,o,u) * gl%nPCSAFTQ(i) * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * gl%QPCSAFTQ(i)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(o))
                    
                        if (u == o) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(o))*0.5d0
                            sum2 = sum2 + (gl%molfractions(i)*gl%J2X1_PCSAFTQ(1,i,o,o) + gl%molfractions(i)*gl%molfractions(o)*gl%J2X2_PCSAFTQ(1,i,o,o,u)) * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * gl%nPCSAFTQ(i)* gl%QPCSAFTQ(i)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(o))
                        else
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(o))*0.5d0
                            sum2 = sum2 + gl%molfractions(o) * gl%molfractions(i) * gl%J2X2_PCSAFTQ(1,i,o,o,u) * gl%QPCSAFTQ(o)**2 * gl%nPCSAFTQ(o) * gl%nPCSAFTQ(i)* gl%QPCSAFTQ(i)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(i) * gl%mPCSAFT(o))
                        end if
                    end if
                    
                    
                    do j = 1, gl%ncomp
                        if (i /= o .and. j /= o) then
                            if (i == u .and. j == u) then 
                                sigma_ij = (gl%sigPCSAFT(u)+gl%sigPCSAFT(u))*0.5d0
                                sum3 = sum3 + gl%molfractions(u)**2 * gl%J2X2_PCSAFTQ(1,u,u,o,u) * gl%nPCSAFTQ(u)*gl%nPCSAFTQ(u) * gl%QPCSAFTQ(u)**2 *gl%QPCSAFTQ(u)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(u)*gl%mPCSAFT(u))
                                sum4 = sum4 + 2.d0*gl%molfractions(u) * gl%J2X1_PCSAFTQ(1,u,u,o) * gl%nPCSAFTQ(u)*gl%nPCSAFTQ(u) * gl%QPCSAFTQ(u)**2 *gl%QPCSAFTQ(u)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(u)*gl%mPCSAFT(u))
                            elseif (i == u .and. j /= u) then 
                                sigma_ij = (gl%sigPCSAFT(u)+gl%sigPCSAFT(j))*0.5d0
                                sum5 = sum5 + gl%molfractions(j) * gl%J2X1_PCSAFTQ(1,u,j,o) * gl%nPCSAFTQ(u)*gl%nPCSAFTQ(j) * gl%QPCSAFTQ(u)**2 *gl%QPCSAFTQ(j)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(u)*gl%mPCSAFT(j))
                                sum6 = sum6 + gl%molfractions(u)*gl%molfractions(j) * gl%J2X2_PCSAFTQ(1,u,j,o,u) * gl%nPCSAFTQ(u)*gl%nPCSAFTQ(j) * gl%QPCSAFTQ(u)**2 *gl%QPCSAFTQ(j)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(u)*gl%mPCSAFT(j))
                            elseif (i /= u .and. j == u) then 
                                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(u))*0.5d0
                                sum8 = sum8 + gl%molfractions(i) * gl%J2X1_PCSAFTQ(1,u,i,o) * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(u) * gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(u)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(i)*gl%mPCSAFT(u))
                                sum9 = sum9 + gl%molfractions(u)*gl%molfractions(i) * gl%J2X2_PCSAFTQ(1,u,i,o,u) * gl%nPCSAFTQ(u)*gl%nPCSAFTQ(i) * gl%QPCSAFTQ(u)**2 *gl%QPCSAFTQ(i)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(u)*gl%mPCSAFT(i))
                            elseif (i /= u .and. j /= u) then 
                                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                                sum7 = sum7 + gl%molfractions(i)*gl%molfractions(j) * gl%J2X2_PCSAFTQ(1,i,j,o,u) * gl%nPCSAFTQ(i)*gl%nPCSAFTQ(j) * gl%QPCSAFTQ(i)**2 *gl%QPCSAFTQ(j)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(i)*gl%mPCSAFT(j))
                            end if
                        end if
                    end do
                end do
               
                A2X2(o,u) = - 0.5625d0 * piPCSAFT * D * 1.d-38 * (2.d0 * (Sum1 + Sum2) + 2.d0 * (sum10 + sum11) + part1 + Sum3 + Sum4 + Sum5 + sum6 + sum7 + sum8 + sum9)
                ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do 
        end do 
        
!DEC$ END IF
end subroutine


subroutine base_A2X2_D(gl,T,D,A2X2_D)





    

implicit none

    type(type_gl) :: gl

    double precision :: sigma_ij, sigma_kj
    double precision :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum, part1
    integer :: i, j, xi, k, o, u
    double precision :: T,D
    double precision, dimension(gl%ncomp,gl%ncomp) :: A2X2_D
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
    call calculate_PCSAFT_functionparts_x2(gl,T,D,getprevious) !Danger: Hier darf A2X2_D nicht selber aufgerufen werden --> Endlosschleife... 
    call J2X2DERIVS_D(gl,T,D,getprevious)
    
      
        do u = 1, gl%ncomp
            do o = 1, gl%ncomp 
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
            
                sigma_ij = (gl%sigPCSAFT(o)+gl%sigPCSAFT(u))*0.5d0
                sum = gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * gl%MyPCSAFTD(u)**2 * gl%nPCSAFTD(u)  / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(o) * gl%mPCSAFT(u))
                sum10 = sum * (gl%molfractions(u) * gl%J2X1_PCSAFTD(1,u,o,u) + gl%J2_PCSAFTD(1,u,o))
                
                if (u == o) then
                    part1 = (- 2.d0 * gl%molfractions(o) * gl%J2X1_PCSAFTD(1,o,o,o) - gl%molfractions(o)**2 * gl%J2X2_PCSAFTD(1,o,o,o,o)) * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2  / (T**2 * kbol**2 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2)  
                   
                    sigma_ij = (gl%sigPCSAFT(o)+gl%sigPCSAFT(u))*0.5d0
                    sum = gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * gl%MyPCSAFTD(u)**2 * gl%nPCSAFTD(u)  / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(o) * gl%mPCSAFT(u))
                    !sum11 = sum * (molfractions(u) * J2X2_PCSAFTD(1,u,o,o,u) + J2X1_PCSAFTD(1,u,o,o))
                    sum11 = sum * (gl%molfractions(u)**2 * gl%J2X2_PCSAFTD(1,u,o,o,u) + 2.d0*gl%molfractions(u)*gl%J2X1_PCSAFTD(1,u,o,o))
                else
                    part1 = - gl%molfractions(o)**2 * gl%J2X2_PCSAFTD(1,o,o,o,u) * gl%MyPCSAFTD(o)**4 * gl%nPCSAFTD(o)**2 / (T**2 * kbol**2 * gl%sigPCSAFT(o)**3 * gl%mPCSAFT(o)**2)
                    
                    sigma_ij = (gl%sigPCSAFT(o)+gl%sigPCSAFT(u))*0.5d0
                    sum = gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * gl%MyPCSAFTD(u)**2 * gl%nPCSAFTD(u)  / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(o) * gl%mPCSAFT(u))
                    sum11 = sum * gl%molfractions(o) * (gl%molfractions(u) * gl%J2X2_PCSAFTD(1,u,o,o,u) + gl%J2X1_PCSAFTD(1,u,o,o))
                end if
                
                do i = 1, gl%ncomp
                    if (i /= u) then 
                        sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(o))*0.5d0
                        sum1 = sum1 + gl%molfractions(i) * gl%J2X1_PCSAFTD(1,i,o,u) * gl%nPCSAFTD(i) * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * gl%MyPCSAFTD(i)**2 / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(o))
                    
                        if (u == o) then 
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(o))*0.5d0
                            sum2 = sum2 + (gl%molfractions(i)*gl%J2X1_PCSAFTD(1,i,o,o) + gl%molfractions(i)*gl%molfractions(o)*gl%J2X2_PCSAFTD(1,i,o,o,u)) * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * gl%nPCSAFTD(i)* gl%MyPCSAFTD(i)**2 / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(o))
                        else
                            sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(o))*0.5d0
                            sum2 = sum2 + gl%molfractions(o) * gl%molfractions(i) * gl%J2X2_PCSAFTD(1,i,o,o,u) * gl%MyPCSAFTD(o)**2 * gl%nPCSAFTD(o) * gl%nPCSAFTD(i)* gl%MyPCSAFTD(i)**2 / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(i) * gl%mPCSAFT(o))
                        end if
                    end if
                    
                    
                    do j = 1, gl%ncomp
                        if (i /= o .and. j /= o) then
                            if (i == u .and. j == u) then 
                                sigma_ij = (gl%sigPCSAFT(u)+gl%sigPCSAFT(u))*0.5d0
                                sum3 = sum3 + gl%molfractions(u)**2 * gl%J2X2_PCSAFTD(1,u,u,o,u) * gl%nPCSAFTD(u)*gl%nPCSAFTD(u) * gl%MyPCSAFTD(u)**2 *gl%MyPCSAFTD(u)**2 / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(u)*gl%mPCSAFT(u))
                                sum4 = sum4 + 2.d0*gl%molfractions(u) * gl%J2X1_PCSAFTD(1,u,u,o) * gl%nPCSAFTD(u)*gl%nPCSAFTD(u) * gl%MyPCSAFTD(u)**2 *gl%MyPCSAFTD(u)**2 / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(u)*gl%mPCSAFT(u))
                            elseif (i == u .and. j /= u) then 
                                sigma_ij = (gl%sigPCSAFT(u)+gl%sigPCSAFT(j))*0.5d0
                                sum5 = sum5 + gl%molfractions(j) * gl%J2X1_PCSAFTD(1,u,j,o) * gl%nPCSAFTD(u)*gl%nPCSAFTD(j) * gl%MyPCSAFTD(u)**2 *gl%MyPCSAFTD(j)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(u)*gl%mPCSAFT(j))
                                sum6 = sum6 + gl%molfractions(u)*gl%molfractions(j) * gl%J2X2_PCSAFTD(1,u,j,o,u) * gl%nPCSAFTD(u)*gl%nPCSAFTD(j) * gl%MyPCSAFTD(u)**2 *gl%MyPCSAFTD(j)**2 / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(u)*gl%mPCSAFT(j))
                            elseif (i /= u .and. j == u) then 
                                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(u))*0.5d0
                                sum8 = sum8 + gl%molfractions(i) * gl%J2X1_PCSAFTD(1,u,i,o) * gl%nPCSAFTD(i)*gl%nPCSAFTD(u) * gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(u)**2 / (kbol**2 * T**2 * sigma_ij**7 * gl%mPCSAFT(i)*gl%mPCSAFT(u))
                                sum9 = sum9 + gl%molfractions(u)*gl%molfractions(i) * gl%J2X2_PCSAFTD(1,u,i,o,u) * gl%nPCSAFTD(u)*gl%nPCSAFTD(i) * gl%MyPCSAFTD(u)**2 *gl%MyPCSAFTD(i)**2 / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(u)*gl%mPCSAFT(i))
                            elseif (i /= u .and. j /= u) then 
                                sigma_ij = (gl%sigPCSAFT(i)+gl%sigPCSAFT(j))*0.5d0
                                sum7 = sum7 + gl%molfractions(i)*gl%molfractions(j) * gl%J2X2_PCSAFTD(1,i,j,o,u) * gl%nPCSAFTD(i)*gl%nPCSAFTD(j) * gl%MyPCSAFTD(i)**2 *gl%MyPCSAFTD(j)**2 / (kbol**2 * T**2 * sigma_ij**3 * gl%mPCSAFT(i)*gl%mPCSAFT(j))
                            end if
                        end if
                    end do
                end do
               
                A2X2_D(o,u) = - piPCSAFT * D * 1.d-38 * (2.d0 * (Sum1 + Sum2) + 2.d0 * (sum10 + sum11) + part1 + Sum3 + Sum4 + Sum5 + sum6 + sum7 + sum8 + sum9)
                ! 1.d-38: comes from units in dimensionless Qi*^2 = 1.d-19  
            end do 
        end do 

        
!DEC$ END IF
end subroutine




    end module pc_saft_A2X_derivs_module
