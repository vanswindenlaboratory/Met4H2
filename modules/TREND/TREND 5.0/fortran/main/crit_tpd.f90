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

    ! module for file crit_tpd.f90
    module crit_tpd_module
    !global use inclusion
    use module_all_types
    use calc_functions


    contains




    subroutine find_crit_tpd(gl,tc_est_inp, rhoc_est_inp, tc_eos, rhoc_eos, pc_eos, iter, errval, nrsubst)
    !M. Thol, 2015/02/18
    !Subroutine to find the saddle point of the EOS, which is the critical point






    implicit none

    type(type_gl) :: gl
    double precision:: tc_est, rhoc_est,tc_est_inp, rhoc_est_inp, tc_est_new, rhoc_est_new
    double precision:: tc_eos, rhoc_eos, pc_eos!, p_calc
    integer:: nrsubst, nEQN, k, i

    integer :: errval, iter

    double precision, dimension(2) :: var, var_est
    double precision, dimension(2) :: EQN
    double precision, dimension(2) :: Delta_D
    double precision, dimension(2,2) :: JacMatrix_crit

    !Exit criterion
    double precision:: eps_eqn, eps_eqn2, eps_eqn3,eps_eqn4

    integer :: irow,icol
    double precision :: det_j, mat_j_inv(2,2), mat(2,2), mat_diag(2,2)

    !define the standard exit criterion
    eps_eqn = 1.D-10
    !define an alternate (softer) exit criterion in case the first criterion cannot be reached out of numeric reasons
    eps_eqn2 = 1.D-8

    mat_diag= reshape((/1.d0,0.d0,0.d0,1.d0/),(/2,2/))

    nEQN = 2
    EQN = 0.d0
    errval = 0
    Delta_D = 0.d0
    mat_j_inv = 0.d0

    tc_est = tc_est_inp
    rhoc_est = rhoc_est_inp

    do iter = 1, 50 ! a maximum of 50 iterations allowed

        var_est = 0.d0
        var_est(1) = tc_est
        var_est(2) = rhoc_est

        !SysOfEqs: calculate system of equations (dpdd=0, d2pdd2=0 => saddle point at the critical point)
        call SysOfEqs_crit(gl,tc_est, rhoc_est, EQN, errval,nrsubst)


       if (errval /= 0) then
            return
        end if

        ! end calculation if a certain accuracy is reached
        if (maxval(dabs(EQN)) < eps_eqn) then
            if (tc_est < 0.d0) then
                errval = -8889
            else
                tc_eos = tc_est
                rhoc_eos = rhoc_est
                pc_eos = p_calc(gl,tc_est, rhoc_est, nrsubst)
            end if
            return
        end if


        Delta_D = -EQN

        ! calculation of the derivatives needed for the solution of the system of equations
        call Jacobi_crit (gl,tc_est, rhoc_est, JacMatrix_crit, errval,nrsubst)

        if (errval /= 0) then
            return
        end if

        det_j = DET2X2_crit(gl,JacMatrix_crit)

        do icol=1,2
            do irow=1,2
                mat(:,:)=JacMatrix_crit(:,:)
                mat(:,irow)=mat_diag(:,icol)
                mat_j_inv(irow,icol)=DET2X2_crit(gl,mat)/det_j
            enddo
        enddo

        Delta_D=matmul(Delta_D,mat_j_inv)

        ! set the new values
        tc_est_new = tc_est + Delta_D(1)
        rhoc_est_new = rhoc_est + Delta_D(2)

        var = 0.d0
        var(1) = tc_est_new
        var(2) = rhoc_est_new

        tc_est = var(1)
        rhoc_est = var(2)

    end do
    errval = -8889
    tc_eos = tc_est
    rhoc_eos = rhoc_est
    pc_eos = p_calc(gl,tc_est, rhoc_est, nrsubst)

    end subroutine


    subroutine SysOfEqs_crit(gl,tc_est, rhoc_est, EQN, errval,nrsubst)

    ! SysOfEqs: (dpdd=0, d2pdd2=0 => saddle point at the critical point)




    implicit none

    type(type_gl) :: gl


    double precision :: tc_est, rhoc_est
    integer :: errval

    integer :: nrsubst
    double precision, dimension(2) :: EQN
    double precision :: dpdd, d2pdd2!, dpdd_calc, d2pdd2_calc

    dpdd = dpdd_calc(gl,tc_est, rhoc_est,nrsubst)
    d2pdd2 = d2pdd2_calc(gl,tc_est, rhoc_est,nrsubst)

    if (isnan(dpdd) .or. isnan(d2pdd2)) errval = -8889

    EQN(1) = dpdd
    EQN(2) = d2pdd2


    end subroutine SysOfEqs_crit




    subroutine Jacobi_crit (gl,tc_est, rhoc_est, JacMatrix_crit, errval,nrsubst)








    implicit none

    type(type_gl) :: gl


    double precision :: tc_est, rhoc_est
    integer :: errval

    integer :: nrsubst
    double precision, dimension(2,2) :: JacMatrix_crit
    !double precision :: D2PDD2_CALC, D2PDDT_CALC, D3PDDDT_CALC, D3PDD3_CALC
    double precision :: D2PDD2, D2PDDT, D3PDDDT, D3PDD3

    double precision:: diff
    double precision:: D2PDD2_p, D2PDD2_m
    double precision:: tc_est_p, tc_est_m, rhoc_est_p, rhoc_est_m

    JacMatrix_crit = 0.d0
    errval = 0
    D2PDD2 = 0.d0
    D2PDDT = 0.d0
    D3PDDDT = 0.d0
    D3PDD3 = 0.d0

    diff = 1.d-4

    !analytic derivatives
    D2PDD2 = D2PDD2_CALC (gl,tc_est, rhoc_est, nrsubst)
    D2PDDT = D2PDDT_CALC(gl,tc_est, rhoc_est, nrsubst)/gl%factortrans     !For real fluids, this derivative must be devided by 10^6 so that it matches the other derivatives, which are given in MPa instead of Pa

    !numerical derivatives
    tc_est_p = tc_est * (1.d0 + diff)
    D2PDD2_p = D2PDD2_CALC(gl,tc_est_p, rhoc_est, nrsubst)
    tc_est_m = tc_est * (1.d0 - diff)
    D2PDD2_m = D2PDD2_CALC(gl,tc_est_m, rhoc_est, nrsubst)
    D3PDDDT = (D2PDD2_p-D2PDD2_m)/(2.d0*diff*tc_est)

    rhoc_est_p = rhoc_est * (1.d0 + diff)
    D2PDD2_p = D2PDD2_CALC(gl,tc_est, rhoc_est_p, nrsubst)
    rhoc_est_m = rhoc_est * (1.d0 - diff)
    D2PDD2_m = D2PDD2_CALC(gl,tc_est, rhoc_est_m, nrsubst)
    D3PDD3 = (D2PDD2_p-D2PDD2_m)/(2.d0*diff*rhoc_est)

    ! set Jacobi Matrix with the derivatives needed
    JacMatrix_crit(1,1) = D2PDDT
    JacMatrix_crit(2,1) = D2PDD2
    JacMatrix_crit(1,2) = D3PDDDT
    JacMatrix_crit(2,2) = D3PDD3

    end subroutine Jacobi_crit


    double precision function DET2X2_crit(gl,mat)

    implicit none

    type(type_gl) :: gl


    double precision:: mat(2,2)

    DET2X2_crit=mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2)

    end function DET2X2_crit



    end module crit_tpd_module
