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


    !Mutke
    module ljf_tansport

    implicit none


    double precision, dimension(7), parameter :: a = (/0D0,3.108096D-1,-1.712113D-1,-7.158049D-1,2.486777D0,-1.783167D0,3.944048D-1/)
    double precision :: Vis_calc
    double precision :: rhostar, T_star, sigmastarB, rhostar_LJ, fstar, zeta, gsigma, denom_
    double precision,  parameter :: pi = 3.14159265358979d0

    double precision, dimension(6), parameter :: omega = (/0d0, 2.8745d0, -2.0265d0, 0.9158d0, -0.1960d0, 0.0160d0/)
    double precision :: summer, summer_sum


    contains


    !VS3 model for lennard jones fluid
    double precision function Ruckenstein_IECR_1997(T,D)

    double precision :: T,D

    T_star = T
    rhostar = D


    ! def get_rhostarDstar(self, *, T_star, rhostar):
    sigmastarB = 1.1532d0*(1d0 + (T_star/0.527d0)**(1d0/2d0))**(-1d0/6d0)
    ! Ruckenstein and Liu define a different density, which we call here rhostar_LJ
    ! See Eq. 40.  "Normal" rho^* is rho*sigma_{LJ}^3
    rhostar_LJ = rhostar*sigmastarB**3d0
    fstar = 1d0 + 0.94605d0*rhostar_LJ**1.5d0 + 1.4022d0*rhostar_LJ**3d0 - 5.6898d0*rhostar_LJ**5d0 + 2.6626d0*rhostar_LJ**7d0
    zeta = rhostar_LJ*(pi/6d0)
    gsigma = (1d0-0.5d0*zeta)/(1d0-zeta)**3d0
    denom_ = 8d0*pi**0.5d0/3d0*sigmastarB**2d0*(gsigma/fstar + 0.4d0/T_star**1.5d0)

    Vis_calc = T_star**0.5d0/denom_

    Ruckenstein_IECR_1997 = Vis_calc

    end function


    !VS4 model for lennard jones fluid
    double precision function Rowley_IJT_1997(T,D)

    integer :: i,j
    double precision :: T,D
    T_star = T
    rhostar = D

    ! def get_etastar_dilute(self, T_star):
    ! """ Rowley, Eq. 10 (b) """ ! padded with one zero for index consistency with paper
    summer = 0d0
    summer_sum = 0d0
    do j = 1, 6
        summer = omega(j)*T_star**(j-1d0)
        summer_sum = summer_sum + summer
    end do
    Vis_calc = 5d0/16d0*(T_star/pi)**(0.5d0)/(summer_sum)

    Rowley_IJT_1997 = Vis_calc

    end function

    !VS5 model for lennard jones fluid
    !double precision function Rowley_IJT_1997(T,D)
    !
    !integer :: i,j
    !double precision :: T,D
    !T_star = T
    !rhostar = D
    !
    !! def get_etastar_dilute(self, T_star):
    !! """ Rowley, Eq. 10 (b) """ ! padded with one zero for index consistency with paper
    !summer = 0d0
    !summer_sum = 0d0
    !do j = 1, 6
    !    summer = omega(j)*T_star**(j-1d0)
    !    summer_sum = summer_sum + summer
    !end do
    !Vis_calc = 5d0/16d0*(T_star/pi)**(0.5d0)/(summer_sum)
    !
    !Rowley_IJT_1997 = Vis_calc
    !
    !end function

    !!VS4 model for lennard jones fluid
    !double precision function Rowley_IJT_1997(T,D)
    !
    !
    !Rowley_IJT_1997 = 1d0
    !
    !end function



    end module






    ! module for file transport.f90
    module transport_module
    !global use inclusion
    use module_all_types
    use calc_functions
    use setup_module
    use ancillary_equations_mix_module
    contains




    !----------------------------------------------------------------------------------
    !Function for calculating viscosity [Pa-s*10^6]
    !----------------------------------------------------------------------------------
    ! M. Thol, 2013

    DOUBLE PRECISION FUNCTION VISDYN_CALC(gl,T,D,NRSUBST)

    use ljf_tansport

    implicit none

    type(type_gl) :: gl


    !---------------------------------------------------------------------------------------------------------
    !Declaration of variables used in DOUBLE PRECISION FUNCTION VISDYN_CALC
    DOUBLE PRECISION:: mue, mue0, mue1,mue2 !Declaration for VS0 variables
    DOUBLE PRECISION:: nue, nue0, nuer, P,etagas,etah  !
    DOUBLE PRECISION:: Xi_mue2,delchi_mue2, Y_mue2, lw_mue2, w_mue2
    DOUBLE PRECISION:: tdiml, T, D,dpdd_state,dpdd_ref, dddp_state,dddp_ref, psid_mue2, tref_mue
    DOUBLE PRECISION:: rhodiml, tau, delta
    DOUBLE PRECISION:: Hi, sum1, sum_new, sum1_mue1, sum2_mue1, Tmelt_eta
    DOUBLE PRECISION:: theta,beta_lc,del0
    double PRECISION::  eta0_lc, eta1, deleta, deletacrit,del,betastar,sp
    DOUBLE PRECISION:: norp,dorp,nexpeta,dexpeta, rho, etal, etadg
    DOUBLE PRECISION:: etadash, etacrit, DSI, TT     !variables for ethylene
    DOUBLE PRECISION:: fcross, etaliq,ciLJ,ciDEL,ciSM,B_eta, C_eta, B_etar,C_etar,etagasred  !variables for methanol
    DOUBLE PRECISION:: Tr, Dr, sigHS, bvol, ksi, gHS, Tstar, sigceta, sigHSr                 !variables for methanol
    DOUBLE PRECISION:: etaexess, lnsstar, sstar, rhoscmol, rhoredH2,TredH2      !variables for hydrogen
    DOUBLE PRECISION,DIMENSION(5):: aiH2      !hydrogen
    DOUBLE PRECISION,DIMENSION(7):: biH2      !hydrogen
    DOUBLE PRECISION,DIMENSION(6):: ciH2      !hydrogen
    DOUBLE PRECISION:: collsum, ci, Rmix
    INTEGER:: n, m, i, j, k, nrsubst,l,o,pcount


    !VS2
    !parameters for eq. (18)
    double precision:: eta_VS2
    double precision:: eta0_VS2
    double precision:: eta1_VS2
    double precision:: eta2_VS2
    !parameters for Eq. (20)
    double precision:: omega_VS2
    !parameters for Eq. (22) - (25)
    double precision:: Frhot_VS2
    double precision:: Gt_VS2
    double precision:: Hrho_VS2

    !VS4
    DOUBLE PRECISION:: tau3
    DOUBLE PRECISION:: ka,kaa,kr,krr,ki,kii, kttt, kaaa,eta4rs
    DOUBLE PRECISION:: phi1,phi2
    DOUBLE PRECISION:: pid, pa, pr, DPDT
    logical :: within_limits

    !VS5
    !parameters for Eq. (1)
    double precision:: eta0_VS5
    double precision:: omegastar
    !parameters for eq.(2)
    double precision:: A_VS5
    double precision:: B_VS5
    double precision:: C_VS5
    double precision:: D_VS5
    double precision:: E_VS5
    double precision:: F_VS5
    double precision:: G_VS5
    double precision:: H_VS5
    double precision:: S_VS5
    double precision:: W_VS5
    !parameter for Eq. (6)
    double precision:: vc_VS5
    double precision:: fc_VS5
    !parameter for Eq. (10)
    double precision:: etak
    double precision:: etap
    double precision:: Y_VS5
    double precision:: G1_VS5
    double precision:: G2_VS5
    double precision:: eta_VS5
    !parameters for eq.(11)
    !a0_vs1 and a1_vs1 were determined by regression of the viscosity data for nonpolar fluids
    !a2_vs1 and a3 were determinde from the data of polar and associating fluid
    double precision, dimension(10):: Ai_VS5
    double precision, dimension(10):: a0_VS5
    double precision, dimension(10):: a1_VS5
    double precision, dimension(10):: a2_VS5
    double precision, dimension(10):: a3_VS5

    !VS9
    double precision:: etacorr_VS9, eta_1, tred, rhored
    double precision:: sum1_VS9, sum2_VS9, sum3_VS9, e_VS9
    double precision:: sigma_VS9, tau_VS9
    double precision:: sum3par1_VS9, sum3par2_VS9, ini_dens, sum_ini_dens
    double precision:: sum_comp_exp, sum_poly, sum_exp, sum_crit
    double precision:: term1, term2, term3, term4, term7
    integer:: nterm_normal

    !VS hydrogen peroxide
    double precision:: A_h2o2_liq, B_h2o2_liq, C_h2o2_liq, D_h2o2_liq, d_crit_h2o2
    double precision:: A_h2o2_vap, B_h2o2_vap, C_h2o2_vap

    !----------------------------------- VS0 ---------------------------------------------
    !parameters of hydrogen: C. D. Muzny, M. L. Huber, A. F. Kazakov, J. Chem. Eng. Data, 58:969, 2013
    double precision :: mwH2
    double precision :: epskH2
    double precision :: sigmaH2
    double precision :: rhosc, rhospec
    !-----------------------------------------------------------------------------------
    !parameters of helium
    double precision :: eta0dash, etaedash, eta0m
    double precision :: BHel, CHel, DHel, dspec

    !parameters of LJF

    !Leslie V. Woodcock, Equation of State for the Viscosity of Lennard-Jones Fluids

    double precision, dimension(9) :: bi_LJF !Table 2
    double precision, dimension(9) :: ti_LJF

    double precision, dimension(6) :: ai22_LJF !Table 1

    double precision :: B_etaLJF_star !eq. 18
    double precision :: C_AH !Ashurst-Hoover constant S.441

    double precision :: omega_22 !collision integrals eq. 13
    double precision :: diff_omega_22, diff_omega_23, omega_22_exp
    double precision :: omega_23
    double precision :: omega_24

    double precision :: f_11
    double precision :: f_12
    double precision :: f_22

    double precision :: f_eta

    double precision :: eta_ideal

    double precision :: sum_ai22, diff_sum !part of eq. 13
    double precision :: exp_sum_ai22


    VISDYN_CALC=0.d0
    mue=0.d0
    mue0=0.d0
    mue1=0.d0
    mue2=0.d0
    tdiml=0.d0
    rhodiml=0.d0
    Hi=0.d0
    sum1=0.d0
    sum_new=0.d0
    nue=0.d0
    nue0=0.d0
    nuer=0.d0
    n=0
    m=0
    i=0
    j=0
    k=0
    Y_mue2=0.d0
    etaexess=0.d0
    sstar=0.d0
    lnsstar=0.d0
    etal=0.d0
    etadg=0.d0
    collsum=0.d0
    ci=0.d0

    !Parameters VS1 model
    tstar=0.d0
    eta0_lc=0.d0
    betastar=0.d0
    beta_lc=0.d0
    del0=0.d0
    sp=0.d0
    etah=0.d0
    norp=0.d0
    dorp=0.d0
    nexpeta=0.d0
    dexpeta=0.d0

    If ((T >= gl%visco%upp_temp_eta(nrsubst)) .OR. (T <= gl%visco%low_temp_eta(nrsubst))) then
        ! error ILIMIT TODO
    endif


    !--------------------------------------------------------------------------------------------------------
    ! Calculating viscosity for [V]iscosity [S]pecification model [0]
    !--------------------------------------------------------------------------------------------------------
    If ((gl%visco%etamodel(nrsubst) == 'VS0') .and. (gl%visco%pointer_lambda(nrsubst) == 'H2O')) then       !begin loop etamodel water

        within_limits = .not. gl%hold_limits
        if (gl%hold_limits) then           ! Range of validity:
            p = P_CALC(gl,T,D,nrsubst)
            ! Pressure range for melting temperature iteration
            if ((p <=  gl%visco%pmaxeta(nrsubst)).and.(p > 0.d0)) then
                ! Melting temperature Iteration:
                Tmelt_eta = tmelt_eq(gl,p, nrsubst)

                !Melting temperature dependent range of validity :

                if  (((p <= 1000.d0).and.(p > 500.d0) .and. (Tmelt_eta <= T) .and. (T <= 373.15d0)) .or. &
                    & ((p <= 500.d0) .and.(p > 350.d0) .and. (Tmelt_eta <= T) .and. (T <= 433.15d0)) .or. &
                    & ((p <= 350.d0) .and.(p > 300.d0) .and. (Tmelt_eta <= T)  .and. (T <= 873.15d0)) .or. &
                    & ((p <= 300.d0) .and.(p >= gl%ptp(nrsubst)) .and. (Tmelt_eta <= T)  .and. (T <= 1173.15d0)) .or. &
                    & ((p < gl%ptp(nrsubst)) .and.(p > 0.d0).and. (273.16 <= T).and. (T <= 1173.15d0))) then
                    within_limits = .true.
                else
                    VISDYN_CALC=-9916.d0
                end if
            else
                VISDYN_CALC=-9916.d0
            end if
        endif
        if (within_limits) then
            tdiml = t / gl%visco%red_temp_eta(nrsubst) ! has to be Input temperature / Reference Temperature
            rhodiml = d / gl%visco%red_dens_eta(nrsubst) ! has to be Input density / Reference Density
            tau = 1.d0/tdiml-1.d0
            delta = rhodiml-1.d0

            ! mue0  < --Viscosity in the dilute-gas limit
            DO n = 1, 4
                sum1 = sum1 + gl%visco%coeff_hieta(nrsubst,n) / tdiml**(n-1)
            END DO

            mue0 = tdiml**0.5d0/sum1

            ! mue1  < --Contribution to viscosity due to finite density
            sum1_mue1=0.d0
            DO i=1, 6
                sum2_mue1=0.d0
                DO j=1, 7
                    sum2_mue1 = sum2_mue1 + gl%visco%hieta(nrsubst,i,j) * (delta)**(j-1)
                END DO
                sum1_mue1 = sum1_mue1 + ((tau)**(i-1)) * sum2_mue1
            END DO


            mue1 = dexp(rhodiml * sum1_mue1)

            !mue2  < --represents the critical enhancement of viscosity

            !IF ((t > lowbo_temp_mue2(nrsubst)) .AND. (t < uppbo_temp_mue2(nrsubst)) &      ! Case selection if close enough to critical Temperature and Density
            !& .AND. (lowbo_dens_mue2(nrsubst) < d) .AND.(d < uppbo_dens_mue2(nrsubst))) THEN

            !Berechnung von delchi_mue2
            dpdd_state=DPDD_CALC(gl,T,D, NRSUBST)
            dddp_state=1.0d0/dpdd_state/D
            tref_mue=gl%visco%T_mue2(nrsubst)*gl%Tc(nrsubst)
            dpdd_ref=DPDD_CALC(gl,tref_mue,D, nrsubst)
            dddp_ref=1.d0/dpdd_ref/D
            delchi_mue2 = rhodiml**2*gl%pc(nrsubst) * (dddp_state - dddp_ref * gl%visco%T_mue2(nrsubst)/tdiml)                         !Eq. (21)

            !Berechnung von Xi_mue2
            Xi_mue2 = gl%visco%Xi0_mue2(nrsubst) * (delchi_mue2/gl%visco%gamma0_mue2(nrsubst))**(gl%visco%nue_mue2(nrsubst)/gl%visco%gamma_mue2(nrsubst))      !Eq. (20)

            if (Xi_mue2 < 0.d0) Xi_mue2=0.d0

            if(gl%visco%pointercrit_eta(nrsubst) == 'I08') then

                IF ((Xi_mue2 >= 0) .AND. (Xi_mue2 <= gl%visco%corr_ref_mue2(nrsubst))) THEN             !Case selection of Y_mue2: Eq. (15)

                    Y_mue2 = 0.2d0 * gl%visco%qc_mue2(nrsubst) * Xi_mue2 * ((gl%visco%qd_mue2(nrsubst) * Xi_mue2)**5) * &
                        & (1.d0 - gl%visco%qc_mue2(nrsubst) * Xi_mue2 + (gl%visco%qc_mue2(nrsubst) * Xi_mue2)**2 &
                        & - 765.d0/504.d0 * (gl%visco%qd_mue2(nrsubst) * Xi_mue2)**2)

                ELSEIF (Xi_mue2 > gl%visco%corr_ref_mue2(nrsubst)) THEN                                 !Case selection of Y_mue2: Eq. (16)
                    !Berechnung von psid_mue2
                    psid_mue2 = acos((1.d0 + (gl%visco%qd_mue2(nrsubst)**2) * (Xi_mue2)**2)**(-1.d0/2.d0))
                    !Berechnung von w
                    w_mue2 = (abs((gl%visco%qc_mue2(nrsubst)*Xi_mue2 - 1.d0) / (gl%visco%qc_mue2(nrsubst)*Xi_mue2 + 1.d0)))**(1.d0/2.d0) &
                        &  * tan(psid_mue2 / 2.d0)
                    !Berechnung von L(w)-Funktion
                    IF ((gl%visco%qc_mue2(nrsubst) * Xi_mue2) > 1) THEN                                 !case selection of lw-function
                        lw_mue2 = log((1.d0 + w_mue2) / (1.d0 - w_mue2))
                    ELSEIF ((gl%visco%qc_mue2(nrsubst) * Xi_mue2) <= 1) THEN                            !case selection of lw-function
                        lw_mue2 = 2.d0 * atan(abs(w_mue2))
                    ENDIF                                                                      !case selection of lw-function


                    Y_mue2 = (1.d0 / 12.d0) * sin(3.d0 * psid_mue2) &
                        & - (1.d0 / (4.d0 * gl%visco%qc_mue2(nrsubst) * Xi_mue2)) * sin(2.d0 * psid_mue2) &
                        & + (1.d0 / (gl%visco%qc_mue2(nrsubst) * Xi_mue2)**2) * (1.d0 - 1.25d0 * (gl%visco%qc_mue2(nrsubst) * Xi_mue2)**2) &
                        & * sin(psid_mue2) - (1.d0 / (gl%visco%qc_mue2(nrsubst) * Xi_mue2)**3) * &
                        & ((1.d0 - 1.5d0 * (gl%visco%qc_mue2(nrsubst) * Xi_mue2)**2) * psid_mue2 - &
                        & ((gl%visco%qc_mue2(nrsubst) * psid_mue2)**2 - 1.d0)**(3.d0/2.d0) * lw_mue2)
                ENDIF                                                                          !end case selection of Y_mue2
                mue2 = dexp(gl%visco%X_mue2(nrsubst)*Y_mue2)
            ELSE
                mue2 = 1.d0
            end if
            VISDYN_CALC = gl%visco%red_vis_eta(nrsubst) * mue0 * mue1 * mue2
        end if
        !------------------------------------------------------------------------------------------------------------------------
        !------------------------------------------------------------------------------------------------------------------------
    ELSEIf ((gl%visco%etamodel(nrsubst) == 'VS0') .and. (gl%visco%pointer_lambda(nrsubst) == 'MEO')) then
        !M. Thol, November 2013
        !viscosity model for methanol of H. W. Xiang, A. Laesecke, M. L. Huber, J. Phys. Chem. Ref. Data, 35(4): 1597-1, 2006
        !------------------------------------------------------------------------------------------------------------------------
        fcross=0.d0
        etagas=0.d0
        etaliq=0.d0
        eta0_lc=0.d0
        B_eta=0.d0
        B_etar=0.d0
        C_etar=0.d0
        C_eta=0.d0
        etagasred=0.d0
        sigceta=0.d0

        !Parameters for the collision integrals
        !------------------------
        !Lennard_Jones: eq. (4)
        gl%visco%a_OH1(1)=1.16145d0
        gl%visco%a_OH1(2)=-0.14874d0
        gl%visco%a_OH1(3)=0.52487d0
        gl%visco%a_OH1(4)=-0.7732d0
        gl%visco%a_OH1(5)=2.16178d0
        gl%visco%a_OH1(6)=-2.43787d0
        !-------------------------
        !Stockmayer: eq. (8)
        gl%visco%a_OH1(7)=0.95976d-3
        !-------------------------
        !residual reduced collision integral: eq. (9)
        gl%visco%a_OH1(8)= 0.10225d0
        gl%visco%a_OH1(9)=-0.97346d0
        gl%visco%a_OH1(10)= 0.10657d0
        gl%visco%a_OH1(11)=-0.34528d0
        gl%visco%a_OH1(12)=-0.44557d0
        gl%visco%a_OH1(13)=-2.58055d0
        !-------------------------

        !Parameters for the second viscosity virial coefficient: eq. (13)
        gl%visco%b_OH1(1)=-0.19572881d2
        gl%visco%b_OH1(2)= 0.21973999d3
        gl%visco%b_OH1(3)=-0.10153226d4
        !b_OH1(4)= 0.24710125d4  !Paper
        gl%visco%b_OH1(4)= 0.247101251d4  !REFPROP
        gl%visco%b_OH1(5)=-0.33751717d4
        gl%visco%b_OH1(6)= 0.24916597d4
        gl%visco%b_OH1(7)=-0.78726086d3
        gl%visco%b_OH1(8)= 0.14085455d2
        gl%visco%b_OH1(9)=-0.34664158d0

        !Parameters for the third viscosity virial coefficient: eq. (14)
        gl%visco%c_OH1(1)=0.186222085d-2
        gl%visco%c_OH1(2)=0.9990338d1

        !Parameters for the temperature and density dependency of sigHS: eq. (17)
        gl%visco%d_OH1(1)=-0.1181909d1
        gl%visco%d_OH1(2)= 0.5031030d0
        gl%visco%d_OH1(3)=-0.6268461d0
        gl%visco%d_OH1(4)= 0.5169312d0
        gl%visco%d_OH1(5)=-0.2351349d0
        gl%visco%d_OH1(6)= 0.53980235d-1
        gl%visco%d_OH1(7)=-0.49069617d-2
        !-------------------------
        gl%visco%e_OH1(1)= 0.4018368d1
        gl%visco%e_OH1(2)=-0.4239180d1
        gl%visco%e_OH1(3)= 0.2245110d1
        gl%visco%e_OH1(4)=-0.5750698d0
        gl%visco%e_OH1(5)= 0.23021026d-1
        gl%visco%e_OH1(6)= 0.25696775d-1
        gl%visco%e_OH1(7)=-0.68372749d-2
        gl%visco%e_OH1(8)= 0.72707189d-3
        gl%visco%e_OH1(9)=-0.29255711d-4

        Tstar=T/epsk
        Tr=T/Tcrit
        DSI=D*mw_oh1/gl%factor
        Dr=DSI/Dcrit

        !collision integral, Lennard-Jones potential: eq. (4)
        ciLJ=gl%visco%a_OH1(1)*Tstar**gl%visco%a_OH1(2)+gl%visco%a_OH1(3)*dexp(gl%visco%a_OH1(4)*Tstar)+gl%visco%a_OH1(5)*dexp(gl%visco%a_OH1(6)*Tstar)

        !residual reduced collision integral: eq. (9)
        ciDEL=gl%visco%a_OH1(8)*Tstar**gl%visco%a_OH1(9)+gl%visco%a_OH1(10)*dexp(gl%visco%a_OH1(11)*Tstar)+gl%visco%a_OH1(12)*dexp(gl%visco%a_OH1(13)*Tstar)

        !collision integral, Stockmayer potential: eq. (8)
        ciSM=ciLJ*(1.d0+delst**2/(1.d0+gl%visco%a_OH1(7)*delst**6)*ciDEL)

        !viscosity of diluted gas: eq. (3)
        eta0_lc=(5.d0/16.d0)*(mw_oh1/avo*1.d-3*kbol*T/PI)**0.5d0/(sig0**2*ciSM)

        !Second viscosity virial coefficient: eq. (13)
        do i=1,7
            B_etar=B_etar+gl%visco%b_OH1(i)/Tstar**(0.25d0*(i-1.d0))
        end do
        B_etar=B_etar+gl%visco%b_OH1(8)/Tstar**2.5d0+gl%visco%b_OH1(9)/Tstar**5.5d0
        B_eta=B_etar*avo*sig0**3

        !Third viscosity virial coefficient: eq. (14)
        C_etar=gl%visco%c_OH1(1)*Tstar**3*dexp(gl%visco%c_OH1(2)/Tstar**0.5d0)
        C_eta=C_etar*(avo*sig0**3)**2

        !viscosity at low densities: eq. (2)
        etagasred=1.d0+B_eta*D+C_eta*D**2
        etagas=etagasred*eta0_lc

        !radial distribution function: eq. (16)
        sigceta=(6.d0*gl%wm(nrsubst)/PI/dcrit/avo*1.d-3)**(1.d0/3.d0)    !Pay attention: Equation is wrong in the paper

        sigHSr=0.d0
        do i=1,7
            sigHSr=sigHSr+gl%visco%d_OH1(i)/Tr**(i-1)
        end do
        do i=1,9
            sigHSr=sigHSr+gl%visco%e_OH1(i)*Dr**i
        end do
        sigHS=sigHSr*sigceta

        bvol=2.d0*PI*avo*sigHS**3/3.d0
        ksi=bvol*D/4.d0
        gHS=(1.d0-0.5d0*ksi)/(1.d0-ksi)**3

        !high-density region: eq. (15)
        etaliq=1.d0/gHS+0.8d0*bvol*D+0.761d0*gHS*(bvol*D)**2

        !transition function from low density to high density
        fcross=1.d0/(1.d0+dexp(5.d0*(Dr-1.d0)))

        !summation of the different contributions: eq. (1)
        VISDYN_CALC=eta0_lc*(fcross*etagasred+(1.d0-fcross)*etaliq)*1.d6

        !--------------------------------------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------------------------
    ELSEIf ((gl%visco%etamodel(nrsubst) == 'VS0') .and. (gl%visco%pointer_lambda(nrsubst) == 'ETY')) then
        !M. Thol, November 2013
        !viscosity model for ethylene of P. M. Holland, B. E. Eaton, H. J. M. Hanley, J. Phys. Chem. Ref. Data, 12(4):917, 1983
        !--------------------------------------------------------------------------------------------------------------------------
        eta0_lc=0.d0
        eta1=0.d0
        etadash=0.d0
        deleta=0.d0
        deletacrit=0.d0

        DSI=D*gl%wm(nrsubst)*1.d-3    !g/cm³
        THETA=(DSI-rhoc_ETHY)/rhoc_ETHY

        !Parameters for eq. (7)
        gl%visco%GV(1)=-3.5098225018D6
        gl%visco%GV(2)= 2.5008406184D6
        gl%visco%GV(3)=-5.8365540744D5
        gl%visco%GV(4)= 4.5549146583D3
        gl%visco%GV(5)= 2.2881683403D4
        gl%visco%GV(6)=-4.7318682077D3
        gl%visco%GV(7)= 4.5022249258D2
        gl%visco%GV(8)=-2.1490688088D1
        gl%visco%GV(9)= 4.1649263233D-1

        !Parameters for eq. (3)
        gl%visco%A_ETHY=0.d0
        gl%visco%B_ETHY=0.d0
        gl%visco%C_ETHY=0.d0
        gl%visco%F_ETHY=358.9d0

        !Parameter for eq. (4)
        gl%visco%E_ETHY=1.d0

        gl%visco%J_ETHY(1)=-4.8544486732d0
        gl%visco%J_ETHY(2)= 1.3033585236d1
        gl%visco%J_ETHY(3)= 2.7808928908d4
        gl%visco%J_ETHY(4)=-1.8241971308d3
        gl%visco%J_ETHY(5)= 1.5913024509d0
        gl%visco%J_ETHY(6)=-2.0513573927d2
        gl%visco%J_ETHY(7)=-3.9478454708d4

        !diluted gas: eq. (7)
        TT = T**(1.d0/3.d0)
        eta0_lc= gl%visco%gv(1)/t+gl%visco%gv(2)/tt**2+gl%visco%gv(3)/tt+gl%visco%gv(4)+gl%visco%gv(5)*tt+gl%visco%gv(6)*tt**2+gl%visco%gv(7)*t+gl%visco%gv(8)*tt**4+gl%visco%gv(9)*tt**5

        !1st density correction for the moderately dense gas: eq. (3)
        !actually that is: eta1=A_ETHY+B_ETHY*(C_ETHY-DLOG(T/F_ETHY))**2
        !as the coefficients are 0 = >  eta1=0
        eta1=0.d0

        !remainder: eq. (4)
        etadash=dexp(gl%visco%J_ETHY(1)+gl%visco%J_ETHY(4)/t)*(dexp(DSI**0.1d0*(gl%visco%J_ETHY(2)+gl%visco%J_ETHY(3)/t**1.5d0) &
            &     +theta*DSI**0.5d0*(gl%visco%J_ETHY(5)+gl%visco%J_ETHY(6)/t+gl%visco%J_ETHY(7)/t**2))-1.0d0)

        !critical enhancement
        etacrit=0.d0

        !summation over the 4 different terms: eq. (1)
        VISDYN_CALC=(eta0_lc + eta1 + etadash + etacrit)/10.d0

        !-------------------------------------------------------------------------------------------------------------
        !-------------------------------------------------------------------------------------------------------------
    ELSEIf ((gl%visco%etamodel(nrsubst) == 'VS0') .and. (gl%visco%pointer_lambda(nrsubst) == 'H2A')) then
        !M. Thol, July 2014
        !viscosity model for hydrogen of C. D. Muzny, M. L. Huber, A. F. Kazakov, J. Chem. Eng. Data, 58:969, 2013
        !CAUTION: The source code does not fit to the publication. Modifications are taken from REFPROP,
        !         errors in the paper are suspected. M. Thol, Feb 2015
        !-------------------------------------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------------------------------
        !Parameters for the hydrogen model
        mwH2=2.01588d0
        epskH2=30.41d0      !K
        sigmaH2=0.297d0     !nm
        rhosc=90.5d0        !kg/m^3

        rhoscmol=rhosc/gl%wm(nrsubst)
        rhospec=D*gl%wm(nrsubst)

        aiH2(1)= 2.09630d-1
        aiH2(2)=-4.55274d-1
        aiH2(3)= 1.43602d-1
        aiH2(4)=-3.35325d-2
        aiH2(5)= 2.76981d-3

        biH2(1)=-0.1870d0
        biH2(2)= 2.4871d0
        biH2(3)= 3.7151d0
        biH2(4)=-11.0972d0
        biH2(5)= 9.0965d0
        biH2(6)=-3.8292d0
        biH2(7)= 0.5166d0

        ciH2(1)= 6.43449673d-3
        ciH2(2)= 4.56334068d-2
        ciH2(3)= 2.32797868d-1
        ciH2(4)= 9.58326120d-1
        ciH2(5)= 1.27941189d-1
        ciH2(6)= 3.63576595d-1

        Tstar=T/epskH2
        rhoredH2=d/rhoscmol
        rhoredH2=rhospec*0.011d0
        TredH2=T/gl%Tc(nrsubst)

        !eq. (4)
        do i=1,5
            lnsstar = lnsstar + aiH2(i)*(dlog(tstar))**(i-1)
        end do
        sstar = dexp(lnsstar)

        !eq. (3)
        eta0_lc = 0.021357d0*(mwH2*T)**0.5d0/(sigmaH2**2*sstar)

        !eq. (7)
        do i=1,7
            Betastar = betastar + biH2(i)/Tstar**(i-1)
        end do

        Betastar=Betastar*0.6022137d0

        !eq. (6)
        beta_lc = sigmaH2**3*Betastar

        !eq. (5)
        eta1 = beta_lc*eta0_lc*D*1.d-3

        !eq. (9) - first part
        etacrit = ciH2(1)*rhoredH2**2*dexp(ciH2(2)*TredH2+ciH2(3)/TredH2+ciH2(4)*rhoredH2**2/(ciH2(5)+TredH2)+ciH2(6)*rhoredH2**6)

        !eq. (9) - second part
        VISDYN_CALC = (eta0_lc + eta1 + etacrit*1.d3)!/10.d0

        !--------------------------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------------
        !M. Thol, Feb 2015
        !viscosity model for helium of V.D. Arp, R.D. McCarty, and D.G. Friend, NIST Technical Note 1334 (revised), 1998
        !--------------------------------------------------------------------------------------------------------------
    ELSEIf ((gl%visco%etamodel(nrsubst) == 'VS0') .and. (gl%visco%pointer_lambda(nrsubst) == 'HE')) then

        !convert into specific unit
        dspec = d * gl%wm(nrsubst) * 1.d-3

        if (T .le. 300.d0) then
            tau = dlog(T)
        else
            tau = dlog(300.d0)
            !Citation from paper:
            !Since Steward's analysis did not include any dense-gas data for temperatures above 300 K, either calculated or
            !experimental, the temperature-dependent excess function given by eq (9) was evaluated at 300 K.
        end if

        !eq. (2)
        eta0dash = -0.135311743d0 / tau + 1.00347841d0 + 1.20654649d0 * tau - 0.149564551d0 * tau**2 + 0.0125208416d0 * tau**3

        !eq. (4)
        BHel = -47.5295259d0 / tau + 87.6799309d0 - 42.0741589d0 * tau + 8.33128289d0 * tau**2 - 0.589252385d0 * tau**3

        !eq. (5)
        CHel = 547.309267d0 / tau - 904.870586d0 + 431.404928d0 * tau - 81.4504854d0 * tau**2 + 5.37008433d0 * tau**3

        !eq. (6)
        DHel = -1684.39324d0 / tau + 3331.0863d0 - 1632.19172d0 * tau + 308.804413d0 * tau**2 - 20.2936367d0 * tau**3

        !eq. (3)
        etaedash = dspec * BHel + dspec**2 * CHel + dspec**3 * DHel

        if (T <= 100.d0) then
            !eq. (1)
            VISDYN_CALC = dexp(eta0dash + etaedash) * 1.d-1
            return
        elseif (T > 100.d0) then
            eta0m = 196.d0 * T ** 0.71938d0 * dexp(12.451d0/T - 295.67/T**2 - 4.1249d0)
            eta0dash = dexp(eta0dash)
            etaedash = dexp(etaedash)
            if (T < 110.d0) then
                eta0_lc = eta0dash + (eta0m - eta0dash) * (T - 100.d0)*1.d-1
            else
                eta0_lc = eta0m
            end if
            VISDYN_CALC = (eta0dash*etaedash+eta0_lc-eta0dash) * 1.d-1
            return
        end if

        !--------------------------------------------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------------------------------
        !M. Thol, Feb 2015
        !viscosity model for neon of  V.A. Rabinovich, A.A. Vasserman, V.I. Nedostup, and L.S. Veksler, Hemisphere Publishing Corp., 1988
        !--------------------------------------------------------------------------------------------------------------------------------
    ELSEIf ((gl%visco%etamodel(nrsubst) == 'VS0') .and. (gl%visco%pointer_lambda(nrsubst) == 'NEO')) then

        VISDYN_CALC = -5243.D0
        return

        !--------------------------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------------
    ELSEIf ((gl%visco%etamodel(nrsubst) == 'VS0') .and. (gl%visco%pointer_lambda(nrsubst) == 'LJF')) then

        bi_LJF(1)=-18.710845605d0
        bi_LJF(2)=208.37732997d0
        bi_LJF(3)=-949.49804258d0
        bi_LJF(4)=2247.6240206d0
        bi_LJF(5)=-2878.5655072d0
        bi_LJF(6)=1735.2192655d0
        bi_LJF(7)=-523.67304508d0
        bi_LJF(8)=179.42313652d0
        bi_LJF(9)=-0.31888014252d0

        ti_LJF(1)=0.d0
        ti_LJF(2)=0.25d0
        ti_LJF(3)=0.5d0
        ti_LJF(4)=0.75d0
        ti_LJF(5)=1.d0
        ti_LJF(6)=1.25d0
        ti_LJF(7)=1.75d0
        ti_LJF(8)=2.d0
        ti_LJF(9)=5.5d0

        ai22_LJF(1)=3.108096d-1
        ai22_LJF(2)=-1.712113d-1
        ai22_LJF(3)=-7.158049d-1
        ai22_LJF(4)=2.486777d0
        ai22_LJF(5)=-1.783167d0
        ai22_LJF(6)=3.944048d-1

        B_etaLJF_star=0.d0
        C_AH=3.025d0 !Ashurst and Hoover

        do i=1,9
            B_etaLJF_star=B_etaLJF_star+bi_LJF(i)*(1/T)**(ti_LJF(i)) !Rainwater and Friend
        end do

        sum_ai22=0.d0

        do i=1,6
            exp_sum_ai22=(i-1)*0.5d0
            sum_ai22=sum_ai22+ai22_LJF(i)*(1/T)**exp_sum_ai22
        end do

        omega_22=exp((1/6.0d0)*dlog(1.d0/T)+dlog(17.d0/18.d0)+sum_ai22)

        omega_22_exp=exp(T**(-2.5d0)*ai22_LJF(6) + T**(-1.5d0)*ai22_LJF(4) + T**(-0.5d0)*ai22_LJF(2) + ai22_LJF(1) + ai22_LJF(3)/T + ai22_LJF(5)/(T**2.0d0))
        !------------------------------------------------------------------------------------------------------------------!richtig bis hier
        diff_omega_22=-(17.d0/108.d0)*T**(-7.d0/6.d0)*omega_22_exp + (17.d0/18.d0)*T**(-1.d0/6.d0)*(-2.5d0*T**(-3.5d0)*ai22_LJF(6) - 1.5d0*T**(-2.5)*ai22_LJF(4) - 0.5d0*T**(-1.5d0)*ai22_LJF(2) - ai22_LJF(3)/T**2.d0 - 2.d0*ai22_LJF(5)/T**3.d0)*omega_22_exp


        omega_23=omega_22+(1/4.d0)*T*diff_omega_22



        !omega_23 mit (1/5)
        !diff_omega_23=-(17.d0/90.d0)*T**(-7.d0/6.d0)*omega_22_exp + (17.d0/15.d0)*T**(-1.d0/6.d0)*(-2.5d0*T**(-3.5d0)*ai22_LJF(6) - 1.5d0*T**(-2.5d0)*ai22_LJF(4) - 0.5d0*T**(-1.5d0)*ai22_LJF(2) - ai22_LJF(3)/T**2.d0 - 2.d0*ai22_LJF(5)/T**3.d0)*omega_22_exp + 0.2d0*T*(0.183641975308642*T**(-13.d0/6.d0)*omega_22_exp - (17.d0/54.d0)*T**(-7.d0/6.d0)*(-2.5d0*T**(-3.5d0)*ai22_LJF(6) - 1.5d0*T**(-2.5d0)*ai22_LJF(4) - 0.5d0*T**(-1.5d0)*ai22_LJF(2) - ai22_LJF(3)/T**2.d0 - 2.d0*ai22_LJF(5)/T**3.d0)*omega_22_exp + (17.d0/18.d0)*T**(-1.d0/6.d0)*(8.75d0*T**(-4.5d0)*ai22_LJF(6) + 3.75d0*T**(-3.5d0)*ai22_LJF(4) + 0.75d0*T**(-2.5d0)*ai22_LJF(2) + 2.d0*ai22_LJF(3)/T**3.d0 + 6.d0*ai22_LJF(5)/T**4.d0)*omega_22_exp + (17.d0/18.d0)*T**(-1.d0/6.d0)*(-2.5d0*T**(-3.5d0)*ai22_LJF(6) - 1.5d0*T**(-2.5d0)*ai22_LJF(4) - 0.5d0*T**(-1.5d0)*ai22_LJF(2) - ai22_LJF(3)/T**2.d0 - 2.d0*ai22_LJF(5)/T**3.d0)**2.d0*omega_22_exp)

        !omega_23 mit (1/4)
        diff_omega_23=-(17.d0/108.d0)*T**(-7.d0/6.d0)*omega_22_exp + (85.d0/72.d0)*T**(-1.d0/6.d0)*(-2.5d0*T**(-3.5d0)*ai22_LJF(6) - 1.5d0*T**(-2.5d0)*ai22_LJF(4) - 0.5d0*T**(-1.5d0)*ai22_LJF(2) - ai22_LJF(3)/T**2.d0 - 2.d0*ai22_LJF(5)/T**3.d0)*omega_22_exp + 0.2d0*T*(0.183641975308642*T**(-13.d0/6.d0)*omega_22_exp - (17.d0/54.d0)*T**(-7.d0/6.d0)*(-2.5d0*T**(-3.5d0)*ai22_LJF(6) - 1.5d0*T**(-2.5d0)*ai22_LJF(4) - 0.5d0*T**(-1.5d0)*ai22_LJF(2) - ai22_LJF(3)/T**2.d0 - 2.d0*ai22_LJF(5)/T**3.d0)*omega_22_exp + (17.d0/18.d0)*T**(-1.d0/6.d0)*(8.75d0*T**(-4.5d0)*ai22_LJF(6) + 3.75d0*T**(-3.5d0)*ai22_LJF(4) + 0.75d0*T**(-2.5d0)*ai22_LJF(2) + 2.d0*ai22_LJF(3)/T**3.d0 + 6.d0*ai22_LJF(5)/T**4.d0)*omega_22_exp + (17.d0/18.d0)*T**(-1.d0/6.d0)*(-2.5d0*T**(-3.5d0)*ai22_LJF(6) - 1.5d0*T**(-2.5d0)*ai22_LJF(4) - 0.5d0*T**(-1.5d0)*ai22_LJF(2) - ai22_LJF(3)/T**2.d0 - 2.d0*ai22_LJF(5)/T**3.d0)**2.d0*omega_22_exp)

        omega_24=omega_23+(1/5.0d0)*T*diff_omega_23

        f_11=4.d0*omega_22
        f_12=7.d0*omega_22-8.d0*omega_23
        f_22=(301.d0/12.d0)*omega_22-28.d0*omega_23+20.d0*omega_24

        f_eta=1.d0+((f_12**2.d0)/(f_11*f_22-(f_12**2.d0)))

        eta_ideal=(5.d0/(16.d0*pi**0.5d0))*(T**0.5d0)*(f_eta/omega_22)

        VISDYN_CALC = eta_ideal*(1.d0+B_etaLJF_star*d+C_AH*(1.d0/T)**(1.d0/3.d0)*d**4.d0) !C_AH falsch ?
        return



    elseif((gl%visco%etamodel(nrsubst) == 'VS0') .and. (gl%visco%pointer_lambda(nrsubst) == 'LJ1')) then

        VISDYN_CALC = Ruckenstein_IECR_1997(T,D)

        !elseif((gl%visco%etamodel(nrsubst) == 'VS4') .and. (gl%visco%pointer_lambda(nrsubst) == 'LJF')) then
        !
        !    VISDYN_CALC = Rowley_IJT_1997(T,D)

        !--------------------------------------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------------------------------------
    ELSEIf ((gl%visco%etamodel(nrsubst) == 'VS0') .and. (gl%visco%pointer_lambda(nrsubst) == 'R23')) then
        !M. Thol, July 2014
        !viscosity model for R23 of Z. Shan, S. G. Penoncello, R. T. Jacobsen, ASHRAE Transactions, 106:1-11 (2000)
        !--------------------------------------------------------------------------------------------------------------
        Tstar = T/gl%visco%lejo_epka(nrsubst)

        !eq. (5)
        DO i = 1,gl%visco%colltermnr(nrsubst)
            collsum=collsum+gl%visco%coeffcoll(nrsubst,i)*dlog(tstar)**gl%visco%powercoll(nrsubst,i)
        END DO
        ci=dexp(collsum)

        !eq. (4)
        etadg=gl%visco%Chapman_Enskog1(nrsubst)*T**0.5d0/(gl%visco%lejo_sigma(nrsubst)**2*ci)

        call R_mix_calc(gl,Rmix)
        !eq. (7)
        etal=gl%visco%c2_R23*(gl%visco%rhol_R23*1.d-3)**2/((gl%visco%rhol_R23-D)*1.d-3)*T**0.5d0*dexp(D/(gl%visco%rhol_R23-D)*gl%visco%delg_R23/Rmix/T)

        ksi=(d-gl%rhoc(nrsubst))*1.d-3
        tau=t-gl%tc(nrsubst)
        !eq.(8a)
        etacrit=4.d0*gl%visco%etamax_R23/((dexp(ksi)+dexp(-ksi))*(dexp(tau)+dexp(-tau)))

        !eq. (3)
        !test1=((rhol_R23-D)/rhol_R23)**c1_R23*etadg
        !test2=(D/rhol_R23)**c1_R23*etal
        VISDYN_CALC = ((gl%visco%rhol_R23-D)/gl%visco%rhol_R23)**gl%visco%c1_R23*etadg +(D/gl%visco%rhol_R23)**gl%visco%c1_R23*etal+etacrit

        return


        !--------------------------------------------------------------------------------------------------------
        ! Calculating viscosity for [V]iscosity [S]pecification model [1]
        !--------------------------------------------------------------------------------------------------------

    ELSE IF (gl%visco%etamodel(nrsubst) == 'VS1') THEN
        !tstar=t/eklj(nrsubst)   Monika: changed April 2017
        tstar=t/gl%visco%treddens(nrsubst)
        !tau=tredeta(nrsubst)*t**(-1.d0)
        tau=gl%visco%tredeta(nrsubst)/t
        del=d/gl%visco%rhoredeta(nrsubst)

        !eta0  < --dilute gas viscosity
        eta0_lc = eta0 (gl,T,nrsubst)

        !beta  < -- term for initial density dependence
        IF ((gl%visco%denstermnr(nrsubst) > 0) ) THEN  !.AND. (0.5 < tstar < 100)
            betastar=0.d0
            DO j=1,gl%visco%denstermnr(nrsubst)
                betastar=betastar+gl%visco%coeffbetastar(nrsubst,j)*tstar**gl%visco%powerbetastar(nrsubst,j)        !for R1234yf and R1234zeE coeffbetastar(1,4) is different in REFPROP than in the Paper
                !write(*,*) 'coeff (',j,') =',coeffbetastar(nrsubst,j)
                !write(*,*) 'expo (',j,') =',powerbetastar(nrsubst,j)
                !write(*,*) 'term (',j,') =', coeffbetastar(nrsubst,j)*tstar**powerbetastar(nrsubst,j)
                !write(*,*) 'eta1B2 =',betastar
                !pause
            END DO
            beta_lc=gl%visco%etaB2(nrsubst)*betastar
        ELSE
            beta_lc=0.d0
        ENDIF

        !etah  < -- residual viscosity

        !del0  < -- close packed density
        IF (gl%visco%term_num1_eta(nrsubst) > 0) THEN
            del0=0.d0
            DO k=1,gl%visco%term_num1_eta(nrsubst),1
                del0=del0+gl%visco%a_vs1(nrsubst,k)/tau**gl%visco%powertau(nrsubst,k)
            END DO

        ELSE IF (gl%visco%term_num1_eta(nrsubst) < 0 ) THEN
            del0=1.d0
            DO k=2,abs(gl%visco%term_num1_eta(nrsubst)),1
                del0=del0+gl%visco%a_vs1(nrsubst,k)/tau**gl%visco%powertau(nrsubst,k)
            END DO
            del0=gl%visco%a_vs1(nrsubst,1)/del0
        ELSE
            del0=1.d0
        END IF

        !sp  < -- simple polynomial terms
        IF (gl%visco%term_num2_eta(nrsubst) /= 0 ) THEN
            DO m=1,gl%visco%term_num2_eta(nrsubst)
                sp=gl%visco%b_vs1(nrsubst,m)/tau**gl%visco%ptausp(nrsubst,m)*del**gl%visco%pdelsp(nrsubst,m) &
                    & *del0**gl%visco%pdel0sp(nrsubst,m)
                IF (gl%visco%expdel(nrsubst,m) >= 1) THEN
                    sp=sp*dexp(-del**gl%visco%expdel(nrsubst,m))
                END IF
                etah=etah+sp
            END DO
        END IF

        !norp  < --numerator of rational polynomial term
        IF (gl%visco%term_num3_eta(nrsubst) /= 0) THEN
            norp=0.d0
            DO n=gl%visco%term_num2_eta(nrsubst)+1,gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst),1
                !original: norp=norp+b(nrsubst,n)/tau**ptausp(nrsubst,n)*del**pdelsp(nrsubst,n) &
                !Monika November 2013: in the fluid file the exponent ptausp is given with a negative sign.
                !For ethane, the sum is divided by "tau**ptausp(nrsubst,n)" instead of multiplied = >  could be a problem for other fluids
                norp=norp+gl%visco%b_vs1(nrsubst,n)/tau**gl%visco%ptausp(nrsubst,n)*del**gl%visco%pdelsp(nrsubst,n) &
                    & *del0**gl%visco%pdel0sp(nrsubst,n)
                IF (gl%visco%expdel(nrsubst,n) >= 1) THEN
                    norp=norp*dexp(-del**gl%visco%expdel(nrsubst,n))
                END IF
            END DO
        ELSE
            norp=1.d0
        END IF


        !dorp  < --denominator of rational polynomial term
        IF (gl%visco%term_num4_eta(nrsubst) /=  0) THEN
            dorp=0.d0
            DO l=gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+1, &
                & gl%visco%term_num2_eta(nrsubst)+gl%visco%term_num3_eta(nrsubst)+gl%visco%term_num4_eta(nrsubst),1
                !original: dorp=dorp+b(nrsubst,l)*tau**ptausp(nrsubst,l)*del**pdelsp(nrsubst,l) &
                !Monika November 2013: in the fluid file the exponent ptausp is given with a negative sign.
                !For ethane, the sum is divided by "tau**ptausp(nrsubst,n)" instead of multiplied = >  could be a problem for other fluids
                dorp=dorp+gl%visco%b_vs1(nrsubst,l)/tau**gl%visco%ptausp(nrsubst,l)*del**gl%visco%pdelsp(nrsubst,l) &
                    & *del0**gl%visco%pdel0sp(nrsubst,l)
                IF (gl%visco%expdel(nrsubst,l) >= 1) THEN
                    dorp=dorp*dexp(-del**gl%visco%expdel(nrsubst,l))
                END IF
            END DO
        ELSE
            dorp=1.d0
        END IF

        !combination of rational polynomial terms
        IF ((gl%visco%term_num3_eta(nrsubst) /= 0) .OR. (gl%visco%term_num4_eta(nrsubst) /= 0)) THEN
            etah=etah+norp/dorp
        END IF

        !nexpeta  < -- numerator of exponential term
        IF (gl%visco%term_num5_eta(nrsubst) /= 0) THEN
            nexpeta=0.d0
            DO o=l,l-1+gl%visco%term_num5_eta(nrsubst),1
                nexpeta=nexpeta+gl%visco%b_vs1(nrsubst,o)*tau**gl%visco%ptausp(nrsubst,o)*del**gl%visco%pdelsp(nrsubst,o) &
                    & *del0**gl%visco%pdel0sp(nrsubst,o)
            END DO
        ELSE
            nexpeta=1.d0
        END IF

        !dexpeta  < -- denominator of exponential term
        IF (gl%visco%term_num6_eta(nrsubst) /= 0) THEN
            dexpeta=0.d0
            DO pcount=o,o-1+gl%visco%term_num6_eta(nrsubst),1
                dexpeta=dexpeta+gl%visco%b_vs1(nrsubst,pcount)*tau**gl%visco%ptausp(nrsubst,pcount)*del**gl%visco%pdelsp(nrsubst,pcount) &
                    & *del0**gl%visco%pdel0sp(nrsubst,pcount)
            END DO
        ELSE
            dexpeta=1.d0
        END IF

        !combination of exponential terms
        IF (gl%visco%term_num5_eta(nrsubst) /= 0 .OR. gl%visco%term_num6_eta(nrsubst) /= 0) THEN
            etah=etah+dexp(nexpeta/dexpeta)
        END IF

        !Unit conversion
        etah=etah*gl%visco%visredeta(nrsubst)

        rho=D/gl%factor
        VISDYN_CALC=eta0_lc*(1.d0+beta_lc*rho)+etah

        !Lars / Andreas March 2014
        !Implementation of ECS models

        !###############################
        gl%residual_ecs = eta0_lc* beta_lc *rho + etah
        !###############################

        !Andreas, October 2013
        !If the viscosity is 0, give back an error
        if (VISDYN_CALC == 0.D0) then
            VISDYN_CALC = -5242.D0
        end if
        RETURN


        !--------------------------------------------------------------------------------------------------------
        ! Calculating viscosity for [V]iscosity [S]pecification model [2]
        ! M. Thol and M. Schaden, July 2015
        !--------------------------------------------------------------------------------------------------------

    ELSE IF (gl%visco%etamodel(nrsubst) == 'VS2') THEN
        ! B.A. Younglove and J.F. Ely, J. Phys. Chem. Ref. Data, 16:577-798, 1987

        ! calculation of the dilute gas contribution
        ! eq. (20)
        ! collision integral [-]

        omega_VS2 = 0.d0
        DO k=1, gl%visco%colltermnr(nrsubst)
            omega_VS2 = omega_VS2 + gl%visco%coeffcoll(nrsubst,k)*(gl%visco%eklj(nrsubst)*t**(-1.d0))**((4.d0-k)/3.d0)
        END DO
        omega_VS2 = 1.d0/omega_VS2

        ! eq. (19)
        ! dilute gas contribution
        eta0_VS2 = gl%visco%const19_VS2(nrsubst)*t**gl%visco%exp19_VS2(nrsubst)/(gl%visco%slj(nrsubst)**2.d0*omega_VS2)


        ! calculation of the moderately dense fluid contribution

        ! eq. (21)
        ! moderately dense fluid contribution
        ! unit conversion of d from [mol/m³] to [mol/dm³]
        eta1_VS2 = (d/gl%factor)*(gl%visco%Fv_VS2(nrsubst,1)+gl%visco%Fv_VS2(nrsubst,2)*(gl%visco%Fv_VS2(nrsubst,3)-log(t/gl%visco%Fv_VS2(nrsubst,4)))**2d0)


        ! calculations of the dense gas contribution

        ! eq. (25)
        ! unit conversion of d from [mol/m³] to [mol/dm³]
        Hrho_VS2 = sqrt(d/gl%factor)*(d/gl%factor-gl%visco%Ev_VS2(nrsubst,8))/gl%visco%Ev_VS2(nrsubst,8)

        ! eq. (23)
        Gt_VS2 = gl%visco%Ev_VS2(nrsubst,1)+gl%visco%Ev_VS2(nrsubst,2)/t

        ! eq. (24)
        ! unit conversion of d from [mol/m³] to [mol/dm³]
        Frhot_VS2 = Gt_VS2+(gl%visco%Ev_VS2(nrsubst,3)+gl%visco%Ev_VS2(nrsubst,4)*t**(-1.5d0))*(d/gl%factor)**0.1d0&
            &           +(gl%visco%Ev_VS2(nrsubst,5)+(gl%visco%Ev_VS2(nrsubst,6)*t**(-1.d0))+(gl%visco%Ev_VS2(nrsubst,7)*t**(-2.d0)))*Hrho_VS2

        ! eq. (22)
        ! dense gas contribution
        eta2_VS2 = exp(Frhot_VS2)-exp(Gt_VS2)

        ! eq. (18)
        ! sum of the dilute gas, moderately dense fluid and dense gas contribution
        eta_VS2 = eta0_VS2 + eta1_VS2 + eta2_VS2

        VISDYN_CALC = eta_VS2
        return

        !--------------------------------------------------------------------------------------------------------
        ! Calculating viscosity for [V]iscosity [S]pecification model [4]
        !--------------------------------------------------------------------------------------------------------
    ELSE IF (gl%visco%etamodel(nrsubst) == 'VS4') THEN
        ! does not completly match the Refprop values due to tau and Tr

        tau = gl%tc(nrsubst) / T   ! Refprop uses critical temperature from EOS
        Tr = T / gl%tc(nrsubst)    ! Refprop uses reducing temperature from viscosity model

        tau3 = tau*tau*tau

        eta0_lc = gl%visco%d0_vs4(nrsubst,1)  + gl%visco%d0_vs4(nrsubst,2) *Tr**gl%visco%d0exp(nrsubst,2) + gl%visco%d0_vs4(nrsubst,3) *Tr**gl%visco%d0exp(nrsubst,3) &
            + gl%visco%d0_vs4(nrsubst,4) *Tr**gl%visco%d0exp(nrsubst,4) + gl%visco%d0_vs4(nrsubst,5) *Tr**gl%visco%d0exp(nrsubst,5)

        if (gl%substcasnr(nrsubst) == '2551-62-4') eta0_lc=eta0_lc*1.d-3   !SF6 => different unit


        phi1 = dexp(tau) - 1.d0
        phi2 = dexp(tau*tau) - 1.d0
        ka = (gl%visco%a_0(nrsubst)  + gl%visco%a_1(nrsubst) *phi1 + gl%visco%a_2(nrsubst) *phi2) * tau
        kaa = (gl%visco%aa0(nrsubst)  + gl%visco%aa1(nrsubst) *phi1 + gl%visco%aa2(nrsubst) *phi2) * tau3
        kr = (gl%visco%b0_vs4(nrsubst)  + gl%visco%b1_vs4(nrsubst) *phi1 + gl%visco%b2_vs4(nrsubst) *phi2) * tau
        krr = (gl%visco%bb0(nrsubst)  + gl%visco%bb1(nrsubst) *phi1 + gl%visco%bb2(nrsubst) *phi2) * tau3
        ki = (gl%visco%c0(nrsubst)  + gl%visco%c_1(nrsubst) *phi1 + gl%visco%c_2(nrsubst) *phi2) * tau
        kii = (gl%visco%cc0(nrsubst)  + gl%visco%cc1(nrsubst) *phi1 + gl%visco%cc2(nrsubst) *phi2) * tau3
        kttt = (gl%visco%dd_0(nrsubst)  + gl%visco%dd_1(nrsubst) *phi1 + gl%visco%dd_2(nrsubst) *phi2) * tau
        kaaa = (gl%visco%e_0(nrsubst)  + gl%visco%e_1(nrsubst) *phi1 + gl%visco%e_2(nrsubst) *phi2) * tau

        p = P_CALC(gl,T,D,nrsubst) * 10.d0      !pressure in bar
        pid = D*gl%REQ(nrsubst)*T*1.d-6 * 10.d0
        DPDT = DPDT_CALC(gl,T,D,nrsubst) * 10.d0
        pa = p - T * DPDT
        pr = p - pa - pid

        eta4rs=(ki*pid + kr*pr + ka*pa + kii*pid*pid + krr*pr*pr + kaa*pa*pa + kttt*(pr+pid)**3 + kaaa * pa**3) * gl%factor
        VISDYN_CALC = eta0_lc + eta4rs

        return

        !--------------------------------------------------------------------------------------------------------
        ! Calculating viscosity for [V]iscosity [S]pecification model [5]
        ! M. Thol and M. Schaden, July 2015
        !--------------------------------------------------------------------------------------------------------
    ELSE IF (gl%visco%etamodel(nrsubst) == 'VS5') THEN

        !Parameters VS5 model
        G_VS5=-6.435d-4
        H_VS5=7.27371d0
        S_VS5=18.0323d0
        W_VS5=-0.7683d0

        a0_VS5(1)=6.32402d0
        a0_VS5(2)=0.12102d-2
        a0_VS5(3)=5.28346d0
        a0_VS5(4)=6.62263d0
        a0_VS5(5)=19.74540d0
        a0_VS5(6)=-1.89992d0
        a0_VS5(7)=24.27450d0
        a0_VS5(8)=0.79716d0
        a0_VS5(9)=-0.23816d0
        a0_VS5(10)=0.68629d-1

        a1_VS5(1)=50.41190d0
        a1_VS5(2)=-0.11536d-2
        a1_VS5(3)=254.209d0
        a1_VS5(4)=38.0957d0
        a1_VS5(5)=7.63034d0
        a1_VS5(6)=-12.5367d0
        a1_VS5(7)=3.44945d0
        a1_VS5(8)=1.11764d0
        a1_VS5(9)=0.67695d-1
        a1_VS5(10)=0.34793d0

        a2_VS5(1)=-51.6801d0
        a2_VS5(2)=-0.62571d-2
        a2_VS5(3)=-168.481d0
        a2_VS5(4)=-8.46414d0
        a2_VS5(5)=-14.3544d0
        a2_VS5(6)=4.98529d0
        a2_VS5(7)=-11.2913d0
        a2_VS5(8)=0.12348d-1
        a2_VS5(9)=-0.8163d0
        a2_VS5(10)=0.59256d0

        a3_VS5(1)=1189.02d0
        a3_VS5(2)=0.37283d-1
        a3_VS5(3)=3898.27d0
        a3_VS5(4)=31.4178d0
        a3_VS5(5)=31.5267d0
        a3_VS5(6)=-18.1507d0
        a3_VS5(7)=69.3466d0
        a3_VS5(8)=-4.11661d0
        a3_VS5(9)=4.02528d0
        a3_VS5(10)=-0.72663d0

        !! T-H. Chung, M. Ajlan, L.L. Lee and K.E. Starling, Ind. Eng. Chem. Res. 1998, 27, 671-679.

        ! eq. (3)
        ! dimensionless temperature [-]
        tstar = t/gl%visco%eklj(nrsubst)

        ! eq. (4) solved to critical volume for eq. (6) [cm³/g]
        vc_VS5 = (gl%visco%slj_VS5(nrsubst)/0.809d0)**3

        eta0_VS5 = eta0 (gl,T,nrsubst)

        ! eq. (11)
        ! nonpolar substances only the first two terms
        ! for polar substances the third term is included
        ! for hydrogen bonding substances the last term is included
        DO i=1, 10
            Ai_VS5(i) = a0_VS5(i) + a1_VS5(i)*gl%visco%accen_VS5(nrsubst) + a2_VS5(i)*gl%visco%dipolered_VS5(nrsubst)**4.d0 + a3_VS5(i)*gl%visco%kappa_VS5(nrsubst)
        END DO

        ! parameters for eq. (10) [-]
        Y_VS5 = d/gl%factortrans*vc_VS5/6.d0      !unit conversion of density d from [mol/m³] to [mol/cm³]
        G1_VS5 = (1.d0-0.5d0*Y_VS5)/(1.d0-Y_VS5)**3.d0
        G2_VS5 = (Ai_VS5(1)*(1.d0-dexp(-Ai_VS5(4)*Y_VS5))/Y_VS5+Ai_VS5(2)*G1_VS5*dexp(Ai_VS5(5)*Y_VS5)+Ai_VS5(3)*G1_VS5)/(Ai_VS5(1)*Ai_VS5(4)+Ai_VS5(2)+Ai_VS5(3))

        ! dense fluid contributions of the viscosity [g/(cm*s)]
        etak = eta0_VS5*(1.d0/G2_VS5+Ai_VS5(6)*Y_VS5)
        etap = ((36.344d-6)*dsqrt(gl%wm(nrsubst)*gl%factor*gl%visco%tc_vs5(nrsubst))/vc_VS5**(2.d0/3.d0))*Ai_VS5(7)*Y_VS5**2.d0*G2_VS5*dexp(Ai_VS5(8)+Ai_VS5(9)/tstar+Ai_VS5(10)/tstar**2.d0)
        ! at the limiting value of low density, Y_VS5 approaches zero and G2_VS5
        ! aproaches unity, so that etap is negligible and eq.(10) reduces to eq.(6)
        ! at extremly high density etap is the major contributing term

        ! eq. (10)
        eta_VS5 = (etak + etap)*gl%factor_VS5eta  !unit conversion from [g/(cm*s)] to [uPa*s]

        VISDYN_CALC = eta_VS5

        return

    ELSE IF(gl%visco%etamodel(nrsubst) == 'VS6') THEN
        ! Calculating viscosity for [V]iscosity [S]pecification model [6]
        ! S. Pohl, Nov. 2019
        !--------------------------------------------------------------------------------------------------------

        ! REF FOR PROGRAMMING (CYCLOHEXANE)
        !Tariq, U., Jusoh, A.R.B., Riesco, N., and Vesovic, V.,
        ! "Reference Correlation of the Viscosity of Cyclohexane from the Triple Point to 700 K and up to 110 MPa,"
        !J. Phys. Chem. Ref. Data, 43, 033101, 2014.




























        !--------------------------------------------------------------------------------------------------------
        ! Calculating viscosity for [V]iscosity [S]pecification model [9]
        ! T. Eisenbach and M. Thol, March 2017; February 2020
        !--------------------------------------------------------------------------------------------------------
    ELSE IF (gl%visco%etamodel(nrsubst) == 'VS9') THEN
        !cf. Vogel, Span, Herrmann, J. Phys. Chem. Ref. Data 44, 043101 (2015)

        !reduced temperature
        tau_VS9=gl%visco%tredetadg(nrsubst)/T
        !reduced density
        sigma_VS9= D / gl%visco%rhoredeta(nrsubst)

        sum1_VS9 = 0.d0
        sum2_VS9 = 0.d0
        sum3_VS9 = 0.d0
        sum3par1_VS9=0.d0
        sum3par2_VS9=0.d0
        ini_dens = 0.d0
        e_VS9 = 0.d0
        eta1 = 0.d0
        tred = 0.d0
        rhored = 0.d0
        etacrit = 0.d0

        sum_comp_exp = 0.d0
        sum_poly = 0.d0
        sum_exp = 0.d0
        sum_crit = 0.d0
        nterm_normal = 0

        eta0_lc = eta0(gl,T,nrsubst)

        !initial density part
        if (gl%visco%ndens_initial(nrsubst) .ne. 0) then
            sum_ini_dens = 0.d0
            if (gl%substcasnr(nrsubst)=="106-97-8") then !-----------------butane---------------------
                do i=2, gl%visco%ndens_initial(nrsubst)+1
                    sum_ini_dens = sum_ini_dens + gl%visco%cieta(nrsubst,i)*tau_VS9**gl%visco%tieta(nrsubst,i)
                end do
                ini_dens=1.d0+gl%visco%cieta(nrsubst,1) * sum_ini_dens * sigma_VS9
                eta0_lc = eta0_lc * ini_dens
            elseif (gl%substcasnr(nrsubst)=="124-38-9") then    !-------------------CO2-----------------------
                tstar = T / gl%visco%tinit_red(nrsubst)
                ini_dens = 0.d0
                do i = 1, gl%visco%ndens_initial(nrsubst) - 1
                    ini_dens = ini_dens + gl%visco%cieta(nrsubst,i) * tstar ** gl%visco%tieta(nrsubst,i)
                end do
                !Eq. (5), page 11
                !Monika: Hier eventuell noch etwas falsch??
                eta1 = eta0_lc * ini_dens * gl%visco%cieta(nrsubst,10) / gl%wm(nrsubst)
            elseif ((gl%substcasnr(nrsubst)=="1333-74-0") .or. (gl%substcasnr(nrsubst)=="7782-39-0")) then !------------hydrogen or deuterium---------------
                tstar = T / gl%visco%tinit_red(nrsubst)
                ini_dens = 0.d0
                do i = 1, gl%visco%ndens_initial(nrsubst) - 1
                    ini_dens = ini_dens + gl%visco%cieta(nrsubst,i) * tstar ** gl%visco%tieta(nrsubst,i)
                end do
                eta1 = ini_dens * gl%visco%etainit_red(nrsubst) * gl%visco%cieta(nrsubst,gl%visco%ndens_initial(nrsubst))**3 * D/gl%factor * eta0_lc
            end if
        end if


        !calculation of the residual part
        if (gl%substcasnr(nrsubst)=="1120-21-4") then           !-------------C11----------------------------
            etacorr_VS9=0.d0
            etacorr_VS9=eta0_lc+(sigma_VS9**(2.d0/3.d0))*((1/tau_VS9)**(1.d0/2.d0))*gl%visco%coeff(nrsubst,1)/(gl%visco%coeff(nrsubst,2)+gl%visco%coeff(nrsubst,3)*(1/tau_VS9)+sigma_VS9**2+(1/tau_VS9)**2+gl%visco%coeff(nrsubst,4)*sigma_VS9*(1/tau_VS9)+gl%visco%coeff(nrsubst,5)*sigma_VS9)
        elseif (gl%substcasnr(nrsubst)=="124-38-9") then        !------------co2--------------------
            !eq.(8), Laesecke et al. (2017), page 14
            tr = T/gl%visco%tredeta(nrsubst)
            sum_poly = (gl%visco%visredeta(nrsubst) /gl%factor * (gl%visco%coeff(nrsubst,5) * Tr * sigma_VS9 **3.d0 + (sigma_VS9 ** 2.d0 + sigma_VS9 ** gl%visco%dexpo(nrsubst,2)) / (Tr + gl%visco%coeff(nrsubst,4))))
            sum_exp = sigma_VS9 ** 2.d0 + sigma_VS9 ** gl%visco%dexpo(nrsubst,2)
            !Eq. (1), Laesecke et al. (2017), page 8
            DSI = D * gl%wm(nrsubst)  !specific density in kg/m³
            etacorr_VS9 = (eta0_lc + DSI*eta1 + sum_poly)*gl%factor
        elseif ((gl%substcasnr(nrsubst)=="1333-74-0") .or. (gl%substcasnr(nrsubst)=="7782-39-0")) then !-----------hydrogen or deuterium--------------
            !Eq. (9)
            tred = T/gl%visco%tredeta(nrsubst)
            rhored = D / gl%visco%rhoredeta(nrsubst)
            term1 = gl%visco%coeff(nrsubst,1)*tred**gl%visco%texpo(nrsubst,1) * rhored**gl%visco%dexpo(nrsubst,1)
            term2 = gl%visco%coeff(nrsubst,2)*tred**gl%visco%texpo(nrsubst,2) * rhored**gl%visco%dexpo(nrsubst,2)
            term3 = gl%visco%coeff(nrsubst,3)*tred**gl%visco%texpo(nrsubst,3) * rhored**gl%visco%dexpo(nrsubst,3)
            term4 = gl%visco%coeff(nrsubst,4)*tred**gl%visco%texpo(nrsubst,4) * rhored**gl%visco%dexpo(nrsubst,4)
            term7 = gl%visco%coeff(nrsubst,7)*tred**gl%visco%texpo(nrsubst,7) * rhored**gl%visco%dexpo(nrsubst,7)
            etacrit = term7 * dexp(term1 + term2 + term4/(gl%visco%coeff(nrsubst,5) + tred) + term3)*gl%factor
            if (gl%substcasnr(nrsubst)=="7782-39-0") etacrit = etacrit * 1.18d0

            etacorr_VS9 = eta0_lc + eta1 + etacrit

        else

            !first sum in final viscosity correlation, polynomial terms in eq. (45)
            do i = 1,gl%visco%term_num1_eta(nrsubst) - gl%visco%term_expo_eta(nrsubst) - gl%visco%term_com_expo_eta(nrsubst)
                sum_poly=sum_poly + gl%visco%coeff(nrsubst,i)*tau_VS9**(gl%visco%texpo(nrsubst,i))*sigma_VS9**(gl%visco%dexpo(nrsubst,i))
            end do

            !second sum in final viscosity correlation, exponential terms in eq. (45)
            do i =gl%visco%term_num1_eta(nrsubst)-gl%visco%term_expo_eta(nrsubst) + 1,gl%visco%term_num1_eta(nrsubst)
                sum_exp=sum_exp+gl%visco%coeff(nrsubst,i)*tau_VS9**(gl%visco%texpo(nrsubst,i))*sigma_VS9**(gl%visco%dexpo(nrsubst,i))*dexp(-sigma_VS9)
            end do

            nterm_normal =  gl%visco%term_num1_eta(nrsubst) - gl%visco%term_com_expo_eta(nrsubst)
            !third sum in final viscosity correlation, complex exponential terms in eq. (45)
            do i =1,gl%visco%term_com_expo_eta(nrsubst)
                sum_comp_exp=sum_comp_exp+gl%visco%coeff(nrsubst,i+nterm_normal)*tau_VS9**(gl%visco%texpo(nrsubst,i+nterm_normal))*sigma_VS9**(gl%visco%dexpo(nrsubst,i+nterm_normal))*(sigma_VS9**gl%visco%li_VS9(nrsubst,i+nterm_normal) * tau_VS9) ** gl%visco%pi_VS9(nrsubst,i+nterm_normal)
            end do

            sum3par1_VS9=(sigma_VS9-1.d0)**2
            sum3par2_VS9=abs(tau_VS9-1.d0)

            do i=gl%visco%term_num1_eta(nrsubst)+1,gl%visco%term_num1_eta(nrsubst)+gl%visco%term_num2_eta(nrsubst)
                sum_crit=sum_crit+(gl%visco%coeff(nrsubst,i)*tau_VS9*sigma_VS9*dexp((1.d0*gl%visco%beta_VS9(nrsubst,i)*sum3par1_VS9)+ (1.d0*gl%visco%eta_VS9(nrsubst,i)*sum3par2_VS9)))
            end do

            !final viscosity correlation eq. (45)
            etacorr_VS9=eta0_lc+eta1+sum_poly + sum_exp + sum_comp_exp + sum_crit
        end if

        VISDYN_CALC=etacorr_VS9!*gl%factor

        return

    ELSEIF (gl%substcasnr(nrsubst) =="7722-84-1") THEN
        !Viscosity model for hydrogen peroxide according to Yaws_b1999_Chemical Properties Handbook

        if (trim(gl%inptorig) .ne. "td") then
            VISDYN_CALC = -9955.d0
            return
        end if

        !Parameter for liquid model
        A_h2o2_liq = -1.615d0
        B_h2o2_liq = 503.8d0
        C_h2o2_liq = 0.0003501d0
        D_h2o2_liq = -0.000001168d0

        !Parameter for vapor model
        A_h2o2_vap = 8.039d0
        B_h2o2_vap = 0.27d0
        C_h2o2_vap = 0.0000829d0

        d_crit_h2o2 = 437.8d0/gl%wm(1)

        if (D .ge. d_crit_h2o2) then    !liquid phase
            if (t .lt. 273.d0) then
                VISDYN_CALC = -9912.d0
                return
            elseif (t .gt. 728.d0) then
                VISDYN_CALC = -9913.d0
                return
            else
                VISDYN_CALC = 10.d0 ** (A_h2o2_liq + B_h2o2_liq / T + C_h2o2_liq * T + D_h2o2_liq * T**2) * gl%factor
            end if
        else
            if (t .lt. 373.d0) then
                VISDYN_CALC = -9912.d0
                return
            elseif (t .gt. 600.d0) then
                VISDYN_CALC = -9913.d0
                return
            else
                VISDYN_CALC = A_h2o2_vap + B_h2o2_vap * T + C_h2o2_vap * T**2
            end if
        end if
    ELSE
        VISDYN_CALC = -5234.d0
    ENDIF

    !Andreas, October 2013, ERRORHANDLING
    !If the viscosity is 0, give back an error
    if (VISDYN_CALC .lt. 1.D-14) then
        VISDYN_CALC = -5243.D0
    end if

    end function




    !-------------------------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------------------------
    double precision function eta0 (gl,T,nrsubst)

    implicit none

    type(type_gl) :: gl

    double precision:: eta0_lc, t, taut, taut12, taut13, denom, tstar, tau_vs9, e_vs9, eta0_vs9, s_star, ci, collsum
    integer:: i, nrsubst
    !parameter for VS5
    double precision:: vc_VS5, omegastar, fc_VS5, A_VS5, B_VS5, C_VS5, D_VS5, E_VS5, F_VS5

    !-------------------------------------------------------------------------------------
    !-------------------------------VS1 model---------------------------------------------
    !-------------------------------------------------------------------------------------
    IF (gl%visco%etamodel(nrsubst) == 'VS1') THEN

        ci=0.d0
        collsum=0.d0
        eta0 = 0.d0

        tstar=t/gl%visco%eklj(nrsubst)

        !eta0vs1  < --dilute gas viscosity
        IF (((gl%visco%pointereta(nrsubst) == 'CI1') .AND. (gl%substcasnr(nrsubst) == '64-17-5')) &              !ethanol
            & .OR.((gl%visco%pointereta(nrsubst) == 'NUL') .AND. (gl%substcasnr(nrsubst) == '306-83-2'))) THEN   !R123
            eta0=gl%visco%a0_vs1(nrsubst)+gl%visco%a1_vs1(nrsubst)*t+gl%visco%a2_vs1(nrsubst)*t**2

        ELSEIF ((gl%visco%pointereta(nrsubst) == 'CI0')) THEN
            IF ((gl%substcasnr(nrsubst) == '354-33-6').or.(gl%substcasnr(nrsubst) == '616-38-6')) THEN      !R125 or DMC
                ci=0.d0
                !ci=1.16145d0/tstar**0.14874d0+0.52487d0*dexp(-0.77320d0*tstar)+2.16178d0*dexp(-2.43780d0*tstar)   !Paper: M. L. Huber, A. Laesecke, Ind. Eng. Chem. Res., 45(12):4447, 2006
                ci=1.16145d0/tstar**0.14874d0+0.52487d0*DEXP(-0.77320d0*tstar)+2.16178d0*DEXP(-2.43787d0*tstar)   !REFPROP 9.1, Source Code
            ELSEIF (gl%substcasnr(nrsubst) == '754-12-1') THEN      !R1234yf approximation of collision integral
                eta0 = (-836950.d0+6336.28d0*t-2.3547d0*t**2+0.0395563d0*t**3)/(39509.1d0+121.018d0*t+t**2)
            ELSEIF (gl%substcasnr(nrsubst) == '29118-24-9') THEN      !R1234zeE  approximation of collision integral
                eta0 = (-963382.d0+9614.09d0*t-13.233*t**2+0.0360562d0*t**3)/(122059.d0-224.741d0*t+t**2)
            ELSE
                ci=0.d0
                !Paper: Neufeld, P.D., Janzen, A.R., Aziz, R.A., 1972. Empirical equations to calculate 16 of the transport collision integrals omega(l.s)* for the Lennard-Jones (12-6) potential. J. Chem. Phys. 57, 1100–1102.
                !----------------------
                !Difference to REFPROP: this part is commented in REFPROP and, therefore, not evaluated:             - (6.435d-4*tstar**0.14874d0) * (sin(18.0323d0*tstar**-0.7683d0-7.27371d0))
                !----------------------
                ci=1.16145d0/tstar**0.14874d0 + 0.52487d0*DEXP(-0.77320d0*tstar) + 2.16178d0*DEXP(-2.43787d0*tstar) ! - (6.435d-4*tstar**0.14874d0) * (sin(18.0323d0*tstar**-0.7683d0-7.27371d0))
            ENDIF

        ELSEIF ((gl%visco%pointereta(nrsubst) == 'CI1')) THEN

            collsum=0.d0
            DO i = 1,gl%visco%colltermnr(nrsubst)
                collsum=collsum+gl%visco%coeffcoll(nrsubst,i)*dlog(tstar)**gl%visco%powercoll(nrsubst,i)
            END DO
            ci=dexp(collsum)

        ELSEIF (gl%visco%pointereta(nrsubst) == 'CI2') THEN

            collsum=0.d0
            DO i = 1,gl%visco%colltermnr(nrsubst),1
                collsum=collsum+gl%visco%coeffcoll(nrsubst,i)*tstar**((i-1.d0)/3.d0-1.d0)
            END DO
            ci=1.d0/collsum
        ENDIF

        !if((trim(gl%visco%pointereta(nrsubst)) /= 'NUL').and.((trim(gl%substcasnr(nrsubst)) /= '64-17-5') .or. (trim(gl%substcasnr(nrsubst)) == '102687-65-0'))) then  !ethanol and r1233zd(E)
        if (eta0 .lt. 1.d-14) then
            eta0=gl%visco%ce(nrsubst)*t**gl%visco%cep(nrsubst)/gl%visco%slj(nrsubst)**2/ci
        end if

        eta0=eta0*gl%visco%visredetadg(nrsubst)

        !-------------------------------------------------------------------------------------
        !-------------------------------VS5 model---------------------------------------------
        !-------------------------------------------------------------------------------------
    ELSEIF (gl%visco%etamodel(nrsubst) == 'VS5') THEN

        ! T-H. Chung, M. Ajlan, L.L. Lee and K.E. Starling, Ind. Eng. Chem. Res. 1998, 27, 671-679.

        A_VS5=1.16145d0
        B_VS5=0.14874d0
        C_VS5=0.52487d0
        D_VS5=0.7732d0
        E_VS5=2.16178d0
        F_VS5=2.43787d0

        ! eq. (3)
        ! dimensionless temperature [-]
        tstar = t/gl%visco%eklj(nrsubst)

        ! eq. (2)
        ! reduced collision integral for eq. (6) [-]
        omegastar = A_VS5/(tstar**B_VS5)+C_VS5*dexp(-D_VS5*tstar)+E_VS5*dexp(-F_VS5*tstar) !+G_VS5*(tstar**B_VS5)*sin(S_VS5*(tstar**W_VS5)-H_VS5)

        ! eq. (4) solved to critical volume for eq. (6) [cm³/g]
        vc_VS5 = (gl%visco%slj_VS5(nrsubst)/0.809d0)**3

        ! eq. (7)
        ! factor Fc to account for molecular structure and polar effects [-]
        fc_VS5 = 1.d0-0.2756d0*gl%visco%accen_VS5(nrsubst)+0.059035d0*gl%visco%dipolered_VS5(nrsubst)**4.d0+gl%visco%kappa_VS5(nrsubst)

        ! eq. (6)
        ! Dilute gas contribution of the viscosity [g/(cm*s)]
        ! unit conversion of molar weight from [kg/mol] to [g/mol]
        eta0 = 4.0785d-5*dsqrt(gl%wm(nrsubst)*gl%factor*t)*fc_VS5/(vc_VS5**(2.d0/3.d0)*omegastar)     !This line causes the tiny difference when comparing to REFPROP. In REFPROP, a parameter is used instead of calculating the entire equation.

        continue


        !-------------------------------------------------------------------------------------
        !-------------------------------VS9 model---------------------------------------------
        !-------------------------------------------------------------------------------------
    ELSEIF (gl%visco%etamodel(nrsubst) == 'VS9') THEN

        eta0_lc = 0.d0

        !reduced temperature
        tau_VS9=gl%visco%tredetadg(nrsubst)/T

        !calculate dilute-gas contribution
        if (gl%substcasnr(nrsubst)=="74-84-0") then !ethane

            if (gl%visco%dilgas_term(nrsubst) .ne. 0) then
                taut = t
                taut12 = dsqrt(t)
                taut13 = t**(1.d0/3.d0)
                denom = gl%visco%coeff_hi(nrsubst,1) + gl%visco%coeff_hi(nrsubst,2) / dexp(taut13)   &
                    &               + (gl%visco%coeff_hi(nrsubst,3) + gl%visco%coeff_hi(nrsubst,4) / dexp(taut13)) / taut12   &
                    &               + gl%visco%coeff_hi(nrsubst,5) * t / dexp(2.d0 * taut13)
                eta0_lc = taut12 / denom
            end if

        else if (gl%substcasnr(nrsubst)=="1120-21-4") then      !-----------------C11---------------------
            eta0_lc = (gl%visco%coeff_hi(nrsubst,1)+(gl%visco%coeff_hi(nrsubst,2)*(1/tau_VS9)) + (gl%visco%coeff_hi(nrsubst,3)*((1/tau_VS9)**2)) + (gl%visco%coeff_hi(nrsubst,4)*((1/tau_VS9)**3)) + (gl%visco%coeff_hi(nrsubst,5)*((1/tau_VS9)**4)) + (gl%visco%coeff_hi(nrsubst,6)*((1/tau_VS9)**5)))/(gl%visco%coeff_hi(nrsubst,7)+(gl%visco%coeff_hi(nrsubst,8)*(1/tau_VS9)))

        else if (gl%substcasnr(nrsubst)=="106-97-8") then       !-----------------BUTANE------------------
            e_VS9 = gl%visco%coeff_hi(nrsubst,2) + gl%visco%coeff_hi(nrsubst,3) * dlog(tau_VS9) + gl%visco%coeff_hi(nrsubst,4) * (dlog(tau_VS9)) ** 2
            eta0_lc = gl%visco%coeff_hi(nrsubst,1)/(tau_VS9**0.5d0 * dexp(e_VS9))

        elseif (gl%substcasnr(nrsubst)=="124-38-9") then        !------------------CO2--------------------
            !Eq. (4) page 9
            denom = gl%visco%coeff_hi(nrsubst,1) + gl%visco%coeff_hi(nrsubst,2) * T ** (1.d0/6.d0) + gl%visco%coeff_hi(nrsubst,3) * dexp(gl%visco%coeff_hi(nrsubst,4) * T ** (1.d0/3.d0)) + (gl%visco%coeff_hi(nrsubst,5) + gl%visco%coeff_hi(nrsubst,6) * T ** (1.d0/3.d0)) / dexp(T ** (1.d0/3.d0)) + gl%visco%coeff_hi(nrsubst,7) * T ** 0.5d0
            eta0_lc = gl%visco%coeff_hi(nrsubst,8) * T ** 0.5d0 / denom

        elseif ((gl%substcasnr(nrsubst)=="1333-74-0") .or. (gl%substcasnr(nrsubst)=="7782-39-0")) then   !-------------hydrogen or deuterium----------------
            s_star = 0.d0
            tstar = T/gl%visco%coeff_hi(nrsubst,gl%visco%dilgas_term(nrsubst))
            do i=1,(gl%visco%dilgas_term(nrsubst)-3)
                s_star =s_star + gl%visco%coeff_hi(nrsubst,i)*dlog(tstar)**(i-1)
            end do
            s_star = dexp(s_star)
            eta0_lc = gl%visco%visredetadg(nrsubst) * dsqrt(T*gl%visco%coeff_hi(nrsubst,(gl%visco%dilgas_term(nrsubst)-1))) / gl%visco%coeff_hi(nrsubst,(gl%visco%dilgas_term(nrsubst)-2))**2/s_star
            if (gl%substcasnr(nrsubst)=="7782-39-0") eta0_lc = eta0_lc * dsqrt(2.d0)

        else
            do i=1,gl%visco%dilgas_term(nrsubst)
                eta0_lc =eta0_lc + gl%visco%coeff_hi(nrsubst,i)*tau_VS9**(-1.d0*i)
            end do
        end if
        eta0 = eta0_lc
    ENDIF

    end function



    !----------------------------------------------------------------------------------
    !Function for calculating thermal conductivity [W/m/K]
    !----------------------------------------------------------------------------------
    ! M. Thol, June 2014

    double precision function TCX_CALC(gl,T, D, nrsubst)

    implicit none

    type(type_gl) :: gl


    double precision:: T, D, lam0, lam1, lamcrit
    integer:: nrsubst
    double precision:: A_h2o2, B_h2o2, C_h2o2, d_crit_h2o2

    lam0 = 0.d0
    lam1 = 0.d0
    lamcrit = 0.d0
    d_crit_h2o2 = 0.d0

    TCX_CALC = 0.d0

    if (gl%tcx%tcmodel(nrsubst) == 'TC0') then
        if(gl%substcasnr(nrsubst) == '75-46-7') then    !R23
            TCX_CALC=tc0R23(gl,T,D,nrsubst)
            return
        elseif(gl%substcasnr(nrsubst) == '7440-59-7') then    !Helium
            TCX_CALC=tc0helium(gl,T,D,nrsubst)
            return
        else
            lam0=tcxtc0(gl,T,D,nrsubst)
            if (gl%tcx%tkmodel(nrsubst) == 'TK1') then
                lamcrit=dlamctk1(gl,T,D,nrsubst)
            elseif (gl%tcx%tkmodel(nrsubst) == 'TK3') then
                lamcrit=dlamctk3(gl,T,D,nrsubst)
            end if
            TCX_CALC=lam0 + lamcrit
            return
        end if
    elseif (gl%tcx%tcmodel(nrsubst) == 'TC1') then
        lam0=lam0tc1(gl,T,nrsubst)
        lam1=dlamtc1(gl,T,D,nrsubst)
        if     (gl%tcx%tkmodel(nrsubst) == 'TK1') then
            lamcrit=dlamctk1(gl,T,D,nrsubst)
        elseif (gl%tcx%tkmodel(nrsubst) == 'TK3') then
            lamcrit=dlamctk3(gl,T,D,nrsubst)
            if (abs(lamcrit + 5245.d0) .lt. 1.d-12) then   !no viscosity model available to calculate the critical enhancement
                TCX_CALC = -5245.d0
                return
            end if
            !CAUTION: Slight differences to Refprop can be caused by a calculated critical pressure (Trend) instead of using the one given in the fluid file (Refprop)
            !elseif (tkmodel(nrsubst) == 'TK6') then
            !    lamcrit=dlamctk6(gl,T,D,nrsubst)
        elseif (gl%tcx%tkmodel(nrsubst) == 'CH4') then
            lamcrit=lamcrch4(gl,T,D,nrsubst)
        end if
        !MT 2019/05/15, deleted the if query, does not make sense
        !if (lamcrit .lt. 1.d-14) then
        !    TCX_CALC = lamcrit
        !else
        TCX_CALC=lam0 + lam1 + lamcrit
        !end if
        return
    elseif (gl%tcx%tcmodel(nrsubst) == 'TC2') then
        TCX_CALC = -5244.d0
        return
    elseif (gl%tcx%tcmodel(nrsubst) == 'TC3') then
        TCX_CALC=lam0tc3(gl,T,nrsubst) + lam123tc3(gl,T,D,nrsubst)
        return
    elseif (gl%tcx%tcmodel(nrsubst) == 'TC4') then
        TCX_CALC = -5244.d0
        return
    elseif (gl%tcx%tcmodel(nrsubst) == 'TC5') then
        TCX_CALC = lambda_TC5(gl,T,D,nrsubst)
        return
    elseif (gl%tcx%tcmodel(nrsubst) == 'TC6') then
        TCX_CALC = -5244.d0
        return
    elseif (gl%substcasnr(1)=="7722-84-1") then
        !thermal conductivity model for hydrogen peroxide according to Yaws_b1999_Chemical Properties Handbook

        if (trim(gl%inptorig) .ne. "td") then
            TCX_CALC = -9955.d0
            return
        end if

        d_crit_h2o2 = 437.8d0/gl%wm(1)

        if (D .ge. d_crit_h2o2) then   !liquid model parameter
            if (t .lt. 273.d0) then
                TCX_CALC = -9912.d0
                return
            elseif (t .gt. 657.d0) then
                TCX_CALC = -9913.d0
                return
            else
                A_h2o2 = 0.4425d0
                B_h2o2 = -0.00018406d0
                C_h2o2 = -0.00000038824d0
            end if
        else   !vapor model parameter
            if (t .lt. 275.d0) then
                TCX_CALC = -9912.d0
                return
            elseif (t .gt. 1200.d0) then
                TCX_CALC = -9913.d0
                return
            else
                A_h2o2 = -0.00858d0
                B_h2o2 = 0.000086933d0
                C_h2o2 = -0.000000006297d0
            end if
        end if

        TCX_CALC = (A_h2o2 + B_h2o2 * T + C_h2o2 * T**2)! / gl%factor
        return

    else
        TCX_CALC = -5234.d0
    end if

    end function TCX_CALC


    !---------------------------------------------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------------------------------------------
    !TCX, TC0
    double precision function tcxtc0(gl,T, D, nrsubst)

    implicit none

    type(type_gl) :: gl


    double precision:: tr, tau, dr, T, D, d_kg, th, tt
    double precision:: tcxlam0, tcxdellam, tcxdellamcr, tcxdellaml
    double precision:: f1, f2, f3, f4
    double precision:: be,b0_lc,c1_lc,c2_lc,ct1,ct2,cr1,cr2,cr3,dr1,d1   !D2O
    double precision, Dimension(9):: GT !ETY
    double precision, Dimension(7):: k  !ETY
    integer:: i
    integer:: nrsubst

    if (gl%tcx%pointer_hardtc(nrsubst) == 'D2O') then
        !Reference:
        !International Association for the Properties of Water and Steam,
        !"Viscosity and Thermal Conductivity of Heavy Water Substance,"
        !Physical Chemistry of Aqueous Systems: Proceedings of the 12th International Conference on the Properties of Water and Steam,
        !Orlando, Florida, September 11-16, A107-A138, 1994.

        tcxlam0 = 0.0d0
        tcxdellam = 0.0d0
        tcxdellamcr = 0.0d0
        tcxdellaml = 0.0d0

        be = gl%tcx%dellamcr_coeff(nrsubst,1)
        b0_lc = gl%tcx%dellamcr_coeff(nrsubst,2)
        c1_lc = gl%tcx%dellamcr_coeff(nrsubst,3)
        c2_lc = gl%tcx%dellamcr_coeff(nrsubst,4)
        ct1 = gl%tcx%dellamcr_coeff(nrsubst,5)
        ct2 = gl%tcx%dellamcr_coeff(nrsubst,6)
        cr1 = gl%tcx%dellamcr_coeff(nrsubst,7)
        cr2 = gl%tcx%dellamcr_coeff(nrsubst,8)
        cr3 = gl%tcx%dellamcr_coeff(nrsubst,9)
        dr1 = gl%tcx%dellamcr_coeff(nrsubst,10)
        d1 = gl%tcx%dellamcr_coeff(nrsubst,11)

        tr = T/gl%tcx%tred_tc(nrsubst)              !Eq. (B4)
        dr = D/gl%tcx%rhored_tc(nrsubst)            !Eq. (B5)
        tau = tr/(dabs(tr-1.1d0)+1.1d0)

        do i=1, gl%tcx%tcxlam0_num(nrsubst)         !Eq. (B8)
            tcxlam0 = tcxlam0 + gl%tcx%lam0_coeff(nrsubst,i)*tr**gl%tcx%lam0_exp(nrsubst,i)
        enddo

        do i=1, gl%tcx%dellam(nrsubst)              !Eq. (B9), part 1
            tcxdellam = tcxdellam + gl%tcx%dellam_coeff(nrsubst,i)*dr**gl%tcx%dellam_exp(nrsubst,i)
        enddo
        tcxdellam = tcxdellam + b0_lc*(1.d0-dexp(be*dr))   !Eq. (B9), part 2

        f1 = dexp(ct1*tr + ct2*tr**2)               !Eq. (B12)
        f2 = dexp(cr1*(dr - 1.d0)**2) + cr2*dexp(cr3*(dr - dr1)**2)     !Eq. (B13)
        f3 = 1.d0 + dexp(60.d0*(tau - 1.d0) + 20.d0)                       !Eq. (B14)
        f4 = 1.d0 + dexp(100.d0*(tau - 1.d0) + 15.d0)                      !Eq. (B15)

        tcxdellamcr = c1_lc*f1*f2*(1.d0 + f2**2*((c2_lc*f1**4)/f3 + 3.5d0*f2/f4))   !Eq. (B10)

        tcxdellaml = d1*f1**1.2d0*(1.d0 - dexp(-(dr/2.5d0)**10))

        tcxtc0 = gl%tcx%etared_tc(nrsubst) * (tcxlam0 + tcxdellam + tcxdellamcr + tcxdellaml)
        return

    elseif (gl%tcx%pointer_hardtc(nrsubst) == 'ETY') then
        GT(1) = -2.903423528d5
        GT(2) =  4.680624952d5
        GT(3) = -1.8954783215d5
        GT(4) = -4.8262235392d3
        GT(5) =  2.243409372d4
        GT(6) = -6.6206354818d3
        GT(7) =  8.9937717078d2
        GT(8) = -6.0559143718d1
        GT(9) =  1.6370306422d0
        k(1) = -1.304503323d1
        k(2) =  1.8214616599d1
        k(3) = -9.903022496d3
        k(4) =  7.420521631d2
        k(5) = -3.0083271933d-1
        k(6) =  9.6456068829d1
        k(7) =  1.350256962d4
        tcxlam0 = 0.d0
        d_kg = D*gl%wm(nrsubst)*1.d-3    ! kg/L
        th = (d_kg-0.221d0)/0.221d0 !(D-rhoc(nrsubst))/rhoc(nrsubst)    ! refprop rhoc = 0.221 kg/L  ?? - >  rhoc = 0.214 kg/L

        !!nach refprop
        !tt = T**(1.0d0/3.0d0)
        !tcxlam0 = GT(1)/T + GT(2)/tt**2 + GT(3)/tt + GT(4) + GT(5)*tt + GT(6)*tt**2 + GT(7)*T + GT(8)*tt**4 + GT(9)*tt**5

        do i=1, 9
            tcxlam0 = tcxlam0 + GT(i)*T**(-4.d0/3.d0 + i/3.d0)
        enddo

        tcxdellam = dexp(k(1)+k(4)/T)*(dexp(d_kg**0.1d0*(k(2)+k(3)/T**1.5d0) + th*d_kg**0.5d0*(k(5)+k(6)/T + k(7)/T**2)) - 1.d0)
        ! rhoc = 0.221d0 - >  tcxdellam(180K,19685mol/m³) = 164.628436040869, rhoc = rhoc(nrusbst) - >  tcxdellam(180K,19685mol/m³) = 171.037125069478
    elseif (gl%tcx%pointer_hardtc(nrsubst) == 'H2O') then
        tcxlam0 = 0.d0
        tcxdellam = 0.d0
        tcxdellaml = 0.d0

        tr = T/gl%tcx%tred_tc(nrsubst)
        dr = D/gl%tcx%rhored_tc(nrsubst) !*1.d-3

        do i=1, gl%tcx%tcxlam0_num(nrsubst)
            tcxlam0 = tcxlam0 + gl%tcx%lam0_coeff(nrsubst,i)/(tr**gl%tcx%lam0_exp(nrsubst,i))
        enddo
        tcxlam0 = tr**0.5d0/tcxlam0


        do i=1, gl%tcx%dellam(nrsubst)
            tcxdellaml = tcxdellaml + (1.0d0/tr - 1.0d0)**gl%tcx%dellam_exp(nrsubst,i)*gl%tcx%dellam_coeff(nrsubst,i)*(dr - 1.0d0)**gl%tcx%dellam_exp2(nrsubst,i)
        enddo
        tcxdellaml = exp(dr*tcxdellaml)

        tcxlam0 = tcxlam0*tcxdellaml   !*tcx_tk(nrsubst)
    else
        tcxtc0 = -5234.d0
        return
    endif
    tcxtc0 = (tcxlam0 + tcxdellam)*1.d-3 ! W m^-1 K^-1

    RETURN
    END

    !--------------------------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------------------------------------
    !lamda0, TC1
    double precision function lam0tc1(gl,T, nrsubst)

    implicit none

    type(type_gl) :: gl


    double precision:: lam0_num, lam0_den, tau, Rmix, cint, cp01, cp0_tc
    double precision:: T, vis0
    integer:: nrsubst, i

    lam0_num = 0.d0
    lam0_den = 0.d0
    cint = 0.d0
    cp01 = 0.d0
    vis0 = 0.d0
    tau = T/gl%tcx%tred_dil(nrsubst)
    cp0_tc = CP0_CALC(gl,T,nrsubst)    ! J/(mol K)
    call R_mix_calc(gl,Rmix)

    do i=1, gl%tcx%num_dil_tc(nrsubst)
        if (gl%tcx%num_tc_powerT_dil(nrsubst,i) >= -90.d0) then
            lam0_num = lam0_num + gl%tcx%num_tc_coeff_dil(nrsubst,i)*tau**gl%tcx%num_tc_powerT_dil(nrsubst,i)
        else
            if (gl%tcx%num_tc_powerT_dil(nrsubst,i) ==  -99.d0) then
                cint = cp0_tc - 2.5d0*Rmix
                cp01 = 1.d0 + gl%tcx%num_tc_coeff_dil(nrsubst,i)*cint
                if (lam0_num <= 0.d0) then
                    lam0_num = cp01
                else
                    lam0_num = lam0_num*cp01
                endif
            elseif (gl%tcx%num_tc_powerT_dil(nrsubst,i) == -97.d0) then
                vis0 = eta0 (gl,T,nrsubst)
                lam0_num = lam0_num + gl%tcx%num_tc_coeff_dil(nrsubst,i)*vis0
            elseif (gl%tcx%num_tc_powerT_dil(nrsubst,i) == -96.d0) then  !quelle?! - >  refprop
                cint = cp0_tc/Rmix-2.5d0
                vis0 = eta0 (gl,T,nrsubst)
                lam0_num = (lam0_num*cint+3.75d0)*Rmix*vis0/gl%wm(nrsubst)*1.d-3
            endif
        endif
    enddo
    !MT: Juni 2014, wegen CO2 auskommentiert
    !if ((cint > 0.d0) .and. (cp01 > 0.d0)) then
    !    lam0_num = lam0_num*cp01
    !endif

    do i=1, gl%tcx%den_dil_tc(nrsubst)
        lam0_den = lam0_den + gl%tcx%den_tc_coeff_dil(nrsubst,i)*(tau**gl%tcx%den_tc_exp_dil(nrsubst,i))
    enddo

    if (lam0_den > 0.d0) then
        lam0tc1 = (lam0_num/lam0_den)*gl%tcx%tcx_dil(nrsubst)              ! W m^-1 K^-1
    else
        lam0tc1 = lam0_num*gl%tcx%tcx_dil(nrsubst)                         ! W m^-1 K^-1
    endif

    RETURN
    END

    double precision function lam0tc3(gl,T,nrsubst)





    implicit none

    type(type_gl) :: gl


    double precision:: lam0
    double precision:: eta00
    double precision:: T
    integer:: nrsubst, i

    lam0 = 0.d0
    ! nach lit
    !do i=1, 9
    !    lam0 = lam0 + eta0(nrsubst,i)*T**((4.d0-i)/3.d0)
    !enddo

    ! nach refprop
    eta00 = (gl%tcx%lj_epskap(nrsubst)/T)**(1.d0/3.d0)
    do i=1, 9
        lam0 = lam0 + gl%tcx%eta0(nrsubst,i)*eta00**(4.d0-i)
    enddo

    lam0tc3 = gl%tcx%const_eq20(nrsubst) * T**gl%tcx%Texp_eq20(nrsubst) * lam0/(gl%tcx%lj_sigma(nrsubst)**2)

    RETURN
    END

    double precision function lam123tc3(gl,T,D,nrsubst)


    implicit none

    type(type_gl) :: gl


    double precision:: lam1, lam2, dellam1, dellamcrit
    double precision:: T, D, del, G, F, H, l, y, vis, del_kg, na, PI_lc, bk
    double precision:: dpdd, dpdt
    integer:: nrsubst

    lam1 = 0.d0
    lam2 = 0.d0
    del = D*1.d-3
    del_kg = D*gl%wm(nrsubst)*1.d-3
    vis = VISDYN_CALC(gl,T,D,nrsubst)
    na = 6.0221367d23
    PI_lc = 3.14159265358979d0
    bk = 1.38054d-16    !?? - >  refprop
    dpdt = DPDT_CALC(gl,T,D,nrsubst)   ! MPa/K
    dpdt = dpdt*1.d7    ! Pa/K
    dpdd = DPDD_CALC(gl,T,D,nrsubst)   ! MPa*m^3/mol
    dpdd = dpdd/gl%wm(nrsubst)*1.d10    ! mPa*m^3/kg ! refprop

    ! lam1
    lam1 = del*(gl%tcx%Fvi(nrsubst,1) + gl%tcx%Fvi(nrsubst,2)*(gl%tcx%Fvi(nrsubst,3) - dlog(T/gl%tcx%Fvi(nrsubst,4)))**2)

    ! lam2
    G = gl%tcx%Evi(nrsubst,1) + gl%tcx%Evi(nrsubst,2)/T
    H = del**0.5d0*(del - gl%tcx%Evi(nrsubst,8))/gl%tcx%Evi(nrsubst,8)
    F = G + (gl%tcx%Evi(nrsubst,3) + gl%tcx%Evi(nrsubst,4)*T**(-1.5d0))*del**0.1d0 + (gl%tcx%Evi(nrsubst,5) + gl%tcx%Evi(nrsubst,6)/T + gl%tcx%Evi(nrsubst,7)/(T**2))*H
    lam2 = dexp(F) - dexp(G)

    ! dellamcrit
    l = gl%tcx%F_tc3(nrsubst)*(gl%tcx%rm(nrsubst)**5.d0*del_kg*na/(gl%wm(nrsubst)*1000.d0)*gl%tcx%lj_epskap(nrsubst)/T)**0.5d0
    ! l = F_tc3(nrsubst)*rm**5.d0*del_kg*(na/wm(nrsubst)*lf_epskap(nrsubst)/T**0.5d0  ! lit (aber eine geschlossene Klammer fehlt)
    y = 6.d0*PI_lc*vis*1.d-5*l*(bk*T*del_kg*na/gl%wm(nrsubst)*0.1d-3)**0.5d0  ! bk =! boltz_lc?
    dellam1 = 0.d0
    if (dpdd >= 0) then
        dellam1 = bk*(T*dpdt)**2/(del_kg*dpdd)**0.5d0/y  !refprop
        ! dellam1 = K*T**2.d0*dpdt*Kt**0.5d0*y   ! lit, K=?, Kt berechnen
    endif
    dellamcrit = dellam1*dexp(-4.25d0*((D-gl%rhoc(nrsubst))/gl%rhoc(nrsubst))**4 - 18.66d0*((T-gl%tc(nrsubst))/gl%tc(nrsubst))**2)*1.0d-5 !refprop
    ! dellamcrit = dellam1*exp(-18.66*((D-rhoc(nrsubst))/rhoc(nrsubst))**4.d0 - 4.25d0*((T-tc(nrsubst))/tc(nrsubst))**2.d0)/1.0d5 ! lit

    lam123tc3 = lam1 + lam2 +  dellamcrit
    RETURN
    END

    ! DeltaLamda, TC1
    double precision function dlamtc1(gl,T,D,nrsubst)

    implicit none

    type(type_gl) :: gl


    double precision:: Deltalam_num, del, tau, ts, Rmix
    double precision:: T, D, rhoexp, delstar
    double precision:: rhosat_vap, rhosat_vap2, psat, psats, zc, beta_lc
    double precision, dimension(5), parameter:: j = (/-0.7377483d0,-1.241532d0,-1.649972d0,2.281949d0,1.43957d0/)
    double precision, dimension(5), parameter:: h = (/-6.589879d0,0.6355175d0,11.31028d0,-10.3872d0,3.393075d0/)
    integer:: i
    integer:: nrsubst     !current substance (number)

    Deltalam_num = 0.d0
    rhosat_vap = 0.d0
    rhosat_vap2 = 0.d0
    psat = 0.d0
    psats = 0.d0
    zc = 0.d0
    beta_lc = 0.355d0
    delstar = 11.d0
    del = D/(1000.d0*gl%tcx%rhored_bgrd(nrsubst))     ! mol/m³ - >  mol/L - >  kg/m³ or [-]
    tau = T/gl%tcx%tred_bgrd(nrsubst)
    ts = 1 - T/gl%tc(nrsubst)
    call R_mix_calc(gl,Rmix)

    do i=1, gl%tcx%num_bgrd_tc(nrsubst)
        if (gl%tcx%rho_exp_bgrd(nrsubst,i) == -99) then
            if ((T < gl%tc(nrsubst)) .and. (D < gl%rhoc(nrsubst))) then
                psat = gl%pc(nrsubst)*exp(h(1)*ts/(1-ts) + h(2)*ts + h(3)*ts**(1.9) + h(4)*ts**2 + h(5)*ts**3)
                psats = psat/gl%pc(nrsubst)

                zc = gl%pc(nrsubst)*10.d6/(Rmix*gl%tc(nrsubst)*gl%rhoc(nrsubst))

                rhosat_vap = psat*10.d6/(Rmix*T)*(1.d0+psat*tau**8*(zc-1)/gl%pc(nrsubst)*(1.d0+(j(1)*ts**beta_lc+j(2)*ts**(2*beta_lc)+j(3)*(ts+ts**4)+j(4)*ts**2)/(1.d0+j(5)*ts)))**(-1)
                rhosat_vap2 = gl%rhoc(nrsubst)*(1.d0-ts)**7*(1.d0-(1.d0-(1.d0-ts)**8/psats)/zc+(1-zc**(-1))*(j(1)*ts**beta_lc+j(2)*ts**(2*beta_lc)+j(3)*(ts+ts**4)+j(4)*ts**2)/(1.d0+j(5)*ts))**(-1)

                delstar = rhosat_vap/gl%rhoc(nrsubst)
                !delstar = rhosat_vap2/rhoc(nrsubst)
            else
                delstar = 11.d0     !refprop = 1?!
            endif
            Deltalam_num = Deltalam_num/delstar
        else
            if (gl%tcx%rho_exp_bgrd(nrsubst,i) > 0.d0) then
                rhoexp = exp(-del**(gl%tcx%rho_exp_bgrd(nrsubst,i)))
            else
                rhoexp = 1.d0
            endif
            Deltalam_num = Deltalam_num + gl%tcx%num_tc_coeff_bgrd(nrsubst,i)*tau**gl%tcx%num_tc_powerT_bgrd(nrsubst,i)*del**gl%tcx%num_tc_rho_bgrd(nrsubst,i)*rhoexp
        endif
    enddo

    dlamtc1 = Deltalam_num *gl%tcx%tcx_bgrd(nrsubst)                 ! W m^-1 K^-1
    RETURN
    END

    ! Lamda Crit, CH4
    double precision function lamcrch4(gl,T,D,nrsubst)

    implicit none

    type(type_gl) :: gl


    double precision:: T, D, del, tau, ts, rhos, Rmix
    double precision:: xt,f,vis, th, omega_lc
    double precision::  dpdt
    double precision:: dpdd

    integer:: nrsubst     !current substance (number)

    del = D/gl%rhoc(nrsubst)
    tau = gl%tc(nrsubst)/T
    rhos = 1.d0 - D/gl%rhoc(nrsubst)
    ts = 1.d0 - T/gl%tc(nrsubst)
    dpdt = 0.d0
    dpdd = 0.d0
    th = 0.d0
    f = 0.d0
    omega_lc = 0.d0
    vis = VISDYN_CALC(gl,T,D,nrsubst)

    call R_mix_calc(gl,Rmix)

    f = exp(-(2.646d0*abs(ts)**0.5d0+2.678d0*rhos**2-0.637d0*rhos))

    if ((abs(ts) < 0.03d0) .and. (abs(rhos) < 0.25d0)) then   ! 185K <= T <= 196K, 7.6 mol dm^-3 <= rho <= 12.7 mol dm^-3
        if (ts < (-abs(rhos)**(1.d0/0.355d0)/(-6.098d0))) then
            th = 1.d0 + 0.287d0*(1 + (-6.098d0)*ts*abs(rhos)**(-1.d0/0.355d0))**(2.d0*0.355d0)
        else
            th = 1.d0
        endif
        omega_lc = -1.401d0*ts*abs(rhos)**(-1.d0/0.355d0)
        xt = 0.1133d0*abs(rhos)**3.352d0*th**0.732d0*(th+omega_lc*(th+0.535d0))**(-1)     !not tested
    elseif ((abs(ts) < 0.03d0) .and. (rhos == 0.0d0)) then
        xt = 0.0801d0*abs(ts)**(-1.19d0)              !not tested
    else
        dpdd = DPDD_CALC(gl,T,D,nrsubst)
        dpdd = dpdd*1.d6/(Rmix*T)        !1.D0 + 2.D0 * D_ARD + DD_ARDD

        xt = 0.28631d0*del*tau/dpdd
        !xt2 = (Pc(nrsubst)*Factor*Factor/(rhoc(nrsubst)*Rmix*tc(nrsubst)))*del*tau/dpdd       !same as xt, but with calculated constant
    endif

    dpdt = DPDT_CALC(gl,T,D,nrsubst)
    dpdt = dpdt*1.0d6/(Rmix*D)        ! 1+D_ARD-DT_ARDT

    if (xt < 0.d0) xt = 1.d5     ! from refprop Quellcode??! xt < 0 - >  math err. why 1.d5??

    lamcrch4 = 91.855d0/(vis*tau**2)*dpdt**2*xt**0.4681d0*f*1.0d-3       ! exponent 0.4681 = (gamma-ni)/gamma = (1.19-0.633)/1.19

    RETURN
    END


    double precision function dlamctk1(gl,T,D,nrsubst)

    implicit none

    type(type_gl) :: gl


    integer:: nrsubst, i
    double precision:: T, D, tau, del

    tau=abs(gl%tcx%tred_tk(nrsubst))/T
    del=D/gl%tcx%rhored_tk(nrsubst)

    dlamctk1=0.d0

    dlamctk1=gl%tcx%a_tk1(nrsubst,1)*dexp(gl%tcx%a_tk1(nrsubst,2)*(tau+gl%tcx%tsum_tk1(nrsubst,2))**gl%tcx%texp_tk1(nrsubst,2) &
        &  +gl%tcx%a_tk1(nrsubst,3)*(del+gl%tcx%dsum_tk1(nrsubst,3))**gl%tcx%dexp_tk1(nrsubst,3))

    return
    end

    double precision function lambda_TC5(gl, T, D, nrsubst)

    !Programmed according to "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties"; T. H. Chung, M. A. Ajlan, L. L. Lee, and K. E. Starling;
    !Ind. Eng. Chem. Res. 1988, 27, 671-679

    implicit none
    type(type_gl) :: gl
    double precision :: T, Tr, D, vc_TC5, Y_TC5, G1_TC5, H2_TC5, lambdak, lambda0, lambdap, Z_TC5, alpha_TC5, beta_TC5, psi_TC5, etadg, CALtoJ
    double precision, dimension(4,7) :: bji
    double precision, dimension(7) :: B
    integer :: nrsubst, i

    CALtoJ = 4.184D0


    bji = 0.d0
    bji(:,1) = (/ 2.41657d0,     0.74824d0,     -0.91858d0,   121.72100d0/)
    bji(:,2) = (/-0.50924d0,    -1.50936d0,    -49.99120d0,    69.98340d0/)
    bji(:,3) = (/ 6.61069d0,     5.62073d0,     64.75990d0,    27.03890d0/)
    bji(:,4) = (/14.54250d0,    -8.91387d0,     -5.63794d0,    74.34350d0/)
    bji(:,5) = (/ 0.79274d0,     0.82019d0,     -0.69369d0,     6.31734d0/)
    bji(:,6) = (/-5.86340d0,    12.80050d0,      9.58926d0,   -65.52920d0/)
    bji(:,7) = (/81.17100d0,   114.15800d0,    -60.84100d0,   466.77500d0/)

    do i = 1,7
        B(i) = bji(1,i) + bji(2,i)*gl%tcx%accen_TC5(nrsubst) + bji(3,i)*gl%tcx%dipolered_TC5(nrsubst)**4.d0 + bji(4,i)*gl%tcx%kappa_TC5(nrsubst)
    end do

    vc_TC5 = (gl%tcx%lj_sigma(nrsubst)/0.809d0)**3 * gl%factor ![cm^3 / mol]

    Y_TC5 = d/gl%factortrans*vc_TC5/6.d0            !unit conversion of density d from [mol/m³] to [mol/cm³] -> Y => dimensionless?

    G1_TC5 = (1.d0-0.5d0*Y_TC5)/(1.d0-Y_TC5)**3.d0  !dimensionless?

    H2_TC5 = (B(1)*(1.d0 - dexp(-B(4)*Y_TC5))/Y_TC5 + B(2)*G1_TC5 * dexp(B(5)*Y_TC5) + B(3)*G1_TC5)/(B(1)*B(4) + B(2) + B(3)) !dimensionless?

    Tr = T/gl%tcx%tc_TC5(nrsubst)                   !dimensionless

    alpha_TC5 = (CP0_CALC(gl, T, nrsubst) - gl%Req(nrsubst))/(gl%Req(nrsubst)) - 1.5d0 ! dimensionless

    beta_TC5 = 0.7862d0 - 0.7109d0 * gl%tcx%accen_TC5(nrsubst) + 1.3168d0 * gl%tcx%accen_TC5(nrsubst)**2.d0 !dimensionless

    Z_TC5 = 2.d0 + 10.5d0 * (Tr)**2.d0  !dimensionless

    psi_TC5 = 1.d0 + alpha_TC5 * (0.215d0 + 0.28288d0 * alpha_TC5 - 1.061d0 * beta_TC5 + 0.26665d0 * Z_TC5 )/(0.6366d0 + beta_TC5 * Z_TC5 + 1.061d0 * alpha_TC5 * beta_TC5) !dimensionless

    etadg = eta0(gl,T,nrsubst)*gl%factor_VS5eta
    ! eta0 unit: Dilute gas contribution of the viscosity [g/(cm*s)]; unit conversion of molar weight from [kg/mol] to [g/mol] and from calorie to Joule is neccessary
    lambda0 = 7.452d0 * CALtoJ * (etadg/(gl%wm(nrsubst)*gl%factor)) * psi_TC5! [mol / (cm s)]


    lambdak = lambda0 * ((1.d0/H2_TC5) + B(6)*Y_TC5) ! [mol / (cm s)]

    lambdap = CALtoJ *(3.039d-4 *dsqrt(gl%tcx%tc_TC5(nrsubst)/(gl%wm(nrsubst))) / (vc_TC5 **(2.d0/3.d0)))*B(7)*Y_TC5**2.d0*H2_TC5* dsqrt(Tr) * gl%factor !conversion from calorie to Joule is neccessary

    !REFPROP USES DIFFERENT FACTOR
    !lambdap = (3.586d-3 *dsqrt(gl%tcx%tc_TC5(nrsubst)/(gl%wm(nrsubst))) / (vc_TC5 **(2.d0/3.d0)))*B(7)*Y_TC5**2.d0*H2_TC5* dsqrt(Tr) * gl%factor

    lambda_TC5 = (lambdak + lambdap) / gl%factor    !unit: W/m-K

    end function lambda_TC5

    ! DeltaLamda Crit, TK3
    double precision function dlamctk3(gl,T,D,nrsubst)

    implicit none

    type(type_gl) :: gl


    double precision:: omega_lc, Omega0, x, xi, deltachi, chi, chir,qd
    double precision:: y, td, cpd, rhod, visd, z_y, kappainv    !water
    double precision:: T, D, del, boltz_lc, avogad_lc, Rmix, cp, cv, T_445, dpdd, vis, PI_lc, dpddref
    double precision:: delT, delD  ! toluene,sf6
    double precision:: m_ety, gamma_ety, beta_ety, gnu_ety, x0_ety, e1_ety, e2_ety, gam0_ety, x_ety, h_ety,dhdx_ety, B_ety, dpdt, xi0_ety,Kt_ety, F_ety ! holland: ethylene
    integer:: nrsubst     !current substance (number)

    omega_lc = 0.d0
    Omega0 = 0.d0
    x = 0.d0
    xi = 0.d0
    deltachi = 0.d0
    chi = 0.d0
    chir = 0.d0
    cp = 0.d0
    cv = 0.d0
    T_445 = 445.d0
    dpddref = 0.d0
    dpdd = 0.d0

    vis = visdyn_calc(gl,T, D, nrsubst)
    if (vis .lt. 1.d-14) then
        dlamctk3 = -5245.d0
        return
    end if

    !old value PI_lc = 3.14159265358979d0
    PI_lc = 3.14159265358979323846d0

    cp = CP_CALC(gl,T,D,nrsubst)   ! J/(mol K)
    cv = CV_CALC(gl,T,D,nrsubst)    ! J/(mol K)

    del = D/gl%factor                ! mol/m³ - >  mol/L

    call R_mix_calc(gl,Rmix)
    avogad_lc = 6.02214076d23
    boltz_lc = 1.380649d-23

    ! prüfen, ob Temperatur größer ist als Tref und abbruch? - >  Refprop

    if (gl%components(nrsubst) == 'CO2') then    ! vesovic
        if (T <= 445d0) then     !vesovic

            dpdd = DPDD_CALC(gl,T,D,nrsubst)   ! MPa/mol*m^3
            chi = (gl%Pc(nrsubst)/((gl%rhoc(nrsubst)**2)*gl%tc(nrsubst)))*D*T/dpdd

            dpdd = DPDD_CALC(gl,gl%tcx%tref_tk(nrsubst),D,nrsubst)    ! MPa/mol*m^3            chir = (gl%Pc(nrsubst)/((gl%rhoc(nrsubst)**2)*gl%tc(nrsubst)))*D*gl%tref_tk(nrsubst)/dpdd

            deltachi = chi - chir*gl%tcx%tref_tk(nrsubst)/T             ! m^3

            x = gl%tcx%xi0_tk(nrsubst)*((deltachi/gl%tcx%gam0_tk(nrsubst))**(gl%tcx%gnu_tk(nrsubst)/gl%tcx%gamma_tk(nrsubst)))  ! m
        else
            dpdd = DPDD_CALC(gl,T_445,D,nrsubst)   ! Mpa/mol*m^3
            chi = (gl%Pc(nrsubst)/((gl%rhoc(nrsubst)**2)*gl%tc(nrsubst)))*D*T_445/dpdd       !vesovic

            dpdd = DPDD_CALC(gl,gl%tcx%tref_tk(nrsubst),D,nrsubst)    ! Mpa/mol*m^3
            chir = (gl%Pc(nrsubst)/((gl%rhoc(nrsubst)**2)*gl%tc(nrsubst)))*D*gl%tcx%tref_tk(nrsubst)/dpdd    !vesovic

            deltachi = chi - chir*gl%tcx%tref_tk(nrsubst)/T_445

            !------------------------------
            !------------------------------
            !2020/02/19, Monika: copied from Refprop, I did not find it in any paper. It is included here to compare to Refprop, but will not be activated for the Trend release
            if (deltachi.le.1d-2.and.d.gt.gl%rhoc(nrsubst)*1.5d0) deltachi=1.d0-2.d0*2.d0**(deltachi-1.d-2)   !deltachi can go negative far from critical.
            !------------------------------
            !------------------------------


            if (deltachi <= 0.d0) then    ! sqrt(-...) - >  math err (e.g. T = 600K, D = 1000kg/m³)
                dlamctk3 = 0.d0
                RETURN
            endif

            xi = gl%tcx%xi0_tk(nrsubst)*((deltachi/gl%tcx%gam0_tk(nrsubst))**(gl%tcx%gnu_tk(nrsubst)/gl%tcx%gamma_tk(nrsubst)))  ! m

            x = xi*exp(-(T-T_445)/10)   ! for T > 445 (vesovic)   ... T = 600, D = 1000 kg/m³ - >  delchi < 0 - >  NaN
        endif
        !else if ((components(nrsubst) == 'TOLUENE') .and. (T > (tc(nrsubst)+10))) then  ! Assael T > Tc+10...15  ; auch SF6 mit anderen koeff; propane ähnlich mit anderen koeff (nur relativ weit entfernt vom kritischen punkt?)
        !    delT = abs(T/tc(nrsubst)-1)
        !    delD = D/rhoc(nrsubst)-1
        !    C1 = 0.2d-3
        !    C2 = 4.5d-2
        !    C3 = 0.09
        !    dlamctk3 = C1/(C2+delT)*exp(-(C3*delD)**2)
        !    RETURN
    elseif (gl%components(nrsubst) == 'ETHYLENE2') then ! Holland ETHYLENE. not working
        delT = (T - gl%tc(nrsubst))/gl%tc(nrsubst)
        delD = (D - gl%rhoc(nrsubst))/gl%rhoc(nrsubst)
        m_ety = 11.0d0
        gamma_ety = 1.19d0
        beta_ety = 0.355d0
        gnu_ety = (2*beta_ety + gamma_ety)/3
        x0_ety = 0.168d0
        e1_ety = 2.17d0
        e2_ety = 0.287d0
        gam0_ety = x0_ety**gamma_ety/(e1_ety*e2_ety**((gamma_ety - 1.0d0)/(2.0d0*beta_ety))) ! = 0.0532?
        x_ety = abs(delT)/(abs(delD)**(1.0d0/beta_ety))     ! & x=-x0_ety??

        !chi
        if ((abs(delD) <= 0.25d0) .and. (abs(delT) <= 0.025d0)) then    ! not tested
            h_ety = e1_ety*((x_ety + x0_ety)/x0_ety)*(1.0d0 + e2_ety*((x_ety + x0_ety)/x0_ety)**(2.d0*beta_ety))**((gamma_ety - 1.0d0)/(2*beta_ety))
            dhdx_ety = (e1_ety/x0_ety)**((gamma_ety - 1.0d0)/(2.0d0*beta_ety))*((gamma_ety - 1.0d0)/(2.0d0*beta_ety))*x_ety**((gamma_ety - 1.0d0)/(2.0d0*beta_ety) - 1.0d0) + e1_ety**((gamma_ety - 1.0d0)/(2.0d0*beta_ety)) + (1.0d0/(x0_ety**gamma_ety))*gamma_ety*x_ety**(gamma_ety - 1.0d0) + (e1_ety*e2_ety)**((gamma_ety - 1.0d0)/(2.0d0*beta_ety))
            chi = abs(delD)**(gamma_ety/beta_ety)*(h_ety - x_ety*dhdx_ety/beta_ety)
            chi = chi**(-1)
        else
            chi = gam0_ety*abs(delT)**(-gamma_ety)  ! =? D**2*Kt_ety**0.5d0*pc(nrsubst)/(rhoc(nrsubst)**2.0d0)
        endif

        B_ety = x0_ety**(-beta_ety) !or B_ety = abs(delD)/(abs(delT)**(beta_ety))
        xi0_ety = 0.69d0/((B_ety**2.0d0*gl%pc(nrsubst)/gam0_ety*boltz_lc*gl%tc(nrsubst))**0.5d0)
        x = xi0_ety*(chi/gam0_ety)**(gnu_ety/gamma_ety)    !or x = xi0_ety*abs(delT)**(-gnu_ety)
        Kt_ety = chi*gl%rhoc(nrsubst)**2/(D**2*gl%pc(nrsubst))

        dpdt = DPDT_CALC(gl,T,D,nrsubst)
        F_ety = exp(-18.66d0*delT**2)*exp(-4.25d0*delD**4)

        dlamctk3 = (m_ety/(D*avogad_lc*boltz_lc*T))**0.5d0*boltz_lc*T**2*dpdt**2*Kt_ety**0.5d0*F_ety/(6.0d0*PI_lc*vis*x)
        RETURN
    elseif (trim(gl%substcasnr(nrsubst)) == '7732-18-5') then   !water model of Huber (IAPWS release)
        td = T/gl%tc(nrsubst)
        rhod = D/gl%rhoc(nrsubst)
        cpd = (cp/(gl%wm(nrsubst)*gl%factor))/0.46151805d0        ! cp/R, cp [J/(mol K)], R [kJ/(kg K)], wm(nrsubst) [kg/mol]

        visd = vis !vis_reduced visd = vis/(1.d-6) - >  wenn vis(water) berechnet wird

        dpdd = DPDD_CALC(gl,T,D,nrsubst)   ! MPa/mol*m^3
        chi = 1.d0/(D*dpdd)        ! 1/MPa

        dpdd = DPDD_CALC(gl,gl%tcx%tref_tk(nrsubst),D,nrsubst)    ! MPa/mol*m^3
        chir = 1.d0/(D*dpdd)       ! 1/MPa

        deltachi = gl%pc(nrsubst)*rhod**2*(chi - chir*gl%tcx%tref_tk(nrsubst)/T)  ! [-] - REFPROP QT rhod**2?? !! - >  lit: rhod**1

        if (deltachi <= 0.0d0) then
            dlamctk3 = 0.0d0
            RETURN
        endif

        x = gl%tcx%xi0_tk(nrsubst)*((deltachi/gl%tcx%gam0_tk(nrsubst))**(gl%tcx%gnu_tk(nrsubst)/gl%tcx%gamma_tk(nrsubst)))

        y = gl%tcx%qd_inverse_tk(nrsubst)**(-1)*x

        kappainv = cv/cp

        if (y < 1.2d-7) then
            z_y = 0.0d0
        else
            z_y = 2.0d0/(PI_lc*y)*(((1.0d0-kappainv)*atan(y) + kappainv*y) - (1.0d0 - exp(-1.0d0/(y**(-1) + y**2/(3.0d0*rhod**2)))))
        endif

        dlamctk3 = 177.8514d0*rhod*cpd*td*z_y/visd
        dlamctk3 = dlamctk3*gl%tcx%tcx_tk(nrsubst)*1.d-3     ! w m^-1 K^-1
        RETURN
    else    ! not vesovic
        dpdd = DPDD_CALC(gl,T,D,nrsubst)   !MPa/mol*m³
        chi = (gl%Pc(nrsubst)/(gl%rhoc(nrsubst)**2))*D/dpdd

        dpddref = DPDD_CALC(gl,gl%tcx%tref_tk(nrsubst),D,nrsubst)
        chir = (gl%Pc(nrsubst)/(gl%rhoc(nrsubst)**2))*D/dpddref*gl%tcx%tref_tk(nrsubst)/T

        deltachi = chi - chir

        !------------------------------
        !------------------------------
        !2020/02/19, Monika: copied from Refprop, I did not find it in any paper. It is included here to compare to Refprop, but will not be activated for the Trend release
        !Differences between Refprop and Trend of up to 0.1% in terms of thermal conductivity may occur
        if (deltachi.le.1d-2.and.d.gt.gl%rhoc(nrsubst)*1.5d0) deltachi=1.d-2*2.d0**(deltachi-1.d-2)   !deltachi can go negative far from critical.
        !------------------------------
        !------------------------------

        if (deltachi <= 0.d0) then    ! sqrt(-...) - >  math err
            dlamctk3 = 0.d0
            RETURN
        endif

        x = gl%tcx%xi0_tk(nrsubst)*((deltachi/gl%tcx%gam0_tk(nrsubst))**(gl%tcx%gnu_tk(nrsubst)/gl%tcx%gamma_tk(nrsubst)))
        !1.080812064326116D-010
        !7.421939014432470D-011
    endif
    !Implementation according to
    ! R.A. Perkins, J.V. Sengers, I.M. Abdulagatov, M.L. Huber
    ! "Simplified Model for the Critical Thermal-Conductivity Enhancement in Molecular Fluids"
    ! Int. J. Thermophys., 34: 191–212 (2013)
    ! DOI: 10.1007/s10765-013-1409-z

    qd=1.d0/gl%tcx%qd_inverse_tk(nrsubst)
    omega_lc = 2.d0/PI_lc * ((cp-cv)/cp*atan(qd*x) + cv*qd*x/cp)

    Omega0 = 2.d0/PI_lc * (1-dexp((-1.d0)/(gl%tcx%qd_inverse_tk(nrsubst)*(x**(-1.d0)) + ((gl%tcx%qd_inverse_tk(nrsubst)**(-1.d0))*x*gl%rhoc(nrsubst)/D)**2/3.d0)))

    dlamctk3 = (gl%tcx%R0_tk(nrsubst)*boltz_lc*T*(omega_lc - Omega0)*del*cp*1.d9)/(6.d0*PI_lc*vis*x)       ! W m^-1 K^-1

    dlamctk3 = dlamctk3*gl%tcx%tcx_tk(nrsubst)     ! W m^-1 K^-1, before: rhored_tk instead of tcx_tk

    RETURN
    END



    double precision function dlamctk6(gl,T,D,nrsubst)
    !MONIKA, July 2014:THIS FUNCTION HAS NOT BEEN TESTED, YET!!!

    implicit none

    type(type_gl) :: gl


    double precision::omega_lc, Omega0, x, xi, deltachi, chi, chir, cp, cv, T_445, vis, PI_lc, avogad_lc, boltz_lc, dpdd
    double precision:: T, D, del, Rmix
    integer:: nrsubst

    omega_lc = 0.d0
    Omega0 = 0.d0
    x = 0.d0
    xi = 0.d0
    deltachi = 0.d0
    chi = 0.d0
    chir = 0.d0
    cp = 0.d0
    cv = 0.d0
    T_445 = 445.d0

    vis = visdyn_calc(gl,T, D, nrsubst)
    PI_lc = 3.14159265358979d0

    cp = CP_CALC(gl,T,D,nrsubst)   ! J/(mol K)
    cv = CV_CALC(gl,T,D,nrsubst)    ! J/(mol K)

    del = D*1.d-3               ! mol/m³ - >  mol/L

    call R_mix_calc(gl,Rmix)
    avogad_lc = 6.0221367d23
    boltz_lc = Rmix/6.0221367d23

    dpdd = DPDD_CALC(gl,T,D,nrsubst)   !MPa/mol*m³
    chi = (gl%Pc(nrsubst)/(gl%rhoc(nrsubst)**2))*D/dpdd

    dpdd = DPDD_CALC(gl,gl%tcx%tref_tk(nrsubst),D,nrsubst)
    chir = (gl%Pc(nrsubst)/(gl%rhoc(nrsubst)**2))*D/dpdd

    deltachi = chi - chir*gl%tcx%tref_tk(nrsubst)/T

    if (deltachi <= 0.d0) then    ! sqrt(-...) - >  math err
        dlamctk6 = 0.d0
        RETURN
    endif

    x = gl%tcx%xi0_tk(nrsubst)*((deltachi/gl%tcx%gam0_tk(nrsubst))**(gl%tcx%gnu_tk(nrsubst)/gl%tcx%gamma_tk(nrsubst)))

    omega_lc = 2/PI_lc * (((cp-cv)/cp)*atan((gl%tcx%qd_inverse_tk(nrsubst)**(-1))*x) + (cv*(gl%tcx%qd_inverse_tk(nrsubst)**(-1))*x)/cp)

    Omega0 = 2/PI_lc * (1-exp((-1)/(gl%tcx%qd_inverse_tk(nrsubst)*(x**(-1)) + ((gl%tcx%qd_inverse_tk(nrsubst)**(-1))*x*gl%rhoc(nrsubst)/D)**2/3)))

    dlamctk6 = (gl%tcx%R0_tk(nrsubst)*boltz_lc*T*(omega_lc - Omega0)*del*cp*1.d9)/(6*PI_lc*vis*x)       ! W m^-1 K^-1

    dlamctk6 = dlamctk6*gl%tcx%tcx_tk(nrsubst)     ! W m^-1 K^-1, before: rhored_tk instead of tcx_tk

    return
    end


    double precision function tc0R23(gl,T,D,nrsubst)

    implicit none

    type(type_gl) :: gl


    double precision:: T, D, tstar, lamdg, laml, lamcrit, tau, ksi, Rmix
    integer:: nrsubst

    Tstar = T/epskR23

    !eq. (6)
    lamdg=gl%tcx%b1_r23 + gl%tcx%b2_r23*T

    call R_mix_calc(gl,Rmix)
    !eq. (7)
    laml=gl%tcx%c2l_R23*(gl%tcx%rholl_R23*1.d-3)**2/((gl%tcx%rholl_R23-D)*1.d-3)*T**0.5d0*dexp(D/(gl%tcx%rholl_R23-D)*gl%tcx%delgl_R23/Rmix/T)

    ksi=(d-gl%rhoc(nrsubst))*1.d-3
    tau=t-gl%tc(nrsubst)
    !eq.(8a)
    lamcrit=4.d0*gl%tcx%lmax_r23/((dexp(ksi)+dexp(-ksi))*(dexp(tau)+dexp(-tau)))

    !eq. (3)
    tc0R23 = (((gl%tcx%rholl_R23-D)/gl%tcx%rholl_R23)**gl%tcx%c1l_R23*lamdg +(D/gl%tcx%rholl_R23)**gl%tcx%c1l_R23*laml+lamcrit)/gl%factor

    end function



    double precision function tc0helium(gl,T,D,nrsubst)

    implicit none

    type(type_gl) :: gl


    double precision :: T, D, lam0, lamc, lame, tau, del
    double precision :: x0, e1, e2, beta, gamm, delta, eta, dcc, bkt, W, x1, x2b, x2be, dhdx, d2kt, deld, delt, DPDD, DPDT, xx
    double precision, dimension(nrsubst) :: x

    integer :: nrsubst, i, j

    tau = (T/gl%tcx%tred_tc(nrsubst))**(1.d0/3.d0)
    del = D/gl%tcx%rhored_tc(nrsubst)

    !Dilute gas conductivity lambda 0 according eq. 3
    lam0 = gl%tcx%lam0_coeff(nrsubst,1) * (T/gl%tcx%tred_tc(nrsubst))**gl%tcx%lam0_exp(nrsubst,1) * dexp(sum(gl%tcx%lam0_coeff(nrsubst,2:gl%tcx%tcxlam0_num(nrsubst)) * (T/gl%tcx%tred_tc(nrsubst))**gl%tcx%lam0_exp(nrsubst,2:gl%tcx%tcxlam0_num(nrsubst))))

    !Excess gas conductivity lambda e according eq. 17
    lame = 0.d0
    do j=1,gl%tcx%dellam(nrsubst)
        if (ABS(gl%tcx%dellam_ln(nrsubst,j)).lt.1d-12 .or. ABS(del).lt.1d-12) then
            lame = lame + gl%tcx%dellam_coeff(nrsubst,j) * tau ** gl%tcx%dellam_exp(nrsubst,j) * del ** gl%tcx%dellam_exp2(nrsubst,j)
        else
            lame = lame + gl%tcx%dellam_coeff(nrsubst,j) * tau ** gl%tcx%dellam_exp(nrsubst,j) * del ** gl%tcx%dellam_exp2(nrsubst,j) * log(del**gl%tcx%dellam_ln(nrsubst,j))
        end if
    end do

    !Critical region gas conductivity lambda c according eq. 7
    lamc = 0.d0
    if (T.ge.3.5d0 .and. T.le.12.) then
        x0=0.392d0
        e1=2.8461d0
        e2=0.27156d0
        beta=1d0/0.3554d0
        gamm=0.1743d0
        delta=4.304d0
        dcc=69.158d0/gl%wm(nrsubst)
        bkt=0d0

        eta = VISDYN_CALC(gl,T,D,NRSUBST)
        DPDD = DPDD_CALC(gl,T,D,nrsubst)
        DPDT = DPDT_CALC(gl,T,D, nrsubst)
        if (D.gt.0.) bkt=1d0/dPdD/D/10d0**6
        deld = abs(d/dcc-1d0)
        delt = abs(T/5.18992d0-1d0)
        W = (delt/0.2d0)**2 + (deld/0.25d0)**2

        if (W.lt.1d0 .and. D.gt.0d0) then
            xx=delt/deld**beta
            x1=(xx+x0)/x0
            x2b=x1**(2d0/beta)
            x2be=(1d0+e2*x2b)**(gamm/2d0*beta)
            dhdx=e1*x2be/x0+e1*e2/x0*x2b*x2be/(1d0+e2*x2b)*gamm
            d2kt=(delta*e1*x1*x2be-xx*dhdx*beta)*deld**(delta-1d0)
            bkt=W*bkt+(1d0-W)*(dcc/D)**2/d2kt/227460d0
        end if

        if (bkt.ge.0. .and. D.gt.0.) then
            lamc = T**2*sqrt(bkt)/D/gl%wm(nrsubst)/eta*(DPDT*10d2)**2*3.726229668d0*dexp(-18.66d0*delt**2-4.25d0*deld**4)*3.4685233d-5
        end if
        
    end if

    tc0helium = lam0 + lamc + lame

    end function



    double precision function ST_CALC(gl,Temp_st, nrsubst)
    ! M. Thol, 2013

    implicit none

    type(type_gl) :: gl


    double precision:: Temp_st, Teta
    integer :: nrsubst, kst
    double precision:: A_h2o2, N_h2o2

    Teta=0.d0
    ST_CALC=0.d0

    if (gl%stn(nrsubst)%stmodel == 'ST1') then
        !Teta with crictical Temp tred
        Teta = (1-(Temp_st/gl%tred(nrsubst)))
        !sigma0*Factor for ST_CALC in mN/m
        Do kst=1, gl%stn(nrsubst)%term_num_st
            ST_CALC = ST_CALC + (gl%stn(nrsubst)%sigma0_st(kst)*gl%Factor * (Teta**gl%stn(nrsubst)%n_exp_st(kst)))
        enddo
        ST_CALC = ST_CALC/gl%factor
    elseif (gl%substcasnr(1)=="7722-84-1") then
        !thermal conductivity model for hydrogen peroxide according to Yaws_b1999_Chemical Properties Handbook

        if ((trim(gl%inptorig) == "tliq") .or. (trim(gl%inptorig) == "tvap") .or. (trim(gl%inptorig) == "tp") .or. (trim(gl%inptorig) == "td")) then
            if (Temp_st .lt. 272.72d0) then
                ST_CALC = -9912.d0
                return
            elseif (Temp_st .gt. 730.15d0) then
                ST_CALC = -9913.d0
                return
            else
                A_h2o2 = 141.031d0
                N_h2o2 = 1.2222d0
                ST_CALC = A_h2o2 * (1.d0 - Temp_st / gl%tc(nrsubst)) ** N_h2o2
                return
            end if
        else
            ST_CALC = -9955.d0
            return
        end if
    else
        ST_CALC = -5234.d0
        return
    endif
    end function



    double precision function DE_CALC(gl,Temp_de, Dens_de, nrsubst)
    ! M. Thol, 2013





    implicit none

    type(type_gl) :: gl


    double precision, dimension(30) :: polar_de
    double precision:: Temp_de, Dens_de, B_de, A_de, C_de, polar_de_ges, eps_lc
    double precision:: Dens_de_new, dielec_one, dielec_two, g_de, Sum_one, Dens_de_new2
    integer :: jde, kde,  nrsubst,last_position, i

    dielec_one=0.d0
    dielec_two=0.d0
    polar_de=0.d0
    polar_de_ges=0.d0
    g_de=0.d0
    B_de=0.d0
    A_de=0.d0
    C_de=0.d0
    Sum_one=0.d0
    Dens_de_new2=0.d0
    Dens_de_new=0.d0
    eps_lc=0.d0


    !checking Mixture Models
    If (nrsubst == 0) then
        Do i=1, gl%ncomp
            if ((gl%de%demodel(i) == 'DE3') .or. (gl%de%demodel(i) == 'DE4')) then
                continue
            else
                DE_CALC = -5235 !Not possible to mix with DE2
                return
            endif
        enddo
    endif
    !
    !First Possible DE Model DE3 (similar to DE4)
    if ((gl%de%demodel(1) == 'DE3') .or. (gl%de%demodel(1) == 'DE4')  ) then

        Dens_de_new = Dens_de/gl%factortrans   !mol/m³ into mol/cm³

        !Dens_mix calculation
        Do i=1,gl%ncomp !number of fluids
            Sum_one = Sum_one + (gl%molfractions(i)/(gl%rhoc(i)))
        enddo
        Dens_de_new2 = Dens_de_new*Sum_one

        DO i=1, gl%ncomp

            !Amue coefficient
            if (gl%de%term_num1_de(i) /= 0) then
                polar_de(i) = gl%de%coeffde(i,1)*gl%de%ref_Temp_de(i)* &
                    ((Dens_de_new2*gl%rhoc(i))**gl%de%dexpde(i,1))*(Temp_de**gl%de%texpde(i,1))
            endif

            !Aepsilon coefficient
            Do kde=(1+gl%de%term_num1_de(i)),(gl%de%term_num1_de(i)+gl%de%term_num2_de(i))
                polar_de(i) = polar_de(i) + gl%de%coeffde(i,kde)*(Dens_de_new2*gl%rhoc(i))**gl%de%dexpde(i,kde)* &
                    ((Temp_de/gl%de%ref_Temp_de(i))-1.d0)**int(gl%de%texpde(i,kde))
            enddo

            !Bepsilon and C coefficient
            Do jde=(1+gl%de%term_num1_de(i)+gl%de%term_num2_de(i)), &
                (gl%de%term_num1_de(i)+gl%de%term_num2_de(i)+gl%de%term_num3_de(i))

                polar_de(i) = polar_de(i) + gl%de%coeffde(i,jde)*(Dens_de_new2*gl%rhoc(i))** &
                    gl%de%dexpde(i,jde)*((gl%de%ref_Temp_de(i)/Temp_de)-1.d0)**int(gl%de%texpde(i,jde))
            enddo

        enddo

        !changing from P to dielectric constant
        if (gl%ncomp == 1) then
            if (gl%de%demodel(1) == 'DE3')  then
                !PCM expression
                DE_CALC = (-1.d0-(2.d0*polar_de(1)))/(polar_de(1)-1.d0)
            else
                !PK eypression
                dielec_one = ((1.d0+(9.d0*polar_de(1)))*0.25d0) + ( ((1.d0+(9.d0*polar_de(1)))*0.25d0)**2 + 0.5d0)**0.5d0
                dielec_two = ((1.d0+(9.d0*polar_de(1)))*0.25d0) - ( ((1.d0+(9.d0*polar_de(1)))*0.25d0)**2 + 0.5d0)**0.5d0
                if (dielec_one >= dielec_two) then
                    DE_CALC = dielec_one
                else
                    DE_CALC = dielec_two
                endif

            endif
            !mixture calculation
        else
            Do i=1, gl%ncomp
                !Mixture Calculation for DE3 Models. First epsilon, then electric polarisation with PK
                if (gl%de%demodel(i) == 'DE3')  then
                    eps_lc =  (-1.d0-(2.d0*polar_de(i)))/(polar_de(i)-1.d0)
                    polar_de(i) = ((eps_lc-1.d0)*(2.d0*eps_lc + 1.d0))/(9.d0*eps_lc)

                endif
            enddo

            !electric Polarisation for whole Mixture
            Do i=1, gl%ncomp
                polar_de_ges = polar_de_ges + (((gl%molfractions(i)/gl%rhoc(i))/Sum_one)*polar_de(i))
            enddo


            !changing from P to dielectric constant with PCM
            dielec_one = ((1.d0+(9.d0*polar_de_ges))*0.25d0) + ( ((1.d0+(9.d0*polar_de_ges))*0.25d0)**2 + 0.5d0)**0.5d0
            dielec_two = ((1.d0+(9.d0*polar_de_ges))*0.25d0) - ( ((1.d0+(9.d0*polar_de_ges))*0.25d0)**2 + 0.5d0)**0.5d0
            if (dielec_one >= dielec_two) then
                DE_CALC = dielec_one
            else
                DE_CALC = dielec_two
            endif
        endif

        !Next possible DE Model DE2
    elseif (gl%de%demodel(1) == 'DE2') then

        Dens_de_new = Dens_de/gl%Factor   !mol/m³ into mol/dm³


        !Sum for all coefficients except the last one (special calculation)
        Do kde=1, gl%de%term_num1_de(1)
            g_de = g_de + gl%de%coeffde(1,kde)*((Dens_de_new/gl%de%ref_dens_de(1))**gl%de%dexpde(1,kde)) &
                *((gl%de%ref_temp_de(1)/Temp_de)**gl%de%texpde(1,kde))
        enddo


        !the last coefficient
        last_position = gl%de%term_num1_de(1)+gl%de%term_num2_de(1)
        !g coefficient
        g_de = g_de + 1.d0 + gl%de%coeffde(1,last_position) &
            *(Dens_de_new/gl%de%ref_dens_de(1))*((Temp_de/gl%de%texpde(1,last_position))-1.d0) &
            **(-gl%de%pexpde(1,last_position))

        !A coefficient
        A_de = ((avogad*((gl%dipole(1)*Debye_change)**2)*Dens_de_new*g_de)/(permit*boltz*Temp_de))*1000.d0

        !B coefficient
        B_de = ((avogad*molec_polar_water*Dens_de_new)/(3.d0*permit))*1000.d0

        !root in numerator with A and B
        C_de = (9.d0 + 2.d0*A_de + 18.d0*B_de + A_de**2 + 10.d0*A_de*B_de + 9.d0*(B_de**2))**0.5d0

        !dielectric constant
        DE_CALC = (1.d0 + A_de + 5.d0*B_de + C_de)/(4.d0 - 4.d0*B_de)

    else ! added by Theresa
        DE_CALC = -5234.d0
    endif

    end function


    end module transport_module
