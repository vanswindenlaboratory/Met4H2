
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
    !Cite as: Span, R.; Beckmuller, R.; Eckermann, T.; Herrig, S.; Hielscher, S.;
    !          Jager, A.; Mickoleit, E.; Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (2020):
    !          TREND. Thermodynamic Reference and Engineering Data 5.0.
    !          Lehrstuhl fur Thermodynamik, Ruhr-Universitat Bochum.

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



    module module_trend_info

    !Andreas November 2013
    !Release information for TREND.
    character(255), parameter:: version_parameter = 'TREND 5.0 Prof. Dr.-Ing. Roland Span, Thermodynamics, Ruhr-Universit酹 Bochum, 2020.'

    end module



    module module_parameters

    !implicit none
    !save

    integer,parameter::bad_nrmixcomp = 1

    !dimension of arrays for tau, delta- (or T, rho-) derivatives (e.g. GETDER, FNRDER)
    integer, parameter :: nderivs = 15
    integer, parameter :: nderivsi = 10 !dimension parameter for ideal part
    integer, parameter :: nderivsdep = 10 !dimension parameter for departure function
    integer, parameter :: nderivshs = 10 !dimension parameter for hard sphere terms

    !dimension of arrays for the first composition derivatives (e.g. in PC SAFT subroutines)
    integer, parameter :: nx1derivs = 10

    !dimension of arrays for the second composition derivatives (e.g. in PC SAFT subroutines)
    integer, parameter :: nx2derivs = 6

    !dimension of arrays for the third composition derivatives (e.g. in PC SAFT subroutines)
    integer, parameter :: nx3derivs = 3

    !dimension of character array of flash properties and units
    integer, parameter :: flash_chardim = 37


    !Boltzmann's constant
    double precision, parameter :: k_B = 1.d0/6.022045d+23*8.31441d0
    !Avogadro constant
    double precision, parameter :: N_A = 6.022045d+23
    !for Phase-Equilibrium calculations
    integer, parameter:: imax = 800

    !Avogadro-const
    double precision, parameter :: Avogad = 6.0221367d0*10.d0**23
    !permittivity of free space
    double precision, parameter :: permit = 8.854187818d0*10.d0**(-12)
    !Boltzman-const
    double precision, parameter :: Boltz = 1.380658d0*10.d0**(-23)
    !molecular polarity of water
    double precision, parameter :: molec_polar_water = 1.636d0*10.d0**(-40)
    !Debye conversion
    double precision, parameter :: Debye_change = 3.308894879d0*10.d0**(-30)

    !parameters of ethylene: P. M. Holland, B. E. Eaton, H. J. M. Hanley, J. Phys. Chem. Ref. Data, 12(4):917, 1983
    double precision,parameter :: rhoc_ETHY=0.221D0   !g/cm?
    double precision,parameter :: M_ETHY=28.054       !g/mol

    !parameters of methanol: H. W. Xiang, A. Laesecke, M. L. Huber, J. Phys. Chem. Ref. Data, 35(4): 1597-1, 2006
    double precision,parameter :: epsk=577.87d0        !LJ parameter: epsilon/k [K] for Stockmayer potential
    double precision,parameter :: sig0=0.3408d-9       !collision diameter [m] for Stockmayer potential
    double precision,parameter :: delst=0.4575d0       !delta for Stockmayer potential
    double precision,parameter :: kbol=1.3806505d-23    !Boltzmann constant [J/K]
    double precision,parameter :: avo=6.0221415d23     !Avogadro constant [1/mol]
    double precision,parameter :: mw_oh1=32.04216d0    !Molar mass [g/mol] used for the vsicosity model
    double precision,parameter :: pi=3.14159265d0
    double precision,parameter :: Tcrit=512.6d0        !critical temperature [K]
    double precision,parameter :: Dcrit=273.d0         !critical density [kg/m設

    !R23
    double precision, parameter:: sigR23=0.4278d0
    double precision, parameter:: epskR23=243.91d0

    !Supporting variables for PR
    double precision, parameter:: sqrt_2 = 2.D0**0.5
    double precision, parameter:: consA_PR = 2.D0**0.5 + 1.D0   !consA_PR = 2^0.5 + 1
    double precision, parameter:: consB_PR = 2.D0**0.5 - 1.D0   !consB_PR = 2^0.5 - 1

    !Andreas J輍er, January 2017
    !PSRK variables
    double precision, parameter :: A1_PSRK = -0.64663D0     !The parameter A1 for the SRK. A1 = -ln((u+1)/u) with u = 1.1

    double precision, parameter :: R_HelmgE = 8.3144598D0    !Universal gas constant

    !costald
    !parameters for tait equation: saturated liquid density
    double precision,parameter:: ataitold=-1.52816d0
    double precision,parameter:: btaitold=1.43907d0
    double precision,parameter:: ctaitold=-0.81446d0
    double precision,parameter:: dtaitold=0.190454d0
    double precision,parameter:: etaitold=-0.296123d0
    double precision,parameter:: ftaitold=0.386914d0
    double precision,parameter:: gtaitold=-0.0427258d0
    double precision,parameter:: htaitold=-0.0480645d0

    !parameters for tait equation: homogeneous density
    double precision,parameter:: atait=-9.070217d0
    double precision,parameter:: btait=62.45326d0
    double precision,parameter:: dtait=-135.1102d0
    double precision,parameter:: ftait=4.79594d0
    double precision,parameter:: gtait=0.250047d0
    double precision,parameter:: htait=1.14188d0
    double precision,parameter:: jtait=0.0861488d0
    double precision,parameter:: ktait=0.0344483d0

    double precision, parameter :: piPCSAFT = 3.14159265359d0

    !PC SAFT
    !universal model constants for equation 18 and 19, Table 1, Gross, Sadowski 2001
    !a_0i, a_1i, a_2i, b_0i, b_1i, b_2i with i element of {0, 6}
    double precision, dimension (0:6), parameter :: a0PCSAFT = (/ 0.9105631445d0, 0.6361281449d0, 2.6861347891d0, -26.547362491d0, 97.759208784d0, -159.59154087d0, 91.297774084d0 /)
    double precision, dimension (0:6), parameter :: a1PCSAFT = (/ -0.3084016918d0, 0.1860531159d0, -2.5030047259d0, 21.419793629d0, -65.255885330d0, 83.318680481d0, -33.746922930d0 /)
    double precision, dimension (0:6), parameter :: a2PCSAFT = (/ -0.0906148351d0, 0.4527842806d0, 0.5962700728d0, -1.7241829131d0, -4.1302112531d0, 13.776631870d0, -8.6728470368d0 /)
    double precision, dimension (0:6), parameter :: b0PCSAFT = (/ 0.7240946941d0, 2.2382791861d0, -4.0025849485d0, -21.003576815d0, 26.855641363d0, 206.55133841d0, -355.60235612d0 /)
    double precision, dimension (0:6), parameter :: b1PCSAFT = (/ -0.5755498075d0, 0.6995095521d0, 3.8925673390d0, -17.215471648d0, 192.67226447d0, -161.82646165d0, -165.20769346d0 /)
    double precision, dimension (0:6), parameter :: b2PCSAFT = (/ 0.0976883116d0, -0.2557574982d0, -9.1558561530d0, 20.642075974d0, -38.804430052d0, 93.626774077d0, -29.666905585d0 /)

    !universal model constants for equation 15, 16 and 17, Table 2, Gross 2005
    !a_0i, a_1i, a_2i, b_0i, b_1i, b_2i, c_0i, c_1i, c_2i with i element of {0, 4}
    double precision, dimension (0:4), parameter :: a0PCSAFTQ = (/ 1.2378308d0, 2.4355031d0, 1.6330905d0, -1.6118152d0, 6.9771185d0 /)
    double precision, dimension (0:4), parameter :: a1PCSAFTQ = (/ 1.2854109d0, -11.465615d0, 22.086893d0, 7.4691383d0, -17.197772d0 /)
    double precision, dimension (0:4), parameter :: a2PCSAFTQ = (/ 1.7942954d0, 0.7695103d0, 7.2647923d0, 94.486699d0, -77.148458d0 /)
    double precision, dimension (0:4), parameter :: b0PCSAFTQ = (/ 0.4542718d0, -4.5016264d0, 3.5858868d0, 0.d0, 0.d0 /)
    double precision, dimension (0:4), parameter :: b1PCSAFTQ = (/ -0.8137340d0, 10.064030d0, -10.876631d0, 0.d0, 0.d0 /)
    double precision, dimension (0:4), parameter :: b2PCSAFTQ = (/ 6.8682675d0, -5.1732238d0, -17.240207d0, 0.d0, 0.d0 /)
    double precision, dimension (0:4), parameter :: c0PCSAFTQ = (/ -0.5000437d0, 6.5318692d0, -16.014780d0, 14.425970d0, 0.d0 /)
    double precision, dimension (0:4), parameter :: c1PCSAFTQ = (/ 2.0002094d0, -6.7838658d0, 20.383246d0, -10.895984d0, 0.d0 /)
    double precision, dimension (0:4), parameter :: c2PCSAFTQ = (/ 3.1358271d0, 7.2475888d0, 3.0759478d0, 0.d0, 0.d0 /)


    !universal model constants for equation 12, 13 and 14, Table 1, Gross 2006
    !a_0i, a_1i, a_2i, b_0i, b_1i, b_2i, c_0i, c_1i, c_2i with i element of {0, 4}
    double precision, dimension (0:4), parameter :: a0PCSAFTD = (/ 0.3043504d0, -0.1358588d0, 1.4493329d0, 0.3556977d0, -2.0653308d0/)
    double precision, dimension (0:4), parameter :: a1PCSAFTD = (/ 0.9534641d0, -1.8396383d0, 2.0131180d0, -7.3724958d0, 8.2374135d0 /)
    double precision, dimension (0:4), parameter :: a2PCSAFTD = (/ -1.1610080d0, 4.5258607d0, 0.9751222d0, -12.281038d0, 5.9397575d0 /)
    double precision, dimension (0:4), parameter :: b0PCSAFTD = (/ 0.2187939d0, -1.1896431d0, 1.1626889d0, 0.d0, 0.d0 /)
    double precision, dimension (0:4), parameter :: b1PCSAFTD = (/ -0.58731640d0, 1.2489132d0, -0.5085280d0, 0.d0, 0.d0 /)
    double precision, dimension (0:4), parameter :: b2PCSAFTD = (/ 3.4869576d0, -14.915974d0, 15.372022d0, 0.d0, 0.d0 /)
    double precision, dimension (0:4), parameter :: c0PCSAFTD = (/ -0.0646774d0, 0.1975882d0, -0.8087562d0, 0.6902849d0, 0.d0 /)
    double precision, dimension (0:4), parameter :: c1PCSAFTD = (/ -0.9520876d0, 2.9924258d0, -2.3802636d0, -0.2701261d0, 0.d0 /)
    double precision, dimension (0:4), parameter :: c2PCSAFTD = (/ -0.6260979d0, 1.2924686d0, 1.6542783d0,-3.4396744d0, 0.d0 /)

    !fit_helmge
    double precision, parameter:: resize_fac_abc = 1000.D0 !Factor to bring a,b,c for fitting to approximately the same order of magnitude
    integer, parameter :: nr_max_params = 100   !Specifies the maximum number of adjustable parameters



    integer, parameter :: coeff_structure_all_zeros = 0
    integer, parameter :: coeff_structure_int = 1
    integer, parameter :: coeff_structure_real = 2

    ! AGA8
    ! Constants
    double precision, parameter :: Rag8 = 8.31451d0    ! [J.mol^-1.K^-1] Universal gas constant

    ! Coefficients and Parameters
    ! a_n, b_n, c_n, k_n, u_n, g_n, q_n, f_n, s_n, w_n with n element of {1:58}
    double precision, dimension(1:58), parameter :: a_n =  (/    0.1538326d0,    1.341953d0,        -2.998583d0, -0.04831228d0,    0.3757965d0,   -1.589575d0,  -0.05358847d0,   0.88659463d0, -0.71023704d0,   -1.471722d0, &
        1.32185035d0, -0.78665925d0, 2.29129d-9,   0.1576724d0,   -0.4363864d0, -0.04408159d0, -0.003433888d0,   0.03205905d0,  0.02487355d0,  0.07332279d0, &
        -0.001600573d0,   0.6424706d0,       -0.4162601d0, -0.06689957d0,    0.2791795d0, - 0.6966051d0, -0.002860589d0, -0.008098836d0,    3.150547d0, 0.007224479d0, &
        -0.7057529d0,   0.5349792d0,      -0.07931491d0,   -1.418465d0,   -5.99905d-17,   0.1058402d0,   0.03431729d0, -0.007022847d0,  0.02495587d0,  0.04296818d0, &
        0.7465453d0,  -0.2919613d0,         7.294616d0,   -9.936757d0, -0.005399808d0,  -0.2432567d0,   0.04987016d0,  0.003733797d0,    1.874951d0, 0.002168144d0, &
        -0.6587164d0, 0.000205518d0,      0.009776195d0, -0.02048708d0,   0.01557322d0, 0.006862415d0, -0.001226752d0,  0.002850908d0 /)
    double precision, dimension(1:58), parameter :: b_n =  (/  1.d0,   1.d0,  1.d0,  1.d0,   1.d0,  1.d0,  1.d0,   1.d0,   1.d0,  1.d0, &
        1.d0,   1.d0,  1.d0,  1.d0,   1.d0,  1.d0,  1.d0,   1.d0,   2.d0,  2.d0, &
        2.d0,   2.d0,  2.d0,  2.d0,   2.d0,  2.d0,  2.d0,   3.d0,   3.d0,  3.d0, &
        3.d0,   3.d0,  3.d0,  3.d0,   3.d0,  3.d0,  3.d0,   4.d0,   4.d0,  4.d0, &
        4.d0,   4.d0,  4.d0,  4.d0,   5.d0,  5.d0,  5.d0,   5.d0,   5.d0,  6.d0, &
        6.d0,   7.d0,  7.d0,  8.d0,   8.d0,  8.d0,  9.d0,   9.d0 /)
    double precision, dimension(1:58), parameter :: c_n =  (/  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  1.d0,  1.d0,   1.d0,  1.d0,  1.d0,   1.d0,   0.d0,  0.d0, &
        1.d0,   1.d0,  1.d0,  1.d0,   1.d0,  1.d0,  1.d0,   0.d0,   1.d0,  1.d0, &
        1.d0,   1.d0,  1.d0,  1.d0,   1.d0,  1.d0,  1.d0,   0.d0,   0.d0,  1.d0, &
        1.d0,   1.d0,  1.d0,  1.d0,   0.d0,  1.d0,  1.d0,   1.d0,   1.d0,  0.d0, &
        1.d0,   0.d0,  1.d0,  1.d0,   1.d0,  1.d0,  1.d0,   1.d0 /)
    double precision, dimension(1:58), parameter :: k_n =  (/  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  3.d0,  2.d0,   2.d0,  2.d0,  4.d0,   4.d0,   0.d0,  0.d0, &
        2.d0,   2.d0,  2.d0,  4.d0,   4.d0,  4.d0,  4.d0,   0.d0,   1.d0,  1.d0, &
        2.d0,   2.d0,  3.d0,  3.d0,   4.d0,  4.d0,  4.d0,   0.d0,   0.d0,  2.d0, &
        2.d0,   2.d0,  4.d0,  4.d0,   0.d0,  2.d0,  2.d0,   4.d0,   4.d0,  0.d0, &
        2.d0,   0.d0,  2.d0,  1.d0,   2.d0,  2.d0,  2.d0,   2.d0 /)
    double precision, dimension(1:58), parameter :: u_n =  (/  0.d0,  0.5d0,  1.d0, 3.5d0, -0.5d0, 4.5d0, 0.5d0,  7.5d0,  9.5d0,  6.d0, &
        12.d0, 12.5d0, -6.d0,  2.d0,   3.d0,  2.d0,  2.d0,  11.d0, -0.5d0, 0.5d0, &
        0.d0,   4.d0,  6.d0, 21.d0,  23.d0, 22.d0, -1.d0, -0.5d0,   7.d0, -1.d0, &
        6.d0,   4.d0,  1.d0,  9.d0, -13.d0, 21.d0,  8.d0, -0.5d0,   0.d0,  2.d0, &
        7.d0,   9.d0, 22.d0, 23.d0,   1.d0,  9.d0,  3.d0,   8.d0,  23.d0, 1.5d0, &
        5.d0, -0.5d0,  4.d0,  7.d0,   3.d0,  0.d0,  1.d0,   0.d0 /)
    double precision, dimension(1:58), parameter :: g_n =  (/  0.d0,   0.d0,  0.d0,  0.d0,   1.d0,  1.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   1.d0,  0.d0,  0.d0,   0.d0,   1.d0,  0.d0, &
        0.d0,   1.d0,  1.d0,  1.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        1.d0,   0.d0,  0.d0,  1.d0,   0.d0,  1.d0,  0.d0,   0.d0 /)
    double precision, dimension(1:58), parameter :: q_n =  (/  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  1.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  1.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  1.d0,  0.d0,   1.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  1.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   1.d0,  0.d0,  0.d0,   0.d0,  0.d0,  1.d0,   0.d0,   1.d0,  0.d0, &
        0.d0,   1.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   1.d0 /)
    double precision, dimension(1:58), parameter :: f_n =  (/  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  1.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  1.d0,   0.d0,   0.d0,  1.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   1.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0 /)
    double precision, dimension(1:58), parameter :: s_n =  (/  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   1.d0,   1.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0 /)
    double precision, dimension(1:58), parameter :: w_n =  (/  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  1.d0, &
        1.d0,   1.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0,   0.d0,  0.d0, &
        0.d0,   0.d0,  0.d0,  0.d0,   0.d0,  0.d0,  0.d0,   0.d0 /)


    ! Characterization Parameters
    ! M_i, E_i, K_i, G_i, Q_i, F_i, S_i, W_i with i components of {1:21} (Element Compunds)
    ! Element Compounds i -      1: Methane
    !                            2: Nitrogen
    !                            3: Carbon dioxide
    !                            4: Ethane
    !                            5: Propane
    !                            6: i-Butane
    !                            7: n-Butane
    !                            8: i-Pentane
    !                            9: n-Pentane
    !                           10: n-Hexane
    !                           11: n-Heptane
    !                           12: n-Octane
    !                           13: n-Nonane
    !                           14: n-Decane
    !                           15: Hydrogen
    !                           16: Oxygen
    !                           17: Carbon monoxide
    !                           18: Water
    !                           19: Hydrogen sulfide
    !                           20: Helium
    !                           21: Argon

    ! Molar mass M_i (g/mol)
    ! Writer's note: These are not the current internationally accepted values of the molar masses,
    ! but are identical to those reported in the 1994 edition of AGA 8 to maintain consistency.
    ! For other properties not related to the DETAIL or GROSS equations of state, the molar masses in AGA 5 or GPA 2145 should be used.
    double precision, dimension(1:21), parameter :: M_i_temp =  (/  16.043d0,   28.0135d0,   44.01d0,   30.07d0, 44.097d0,  58.123d0, 58.123d0,   72.15d0,  72.15d0, 86.177d0, &
        100.204d0,   114.231d0, 128.258d0, 142.285d0, 2.0159d0, 31.9988d0,  28.01d0, 18.0153d0, 34.082d0, 4.0026d0, &
        39.948d0 /)
    ! Energy parameter E_i (K)
    double precision, dimension(1:21), parameter :: E_i_temp =  (/  151.3183d0,   99.73778d0,   241.9606d0,   244.1667d0, 298.1183d0, 324.0689d0, 337.6389d0, 365.5999d0, 370.6823d0, 402.636293d0, &
        427.72263d0, 450.325022d0, 470.840891d0, 489.558373d0, 26.95794d0, 122.7667d0, 105.5348d0, 514.0156d0,  296.355d0,   2.610111d0, &
        119.6299d0 /)
    ! Size parameter K_i (m^3/kmol)^(1/3)
    double precision, dimension(1:21), parameter :: K_i_temp =  (/ 0.4619255d0, 0.4479153d0, 0.4557489d0, 0.5279209d0,  0.583749d0, 0.6406937d0, 0.6341423d0, 0.6738577d0, 0.6798307d0, 0.7175118d0, &
        0.7525189d0, 0.784955d0,  0.8152731d0, 0.8437826d0, 0.3514916d0, 0.4186954d0, 0.4533894d0, 0.3825868d0, 0.4618263d0, 0.3589888d0, &
        0.4216551d0 /)
    ! Orientation parameter G_i
    double precision, dimension(1:21), parameter :: G_i_temp =  (/  0.d0, 0.027815d0, 0.189065d0,   0.0793d0, 0.141239d0, 0.256692d0, 0.281835d0, 0.332267d0, 0.366911d0, 0.289731d0, &
        0.337542d0, 0.383381d0, 0.427354d0, 0.469659d0, 0.034369d0,    0.021d0, 0.038953d0,   0.3325d0,   0.0885d0,       0.d0, &
        0.d0 /)
    ! Quadruple parameter Q_i
    double precision, dimension(1:21), parameter :: Q_i_temp =  (/ 0.d0, 0.d0, 0.69d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
        0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.06775d0, 0.633276d0, 0.d0, &
        0.d0 /)
    ! High temperature parameter F_i
    double precision, dimension(1:21), parameter :: F_i_temp =  (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
        0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
        0.d0 /)
    ! Dipole parameter S_i
    double precision, dimension(1:21), parameter :: S_i_temp =  (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,     0.d0,   0.d0, 0.d0, &
        0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.5822d0, 0.39d0, 0.d0, &
        0.d0 /)
    ! Association parameter W_i
    double precision, dimension(1:21), parameter :: W_i_temp =  (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
        0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
        0.d0 /)


    ! Binary Interaction Parameter Values
    ! E_ij, U_ij, K_ij, G_ij
    ! Writer's note: Values for the j,i pair are equal to the i,j pair, that is, Eji =Eij, Uji =Uij, Kji =Kij, and Gji =Gij.
    ! Values of 1.0 should be used for all binary interaction parameters except for the entries in this table.

    !                                                                        1: Methane  2: Nitrogen  3: Carbon dioxide    4: Ethane   5: Propane  6: i-Butane  7: n-Butane  8: i-Pentane  9: n-Pentane  10: n-Hexane  11: n-Heptane  12: n-Octane  13: n-Nonane  14: n-Decane  15: Hydrogen  16: Oxygen  17: Carbon monoxide    18: Water  19: Hydrogen sulfide  20: Helium  21: Argon
    double precision, dimension(21,21), parameter :: E_ij_temp = reshape( (/       1.d0,   0.97164d0,        0.960644d0,        1.d0,  0.994635d0,   1.01953d0,  0.989844d0,    1.00235d0,   0.999268d0,   1.107274d0,     0.88088d0,   0.880973d0,   0.881067d0,   0.881161d0,    1.17052d0,       1.d0,          0.990126d0,  0.708218d0,           0.931484d0,       1.d0,      1.d0, &         !  1: Methane
        0.97164d0,        1.d0,         1.02274d0,   0.97012d0,  0.945939d0,  0.946914d0,  0.973384d0,    0.95934d0,    0.94552d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,    1.08632d0,    1.021d0,           1.00571d0,  0.746954d0,           0.902271d0,       1.d0,      1.d0, &         !  2: Nitrogen
        0.960644d0,   1.02274d0,              1.d0,  0.925053d0,  0.960237d0,  0.906849d0,  0.897362d0,   0.726255d0,   0.859764d0,   0.855134d0,    0.831229d0,    0.80831d0,   0.786323d0,   0.765171d0,    1.28179d0,       1.d0,               1.5d0,  0.849408d0,           0.955052d0,       1.d0,      1.d0, &         !  3: Carbon dioxide
        1.d0,   0.97012d0,        0.925053d0,        1.d0,   1.02256d0,        1.d0,   1.01306d0,         1.d0,    1.00532d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,    1.16446d0,       1.d0,                1.d0,  0.693168d0,           0.946871d0,       1.d0,      1.d0, &         !  4: Ethane
        0.994635d0,  0.945939d0,        0.960237d0,   1.02256d0,        1.d0,        1.d0,    1.0049d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,   1.034787d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  5: Propane
        1.01953d0,  0.946914d0,        0.906849d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,        1.3d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  6: i-Butane
        0.989844d0,  0.973384d0,        0.897362d0,   1.01306d0,    1.0049d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,        1.3d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  7: n-Butane
        1.00235d0,   0.95934d0,        0.726255d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  8: i-Pentane
        0.999268d0,   0.94552d0,        0.859764d0,   1.00532d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  9: n-Pentane
        1.107274d0,        1.d0,        0.855134d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           1.008692d0,       1.d0,      1.d0, &         ! 10: n-Hexane
        0.88088d0,        1.d0,        0.831229d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           1.010126d0,       1.d0,      1.d0, &         ! 11: n-Heptane
        0.880973d0,        1.d0,         0.80831d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           1.011501d0,       1.d0,      1.d0, &         ! 12: n-Octane
        0.881067d0,        1.d0,        0.786323d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           1.012821d0,       1.d0,      1.d0, &         ! 13: n-Nonane
        0.881161d0,        1.d0,        0.765171d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           1.014089d0,       1.d0,      1.d0, &         ! 14: n-Decane
        1.17052d0,   1.08632d0,         1.28179d0,   1.16446d0,  1.034787d0,       1.3d0,       1.3d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,               1.1d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 15: Hydrogen
        1.d0,     1.021d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 16: Oxygen
        0.990126d0,   1.00571d0,             1.5d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,        1.1d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 17: Carbon monoxide
        0.708218d0,  0.746954d0,        0.849408d0,  0.693168d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 18: Water
        0.931484d0,  0.902271d0,        0.955052d0,  0.946871d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,   1.008692d0,    1.010126d0,   1.011501d0,   1.012821d0,   1.014089d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 19: Hydrogen sulfide
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 20: Helium
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0  /), &     ! 21: Argon
        (/ 21, 21 /) )

    !                                                                        1: Methane  2: Nitrogen  3: Carbon dioxide    4: Ethane   5: Propane  6: i-Butane  7: n-Butane  8: i-Pentane  9: n-Pentane  10: n-Hexane  11: n-Heptane  12: n-Octane  13: n-Nonane  14: n-Decane  15: Hydrogen  16: Oxygen  17: Carbon monoxide    18: Water  19: Hydrogen sulfide  20: Helium  21: Argon
    double precision, dimension(21,21), parameter :: U_ij_temp = reshape( (/       1.d0,  0.886106d0,        0.963827d0,        1.d0,  0.990877d0,        1.d0,  0.992291d0,         1.d0,    1.00367d0,   1.302576d0,    1.191904d0,   1.205769d0,   1.219634d0,   1.233498d0,    1.15639d0,       1.d0,                1.d0,        1.d0,           0.736833d0,       1.d0,      1.d0, &         !  1: Methane
        0.886106d0,        1.d0,        0.835058d0,  0.816431d0,  0.915502d0,        1.d0,  0.993556d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,   0.408838d0,       1.d0,                1.d0,        1.d0,           0.993476d0,       1.d0,      1.d0, &         !  2: Nitrogen
        0.963827d0,  0.835058d0,              1.d0,   0.96987d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,   1.066638d0,    1.077634d0,   1.088178d0,   1.098291d0,   1.108021d0,         1.d0,       1.d0,               0.9d0,        1.d0,            1.04529d0,       1.d0,      1.d0, &         !  3: Carbon dioxide
        1.d0,  0.816431d0,         0.96987d0,        1.d0,  1.065173d0,      1.25d0,      1.25d0,       1.25d0,       1.25d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,    1.61666d0,       1.d0,                1.d0,        1.d0,           0.971926d0,       1.d0,      1.d0, &         !  4: Ethane
        0.990877d0,  0.915502d0,              1.d0,  1.065173d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  5: Propane
        1.d0,        1.d0,              1.d0,      1.25d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  6: i-Butane
        0.992291d0,  0.993556d0,              1.d0,      1.25d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  7: n-Butane
        1.d0,        1.d0,              1.d0,      1.25d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  8: i-Pentane
        1.00367d0,        1.d0,              1.d0,      1.25d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  9: n-Pentane
        1.302576d0,        1.d0,        1.066638d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           1.028973d0,       1.d0,      1.d0, &         ! 10: n-Hexane
        1.191904d0,        1.d0,        1.077634d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           1.033754d0,       1.d0,      1.d0, &         ! 11: n-Heptane
        1.205769d0,        1.d0,        1.088178d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           1.038338d0,       1.d0,      1.d0, &         ! 12: n-Octane
        1.219634d0,        1.d0,        1.098291d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           1.042735d0,       1.d0,      1.d0, &         ! 13: n-Nonane
        1.233498d0,        1.d0,        1.108021d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           1.046966d0,       1.d0,      1.d0, &         ! 14: n-Decane
        1.15639d0,  0.408838d0,              1.d0,   1.61666d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 15: Hydrogen
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 16: Oxygen
        1.d0,        1.d0,             0.9d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 17: Carbon monoxide
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 18: Water
        0.736833d0,  0.993476d0,         1.04529d0,  0.971926d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,   1.028973d0,    1.033754d0,   1.038338d0,   1.042735d0,   1.046966d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 19: Hydrogen sulfide
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 20: Helium
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0  /), &     ! 21: Argon
        (/ 21, 21 /) )

    !                                                                        1: Methane  2: Nitrogen  3: Carbon dioxide    4: Ethane   5: Propane  6: i-Butane  7: n-Butane  8: i-Pentane  9: n-Pentane  10: n-Hexane  11: n-Heptane  12: n-Octane  13: n-Nonane  14: n-Decane  15: Hydrogen  16: Oxygen  17: Carbon monoxide    18: Water  19: Hydrogen sulfide  20: Helium  21: Argon
    double precision, dimension(21,21), parameter :: K_ij_temp = reshape( (/       1.d0,   1.00363d0,        0.995933d0,        1.d0,  1.007619d0,        1.d0,  0.997596d0,         1.d0,   1.002529d0,   0.982962d0,    0.983565d0,   0.982707d0,   0.981849d0,   0.980991d0,    1.02326d0,       1.d0,                1.d0,        1.d0,            1.00008d0,       1.d0,      1.d0, &         !  1: Methane
        1.00363d0,        1.d0,        0.982361d0,   1.00796d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,    1.03227d0,       1.d0,                1.d0,        1.d0,           0.942596d0,       1.d0,      1.d0, &         !  2: Nitrogen
        0.995933d0,  0.982361d0,              1.d0,   1.00851d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,   0.910183d0,    0.895362d0,   0.881152d0,    0.86752d0,   0.854406d0,         1.d0,       1.d0,                1.d0,        1.d0,            1.00779d0,       1.d0,      1.d0, &         !  3: Carbon dioxide
        1.d0,   1.00796d0,         1.00851d0,        1.d0,  0.986893d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,    1.02034d0,       1.d0,                1.d0,        1.d0,           0.999969d0,       1.d0,      1.d0, &         !  4: Ethane
        1.007619d0,        1.d0,              1.d0,  0.986893d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  5: Propane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  6: i-Butane
        0.997596d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  7: n-Butane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  8: i-Pentane
        1.002529d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  9: n-Pentane
        0.982962d0,        1.d0,        0.910183d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           0.96813d0,       1.d0,      1.d0, &         ! 10: n-Hexane
        0.983565d0,        1.d0,        0.895362d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           0.96287d0,       1.d0,      1.d0, &         ! 11: n-Heptane
        0.982707d0,        1.d0,        0.881152d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           0.957828d0,       1.d0,      1.d0, &         ! 12: n-Octane
        0.981849d0,        1.d0,         0.86752d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           0.952441d0,       1.d0,      1.d0, &         ! 13: n-Nonane
        0.980991d0,        1.d0,        0.854406d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,           0.948338d0,       1.d0,      1.d0, &         ! 14: n-Decane
        1.02326d0,   1.03227d0,              1.d0,   1.02034d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 15: Hydrogen
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 16: Oxygen
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 17: Carbon monoxide
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 18: Water
        1.00008d0,  0.942596d0,         1.00779d0,  0.999969d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,    0.96813d0,     0.96287d0,   0.957828d0,   0.952441d0,   0.948338d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 19: Hydrogen sulfide
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 20: Helium
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0  /), &     ! 21: Argon
        (/ 21, 21 /) )

    !                                                                        1: Methane  2: Nitrogen  3: Carbon dioxide    4: Ethane   5: Propane  6: i-Butane  7: n-Butane  8: i-Pentane  9: n-Pentane  10: n-Hexane  11: n-Heptane  12: n-Octane  13: n-Nonane  14: n-Decane  15: Hydrogen  16: Oxygen  17: Carbon monoxide    18: Water  19: Hydrogen sulfide  20: Helium  21: Argon
    double precision, dimension(21,21), parameter :: G_ij_temp = reshape( (/       1.d0,        1.d0,        0.807653d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,    1.95731d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  1: Methane
        1.d0,        1.d0,        0.982746d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  2: Nitrogen
        0.807653d0,  0.982746d0,              1.d0,  0.370296d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,   1.67309d0,                 1.d0,       1.d0,      1.d0, &         !  3: Carbon dioxide
        1.d0,        1.d0,        0.370296d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  4: Ethane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  5: Propane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  6: i-Butane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  7: n-Butane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  8: i-Pentane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         !  9: n-Pentane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 10: n-Hexane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 11: n-Heptane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 12: n-Octane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 13: n-Nonane
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 14: n-Decane
        1.95731d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 15: Hydrogen
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 16: Oxygen
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 17: Carbon monoxide
        1.d0,        1.d0,         1.67309d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 18: Water
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 19: Hydrogen sulfide
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0, &         ! 20: Helium
        1.d0,        1.d0,              1.d0,        1.d0,        1.d0,        1.d0,        1.d0,         1.d0,         1.d0,         1.d0,          1.d0,         1.d0,         1.d0,         1.d0,         1.d0,       1.d0,                1.d0,        1.d0,                 1.d0,       1.d0,      1.d0  /), &     ! 21: Argon
        (/ 21, 21 /) )


    ! Coefficients and Paramesters of the Ideal-gas Heat Capacity Equations
    ! n0_ik with i element of {1:21} (Element Compounds), k element of {1:7}
    double precision, dimension(7,21), parameter :: n0_ik_temp = reshape( (/ 29.83843397d0,  -15999.69151d0,     4.00088d0,  0.76315d0,   0.00460d0,     8.74432d0,  -4.46921d0, &
        17.56770785d0,  -2801.729072d0,     3.50031d0,  0.13732d0,   -0.1466d0,     0.90066d0,        0.d0, &
        20.65844696d0,  -4902.171516d0,     3.50002d0,  2.04452d0,  -1.06044d0,     2.03366d0,   0.01393d0, &
        36.73005938d0,  -23639.65301d0,     4.00263d0,  4.33939d0,   1.23722d0,     13.1974d0,  -6.01989d0, &
        44.70909619d0,  -31236.63551d0,     4.02939d0,  6.60569d0,   3.19700d0,     19.1921d0,  -8.37267d0, &
        34.30180349d0,  -38525.50276d0,     4.06714d0,  8.97575d0,   5.25156d0,     25.1423d0,   16.1388d0, &
        36.53237783d0,  -38957.80933d0,     4.33944d0,  9.44893d0,   6.89406d0,     24.4618d0,   14.7824d0, &
        43.17218626d0,  -51198.30946d0,          4.d0,  11.7618d0,   20.1101d0,     33.1688d0,        0.d0, &
        42.67837089d0,  -45215.83000d0,          4.d0,  8.95043d0,   21.8360d0,     33.4032d0,        0.d0, &
        46.99717188d0,  -52746.83318d0,          4.d0,  11.6977d0,   26.8142d0,     38.6164d0,        0.d0, &
        52.07631631d0,  -57104.81056d0,          4.d0,  13.7266d0,   30.4707d0,     43.5561d0,        0.d0, &
        57.25830934d0,  -60546.76385d0,          4.d0,  15.6865d0,   33.8029d0,     48.1731d0,        0.d0, &
        62.09646901d0,  -66600.12837d0,          4.d0,  18.0241d0,   38.1235d0,     53.3415d0,        0.d0, &
        65.93909154d0,  -74131.45483d0,          4.d0,  21.0069d0,   43.4931d0,     58.3657d0,        0.d0, &
        13.07520288d0,  -5836.943696d0,     2.47906d0,  0.95806d0,   0.45444d0,     1.56039d0,  -1.37560d0, &
        16.80171730d0,  -2318.322690d0,     3.50146d0,  1.07558d0,   1.01334d0,          0.d0,        0.d0, &
        17.45786899d0,  -2635.244116d0,     3.50055d0,  1.02865d0,   0.00493d0,          0.d0,        0.d0, &
        21.57882705d0,  -7766.733078d0,     4.00392d0,  0.01059d0,   0.98763d0,     3.06904d0,        0.d0, &
        21.58309440d0,  -6069.035869d0,          4.d0,  3.11942d0,   1.00243d0,          0.d0,        0.d0, &
        10.04639507d0,      -745.375d0,         2.5d0,       0.d0,        0.d0,          0.d0,        0.d0, &
        10.04639507d0,      -745.375d0,         2.5d0,       0.d0,        0.d0,          0.d0,        0.d0 /), &
        (/ 7, 21 /) )

    ! theta0_ik with i element of {1:21} (Element Compounds); k element of {4:7}
    double precision, dimension(4:7,21), parameter :: theta0_ik_temp = reshape( (/ 820.659d0,  178.410d0,  1062.82d0,  1090.53d0, &
        662.738d0,  680.562d0,  1740.06d0,       0.d0, &
        919.306d0,  865.070d0,  483.553d0,  341.109d0, &
        559.314d0,  223.284d0,  1031.38d0,  1071.29d0, &
        479.856d0,  200.893d0,  955.312d0,  1027.29d0, &
        438.270d0,  198.018d0,  1905.02d0,  893.765d0, &
        468.270d0,  183.636d0,  1914.10d0,  903.185d0, &
        292.503d0,  910.237d0,  1919.37d0,       0.d0, &
        178.670d0,  840.538d0,  1774.25d0,       0.d0, &
        182.326d0,  859.207d0,  1826.59d0,       0.d0, &
        169.789d0,  836.195d0,  1760.46d0,       0.d0, &
        158.922d0,  815.064d0,  1693.07d0,       0.d0, &
        156.854d0,  814.882d0,  1693.79d0,       0.d0, &
        164.947d0,  836.264d0,  1750.24d0,       0.d0, &
        228.734d0,  326.843d0,  1651.71d0,  1671.69d0, &
        2235.71d0,  1116.69d0,       0.d0,       0.d0, &
        1550.45d0,  704.525d0,       0.d0,       0.d0, &
        268.795d0,  1141.41d0,  2507.37d0,       0.d0, &
        1833.63d0,  847.181d0,       0.d0,       0.d0, &
        0.d0,       0.d0,       0.d0,       0.d0, &
        0.d0,       0.d0,       0.d0,       0.d0 /), &
        (/ 4, 21 /) )

    end module module_parameters

    module module_hdrt_parameters
    integer, parameter :: np_4ph = 200
    ! Fundamental constants:
    !data Rgas/8.314462d0/    ! [J.mol^-1.K^-1] universal gas constant (RUB value 13.12.2011)
    !data kB/1.38064880d-23/   ! [J/K] Boltzmann constant
    !data Nav/6.02214129d23/   ! [mol^-1] Avogadro constant (Rgas/kB)
    ! CODATA, 4.3.2017, (http://physics.nist.gov/cuu/Constants/index.html)
    double precision, parameter :: Rgas = 8.3144598d0       ! [J.mol^-1.K^-1] universal gas constant
    double precision, parameter :: Nav = 6.022140857d23     ! [mol^-1] Avogadro constant (Rgas/kB)
    double precision, parameter :: kB = 1.38064852d-23      ! [J/K] Boltzmann constant
    double precision, parameter :: kcal = 4186.8d0          ! [J/kcal] - conversion

    end module module_hdrt_parameters








    module module_all_types

    use, intrinsic :: iso_c_binding
    use, intrinsic :: ieee_arithmetic
    use module_parameters
    use module_hdrt_parameters
    !use functionals

    type type_model_info
        integer,dimension(13):: eq_type_list=(/1,2,3,4,6,7,8,9,15,51,52,53,81/)
    end type
    
    type type_state_point

        logical :: fluid_present
        character(12) :: input
        double precision, dimension(2) :: prop
        double precision, dimension(:), allocatable :: molfractions
        !handle for gl copy and phase property type
        type(c_ptr) :: gl_check_handle
        type(c_ptr) :: ph_prop_handle
        integer :: nr_checks !dimension of ph_prop_handle
    end type

    type type_uncertainty

        !implicit none
        !save

        !in some cases of uncertainty estimation it is not wanted to have a twophase check | default=.true.
        LOGICAL :: blnCheck2Phase
        !path of uncertainty file
        !CHARACTER(255) :: UNCTY_FILE = ''
        !variables contain known uncertainty

        double precision :: uncty

        DOUBLE PRECISION :: dblUncty_CP
        DOUBLE PRECISION :: dblUncty_D
        DOUBLE PRECISION :: dblUncty_WS
        DOUBLE PRECISION :: dblUncty_P
        DOUBLE PRECISION :: dblUncty_CV
        DOUBLE PRECISION :: dblUncty_CP_vle
        DOUBLE PRECISION :: dblUncty_D_vle
        DOUBLE PRECISION ::	dblUncty_WS_vle
        DOUBLE PRECISION :: dblUncty_P_vle
        DOUBLE PRECISION :: dblUncty_CV_vle
        !variable will contain biggest estimated uncertainty at the end of calculation
        DOUBLE PRECISION :: estUncty_out

        !Andreas Aprip 2013 OUTPUT VARIABLES FOR UNCERTAINTIES
        !-------------------------------
        Double precision:: Uncty_cp
        Double precision:: Uncty_d
        Double precision:: Uncty_w
        Double precision:: Uncty_p
        Double precision:: Uncty_cv
        Double precision:: Uncty_h
        Double precision:: Uncty_s
        Double precision:: Uncty_u
        !-------------------------------
        ! New check variables if same temperature and density was called last time
        double precision:: temp_last, dens_last
        !-------------------------------

        double precision :: rho_limit_uncty
        double precision :: cp_limit_uncty
        double precision :: cv_limit_uncty
        double precision :: ws_limit_uncty

        double precision, dimension(100) :: ni_varied_pos
        double precision, dimension(100) :: ni_varied_neg
        double precision, dimension(20) :: ni_varied_pos_ideal
        double precision, dimension(20) :: ni_varied_neg_ideal
        double precision, dimension(100) :: used_coeff_variance_pos
        double precision, dimension(100) :: used_coeff_variance_neg
        integer, dimension(100) :: used_base_pos
        integer, dimension(100) :: used_base_neg

    end type type_uncertainty

    !Erik, April 2018
    type type_COSMO_SAC
        !--------------------------------------------------------------------------------------
        ! Variables for the COSMO-SAC model
        !--------------------------------------------------------------------------------------

        !integer :: interval !, version
        double precision, dimension(30) :: Acosmo, Vcosmo
        double precision, dimension(51) :: counter
        double precision, dimension(51,30) :: sigma
        !Variables COSMO-SAC Version 2 and 3
        !double precision, dimension(30,3) :: Acosmo_ver23
        double precision, dimension(51,30,3) :: sigma_v23
        double precision, dimension(51,3) :: counter_v23
        integer :: COSMO_ver                !Choose COSMO-SAC version with this variable:
        !1 COSMO-SAC Version by Lin, Sandler 2002;
        !2 COSMO-SAC improvements by Hsieh, Sandler, Lin 2010;
        !3 COSMO-SAC considering dispersive interactions by Hsieh, Lin, Vrabec 2014
        !4 COSMO-SAC Xiong et al 2014
        !november 2018, Erik
        logical :: solv_method
        integer, dimension(3) :: interval_min, interval_max
        double precision, dimension(51,51,3,3) :: sum_square_counter, hb_inter

        double precision :: temp_cosmo_prev
        double precision, dimension(30) :: molefraction_cosmo_prev
        double precision, dimension(51,30,3) :: seggamma_pure_prev
        double precision, dimension(51,30) :: seggamma_pure_prev_v1
        double precision, dimension(51,3) :: seggamma_prev
        double precision, dimension(51) :: seggamma_prev_v1
        double precision, dimension(15,30) :: ln_gamma_C_cosmo_prev, ln_gamma_R_cosmo_prev
        double precision, dimension(15) :: gE_C_cosmo_prev, gE_R_cosmo_prev
        integer, dimension(30) :: molecule_type    !for dispersive interactions
        double precision, dimension(30) :: eps_molecule    !for dispersive interactions
        double precision, dimension(30,17) :: m_tau
        double precision, dimension(30) :: NBT_cosmo    !Normal boiling temperature to calculate molar Volume in Version 4 (Xiong)
        double precision, dimension(30) :: ln_gamma_dsp
        !integer :: int

        !for analytical derivations of ln_gamma_res
        logical :: analytical
        double precision, dimension(51) :: seggamma_gl
        double precision, dimension(51,30) :: seggamma_pure_gl
        double precision, dimension(51,3) :: seggamma_v23_gl
        double precision, dimension(51,30,3) :: seggamma_pure_v23_gl
        double precision, dimension(51,51) :: dFi_mix_v1_gl
        double precision, dimension(51,51,3,3) :: dFi_mix_gl
        double precision, dimension(51,51,30) :: dFi_pure_v1_gl
        double precision, dimension(51,51,30,3,3) :: dFi_pure_gl
        double precision, dimension(51,51) :: d2Fi_dseggamma_dT_mix
        double precision, dimension(51) :: dseggamma_dT_mix
        !numerical test of Fi0_mix_v1_gl
        double precision, dimension(51) :: Fi0_mix_v1_gl
        integer, dimension(3) :: start_sigma, end_sigma

        double precision, dimension(51) :: sigma_profile_mix_gl        !Sigma-profile of the mixture
        double precision, dimension(51,3) :: sigma_profile_mix_v23_gl
        double precision, dimension(51,51) :: delta_w_gl                         !Electrostatic interactions
        double precision, dimension(51,51,3,3) :: delta_w_v23_gl
        double precision :: aeff_gl

        !for second derivative of segment activity coefficient in respect to xa and xb
        double precision, dimension(51,30) :: dseggamma_dxa_mix
        double precision, dimension (51,30) :: dFi_dxa

        !numerical test
        double precision, dimension(51) :: dDetBdt_mix_numeric
        double precision, dimension(51,30) :: d2seggamma_dxa_dxb_numeric, dDetAdxb_mix_numeric, d2seggamma_dxa_dT_mix_test
        double precision, dimension(15,30,30) :: ln_gamma_R_Trho_dxa_numeric, ln_gamma_R_Trho_dxa_test
        double precision, dimension(60,60,30) :: dMatrixB_dxa
        double precision, dimension(60,60,51,30) :: dMatrixA_dxa_num
        double precision, dimension(51,51,30) :: d2F_dseggammadxa_num
        double precision, dimension(51,30) :: d2F_dxadxb_num, d2lnseggammadxadT_mix_numeric, detMatrixA_mix_num
        double precision, dimension(30) :: gE_R_dxadxb_num
        double precision, dimension(15,30) :: gE_R_Trho_dxa_num, gE_R_Trho_dxa_test
        double precision, dimension(60,60,51,30) :: MatrixA
        double precision, dimension(15) :: gE_R_Trho_num
        double precision, dimension(15,30) :: ln_gamma_R_Trho_test

    end type type_COSMO_SAC


    type type_viscosity

        logical:: eta_read      !true == viscosity model has been read
        logical:: h2o_read
        logical:: coll_read

        character*255, dimension(30):: etamodel
        character*255, dimension(30)::pointereta
        character*255, dimension(30)::pointerceaux

        double precision, dimension(30):: tmineta
        double precision, dimension(30):: tmaxeta
        double precision, dimension(30):: pmaxeta
        double precision, dimension(30):: rhomaxeta

        integer,          dimension(30):: dilgas_term


        !double precision, dimension(30):: hij_terms_eta
        character*255, dimension(30):: pointercrit_eta, pointer_lambda, pointer_collmod
        double precision, dimension(30):: low_temp_eta                      !lower temperature limit of the viscosity model [K]
        double precision, dimension(30):: upp_temp_eta                      !upper temperature limit of the viscosity model [K]
        double precision, dimension(30):: max_dens_eta                      !upper density limit of the viscosity model [mol/l]
        double precision, dimension(30):: red_temp_eta                      !reducing temperature of the viscosity model [K]
        double precision, dimension(30):: red_dens_eta                      !reducing density of the viscosity model [mol/l]
        double precision, dimension(30):: red_vis_eta                       !reducing viscosity of the viscosity model [Pa*s]
        double precision, dimension(30):: Chapman_Enskog1, Chapman_Enskog2
        integer, dimension(30):: nterm1_eta, nterm2_eta, nterm3_eta, nterm4_eta
        integer, dimension(30):: nterm5_eta, nterm6_eta, nterm7_eta, nterm8_eta
        integer, dimension(30):: nterm9_eta, nterm_eta
        double precision, dimension(30,100):: coeff_dil
        double precision, dimension(30,100):: exp_dil
        double precision, dimension(30,100):: coeff_hieta
        double precision, dimension(30,100,100):: hieta
        integer,          dimension(30,100):: iexpeta
        integer,          dimension(30,100):: jexpeta
        double precision, dimension(30,100):: exp1eta
        double precision, dimension(30,100):: exp2eta
        double precision, dimension(30,100):: coeff_ai                      !viscosity in the zero density limit
        double precision, dimension(30,100):: exp_ai
        double precision, dimension(30,100):: coeff_bi                      !second viscosity virial coefficient
        double precision, dimension(30,100):: exp_bi
        double precision, dimension(30,100):: coeff_ci                      !third viscosity virial coefficient
        double precision, dimension(30,100):: exp_ci
        double precision, dimension(30,100):: coeff_di                      !reduced viscosity at high density according to the Enskog theory
        double precision, dimension(30,100):: exp_di
        double precision, dimension(30,100):: coeff_ei                      !reduced viscosity at high density according to the Enskog theory
        double precision, dimension(30,100):: exp_ei
        double precision, dimension(30,100):: bieta                         !coefficient
        double precision, dimension(30,100):: tieta                         !temperature exponent
        double precision, dimension(30,100):: cieta                         !coefficient
        double precision, dimension(30,100):: ti_eta                        !temperature exponent
        double precision, dimension(30,100):: alpha_eta                     !coefficient
        double precision, dimension(30,100):: exp_tau_eta                   !temperature exponent
        double precision, dimension(30,100):: exp_del_eta                   !density exponent
        integer, dimension(30):: ndilgas_term, ndens_initial
        double precision, dimension(30,100):: coeff_in,texp_in       !coefficients and temperature exponents for the initial density contribution
        double precision, dimension(30):: tinit_red, etainit_red !reducing parameters for the initial density contribution
        character*255, dimension(30):: pointer_eta
        double precision, dimension(30):: lejo_sigma                        !Lennard-Jones coefficient sigma [nm]
        double precision, dimension(30):: lejo_epka                         !Lennard-Jones coefficient epsilon/kappa [K]
        double precision, dimension(30):: c1_eta
        double precision, dimension(30):: tred_eta              !reducing parameter
        double precision, dimension(30):: etared                !reducing parameter
        double precision, dimension(30,100):: Ninuer            !coefficient of the residual part
        double precision, dimension(30,100):: tinuer            !temperature exponent of the residual part
        double precision, dimension(30,100):: dinuer            !density exponent of the residual part
        double precision, dimension(30,100):: di0nuer           !density exponent of del0 of the residual part
        double precision, dimension(30,100):: linuer            !density exponent in the exponential part of the residual part

        double precision:: rhol_R23, c1_R23, c2_R23, delg_R23, etamax_R23 !parameters for the VS0 model of R23

        ! parameters for the critical enhancement of water
        double precision, dimension(30):: lowbo_temp_mue2                   !lower boundary [K] of the region around the crit. temp.
        double precision, dimension(30):: uppbo_temp_mue2                   !upper boundary [K] of the region around the crit. temp.
        double precision, dimension(30):: lowbo_dens_mue2                   !lower boundary	[kg/m^3] of the region around the crit. dens.
        double precision, dimension(30):: uppbo_dens_mue2                   !upper boundary	[kg/m^3] of the region around the crit. Dens.
        double precision, dimension(30):: rhoc_eta
        double precision, dimension(30):: Tc_eta
        double precision, dimension(30):: corr_ref_mue2                     !reference correlation length
        double precision, dimension(30):: qc_mue2
        double precision, dimension(30):: qd_mue2
        double precision, dimension(30):: nue_mue2
        double precision, dimension(30):: gamma_mue2
        double precision, dimension(30):: Xi0_mue2
        double precision, dimension(30):: gamma0_mue2
        double precision, dimension(30):: T_mue2
        double precision, dimension(30):: X_mue2
        !----------------------------------------------------------------------------------------------------------------------------------------
        !!parameters of methanol: H. W. Xiang, A. Laesecke, M. L. Huber, J. Phys. Chem. Ref. Data, 35(4): 1597-1, 2006
        !double precision,parameter :: epsk=577.87d0        !LJ parameter: epsilon/k [K] for Stockmayer potential
        !double precision,parameter :: sig0=0.3408d-9       !collision diameter [m] for Stockmayer potential
        !double precision,parameter :: delst=0.4575d0       !delta for Stockmayer potential
        !double precision,parameter :: kbol=1.3806505d-23    !Boltzmann constant [J/K]
        !double precision,parameter :: avo=6.0221415d23     !Avogadro constant [1/mol]
        !double precision,parameter :: mw_oh1=32.04216d0    !Molar mass [g/mol] used for the vsicosity model
        !double precision,parameter :: pi=3.14159265d0
        !double precision,parameter :: Tcrit=512.6d0        !critical temperature [K]
        !double precision,parameter :: Dcrit=273.d0         !critical density [kg/m設

        double precision, dimension(13) :: a_OH1    !parameters for collision integrals: eq. (4), eq. (8), eq. (9)
        double precision, dimension(9) :: b_OH1     !parameters for second viscosity virial coefficient: eq. (13)
        double precision, dimension(2) :: c_OH1     !Parameters for the third viscosity virial coefficient: eq. (14)
        double precision, dimension(7) :: d_OH1     !Parameters for the temperature and density dependency of sigHS: eq. (17)
        double precision, dimension(9) :: e_OH1     !Parameters for the temperature and density dependency of sigHS: eq. (17)
        !----------------------------------------------------------------------------------------------------------------------------------------
        !!parameters of ethylene: P. M. Holland, B. E. Eaton, H. J. M. Hanley, J. Phys. Chem. Ref. Data, 12(4):917, 1983
        !double precision,parameter :: rhoc_ETHY=0.221D0   !g/cm?
        !double precision,parameter :: M_ETHY=28.054       !g/mol
        double precision, dimension(9) :: GV             !coefficients and exponents for the diluted gas, eq. (7)
        double precision :: A_ETHY, B_ETHY, C_ETHY, F_ETHY           !parameters for 1st density correction for the moderately dense gas: eq. (3)
        double precision :: E_ETHY                                   !Parameter for remainder, eq. (4)
        double precision, dimension(7) :: J_ETHY                     !Parameters for remainder, eq. (4)
        !----------------------------------------------------------------------------------------------------------------------------------------
        !parameters collision model
        integer, dimension(30):: nterm_coll
        double precision, dimension(30,100):: ticoll            !coeff, power of Tstar

        !
        double precision, dimension(30):: tredeta
        double precision, dimension(30):: rhoredeta
        double precision, dimension(30):: visredeta

        integer, dimension(30):: term_num1_eta, term_num2_eta, term_num3_eta, term_num4_eta
        integer, dimension(30):: term_num5_eta, term_num6_eta, term_num7_eta, term_num8_eta


        !VS1-Model
        double precision, dimension(30):: low_tempcoll
        double precision, dimension(30):: upp_tempcoll
        double precision, dimension(30):: low_prescoll
        double precision, dimension(30):: upp_denscoll
        double precision, dimension(30):: eklj
        double precision, dimension(30):: slj
        double precision, dimension(30):: ce
        double precision, dimension(30):: cep
        double precision, dimension(30):: a0_vs1
        double precision, dimension(30):: a1_vs1
        double precision, dimension(30):: a2_vs1
        double precision, dimension(30):: tredetadg
        double precision, dimension(30):: visredetadg

        double precision, dimension(30):: etaB2
        double precision, dimension(30):: treddens
        integer, dimension(30):: colltermnr
        integer, dimension(30):: denstermnr


        double precision, dimension(30,100):: coeffcoll
        double precision, dimension(30,100):: powercoll

        double precision, dimension(30,100):: coeffbetastar
        double precision, dimension(30,100):: powerbetastar

        double precision, dimension(30,100):: a_vs1
        double precision, dimension(30,100):: powertau

        double precision, dimension(30,100):: b_vs1
        double precision, dimension(30,100):: ptausp
        double precision, dimension(30,100):: pdelsp
        double precision, dimension(30,100):: pdel0sp
        double precision, dimension(30,100):: expdel

        !VS2
        !B. A. Younglove, J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982
        double precision, dimension(30):: const19_VS2
        double precision, dimension(30):: exp19_VS2
        double precision, dimension(30,4):: Fv_VS2
        double precision, dimension(30,8):: Ev_VS2

        !VS4
        double precision, dimension(30):: a_0,a_1,a_2,aa0,aa1,aa2,b0_vs4,b1_vs4,b2_vs4,bb0,bb1,bb2,c0,c_1,c_2,cc0,cc1,cc2,dd_0,dd_1,dd_2,e_0,e_1,e_2
        double precision, dimension(30,6):: d0_vs4, d0exp

        !VS5
        !T-H. Chung, M. Ajlan, L.L. Lee and K.E. Starling, Ind. Eng. Chem. Res. 1998, 27, 671-679.
        double precision, dimension(30):: slj_VS5
        double precision, dimension(30):: tc_VS5
        double precision, dimension(30):: accen_VS5
        double precision, dimension(30):: dipolered_VS5
        double precision, dimension(30):: kappa_VS5
        double precision, dimension(30):: addchung_VS5

        !VS6





        !VS9
        !E. Vogel, R. Span, and S. Herrmann, Journal of Physical and Chemical Reference Data 44, 043101 (2015)
        double precision, dimension (30,20)::coeff
        double precision, dimension (30,20)::texpo
        double precision, dimension (30,20)::dexpo
        double precision, dimension (30,20)::beta_VS9
        double precision, dimension (30,20)::eta_VS9
        double precision, dimension (30,20)::pi_VS9
        double precision, dimension (30,20)::li_VS9
        double precision, dimension (30,20)::par1
        double precision, dimension (30,20)::par2

        character*10, dimension (30):: omega_model

        !
        !additional coefficient for the dilute-gas contribution in VS9 according to Herrmann et al. (2018)
        double precision, dimension (30,20):: coeff_hi

        integer, dimension(30):: term_expo_eta
        integer, dimension(30):: term_com_expo_eta

    end type type_viscosity

    type surface_tension_t

        character*255::  stmodel
        double precision:: low_temp_st = 0d0
        double precision:: upp_temp_st
        double precision:: upp_press_st
        double precision:: max_dens_st
        integer      :: term_num_st
        double precision:: Temp_crit_st
        double precision, dimension(100) :: sigma0_st,n_exp_st

    end type surface_tension_t

    type type_thermal_conductivity

        !implicit none
        !save

        logical:: tcx_read      !true == thermal conductivity model has been read

        character*10, dimension(30):: tcmodel, tkmodel, pointer_hardtc

        integer, dimension(30):: num_dil_tc, den_dil_tc, num_bgrd_tc, den_bgrd_tc

        double precision, dimension(30):: tcx_dil, tcx_bgrd
        double precision, dimension(30):: tmin_tc
        double precision, dimension(30):: tmax_tc
        double precision, dimension(30):: pmax_tc
        double precision, dimension(30):: rhomax_tc
        double precision, dimension(30):: tred_dil, tred_bgrd, rhored_bgrd
        double precision, dimension(30):: lj_sigma, lj_epskap, const_eq20, Texp_eq20, F_tc3, rm ! TC3
        double precision, dimension(30,30):: eta0, Fvi, Evi ! TC3
        integer, dimension(30):: tcxlam0_num,tcxlam0_den, dellam, dellamcr
        double precision, dimension(30):: tred_tc, rhored_tc, etared_tc      !TC0 D2O
        double precision, dimension(30,30):: lam0_coeff, dellam_coeff, dellamcr_coeff, lam0_exp, dellam_exp, dellam_exp2, dellamcr_exp, dellam_ln   !TC0 D2O
        double precision, dimension(30,30):: num_tc_coeff_dil, den_tc_coeff_dil, den_tc_exp_dil, num_tc_powerT_dil
        double precision, dimension(30,30):: num_tc_coeff_bgrd, num_tc_rho_bgrd, num_tc_powerT_bgrd, rho_exp_bgrd

        double precision, dimension(30):: tmin_tk
        double precision, dimension(30):: tmax_tk
        double precision, dimension(30):: pmax_tk
        double precision, dimension(30):: rhomax_tk
        integer, dimension(30):: term_tk
        double precision, dimension(30):: tred_tk, rhored_tk, tcx_tk
        double precision, dimension(30):: gnu_tk, gamma_tk, R0_tk, z_visonly_tk, c_tk, xi0_tk, gam0_tk, qd_inverse_tk, tref_tk
        logical:: tcx_crit_read

        !TC5
        !T-H. Chung, M. Ajlan, L.L. Lee and K.E. Starling, Ind. Eng. Chem. Res. 1998, 27, 671-679.
        !double precision, dimension(30):: slj_TC5
        double precision, dimension(:), allocatable:: tc_TC5
        double precision, dimension(:), allocatable:: accen_TC5
        double precision, dimension(:), allocatable:: dipolered_TC5
        double precision, dimension(:), allocatable:: kappa_TC5
        double precision, dimension(:), allocatable:: addchung_TC5

        !R23
        double precision:: rholl_r23, b1_r23, b2_r23, c1l_r23, c2l_r23, delgl_r23, lmax_r23       !parameters for therm. cond. of R23

        !TK1 model
        Integer, dimension(30):: npnum_tk1,npdenom_tk1,nexp_tk1,nspare_tk1
        double precision, dimension(30,100):: a_tk1, tsum_tk1, texp_tk1, dsum_tk1, dexp_tk1

    end type type_thermal_conductivity

    type type_dielectric_constant

        logical:: de_read

        character*255, dimension(30):: demodel
        double precision, dimension(30):: tmin_de
        double precision, dimension(30):: tmax_de
        double precision, dimension(30):: pmax_de
        double precision, dimension(30):: rhomax_de
        double precision, dimension(30):: ref_temp_de
        double precision, dimension(30):: ref_dens_de
        double precision, dimension(30):: ref_press_de
        integer,          dimension(30):: term_num1_de
        integer,          dimension(30):: term_num2_de
        integer,          dimension(30):: term_num3_de
        integer,          dimension(30):: term_num4_de
        integer,          dimension(30):: term_num5_de
        integer,          dimension(30):: term_num6_de

        double precision, dimension(30,100)::coeffde,texpde,dexpde,pexpde

    end type type_dielectric_constant



    type type_lit_ref

        character (len=2000), dimension(30) :: lit_ref_res
        character (len=2000), dimension(30) :: lit_ref_de
        character (len=2000), dimension(30) :: lit_ref_tcx
        character (len=2000), dimension(30) :: lit_ref_tcx_aux
        character (len=2000), dimension(30) :: lit_ref_stn
        character (len=2000), dimension(30) :: lit_ref_eta
        character (len=2000), dimension(30) :: lit_ref_eta_aux
        character (len=2000), dimension(30) :: lit_ref_mlt
        character (len=2000), dimension(30) :: lit_ref_sbl
        character (len=2000) :: lit_ref_mix
        character (len=2000), dimension(30) :: lit_ref_ideal
        character (len=2000), dimension(30) :: lit_ref_pv
        character (len=2000), dimension(30) :: lit_ref_dv
        character (len=2000), dimension(30) :: lit_ref_dl
        character (len=2000), dimension(30) :: lit_ref_saft

    end type type_lit_ref

    type eos_coefficients
        !type module_eos_coefficients

        !implicit none
        !save

        !Number of coefficients and exponents
        integer, dimension(30) :: nreg, ncrt, nna     !number of regular terms, of gaussian terms,
        !  of nonanalytic terms
        integer, dimension(30) :: nlreg, nlcrt, nlna, nr_nf_col  !number of columns of data
        character(len=3), dimension(30) :: cppcheck   !used to compare if CPP is described in file
        integer, dimension(30) :: I_pol, I_exp, I_GBS, I_NA, nr_nf! number of nonanalytic terms

        !Data Arrays:
        double precision, dimension(:,:) , allocatable :: ni     !coefficients
        double precision, dimension(:,:) , allocatable :: ti    !temp exponents
        double precision, dimension(:,:) , allocatable :: di    !density exponents
        integer, dimension(:,:) , allocatable :: p_i    !density exponents in exponential term
        double precision, dimension(:,:) , allocatable :: gama   !factor in exponential term, mostly 1

        double precision, dimension(:,:) , allocatable :: eps    !epsilon in gaussian term
        double precision, dimension(:,:) , allocatable :: beta   !beta in gaussian term
        double precision, dimension(:,:) , allocatable :: gam    !gamma in gaussian term
        double precision, dimension(:,:) , allocatable :: eta    !eta in gaussian term
        double precision, dimension(:,:) , allocatable :: tli    !exponent always 2 (not needed?)
        double precision, dimension(:,:) , allocatable :: pli    !exponent always 2 (not needed?)

        double precision, dimension(:,:) , allocatable :: etana  !nonanalyticterm 1
        double precision, dimension(:,:) , allocatable :: eidna  !nonanalyticterm 2
        double precision, dimension(:,:) , allocatable :: eitna  !nonanalyticterm 3

        logical:: hard_sphere(30)
        double precision, dimension(:), allocatable:: h1, h2, h3, h4, dof2 !packing fraction parameters for
        !hard sphere terms


        !new term vector
        double precision, dimension(:) , allocatable:: x_p

        !end type module_eos_coefficients
    end type eos_coefficients

    type type_ge
        !*****************************************************************************************************
        !PRELIMINARY CODE
        !CAN LATER BE MOVED TO OTHER FILES!!!
        !ANDREAS J鬍ER, SEPTEMBER 2016
        !*****************************************************************************************************
        !type module_gE
        double precision, dimension(nderivs) :: gE_C        !Combinatorial part of gE and derivatives with respect to tau and delta
        double precision, dimension(nderivs) :: gE_R        !Residual part of gE and derivatives with respect to tau and delta
        double precision, dimension(nderivs, 30) :: ln_gamma_C       !Combinatorial activity coefficients and derivatives with respect to delta and tau
        double precision, dimension(nderivs, 30) :: ln_gamma_R       !Residual activity coefficients and derivatives with respect to delta and tau
        double precision, dimension(nderivs, 30) :: gE_C_dxa        !Combinatorial part of gE and derivatives with respect to tau and delta and xa
        double precision, dimension(nderivs, 30) :: gE_R_dxa        !Residual part of gE and derivatives with respect to tau and delta and xa
        double precision, allocatable:: ln_gamma_C_dxa(:,:,:)       !Combinatorial activity coefficients and derivatives with respect to delta and tau and xa
        double precision, allocatable:: ln_gamma_R_dxa(:,:,:)
        double precision, allocatable:: gE_C_dxadxb(:,:,:)                 !Combinatorial part of gE and derivatives with respect to tau and delta and xa and xb
        double precision, allocatable:: gE_R_dxadxb(:,:,:)                 !Residual part of gE and derivatives with respect to tau and delta and xa and xb
        double precision, allocatable:: ln_gamma_C_dxadxb(:,:,:,:)       !Combinatorial activity coefficients and derivatives with respect to delta and tau and xa and xb
        double precision, allocatable:: ln_gamma_R_dxadxb(:,:,:,:)       !Residual activity coefficients and derivatives with respect to delta and tau and xa and xb
        !end type module_gE

        !Variables to safe results in order to !save calculation time
        double precision :: Temp_prev!, Dens_prev                   !Temperature and Density that the routine was last called with
        integer, dimension(15):: GETDER_prev                        !array specifier to indicate, which derivatives were calculated last time the routine was called
        double precision, dimension(30):: molfrac_prev              !composition the routine was last called with
        integer, dimension(30):: Eq_type_prev                       !Equation type that the routine was last called with
        integer :: mixtype_prev                                     !mixtype the routine was last called with
        integer :: C_or_R_prev
        double precision, dimension(15) :: gE_C_prev                 !Combinatorial part of gE and derivatives with respect to tau and delta
        double precision, dimension(15) :: gE_R_prev                 !Residual part of gE and derivatives with respect to tau and delta
        double precision, dimension(15, 30) :: ln_gamma_C_prev       !Combinatorial activity coefficients and derivatives with respect to delta and tau
        double precision, dimension(15, 30) :: ln_gamma_R_prev       !Residual activity coefficients and derivatives with respect to delta and tau

        double precision :: Temp_dxa_prev!, Dens_dxa_prev               !Temperature and Density that the routine was last called with
        integer, dimension(15):: GETDER_dxa_prev                        !array specifier to indicate, which derivatives were calculated last time the routine was called
        double precision, dimension(30):: molfrac_dxa_prev              !composition the routine was last called with
        integer, dimension(30):: Eq_type_dxa_prev                       !Equation type that the routine was last called with
        integer :: mixtype_dxa_prev                                     !mixtype the routine was last called with
        integer :: C_or_R_dxa_prev
        double precision, dimension(15,30) :: gE_C_dxa_prev             !Combinatorial part of gE and derivatives with respect to tau and delta
        double precision, dimension(15,30) :: gE_R_dxa_prev             !Residual part of gE and derivatives with respect to tau and delta
        double precision, dimension(15,30,30) :: ln_gamma_C_dxa_prev    !Combinatorial activity coefficients and derivatives with respect to delta and tau
        double precision, dimension(15,30,30) :: ln_gamma_R_dxa_prev    !Residual activity coefficients and derivatives with respect to delta and tau

        double precision :: Temp_dxadxb_prev!, Dens_dxadxb_prev             !Temperature and Density that the routine was last called with
        integer, dimension(15):: GETDER_dxadxb_prev                         !array specifier to indicate, which derivatives were calculated last time the routine was called
        double precision, dimension(30):: molfrac_dxadxb_prev               !composition the routine was last called with
        integer, dimension(30):: Eq_type_dxadxb_prev                        !Equation type that the routine was last called with
        integer :: mixtype_dxadxb_prev                                      !mixtype the routine was last called with
        integer :: C_or_R_dxadxb_prev
        double precision, dimension(15,30,30) :: gE_C_dxadxb_prev              !Combinatorial part of gE and derivatives with respect to tau and delta
        double precision, dimension(15,30,30) :: gE_R_dxadxb_prev              !Residual part of gE and derivatives with respect to tau and delta
        double precision, dimension(15,30,30,30) :: ln_gamma_C_dxadxb_prev     !Combinatorial activity coefficients and derivatives with respect to delta and tau
        double precision, dimension(15,30,30,30) :: ln_gamma_R_dxadxb_prev     !Residual activity coefficients and derivatives with respect to delta and tau
    end type type_ge


    type type_sea

        !The Array with the constants
        double precision, dimension(8,8,7) :: GIJK_list    !done
        !The range of Validity for the Gibbs Equation for Seawater
        double precision :: tmin_seawater
        double precision :: tmax_seawater
        double precision :: pmin_seawater
        double precision :: pmax_seawater
        double precision :: salmin_seawater
        double precision :: salmax_seawater!delta_testMB
        !For Developer application ;)
        !double precision :: delta_D
        !double precision :: delta_T
        double precision :: Regula_temp
        double precision :: regula_press
        double precision :: regula_rho
        double precision :: salinity
        double precision :: seap
        double precision :: seap_save
        double precision :: wm_sea
        double precision :: wm_salt
        double precision :: x_s
        double precision :: x_w
        double precision :: sal_m
        double precision :: molcoef
        double precision :: seapmin_seawater
        double precision, dimension(30,5) :: dummyarray, dummyarray1
        !double precision :: watermix_roh_liq
        !double precision :: vapfrac_reac
        !double precision :: seap_save!, dpdxiv
        !double precision, dimension(30) :: wm_reac
        !double precision, dimension(30) ::wm_save!, seastartcomp, watermix_x_vap, watermix_x_liq
        !logical:: seawater, sal, seaphase, seacalc,  seawatercalled, seawaterflash
        !integer :: salpos, seapos

    end type type_sea

    type type_el
        !subtype for electrolytes, Benedikt 2019

        integer :: solpos
        integer :: n_salts, n_fluids
        integer :: salt
        integer :: saltpos
        integer :: modelflag
        integer :: functionflag
        integer :: betaflag
        integer, dimension(30) :: z_salt, mapping, ref_molality, eq_type_el, eq_types_in
        double precision :: tred, tred2, pred, pred2, pi, Na, e_const, eps0, es, kb, vm, vx, alpha1, alpha2, zm, t0, p0, el_permit, zx, a_part, bc_part, press, deybe, m_ref, x_salt, mol_in, wm_salt, n2, g_Exold ,mix_old
        double precision, dimension(30)::mw , m, pr, tr, b, molality, moles_input, cp0_param
        double precision, dimension(4,3,7) :: bijk
        double precision, dimension(4,30) :: pitzer_param, g0_param, lit_param
        double precision, dimension(4,10) :: ref_param
        character(30), dimension(30) :: salts, fluids_input, fluids_wo_salts
        logical :: resorted, mixedelectrolytes, fitting, brine

        logical :: strongel

    end type type_el

    type type_dew
        logical:: de_read

        character*255, dimension(30):: demodel
        double precision, dimension(30):: tmin_de
        double precision, dimension(30):: tmax_de
        double precision, dimension(30):: pmax_de
        double precision, dimension(30):: rhomax_de
        double precision, dimension(30):: ref_temp_de
        double precision, dimension(30):: ref_dens_de
        double precision, dimension(30):: ref_press_de
        integer,          dimension(30):: term_num1_de
        integer,          dimension(30):: term_num2_de
        integer,          dimension(30):: term_num3_de
        integer,          dimension(30):: term_num4_de
        integer,          dimension(30):: term_num5_de
        integer,          dimension(30):: term_num6_de

        double precision, dimension(30,100)::coeffde,texpde,dexpde,pexpde
    end type type_dew

    type type_des
        logical:: de_read

        character*255, dimension(30):: demodel
        double precision, dimension(30):: tmin_de
        double precision, dimension(30):: tmax_de
        double precision, dimension(30):: pmax_de
        double precision, dimension(30):: rhomax_de
        double precision, dimension(30):: ref_temp_de
        double precision, dimension(30):: ref_dens_de
        double precision, dimension(30):: ref_press_de
        integer,          dimension(30):: term_num1_de
        integer,          dimension(30):: term_num2_de
        integer,          dimension(30):: term_num3_de
        integer,          dimension(30):: term_num4_de
        integer,          dimension(30):: term_num5_de
        integer,          dimension(30):: term_num6_de

        double precision, dimension(30,100)::coeffde,texpde,dexpde,pexpde
    end type type_des




    double precision, PARAMETER :: FOR_D_QNAN = TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)

    type fnrderivs_cache_t
        ! simple cache for storing the last values for tau and del calculated in fnrderivs

        ! flag and old value for del
        !logical :: del_ok
        double precision :: del_old = FOR_D_QNAN

        !cached data
        double precision, dimension(0:6) :: delp_i
        double precision :: logdel


        ! flag and old value for tau
        !logical :: tau_ok
        double precision :: tau_old = FOR_D_QNAN

        !cached data
        double precision :: tau2, tau3
        double precision :: logtau



    end type fnrderivs_cache_t


    type depfuncfnr_cache_t
        ! simple cache for storing the last values for tau and del calculated in fnrderivs

        ! full-cache for mixtures with up to 5 components (needs 10 entries for the upper triangle of a diag-matrix (4+3+2+1))
        ! index 0 is used for non-cached data to have uniform access

        ! flag and old value for del
        !logical :: del_ok
        double precision, dimension(0:10) :: del_old = FOR_D_QNAN


        !cached data
        double precision, dimension(25,0:10) :: del_pow_di_depf


        ! flag and old value for tau
        !logical :: tau_ok
        double precision, dimension(0:10) :: tau_old = FOR_D_QNAN

        !cached data
        double precision, dimension(25,0:10) :: tau_pow_ti_depf



    end type depfuncfnr_cache_t


    type cross
        double precision :: delta_1,alpha,beta,gama,b0,b2,del_1
        double precision , allocatable :: param(:)

        integer,dimension(15) ::  dx = (/ 0,0,0,1,2,1,2,0,3,1,0,1,2,3,4 /)
        integer,dimension(15) ::  dy = (/ 0,1,2,0,0,1,1,3,0,2,4,3,2,1,0 /)
    end type



    type type_gl
        
        type(type_model_info):: m_info
        
        type(type_state_point) :: s_p
        !type module_ref

        type(fnrderivs_cache_t) :: fnrderivs_cache
        type(depfuncfnr_cache_t) :: depfuncfnr_cache
        !type(func):: ff

        !if molfractions of a mixture are zero this array contains the positions of the original mixture of the fluids that are present
        integer,allocatable :: comp_map(:)
        !zero_comp: false: all specified molfractions are > 1d-14, true: one or more specified molfraction is < 1.d-14
        logical :: zero_comp, zero_spec
        !storage for original input
        character (30), allocatable :: fluids_zero_orig(:)
        double precision, allocatable :: moles_zero_orig(:), moles_zero_orig_spec(:)
        integer, allocatable :: EOS_indicator_zero_orig(:)

        !implicit none
        !save
        !boolean whether to calculate reference point
        logical :: calc_ref
        !boolean whether it has been set
        logical :: ref_set
        !integer to
        ! - set the range of validity of transport properties (1: viscosity, 2: thermal conductivity, 3: dielectric constant)
        ! - allocate the variable for
        integer:: transport !(1: viscosity, 2: thermal conductivity, 3: dielectric constant)

        logical :: VLE_needed !set to false if only critical/triple parameters or acentric factor or MW are needed

        !end type module_ref
        logical :: gemix

        !type module_eos_coefficients
        type(eos_coefficients), allocatable:: eos_coeff
        !end type


        !type module_fluid_parameters

        !implicit none
        !save

        !flag variable to check whether a fluid has already been loaded or not
        !on 1st call string is empty
        !on later calls the variable contains the string of the previous call
        character(30), dimension(30) :: components_old
        !integer :: write(*,*)
        logical :: same_components
        character(255) :: path_arg

        double precision, dimension(2) :: inputprops
        double precision, dimension(30) :: moles_read

        !mixture specification:
        character(30), dimension(30) :: components      !Array of components in the mixture
        double precision, dimension(30) :: molfractions     !Array of molefractions of the components
        integer :: ncomp , ncomp_old, nphases                                 !Number of Components and phases in the Mixture

        !DEC$ IF DEFINED(HCgen)
        !Name of Substance:
        character(len=255), dimension(30) :: substfullname
        character(len=12),  dimension(30) :: substshortname
        character(len=255), dimension(30) :: substchemname
        character(len=255), dimension(30) :: substsynonym

        character(len=12), dimension(30) :: substcasnr
        !DEC$ ELSE
        character(len=255), dimension(:), allocatable :: substfullname
        character(len=12),  dimension(:), allocatable :: substshortname
        character(len=255), dimension(:), allocatable :: substchemname
        character(len=255), dimension(:), allocatable :: substsynonym

        character(len=12), dimension(:), allocatable :: substcasnr
        !DEC$ ENDIF

        !Critical Point data:
        double precision, dimension(30) :: tc      !temperature at critical point [K]
        double precision, dimension(30) :: pc      !pressure at critical point [MPa]
        double precision, dimension(30) :: rhoc    !density at critical point [mol/m設

        !Triple Point data:
        double precision, dimension(30) :: ttp     !temperature at triple point [K]
        double precision, dimension(30) :: ptp     !pressure at triple point [MPa], calculated from
        !    vapor pressure equation
        double precision, dimension(30) :: ptpmod  !pressure at triple point [MPa] * 0.1
        double precision, dimension(30) :: rhotp   !density at triple point [mol/m設

        !Normal Boiling Point data:
        double precision, dimension(30) :: tnbp    !temperature at normal boiling point [K]

        !Further data:
        double precision, dimension(30) :: wm      !molecular weight [kg/mol]
        double precision, dimension(30) :: accen   !acentric factor [-]
        double precision, dimension(30) :: dipole  !dipole moment [debye]

        !Limits:
        double precision, dimension(30) :: tminfluid      !lower temerature limit [K]
        double precision, dimension(30) :: tmaxfluid      !upper temerature limit [K]
        double precision, dimension(30) :: pmaxfluid      !upper pressure limit [MPa]
        double precision, dimension(30) :: rhomaxfluid    !maximum denisity [mol/m設
        logical :: hold_limits                            !whether to ignore EOS limits or not

        !Gaskonstant:
        double precision, dimension(30) :: Req          !Gasconstant used in RefProp Equation [J/mol-K]

        !Factors for switching between SI units and reduced units
        double precision :: Factor
        double precision:: factorpress
        double precision:: factortrans
        double precision:: factorrbwr
        double precision:: factor_VS5eta
        double precision:: factor_VS5slj

        character(len=12) :: inptorig

        logical:: phaseboundary

        logical, dimension(30) :: pvexist           ! ancillary equation for vapor pressure available
        logical, dimension(30) :: meltexist         ! ancillary equation for melting pressure available

        !New variables to distingiush between different EOS used
        integer, dimension(30):: Eq_type     !This variable indicates which EOS type is used
        !  for the specific fluid:
        !1: Helmholtz-EOS
        !2: SRK
        integer:: Mix_type                   !1: Lorentz-Berthelot or modified Helmholtz mixing rules
        !  (All types of EOS might be integrated to Helmholtz EOS
        !  and used with these mixing rules!!)
        !2: SRK-Mixing rules (Only SRK equations for ALL fluids!!!)

        logical :: twophasecalc              !flag for density iteration assigning a twophase iteration,
        !  which means that the solution doesn't have to meet all
        !  physical criteria

        logical :: bad_mixmodel             !flag for density iteration: physical solution doesn't have to meet all physical criteria

        double precision :: factor_decide4vap

        integer, dimension(5) :: phase_id

        !Parameters of the ideal isochoric heat capacity
        !Thu Febuary 2013
        double precision, dimension (30):: A_cv0
        double precision, dimension (30):: B_cv0
        double precision, dimension (30):: C_cv0
        double precision, dimension (30):: D_cv0
        double precision, dimension (30):: E_cv0
        double precision, dimension (30):: F_cv0
        double precision, dimension (30):: G_cv0


        !Parameters of the ideal isobaric heat capacity for the cp0-modell Joback (1986)
        !The parameters are the sum of the influences and number of the molecules of the substances
        !Thu 2014/05
        double precision:: cp0_A
        double precision:: cp0_B
        double precision:: cp0_C
        double precision:: cp0_D

        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !! For some fluids the ideal part and the residual part are reduced with different R
        !! When a property is calculated, the ideal part is multiplied with the same R as the residual
        !! part in REFPROP (e.g. Prop=(ai+ar)*R_res). When it comes to mixtures this problem gets even
        !! worse, since fluids in the mixture might have different Rs for residual parts and ideal parts
        !! Thus (in agreement with the GERG-software) the ideal part and the residual part are
        !! transformed to the same R in the sub_file_input routine by multiplying all linear
        !! coefficients with the equation specific R and dividing by the universal R
        !! Andreas, November 2013
        !! NOT YET IMPLEMENTED!!!!! (REASON FOR NOT IMPLEMENTING THIS AT THE MOMENT:
        !! NOT ONLY THE COEFFICIENTS NEED TO BE CORRECTED, BUT ALSO SOME "1"s
        !! (FOR EXAMPLE: cp = 1 +... = R_old/R_new + ...)
        !logical :: test_mode = .false.
        !  test_mode = true  -- >  Pure fluids are calculated consistently to REFPROP
        !  (ignoring the "R"-problem). USE FOR TESTING ONLY
        !  test_mode = false -- >  The "R"-correction is applied. In this case the results are consistent
        !  to the GERG-software
        !double precision, parameter :: R_univ = 8.314462    !J/molK
        ! -- >  Universal gas constant used for all fluids
        !!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !end type module_fluid_parameters

        !type module_ancillary_parameters  ! alle Indizes von Theresa getauscht!

        !implicit none
        !save

        !Parameters for temperature iteration using the melting pressure equation

        double precision :: Tsatmelt_old, pmelt_high, melt_p
        logical :: melt_it, p_fake

        !Temperature variables:

        double precision, dimension(30) :: pmeltmintemp !melting pressure lower temperature limit
        double precision, dimension(30) :: pmeltmaxtemp  !melting pressure upper temperature limit

        !Temperature and Density for Reducing:
        double precision, dimension(30) :: vpred     !reducing (red.) parameter (par.) for vapor pressure
        double precision, dimension(30) :: tvpred    !temperature red. par. for vapour pressure
        double precision, dimension(30) :: dlred     !red. par. for saturated liquid density [mol/m設
        double precision, dimension(30) :: tdlred    !temperature red. par. for saturated liquid density
        double precision, dimension(30) :: dvred     !red. par. for saturated vapour density [mol/m設
        double precision, dimension(30) :: tdvred    !temperature red. par. for saturated vapour density
        double precision, dimension(30) :: pmeltred  !red. par. for melting pressure
        double precision, dimension(30) :: tpmeltred !temperature red. par. for melting pressure
        double precision, dimension(30) :: psubred   !red. parameter for sublimation pressure
        double precision, dimension(30) :: tpsubred  !temperature red. par. for sublimation pressure

        !Number of coefficients and exponents
        integer, dimension(30) :: nvpcoeff           !number of vapour pressure coefficients
        integer, dimension(30) :: ndlcoeff           !number of saturated liquid density coefficients
        integer, dimension(30) :: ndvcoeff           !number of saturated vapour density coefficients
        !number of melting pressure coefficients
        integer, dimension(30) :: npmeltcoeff1, npmeltcoeff2, npmeltcoeff3,npmeltcoeff4,npmeltcoeff5
        !number of sublimation pressure coefficients
        integer, dimension(30) :: npsubcoeff1, npsubcoeff2, npsubcoeff3

        !Type of Equation
        integer, dimension(30) :: vptype                     !vapour pressure equation type
        integer, dimension(30) :: dltype                     !saturated liquid density equation type
        integer, dimension(30) :: dvtype                     !saturated vapour density equation type
        integer, dimension(30) :: pmelttype                  !meling pressure equation type
        integer, dimension(30) :: psubtype                   !sublimation pressure equation type

        !Data Arrays:
        double precision, dimension(:,:),allocatable :: vpcoeff    !vapour pressure coefficients
        double precision, dimension(:,:),allocatable :: vpexp      !vapour pressure exponents
        double precision, dimension(:,:),allocatable :: dlcoeff    !saturated liquid density coefficients
        double precision, dimension(:,:),allocatable :: dlexp      !saturated liquid density exponents
        double precision, dimension(:,:),allocatable :: dvcoeff    !saturated vapour density coefficients
        double precision, dimension(:,:),allocatable :: dvexp      !saturated vapour density exponents
        double precision, dimension(:,:),allocatable :: pmeltcoeff !melting pressure coefficients
        double precision, dimension(:,:),allocatable :: pmeltexp   !melting pressure exponents
        double precision, dimension(:,:),allocatable :: psubcoeff  !sublimation pressure coefficients
        double precision, dimension(:,:),allocatable :: psubexp    !sublimation pressure exponents

        !end type module_ancillary_parameters

        !type module_general_eos_parameters ! alle Indizes von Theresa getauscht!

        !implicit none
        !save

        !Temperature and Density for Reducing:
        double precision, dimension(30) :: tred   !temperature for reducing [K]
        double precision, dimension(30) :: rhored !density for reducing [mol/m設
        double precision :: tredmix               !calculated temperature for reducing in mixtures [K]
        double precision :: rhoredmix             !calculated density for reducing in mixtures [mol/m設

        !Reference State:
        character(len=3), dimension(30):: refstate !can be others, ...
        double precision, dimension(30) :: tref    !temperature choosen for reference point [K]
        double precision, dimension(30) :: pref    !pressure choosen for reference point [MPa]
        double precision, dimension(30) :: rhoref  !density choosen for reference point [mol/m設
        !(not read from refprop file, instead calculated in subroutine)
        double precision, dimension(30) :: href    !enthalpy at reference point [J/mol]
        double precision, dimension(30) :: sref    !entropy at reference point [J/(mol K)]
        logical :: ref
        !end type module_general_eos_parameters



        !type module_ideal_gas_coefficients ! alle Indizes von Theresa getauscht!

        !implicit none
        !save

        !c1 and c2 not from Refprop file, they will be calculated at the end of the read_file routine
        double precision, dimension(30) :: c1, c2

        !Temperature and Density for Reducing:
        double precision, dimension(30) :: tcp0red     !temperature reducing parameter for cp0
        double precision, dimension(30) :: cp0red      !reducing parameter for cp0

        !Number of coefficients and exponents
        integer, dimension(30) :: ncp0cosh,ncp0sinh   !number of cosh and sinh coefficients of cp0 data
        !number of polynom and planck-einstein coefficients of cp0 data
        integer, dimension(30) :: ncp0poly,ncp0pl
        integer, dimension(30) :: ncp0el              !special term for the the ideal gas of air.ppf
        integer, dimension(30) :: ncp0c           !number of constant coefficients of alfa0 data
        !number of polynom and planck-einstein coefficients of cp0 data
        integer, dimension(30) :: ncp0log,ncp0,ncp0logexp


        !DEC$ IF DEFINED(HCgen)
        !Data Arrays:
        double precision, dimension(20,30) :: cp0coeff    !cp0 coefficients
        !cp0 exponents = temperature exponents ti* for cp0 equation
        double precision, dimension(20,30) :: cp0exp
        double precision, dimension(20,30) :: cp0hyp1      !cp0 parameter for sinh, cosh of cpp models
        double precision, dimension(20,30) :: cp0hyp2      !cp0 parameter for sinh, cosh of cpp models
        double precision, dimension(20,30) :: cp0hyp3      !cp0 parameter for sinh, cosh of cpp models
        !DEC$ ELSE
        !Data Arrays:
        double precision, dimension(:,:) ,allocatable:: cp0coeff    !cp0 coefficients
        !cp0 exponents = temperature exponents ti* for cp0 equation
        double precision, dimension(:,:),allocatable :: cp0exp
        double precision, dimension(:,:),allocatable :: cp0hyp1      !cp0 parameter for sinh, cosh of cpp models
        double precision, dimension(:,:),allocatable :: cp0hyp2      !cp0 parameter for sinh, cosh of cpp models
        double precision, dimension(:,:),allocatable :: cp0hyp3      !cp0 parameter for sinh, cosh of cpp models
        !DEC$ ENDIF

        logical, dimension(30) :: cpmodel, phkmodel, phmodel               !variables to check which ideal part model is used

        !end type module_ideal_gas_coefficients

        !type module_mixture_parameters

        !implicit none
        !save

        !arrays for reducing functions coefficients for mixtures
        double precision, dimension(30,30) :: rfbetat
        double precision, dimension(30,30) :: rfbetarho
        double precision, dimension(30,30) :: rfgammat
        double precision, dimension(30,30) :: rfgammarho

        !----------------------------------------------------------------------------------------------
        !DEPARTURE FUNCTION COEFFICIENT STORAGE VARIABLES CHANGED!!! Andreas Aug. 2010
        double precision, dimension(:,:,:), allocatable:: dfn
        double precision, dimension(:,:,:), allocatable:: dfd
        double precision, dimension(:,:,:), allocatable:: dft
        double precision, dimension(:,:,:), allocatable:: dfl
        double precision, dimension(:,:,:), allocatable:: dfp
        double precision, dimension(:,:,:), allocatable:: dfeta
        double precision, dimension(:,:,:), allocatable:: dfeps
        double precision, dimension(:,:,:), allocatable:: dfbeta
        double precision, dimension(:,:,:), allocatable:: dfgamma
        double precision, dimension(:,:,:), allocatable:: dfgeta    !eta in gaussian term
        double precision, dimension(:,:,:), allocatable:: dfgbeta   !beta in gaussian term
        double precision, dimension(:,:,:), allocatable:: dfggam    !gamma in gaussian term
        double precision, dimension(:,:,:), allocatable:: dfgeps    !epsilon in gaussian term

        integer, dimension(:,:), allocatable :: dfd_coeff_structure
        integer, dimension(:,:), allocatable :: dfl_coeff_structure


        double precision, dimension(30,30):: Fij
        integer, dimension(30,30):: dfpol
        integer, dimension(30,30):: dfexp
        integer, dimension(30,30):: dfgau

        !Binary mixing parameter for quadratic mixing rules for the residual reduced Helmholtz energy

        !double precision, dimension(30,30):: kij_Helm

        !Variables for cross coefficients of new mixture model, Andreas J輍er July 2016
        !double precision, dimension(30,30,50) :: equal_terms    !Go through all terms and identify equal terms for all components
        !The information is always stored in the following way: equal_terms(fluid1, fluid2, nr_of_equal_terms) = termnr_fluid2
        !Example:
        !equal_terms(1,2,1) = 3 -> The third term of fluid 2 is also used in fluid 1 (First equal term found)
        !equal_terms(2,3,10) = 8 -> The eights term of fluid 3 is also used in fluid 2 (Tenth equal term found)
        !For all equal terms, only one binary interaction parameter Helm_k_nij is needed
        !integer:: nr_of_equal_terms(30,30)                      !Stores how many equal terms have been found per binary mixture

        !double precision, dimension(:,:,:), allocatable :: Helm_k_nij    !The binary interaction parameter per term for fluid i and j: Helm_k_nij(fluid1, fluid2, termnr)

        !Variables for the new mixture model
        double precision, dimension(:,:,:,:), allocatable :: mix_nij    !coefficients
        double precision, dimension(:,:,:), allocatable :: mix_tij    !temp exponents
        double precision, dimension(:,:,:), allocatable :: mix_dij    !density exponents
        integer, dimension(:,:,:), allocatable :: mix_p_ij    !density exponents in exponential term
        double precision, dimension(:,:,:), allocatable :: mix_gama   !factor in exponential term, mostly 1

        !double precision, dimension(30,30,40) :: mix_eps    !epsilon in gaussian term
        !double precision, dimension(30,30,40) :: mix_beta   !beta in gaussian term
        !double precision, dimension(30,30,40) :: mix_gam    !gamma in gaussian term
        !double precision, dimension(30,30,40) :: mix_eta    !eta in gaussian term
        !double precision, dimension(30,30,40) :: mix_tli    !exponent always 2 (not needed?)
        !double precision, dimension(30,30,40) :: mix_pli    !exponent always 2 (not needed?)

        integer, dimension(:,:), allocatable :: mix_nreg               !number of regular terms (of binary mixture interaction terms)


        !end type module_mixture_parameters

        type(type_uncertainty), allocatable :: uncty

        !type module_psrk

        !implicit none
        !save

        !check if there is enough data for calculations:
        logical :: psrk_ok



        !double precision, dimension(100,33) :: ufc
        !   up to 100 lines: group X, group Y, ...
        !   up to 33 columns: groupnr, subgroupnr, Qk, amount of group X in comp A, B, C, ...

        !number of different groups:
        integer :: ngroup
        !number of different groups including different subgroups:
        integer :: ngroup_nsubgroup

        !coefficients:
        double precision, dimension(30,3) :: ccoeff !30 lines: component, 3 columns: c1, c2,c3

        !double precision, dimension(100,100) :: atcoeff
        !double precision, dimension(100,100) :: btcoeff
        !double precision, dimension(100,100) :: ctcoeff

        !end type module_psrk

        !type module_asso

        !implicit none
        !save

        !check if association terms exsists
        logical, dimension(30) :: assoexist

        !model:
        integer, dimension(30) :: assomodel

        !coefficients:
        double precision, dimension(30) :: assoc      !multiplicative coefficient
        double precision, dimension(30) :: assovolred !reduced normalization volume
        double precision, dimension(30) :: assovolint !Len.-Jones segment diameter*volume of interaction
        double precision, dimension(30) :: assoenergy !association energy

        !end type module_asso

        !type module_VLE
        !--------------------------------------------------------------------------------
        ! In this type veraiables are stored that are needed in several
        ! routines for the calculation of the phase equilibrium.
        ! Storing variables in this type is important especially for densities
        ! that have to be calculateed iteratively.
        !--------------------------------------------------------------------------------
        ! J. Gernert, Jan. 2011

        double precision, dimension(30):: rho_TP_pure   ! Pure fluid densities at given T, p
        ! Pure fluid fugacity coefficients at given T and a density calculated from T and mixture p
        double precision, dimension(30):: fugco_pure
        double precision:: rho_vap                      ! Mixture density, vapor phase
        double precision:: rho_liq                      ! Mixture density, liquid phase
        ! In case of three phase equilibria, two liquid densities are necessary
        double precision:: rho_liq1, rho_liq2, rho_liq3
        !integer, parameter:: imax = 800
        ! arrays for storing points on the phase boundary, calculated in the routine 'phasenv'
        double precision, dimension(imax):: T_pts, p_pts, rholiq_pts, rhovap_pts
        double precision, dimension(imax, 30):: x_pts
        integer, dimension(imax)::pointID  ! ID parameter for special points:
        ! normal pts: 0, crit. pt: 1, cricondentherm: 2, cricondenbar: 3
        ! specified temp.:4, specified pressure: 5
        integer:: phasenv_pts       ! number of points calculated on the phase boundary

        logical:: savebounds_p
        double precision:: min_factor
        double precision:: pmin_old, pmax_old
        !end type module_VLE

        !type module_phasedet_pure

        !implicit none
        !save

        double precision:: rho_it      ! density internaly iterated
        double precision:: p_it        ! pressure from which density was iterated - >  control parameter

        !end type module_phasedet_pure





        !type module_solids
        !--------------------------------------------------------------------------------
        ! In this type variables are stored that are needed for
        ! several solid routines.
        !--------------------------------------------------------------------------------
        ! A. J輍er, May. 2012

        !implicit none
        !save

        !Contains a list of all substances that form hydrates and are impemented to this software
        !If a known hydrate former is not in this list, calculation of that types of hydrate is not
        !possible at the moment
        character (30), dimension (:), allocatable :: Hydrate_list

        !In this list all hydrate formers that could be identified in the fluidlist are stored
        character (30), dimension (30) :: Hydrate_formers
        !In this list the position of all hydrate formers according to the list "Hydrate_formers"
        !   in the fluid list is stored
        !Example for the use of the hydrate variables:
        !Given mixture:         Fluid = "water";"Decane";Methane";"CO2";"Nonane"
        !Possible Hydrates:     Hydrate_list = "water";"CO2";Methane";
        !Hydrate_formers        Hydrate_formers = "water";"Methane";"CO2" (Order like in the fluid vector)
        !Hydrate_pos            Hydrate_pos = 1;3;4 (Positions in the fluid vector)
        !integer, dimension(30):: Hydrate_pos

        integer:: nrofhydrateformers                            !The number of how many components in the mixture form hydrates
        character (30), dimension (30) :: Fluidlist_hydrate     !This is a changed fluid list for internal use only
        double precision, dimension(30) :: moleslist_hydrate    !Sorted molfractions according to fluidlist_hydrate
        integer, dimension(30):: mapping                        !In this list the position of all substances of the changed list "Fluidlist_hydrate" is stored
        integer:: hdrt_structure_stable                         !The Stable hydrate strucutre type is saved here.

        character (30), dimension (:), allocatable :: Solid_list            !Contains a list of all substances for which solid state equations are available
        character (30), dimension (30) :: solid_subst           !In this list all components for which solid equations are available that could be identified in the fluidlist are stored
        integer, dimension(30):: solid_indicator                !This is a list, containing information on which position(s) substances with solid equations availabe stand in Fluidlist_hydrate

        !Example for the use of the solid variables:
        !Given mixture:         Fluid = "water";"Decane";Methane";"CO2";"Nonane"
        !Possible Hydrates:     Fluidlist_hydrate = "water";"CO2";"methane";"Decane"; Nonane"    < -- This list is used internally as fluid list
        !solid_list             solid_list = "water";"CO2"
        !solid_indicator        solid_indicator = 1;2;0,0,0,0    < -- the numbers refer to solid_list
        integer:: solid_pos                                     !Position of the solid in the vector Fluidlist_hydrate


        integer, dimension(2) :: solidtype                      !This variable might be used to indicate which solids occur in the phase equilibrium, this has to be figured out in the routine PhaseDet
        !solidtype(1) -- >  solidtype(1) = 0 -- >  No pure solid phase
        !                 solidtype(1) = solid nr according to solid_list
        !solidtype(2) -- >  solidtype(2) = 0 -- >  No hydrate phase
        !                 solidtype(2) = 1 -- >  Hydrate phase

        integer :: solidtype_akt_phase                          !the variable "solidtype" is the working variable that determines the solid flash calculation
        !In the course of iteration it becomes necessary to !save the information if the solid phase in the
        !phase equilibrium found is solid water or solid co2 (other solids not yet implemented). This information is
        !stored in this variable. It is: SOLID WATER: solidtype_akt_phase = 1 ; SOLID CO2: solidtype_akt_phase = 2
        integer :: solidpos_akt_phase                           !The same as solidtype_akt_phase, but this variable saves the position in the fluid vector where the pure solid
        !former stands


        integer:: nrofsolids                                    !The number of how many components in the mixture form solids

        character (25), dimension (30) ::Fluidlist_solids       !This is a changed fluid list for internal use only

        logical :: check_solid                                     !True -- >  the code checks for the formation of solid phases
        !False -- >  No handling of solid phases in the code

        double precision, dimension(6) :: p_Q_Hyd_2C            !Variable to !save pressures at quadruple points for binary hydrate forming system.
        !This information makes ph and ps flash calculations for binary systems easier (and a little faster)
        !Q(1)       : VLwHIw
        !Q(2)       : VLwLxH  !If system forms LLE
        !Q(3)       : VLxHIx  !If system forms LLE and Solid EOS is available
        !Q(4)       : LwLxHIx !If system forms LLE and solid EOS is available
        !Q(5),Q(6)  : Backup
        double precision, dimension(6) :: T_Q_Hyd_2C           !Variable to !save temperatures at quadruple points for binary hydrate forming system.
        !This information makes pMLTL calculations close to quaduple points safer and faster
        !Q(1)       : VLwHIw
        !Q(2)       : VLwLxH  !If system forms LLE
        !Q(3)       : VLxHIx  !If system forms LLE and Solid EOS is available
        !Q(4)       : LwLxHIx !If system forms LLE and solid EOS is available
        !Q(5),Q(6)  : Backup


        !Limits for solid equations
        !Andreas Feb 2014
        double precision :: pmax_Dryice
        double precision :: Tmin_Dryice
        double precision :: pmax_Waterice
        double precision :: Tmin_Waterice




        !end type module_solids


        !
        !type module_droplet_VLE
        !
        !!implicit none
        !!save
        !
        !logical :: droplet          ! indicator for different VLE algorithm, mueL(pL) = mueV(pV), pL != pV, Laplace equation: pL = pV + 2s/r  (surface tension s, droplet radius r)
        !double precision :: pl      ! pressure inside of the droplet ( >  pv)
        !
        !end type module_droplet_VLE


        type(type_dielectric_constant), allocatable :: de
        type(type_dew), allocatable :: dew
        type(type_des), allocatable :: des


        logical:: stn_read
        type(surface_tension_t), allocatable, dimension(:):: stn
        type(type_viscosity), allocatable :: visco

        !type module_ideal

        !implicit none
        !save

        double precision:: rhoideal

        !end type

        !type module_bwr

        !implicit none
        !save

        logical:: bwr_read
        Logical :: ben_read   !Bender EOS has been read
        Logical :: STA_read   !Starling EOS has been read

        double precision, dimension(30,100):: ncoeff
        double precision, dimension(30):: gama_bwr

        !end type



        type(type_thermal_conductivity), allocatable :: tcx
        !type module_trans_limits

        !implicit none
        !save

        double precision, dimension(30):: tmintrans
        double precision, dimension(30):: tmaxtrans
        double precision, dimension(30):: pmaxtrans
        double precision, dimension(30):: rhomaxtrans

        !end type


        !type module_FNR


        !implicit none
        !save

        !double precision :: TEMP_FNR_OLD, DENS_FNR_OLD, tredmix_OLD, rhoredmix_old
        integer :: nrsubst_old, mix_type_old
        integer, dimension(30):: Eq_type_old
        !!character(len=12), dimension(30) :: comp_FNR_old
        !!double precision, dimension(30) :: mol_FNR_old
        !double precision :: FNR_OLD
        !double precision, dimension(15) :: FNRDER_OLD
        !integer :: zaehler, fnr_c
        !double precision, dimension(70,30) :: delpi_old
        !double precision, dimension (20) :: gauss_term_old, deleps_old, taugam_old
        !double precision, dimension (70) :: reg_term_old
        !double precision :: nacalc_old

        !end type

        !type module_FNR_MIX


        !implicit none
        !save

        !double precision :: TEMP_FNR_MIX_OLD, DENS_FNR_MIX_OLD, tredmix_mix_OLD, rhoredmix_mix_old
        !character(len=12), dimension(30) :: comp_FNR_MIX_old
        !double precision, dimension(30) :: mol_FNR_MIX_old
        !double precision :: FNR_MIX_OLD
        !double precision, dimension(15) :: FNRDER_MIX_OLD
        !integer :: zaehler_mix, cn_fnm, za_fnm, uncomp

        !end type

        !type module_cubic

        !implicit none
        !save

        double precision, dimension(30):: ac_SRK, ac_PR
        double precision, dimension(30):: m_SRK

        double precision, dimension(30):: ai_SRK, bi_SRK, ai_PR, bi_PR        !Constants a and b for every fluid used of the cubic EOS

        double precision, dimension(:,:), allocatable:: kij_SRK, aij_SRK, kij_PR, aij_PR !Mixing rules kij for every binary combination of fluids used
        double precision, dimension(:,:), allocatable:: lij_SRK, bij_SRK, lij_PR, bij_PR !If quadratic mixing rules are used for b, these variables are needed

        double precision:: a_SRK, b_SRK, a_PR, b_PR                           !Mixture parameters for a and b

        double precision:: rhored_SRK, pred_SRK, Tred_SRK                     !Reducing parameters for integrated SRK -- >  CAUTION!! SET TO 1 FOR NOW
        double precision:: rhored_PR, pred_PR, Tred_PR                        !Reducing parameters for integrated PR -- >  CAUTION!! SET TO 1 FOR NOW
        double precision:: R_SRK                                              !Universelle Gaskonstante: R = 8.3144598 J/mol K
        double precision:: R_PR                                               !Universelle Gaskonstante: R = 8.3144598 J/mol K
        double precision, dimension(30):: kappa

        double precision:: Temp_init                                          !This is the last temperature the cubic EOS was initialized with. This is an important variable, which is needed to check whether
        !the EOS has to be initialized or not by comparing this temperature to the actual one

        !!Supporting variables for PR
        !double precision, parameter:: sqrt_2 = 2.D0**0.5
        !double precision, parameter:: consA_PR = 2.D0**0.5 + 1.D0   !consA_PR = 2^0.5 + 1
        !double precision, parameter:: consB_PR = 2.D0**0.5 - 1.D0   !consB_PR = 2^0.5 - 1

        !Andreas J輍er, January 2017
        !Mathias and Copeman function for better representation of the vapor-pressure curve
        double precision, dimension(3,30) :: Cji_cubic          !For the parameter a of the SRK or PR, it is: a = a0 * (1 + C1i*(1-(T/Tci)**0.5) + C1i*(1-(T/Tci)**0.5)**2 + C1i*(1-(T/Tci)**0.5)**3)**2
        !j -> parameter number, i-> fluidnumber

        !!Andreas J輍er, January 2017
        !!PSRK variables
        !double precision, parameter :: A1_PSRK = -0.64663D0     !The parameter A1 for the SRK. A1 = -ln((u+1)/u) with u = 1.1

        !end type module_cubic


        !type module_HelmgE

        !implicit none
        !save

        !Packing fraction u = v/b = ui = vi/bi, where v is the molar volume of the liquid calculated at a reference pressure p0 and b is the covolume of a cubic equation of state
        !Here (like for PSRK) the reference pressure p0 = 101325 Pa = 1 atm was chosen.
        !At these conditions, the average packing fraction on the bubble line for all fluids in REFPROP or TREND is about 1.17
        double precision:: u_pack

        !Covolume
        double precision, dimension(30):: bi_HelmgE        !Covolume for all fluids used in the Helm+gE model
        double precision:: b_HelmgE                        !Covolume of the mixture

        double precision:: rho_mix_ref                      !Reference density of the liquid in the mixture at reference pressure
        double precision, dimension(30):: rho_i_ref         !Reference densities of all fluids in the mixture at reference pressure

        !double precision, parameter :: R_HelmgE = 8.3144598D0    !Universal gas constant

        integer:: LIN_or_LB         !1:LB mixing rules for the reducing parameters, 2: Linear mixing rules for the reducing parameters

        !end type


        !type module_LKP

        !implicit none
        !save

        !LKP-constants
        double precision :: lkp_b1_0
        double precision :: lkp_b2_0
        double precision :: lkp_b3_0
        double precision :: lkp_b4_0
        double precision :: lkp_c1_0
        double precision :: lkp_c2_0
        double precision :: lkp_c3_0
        double precision :: lkp_c4_0
        double precision :: lkp_d1_0
        double precision :: lkp_d2_0
        double precision :: lkp_beta_0
        double precision :: lkp_gamma_0
        double precision :: lkp_w_ref
        double precision :: lkp_b1_ref
        double precision :: lkp_b2_ref
        double precision :: lkp_b3_ref
        double precision :: lkp_b4_ref
        double precision :: lkp_c1_ref
        double precision :: lkp_c2_ref
        double precision :: lkp_c3_ref
        double precision :: lkp_c4_ref
        double precision :: lkp_d1_ref
        double precision :: lkp_d2_ref
        double precision :: lkp_beta_ref
        double precision :: lkp_gamma_ref
        double precision, dimension(30) :: vc_LKP

        !Sonstige Variablen (sp酹er durch in TREND verwendete Modulvariablen ersetzen)
        double precision, dimension(30, 30) :: kij_LKP, tc_ij_LKP, vc_ij_LKP

        !Sonstige Konstanten (sp酹er durch in TREND verwendete Konstanten ersetzen)
        double precision :: help_tau
        double precision :: help_delta
        double precision ::accenLKPMix, zcLKP, wLKP

        !end type module_LKP



        !type module_costald_eq

        !implicit none
        !save

        !!parameters for tait equation: saturated liquid density
        !    double precision,parameter:: ataitold=-1.52816d0
        !    double precision,parameter:: btaitold=1.43907d0
        !    double precision,parameter:: ctaitold=-0.81446d0
        !    double precision,parameter:: dtaitold=0.190454d0
        !    double precision,parameter:: etaitold=-0.296123d0
        !    double precision,parameter:: ftaitold=0.386914d0
        !    double precision,parameter:: gtaitold=-0.0427258d0
        !    double precision,parameter:: htaitold=-0.0480645d0
        !
        !    !parameters for tait equation: homogeneous density
        !    double precision,parameter:: atait=-9.070217d0
        !    double precision,parameter:: btait=62.45326d0
        !    double precision,parameter:: dtait=-135.1102d0
        !    double precision,parameter:: ftait=4.79594d0
        !    double precision,parameter:: gtait=0.250047d0
        !    double precision,parameter:: htait=1.14188d0
        !    double precision,parameter:: jtait=0.0861488d0
        !    double precision,parameter:: ktait=0.0344483d0
        !
        double precision,dimension(50):: omega, vstar

        !mixing parameters
        double precision :: omegaCOS, vstarCOS, TcCOS

        !parameters for RKM equation
        double precision::k1(0:10,0:10)
        double precision::k2(0:10,0:10)
        double precision::RKM_V(0:100,0:8)
        double precision::RKM_ps_me(0:100,0:1)
        double precision::RKM_x(0:8,0:9)
        double precision::RKM_M(1:8)
        double precision::RKM_Tc(1:8)
        character (12), dimension (30,3) :: RKM_fluids

        !end type

        !type module_flash

        !implicit none
        !save

        logical:: startvaluespin

        !end type

        !type module_PCSAFT

        !implicit none
        !save

        integer :: n_start, n_end
        double precision :: molfractions_save
        logical :: mol_save

        !number Pi
        !double precision, parameter :: piPCSAFT = 3.14159265359d0
        !double precision, parameter :: piPCSAFT = 3.14159265358979d0

        !mixing parameter for binary mixtures
        double precision, dimension (30,30) :: kij_PCSAFT

        !parameters for PC SAFT m, sigma, epsilon/k
        double precision, dimension (30):: mPCSAFT
        double precision, dimension (30):: sigPCSAFT
        double precision, dimension (30):: epskPCSAFT
        double precision, dimension (30):: kabPCSAFT
        double precision, dimension (30):: epsabkPCSAFT
        integer, dimension (30) :: n_sites
        character(1), dimension(30):: SAFTmodel_assoc

        !parameters for PCP SAFT n, Q
        double precision, dimension (30):: nPCSAFTQ
        double precision, dimension (30):: QPCSAFTQ
        !parameters for PCP SAFT n, My
        double precision, dimension (30):: nPCSAFTD
        double precision, dimension (30):: MyPCSAFTD

        !!universal model constants for equation 18 and 19, Table 1, Gross, Sadowski 2001
        !!a_0i, a_1i, a_2i, b_0i, b_1i, b_2i with i element of {0, 6}
        !double precision, dimension (0:6), parameter :: a0PCSAFT = (/ 0.9105631445d0, 0.6361281449d0, 2.6861347891d0, -26.547362491d0, 97.759208784d0, -159.59154087d0, 91.297774084d0 /)
        !double precision, dimension (0:6), parameter :: a1PCSAFT = (/ -0.3084016918d0, 0.1860531159d0, -2.5030047259d0, 21.419793629d0, -65.255885330d0, 83.318680481d0, -33.746922930d0 /)
        !double precision, dimension (0:6), parameter :: a2PCSAFT = (/ -0.0906148351d0, 0.4527842806d0, 0.5962700728d0, -1.7241829131d0, -4.1302112531d0, 13.776631870d0, -8.6728470368d0 /)
        !double precision, dimension (0:6), parameter :: b0PCSAFT = (/ 0.7240946941d0, 2.2382791861d0, -4.0025849485d0, -21.003576815d0, 26.855641363d0, 206.55133841d0, -355.60235612d0 /)
        !double precision, dimension (0:6), parameter :: b1PCSAFT = (/ -0.5755498075d0, 0.6995095521d0, 3.8925673390d0, -17.215471648d0, 192.67226447d0, -161.82646165d0, -165.20769346d0 /)
        !double precision, dimension (0:6), parameter :: b2PCSAFT = (/ 0.0976883116d0, -0.2557574982d0, -9.1558561530d0, 20.642075974d0, -38.804430052d0, 93.626774077d0, -29.666905585d0 /)
        !
        !!universal model constants for equation 15, 16 and 17, Table 2, Gross 2005
        !!a_0i, a_1i, a_2i, b_0i, b_1i, b_2i, c_0i, c_1i, c_2i with i element of {0, 4}
        !double precision, dimension (0:4), parameter :: a0PCSAFTQ = (/ 1.2378308d0, 2.4355031d0, 1.6330905d0, -1.6118152d0, 6.9771185d0 /)
        !double precision, dimension (0:4), parameter :: a1PCSAFTQ = (/ 1.2854109d0, -11.465615d0, 22.086893d0, 7.4691383d0, -17.197772d0 /)
        !double precision, dimension (0:4), parameter :: a2PCSAFTQ = (/ 1.7942954d0, 0.7695103d0, 7.2647923d0, 94.486699d0, -77.148458d0 /)
        !double precision, dimension (0:4), parameter :: b0PCSAFTQ = (/ 0.4542718d0, -4.5016264d0, 3.5858868d0, 0.d0, 0.d0 /)
        !double precision, dimension (0:4), parameter :: b1PCSAFTQ = (/ -0.8137340d0, 10.064030d0, -10.876631d0, 0.d0, 0.d0 /)
        !double precision, dimension (0:4), parameter :: b2PCSAFTQ = (/ 6.8682675d0, -5.1732238d0, -17.240207d0, 0.d0, 0.d0 /)
        !double precision, dimension (0:4), parameter :: c0PCSAFTQ = (/ -0.5000437d0, 6.5318692d0, -16.014780d0, 14.425970d0, 0.d0 /)
        !double precision, dimension (0:4), parameter :: c1PCSAFTQ = (/ 2.0002094d0, -6.7838658d0, 20.383246d0, -10.895984d0, 0.d0 /)
        !double precision, dimension (0:4), parameter :: c2PCSAFTQ = (/ 3.1358271d0, 7.2475888d0, 3.0759478d0, 0.d0, 0.d0 /)
        !
        !
        !!universal model constants for equation 12, 13 and 14, Table 1, Gross 2006
        !!a_0i, a_1i, a_2i, b_0i, b_1i, b_2i, c_0i, c_1i, c_2i with i element of {0, 4}
        !double precision, dimension (0:4), parameter :: a0PCSAFTD = (/ 0.3043504d0, -0.1358588d0, 1.4493329d0, 0.3556977d0, -2.0653308d0/)
        !double precision, dimension (0:4), parameter :: a1PCSAFTD = (/ 0.9534641d0, -1.8396383d0, 2.0131180d0, -7.3724958d0, 8.2374135d0 /)
        !double precision, dimension (0:4), parameter :: a2PCSAFTD = (/ -1.1610080d0, 4.5258607d0, 0.9751222d0, -12.281038d0, 5.9397575d0 /)
        !double precision, dimension (0:4), parameter :: b0PCSAFTD = (/ 0.2187939d0, -1.1896431d0, 1.1626889d0, 0.d0, 0.d0 /)
        !double precision, dimension (0:4), parameter :: b1PCSAFTD = (/ -0.58731640d0, 1.2489132d0, -0.5085280d0, 0.d0, 0.d0 /)
        !double precision, dimension (0:4), parameter :: b2PCSAFTD = (/ 3.4869576d0, -14.915974d0, 15.372022d0, 0.d0, 0.d0 /)
        !double precision, dimension (0:4), parameter :: c0PCSAFTD = (/ -0.0646774d0, 0.1975882d0, -0.8087562d0, 0.6902849d0, 0.d0 /)
        !double precision, dimension (0:4), parameter :: c1PCSAFTD = (/ -0.9520876d0, 2.9924258d0, -2.3802636d0, -0.2701261d0, 0.d0 /)
        !double precision, dimension (0:4), parameter :: c2PCSAFTD = (/ -0.6260979d0, 1.2924686d0, 1.6542783d0,-3.4396744d0, 0.d0 /)

        !----------------------------------------------------
        ! Avoiding double calculations:
        !----------------------------------------------------
        ! type variables for all PC SAFT function parts

        ! logical for determining if the arrays have already been allocated:
        logical :: already_allocated_PCSAFT

        ! !save the independent variables:
        double precision :: T_PCSAFT, D_PCSAFT
        double precision, dimension (30) :: x_PCSAFT ! note: this is hardcoded because molfractions is also hardcoded with 30
        integer :: nrsubst_PCSAFT

        ! saving which derivatives have been calculated for all function parts in the following type variables:
        ! the T, rho derivatives:
        integer, dimension (:), allocatable :: ar_calculated_PCSAFT
        ! the dxi (x1) derivatives:
        integer, dimension (:), allocatable :: arx1_calculated_PCSAFT
        ! the dxi dxj (x2) derivatives
        integer, dimension (:), allocatable :: arx2_calculated_PCSAFT
        ! the dxi dxj dxk (x3) derivatives
        integer, dimension (:), allocatable :: arx3_calculated_PCSAFT

        ! saving the results of the derivatives for all function parts in the following type variables
        ! mean segment diameter mmean
        double precision :: mmean_PCSAFT
        double precision, dimension (:,:), allocatable :: mmeanQij_PCSAFT
        double precision, dimension (:,:,:), allocatable :: mmeanQijk_PCSAFT
        double precision, dimension (:,:), allocatable :: mmeanDij_PCSAFT
        double precision, dimension (:,:,:), allocatable :: mmeanDijk_PCSAFT
        ! parameters a_j and b_j as ab
        double precision, dimension (2,0:6) :: ab_PCSAFT
        double precision, dimension (:,:,:,:), allocatable :: ab_PCSAFTQ
        double precision, dimension (:,:,:,:), allocatable :: c_PCSAFTQ
        !double precision, dimension (:,:,:,:,:), allocatable :: abx1_PCSAFTQ
        !double precision, dimension (:,:,:,:,:), allocatable :: cx1_PCSAFTQ
        double precision, dimension (:,:,:,:), allocatable :: ab_PCSAFTD
        double precision, dimension (:,:,:,:), allocatable :: c_PCSAFTD
        double precision, dimension (:,:,:), allocatable :: abx1_PCSAFT
        double precision, dimension (:,:,:,:), allocatable :: abx2_PCSAFT
        double precision, dimension (:,:,:,:,:), allocatable :: abx3_PCSAFT
        ! the T, rho derivatives:
        double precision, dimension (:), allocatable :: ar_PCSAFT, adisp_PCSAFT, c_PCSAFT, i1_PCSAFT, i2_PCSAFT, meo1_PCSAFT, meo2_PCSAFT, ahc_PCSAFT, ahs_PCSAFT, z0_PCSAFT, z1_PCSAFT, z2_PCSAFT, z3_PCSAFT
        double precision, dimension (:,:), allocatable :: gii_PCSAFT, di_PCSAFT, gij_PCSAFT
        double precision, dimension (:), allocatable :: A2_PCSAFTQ, A3_PCSAFTQ, AQQ_PCSAFTQ, A2_PCSAFTD, A3_PCSAFTD, ADD_PCSAFTD, AASSOC_PCSAFT
        double precision, dimension (:,:,:), allocatable :: J2_PCSAFTQ, J2_PCSAFTD
        double precision, dimension (:,:,:,:), allocatable :: J3_PCSAFTQ, J3_PCSAFTD
        double precision, dimension (:,:), allocatable :: xA_PCSAFT
        double precision, dimension (:,:,:,:), allocatable :: delta_AB ! fixed to binary associating mixtures and a maximum of 4 sites per molecule
        ! the dxi (x1) derivatives:
        double precision, dimension (:,:), allocatable :: arx1_PCSAFT, adispx1_PCSAFT, cx1_PCSAFT, i1x1_PCSAFT, i2x1_PCSAFT, meo1x1_PCSAFT, meo2x1_PCSAFT, ahcx1_PCSAFT, ahsx1_PCSAFT, z0x1_PCSAFT, z1x1_PCSAFT, z2x1_PCSAFT, z3x1_PCSAFT
        double precision, dimension (:,:,:), allocatable :: giix1_PCSAFT
        double precision, dimension (:,:), allocatable :: A2X1_PCSAFTQ, A3X1_PCSAFTQ, AQQX1_PCSAFTQ, A2X1_PCSAFTD, A3X1_PCSAFTD, ADDX1_PCSAFTD, AASSOCX1_PCSAFT
        double precision, dimension (:,:,:,:), allocatable :: J2X1_PCSAFTQ, J2X1_PCSAFTD
        double precision, dimension (:,:,:,:,:), allocatable :: J3X1_PCSAFTQ, J3X1_PCSAFTD
        ! the dxi dxj (x2) derivatives
        double precision, dimension (:,:,:), allocatable :: arx2_PCSAFT, adispx2_PCSAFT, cx2_PCSAFT, i1x2_PCSAFT, i2x2_PCSAFT, meo1x2_PCSAFT, meo2x2_PCSAFT, ahcx2_PCSAFT, ahsx2_PCSAFT
        double precision, dimension (:,:,:,:), allocatable :: giix2_PCSAFT
        double precision, dimension (:,:,:), allocatable :: A2X2_PCSAFTQ, A3X2_PCSAFTQ, AQQX2_PCSAFTQ, A2X2_PCSAFTD, A3X2_PCSAFTD, ADDX2_PCSAFTD, AASSOCX2_PCSAFT
        double precision, dimension (:,:,:,:,:), allocatable :: J2X2_PCSAFTQ, J2X2_PCSAFTD
        double precision, dimension (:,:,:,:,:,:), allocatable :: J3X2_PCSAFTQ, J3X2_PCSAFTD
        ! the dxi dxj dxk (x3) derivatives
        double precision, dimension (:,:,:,:), allocatable :: arx3_PCSAFT, adispx3_PCSAFT, cx3_PCSAFT, i1x3_PCSAFT, i2x3_PCSAFT, ahcx3_PCSAFT, ahsx3_PCSAFT
        double precision, dimension (:,:,:,:,:), allocatable :: giix3_PCSAFT
        !----------------------------------------------------

        !end type



        !type module_generalized_eos

        double precision:: polfac_eq
        !integer:: gen_i

        !EOS Alexandrov et al.
        double precision, dimension(14, 4) :: coefficientsi
        double precision, dimension(14, 3) :: exponentsi
        integer :: anz_term_igor

        !EOS Span
        double precision, dimension(10, 3) :: coefficientss
        double precision, dimension(10, 4) :: exponentss
        integer :: anz_term_span

        !EOS Sun & Ely
        double precision, dimension(14, 3) :: coefficientssun   !am-Werte
        double precision, dimension(14, 3) :: exponentssun
        integer :: anz_term_sun

        double precision :: acenfactor_1, polarfactor_1, acenfactor_2, polarfactor_2, acenfactor_3, polarfactor_3

        !end type


        !type module_literature_ref


        !implicit none
        !save

        type(type_lit_ref), allocatable :: litref

        !end type  module_literature_ref

        !type for ge model
        type(type_ge), allocatable:: ge

        type(type_sea), allocatable :: sea

        type(type_el), allocatable :: el

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !In this FORTRAN-code different gE-models are programed
        !The following models are available:
        ! -NRTL
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        !type module_NRTL
        !--------------------------------------------------------------------------------------
        ! Variables for the NRTL model
        ! In the original NRTL model, only three adjustable parameters are required per binary
        ! mixture. However, later modifications of the model introduced a temperature
        ! dependence of the parameters, which allow for a representation of experimental data
        ! over a larger temperature range.
        !--------------------------------------------------------------------------------------
        ! Andreas Jaeger, June 2016
        !Never used: commented
        ! also sub subroutine gE_NRTL(gl,Temp, gE, ln_gamma) is commented

        !double precision, dimension(30,30):: alpha_ij_0_NRTL
        !double precision, dimension(30,30):: alpha_ij_1_NRTL
        !
        !double precision, dimension(30,30):: tau_Aij_NRTL
        !double precision, dimension(30,30):: tau_Bij_NRTL
        !double precision, dimension(30,30):: tau_Cij_NRTL
        !double precision, dimension(30,30):: tau_Dij_NRTL
        !double precision, dimension(30,30):: tau_Eij_NRTL
        !double precision, dimension(30,30):: tau_Fij_NRTL
        !
        !double precision, dimension(30,30):: alpha_ij_NRTL
        !double precision, dimension(30,30):: tau_ij_NRTL
        !double precision, dimension(30,30):: G_ij_NRTL

        !end type module_NRTL




        !type module_UNIFAC
        !--------------------------------------------------------------------------------------
        ! Variables for the UNIFAC model
        !--------------------------------------------------------------------------------------
        ! Andreas Jaeger, January 2017

        double precision, dimension(:,:), allocatable :: R_ik     !Group volume parameters for group k of molecule i (i=30, k=100)
        double precision, dimension(:,:), allocatable :: Q_ik     !Group area parameters for group k of molecule i (i=30, k=100)
        double precision, dimension(:,:), allocatable :: v_ik     !Number of groups k in molecules of type i (i=30, k=100)
        double precision, dimension(:,:),allocatable :: a_nm    !Interaction parameter a for group n and m
        double precision, dimension(:,:),allocatable :: b_nm    !Interaction parameter b for group n and m
        double precision, dimension(:,:),allocatable :: c_nm    !Interaction parameter c for group n and m

        integer, dimension(30) :: nr_of_groups_i        !Number of different groups of which molecule i is composed of (example: propane 2x1 + 1x2 -> nr_of_groups(propane) = 2)
        integer, dimension(30) :: group_start_index     !This vector contains the indices where the groups of each component start, when all groups are written in one array (in maingroups_mix_list)
        integer :: nr_of_groups_mix                     !Total number of groups in the mixture (with double groups)
        integer, dimension(30,100) :: subgroups_ik_list     !List that contains the subgroups k (max 100) of every fluid i (max 30) in the mixture
        integer, dimension(3000) :: maingroups_mix_list     !List that contains all maingroups that appear in the mixture

        !Example for use of the variables subgroups_ik_list and maingroups_mix_list:
        !Mixture of ethane (subgroups: 2 x 1, maingroups: 2 x 1) + propane (subgroups: 2 x 1 + 1 x 2, maingroups: 3 x 1) + co2 (subgroup: 117, maingroup: 56)
        !subgroups_ik_list(1,1) = 1     v_ik(1,1) = 2
        !subgroups_ik_list(2,1) = 1     v_ik(2,1) = 2
        !subgroups_ik_list(2,2) = 2     v_ik(2,2) = 1
        !subgroups_ik_list(3,1) = 117   v_ik(3,1) = 1
        !maingroups_mix_list(1) = 1     (2 x 1)
        !maingroups_mix_list(2) = 1     (2 x 1)
        !maingroups_mix_list(3) = 1     (1)
        !maingroups_mix_list(5) = 56    (56)




        !end type module_UNIFAC

        !type with fit variables for Helm+gE model.
        !type module_fit_helmge

        !implicit none
        !save

        !!double precision, parameter:: resize_fac_abc = 1000.D0 !Factor to bring a,b,c for fitting to approximately the same order of magnitude
        !integer:: nr_of_fit_Ru      !Indicates how many group volumes shall be fitted
        !integer:: nr_of_fit_Qu      !Indicates how many group surfaces shall be fitted
        !integer:: nr_of_fit_auv     !Indicates how many group interaction parameters a shall be fitted
        !integer:: nr_of_fit_buv     !Indicates how many group interaction parameters b shall be fitted
        !integer:: nr_of_fit_cuv     !Indicates how many group interaction parameters c shall be fitted
        !!integer, parameter :: nr_max_params = 100   !Specifies the maximum number of adjustable parameters
        !integer, dimension(30,100):: nr_of_Ru_in_i !Indicates how often group Ru is present in molecule i
        !integer, dimension(30,100):: nr_of_Qu_in_i !Indicates how often group Qu is present in molecule i
        !integer, dimension(100) :: mul_fit_par !Indicates the multiplicity of the groups in the mixture.
        !!The parameters are given in the same order as the fitted params, see variable param_ind
        !!Only the fitted variables stand in this variable
        !!Example, see example below for the use of variable pos_fit_par
        !integer, dimension(100,100,2):: pos_fit_par     !Saves the position of the fit parameters as in the variables R_ik, Q_ik, a_nm, b_nm, c_nm
        !The data is stored in the variable as follows
        !First pos: The number i of the fitted parameter as indicated in the fit variable param_ind(i)
        !Second pos: If a group is present in the mixture more than once (see variable mul_fit_par), the different positions are indicated by the second entry
        !Third pos: For Q_u and R_U: 1: nr. of fluid i, 2: nr. of group k of fluid i
        !           For anm, bnm, cm: 1: row of interaction parameter matrix, 2: column of interaction parameter matrix
        !Example:
        !Binary mixture of CO2 and Ethane
        !CO2: subgroup 117, main group 56
        !Ethan: subgroup: 1 and 1, main group: 1 and 1
        !The variable fit_params (see fitter) holds all parameters of the model and the
        !variable mul_fit_par gives the multiplicity of the parameter. For the given example:
        !fit_params(1) = R_117      Rk(1,1)
        !fit_params(2) = R_1        Rk(2,1)
        !fit_params(3) = R_1        Rk(2,2)
        !fit_params(4) = Q_117      Qk(1,1)
        !fit_params(5) = Q_1        Rk(2,1)
        !fit_params(6) = Q_1        Rk(2,2)
        !fit_params(7) = a_56_56    a_nm(1,1)
        !fit_params(8) = a_56_1     a_nm(1,2)
        !fit_params(9) = a_56_1     a_nm(1,3)
        !fit_params(10) = a_1_56    a_nm(2,1)
        !fit_params(11) = a_1_1     a_nm(2,2)
        !fit_params(12) = a_1_1     a_nm(2,3)
        !fit_params(13) = a_1_56    a_nm(3,1)
        !fit_params(14) = a_1_1     a_nm(3,2)
        !fit_params(15) = a_1_1     a_nm(3,3)
        !fit_params(16) = b_56_56   b_nm(1,1)
        !fit_params(17) = b_56_1    b_nm(1,2)
        !fit_params(18) = b_56_1    b_nm(1,3)
        !fit_params(19) = b_1_56    b_nm(2,1)
        !fit_params(20) = b_1_1     b_nm(2,2)
        !fit_params(21) = b_1_1     b_nm(2,3)
        !fit_params(22) = b_1_56    b_nm(3,1)
        !fit_params(23) = b_1_1     b_nm(3,2)
        !fit_params(24) = b_1_1     b_nm(3,3)
        !fit_params(25) = c_56_56   c_nm(1,1)
        !fit_params(26) = c_56_1    c_nm(1,2)
        !fit_params(27) = c_56_1    c_nm(1,3)
        !fit_params(28) = c_1_56    c_nm(2,1)
        !fit_params(29) = c_1_1     c_nm(2,2)
        !fit_params(30) = c_1_1     c_nm(2,3)
        !fit_params(31) = c_1_56    c_nm(3,1)
        !fit_params(32) = c_1_1     c_nm(3,2)
        !fit_params(33) = c_1_1     c_nm(3,3)

        !Assume that we want to fit the parameters: R_1, Q_117, a_56_1, and b_1_56
        !Then:
        !param_ind(1) = 2       R1     -> First fit variable
        !param_ind(2) = 4       Q_117  -> Second fit variable
        !param_ind(3) = 7       a_56_1 -> Third fit variable
        !param_ind(4) = 16      b_1_56 -> Fourth fit variable
        !The position variable would finally be filled as follows:
        !pos_fit_par(1,1,1) = 2     !First variable, first position, fluidnr ->       Rk(2<-, 1)
        !pos_fit_par(1,1,2) = 1     !First variable, first position, groupnr ->       Rk(2  , 1<-)
        !pos_fit_par(1,2,1) = 2     !First variable, second position, fluidnr ->      Rk(2<-, 2)
        !pos_fit_par(1,2,2) = 2     !First variable, second position, groupnr ->      Rk(2  , 2<-)
        !pos_fit_par(2,1,1) = 1     !Second variable, first position, fluidnr ->      Qk(1<-, 1)
        !pos_fit_par(2,1,2) = 1     !Second variable, first position, groupnr ->      Qk(1  , 1<-)
        !pos_fit_par(3,1,1) = 1     !Third variable, first position, row nr. ->      anm(1<-, 2)
        !pos_fit_par(3,1,2) = 2     !Third variable, first position, column nr. ->   anm(1  , 2<-)
        !pos_fit_par(3,2,1) = 1     !Third variable, second position, row nr. ->     anm(1<-, 3)
        !pos_fit_par(3,2,2) = 3     !Third variable, second position, column nr. ->  anm(1  , 3<-)
        !pos_fit_par(4,1,1) = 2     !Fourth variable, first position, row nr. ->     bnm(2<-, 1)
        !pos_fit_par(4,1,2) = 1     !Fourth variable, first position, column nr. ->  bnm(2  , 1<-)
        !pos_fit_par(4,2,1) = 3     !Fourth variable, second position, row nr. ->    bnm(3<-, 1)
        !pos_fit_par(4,2,2) = 1     !Fourth variable, second position, column nr. -> bnm(3  , 1<-)

        !store the information of the multiplicity of the fit variables in mul_fit_par
        !mul_fit_par(1) = 2
        !mul_fit_par(2) = 1
        !mul_fit_par(3) = 2
        !mul_fit_par(4) = 2

        !Check how often FITTED Qu and Ru are present in molecule i
        !For the given example:
        !nr_of_Ru_in_i(i,u)     -> 2 substances, 2 fitted parameters (among Ru+Qu): 1 Ru, 1 Qu
        !nr_of_Ru_in_i(1,1) = 0     no group R_1 in molecule 1
        !nr_of_Ru_in_i(2,1) = 2     two groups R_1 in molecule 2
        !nr_of_Qu_in_i(1,2) = 1     one group Q_117 in molecule 1
        !nr_of_Qu_in_i(2,2) = 0     no group Q_117 in molecule 1



        !end type




        !type hdrt_property_definition
        !******************************************************************************
        !
        !  Module with property constants
        ! October 2012; modified May-October 2015
        ! November 2016:
        ! universal reference state conditions T0 = T0a = 273.15 K, p0 = p0a = 1 Pa
        ! new lattice parameter correlation - conference article EFM 2016
        ! REOS - Reference EoS (not GERG)
        ! BS - Langmuir from - multilayered water shell by Ballard & Sloan
        ! M4 - compressibility from Murnaghan EoS
        ! CONS - constant shell radii for all hydrate formers
        ! constant gw0 and hw0 for both sI and sII hydrates

        !/////////////////////////////////////////////////////////////////////////

        logical :: numderiv




        ! Equations of state for fluid phases:
        character(4) :: EoS_type                      ! 'reos' = reference EoSs + EoS-CG mix rules or GERG mix rules if EoS-CG unavailable
        ! 'gerg' = GERG 2004 with GERG orig. reference state
        ! 'srk' = SRK with quadratic mixing rules for a and linear mixing rules for b used
        !         If available, binary interaction parameter kij used
        ! Shell radius def.
        character*4 :: Ra_def                        ! 'r(a)' = T,p dependent shell radii Ra = R0*a0/a(T,p)
        ! 'tref' = constant shell radii Ra = R0(T_ref) * a(T_ref)/a_ref
        ! 'cons' = constant shell radii Ra, Ra = R_S, R_L ... used since 10.2016

        !character(2) :: hdrt_structure
        character(2) :: Langmuir                      ! 'jh' = John & Holder (1982)    ... three water shells
        ! 'ks' = Klauda & Sandler (2002) ... T-fit
        ! 'bs' = Ballard & Sloan (2002)  ... multilayered cages
        ! 'vj' = Vins & Jaeger :-) ... multilayered cages + two additional layers
        ! 'tr' = Trout's group from MIT ... "Experimental" Langmuir constants
        !                                    --> fit two potential parameters m and r0 to potential derived from Langmuir constant for small and large cavity
        !                                    Parameters available for: Methane, Ethane, Propane, Argon

        character(2) :: Latt_par                      ! 'm4' = correlation based on Murnaghan EoS with dB0 = 4 = const.
        ! 'bs' = Ballard & Sloan = constant compressibility (kapa1)
        ! 'vj' = Vins & Jaeger = pressure dependent compress. (kapa1,2)

        ! Reference conditions @ T0, p0
        double precision :: gw_B0                               ! Reference state Gibbs energy of water in empty lattice [J/mol]
        double precision :: hw_B0                               ! Reference state Enthalpy of empty lattice [J.mol^-1]
        double precision :: gw_B0_a0, hw_B0_a0
        double precision, dimension(2) :: kb_gw                 ! slope and offset of the linear correlation for gw_B0
        double precision, dimension(2) :: kb_hw                 ! slope and offset of the linear correlation for hw_B0

        double precision, dimension(3) :: alfa ! thermal expansion parameters
        double precision :: T0_a, p0_a ! T0_a: [K] reference temperature for a [A]
        ! p0_a: [Pa] reference pressure for a [A] ... p0_a = 1 Pa for 'm4' since October 2016
        double precision, dimension(30) :: a0_hdrt, a_hc, sigma, eps_k, sigmad, epsd_k       !sigmad and epsd_k set "0" for all single occupand hydrate formers! physically correct: Kihara parameters "JJ"
        double precision, dimension(30) :: r_vdw    !half Moleculecular distance (vdW Radii) for double occupancy, if =0, then not known for double occ., fitting parameter?!
        double precision :: a_ref                   ! [A] reference lattice parameter @ T_ref
        double precision :: T_ref                   ! [K] reference temperature for the shell radii
        double precision :: p_ref                   ! [Pa] reference pressure for the shell radii
        double precision :: l_core
        !double precision, dimension(6) :: R_S, R_L, R_M
        !integer, dimension(6) :: z_S, z_L, z_M
        double precision, dimension(3,6) :: R_shell
        integer, dimension(3,6) :: z_shell

        integer :: N_hdrts                              !The number of how many components in the mixture form hydrates (with water)
        integer :: N_guests                             !The number of guest components forming hydrates (WITHOUT water)
        double precision, dimension(30) :: moles_hdrt, Y_yj_tf   !Sorted molfractions equal to moleslist_hydrate from module_solids
        double precision, dimension(2,30) :: B0_hdrt         !Bulk modulus of pure hydrate formers
        double precision, dimension(2) :: B0_m          !Bulk modulus of mixed hydrate
        ! [Pa, -] Bulk modulus and its derivative from Murnaghan EoS - mixed hydrate
        double precision :: a0_m                        !lattice parameter of mixed hydrate at T0_a & p0_a = 0 Pa
        double precision :: a_T_ref                     !lattice parameter of mixed hydrate at T_ref & p_ref
        ! [A] lattice parameter of the gas at T_ref and low p of 1 Mpa

        integer :: occupmax

        double precision, dimension(3,30) :: KS_s, KS_l         ! Klauda&Sandler T-fit for the Langmuir constant - only sI & sII

        double precision, dimension(3) :: kapa
        double precision, dimension(2,30) :: m_const, rs_const
        character(6) :: Latt_par_mix                  !mixing rule for lattice parameter of mixed hydrates
        !Latt_par_mix = 'volume'/   ! 'volume' = ideal mixing rule for volume of unit cell
        !Latt_par_mix = 'vegard'/   ! 'vegard' = Vegard's mixing rule for lattice parameter
        character(3) :: mixdecision                     !mixing rule for Langmuir Constants of double occupied cages
        !///////////////////////////////////////////////////////////////////////////////////////////////////


        ! -------------------------
        ! Structure of the hydrate:
        ! -------------------------
        ! data hdrt_structure/'s1'/  ! 's1' = hydrate structure sI
        ! 's2' = hydrate structure sII
        ! 'sH' = hydrate structure sH (not programmed yet)
        ! value for 'hdrt_structure' needs to be defined outside of the 'hdrt_property_definition' module
        ! ... needed in searching for stable hydrate structure at given conditions

        ! --------------------------------------------------------
        ! Molar heat capacity of Hydrate (initial values - ice Ih)
        ! --------------------------------------------------------
        ! heat capacity of the empty hydrate lattice can be well approximated
        ! by the cubic, respectively hexagonal, ice Ih ... fit from IAPWS ice Ih EoS
        ! for p0 = 2.0 MPa ... T = 140 - 300 K
        !data cpH_a/0.12814d0/     ! [J.mol^-1] constant in cp(T): cp = cp_a*T + cp_b
        !data cpH_b/2.74566d0/     ! [J.mol^-1.K^-1] constant in cp(T) expression
        ! for p0 = 0 MPa ... T = 120 - 300 K
        double precision ::  cpH_a     ! [J.mol^-1] constant in cp(T): cp = cp_a*T + cp_b
        double precision ::  cpH_b     ! [J.mol^-1.K^-1] constant in cp(T) expression
        ! linear correlation for cp(ice Ih) stays the SAME for p0 = 1 Pa!!! - checked November 2016

        ! -----------------------------
        ! Reference conditions @ T0, p0
        ! -----------------------------
        double precision ::  T0        ! [K] reference temperature for Gibbs energy
        !data p0/2.0d6/           ! [Pa] reference pressure for Gibbs energy
        !........................................................................
        !THE REFERENCE POINT OF THE EMPTY BETA LATTICE WAS CHANGED TO p0 = 0 PA
        !THE REASON FOR THIS IS, THAT AT 0 PRESSURE ALL HYDRATES OF ONE TYPE MAY
        !HAVE COMPARABLE REFERENCE CONDITIONS, SINCE THE PRESSURE TERM VANISHES,
        !THE THERMAL EXPANSION TERM IS THE SAME FOR ALL HYDRATES OF THE SAME
        !STRUCTURE. THUS ONLY a0 SHALL MAKE A DIFFERENCE IN THE REFERENCE
        !PARAMETERS. Andreas October 2013
        !data p0/0.D0/           ! [Pa] reference pressure for Gibbs energy
        !........................................................................
        !The reference pressure both for the lattice parameter a(T0_a, p0_a) and
        !the empty beta-lattice was changed to 1 Pa ... October 2016
        !Reason: calculation of entropy of an ideal gas at the ref. pressure 1 Pa
        double precision ::  p0           ! [Pa] reference pressure for Gibbs energy = 1 Pa


        !      INITIATION of other variables

        !       |           |            |
        !       |           |            |
        !       |           |            |
        !       V           V            V

        ! -----------------------------------
        ! # of H2O molecules in one unit cell
        ! -----------------------------------
        double precision :: Nw               ! [-] number of water molecules in unit cell of given hydrate structure

        ! ------------------------------
        ! # of cavities per H2O molecule
        ! ------------------------------
        double precision, dimension(3) :: v_cavi                     ! first small then large cavity (third is for medium cavity in sH structure)
        integer :: N_cavi                                           !Number of cavities in the given hydrate structure sI & sII = 2, sH = 3
        ! number of cavities ... 2 for sI and sII (3 for sH)

        ! --------------------
        ! Lattice parameter a:
        ! --------------------
        !data a0/0.d0/               ! [A] lattice param @ T0, p0
        !data kapa/0.d0, 0.d0, 0.d0/ ! isothermal compressibilty - coefs.

        ! ------------------------------------------------------------
        ! Kihara potential parameters & Reference conditions @ T0, p0:
        ! ------------------------------------------------------------
        !data a_hc/0.d0/     ! [A] hard core radius (considered as constant)
        !data sigma/0.d0/    ! [A] collision diameter
        !data eps_k/0.d0/    ! [K] energy parameter - potential depth


        !!New definition if the reference conditions
        !!The reference condition according to Ballard and Sloan are devided by a0, such that
        !!they are universal for each hydrate former of one structure
        !data gw_B0_a0/0.d0/        ! Gibbs energy of water in empty lattice [J/mol], devided by the lattice parameter [A] @T0 and p0
        !data hw_B0_a0/0.d0/        ! Enthalpy of empty lattice [J.mol^-1], devided by the lattice parameter [A] @T0 and p0

        ! -------------------------------------------
        ! Water shell radii and coordination numbers:
        ! -------------------------------------------
        !data R_S/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/      ! small cavity
        !data z_S/0, 0, 0, 0, 0, 0/
        !data R_L/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/ ! large cavity
        !data z_L/0, 0, 0, 0, 0, 0/
        !data R_M/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/ ! medium cavity - only for sH hydrate
        !data z_M/0, 0, 0, 0, 0, 0/

        ! ----------------------------------
        ! Klauda and Sandler T-fit for C_iJ:
        ! ----------------------------------
        !data KS_s/0.d0, 0.d0, 0.d0/
        !data KS_l/0.d0, 0.d0, 0.d0/

        ! --------------------------------
        ! Trout's group from MIT for C_iJ:
        ! --------------------------------
        !data m_const/0.D0, 0.D0/     !Small & large cavity
        !data rs_const/0.D0, 0.D0/    !Small & large cavity

        !end type hdrt_property_definition


        !type module_DryIce
        !--------------------------------------------------------------------------------
        ! In this module variables and constants are stored that are needed for the
        ! fundamental equation of dry ice
        !--------------------------------------------------------------------------------
        ! A. Jaeger, Dec. 2011

        !Constants g0 - g10, n
        double precision :: g0, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10
        double precision :: exp_n
        !constants g0_alpha to g8_alpha
        double precision :: g0_a, g1_a, g2_a, g3_a, g4_a, g5_a, g6_a, g7_a, g8_a
        !constants g0_kappa to g2_kappa
        double precision :: g0_k, g1_k, g2_k
        !Reference point
        double precision :: T0_dryice, p0_dryice, R_CO2
        !constant pi
        double precision :: const_pi

        !end type module_DryIce


        !type module_WaterIce
        !--------------------------------------------------------------------------------
        ! In this module variables and constants are stored that are needed for the
        ! fundamental equation of dry ice
        !--------------------------------------------------------------------------------
        ! A. Jaeger, Feb. 2011

        !Constants g0, g00 - g04
        double precision :: g0_waterice, g00_waterice, g01_waterice, g02_waterice, g03_waterice, g04_waterice
        !constant s0 (for linking the solid equation to the IAPWS-95)
        double precision :: s0
        !complex number constants
        double complex :: t1, t2, r1, r2, r20, r21, r22
        !Reference point
        double precision :: Ttr_water, ptr_water, p0_water
        !Molar mass
        double precision :: M_H20


        !end type module_WaterIce

        !type module_seawater
        !--------------------------------------------------------------------------------
        ! In this module constants are stored that are needed for the
        ! gibbs equation of Seawater
        !--------------------------------------------------------------------------------
        !  May. 2018

        !these locialcs and integers are needed in flash routines --> not moved to seawater type
        logical:: seawater, sal, seacalc,  seawatercalled, seawaterflash!, seaphase
        integer :: salpos, seapos

        !general variables for ge-carrier in props_derivs_mix and further model handling that can not be stored sensfully in the splitted up types
        !Benedikt 09/2019
        double precision :: gepress
        integer :: modelflag
        !logical :: ge_carrier
        logical :: el_present
        logical :: gecarrier

        character(255):: path



        !end type module_seawater




        !type module_ecs


        !reference fluid variables
        character(12) :: input_0
        double precision :: t_0
        double precision :: d_0
        character(255) :: fluids_0
        character(255) :: moles_0
        double precision, dimension(30) :: molev_0
        character(255) :: eos_indicator_0
        character(255) :: eos_indicator_0_call
        double precision :: tc_0
        double precision :: dc_0
        double precision :: m_0
        double precision :: sigma_0
        double precision :: epsilon_0

        !fluid variables
        double precision, dimension(30):: tc_fluid
        double precision, dimension(30) :: dc_fluid
        double precision, dimension(30) :: w_fluid
        double precision, dimension(30) :: m_fluid
        integer :: ncomp_fluid
        double precision :: alpha_temp
        double precision :: dalpha_ddelta_temp
        double precision :: d2alpha_ddelta2_temp
        double precision :: dalpha_dtau_temp
        double precision :: d2alpha_ddelta_dtau_temp
        double precision :: z_temp
        double precision :: residual_ecs

        !end type module_ecs

        !type module_bell_jaeger
        Double precision, allocatable :: n2d2lnfidnjdnk(:,:,:)
        double precision, allocatable :: ndndlnfidnjdnk(:,:,:)
        double precision, allocatable :: dndndlnfidnjdnkddel(:,:,:)
        double precision, allocatable :: dn2d2lnfidnjdnkddel(:,:,:)
        double precision, allocatable :: dn2d2lnfidnjdnkdtau(:,:,:)
        double precision, allocatable :: dndndlnfidnjdnkdtau(:,:,:)
        double precision, allocatable :: d_n_dlnfi_dnj_dxm(:,:,:)
        double precision, allocatable :: d2_n_dlnfi_dnj_dxmdtau(:,:,:)
        double precision, allocatable :: d2_n_dlnfi_dnj_dxmddel(:,:,:)
        Double precision, allocatable :: dndnardnidnjALL(:,:,:)
        Double precision, allocatable :: dndlnfidnjdxk(:,:,:)
        Double precision, allocatable :: ndardnidxi_all(:,:,:)
        Double precision, allocatable :: dndndardnidnjdxk_ALL(:,:,:,:)
        Double precision, allocatable :: d2ndlnfidnjdxkddel(:,:,:)
        Double precision, allocatable :: d2ndlnfidnjdxkdtau(:,:,:)
        double precision, allocatable :: d2PSI_Y_dxjdxk(:,:,:,:)
        double precision, allocatable :: MIXDERIVFNR_dxidxj(:,:,:)
        double precision, allocatable :: dPSI_Y_dxj(:,:,:)
        double precision, allocatable :: MIXDERIVFNR_dxidxjdxk(:,:,:,:)
        double precision, allocatable :: PSI_Y(:,:)
        Double precision, allocatable :: ndndardnidnj_ALL(:,:,:)
        Double precision, allocatable :: dndardnidxi_ALL(:,:,:)
        Double precision, allocatable :: d2ndardnidxidxj_ALL(:,:,:,:)

        integer :: converttype, unitin, unitstat, unitout, already_converted                      !Wird 鄣nlich dem Calctype_internal f�r die Umrechnung von molaren und spezifischen Gr廲en genutzt. Fallunterscheidungen nach diesem Typen!
        logical :: vir
        character(20) :: unitdefinition
        !end type module_bell_jaeger


        !!Erik, April 2018
        !!--------------------------------------------------------------------------------------
        !! Variables for the COSMO-SAC model
        !!--------------------------------------------------------------------------------------
        !
        !!integer :: interval !, version
        !double precision, dimension(30) :: Acosmo, Vcosmo
        !double precision, dimension(51) :: counter
        !double precision, dimension(51,30) :: sigma
        !!Variables COSMO-SAC Version 2 and 3
        !!double precision, dimension(30,3) :: Acosmo_ver23
        !double precision, dimension(51,30,3) :: sigma_v23
        !double precision, dimension(51,3) :: counter_v23
        !integer :: COSMO_ver                !Choose COSMO-SAC version with this variable:
        !                                        !1 COSMO-SAC Version by Lin, Sandler 2002;
        !                                        !2 COSMO-SAC improvements by Hsieh, Sandler, Lin 2010;
        !                                        !3 COSMO-SAC considering dispersive interactions by Hsieh, Lin, Vrabec 2014
        !                                        !4 COSMO-SAC Xiong et al 2014
        !double precision :: temp_cosmo_prev
        !double precision, dimension(30) :: molefraction_cosmo_prev
        !double precision, dimension(51,2,3) :: seggamma_pure_prev   !nach test 2 auf 30 鄚dern
        !double precision, dimension(51,3) :: seggamma_prev
        !double precision, dimension(15,30) :: ln_gamma_C_cosmo_prev, ln_gamma_R_cosmo_prev
        !double precision, dimension(15) :: gE_C_cosmo_prev, gE_R_cosmo_prev
        !integer, dimension(30) :: molecule_type    !for dispersive interactions
        !double precision, dimension(30) :: eps_molecule    !for dispersive interactions
        !double precision, dimension(30,17) :: m_tau
        !double precision, dimension(30) :: NBT_cosmo    !Normal boiling temperature to calculate molar Volume in Version 4 (Xiong)
        !double precision, dimension(30) :: ln_gamma_dsp
        !
        type(type_COSMO_SAC), allocatable:: cosmo

        type(cross),allocatable:: crs
        logical:: is_crs
        logical:: virial_num

        logical :: trend_calc_called    ! info if trend_calc was already called
        

    end type type_gl




    !interface
    !    !subroutine uppertolower_char(string,n)
    !    !    integer :: n
    !    !    character*n :: string                    ! fluid vector
    !    !end subroutine uppertolower_char
    !
    !    !subroutine  Sorter(x, dim, length)
    !    !    integer :: dim, length
    !    !    character*length, dimension(dim) :: x
    !    !end subroutine Sorter
    !end interface

    contains



    end module module_all_types


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !-------------------------------------------------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------------------------------------------------
    !S. Pohl
    !Module for the data type for fitting equation of state
    module main_Data



    type fitter_data_main


        double precision, allocatable:: pp(:), dd(:), tt(:),  pprop(:), dum(:), properr(:), ww(:), tau(:), delta(:)
        character*255, allocatable:: author(:) , proptype(:)
        integer:: ndata,n_pot,ndata_a_crit,ndata_norm
        double precision, allocatable:: t_red(:), r_red(:), t_trip(:), r_const(:), mol_w(:), fluid_param(:,:)
        integer, allocatable:: npoint(:), authorid(:),n_ptypes(:,:),pos_ptypes(:,:), pot_id(:), next_pot(:)
        integer, allocatable:: n_ptypes_crit(:,:),pos_ptypes_crit(:,:),next_pot_crit(:) !positions and numbers of critical data
        logical, allocatable:: pos_mem(:,:),pos_mem_crit(:,:)
        !data save area for the data for ancillary equation fitting
        double precision, allocatable, dimension(:,:,:):: ancillary_fit_data !(datap,4,npot)


    end type



    !this is a copy of the parameters of each equation saved on the heap
    double precision, allocatable:: eqn_parameter_mem(:,:,:)
    integer, allocatable:: size_vector_mem(:,:)
    integer:: max_parameter


    !for give the data to trend


    integer, allocatable:: identify(:,:)



    !ancillary equations
    double precision, allocatable:: vpcoeff(:),vpexp(:), dlcoeff(:), dlexp(:), dvcoeff(:), dvexp(:)
    integer, dimension(6):: anc_terms !(1) nvpcoeff, (2) ndlcoeff, (3) ndvcoeff


    !bundle to give it to trend
    double precision, allocatable:: anc_param(:,:) !vpcoeff,vexp,dlcoeff,dlexp,dvcoeff,dvexp

    end module

    module constraints

    !type for one constraint
    ! Start for prop and end of prop e.q. t_min, t_max
    ! Values calculated for the prop
    ! Violation of the constraint
    type results
        double precision , allocatable:: res_col(:)
    end type

    type constraint_individ
        double precision, allocatable:: vals(:),violation(:)

        type(results),allocatable:: res_row(:)
    end type

    type constraints_all_mems
        type(constraint_individ), allocatable:: calcs(:)
        double precision , allocatable:: jac_const(:,:)
        double precision ,allocatable:: const_dev(:)
    end type

    type grr
        double precision , allocatable:: grid_1(:)
        double precision , allocatable:: grid_2(:)
    end type

    type constraint
        double precision,allocatable :: start_1(:),ende_1(:),compare_val(:),start_2(:),ende_2(:),weight(:)
        integer, allocatable:: nrs_1(:),nrs_2(:)
        type(grr),allocatable:: grids(:)
        character*12,allocatable:: calctype(:),compare_sign(:),input(:)
        integer:: nr_cst_mem
        type(constraints_all_mems),allocatable:: mem(:)
    end type


    end module

    module module_matrices



    type levenberg


        double precision, allocatable:: r_vec(:),r_vec_for_rho(:)
        double precision, allocatable:: jac(:,:)
        double precision, allocatable:: constr_jac(:)
        double precision, allocatable::  wc(:)

        double precision, allocatable:: N1(:,:)
        double precision, allocatable:: N2(:,:)
        double precision, allocatable:: N3(:)
        double precision, allocatable:: IM(:,:)
        double precision, allocatable:: LPD1(:)
        double precision, allocatable:: RPD1(:)
        double precision, allocatable:: RPD1_S(:)
        double precision, allocatable:: LPD2(:,:)
        double precision:: maxN
        double precision:: product_rpd1
        double precision:: alpha
        double precision:: decide

        !controlling damping of inverse
        double precision:: damp_inverse
        logical:: first_r_mainpulation
        integer:: count_same
        double precision:: SSQ_LAST


        !shows if critical data is added(1) to the fit or not(0)
        integer:: fit_crit,fit_crit_save
        double precision:: ssq_Dev


    end type





    end module module_matrices

    module module_eqn_props

    use, intrinsic :: iso_c_binding
    !use crs_var
    !variables for the lagrange constrain function
    type general_constrains

        double precision:: penalty_func
        double precision, allocatable:: slack_var(:,:,:)
        double precision, allocatable:: mue(:,:,:,:) ,psi(:,:,:,:)
        integer:: nr_constrains !same for all potentials
        integer, allocatable :: nr_temp(:)
        integer, allocatable:: nr_pts(:)
        double precision, allocatable:: temp(:,:,:),dense(:,:,:),press(:,:,:,:),targ(:,:,:,:),tau(:,:,:),delta(:,:,:)
        double precision:: dense_perc
        integer, allocatable:: ctype(:)
        logical:: constr_on !if .true. constraints are calculated, else they were skipped


    end type

    !ancillary equation for each equation
    type anc

        integer, allocatable:: nvpcoeff(:)
        double precision, allocatable:: vpcoeff(:,:),vpexp(:,:)
        integer, allocatable:: ndlcoeff(:)
        double precision, allocatable:: dlcoeff(:,:),dlexp(:,:)
        integer, allocatable:: ndvcoeff(:)
        double precision, allocatable:: dvcoeff(:,:),dvexp(:,:)


    end type

    type log_constraints
        double precision:: mue_fac
    end type



    type props

        type(c_ptr):: handle
        double precision, dimension(50,11):: param, param_old, param_init, param_mem !parameter matrix
        double precision , dimension(:), allocatable:: crs_param,crs_param_old

        !arrays for the general coefficient function
        double precision, allocatable:: general_coeffs(:,:),general_coeffs_old(:,:),general_coeffs_mem(:,:)  !memory for general coefficients(by now maximal 5)
        double precision, allocatable:: general_coeffs_result(:,:,:)    !for calculations in Trend (term_nr,datapoint_nr)
        double precision, allocatable:: general_coeff_factors(:,:)

        integer:: nr_gen_coeffs

        !arrays for the general temperature exponent function
        double precision, allocatable:: general_texp(:,:),general_texp_old(:,:),general_texp_mem(:,:)
        double precision, allocatable:: general_texp_result(:,:,:)
        double precision, allocatable:: general_texp_factors(:,:)

        integer:: nr_gen_texp

        !index of the best general solution (needed for output later)
        integer:: best_general_idx

        !matrix with factores of texp function
        double precision, allocatable:: general_texp_start(:,:)
        double precision, allocatable:: work_array(:)



        integer:: I_P, I_E, I_G, I_L, I_L_D, n_term, n_reg,ndata,n_pot
        integer:: I_K
        double precision:: ssq,ssq_lm,ssq_old,ssq_mem
        double precision, allocatable:: delta_a(:)

        integer, allocatable:: fit_indicator(:)
        integer:: act_fit_comb

        double precision, allocatable:: t_red(:), r_red(:), p_red(:)
        double precision, allocatable:: t_red_old(:), r_red_old(:), p_red_old(:)


        double precision, allocatable, dimension(:,:,:):: red_num_derivs

        integer, allocatable:: pot_id(:)
        !for solving equation system
        integer, allocatable:: ipiv(:)
        integer:: info
        integer:: main_fails

        integer:: fail_step
        integer:: fail_step_count
        !type(anc):: anc_eqn
        type(general_constrains)::  constr
        type(log_constraints):: log_constr

        logical:: kernel_fit
        logical:: kernel_active
        integer:: nr_crs_param
        !type(general_term) :: gterm

    end type

    end module module_eqn_props

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !-------------------------------------------------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------------------------------------------------

    module module_hdrt_phaselines
    use module_all_types
    use module_hdrt_parameters

    !variables for phase equilibrium lines
    type type_hdrt_4ph_lines
        !>Variables for Quadruple Lines
        !>Calculate values for four phase line only for physical reasonable values of phasefraction (0 <= beta_ph <= 1)
        logical :: fourphase_physical, show_progress
        !!
        !!Q_map:
        !!in rows start points of quadruple line "types" are saved:
        !!VLHI, VLLH, etc.
        !!in columns start points for a three phase line from a quadruple line are saved
        !!1: VLwHIw: VLH, LHI
        !!2: VLwLcH: VLL, LwLcH
        !!3: not defined for mixtures now
        !!4: not defined for mixtures now
        integer, dimension(6) :: Q_map
        !>Type of the quadruple point. Dimensions:
        !!(1: steps between pure components 1 is pure guest#1, np_4ph is pure guest#2)
        character(6), dimension(np_4ph) :: Q_point_mix
        !>pressure and Temperature of the corresponding quadruple line
        !!(1: steps between pure components 1 is pure guest#1, np_4ph is pure guest#2)
        double precision, dimension(np_4ph) :: press_Q_mix, Temp_Q_mix
        !>Phasefractions are stored in this variable. Dimensions are:
        !!(1: phase#, 2: steps between pure components 1 is pure guest#1, np_4ph is pure guest#2)
        double precision, dimension(4,np_4ph) :: phasefrac_mix, rho_Q_mix
        !>Positions of emerging phases at roots are saved in this variable
        !!(1: root#, 2: emerging phase#)
        integer, dimension(6,2) :: phaseemerge_mix
        !>Compositions of each phase and each Quadruple Point are saved in these variables. (:,1) belongs to the Quadruple Point of pure hydrate: water - guest#1.
        !!(:,100) belongs to the Quadruple Point of pure hydrate: water - guest#2.
        !!Dimensions are: (1: molfractions of components, 2: steps between pure components 1 is pure guest#1, 100 ist pure guest#2)
        double precision, dimension(30,np_4ph) :: x_Q_ph1_mix, x_Q_ph2_mix, x_Q_ph3_mix, x_Q_ph4_mix
        logical, dimension (np_4ph) :: Qexist_mix

    end type type_hdrt_4ph_lines

    type type_hdrt_3ph_lines

        !Variables for threephase lines
        integer :: np_3ph, nflashs, loops, EqTypeStart
        character(20) :: trline_mix, strFlash
        character(255) :: pathout
        character(6) :: EqType
        double precision, dimension(:), allocatable :: press_tr_mix, temp_tr_mix
        double precision, dimension(:,:), allocatable :: chempot, lnfi
        double precision, dimension(:,:,:), allocatable :: x_ph
        double precision, dimension(:,:), allocatable :: rho_tr_mix
        double precision, dimension(:,:,:), allocatable :: occup
        double precision, dimension(:,:), allocatable :: hyd_nr
        double precision, dimension(:), allocatable :: s_hyd, h_hyd, h_melt, tpd_sol, a0_other_struct

    end type type_hdrt_3ph_lines

    type type_hdrt_pheq
        !>toggle for turning head row in outputdata files on or off
        logical :: print_head, show_progress

        !>Variables for Calculation of (Hydrate) phase equilibirum lines
        integer :: CompositionType !1: Overall Comp   2: Composition of the VaporPhase at Equilibrium (water-free base)

        ! Public (global) arguments
        character(len=255) :: CompsHdrtsAll   !Guests forming hydrates - variable for the naming of files


        !> n_3ph_lines: number of three phase lines without different betas = 0
        !> cur_3ph: three phase line that is calculated at the moment
        integer :: n_3ph_lines
        ! subtype for three phase lines
        ! each dimension is for one three phase equilibirum line with beta pi = 0
        type(type_hdrt_3ph_lines), dimension(:), allocatable :: phl3

        ! subtype for four phase lines
        !4 dimensional for each phase equilibirum type
        type(type_hdrt_4ph_lines), dimension(4) :: phl4

    end type type_hdrt_pheq
    end module module_hdrt_phaselines

    !Module that contains the type for isochoric tracing algorithms
    !Andreas J輍er October 2018
    module module_isoch_therm
    !Type for isochoric thermodynamics (construct px-diagrams using the Helmholtz energy density Psi = A / V = a * rho)
    !Andreas J輍er, October 2018
    !Note: All x-derivatives are taken with all other x constant. That means, that all x are assumed being independent
    type type_isoch_therm

        !The inverse reduced temperature
        double precision:: tau_it
        !Reduced density
        double precision:: delta_it

        !Dimensionless Helmholtz energies and derivatives with respect to delta and tau
        double precision, dimension(:), allocatable:: alpha_it
        !Residual dimensionless Helmholtz energies and derivatives with respect to delta and tau
        double precision, dimension(:), allocatable:: alphar_it
        !Ideal dimensionless Helmholtz energies and derivatives with respect to delta and tau
        double precision, dimension(:), allocatable:: alpha0_it

        !Reducing temperature
        double precision:: Tred_it
        !Derivative of reducing temperature with respect to rho_i at constant rhok
        double precision, dimension(:), allocatable :: dTr_drhoi_it
        !Derivative of reducing temperature with respect to xi at constant xk
        double precision, dimension(:), allocatable :: dTr_dxi_it
        !Second derivative of reducing temperature with respect to rho_i and rho_j at constant rhok
        double precision, dimension(:,:), allocatable :: d2Tr_drhoidrhoj_it
        !Second derivative of reducing temperature with respect to xi and xj at constant xk
        double precision, dimension(:,:), allocatable :: d2Tr_dxidxj_it

        !Reducing density
        double precision:: rhored_it
        !Derivative of reducing density with respect to rho_i at constant rhok
        double precision, dimension(:), allocatable :: drhor_drhoi_it
        !Derivative of reducing density with respect to xi at constant xk
        double precision, dimension(:), allocatable :: drhor_dxi_it
        !Second derivative of reducing density with respect to rho_i and rho_j at constant rhok
        double precision, dimension(:,:), allocatable :: d2rhor_drhoidrhoj_it
        !Second derivative of reducing density with respect to xi and xj at constant xk
        double precision, dimension(:,:), allocatable :: d2rhor_dxidxj_it

        !Derivative of the dimensionless Helmholtz energy with respect to xi at constant tau, delta, and xj
        double precision, dimension(:), allocatable :: dalpha_dxi_it
        !Second derivative of the dimensionless Helmholtz energy with respect to delta and xi at constant tau, xk
        double precision, dimension(:), allocatable :: d2alpha_dxiddel_it
        !Second derivative of the dimensionless Helmholtz energy with respect to tau and xi at constant delta, xk
        double precision, dimension(:), allocatable :: d2alpha_dxidtau_it
        !Second derivative of the dimensionless Helmholtz energy with respect to xi and xj at constant tau, delta, xk
        double precision, dimension(:,:), allocatable :: d2alpha_dxidxj_it

        !Derivative of the residual dimensionless Helmholtz energy with respect to xi at constant tau, delta, and xj
        double precision, dimension(:), allocatable :: dalphar_dxi_it
        !Second derivative of the residual dimensionless Helmholtz energy with respect to delta and xi at constant tau, xk
        double precision, dimension(:), allocatable :: d2alphar_dxiddel_it
        !Second derivative of the residual dimensionless Helmholtz energy with respect to tau and xi at constant delta, xk
        double precision, dimension(:), allocatable :: d2alphar_dxidtau_it
        !Second derivative of the residual dimensionless Helmholtz energy with respect to xi and xj at constant tau, delta, xk
        double precision, dimension(:,:), allocatable :: d2alphar_dxidxj_it

    end type type_isoch_therm
    end module module_isoch_therm

    module help_funcs
    contains
    integer function locate_entry(char, vector) result(pos)

    implicit none
    integer :: iter, e
    character(*), dimension(*) :: vector
    character(*) :: char

    do iter  = 1,len(vector)
        if (trim(vector(iter)) == char) then
            pos = iter
            return
        endif
    enddo
    end function locate_entry
    end module

    !>  F type for alphabetical sort

    ! --------------------------------------------------------------------
    ! INTEGER FUNCTION  FindMinimum():
    !    This function returns the location of the minimum in the section
    ! between Start and End.
    ! --------------------------------------------------------------------

    FUNCTION  FindMinimum(x, Start, End, dim, length)
    IMPLICIT  NONE
    integer::dim, length
    character*length, DIMENSION(dim), INTENT(IN) :: x
    INTEGER, INTENT(IN)                :: Start, End
    character*length                            :: Minimum
    INTEGER                            :: Location
    INTEGER                            :: i,FindMinimum

    Minimum  = x(Start)         ! assume the first is the min
    Location = Start                   ! record its position
    DO i = Start+1, End         ! start with next elements
        IF (x(i) < Minimum) THEN !   if x(i) less than the min?
            Minimum  = x(i)       !      Yes, a new minimum found
            Location = i                !      record its position
        END IF
    END DO
    FindMinimum = Location             ! return the position
    END FUNCTION  FindMinimum

    SUBROUTINE  Swap(a, b, length)
    IMPLICIT  NONE
    integer :: length
    character*length, INTENT(INOUT) :: a, b
    character*length                :: Temp

    Temp = a
    a    = b
    b    = Temp
    END SUBROUTINE  Swap

    SUBROUTINE  Sorter(x, dim, length)
    IMPLICIT  NONE
    integer:: dim, length
    character*length, DIMENSION(dim), INTENT(INOUT) :: x
    INTEGER                               :: i
    INTEGER                               :: Location,FindMinimum

    DO i = 1, dim-1                    ! except for the last
        Location = FindMinimum(x, i, dim, dim, length) ! find min from this to last
        CALL  Swap(x(i), x(Location), length)   ! swap this and the minimum
    END DO
    END SUBROUTINE  Sorter


