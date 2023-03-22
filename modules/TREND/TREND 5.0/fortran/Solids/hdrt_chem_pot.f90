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

    ! module for file hdrt_chem_pot.f90
    module hdrt_chem_pot_module
    !global use inclusion
    use module_all_types
    use calc_functions
    use rhomix_pt_module
    use reduced_parameters_calc_module
    use vle_derivs_module


    contains



    !  Hydrate model - SUBROUTINES
    !  fug_g must be everywhere vector with dimension N_guests, i.e. N_hrts - 1 (water)


    !******************************************************************************
    subroutine hdrt_chem_potent_w(gl,temp,press,fug_g,chpw)
    !******************************************************************************
    !
    !  Chemical potential of water in the hydrate phase [J/mol]
    !  June 2011; Vaclav May 2015
    !  April 2017, Benedikt modified for double ocupancy



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, intent(out) :: chpw
    ! Local variables
    double precision :: gw_B, Cf
    double precision, dimension(3,30) :: CiJ, CiJd
    double precision, dimension(30) :: Y_lang
    double precision, dimension(3) :: CiJ_mix
    integer :: i, J
    !//////////////////////////////////////////////////////////////////

    chpw = 0.d0

    ! Gibbs energy of water in empty lattice [J/mol]
    call hdrt_gibbs_w_B(gl,temp,press,fug_g,gw_B)

    ! Langmuir constant @ T,p [-]
    call hdrt_Langmuir(gl,temp,press,CiJ)


    ! Chemical potential of water in hydrate phase [J/mol]
    Cf = 0.d0
    do i = 1,gl%N_cavi


        if (i==1) then
            Cf = 0.d0
            do J = 1, gl%N_guests
                Cf = Cf + CiJ(i,J)*fug_g(J)
            end do

        elseif (i==2) then
            Cf = 0.d0
            do J = 1,gl%N_guests

                if (dabs(gl%r_vdW(J)) > 1.d-14) then
                    call hdrt_Langmuir_double(gl,temp,press,CiJd)
                    if ((gl%mixdecision=='add' .or. gl%mixdecision=='mul' .or. gl%mixdecision=='mol').and.(gl%N_guests > 1)) then
                        call hdrt_langmuirmix(gl,CiJ, CiJd, CiJ_mix, Y_lang)
                        Cf = Cf + CiJ(i,J)*fug_g(J) + CiJ_mix(i)*fug_g(J)**2.d0        !Langmuir mix
                    else
                        Cf = Cf + CiJ(i,J)*fug_g(J) + CiJ(i,J)*CiJd(i,J)*fug_g(J)**2.d0     !Only one guest	!As in Klauda and Sandler, Asiaee considered (CiJd*fug)**2 Combined from Martin, K&S and Asiaee
                    end if
                else
                    Cf = Cf + CiJ(i,J)*fug_g(J)

                end if

            end do

        end if
        if ((1.d0 + Cf) > 0.d0) then
            chpw = chpw + gl%v_cavi(i)*dlog(1.d0 + Cf)
        else
            chpw = 0.d0
        endif
    end do
    chpw = gw_B - Rgas*temp*chpw

    end subroutine hdrt_chem_potent_w
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    subroutine hdrt_occupancy(gl,temp,press,fug_g,occup,CiJ,occup_single, occup_double)
    !******************************************************************************
    !
    !  Occupancies of hydrate cavities
    !  Vaclav - June 2015, Benedikt April 2017 (doubble occupancy)




    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, dimension(3,30), intent(out) :: occup, CiJ
    ! Local variables
    integer :: i, J, a, b
    double precision, dimension(3) :: CiJ_fug, CiJ_mix
    double precision, dimension(30) :: v_occup, Y_lang, Y_occup, minvec!, occup_ls, occup_ld, occup_sms, occup_smd !occup_ls: single occupied large cavities; occup_ld: double occupied large cavities,
    !occup_sms: single occupied small cavitites, occup_smd: double occupied large cavities
    double precision, dimension(3,30) :: CiJd, occup_double, occup_single
    double precision :: v_occup_all, press_old
    !//////////////////////////////////////////////////////////////////

    occup = 0.d0
    CiJ = 0.d0
    CiJ_fug = 0.d0
    CiJ_mix = 0.d0
    v_occup = 0.d0
    Y_lang = 0.d0
    Y_occup = 0.d0
    v_occup_all = 0.d0
    press_old = 0.d0

    !N_guests=2
    do i = 1,gl%N_guests
        gl%Y_yj_tf(i) = gl%molfractions(i+1)/sum(gl%molfractions(2:))
    end do

    do b=1,200
        minvec=1.d-5
        !! Langmuir constants @ T,p [MPa-1]
        !!call hdrt_Langmuir(temp,press,CiJ)
        !
        !! Cavity occupancies
        !!do i = 1,N_cavi
        !!    CiJ_fug(i) = 0.d0
        !!    do J = 1, N_guests
        !!        CiJ_fug(i) = CiJ_fug(i) + CiJ(i,J)*fug_g(J)
        !!    end do
        !!
        !!    do J = 1, N_guests
        !!        occup(i,J) = CiJ(i,J)*fug_g(J)/( 1.d0 + CiJ_fug(i) )
        !!    end do
        !!
        !!end do
        !
        !!       Adapted for double occupancy Benedikt April 2017
        !
        !occup_ls=0.d0
        !occup_ld=0.d0
        !occup_sms=0.d0
        !occup_smd=0.d0
        !
        !! Langmuir constants @ T,p [MPa-1]
        !call hdrt_Langmuir(temp,press,CiJ)				!Hier oder in i= Schleife?
        !
        !
        !
        !
        !5050 format(f12.5,' ',f16.8,' ',f12.8,' ',f12.8,' ',f12.8,' ',f12.8,' ',f12.8,' ',f12.8,' ',f12.8,' ',f12.8,' ',f12.8)
        !
        !occup_double=0.d0
        !
        !! Cavity occupancies
        !do i = 1,1!N_cavi			!Loop for single occupancy
        !
        !    !if (i==1) then
        !
        !        CiJ_fug(i) = 0.d0
        !
        !        do J = 1, N_guests
        !            CiJ_fug(i) = CiJ_fug(i) + CiJ(i,J)*fug_g(J)
        !        end do
        !
        !        do J = 1, N_guests
        !            occup(i,J) = CiJ(i,J)*fug_g(J)/( 1.d0 + CiJ_fug(i) )
        !        end do
        !
        !    !elseif (i==2) then
        !	   ! CiJ_fug(i) = 0.d0
        !        ! Only considerung the twice occupied catities
        !		! Single occupation defined as before!
        !      !  do J = 1, N_guests
        !		    !if (r_vdW(J) > 1.d-14) then
        !			   ! call hdrt_Langmuir_double(temp,press,CiJd)
        !			  !  CiJ_fug(i) = CiJ_fug(i) + CiJd(i,J)*fug_g(J)**2.d0 !+ CiJ(i,J)*fug_g(J)	!Neuen Langmuir einfügen Double occ, done   + CiJ(i,J)*fug_g(J) +
        !            !else
        !		!	     CiJ_fug(i) = CiJ_fug(i) + CiJ(i,J)*fug_g(J)			!Single occ
        !		    !end if
        !        !end do
        !
        !        !do J = 1, N_guests
        !		    !if (r_vdW(J) > 1.d-14) then
        !			!    call hdrt_Langmuir_double(temp,press,CiJd)    ! ist ja nicht notwendig
        !			!    occup(i,J) = (CiJd(i,J)*(fug_g(J)**2.d0)) /( 1.d0 + CiJ_fug(i))!*fug_g(J) + (CiJ_fug(i)**2.d0 ))  !Neue Occupancy einfügen
        !		    !else
        !		!	    occup(i,J) = CiJ(i,J)*fug_g(J) / ( 1.d0 + CiJ_fug(i) ) ! Da auch große Kavitäten nicht von "Doppelformern" belegt sein müssen.
        !		 !   end if
        !        !end do
        !    !else
        !	!    write(*,*) 'Medium cavity to be considered'
        !
        !    !end if
        !
        !
        !end do
        !
        !do i = 2,N_cavi		! For  double occupated cells.
        !
        !    !if (i==1) then
        !
        !        !CiJ_fug(i)=0
        !
        !   if (i==2) then
        !
        !	    CiJ_fug(i) = 0.d0
        !        ! Only considerung the second occupamd
        !		! Single occupation defined as before! (above and below; depends on cavity number)
        !
        !        do J = 1, N_guests
        !		    if (dabs(r_vdW(J)) > 1.d-14) then
        !			    call hdrt_Langmuir_double(temp,press,CiJd)
        !
        !			      CiJ_fug(i) = CiJ_fug(i) + CiJ(i,J)*CiJd(i,J)*fug_g(J)**2.d0 + CiJ(i,J)*fug_g(J) !Das hier ist es normalerweise
        !                  !CiJ_fug(i) = CiJ_fug(i) + CiJ(i,J)*CiJd(i,J)*press**2.d0 + CiJ(i,J)*press     !using presseure instead of fugacity
        !                !CiJ_fug(i) = CiJ_fug(i) + CiJd(i,J)*CiJd(i,J)*fug_g(J)**2.d0 + CiJ(i,J)*fug_g(J)	!	Klauda & Sandler; Martin	!Neuen Langmuir einfügen Double occ, done   + CiJ(i,J)*fug_g(J) +
        !            else
        !			    CiJ_fug(i) = CiJ_fug(i) + CiJ(i,J)*fug_g(J)
        !		   end if
        !        end do
        !
        !        occup_double=0.d0
        !        do J = 1, N_guests
        !		   if (dabs(r_vdW(J)) > 1.d-14) then
        !			    call hdrt_Langmuir_double(temp,press,CiJd)    ! ist ja nicht notwendig
        !                !occup_single(i,J) = CiJ(i,J)*fug_g(J) / (1.d0 + CiJ_fug(i))
        !                occup_single(i,J) = CiJ(i,J)*fug_g(J) / (1.d0 + CiJ_fug(i))
        !
        !			    occup_double(i,J) = (CiJ(i,J)*CiJd(i,J)*(fug_g(J)**2.d0)) /( 1.d0 + CiJ_fug(i))!*fug_g(J)**2.d0)! + (CiJ_fug(i)**2.d0 ))  !Neue Occupancy einfügen; Normal das hier
        !                !occup_double(i,J) = (CiJd(i,J)*(fug_g(J)**2.d0)) /( 1.d0 + CiJ_fug(i))
        !                occup_ld(J)=occup_double(i,J)
        !		    else
        !
        !                !CiJ_fug(i) = CiJ_fug(i) + CiJ(i,J)*fug_g(J) !passiert ja schon vorher!!!
        !
        !                occup_single(i,J) = CiJ(i,J)*fug_g(J)/( 1.d0 + CiJ_fug(i) )
        !                occup_ls(J)=occup_double(i,J)
        !
        !		    end if
        !        end do
        !    else
        !	    write(*,*) 'Medium cavity to be considered'
        !
        !    end if
        !
        !end do
        !
        !do i=2,N_cavi
        !
        !	do J=1,N_guests
        !
        !		occup(i,J) = occup_single(i,J) + occup_double(i,J)*2.d0
        !
        !	end do
        !
        !end do
        !!a=N_guests
        !if (press<1100) then
        !write(50,5050) press, occup(1,1:N_guests), occup_single(2,1:N_guests), occup_double(2,1:N_guests)
        !end if
        !!temp_occ=temp
        !!press_occ=press
        !!occ_small=occup(1,1)
        !!occ_largesingle=occup_single(2,:)
        !!occ_largedouble=occup_double(2,:)

        !New part of occupancy; no decisions for different hydrate formers; all cavities have the same mathematical approach


        call hdrt_Langmuir(gl,temp, press, CiJ)
        call hdrt_langmuir_double(gl,temp, press, CiJd)
        occup_single=0.d0
        occup_double=0.d0

        do i=1,gl%N_cavi
            CiJ_fug(i)=0.d0

            do J=1,gl%N_guests
                if ((gl%mixdecision=='add' .or. gl%mixdecision=='mul' .or. gl%mixdecision=='mol').and.(gl%numderiv .eqv. .true.)) then
                    !call hdrt_langmuirmix(temp, press,fug_g, CiJ_mix, Y_lang)
                    call hdrt_langmuirmix(gl,CiJ, CiJd, CiJ_mix, Y_lang)
                    CiJ_fug(i) = CiJ_fug(i) + CiJ(i,J)*fug_g(J) + CiJ_mix(i)*fug_g(J)**2.d0            !mixed double occ
                else
                    CiJ_fug(i) = CiJ_fug(i) + CiJ(i,J)*fug_g(J) + CiJ(i,J)*CiJd(i,J)*fug_g(J)**2.d0
                end if
            end do

            do J=1,gl%N_guests		!Not to be placed in fist loop, because of sum CiJ_fug
                if (i==1) then
                    occup_single(i,J) = CiJ(i,J)*fug_g(J) / (1.d0 + CiJ_fug(i))
                    !occup_double(i,J) = CiJ_mix(i)*fug_g(J)**2.d0 / (1.d0+CiJ_fug(i))           !mixed double occ
                    occup_double(i,J) = CiJ(i,J)*CiJd(i,J)*fug_g(J)**2.d0 / (1.d0+CiJ_fug(i))
                else if (i==2) then
                    occup_single(i,J)=CiJ(i,J)*fug_g(J) / (1.d0 + CiJ_fug(i))
                    !occup_double(i,J) = CiJ_mix(i)*fug_g(J)**2.d0 / (1.d0+CiJ_fug(i))           !mixed double occ
                    occup_double(i,J)=CiJ(i,J)*CiJd(i,J)*fug_g(J)**2.d0 / (1.d0 + CiJ_fug(i))
                else
                    !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
                    write(*,*)'Medium cavity of structure H not implementated yet'
                    !DEC$ END IF ! WO_WRITE_TO_CONSOLE

                end if
            end do

            occup(i,:) = occup_single(i,:) + 2.d0*occup_double(i,:)


        end do	! i-loop (cavities)

        v_occup = 0.d0
        v_occup_all = 0.d0
        do J = 1, gl%N_guests
            do i = 1,gl%N_cavi
                v_occup(J) = v_occup(J) + (gl%v_cavi(i)*occup(i,J))
            end do
            v_occup_all = v_occup_all + v_occup(J)
        end do

        do J = 1, gl%N_guests
            Y_occup(J) = v_occup(J)/v_occup_all
        end do



        !call hdrt_YJ_Tf(temp, press, fug_g, Y)


        if ((gl%N_guests==1).or.(gl%numderiv .eqv. .false.)) then
            return
        end if

        if ((dabs((gl%Y_yj_tf(1)-Y_occup(1))/(gl%Y_yj_tf(1)))<1.d-6) .and. (dabs((gl%Y_yj_tf(2)-Y_occup(2))/(gl%Y_yj_tf(2)))<1.d-6) ) then

            return
        else
            gl%Y_yj_tf = Y_occup
        end if



        if(b>190) then
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
            write(*,*)'convergence Problem during calculation of gl%langmuir mix', b
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        end if

    end do !(b loop)



    end subroutine hdrt_occupancy
    !******************************************************************************
    !******************************************************************************



    !******************************************************************************
    subroutine hdrt_mole_fract(gl,temp,press,fug_g,occup,CiJ,xH,occup_single, occup_double)
    !******************************************************************************
    !
    !  Mole fraction of gas (2, ... J) and water (1) in the hydrate phase [-]
    !  June 2011; Vaclav May 2015



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, dimension(30), intent(out) ::  xH
    double precision, dimension(3,30), intent(out) :: occup, CiJ, occup_single, occup_double
    ! Local variables
    integer :: i, J
    double precision, dimension(3) :: CiJ_fug
    double precision, dimension(30) :: v_occup
    double precision :: v_occup_all
    !//////////////////////////////////////////////////////////////////

    xH = 0.d0
    occup = 0.d0
    CiJ = 0.d0
    occup_single = 0.d0
    occup_double = 0.d0

    ! Cavity occupancies
    call hdrt_occupancy(gl,temp,press,fug_g,occup,CiJ,occup_single, occup_double)

    ! Hydrate phase mole fraction
    v_occup_all = 0.d0
    do J = 1, gl%N_guests
        v_occup(J) = 0.d0
        do i = 1,gl%N_cavi
            v_occup(J) = v_occup(J) + gl%v_cavi(i)*occup(i,J)
        end do
        v_occup_all = v_occup_all + v_occup(J)
    end do

    xH(1) = 1.d0/(1.d0 + v_occup_all) ! H2O
    do J = 1, gl%N_guests
        xH(J+1) = v_occup(J)/(1.d0 + v_occup_all)   ! guest (H2O is 1st)
    end do


    end subroutine hdrt_mole_fract
    !******************************************************************************
    !******************************************************************************


    !******************************************************************************
    subroutine hdrt_gibbs_w_B(gl,temp,press,fug_g,gw_B)
    !******************************************************************************
    !
    !  Gibbs energy of water in reference empty hydrate lattice [J/mol]
    !  June 2011; 2015, Vaclav



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, intent(out) :: gw_B
    ! Local variables
    double precision :: int_hw, int_vol, int_vol_num
    double precision :: hw_B, a0_m_YTf, a_mix
    integer :: J, check_B0_const
    !//////////////////////////////////////////////////////////////////


    ! Enthalpy of empty lattice [J.mol^-1]
    hw_B = gl%cpH_a/2.d0*(temp**2-gl%T0**2) + gl%cpH_b*(temp-gl%T0) + gl%hw_B0                   ! <- Old version with reference conditions as defined by Ballard and Sloan, Andreas Oct 2013

    ! Integral (hw_B/temp**2)
    int_hw =  -gl%hw_B0*(1.d0/temp-1.d0/gl%T0) + gl%cpH_a/2.d0*(temp-gl%T0) + &               ! <- Old version with reference conditions as defined by Ballard and Sloan, Andreas Oct 2013
        gl%cpH_a*gl%T0**2/2.d0*(1.d0/temp-1.d0/gl%T0) + gl%cpH_b*dlog(temp/gl%T0) + &
        gl%cpH_b*gl%T0*(1.d0/temp-1.d0/gl%T0)

    ! Mixed Hydrates
    ! Check if bulk modulus B0_hdrt is constant
    check_B0_const = 1;
    do J = 1,gl%N_guests
        if (dabs(gl%B0_hdrt(1,J)-10.d9).gt.1.d3) then
            check_B0_const = 0  ! B0_hdrt is not constant for all hydrate formers
        end if
    end do

    if ((gl%Latt_par == 'm4').and.(check_B0_const.eq.1)) then
        ! Lattice parameter @ T0_a, p0_a of mixed hydrate @ T & p
        call hdrt_latt_par_mix_0_Y_Tf(gl,temp,press,fug_g,a0_m_YTf)
        !call hdrt_latt_par_mix(temp,press,fug_g,a_mix)          !only for testing
        !write(*,*) (a0_m_YTf - a_mix) / a_mix                   !see above


        ! Analytical integration of the molar volume vol_B - correlation based on Murnaghan EoS with dB0 = 4 = const.
        int_vol = ( ((gl%B0_m(1)+gl%B0_m(2)*press)/(gl%B0_m(1)+gl%B0_m(2)*gl%p0_a))**((gl%B0_m(2)-1)/gl%B0_m(2)) - &
            ((gl%B0_m(1)+gl%B0_m(2)*gl%p0)/(gl%B0_m(1)+gl%B0_m(2)*gl%p0_a))**((gl%B0_m(2)-1)/gl%B0_m(2)) ) * (gl%B0_m(1)+gl%B0_m(2)*gl%p0_a)/(gl%B0_m(2)-1)
        int_vol = int_vol * Nav/(gl%Nw*temp)*(a0_m_YTf*1.d-10)**3*dexp( 3.d0*gl%alfa(1)*(temp-gl%T0_a)  &
            + 3.d0*gl%alfa(2)*(temp-gl%T0_a)**2 + 3.d0*gl%alfa(3)*(temp-gl%T0_a)**3 ) ! Integral(vol_B)/temp

    else    ! non-constant bulk modulus or other definition of pressure dependence

        ! Numerical integration of the molar volume vol_B
        ! ... problem with evaluation of  B0_m for mixture, therefore numerical integration required
        !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
        write(*,*) 'Warning - volume integral in hdrt_gibbs_w_B calculated numerically'
        write(*,*) '        - bulk modulus is not the same for all involved hydrates'
        !DEC$ END IF ! WO_WRITE_TO_CONSOLE
        call hdrt_int_volB_num(gl,temp,press,fug_g, int_vol) ! Integral(vol_B)
        int_vol = int_vol/temp  ! Integral(vol_B)/temp

    end if

    !! Analytical integral of volume checked numerically for mixed hydrates - OK - February 2016
    !call hdrt_int_volB_num(temp,press,fug_g, int_vol_num) ! Integral(vol_B)
    !    int_vol_num = int_vol_num/temp  ! Integral(vol_B)/temp

    ! Gibbs energy of water in empty lattice [J/mol]
    gw_B = temp*( gl%gw_B0/gl%T0 - int_hw + int_vol )    !<- Old version with reference conditions as defined by Ballard and Sloan, Andreas Oct 2013

    end subroutine hdrt_gibbs_w_B
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    subroutine hdrt_int_volB_num(gl,temp,press,fug_g,int_vol)
    !******************************************************************************
    !
    !  Numerical integration of the MOLAR VOLUME of empty water lattice - integral(vol_B)
    !  January 2012, Vaclav October 2015
    !  March 2017, Vaclav - fug_g added as an input to hdrt_latt_vol_B



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    integer, parameter :: p_dimm = 350
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, intent(out) :: int_vol
    ! Local variables
    integer :: i
    double precision :: a_latt
    double precision, dimension(p_dimm) :: p, volume
    !//////////////////////////////////////////////////////////////////

    ! Pressure vector
    p(1) = gl%p0 ! [Pa] reference pressure for Gibbs energy
    do i = 2,p_dimm
        p(i) = p(i-1) + (press-p(1))/(p_dimm-1)
    end do

    call hdrt_latt_vol_B(gl,temp,p(1),fug_g,a_latt,volume(1)) ! vol_B @ p0

    ! Trapezoidal numerical integration of molar volume vol_B
    ! -------------------------------------------------------
    int_vol = 0.d0
    do i = 1,p_dimm-1
        call hdrt_latt_vol_B(gl,temp,p(i+1),fug_g,a_latt,volume(i+1))
        int_vol = int_vol + (p(i+1)-p(i))*(volume(i+1)+volume(i))/2.d0
    end do

    end subroutine hdrt_int_volB_num
    !******************************************************************************
    !******************************************************************************



    !******************************************************************************
    subroutine hdrt_latt_par_mix_0_Y_Tf(gl,temp,press,fug_g, a0_m_YTf)
    !******************************************************************************
    !
    !  Mixing rule for the lattice parameter [A] @ T0_a & p0_a based on YJ(T,fug)
    !  a) ZERO pressure (fug_g = 0) and temperature 'temp'
    !  b) defined fugacity 'fug_g' and temperature 'temp'
    !  February 2016, Vaclav
    !  March 2016, Sebastian
    !  March 2017, Vaclav



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, intent(out) :: a0_m_YTf
    ! Local variables
    integer :: J
    double precision, dimension(30) :: Y
    !//////////////////////////////////////////////////////////////////

    ! Pure hydrate
    if (gl%N_guests == 1) then
        a0_m_YTf = gl%a0_hdrt(1)

        ! Mixed hydrates
    else

        ! Composition of guests @ temperature 'temp' and fugacity 'fug_g'
        ! ---------------------------------------------------------------
        !HERE IF FUGACITY == 0
        if (maxval(fug_g) < 1.d-12) then
            call hdrt_YJ_zero_press(gl,temp,Y)
        else
            call hdrt_YJ_Tf(gl,temp,press,fug_g,Y)
        endif
        !HERE IF FUGACITY /= 0
        !new mixing rule

        ! Lattice parameter of mixed hydrate @ T0_a and p0_a
        ! --------------------------------------------------

        if (gl%Latt_par_mix == 'volume') then
            a0_m_YTf = 0.d0
            do J = 1,gl%N_guests
                a0_m_YTf = a0_m_YTf + Y(J)*gl%a0_hdrt(J)**3
            end do
            a0_m_YTf = a0_m_YTf**(1.d0/3.d0)  ! [A] lattice parameter of mixed hydrate @ T0_a and p0_a

        elseif (gl%Latt_par_mix == 'vegard') then
            a0_m_YTf = 0.d0
            do J = 1,gl%N_guests
                a0_m_YTf = a0_m_YTf + Y(J)*gl%a0_hdrt(J) ! [A] lattice parameter of mixed hydrate @ T0_a and p0_a
            end do

        end if
    end if

    end subroutine hdrt_latt_par_mix_0_Y_Tf
    !******************************************************************************
    !******************************************************************************



    !******************************************************************************
    subroutine hdrt_latt_par_mix(gl,temp,press,fug_g,a_mix)
    !******************************************************************************
    !
    !  Lattice parameter of the mixed hydrate
    !  February 2016, Vaclav
    !  March 2016, Vaclav - update for a_M0(T0,p0) = sum( a0J(T0,p0)*YJ(T,fug_g) )



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp,press
    ! Output arguments
    double precision, intent(out) :: a_mix
    ! Local variables
    integer :: J, i, check_B0_const
    double precision :: a_M0
    double precision, dimension(30) :: Y, a_pure, fug_g
    !//////////////////////////////////////////////////////////////////

    ! Lattice parameter of mixed hydrate a0_hdrt @ T0_a & p0_a
    ! ---------------------------------------------------
    call hdrt_latt_par_mix_0_Y_Tf(gl,temp, press, fug_g, a_M0)

    ! Lattice parameter of mixed hydrate @ T,p [A = 1.d-10 m]
    ! -------------------------------------------------------
    ! Check if B0_hdrt is constant
    check_B0_const = 1;
    do J = 1,gl%N_guests
        if (dabs(gl%B0_hdrt(1,J)-10.d9).gt.1.d3) then
            check_B0_const = 0  ! B0_hdrt is not constant for all hydrate formers
        end if
    end do

    if ((gl%Latt_par == 'm4').and.(check_B0_const.eq.1)) then
        ! Our correlation based on Murnaghan EoS with dB0 = 4 = const. ... p0_a = 0.d0
        ! Constant bulk modulus and alfa-s for all hydrates of the same structure
        a_mix = a_M0*dexp( gl%alfa(1)*(temp-gl%T0_a) + gl%alfa(2)*(temp-gl%T0_a)**2 + gl%alfa(3)*(temp-gl%T0_a)**3 )* &
            ( (gl%B0_m(1)+gl%B0_m(2)*press)/(gl%B0_m(1)+gl%B0_m(2)*gl%p0_a) )**( -1.d0/(3.d0*gl%B0_m(2)) )    ! [A] lattice param.

    else    ! non-constant bulk modulus or other definition of pressure dependence

        do i = 1,gl%N_guests
            !if (Latt_par == 'vj') then
            !! New compress. - our formulation for pure CO2 - kapa 1,2
            !a_pure(i) = a0_hdrt(i)*dexp( alfa(1)*(temp-T0_a) + alfa(2)*(temp-T0_a)**2 + alfa(3)*(temp-T0_a)**3 )* &
            !            (1.d0 - kapa(1)*(press-p0_a)/(1.d0+kapa(2)*(press-p0_a)) )    ! [A] lattice param.
            !
            !elseif (Latt_par == 'bs') then
            !! Original correlation of Ballard and Sloan - exponential kapa1
            !a_pure(i) = a0_hdrt(i)*dexp( alfa(1)*(temp-T0_a) + alfa(2)*(temp-T0_a)**2 + alfa(3)*(temp-T0_a)**3 - &
            !            kapa(1)*(press-p0_a) )    ! [A] lattice param.
            !
            !elseif (Latt_par == 'm4') then
            ! Our correlation based on Murnaghan EoS with dB0 = 4 = const. ... p0_a = 0.d0
            a_pure(i) = gl%a0_hdrt(i)*dexp( gl%alfa(1)*(temp-gl%T0_a) + gl%alfa(2)*(temp-gl%T0_a)**2 + gl%alfa(3)*(temp-gl%T0_a)**3 )* &
                ( (gl%B0_hdrt(1,i)+gl%B0_hdrt(2,i)*press)/(gl%B0_hdrt(1,i)+gl%B0_hdrt(2,i)*gl%p0_a) )**( -1.d0/(3.d0*gl%B0_hdrt(2,i)) )    ! [A] lattice param.
            !end if
        end do

        ! Composition of guests @ 'temp' and fug_g
        ! -----------------------------------------------------
        !call hdrt_YJ_zero_press(temp,Y)
        call hdrt_YJ_Tf(gl,temp, press, fug_g, Y)

        ! Mixed lattice param. for non-constant B0_hdrt @ pressure 'press' & temperature 'temp'
        ! --------------------------------------------------------------------------------
        if (gl%Latt_par_mix == 'volume') then
            a_mix = 0.d0
            do J = 1,gl%N_guests
                a_mix = a_mix + Y(J)*a_pure(J)**3
            end do
            a_mix = a_mix**(1.d0/3.d0)  ! [A] lattice parameter of mixed hydrate @ Temp and 0 pressure

        elseif (gl%Latt_par_mix == 'vegard') then
            a_mix = 0.d0
            do J = 1,gl%N_guests
                a_mix = a_mix + Y(J)*a_pure(J) ! [A] lattice parameter of mixed hydrate @ Temp and 0 pressure
            end do

        end if

    end if  ! end of constant bulk modulus loop

    end subroutine hdrt_latt_par_mix
    !******************************************************************************
    !******************************************************************************



    !******************************************************************************
    subroutine hdrt_latt_vol_B(gl,temp,press,fug_g, a_latt,vol_B)
    !******************************************************************************
    !
    !  Lattice parameter of final (mixed) hydrate [A]
    !  Molar volume of the empty hydrate lattice [m3/mol]
    !  February 2016, Vaclav
    !  March 2017, Vaclav - fug_g added as an input to hdrt_latt_par_mix



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, intent(out) :: a_latt, vol_B
    !//////////////////////////////////////////////////////////////////

    ! Lattice parameter of mixed hydrate [A = 1.d-10 m]
    ! -------------------------------------------------
    call hdrt_latt_par_mix(gl,temp,press,fug_g,a_latt)

    ! Molar volume of the empty lattice @ T,p [m3/mol]
    ! ------------------------------------------------
    vol_B = Nav/gl%Nw*(a_latt*1.d-10)**3   ! [m3/mol] molar volume of empty beta-lattice

    end subroutine hdrt_latt_vol_B
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    subroutine hdrt_density(gl,temp,press,fug_g,occup,CiJ,rho_H)
    !******************************************************************************
    !
    !  Molar density of gas hydrate [mol/m^3]
    !  Andreas March 2014; Vaclav May 2015, February 2016
    !  March 2017 - fug_g added as an input to hdrt_latt_vol_B




    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in)  :: temp, press
    double precision, dimension(30), intent(in)  :: fug_g
    !double precision, dimension(30) :: occup_ls, occup_ld, occup_sms, occup_smd
    ! Output arguments
    double precision, intent(out) ::  rho_H
    double precision, dimension(3,30) :: occup, CiJ,occup_single, occup_double
    !Auxiliary variables
    integer :: i, J
    double precision :: a_latt, vol_B, v_G
    !//////////////////////////////////////////////////////////////////

    rho_H = 0.d0
    v_G = 0.d0

    ! Molar volume of the empty lattice @ T,p [m3/mol]
    ! ------------------------------------------------
    call hdrt_latt_vol_B(gl,temp,press,fug_g,a_latt,vol_B)

    ! Occupancies of hydrate cavities
    ! -------------------------------
    call hdrt_occupancy(gl,temp,press,fug_g,occup,CiJ,occup_single, occup_double)

    do i = 1, gl%N_cavi
        do J = 1, gl%N_guests
            v_G = v_G + gl%v_cavi(i) * occup(i,J)
        end do
    end do


    rho_H = (1.D0/vol_B)*(1.d0 + v_G)  ! [mol.m^-3] molar density of hydrate

    end subroutine hdrt_density
    !******************************************************************************
    !******************************************************************************



    subroutine hdrt_YJ_zero_press(gl,temp,Y)
    !******************************************************************************
    !
    !  Guest composition at ZERO pressure and temperature 'temp'
    !  Y(J) from limit p --> 0 Pa (L'Hospital rule)
    !
    !  Vaclav July, October 2015, February 2016
    !  Benedikt July 2017, considering double occupancy
    !


    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp
    ! Output arguments
    double precision, dimension(30), intent(out) :: Y
    ! Local variables
    double precision :: v_CiJ_xJ_all, press
    double precision, dimension(3,30) :: CiJ, CiJd
    integer :: i,J
    double precision, dimension(30) :: v_CiJ_xJ
    !//////////////////////////////////////////////////////////////////

    press = 0.d0    ! 0 Pa

    call hdrt_Langmuir(gl,temp,press,CiJ)
    call hdrt_Langmuir_double(gl,temp, press, CiJd)

    ! Mole fraction of guests Y(J) - limit solution @ zero pressure
    v_CiJ_xJ_all = 0.d0
    do J = 1, gl%N_guests
        v_CiJ_xJ(J) = 0.d0
        do i = 1,gl%N_cavi
            v_CiJ_xJ(J) = v_CiJ_xJ(J) + gl%v_cavi(i)*CiJ(i,J)*gl%moles_hdrt(J+1)		!correct implementation has to be checked! no term for double occ needed. (L' Hospital rule)
            ! in moles_hdrt water is the 1st component
        end do
        v_CiJ_xJ_all = v_CiJ_xJ_all + v_CiJ_xJ(J)
    end do

    do J = 1, gl%N_guests
        Y(J) = v_CiJ_xJ(J)/v_CiJ_xJ_all   ! gas (H2O is 1st)
    end do

    gl%Y_yj_tf=Y       !storing Y in the module hdrt_property_definition to compare it in hdrt_occupancy with Langmiurmix
    end subroutine hdrt_YJ_zero_press
    !******************************************************************************
    !******************************************************************************


    !******************************************************************************
    subroutine hdrt_YJ_Tf(gl,temp, press, fug_g, Y)
    !******************************************************************************
    !
    !  Guest composition at specified temperature 'temp' and fugacity 'fug_g'
    !  Y does not actually depend on pressure
    !
    !  March 2016, Sebastian
    !  March 2017, Vaclav -



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, dimension(30), intent(out) :: Y
    ! Local variables
    double precision :: CiJ_fug_all, v_CiJ_fug_all
    double precision, dimension(3,30) :: CiJ
    integer :: k,i,m,J
    double precision, dimension(30) :: v_CiJ_fug!, occup_ls, occup_ld, occup_sms, occup_smd

    double precision, dimension(3,30) :: occup,occup_single, occup_double
    double precision, dimension(30) :: v_occup, Y2
    double precision :: v_occup_all
    !//////////////////////////////////////////////////////////////////

    ! WRONG to be deleted
    ! *******************
    !call hdrt_Langmuir(temp,press,CiJ)
    !
    !CiJ_fug_all = 0.d0
    !do J = 1, N_guests
    !    do i = 1,N_cavi
    !        CiJ_fug_all = CiJ_fug_all + CiJ(i,J)*fug_g(J)
    !    enddo
    !enddo
    !
    !! Mole fraction of guests Y(J) - limit solution @ zero pressure
    !v_CiJ_fug_all = 0.d0
    !v_CiJ_fug = 0.d0
    !do J = 1, N_guests
    !    do i = 1,N_cavi
    !        v_CiJ_fug(J) = v_CiJ_fug(J) + (v_cavi(i)*CiJ(i,J)*fug_g(J))/(1 + CiJ_fug_all)   ! CiJ_fug_all(i) shall be vector, I am afraid
    !    end do
    !    v_CiJ_fug_all = v_CiJ_fug_all + v_CiJ_fug(J)
    !end do
    !
    !do J = 1, N_guests
    !    Y(J) = v_CiJ_fug(J)/v_CiJ_fug_all   ! gas (H2O is 1st)
    !end do


    ! Cage occupancies
    call hdrt_occupancy(gl,temp,press,fug_g,occup,CiJ,occup_single, occup_double)

    ! Mole fraction of guests Y(J) at given temperature and fugacity
    v_occup = 0.d0
    v_occup_all = 0.d0
    do J = 1, gl%N_guests
        do i = 1,gl%N_cavi
            v_occup(J) = v_occup(J) + (gl%v_cavi(i)*occup(i,J))
        end do
        v_occup_all = v_occup_all + v_occup(J)
    end do

    do J = 1, gl%N_guests
        Y(J) = v_occup(J)/v_occup_all   ! gas (H2O is 1st)
    end do

    end subroutine hdrt_YJ_Tf
    !******************************************************************************
    !******************************************************************************


    !  CALL THIS SUBROUTINE EVERYTIME OVERALL COMPOSITION IS CHANGED!!!
    !******************************************************************************
    subroutine hdrt_ref_a_gw0_hw0_COC(gl)
    !******************************************************************************
    !
    !  Subroutine needs to be called only once for a Constant Overall Composition = COC
    !
    !  Lattice parameter of final (gl,mixed) hydrate a_T_ref [A] @ T_ref and p_ref = 0 MPa
    !  for evaluation of corrected shell radii Ra = R_shell*a_T_ref/a_ref
    !  T_ref(sI) = 248.15 K, T_ref(sII) = 253.15 K
    !
    !  Lattice parameter at zero pressure & T0_a, i.e. a0_m, needed for:
    !  Evaluation of the reference state parameters gw_B0 and hw_B0
    !
    !  Vaclav October 2015, February 2016



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    ! T_ref and p_ref = 0Pa defined in module 'hdrt_property_definition'
    ! T0_a and p0_a defined in module 'hdrt_property_definition'

    ! Output arguments
    ! a_T_ref defined in module 'hdrt_property_definition' @ T_ref and p_ref = 0 MPa
    ! gw_B0 and hw_B0 defined in module 'hdrt_property_definition'

    ! Local variables
    double precision :: temp, press, devi, a_latt
    integer :: J,k
    double precision, dimension(30) :: a_pure, Y, Y_0_press, fug_g
    double precision :: a_0_press
    !//////////////////////////////////////////////////////////////////
    a_0_press = 0.d0
    !SH: 11/2017
    !since the Gibbs Energy and Enthalpy are not dependend on the lattice paramter, the calculation of the lattice
    !parameter at reference conditions is not needed anymore. If the dependency on the lattice parameter will
    !ever be reactivated the block in the following statement is used again
    if ((dabs(gl%kb_gw(1)) > 1.d-14) .and. (dabs(gl%kb_hw(1)) > 1.d-14)) then
        ! Mixed lattice parameter a_T_ref @ T_ref, p_ref = 0 MPa for shell radii correction
        ! ---------------------------------------------------------------------------------
        ! T_ref & p_ref are constants for sI and sII hydrates defined in hdrt_modules

        ! Pure hydrate
        if (gl%N_guests == 1) then
            fug_g(1) = gl%p_ref
            call hdrt_latt_par_mix(gl,gl%T_ref,gl%p_ref,fug_g, gl%a_T_ref)
        else

            ! Mixed hydrate -> iterative solution of a_T_ref
            ! Initial estimate for a_T_ref determined from the overall composition
            !fug_g = 0.d0
            fug_g(1:29) = gl%moles_hdrt(2:30)*gl%p_ref    ! fugacity = partial pressure
            call hdrt_latt_par_mix_0_Y_Tf(gl,gl%T_ref,gl%p_ref,fug_g,gl%a_T_ref)

            devi = 1.d0
            do k = 1,500
                ! Lattice parameter of mixed hydrate [A = 1.d-10 m] @ T_ref & p_ref = 0 MPa
                ! -------------------------------------------------------------------------
                call hdrt_latt_par_mix(gl,gl%T_ref,gl%p_ref,fug_g,a_latt)

                devi = dabs(a_latt - gl%a_T_ref)
                gl%a_T_ref = a_latt

                if (devi.le.1.d-7) then
                    exit
                end if

                if (k.ge.495) then
                    !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)
                    write(*,*) 'Problem with gl%a_T_ref convergence!'
                    !DEC$ END IF ! WO_WRITE_TO_CONSOLE
                end if

            end do

        end if  ! pure/mixed hydrate

        ! Lattice parameter of mixed hydrate
        ! @ conditions for the reference lattice parameter T0_a & p0_a
        ! -------------------------------------------------------------
        !call hdrt_latt_param_0_press(T0_a,a_0_press,Y)
        !fug_g = 0.d0
        fug_g(1:29) = gl%moles_hdrt(2:30)*gl%p0_a    ! fugacity = partial pressure
        call hdrt_latt_par_mix_0_Y_Tf(gl,gl%T0_a,gl%p0_a,fug_g,a_0_press)

        gl%a0_m = a_0_press    ! Reference lattice parameter @ T0_a and p0_a
    endif

    ! Reference state parameters
    ! gw_B0 and hw_B0 from linear function of unit cell volume at T0_a & p0_a
    ! -----------------------------------------------------------------------
    gl%gw_B0 = gl%kb_gw(1)*a_0_press**3 + gl%kb_gw(2)     ! Reference state Gibbs energy of water in empty lattice [J/mol]
    gl%hw_B0 = gl%kb_hw(1)*a_0_press**3 + gl%kb_hw(2)     ! Reference state enthalpy of empty lattice [J.mol^-1]



    end subroutine hdrt_ref_a_gw0_hw0_COC
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    subroutine hdrt_Langmuir(gl,temp,press,CiJ)
    !******************************************************************************
    !
    !  Langmuir constant for all cavities [-]
    !  June 2011; May 2015 Vaclav



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    ! Output arguments
    double precision, dimension(3,30), intent(out) :: CiJ
    ! Local variables
    double precision :: pi_loc, a_latt, vol_B
    double precision, dimension(6) :: Ra  ! [A] shell radius dependent on lattice param. a
    integer, dimension(6) :: z   ! [-] coordination number
    double precision :: int_wr, C_num
    integer :: i, J
    !//////////////////////////////////////////////////////////////////

    CiJ = 0.d0

    do J = 1,gl%N_guests

        if (gl%Langmuir == 'ks') then  ! Langmuir constant from Klauda & Sandler T-fit

            CiJ(1,J) = dexp(gl%KS_s(1,J) + gl%KS_s(2,J)/temp + gl%KS_s(3,J)/temp**2)
            CiJ(2,J) = dexp(gl%KS_l(1,J) + gl%KS_l(2,J)/temp + gl%KS_l(3,J)/temp**2)

            if (gl%KS_s(1,J) == 0.d0) then
                CiJ(1,J) = 0.d0
            elseif (gl%KS_l(1,J) == 0.d0) then
                CiJ(2,J) = 0.d0
            end if

        else ! Langmuir const. from BS or JH

            pi_loc = 3.141592653589793d0

            do i = 1,gl%N_cavi

                ! i-th cavity
                !if (Ra_def == 'r(a)') then
                !Ra = R_shell(i,1:6)*a_latt/a_ref ! [A] shell radius dependent on lattice param. a(T,p)
                if (gl%Ra_def == 'tref') then
                    Ra = gl%R_shell(i,1:6)*gl%a_T_ref/gl%a_ref ! [A] shell radius DEPENDENT only on a(T_ref)
                elseif (gl%Ra_def == 'cons') then
                    Ra = gl%R_shell(i,1:6) ! [A] shell radius INDEPENDENT of hydrate type
                end if

                z = gl%z_shell(i,1:6)      ! [-] coordination number at i-th cavity

                call hdrt_int_pot_func(gl,temp,Ra,z,J,int_wr)

                CiJ(i,J) = 4.d0*pi_loc/kB/temp*int_wr

            end do

        end if

    end do  ! do-loop N_guests (J)

    end subroutine hdrt_Langmuir
    !******************************************************************************
    !******************************************************************************



    !******************************************************************************
    subroutine hdrt_Langmuir_double(gl,temp,press,CiJd)
    !******************************************************************************
    !
    !  Function for second Langmiur constant for double occupancy, only for sII Hydrates, considering only i=2 cavity
    !  April 2017, Benedikt, noch kopiert aus Langmiur single



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    ! Output arguments
    double precision, dimension(3,30), intent(out) :: CiJd
    ! Local variables
    double precision :: pi_loc, a_latt, vol_B
    double precision, dimension(6) :: Ra  ! [A] shell radius dependent on lattice param. a
    !double precision, dimension(3,30) :: CiJd
    integer, dimension(6) :: z   ! [-] coordination number
    double precision :: int_wr, C_num, int_wrd                           !Neues Potentialintegral beachten
    integer :: i, J
    !//////////////////////////////////////////////////////////////////

    CiJd = 0.d0

    do J = 1,gl%N_guests

        if (dabs(gl%r_vdW(J)) > 1.d-14) then

            !if (Langmuir == 'ks') then  ! Langmuir constant from Klauda & Sandler T-fit
            !
            !   CiJ(1,J) = dexp(KS_s(1,J) + KS_s(2,J)/temp + KS_s(3,J)/temp**2)
            !  CiJ(2,J) = dexp(KS_l(1,J) + KS_l(2,J)/temp + KS_l(3,J)/temp**2)

            !  if (KS_s(1,J) == 0.d0) then
            !        CiJ(1,J) = 0.d0
            !   elseif (KS_l(1,J) == 0.d0) then
            !        CiJ(2,J) = 0.d0
            !  end if

            !else ! Langmuir const. from BS or JH

            pi_loc = 3.141592653589793d0

            do i = 1,gl%N_cavi

                if (i==1) then
                    CiJd(i,J)=0.d0

                elseif (i==2) then
                    CiJd(i,J)=0.d0
                    ! i-th cavity
                    !if (Ra_def == 'r(a)') then
                    !Ra = R_shell(i,1:6)*a_latt/a_ref ! [A] shell radius dependent on lattice param. a(T,p)
                    if (gl%Ra_def == 'tref') then
                        Ra = gl%R_shell(i,1:6)*gl%a_T_ref/gl%a_ref ! [A] shell radius DEPENDENT only on a(T_ref)
                    elseif (gl%Ra_def == 'cons') then
                        Ra = gl%R_shell(i,1:6) ! [A] shell radius INDEPENDENT of hydrate type
                    end if

                    z = gl%z_shell(i,1:6)      ! [-] coordination number at i-th cavity

                    call hdrt_int_pot_double_func(gl,temp,Ra,z,J,int_wrd)                      !neue Potentialfunktion einfügen
                    !all hdrt_int_pot_func_double_int(temp,Ra,z,J,int_wrd)

                    CiJd(i,J) = 4.d0*pi_loc/kB/temp*int_wrd       !Integral für neues Potential einfügen

                else
                    CiJd(i,J)=0

                end if

            end do

        else
            continue
        end if

        !end if

    end do  ! do-loop N_guests (J)

    end subroutine hdrt_Langmuir_double
    !******************************************************************************
    !******************************************************************************


    !******************************************************************************
    !subroutine hdrt_langmuirmix(temp, press,fug_g, CiJ_mix, Y_lang)
    subroutine hdrt_langmuirmix(gl,CiJ, CiJd, CiJ_mix, Y_lang)
    !******************************************************************************
    !
    ! Benedikt July 2017; Mixing rule for Cell Potential function of double occupied cells
    ! considering composition of the hydrate to mix the Langmiurconstants for double occupancy (mol)
    ! other switches (add, mul) only for comparison
    ! Routine works also, when only one guest is in the hydrate
    !



    !/////////////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl

    double precision, dimension (3,30), intent(in) :: CiJ, CiJd
    !double precision, intent(in) :: temp, press
    integer :: i,J
    double precision, dimension (3,30) :: CiJ_mix_store
    double precision, dimension (30) :: Y_lang, fug_g
    double precision, dimension (3) :: CiJ_mix
    double precision :: N_guests_double, sum_moles_hdrt, sumY_yj
    !//////////////////////////////////////////////////////////////////

    !N_guests=2
    if (gl%N_guests>1) then

        !call hdrt_Langmuir(temp, press,CiJ)
        !call hdrt_langmuir_double(temp, press, CiJd)

        if (gl%mixdecision=='mul') then
            Cij_mix(1:gl%N_cavi)=1.d0
            Cij_mix_store=0.d0

            do i=1,gl%N_cavi
                do J=1,gl%N_guests
                    CiJ_mix_store(i,J)=CiJ(i,J)*CiJd(i,J)
                end do

                do J=1,gl%N_guests
                    CiJ_mix(i)=CiJ_mix(i)*CiJ_mix_store(i,J)
                end do

                N_guests_double=dble(gl%N_guests)
                CiJ_mix(i)=(Cij_mix(i))**(1.d0/N_guests_double)
            end do

        elseif(gl%mixdecision=='add') then

            Cij_mix(1:gl%N_cavi)=0.d0
            Cij_mix_store=0.d0

            do i=1,gl%N_cavi
                do J=1,gl%N_guests
                    CiJ_mix_store(i,J)=CiJ(i,J)*CiJd(i,J)
                end do

                do J=1,gl%N_guests
                    CiJ_mix(i)=CiJ_mix(i)+CiJ_mix_store(i,J)
                end do

                N_guests_double=dble(gl%N_guests)

                CiJ_mix(i)=(Cij_mix(i)) / N_guests_double
            end do

        elseif(gl%mixdecision=='mol') then

            Cij_mix(1:gl%N_cavi)=0.d0
            Cij_mix_store=0.d0
            sum_moles_hdrt=0.d0
            !do J=1,N_guests
            !    sum_moles_hdrt  = sum_moles_hdrt + moles_hdrt(J+1)
            !end do
            sumY_yj=0.d0
            do J=1,gl%N_guests
                sumY_yj = sumY_yj + gl%Y_yj_tf(j)
            end do

            !do J=1,N_guests
            !    Y_lang(J)= (moles_hdrt(J+1)) / (sum_moles_hdrt)
            !end do

            do J=1,gl%N_guests
                Y_lang(J)= (gl%Y_yj_tf(J)) / (sumY_yj)
            end do

            do i=1,gl%N_cavi
                do J=1,gl%N_guests
                    CiJ_mix_store(i,J)=CiJd(i,J)*CiJ(i,J)
                end do


                do J=1,gl%N_guests
                    CiJ_mix(i)=CiJ_mix(i)+( Y_lang(J)*CiJ_mix_store(i,J) )
                end do

                N_guests_double=dble(gl%N_guests)

                ! CiJ_mix(i)=(Cij_mix(i)) / N_guests_double
            end do

        else
            CiJ_mix(1)=CiJ(1,1)*CiJd(1,1)
            CiJ_mix(2)=CiJ(2,1)*CiJd(2,1)


        end if

    else
        continue
    end if

    !********************************************************************
    end subroutine hdrt_langmuirmix
    !********************************************************************
    !********************************************************************




    !******************************************************************************
    subroutine hdrt_int_pot_func(gl,temp,Ra,z,gue,int_wr)
    !******************************************************************************
    !
    !  Integral of the POTENTIAL FUNCTION of gas-water interaction w(r)
    !  June 2011; Vaclav May 2015



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    integer, parameter :: r_dimm = 250
    double precision, intent(in) :: temp
    double precision, dimension(6), intent(in):: Ra  ! [A] shell radius dependent on lattice param. a
    integer, dimension(6), intent(in):: z            ! [-] coordination number
    integer :: gue                                   ! guest index
    ! Output arguments
    double precision, intent(out) :: int_wr
    ! Local variables
    integer, dimension(4),parameter :: N = (/ 10, 11, 4, 5 /)
    double precision, dimension(4),parameter :: N_inv = 1d0/N
    double precision, dimension(4) :: delta
    integer :: shell, i, j, k
    double precision, dimension(r_dimm) :: r, wr, ex_wr_r2
    !double precision, dimension(6,r_dimm) :: w
    double precision :: w
    !double precision, dimension(6,r_dimm) :: w_kB_T

    !//////////////////////////////////////////////////////////////////
    double precision :: help_res,help_res2
    double precision, dimension(2) :: N_i_pow_minus_1,N_i_pow_minus_2
    double precision, dimension(2,4) :: N_i_pow_minus_all
    ! Cell potential function for given cavity 'wr'
    ! ---------------------------------------------
    double precision :: temp_inv
    double precision :: ra_k_inv,ra_k_inv_5,ra_k_inv_11
    double precision :: sigma_gue_6,sigma_gue_12
    do i=1,6
        if (z(i)/=0) then
            shell = i  ! number of water shells
        end if
    end do

    r(1) = 1.d-10
    do i = 1,r_dimm-1
        r(i+1) = r(i) + (Ra(1)-gl%a_hc(gue)-r(1))/(r_dimm-1)  ! [A] radial distance
    end do
    r(r_dimm) = Ra(1)-gl%a_hc(gue)-1.d-10

    !N = (/ 10, 11, 4, 5 /)     ! exponents for delta
    temp_inv = 1d0/temp

    sigma_gue_6 = gl%sigma(gue)**6
    sigma_gue_12 = sigma_gue_6**2

    wr = 0d0
    
    do k = 1,shell  ! k-loop over shells

        ra_k_inv= 1d0/Ra(k)
        ra_k_inv_5 = ra_k_inv**5
        ra_k_inv_11 = ra_k_inv_5*ra_k_inv_5*ra_k_inv

        do j = 1,r_dimm    ! j-radial coordination
            !delta(1:4) = 1.d0/N(1:4)*( (1.d0-r(j)/Ra(k)-gl%a_hc(gue)/Ra(k))**(-N(1:4)) - (1.d0+r(j)/Ra(k)-gl%a_hc(gue)/Ra(k))**(-N(1:4)) )  ! [-]
            help_res = r(j)*ra_k_inv
            help_res2 = gl%a_hc(gue)*ra_k_inv
            N_i_pow_minus_1(1) = 1d0/(1.d0-help_res-help_res2)
            N_i_pow_minus_1(2) = 1d0/(1.d0+help_res-help_res2)
            N_i_pow_minus_2(:) = N_i_pow_minus_1(:)**2
            N_i_pow_minus_all(:,3) = N_i_pow_minus_2(:)**2
            N_i_pow_minus_all(:,4) = N_i_pow_minus_all(:,3)*N_i_pow_minus_1(:)
            N_i_pow_minus_all(:,1) = N_i_pow_minus_all(:,4)**2
            N_i_pow_minus_all(:,2) = N_i_pow_minus_all(:,1) * N_i_pow_minus_1(:)

            delta = N_inv*(N_i_pow_minus_all(1,:)-N_i_pow_minus_all(2,:))

            w = sigma_gue_12*ra_k_inv_11/r(j)*(delta(1)+help_res2*delta(2)) - &
                sigma_gue_6*ra_k_inv_5/r(j)*(delta(3) + help_res2*delta(4)) ! [-] square bracket of potential function

            !w(k,j) = 2.d0*z(k)*gl%eps_k(gue)*w(k,j)  ! [K] potential function (eps = eps/kB [K])

            wr(j) = wr(j) + w*2.d0*z(k)*gl%eps_k(gue)*temp_inv
    !        w_kB_T(k,j) = w*2.d0*z(k)*gl%eps_k(gue)*temp_inv;  ! [-] exp term in the integrant

        end do  ! j-cycle
    end do  ! k-cycle

   ! w_kB_T(1:shell,:) = w(1:shell,:)*temp_inv;  ! [-] exp term in the integrant

    !do j = 1,r_dimm
    !    wr(j) = 0.d0
    !    do k = 1,shell
    !        wr(j) = wr(j) + w_kB_T(k,j) ! [-] sum of w_kB_T
    !    end do
    !end do
    !wr = sum(w_kB_T(1:shell,:),dim=1)


    ! Integrand ... exp(-wr/kB/T)*r^2
    ex_wr_r2(1) = 0.d0

    do i = 2,r_dimm
        ex_wr_r2(i) = dexp(-wr(i))*r(i)**2 ! [m2]
    end do ! i-cycle


    ! Trapezoidal numerical integration of potential function
    ! -------------------------------------------------------
    int_wr = 0.d0

    do i = 1,r_dimm-1
        int_wr = int_wr +  &
            (r(i+1)-r(i))*( ex_wr_r2(i+1) + ex_wr_r2(i) )/2.d0
    end do ! i-cycle

    int_wr = int_wr*1.d-30                                      !Für Volumenintegral; Maße von Meter zu Anström in drei Dimensionen

    end subroutine hdrt_int_pot_func
    !******************************************************************************
    !******************************************************************************



    !******************************************************************************
    subroutine hdrt_int_pot_double_func(gl,temp,Ra,z,gue,int_wr_double)
    !******************************************************************************
    !
    !  Integral of the POTENTIAL FUNCTION of gas-water interaction w(r) for second Langmuir constant
    !  June 2011; Vaclav May 2015; Benedikt April 2017



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    integer, parameter :: r_dimm = 250
    double precision, intent(in) :: temp
    double precision, dimension(6), intent(in):: Ra  ! [A] shell radius dependent on lattice param. a
    integer, dimension(6), intent(in):: z            ! [-] coordination number
    integer :: gue                                   ! guest index
    ! Output arguments
    double precision, intent(out) :: int_wr_double
    ! Local variables
    integer, dimension(4),parameter :: N = (/ 10, 11, 4, 5 /)
    double precision, dimension(4),parameter :: N_inv = 1d0/N
    double precision :: wr_gg
    double precision, dimension(4) :: delta_r1, delta, delta_r2, delta_gg
    integer :: shell, i, j, k
    double precision, dimension(r_dimm) :: r, wr, ex_wr_r2
!    double precision, dimension(6,r_dimm) ::  w_kB_T
    double precision :: w
    !//////////////////////////////////////////////////////////////////
    double precision :: help_res,help_res2
    double precision, dimension(2) :: N_i_pow_minus_1,N_i_pow_minus_2
    double precision, dimension(2,4) :: N_i_pow_minus_all
    double precision :: temp_inv
    double precision :: ra_k_inv,ra_k_inv_5,ra_k_inv_11
    double precision :: sigma_gue_6,sigma_gue_12

    temp_inv =  1d0/temp
    ! Cell potential function for given cavity 'wr'
    ! ---------------------------------------------
    w = 0
    do i=1,6
        if (z(i)/=0) then
            shell = i  ! number of water shells
        end if
    end do

    !Beginn Gästeschleife
    !do L=1,N_guests
    !    if (r_vdW(L)>1.d-14) then

    r(1) = 1.d-10
    !rd(1)=1.d-10
    do i = 1,r_dimm-1
        r(i+1) = r(i) + (Ra(1)-gl%a_hc(gue)-r(1)-gl%r_vdW(gue))/(r_dimm-1)  ! [A] radial distance
        !   rd(i+1)= rd(i) + (Ra(1)-a_hc(gued)-rd(1))/(r_dimm-1) !Integral for second occupand
    end do
    r(r_dimm) = Ra(1)-gl%a_hc(gue)-gl%r_vdW(gue)-1.d-10
    !  rd(r_dimm) = Ra(1)-a_hc(gued)-1.d-10

    !N =     ! exponents for delta
    !(/ 10, 11, 4, 5 /)
    sigma_gue_6 = gl%sigma(gue)**6
    sigma_gue_12 = sigma_gue_6**2

    wr = 0d0
    do k = 1,shell  ! k-loop over shells

        ra_k_inv= 1d0/Ra(k)
        ra_k_inv_5 = ra_k_inv**5
        ra_k_inv_11 = ra_k_inv_5*ra_k_inv_5*ra_k_inv

        do j = 1,r_dimm    ! j-radial coordination
            help_res = (r(j)-gl%r_vdW(gue))*ra_k_inv
            help_res2 = gl%a_hc(gue)*ra_k_inv
            N_i_pow_minus_1(1) = 1d0/(1.d0-help_res-help_res2)
            N_i_pow_minus_1(2) = 1d0/(1.d0+help_res-help_res2)
            N_i_pow_minus_2(:) = N_i_pow_minus_1(:)**2
            N_i_pow_minus_all(:,3) = N_i_pow_minus_2(:)**2
            N_i_pow_minus_all(:,4) = N_i_pow_minus_all(:,3)*N_i_pow_minus_1(:)
            N_i_pow_minus_all(:,1) = N_i_pow_minus_all(:,4)**2
            N_i_pow_minus_all(:,2) = N_i_pow_minus_all(:,1) * N_i_pow_minus_1(:)

            ! i-over delta N
            delta_r1 = N_inv*(N_i_pow_minus_all(1,:)-N_i_pow_minus_all(2,:))

            help_res = (r(j)+gl%r_vdW(gue))*ra_k_inv
            N_i_pow_minus_1(1) = 1d0/(1.d0-help_res-help_res2)
            N_i_pow_minus_1(2) = 1d0/(1.d0+help_res-help_res2)
            N_i_pow_minus_2(:) = N_i_pow_minus_1(:)**2
            N_i_pow_minus_all(:,3) = N_i_pow_minus_2(:)**2
            N_i_pow_minus_all(:,4) = N_i_pow_minus_all(:,3)*N_i_pow_minus_1(:)
            N_i_pow_minus_all(:,1) = N_i_pow_minus_all(:,4)**2
            N_i_pow_minus_all(:,2) = N_i_pow_minus_all(:,1) * N_i_pow_minus_1(:)

            delta_r2= N_inv*(N_i_pow_minus_all(1,:)-N_i_pow_minus_all(2,:))

            !do i = 1,4      ! i-over delta N
            !    delta_r1(i) = 1.d0/N(i)*( (1.d0-(r(j)-gl%r_vdW(gue))/Ra(k)-gl%a_hc(gue)/Ra(k))**(-N(i)) &
            !        &   - (1.d0+(r(j)-gl%r_vdW(gue))/Ra(k)-gl%a_hc(gue)/Ra(k))**(-N(i)) )  ! [-]     (r-r_vdW)
            !
            !    delta_r2(i) = 1.d0/N(i)*( (1.d0-(r(j)+gl%r_vdW(gue))/Ra(k)-gl%a_hc(gue)/Ra(k))**(-N(i)) &
            !        & - (1.d0+(r(j)+gl%r_vdW(gue))/Ra(k)-gl%a_hc(gue)/Ra(k))**(-N(i)) )  ! [-]     (r+r_vdW)
            !
            !    !delta_gg(i) = 1.d0/N(i)*( (1.d0-(abs(r(j)-r_vdW(gue)))/Ra(k)-a_hc(gued)/Ra(k))**(-N(i)) &       !Vermutlich nicht notwendig
            !    !- (1.d0+(abs(r(j)-r_vdW(gue)))/Ra(k)-a_hc(gued)/Ra(k))**(-N(i)) )
            !end do ! i-cycle


            w = (sigma_gue_12*ra_k_inv_11/(r(j)-gl%r_vdW(gue))*(delta_r1(1)+help_res2*delta_r1(2)) - &
                & sigma_gue_6*ra_k_inv_5/(r(j)-gl%r_vdW(gue))*(delta_r1(3) + help_res2*delta_r1(4))) + & ! [-] square bracket of potential function 2* because w(r1)+w(r2)
                & (sigma_gue_12*ra_k_inv_11/(r(j)+gl%r_vdW(gue))*(delta_r2(1)+help_res2*delta_r2(2)) - &  ! second occupand <--> lattice
                & sigma_gue_6*ra_k_inv_5/(r(j)+gl%r_vdW(gue))*(delta_r2(3) + help_res2*delta_r2(4)))








            !w(j,k) = 2.d0*z(k)*gl%eps_k(gue)*w(j,k)  ! [K] potential function (eps = eps/kB [K])
            !wd(j,k) = 2.d0*z(k)*eps_k(gued)*wd(j,k)

            wr(j) = wr(j) + w*2.d0*z(k)*gl%eps_k(gue)*temp_inv    ! [-] exp term in the integrant
!            w_kB_T(k,j) = w*2.d0*z(k)*gl%eps_k(gue)*temp_inv;    ! [-] exp term in the integrant

            !wd_kb_T(j,k) = wd(j,k)/temp


            !end if
            !Ende Gästeschleife

        end do  ! j-cycle
    end do  ! k-cycle

    !do j = 1,r_dimm
    !    wr(j) = 0.d0
    !    do k = 1,shell
    !        wr(j) = wr(j) + w_kB_T(k,j)! + wd_kb_T(j,k) ! [-] sum of w_kB_T         hier aufsummieren der einzelenen Interaktionen
    !    end do
    !end do

!    wr = sum(w_kB_T(1:shell,:),dim=1)

    ! Integrand ... exp(-wr/kB/T)*r^2
    ex_wr_r2(1) = 0.d0



    do i = 2,r_dimm
        ex_wr_r2(i) = dexp(-wr(i))*r(i)**2.d0 ! [m2]  +
    end do ! i-cycle


    ! Trapezoidal numerical integration of potential function
    ! -------------------------------------------------------
    int_wr_double = 0.d0

    do i = 1,r_dimm-1
        int_wr_double = int_wr_double +  &
            &  (r(i+1)-r(i))*( ex_wr_r2(i+1) + ex_wr_r2(i) )/2.d0
    end do ! i-cycle

    wr_gg= 4.d0*gl%epsd_k(gue)*(((gl%sigmad(gue))/(2.d0*gl%r_vdW(gue)-2*gl%a_hc(gue)))**12.d0 &     !Potentail funktion for guest-guest interaction,
        & -(((gl%sigmad(gue))/(2.d0*gl%r_vdW(gue)-2*gl%a_hc(gue)))**6.d0))     !spherical Kihara used

    wr_gg = dexp(-wr_gg*temp_inv)

    int_wr_double = int_wr_double*1.d-30

    int_wr_double = int_wr_double*wr_gg                                !Umstellen des Integrals ermöglicht vorziehen des e^w_gg Terms


    !Für Einheiten (Angström^3 zu Meter)


    end subroutine hdrt_int_pot_double_func
    !******************************************************************************
    !******************************************************************************



    !!******************************************************************************
    !subroutine hdrt_int_pot_func_double_int(temp,Ra,z,gue,int_wr)
    !!******************************************************************************
    !!
    !!  Integral of the POTENTIAL FUNCTION of gas-gas-water interaction w(r)
    !!  June 2011; Vaclav May 2015
    !!  Benedikt July 2017
    !

    !
    !!//////////////////////////////////////////////////////////////////
    !implicit none
    !
    !! Input arguments
    !integer, parameter :: r_dimm = 250
    !double precision, intent(in) :: temp
    !double precision, dimension(6), intent(in):: Ra  ! [A] shell radius dependent on lattice param. a
    !integer, dimension(6), intent(in):: z            ! [-] coordination number
    !integer :: gue                                   ! guest index
    !! Output arguments
    !double precision, intent(out) :: int_wr
    !! Local variables
    !integer, dimension(4) :: N
    !double precision, dimension(4) :: delta
    !integer :: shell, i, j, k
    !double precision, dimension(r_dimm) :: r, wr, ex_wr_r2, wr_gg
    !double precision, dimension(r_dimm,6) :: w, w_kB_T, w_gg, w_gg_kB_T
    !!//////////////////////////////////////////////////////////////////
    !
    !! Cell potential function for given cavity 'wr'
    !! ---------------------------------------------
    !
    !do i=1,6
    !    if (z(i)/=0) then
    !        shell = i  ! number of water shells
    !    end if
    !end do
    !
    !r(1) = 1.d-10
    !do i = 1,r_dimm-1
    !    r(i+1) = r(i) + (Ra(1)-a_hc(gue)-r(1))/(r_dimm-1)  ! [A] radial distance
    !end do
    !    r(r_dimm) = Ra(1)-a_hc(gue)-1.d-10
    !
    !N = (/ 10, 11, 4, 5 /)     ! exponents for delta
    !
    !do k = 1,shell  ! k-loop over shells
    !    do j = 1,r_dimm    ! j-radial coordination
    !        do i = 1,4      ! i-over delta N
    !        delta(i) = 1.d0/N(i)*( (1.d0-r(j)/Ra(k)-a_hc(gue)/Ra(k))**(-N(i)) &
    !                              - (1.d0+r(j)/Ra(k)-a_hc(gue)/Ra(k))**(-N(i)) )  ! [-]
    !        end do ! i-cycle
    !
    !w(j,k) = sigma(gue)**12/Ra(k)**11/r(j)*(delta(1)+a_hc(gue)/Ra(k)*delta(2)) - &
    !         sigma(gue)**6/Ra(k)**5/r(j)*(delta(3) + a_hc(gue)/Ra(k)*delta(4)) ! [-] square bracket of potential function
    !
    !w_gg(j,1) = 4.d0*epsd_k(gue)*((sigmad(gue)/((r(j)-2*a_hc(gue))))**12.d0 - (sigmad(gue)/((r(j)-2*a_hc(gue))))**6.d0)
    !
    !w(j,k) = 2.d0*z(k)*eps_k(gue)*w(j,k)  ! [K] potential function (eps = eps/kB [K])
    !
    !w_kB_T(j,k) = w(j,k)/temp;  ! [-] exp term in the integrant
    !
    !w_gg_kB_T(j,1) = w_gg(j,k)/temp
    !
    !    end do  ! j-cycle
    !end do  ! k-cycle
    !
    !do j = 1,r_dimm
    !    wr(j) = 0.d0
    !    wr_gg(j)=0.d0
    !    do k = 1,shell
    !        wr(j) = wr(j) + w_kB_T(j,k) ! [-] sum of w_kB_T
    !        wr_gg(j)  =wr_gg(j) + w_gg_kB_T(j,1)
    !    end do
    !end do
    !
    !! Integrand ... exp(-wr/kB/T)*r^2
    !ex_wr_r2(1) = 0.d0
    !
    !do i = 2,r_dimm
    !    ex_wr_r2(i) = dexp(-2.d0*wr(i)-wr_gg(i))*r(i)**2 ! [m2]
    !end do ! i-cycle
    !
    !
    !! Trapezoidal numerical integration of potential function
    !! -------------------------------------------------------
    !int_wr = 0.d0
    !
    !do i = 1,r_dimm-1
    !int_wr = int_wr +  &
    !        (r(i+1)-r(i))*( ex_wr_r2(i+1) + ex_wr_r2(i) )/2.d0
    !end do ! i-cycle
    !
    !int_wr = int_wr*1.d-30                                      !Für Volumenintegral; Maße von Meter zu Anström (10^-10^3)
    !
    !    end subroutine hdrt_int_pot_func_double_int
    !!******************************************************************************
    !!******************************************************************************


    !subroutine hdrt_int_pot_double_func2(temp,Ra,z,gue, int_wr_double)
    !!******************************************************************************
    !!******************************************************************************
    !!
    !!  Integral of the POTENTIAL FUNCTION of gas-water interaction w(r) for second Langmuir constant
    !!  Benedikt, May 2017
    !

    !
    !!//////////////////////////////////////////////////////////////////
    !implicit none
    !
    !! Input arguments
    !integer, parameter :: r_dimm = 250
    !double precision, intent(in) :: temp
    !double precision, dimension(6), intent(in):: Ra  ! [A] shell radius dependent on lattice param. a
    !integer, dimension(6), intent(in):: z            ! [-] coordination number
    !integer :: gue                                   ! guest index
    !! Output arguments
    !double precision, intent(out) :: int_wr_double
    !! Local variables
    !integer, dimension(4) :: N
    !!double precision ::
    !double precision, dimension(4) :: delta_r1, delta, delta_r2, delta_gg
    !integer :: shell, i, j, k
    !double precision, dimension(r_dimm) :: r, rg, wr, ex_wr_r2, ex_wr_gg_r2, int_pot_double, int_wr_gg, wr_gg!, int_wr_double
    !double precision, dimension(r_dimm,6) :: w, w_kB_T
    !!//////////////////////////////////////////////////////////////////
    !
    !! Cell potential function for given cavity 'wr'
    !! ---------------------------------------------
    !do i=1,6
    !    if (z(i)/=0) then
    !        shell = i  ! number of water shells
    !    end if
    !end do
    !
    !!r(1) = 1.d-10
    !!do i = 1,r_dimm-1
    !!    r(i+1) = r(i) + (Ra(1)-a_hc(gue)-r(1))/(r_dimm-1)  ! [A] radial distance
    !!end do
    !!    r(r_dimm) = Ra(1)-a_hc(gue)-1.d-10
    !
    !	r(1)=1.d-10
    !do i=1,r_dimm-1
    !
    !	r(i+1) = r(i) + ((Ra(1)-a_hc(gue)-r_vdW(gue)-r(1))/(r_dimm-1))
    !
    !end do
    !
    !    r(r_dimm) = Ra(1)-a_hc(gue)-1.d-10
    !
    !N = (/ 10, 11, 4, 5 /)     ! exponents for delta
    !
    !do k = 1,shell  ! k-loop over shells
    !    do j = 1,r_dimm    ! j-radial coordination
    !        do i = 1,4      ! i-over delta N
    !        delta_r1(i) = 1.d0/N(i)*( (1.d0-(r(j)-r_vdW(gue))/Ra(k)-a_hc(gue)/Ra(k))**(-N(i)) &
    !                           &   - (1.d0+(r(j)-r_vdW(gue))/Ra(k)-a_hc(gue)/Ra(k))**(-N(i)) )  ! [-]     (r-r_vdW)
    !
    !        delta_r2(i) = 1.d0/N(i)*( (1.d0-(r(j)+r_vdW(gue))/Ra(k)-a_hc(gue)/Ra(k))**(-N(i)) &
    !                             & - (1.d0+(r(j)+r_vdW(gue))/Ra(k)-a_hc(gue)/Ra(k))**(-N(i)) )  ! [-]     (r+r_vdW)
    !
    !
    !        end do ! i-cycle
    !
    !        w(j,k) = (sigma(gue)**12.d0/Ra(k)**11.d0/(r(j)-r_vdW(gue))*(delta_r1(1)+a_hc(gue)/Ra(k)*delta_r1(2)) - &
    !                & sigma(gue)**6.d0/Ra(k)**5.d0/(r(j)-r_vdW(gue))*(delta_r1(3) + a_hc(gue)/Ra(k)*delta_r1(4))) + & ! [-] square bracket of potential function 2* because w(r1)+w(r2)
    !               & (sigma(gue)**12.d0/Ra(k)**11.d0/(r(j)+r_vdW(gue))*(delta_r2(1)+a_hc(gue)/Ra(k)*delta_r2(2)) - &  ! second occupand <--> lattice
    !                & sigma(gue)**6.d0/Ra(k)**5.d0/(r(j)+r_vdW(gue))*(delta_r2(3) + a_hc(gue)/Ra(k)*delta_r2(4)))
    !
    !
    !
    !        w(j,k) = 2.d0*z(k)*eps_k(gue)*w(j,k)  ! [K] potential function (eps = eps/kB [K])
    !
    !
    !        w_kB_T(j,k) = w(j,k)/temp;    ! [-] exp term in the integrant
    !
    !    end do  ! j-cycle
    !end do  ! k-cycle
    !
    !do j = 1,r_dimm
    !    wr(j) = 0.d0
    !    do k = 1,shell
    !        wr(j) = wr(j) + w_kB_T(j,k) ! [-] sum of w_kB_T
    !    end do
    !end do
    !
    !! Integrand ... exp(-wr/kB/T)*r^2
    !ex_wr_r2(1) = 0.d0
    !
    !do i = 2,r_dimm
    !    ex_wr_r2(i) = dexp(-wr(i))*r(i)**2 ! [m^2]
    !end do ! i-cycle
    !
    !
    !
    !! Trapezoidal numerical integration of potential function
    !! -------------------------------------------------------
    !int_wr_double = 0.d0
    !
    !do i = 1,r_dimm-1
    !int_wr_double = int_wr_double +  &
    !      &  (r(i+1)-r(i))*( ex_wr_r2(i+1) + ex_wr_r2(i) )/2.d0
    !end do ! i-cycle
    !
    !!Integral for variable distance of guests in guestpot  - numerical integration
    !
    !
    !	rg(1)=2*a_hc(gue)+1.d-10
    !do i=1,r_dimm-1
    !
    !	rg(i+1) = rg(i) + ((Ra(1)-a_hc(gue)-rg(1))/(r_dimm-1))
    !
    !end do
    !
    !    rg(r_dimm) = Ra(1)-a_hc(gue)-1.d-10
    !
    !do i=1,r_dimm
    !
    !	wr_gg(i)=4.d0*epsd_k(gue)*(((sigmad(gue))/(rg(i)-2*a_hc(gue)))**12.d0 &     !Potentail funktion for guest-guest interaction,
    !                       & -((sigmad(gue))/(rg(i)-2*a_hc(gue)))**6.d0)
    !
    !	wr_gg(i)=wr_gg(i)/temp
    !
    !	ex_wr_gg_r2(i)=	dexp(-wr_gg(i))*rg(i)**2.d0
    !
    !end do
    !
    !int_wr_gg = 0.d0
    !
    !do i = 1,r_dimm-1
    !	int_wr_gg = int_wr_gg +  &
    !      &  (rg(i+1)-rg(i))*( ex_wr_gg_r2(i+1) + ex_wr_gg_r2(i) )/2.d0
    !end do
    !                                  !Umstellen des Integrals ermöglicht vorziehen des e^w_gg Terms
    !
    !
    !int_wr_gg = int_wr_gg*1.d-30
    !int_wr_double= int_wr_double*1d-30                                             !Für Einheiten (Angström^3 zu Meter)
    !
    !int_pot_double=int_wr_gg*int_wr_double
    !
    !!wr_out=(wr_gg/temp)+wr              ! nur für Zellpotentialbilder
    !!boltz=ex_wr_r2*ex_wr_gg_r2
    !!wr_out=dlog(wr_out)
    !!r_out=r
    !
    !    end subroutine hdrt_int_pot_double_func2
    !******************************************************************************


    !------------Chemical potential derivatives------------------------------------



    !******************************************************************************
    subroutine hdrt_Df_chem_potent_w(gl,temp,press,fug_g,Dchpw_Df)
    !******************************************************************************
    !
    !  Fugacity-derivative of the chemical potential
    !  of water in the hydrate phase [J/mol]
    !  June 2011; Vaclav May 2015   April 2017; Benedikt --> numerical derivation



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, dimension(30), intent(out) :: Dchpw_Df    ! dimension (N_hdrts x 1)
    ! Local variables
    double precision, dimension(3,30) :: CiJ
    double precision, dimension(3) :: CiJ_fug
    double precision, dimension(30) :: fug_locp, fug_locm, chempm, chempp, fug_loc
    double precision :: loc_chpw
    integer :: i, J, K
    !//////////////////////////////////////////////////////////////////
    if (gl%numderiv == .false.) then
        ! Langmuir constants @ T,p [MPa-1]
        call hdrt_Langmuir(gl,temp,press,CiJ)

        ! Derivative to fug_g(K)
        Dchpw_Df = 0.d0

        do K = 2,gl%N_hdrts    ! loop through all guests (water at K=1)

            do i = 1,gl%N_cavi
                CiJ_fug(i) = 0.d0
                do J = 1, gl%N_guests
                    CiJ_fug(i) = CiJ_fug(i) + CiJ(i,J)*fug_g(J)
                end do

                Dchpw_Df(K) = Dchpw_Df(K) + gl%v_cavi(i)*CiJ(i,K-1)/( 1.d0 + CiJ_fug(i) )
            end do

            Dchpw_Df(K) = - Rgas*temp*Dchpw_Df(K)

        end do

        ! Dchpw_Df is vector with dimension N_hdrts, i.e. (N_guests+1)
        ! first element is 0 - it corresponds to derivative to fugacity of water

    elseif (gl%numderiv == .true.) then

        !--------------------------------------------------------------------------------
        !  Numeric derivation
        chempp=0.d0
        chempm=0.d0
        fug_locp = fug_g*(1+1.d-6)
        fug_locm = fug_g*(1-1.d-6)
        do J = 1,gl%N_guests

            fug_loc = fug_g
            fug_loc(J) = fug_locp(J)

            call hdrt_chem_potent_w(gl,temp,press,fug_loc,loc_chpw)
            chempp(J) = loc_chpw
        end do

        do J = 1,gl%N_guests

            fug_loc = fug_g
            fug_loc(J) = fug_locm(J)

            call hdrt_chem_potent_w(gl,temp,press,fug_loc,loc_chpw)
            chempm (J) = loc_chpw

        end do

        Dchpw_Df(1) = 0.d0

        do J=1,gl%N_guests

            Dchpw_Df(J+1) = (chempp(J) - chempm(J))/(fug_locp(J) - fug_locm(J))   !"Durchreichen" der Ableitungen von N_guests auf N_hdrts

        end do
    endif

    end subroutine hdrt_Df_chem_potent_w
    !******************************************************************************
    !******************************************************************************



    !******************************************************************************
    subroutine hdrt_Dp_chem_potent_w(gl,temp,press,fug_g,Dchpw_Dp)
    !******************************************************************************
    !
    !  Pressure-derivative of the chemical potential
    !  of water in the hydrate phase [J/mol]
    !  June 2011; Vaclav May 2015



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, intent(out) :: Dchpw_Dp
    ! Local variables
    integer :: i
    double precision, dimension(2) :: chpw, p
    double precision :: gw_B, loc_chpw
    double precision, dimension(3,30) :: CiJ
    !//////////////////////////////////////////////////////////////////

    p(1) = press*(1.d0+1.d-6)
    p(2) = press*(1.d0-1.d-6)

    do i=1,2

        ! Chemical potential of water in hydrate phase [J/mol]
        call hdrt_chem_potent_w(gl,temp,p(i),fug_g,loc_chpw)
        chpw(i) = loc_chpw

    end do

    Dchpw_Dp = ( chpw(1)-chpw(2) )/(p(1)-p(2))

    end subroutine hdrt_Dp_chem_potent_w
    ! Pressure derivative could be solved analytically if CiJ = CiJ(T) and not CiJ(T,p)
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    subroutine hdrt_DT_chem_potent_w(gl,temp,press,fug_g,Dchpw_DT)
    !******************************************************************************
    !
    !  Temperature-derivative of the chemical potential
    !  of water in the hydrate phase [J/mol]
    !  June 2011; Vaclav May 2015



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, intent(out) :: Dchpw_DT
    ! Local variables
    integer :: i
    double precision, dimension(2) :: chpw, T
    double precision :: gw_B, loc_chpw
    double precision, dimension(3,30) :: CiJ
    !//////////////////////////////////////////////////////////////////

    T(1) = temp*(1.d0+1.d-4)
    T(2) = temp*(1.d0-1.d-4)

    do i=1,2

        ! Chemical potential of water in hydrate phase [J/mol]
        call hdrt_chem_potent_w(gl,T(i),press,fug_g,loc_chpw)
        chpw(i) = loc_chpw

    end do

    Dchpw_DT = ( chpw(1)-chpw(2) )/(T(1)-T(2))

    end subroutine hdrt_DT_chem_potent_w
    !******************************************************************************
    !******************************************************************************





    !******************************************************************************
    subroutine hdrt_Dxi_Dfj(gl,temp,press,fug_g,DoccupiK_Df,CiJ,DxJ_Df)
    !******************************************************************************
    !
    !  Derivative of molar composition in the hydrate phase [Pa-1]
    !  water is the first component
    !  Vaclav October 2015
    !  Set to numerical derivation: Benedikt, April 2017



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, dimension(3,30), intent(out) :: CiJ
    double precision, dimension(3,30,30), intent(out) :: DoccupiK_Df
    !In the following matrix the derivatives of the hydrate composition with respect to fugacity is stored
    double precision, dimension(30,30), intent(out):: DxJ_Df
    ! Local variables
    double precision, dimension(3,30):: occup, locp_occup, locm_occup, occup_p, occup_m, occup_single, occup_double
    double precision, dimension(3):: CiK_fK
    double precision, dimension(30):: sum_iK, v_occup, loc_xH, xHm, xHp, fug_locm, fug_locp, fug_loc, num_diff!, occup_ls, occup_ld, occup_sms, occup_smd
    double precision:: v_occup_all
    integer:: i, J, K, M, xJ, fM
    double precision, dimension(30,30) :: Dxguest_Df, sum_i
    !//////////////////////////////////////////////////////////////////

    CiK_fK = 0.d0           ! sum_K(C_iK*f_K)
    v_occup = 0.d0          ! sum_i(v_i*occup_iJ)
    v_occup_all = 0.d0      ! sum_i(v_i*sum_J(occup_iJ))
    DoccupiK_Df = 0.d0      ! D(occup_iK)/D(f_M)
    sum_i = 0.d0            ! sum_i( v_i * D(occup_iK)/D(f_M) )
    sum_iK = 0.d0           ! sum_i( v_i * sum_K(D(occup_iK)/D(f_M)) )

    if (gl%numderiv == .false.) then
        !
        ! Cage occupancies for all components @ T,p [-]
        call hdrt_occupancy(gl,temp,press,fug_g,occup,CiJ,occup_single, occup_double)


        ! Derivative of cage occupancy D(occup_iK)/+D(f_M)
        do i = 1,gl%N_cavi
            CiK_fK(i) = 0.d0
            do K = 1,gl%N_guests
                CiK_fK(i) = CiK_fK(i) + CiJ(i,K)*fug_g(K)
                v_occup_all = v_occup_all + gl%v_cavi(i)*occup(i,K)
            end do

            do K = 1,gl%N_guests
                do M = 1,gl%N_guests
                    if (M.eq.K) then
                        DoccupiK_Df(i,K,M) = CiJ(i,M)*(1.d0+CiK_fK(i)-CiJ(i,M)*fug_g(M)) / (1+CiK_fK(i))**2.d0
                    else
                        DoccupiK_Df(i,K,M) = -CiJ(i,M)*CiJ(i,K)*fug_g(K) / (1+CiK_fK(i))**2.d0
                    end if
                end do
            end do

        end do

        ! Derivative of hydrate phase mole fractions of GUESTS according to the fugacities of the guests
        do J = 1,gl%N_guests
            do M = 1,gl%N_guests
                sum_i(J,M) = 0.d0
                sum_iK(M) = 0.d0
                do K = 1,gl%N_guests
                    sum_i(J,M) = 0.d0
                    v_occup(K) = 0.d0
                    do i = 1,gl%N_cavi
                        sum_iK(M) = sum_iK(M) + gl%v_cavi(i)*DoccupiK_Df(i,K,M)
                        sum_i(J,M) = sum_i(J,M) + gl%v_cavi(i)*DoccupiK_Df(i,J,M)

                        v_occup(K) = v_occup(K) + gl%v_cavi(i)*occup(i,K)
                    end do
                end do

                Dxguest_Df(J,M) = ( sum_i(J,M)*(1.d0+v_occup_all)-v_occup(J)*sum_iK(M) ) / (1.d0+v_occup_all)**2.d0

            end do
        end do
        ! In Dxguest_Df, rows (J) stand for mole fraction of guest J and columns stand for fugacity derivative (M)
        ! Dxguest_Df needs to be transposed according to DxJ_Df !!!


        ! Derivative of hydrate phase mole fractions according to the fugacities of the guests (all comps., i.e. also WATER)
        ! Rows stand for fugacities and columns for the mole fractions
        ! x1 = mole fraction of water, f1 = fugacity fo water ... for all components in hydrate phase: D(xJ)/D(f1) = 0
        DxJ_Df = 0.D0   !The derivatives of the composition wrt the fugacity of H2O in the mixture are 0
        ! so: DxJ_Df(1,:) = 0.D0

        do xJ = 1,gl%N_guests
            do fM = 1,gl%N_guests
                DxJ_Df(fM+1,xJ+1) = Dxguest_Df(xJ,fM)     ! guest components only - Dxguest_Df is transposed
            end do
        end do

        do fM = 1,gl%N_guests
            DxJ_Df(fM+1,1) = -sum_iK(fM) / (1.d0+v_occup_all)**2.d0    ! water as the first component
        end do

        ! DxJ_Df is matrix with dimension N_hdrts x N_hdrts
        ! first line is 0 and corresponds to derivative according to water fugacity
        ! elements in first column are derivatives of composition of water D(x1)/D(f_M)

    elseif (gl%numderiv == .true.) then
        !----------------------------------------------------------------------------------------------------------------
        !Numeric derivative

        fug_locp = fug_g*(1+1.d-4)
        fug_locm = fug_g*(1-1.d-4)
        xHp=0.d0
        xHm=0.d0
        DxJ_Df=0.d0


        do J = 1,gl%N_guests

            fug_loc = fug_g
            fug_loc(J) = fug_locp(J)

            call hdrt_mole_fract(gl,temp,press,fug_loc,occup,CiJ,loc_xH,occup_single, occup_double)    !hdrt_mole_fract(temp,press,fug_g,occup,CiJ,xH)
            xHp  = loc_xH
            occup_p  = occup

            fug_loc = fug_g
            fug_loc(J) = fug_locm(J)

            call hdrt_mole_fract(gl,temp,press,fug_loc,occup,CiJ,loc_xH,occup_single, occup_double)
            xHm  = loc_xH
            occup_m  = occup

            num_diff =  (xHp - xHm)/(fug_locp - fug_locm)

            do i = 1,gl%N_hdrts

                DxJ_Df(J+1,i) = (xHp(i) - xHm(i))/(fug_locp(J) - fug_locm(J))
                do k = 1, gl%N_cavi
                    DoccupiK_Df(k,i,J) = (occup_p(k,i) - occup_m(k,i))/(fug_locp(J) - fug_locm(J))
                enddo
            end do

        end do



    endif



    end subroutine hdrt_Dxi_Dfj
    !******************************************************************************
    !******************************************************************************

    !******************************************************************************
    subroutine hdrt_Dxi_DT(gl,temp,press,fug_g,Dxi_DT)
    !******************************************************************************
    !
    !  Numerical Temperature-derivative of the composition
    !  of water in the hydrate phase
    !  March 2017, Sebastian and Andreas



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, dimension(30), intent(out) :: Dxi_DT
    ! Local variables
    integer :: i
    double precision, dimension(2) :: T
    double precision, dimension(2,30) :: xH
    double precision, dimension(30) :: loc_xH!, occup_ls, occup_ld, occup_sms, occup_smd
    double precision, dimension(3,30) :: CiJ, occup,occup_single, occup_double
    !//////////////////////////////////////////////////////////////////

    T(1) = temp*(1.d0+1.d-4)
    T(2) = temp*(1.d0-1.d-4)

    do i=1,2

        ! Chemical potential of water in hydrate phase [J/mol]
        call hdrt_mole_fract(gl,T(i),press,fug_g,occup,CiJ,loc_xH, occup_single, occup_double)
        xH(i,:) = loc_xH

    end do

    Dxi_DT = ( xH(1,:)-xH(2,:) )/(T(1)-T(2))

    end subroutine hdrt_Dxi_DT
    !******************************************************************************
    !******************************************************************************



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !Andreas January 2014
    !New Routines for the calculation of the enthalpy and entropy of hydrate
    ! December 2015
    ! Vaclav's try for updating the subroutines for mixed hydrates ... Not sure whether it will work ;-)
    ! March 2017
    ! Update of the subroutine 'hdrt_entropy_w_B' for new mixing rule a_M0(T0,p0) = sum( a0J(T0,p0)*YJ(T,fug_g) )
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !******************************************************************************
    subroutine hdrt_entropy_w_B(gl,temp,press,fug_g, sw_B, error)
    !******************************************************************************
    !
    !  Entropy of water in reference empty hydrate lattice [J/molK]
    !  January 2014, Andreas; February 2016, Vaclav - update for mixtures
    !  March 2017, Vaclav - new mixing rule a_M0(T0,p0) = sum( a0J(T0,p0)*YJ(T,fug_g) )
    !  July 2017, Benedikt - adapted for double occupancy



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, intent(out) :: sw_B
    ! Local variables
    double precision :: part_hw, part_vol, tempp, tempm !, part_vol_num
    integer:: error
    ! Mixed hydrates
    integer:: check_B0_const
    integer:: J, i
    double precision :: a_MT0, Da_MT03_DT, press_a_MT0, Da_MT03_DT_num, sw_B_num
    double precision, dimension(3,30) :: CiJ, DCiJ_DT    ! Langmuir constant and its derivative wrt temperature T
    double precision, dimension(3,30) :: occup, Docc_DT,occup_single, occup_double, occupm, occupp  ! cage occupany and its derivative wrt temperature T
    double precision, dimension(30) :: DY_DT, Ym, Yp             ! guest composition derivative wrt temperature T
    ! Sums of the Langmuir constant and fugacity
    double precision, dimension(3) :: sum_CiJ_fugJ, sum_DCiJ_DT_fugJ
    ! Sums of the coordination number and cage occupancy
    double precision, dimension(30) :: sum_vi_DocciJ_DT,  sum_vi_occiJ!, occup_ls, occup_ld, occup_sms, occup_smd
    double precision :: sum_vi_DocciJ_DT_all,  sum_vi_occiJ_all
    ! Numerical derivation of the volume
    double precision, dimension(2) :: T, int_vol
    !//////////////////////////////////////////////////////////////////

    ! Contribution of the enthalpy term
    part_hw =  gl%hw_B0 / gl%T0 + gl%cpH_a*(temp-gl%T0) + gl%cpH_b*dlog(temp/gl%T0)   ! 2nd term in Eq. 4.58 in Andreas' diss.

    ! Contribution of the volume term
    if (gl%Latt_par == 'vj') then
        ! Not yet implemented
        error = -11111

    elseif (gl%Latt_par == 'bs') then
        ! Not yet implemented
        error = -11111

    elseif (gl%Latt_par == 'm4') then

        ! Check if B0_hdrt is constant - Mixed hydrates 2015
        check_B0_const = 1;
        do J = 1,gl%N_guests
            if (dabs(gl%B0_hdrt(1,J)-10.d9).gt.1.d3) then
                check_B0_const = 0  ! B0_hdrt is not constant for all hydrate formers
            end if
        end do

        if (check_B0_const.eq.1) then
            ! Lattice parameter of mixed hydrate a0_hdrt @ T0 and p0 from Y(T,fug)
            call hdrt_latt_par_mix_0_Y_Tf(gl,temp,press,fug_g,a_MT0)

            ! Cage occupancies & Langmuir constant @ temp & fug_g
            call hdrt_occupancy(gl,temp,press,fug_g,occup,CiJ,occup_single, occup_double)
            ! Temperature derivatice of the Langmuir constant
            call hdrt_dC_dT(gl,temp,press, DCiJ_DT, error)         ! (numerically)

            ! Temperature derivative of the cage occupancies Docc_DT(i,J)
            ! -----------------------------------------------------------
            !sum_CiJ_fugJ = 0.d0                                            ! Beginn analytische Ableitung
            !sum_DCiJ_DT_fugJ = 0.d0
            !
            !do i = 1,N_cavi
            !    do J = 1,N_guests
            !        sum_CiJ_fugJ(i) = sum_CiJ_fugJ(i) + CiJ(i,J)*fug_g(J)
            !        sum_DCiJ_DT_fugJ(i) = sum_DCiJ_DT_fugJ(i) + DCiJ_DT(i,J)*fug_g(J)
            !    end do
            !end do
            !
            !do i = 1,N_cavi
            !    do J = 1,N_guests
            !        Docc_DT(i,J) = (1.d0+sum_CiJ_fugJ(i))*fug_g(J)*DCiJ_DT(i,J) - CiJ(i,J)*fug_g(J)*sum_DCiJ_DT_fugJ(i)
            !        Docc_DT(i,J) = Docc_DT(i,J) / (1.d0+sum_CiJ_fugJ(i))**2
            !    end do
            !end do                                 ! Ende analytische Ableitung

            !numerically considering double occupation  Benedikt July 2017			Vermutlich aufgrund der numerischen Ableitung von DY_DT nicht notwendig.

            tempp = temp *(1.d0 + 1.d-4)

            call hdrt_occupancy(gl,tempp,press,fug_g,occupp,CiJ,occup_single, occup_double)

            tempm = temp * (1.d0 - 1.d-4)

            call hdrt_occupancy(gl,tempm,press,fug_g,occupm,CiJ,occup_single, occup_double)

            if (gl%N_guests>1.d0) then
                continue
            end if

            do i=1,gl%N_cavi
                do J=1,gl%N_guests
                    Docc_DT(i,J) = (occupp(i,J) - occupm(i,J))/ (tempp-tempm)
                end do
            end do

            ! ! Temperature derivative of guest composition DY_DT(J)
            ! ! -----------------------------------------------------
            ! sum_vi_occiJ = 0.d0
            ! sum_vi_DocciJ_DT = 0.d0
            ! sum_vi_occiJ_all = 0.d0
            ! sum_vi_DocciJ_DT_all = 0.d0

            ! do J = 1,N_guests
            ! do i = 1,N_cavi
            ! sum_vi_occiJ(J) = sum_vi_occiJ(J) + v_cavi(i)*occup(i,J)
            ! sum_vi_DocciJ_DT(J) = sum_vi_DocciJ_DT(J) + v_cavi(i)*Docc_DT(i,J)

            ! sum_vi_occiJ_all = sum_vi_occiJ_all + v_cavi(i)*occup(i,J)
            ! sum_vi_DocciJ_DT_all = sum_vi_DocciJ_DT_all + v_cavi(i)*Docc_DT(i,J)
            ! ! in moles_hdrt water is the 1st component
            ! end do
            ! end do

            ! do J = 1,N_guests
            ! DY_DT(J) = sum_vi_DocciJ_DT(J)*sum_vi_occiJ_all - sum_vi_occiJ(J)*sum_vi_DocciJ_DT_all
            ! DY_DT(J) = DY_DT(J) / sum_vi_occiJ_all**2
            ! end do

            call hdrt_YJ_Tf(gl,tempp, press, fug_g, Yp)

            call hdrt_YJ_Tf(gl,tempm, press, fug_g, Ym)

            DY_DT = (Yp-Ym) / (tempp-tempm)

            ! Temperature derivative of lattice parameter**3 d(a_MT0**3)/dT
            ! -------------------------------------------------------------
            Da_MT03_DT = 0.d0
            if (gl%Latt_par_mix == 'volume') then
                do J = 1,gl%N_guests
                    Da_MT03_DT = Da_MT03_DT + gl%a0_hdrt(J)**3 * DY_DT(J)
                end do
            elseif (gl%Latt_par_mix == 'vegard') then
                do J = 1,gl%N_guests
                    Da_MT03_DT = Da_MT03_DT + gl%a0_hdrt(J) * DY_DT(J)
                end do
                Da_MT03_DT = 3.d0*a_MT0**2*Da_MT03_DT
            end if

            !! Numerical check of Da_MT03_DT - to be commented later - OK - checked February 2016 & March 2017
            !T(1) = temp + 1.d-2
            !T(2) = temp - 1.d-2
            !    call hdrt_latt_par_mix_0_Y_Tf(T(1),press,fug_g,int_vol(1))
            !    call hdrt_latt_par_mix_0_Y_Tf(T(2),press,fug_g,int_vol(2))
            !Da_MT03_DT_num = (int_vol(1)**3-int_vol(2)**3)/(T(1)-T(2))


            ! Analytical integration of the molar volume vol_B - correlation based on Murnaghan EoS with dB0 = 4 = const.
            !Pressure dependent part
            part_vol = ( ((gl%B0_m(1)+gl%B0_m(2)*press)/(gl%B0_m(1)+gl%B0_m(2)*gl%p0_a))**((gl%B0_m(2)-1)/gl%B0_m(2)) - &
                ((gl%B0_m(1)+gl%B0_m(2)*gl%p0)/(gl%B0_m(1)+gl%B0_m(2)*gl%p0_a))**((gl%B0_m(2)-1)/gl%B0_m(2)) ) * (gl%B0_m(1)+gl%B0_m(2)*gl%p0_a)/(gl%B0_m(2)-1.D0)
            !Temperature dependent part
            part_vol = part_vol * Nav/gl%Nw*dexp( 3.d0*gl%alfa(1)*(temp-gl%T0_a)  &
                + 3.d0*gl%alfa(2)*(temp-gl%T0_a)**2 + 3.d0*gl%alfa(3)*(temp-gl%T0_a)**3 ) * &
                ( a_MT0**3*(3.d0*gl%alfa(1)+ 6.d0*gl%alfa(2)*(temp-gl%T0_a) + 9.d0*gl%alfa(3)*(temp-gl%T0_a)**2) + Da_MT03_DT )*1.d-30

            ! Analytical derivative of integral of volume checked numerically for mixed hydrates - OK - February 2016

        else    ! non-constant bulk modulus or other definition of pressure dependence
            !DEC$ IF .NOT. DEFINED(WO_WRITE_TO_CONSOLE)

            write(*,*) 'Warning - derivative of volume integral in hdrt_entropy_w_B calculated numerically'
            write(*,*) '        - bulk modulus is not the same for all involved hydrates'
            !DEC$ END IF ! WO_WRITE_TO_CONSOLE
            ! Numerical integration of the molar volume vol_B
            ! and subsequent numerical derivation wrt temp
            T(1) = temp + 1.d-3
            T(2) = temp - 1.d-3

            call hdrt_int_volB_num(gl,T(1),press,fug_g,int_vol(1)) ! Integral(vol_B)
            call hdrt_int_volB_num(gl,T(2),press,fug_g,int_vol(2)) ! Integral(vol_B)

            part_vol = (int_vol(1)-int_vol(2))/(T(1)-T(2))

        end if

    end if

    ! Entropy of water in empty lattice [J/molK]
    ! Note that this entropy is the molar entropy of the empty beta lattice and thus divided by the
    ! molar amount OF WATER MOLECULES.
    sw_B = -gl%gw_B0/gl%T0 + part_hw - part_vol

    !!Numerical check - OK - February 2016 & March 2017
    !call hdrt_entropy_w_B_num(temp,press,fug_g, sw_B_num, error)

    end subroutine hdrt_entropy_w_B
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    subroutine hdrt_entropy(gl,temp, press, x_hyd, fug_g, s, error)
    !******************************************************************************
    !
    !  Entropy of hydrate [J/molK]
    !  January 2014, Andreas
    !  December 2015, Vaclav - update for mixed hydrates
    !  July 2017, Benedikt - update for double occupied hydrates





    !//////////////////////////////////////////////////////////////////



    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision :: temp, press, sum !, FUGCOPURE_CALC
    double precision, dimension(30) :: fug_g, fug_g0
    double precision :: dchem_pot_dT_fluid, dchem_pot_dx_fluid, dfug_dT_fluid, dfug_dx_fluid
    ! Output arguments
    double precision :: s
    ! Local variables
    double precision :: s_empty_lattice, s_part_guest, Z
    !double precision, dimension(2)  :: sum_CJifJ, sum_dCJidfJ
    !double precision, dimension(30,2) :: CiJ    !Langmuir constant C(guest molecule, cavity type)
    !double precision, dimension(30,2) :: d_CiJ_d_T    !Derivative of Langmuir constant C(guest molecule, cavity type) wrt the temperatur T
    integer:: error, i, J, n

    double precision :: part1, part2, dxw_H_df, dxw_H_dT, dxw_Fluid_dT_constxH, dxw_fluid_dT, dxw_fluid_dxj
    double precision, dimension(2) :: dchem_pot_dT
    double precision, dimension(30) :: x_hyd, d_fug_g_dT
    !double precision, dimension(2) :: sum_CJi_dfJdT, sum_CJi_dfJdxj

    double precision :: dChemPot_wH_dT  !Partial derivative of the chemical potential of water in hydrate wrt temperature at constant p and fJ
    double precision :: p0_guest    ! ,s0_guest
    double precision ::  rho0_guest, rho0_est, rhoredmixorg, tredmixorg !,  G_CALC, S_CALC, SR_CALC,
    integer:: mixflag, iphase

    ! Mixed hydrates 2015
    double precision, dimension(3,30) :: CiJ, CiJd    !Langmuir constant C(guest molecule, cavity type)
    double precision, dimension(3)  :: sum_CJifJ, sum_dCJidfJ, sum_dCijdouble_dfJ, sum_CJidouble_fJ, sum_dC_dT_double
    double precision, dimension(3,30) :: d_CiJ_d_T, d_CiJd_d_T    !Derivative of Langmuir constant C(guest molecule, cavity type) wrt temperature T
    double precision, dimension(3) :: sum_CJi_dfJdT, sum_CJi_dfJdxj
    double precision, dimension(30) :: s0_guest

    ! For Verification in doube occ test
    double precision :: cphwp, cphwm,s_empty_lattice_div, dChemPot_wH_dT_num, sum_s_guests, tempp, tempm

    !//////////////////////////////////////////////////////////////////


    !Calculate the contribution of the empty lattice to the entropy
    !call hdrt_entropy_w_B(temp,press,s_empty_lattice, error)    ! numerically checked for mixtures - OK ... Eq. 4.58 in Andreas' diss.
    call hdrt_entropy_w_B(gl,temp,press,fug_g,s_empty_lattice, error)    ! fugacity also as an input ... March 2017


    !Calculate the contribution of the guest molecule(s) to the entropy

    ! Pure hydrates: call hdrt_Langmuir(temp,press,CiJ(1,:))
    ! Mixed hydrates 2015 ... CiJ(3,30)
    call hdrt_Langmuir(gl,temp,press,CiJ)
    call hdrt_Langmuir_double(gl,temp,press,CiJd)
    ! Eventually hdrt_mole_fract(temp,press,fug_g,occup,CiJ,x_hyd) can be called here in order to x_hyd

    !Calculate the sum CJi * fJ over all guest molecules J for each cavity type i
    sum_CJifJ = 0.D0
    sum_CJidouble_fJ = 0.d0
    do i = 1, gl%N_cavi
        do J = 1, gl%nrofhydrateformers-1
            ! Pure hydrates: sum_CJifJ(i) = sum_CJifJ(i) + CiJ(j,i) *fug_g(j)
            ! Mixed hydrates 2015
            sum_CJifJ(i) = sum_CJifJ(i) + CiJ(i,J) *fug_g(J)
            !Double occupancy 2017:
            sum_CJidouble_fJ(i) = sum_CJidouble_fJ(i) + CiJ(i,J)*CiJd(i,J)*fug_g(J)**2.d0
        end do
    end do

    call hdrt_dC_dT(gl,temp,press, d_CiJ_d_T, error)
    call hdrt_dCd_dT(gl,temp, press, d_CiJd_d_T)
    !Calculate the sum (d_CJi_d_T * fJ) over all guest molecules J for each cavity type i
    sum_dCJidfJ = 0.D0
    sum_dCiJdouble_dfJ = 0.d0
    do i = 1, gl%N_cavi
        do J = 1, gl%nrofhydrateformers-1
            ! Pure hydrates: sum_dCJidfJ(i) = sum_dCJidfJ(i) + d_CiJ_d_T(j,i) *fug_g(j)  ! beware of flipped indexes by d_CiJ_d_T(guest,cavity)
            ! Mixed hydrates 2015:
            sum_dCJidfJ(i) = sum_dCJidfJ(i) + d_CiJ_d_T(i,J) *fug_g(J)
            !Double occupancy 2017:
            sum_dCiJdouble_dfJ(i) = sum_dCiJdouble_dfJ(i) + d_CiJd_d_T(i,J)*fug_g(J)**2.d0
        end do
    end do

    !further loops for double occupancy:
    !loops for double occupation derivates:
    sum_dC_dT_double = 0.d0
    do i=1,gl%N_cavi
        do J=1, gl%nrofhydrateformers-1
            sum_dC_dT_double(i) = sum_dC_dT_double(i) + d_CiJ_d_T(i,J)*CiJd(i,J)*fug_g(J)**2.d0 + d_CiJd_d_T(i,J)*CiJ(i,J)*fug_g(J)**2.d0
        end do
    end do
    !Calculate the contribution of the guest molecules to the entropy
    s_part_guest = 0.D0
    do i = 1, gl%N_cavi
        s_part_guest = s_part_guest + Temp * gl%v_cavi(i) / (1.D0 + sum_CJifJ(i) + sum_CJidouble_fJ(i)) * (sum_dCJidfJ(i) + sum_dC_dT_double(i) )&  !<-- fehlt nach temp*v_cavi(i) / nicht noch was?, nein, ich glaube nicht!
            &+ gl%v_cavi(i) * dlog(1.D0 + sum_CJifJ(i) + sum_CJidouble_fJ(i)) ! 2 last terms in Eq. 4.57 in Andreas' diss. , addaped for double occupancy
    end do

    s_part_guest = s_part_guest * Rgas


    dChemPot_wH_dT = - (s_empty_lattice + s_part_guest)     ! Eq. 4.57 in Andreas' diss.

    !! Pure hydrates
    !!Add the guest dependent part:
    !!Get the hydration number
    !Z = 0.D0
    !do J = 1, nrofhydrateformers-1
    !    do i = 1, N_cavi
    !        Z = Z + v_cavi(i)*CiJ(i,J) *fug_g(J)/(1.D0+sum_CJifJ(i))
    !    end do
    !end do
    !Z = 1.D0 / Z
    !!hyd_nr = 1.D0 / (v_cavi(1) * occup(1) + v_cavi(2) * occup(2))
    !!Z = Z + 1.D0

    !Get the entropy of the pure guest at the reference pressure p0_guest
    !ORIGINAL START
    p0_guest = 1.D-6  !p0 in MPa ... Eq. 4.55 to 4.56 in Andreas' diss.

    !! Pure hydrates
    !mixflag = 2 !Guest at second position in fluid vector
    !iphase = 0  !Most likely vapor, but also check for liquid, you never know....
    !rho0_est = 0.D0
    !
    !rhoredmixorg = rhoredmix
    !tredmixorg = tredmix
    !rhoredmix = rhored(2)
    !tredmix = tred(2)
    !rho0_guest = rhomix_calc(gl,temp, p0_guest, rho0_est, iphase, mixflag)
    !rho0_guest = p0_guest*1.D6/Rgas/Temp
    !!Ideal gas entropy
    !s0_guest = S_CALC(gl,temp,rho0_guest, mixflag) - SR_CALC(gl,temp,rho0_guest, mixflag)
    !rhoredmix = rhoredmix
    !tredmix = tredmixorg
    !!116.945500422556

    rhoredmixorg = gl%rhoredmix
    tredmixorg = gl%tredmix

    ! Mixed hydrates 2015
    do J = 1,gl%N_guests
        mixflag = J+1 !Guest in the fluid vector - hydrate formers need to be right after water in the fluid list
        iphase = 0    !Most likely vapor, but also check for liquid, you never know....
        rho0_est = 0.D0

        gl%rhoredmix = gl%rhored(J+1)
        gl%tredmix = gl%tred(J+1)
        !rho0_guest = rhomix_calc(gl,temp, p0_guest, rho0_est, iphase, mixflag)
        rho0_guest = p0_guest*1.D6/Rgas/Temp
        !Ideal gas entropy
        s0_guest(J) = S_CALC(gl,temp,rho0_guest, mixflag) - SR_CALC(gl,temp,rho0_guest, mixflag)

        gl%rhoredmix = rhoredmixorg
        gl%tredmix = tredmixorg
    end do


    ! Pure hydrates
    ! s =  x_hyd(1) * (-dChemPot_wH_dT + 1.D0 / Z * (s0_guest - Rgas*dlog(fug_g(1)/(p0_guest*1.D6))) )

    ! Mixed hydrates 2015
    s =  - x_hyd(1)*dChemPot_wH_dT
    do J = 1,gl%N_guests
        s =  s + x_hyd(J+1)*( s0_guest(J) - Rgas*dlog(fug_g(J)/(p0_guest*1.D6)) )    ! Eq. 4.55 in Andreas' diss.
    end do
    !ORIGINAL END

    !!TEST START
    !p0_guest = 1.D-8 - 1.d-7  !p0 in MPa ... Eq. 4.55 to 4.56 in Andreas' diss.
    !open(unit = 999, file='Z:\Semrau\HydMix\test.txt')
    !write(999,*)'p0_guest,sum,entropy'
    !do i = 1,10000
    !p0_guest = p0_guest + 1.d-7  !p0 in MPa ... Eq. 4.55 to 4.56 in Andreas' diss.
    !
    !        !! Pure hydrates
    !        !mixflag = 2 !Guest at second position in fluid vector
    !        !iphase = 0  !Most likely vapor, but also check for liquid, you never know....
    !        !rho0_est = 0.D0
    !        !
    !        !rhoredmixorg = rhoredmix
    !        !tredmixorg = tredmix
    !        !rhoredmix = rhored(2)
    !        !tredmix = tred(2)
    !        !rho0_guest = rhomix_calc(gl,temp, p0_guest, rho0_est, iphase, mixflag)
    !        !rho0_guest = p0_guest*1.D6/Rgas/Temp
    !        !!Ideal gas entropy
    !        !s0_guest = S_CALC(gl,temp,rho0_guest, mixflag) - SR_CALC(gl,temp,rho0_guest, mixflag)
    !        !rhoredmix = rhoredmix
    !        !tredmix = tredmixorg
    !        !!116.945500422556
    !
    !rhoredmixorg = rhoredmix
    !tredmixorg = tredmix
    !
    !! Mixed hydrates 2015
    !do J = 1,N_guests
    !    mixflag = J+1 !Guest in the fluid vector - hydrate formers need to be right after water in the fluid list
    !    iphase = 0    !Most likely vapor, but also check for liquid, you never know....
    !    rho0_est = 0.D0
    !
    !    rhoredmix = rhored(J+1)
    !    tredmix = tred(J+1)
    !    rho0_guest = rhomix_calc(gl,temp, p0_guest, rho0_est, iphase, mixflag)
    !    !rho0_guest = p0_guest*1.D6/Rgas/Temp
    !    !Ideal gas entropy
    !    !s0_guest(J) = S_CALC(gl,temp,rho0_guest, mixflag) - SR_CALC(gl,temp,rho0_guest, mixflag)
    !    s0_guest(J) = G_CALC(gl,temp,rho0_guest, mixflag)
    !    fug_g0(J) = FUGCOPURE_CALC(gl,temp,rho0_guest, mixflag) * p0_guest
    !
    !    rhoredmix = rhoredmixorg
    !    tredmix = tredmixorg
    !end do
    !
    !
    !! Pure hydrates
    !! s =  x_hyd(1) * (-dChemPot_wH_dT + 1.D0 / Z * (s0_guest - Rgas*dlog(fug_g(1)/(p0_guest*1.D6))) )
    !sum = 0.d0
    !! Mixed hydrates 2015
    !s =  - x_hyd(1)*dChemPot_wH_dT
    !do J = 1,N_guests
    !    s =  s + x_hyd(J+1)*( s0_guest(J) - Rgas*dlog(fug_g(J)/(fug_g0(J)*1.D6)) )    ! Eq. 4.55 in Andreas' diss.
    !end do
    !do J = 1,N_guests
    !    sum =  sum + x_hyd(J+1)*( s0_guest(J) + temp*Rgas*dlog(fug_g(J)/(fug_g0(J)*1.D6)) )    ! Eq. 4.55 in Andreas' diss.
    !end do
    !write(999,*)p0_guest,sum,s
    !end do
    !close(999)
    !!TEST END

    !Double occ TEST START          Benedikt July 2017
    ! Verification of adapted eq. 4.57 and lattice routines for double occ

    p0_guest = 1.d-6
    !open(unit = 999, file='Z:\Semrau\HydMix\test.txt')
    !write(999,*)'p0_guest,sum,entropy'
    sum_s_guests=0.d0
    do J=1,gl%N_guests
        sum_s_guests = sum_s_guests + x_hyd(J+1)*(s0_guest(J) - Rgas*dlog(fug_g(J)/(p0_guest*1.d6))) ! eq. 4.55 - 2
    end do

    !check = s - s_part_guest

    !numerical derivative of dChemPot_wH_dT:

    tempp = temp *(1.d0 + 1.d-3)

    call hdrt_chem_potent_w(gl,tempp, press, fug_g, cphwp)

    tempm = temp * (1.d0 - 1.d-3)

    call hdrt_chem_potent_w(gl,tempm, press, fug_g, cphwm)

    dChemPot_wH_dT_num = (cphwp - cphwm) / (tempp - tempm)


    s_empty_lattice_div = - dChemPot_wH_dT_num - sum_s_guests - s_empty_lattice !not sure if it works



    !if (dabs(s_empty_lattice_div/s_empty_lattice)>1.d-3) then
    !	write(999,*) 'Entropy calculation incorrect'
    !	write(999,*) 'ceck lattice routines and derivates'
    !	write(999,*) 's_empty_lattice_div, s_empty_lattice, dChemPot_wH_dT_num, dChemPot_wH_dT'
    !	write(999,*) s_empty_lattice_div, s_empty_lattice, dChemPot_wH_dT_num, dChemPot_wH_dT
    !else
    !   write(999,*) ' '
    !write(999,*) 's_empty_lattice_div, s_empty_lattice, sum_s_guests, dChemPot_wH_dT_num, dChemPot_wH_dT'
    !write(999,*) s_empty_lattice_div, s_empty_lattice, sum_s_guests, dChemPot_wH_dT_num, dChemPot_wH_dT
    !   write(999,*) ' '
    !end if


    !Double occ TEST END



    end subroutine hdrt_entropy
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    subroutine hdrt_dC_dT(gl,temp, press, DCiJDT, error)
    !******************************************************************************
    !
    !  Temperature-derivative of the Langmuir constant
    !  of water in the hydrate phase [J/mol]
    !  June 2011, Andreas
    !  December 2015, Vaclav - updated for mixed hydrates




    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in):: temp, press
    ! Output arguments
    double precision, dimension(3,30), intent(out) :: DCiJDT
    integer :: error
    ! Local variables
    integer :: i
    double precision, dimension(2) :: T
    double precision, dimension(3,30,2) :: CiJ
    !//////////////////////////////////////////////////////////////////

    DCiJDT = 0.d0

    T(1) = temp*(1.d0+1.d-4)
    T(2) = temp*(1.d0-1.d-4)

    do i = 1, 2
        ! Langmuir constant @ T,p [-]
        call hdrt_Langmuir(gl,T(i),press,CiJ(:,:,i))
    end do

    DCiJDT = ( CiJ(:,:,1)-CiJ(:,:,2) )/(T(1)-T(2))

    end subroutine hdrt_dC_dT
    !******************************************************************************
    !******************************************************************************



    !******************************************************************************
    subroutine hdrt_dCd_dT(gl,temp, press, DCiJdDT)
    !******************************************************************************
    !
    !  Double occupancy version, Benedikt July 2017
    !
    !  Temperature-derivative of the Langmuir constant
    !  of water in the hydrate phase [J/mol]
    !  June 2011, Andreas
    !  December 2015, Vaclav - updated for mixed hydrates
    !  Double occupancy version, Benedikt July 2017




    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in):: temp, press
    ! Output arguments
    double precision, dimension(3,30), intent(out) :: DCiJdDT
    integer :: error
    ! Local variables
    integer :: i
    double precision, dimension(2) :: T
    double precision, dimension(3,30,2) :: CiJd
    !//////////////////////////////////////////////////////////////////

    DCiJdDT = 0.d0

    T(1) = temp*(1.d0+1.d-4)
    T(2) = temp*(1.d0-1.d-4)

    do i = 1, 2
        ! Langmuir constant @ T,p [-]
        call hdrt_Langmuir_double(gl,T(i),press,CiJd(:,:,i))
    end do

    DCiJdDT = ( CiJd(:,:,1)-CiJd(:,:,2) )/(T(1)-T(2))

    end subroutine hdrt_dCd_dT
    !******************************************************************************
    !******************************************************************************






    !******************************************************************************
    subroutine hdrt_Gibbs_energy(gl,x_hyd, chem_pot, g, error)
    !******************************************************************************
    !
    !  Gibbs enthalpy of hydrate [J/mol]
    !  January 2014
    !  IMPORTANT NOTE: THE PHASE EQUILIBRIUM WITH HYDRATE NEEDS TO BE CALCULATED BEFORE
    !  THIS ROUTINE IS CALLED!!!! THE COMPOSITION OF THE HYDRATE PHASE AND THE CHEMICAL
    !  POTENTIALS OF WATER AND THE GUEST MOLECULES IN HYDRATE ARE INPUTS TO THIS ROUTINE
    !   chempot(1) = chempot of water in hydrate
    !   chempot(2) = chempot of guest 1 in hydrate
    !   chempot(3) = chempot of guest 2 in hydrate
    !   ...
    !   x_hyd(1) = molfraction of water in hydrate
    !   x_hyd(1) = molfraction of guest 1 in hydrate
    !   x_hyd(1) = molfraction of guest 2 in hydrate
    !   ...





    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, dimension(30), intent(in) :: chem_pot, x_hyd
    ! Output arguments
    double precision, intent(out) :: g
    integer:: error
    ! Local variables
    integer:: n
    !//////////////////////////////////////////////////////////////////

    !Calculate the gibbs energy of the hydrate
    g = 0.D0
    do n = 1, gl%nrofhydrateformers
        g = g + x_hyd(n) * chem_pot(n)
    end do

    end subroutine hdrt_Gibbs_energy
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    subroutine hdrt_enthalpy(gl,Temp, press, x_hyd, chem_pot, fug_g, h, error)
    !******************************************************************************
    !
    !  Enthalpy of hydrate [J/mol]
    !  January 2014; VV December 2015
    !  IMPORTANT NOTE: THE PHASE EQUILIBRIUM WITH HYDRATE NEEDS TO BE CALCULATED BEFORE
    !  THIS ROUTINE IS CALLED!!!! THE COMPOSITION OF THE HYDRATE PHASE AND THE CHEMICAL
    !  POTENTIALS OF WATER AND THE GUEST MOLECULES IN HYDRATE ARE INPUTS TO THIS ROUTINE
    !   chempot(1) = chempot of water in hydrate
    !   chempot(2) = chempot of guest 1 in hydrate
    !   chempot(3) = chempot of guest 2 in hydrate
    !   ...
    !   x_hyd(1) = molfraction of water in hydrate
    !   x_hyd(1) = molfraction of guest 1 in hydrate
    !   x_hyd(1) = molfraction of guest 2 in hydrate
    !   ...





    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, dimension(30) :: chem_pot, x_hyd, fug_g
    double precision :: temp, press
    double precision :: dchem_pot_dT_fluid, dchem_pot_dx_fluid, dfug_dT_fluid, dfug_dx_fluid
    ! Output arguments
    double precision :: h
    ! Local variables
    double precision:: g, s
    integer:: error, n

    !//////////////////////////////////////////////////////////////////

    !Calculate the gibbs energy of hydrate
    call hdrt_Gibbs_energy(gl,x_hyd, chem_pot, g, error)

    !Calculate the entropy of hydrate
    call hdrt_entropy(gl,temp, press, x_hyd, fug_g, s, error)

    !Calculate the enthalpy of hydrate
    h = g + Temp * s

    end subroutine hdrt_enthalpy
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    subroutine hdrt_entropy_w_B_num(gl,temp,press, fug_g, sw_B, error)
    !******************************************************************************
    !
    !  Numerical calculation of the entropy of water in reference empty hydrate lattice [J/molK]
    !  February 2014



    !//////////////////////////////////////////////////////////////////

    implicit none

    type(type_gl) :: gl


    ! Input arguments
    double precision, intent(in) :: temp, press
    double precision, dimension(30), intent(in) :: fug_g
    ! Output arguments
    double precision, intent(out) :: sw_B
    ! Local variables
    integer:: error
    double precision, dimension(2) :: T, g_w_b
    !//////////////////////////////////////////////////////////////////

    T(1) = temp*(1.d0+1.d-4)
    T(2) = temp*(1.d0-1.d-4)

    call hdrt_gibbs_w_B(gl,T(1),press,fug_g,g_w_b(1))
    call hdrt_gibbs_w_B(gl,T(2),press,fug_g,g_w_b(2))

    sw_B = -(g_w_b(2) - g_w_b(1)) / (T(2) - T(1))

    end subroutine hdrt_entropy_w_B_num
    !******************************************************************************
    !******************************************************************************




    !******************************************************************************
    !subroutine hdrt_ancillary_hs(temp, d_fluid, press, x_fluid, chem_pot, dchem_pot_dT_fluid, dchem_pot_dx_fluid, &
    !    &   dfug_dT_fluid, dfug_dx_fluid, fug_g, error)
    subroutine hdrt_ancillary_hs(gl,temp, d_fluid, press, x_fluid, chem_pot, fug_g, error)
    !******************************************************************************
    ! Ancillary function for the calculation of chemical potential and fugacity
    ! derivatives that are needed for the enthalpy calculation of hydrates
    ! Andreas, March 2014; VV December 2015
    !INPUT: Temp, d_fluid, x_fluid, press                                !Temperature in K, density in mol/m³, composition of fluid phase, and pressure in PA!!!!!
    !OUTPUT: Chem_pot, d_chempot_d_T, d_chempot_d_x, d_f_dT, d_f_d_x, f  !Partial derivatives of the chemical potential and the fugacity w.r.t. T and x and the fugacity f. All pressures / fugacities and respective derivatives in Pa

    !//////////////////////////////////////////////////////////////////







    implicit none

    type(type_gl) :: gl


    double precision :: temp, d_fluid, press
    double precision, dimension(30) :: x_fluid

    double precision, dimension(30) :: chem_pot, fug_g
    double precision :: dchem_pot_dT_fluid, dchem_pot_dx_fluid, dfug_dT_fluid, dfug_dx_fluid

    double precision, dimension(30) :: dChempoti_dT_fluid, d2nadnidT, dnadni, dlnphiidT_fluid, dfi_dT_fluid
    double precision, dimension(30,30) :: dChempoti_dxj_fluid, dlnfidXj_fluid

    double precision :: p, Rmix,  rhoredmix_orig, tredmix_orig
    double precision, dimension(30) :: lnfi_fluid, molfrac_orig
    integer:: j, k, error

    !//////////////////////////////////////////////////////////////////

    chem_pot = 0.D0
    fug_g = 0.D0
    dchem_pot_dT_fluid = 0.D0
    dchem_pot_dx_fluid = 0.D0
    dfug_dT_fluid = 0.D0
    dfug_dx_fluid = 0.d0
    error = 0

    !Save the overall composition and respective reducing parameters
    rhoredmix_orig = gl%rhoredmix
    tredmix_orig = gl%tredmix
    molfrac_orig = gl%molfractions

    !recalculate the reducing parameters with the fluid phase
    gl%molfractions = x_fluid
    call reduced_parameters_calc(gl,temp)

    p = press / 1.D6    !MPa
    call lnf_mix(gl,Temp, d_fluid, p, lnfi_fluid)
    call Chempot_CALC(gl,Temp, d_fluid, chem_pot, 0)
    !Get the composition of the hydrate phase
    ! Pure hydrates: fug_g(1) = dexp(lnfi_fluid(2)) * 1.D6
    ! Mixed hydrates 2015
    fug_g(1:gl%N_guests) = dexp(lnfi_fluid(2:gl%N_hdrts)) * 1.D6

    !Set back module variables
    gl%rhoredmix = rhoredmix_orig
    gl%tredmix = tredmix_orig
    gl%molfractions = molfrac_orig

    end subroutine hdrt_ancillary_hs
    !******************************************************************************
    !******************************************************************************



    end module hdrt_chem_pot_module
