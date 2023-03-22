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
    
    module other_reader
	use module_all_types
	use file
	use String_Functions
	use calc_functions
	use cubic_eos_module
	use initialize_module
	use crit_tpd_module
	use rhomix_pt_module
	implicit none

	contains

	subroutine read_to_gl(gl,sub_part,nrsubst,read_flag,errorflag)
	type(type_gl):: gl
	type(return_pos_t):: sub_part
	integer :: nrsubst,read_flag,errorflag,errorflag_loop
	errorflag_loop = -1
	do while (errorflag_loop /= 0)
		if (read_flag == 2) then
			call read_simples_to_gl(gl,sub_part,nrsubst,errorflag)
		elseif (read_flag == 5) then
			call read_gen_to_gl(gl,sub_part,nrsubst,errorflag)
		elseif (read_flag == 6) then
			call read_saft_to_gl(gl,sub_part,nrsubst,errorflag)
		endif
		if (errorflag /= 0) then
			errorflag = -8878
			return
		endif
		if ((trim(gl%substfullname(nrsubst)) == trim(gl%components(nrsubst))) .OR. (trim(gl%substshortname(nrsubst)) == trim(gl%components(nrsubst)))  &
			& .OR. (trim(gl%substcasnr(nrsubst)) == trim(gl%components(nrsubst)))) then
			errorflag_loop = 0
		else
			errorflag_loop = -1
			sub_part%p = sub_part%p(2:size(sub_part%p))
			sub_part%pos_id = sub_part%pos_id(2:size(sub_part%pos_id))
		endif
	enddo
	end subroutine read_to_gl


	subroutine read_simples_to_gl(gl,sub_part,nrsubst,errorflag)

	type(type_gl):: gl
	type(return_pos_t):: sub_part
	integer :: nrsubst,errorflag
	!SRK: WARNING: AT THE MOMENT ONLY 30 GROUPS ARE READ IN. IF MORE GROUPS ARE NEEDED, THE FOLLOWING READ STATEMENT AS WELL AS THE FILE srk.fld NEED TO BE ADJUSTED
	read(sub_part%p(1),*,iostat= errorflag) gl%substfullname(nrsubst), gl%substshortname(nrsubst), gl%substcasnr(nrsubst), gl%wm(nrsubst), gl%accen(nrsubst), &
		& gl%pc(nrsubst), gl%tc(nrsubst), gl%ptp(nrsubst), gl%ttp(nrsubst), gl%A_cv0(nrsubst),&
		& gl%B_cv0(nrsubst), gl%C_cv0(nrsubst), gl%D_cv0(nrsubst),gl%E_cv0(nrsubst),gl%F_cv0(nrsubst),gl%G_cv0(nrsubst), gl%Cji_cubic(:,nrsubst), &
		& gl%nr_of_groups_i(nrsubst), gl%subgroups_ik_list(nrsubst,1:30)
	if (errorflag /= 0) then
		errorflag = 0
		!read of SRK without GROUPS, PR, and LKP
		read(sub_part%p(1),*,iostat= errorflag) gl%substfullname(nrsubst), gl%substshortname(nrsubst), gl%substcasnr(nrsubst), gl%wm(nrsubst), gl%accen(nrsubst), &
			& gl%pc(nrsubst), gl%tc(nrsubst), gl%ptp(nrsubst), gl%ttp(nrsubst), gl%A_cv0(nrsubst),&
			& gl%B_cv0(nrsubst), gl%C_cv0(nrsubst), gl%D_cv0(nrsubst),gl%E_cv0(nrsubst),gl%F_cv0(nrsubst),gl%G_cv0(nrsubst)
		if (errorflag /= 0) then
			errorflag = 0
			read(sub_part%p(1),*,iostat= errorflag) gl%substfullname(nrsubst), gl%substshortname(nrsubst), gl%substcasnr(nrsubst), gl%wm(nrsubst), gl%accen(nrsubst), &
				& gl%pc(nrsubst), gl%tc(nrsubst)
		endif
	endif
	end subroutine read_simples_to_gl



	subroutine read_simples(gl,path,folder,nrsubst,is_hc,errorflag)
	character(*) :: path,folder
	integer:: nrsubst,errorflag,errorflag_loop
	type(type_gl) :: gl
	type(file_t) :: SIMPLE_file
	type(return_pos_t):: sub_part
	character(256),dimension(:),allocatable:: get_lines
	character(10), dimension(1):: SIMPLE_fld_name
	!Variables needed for reading cp0-model of Joback
	integer :: i
	double precision, dimension (4) :: dA, dB, dC, dD   !dummy variables
	integer,dimension(30) :: N_m     !number of molecules of the type
	logical :: is_hc ! info if file is hc or normal path file

	! Body of read_SIMPLES

	! Init the file variable
	SIMPLE_file = file_t(is_hc)

	call get_content(SIMPLE_file,path,gl%components,folder)
	if (SIMPLE_file%errorflag /= 0) then
		errorflag = SIMPLE_file%errorflag
		return
	endif
	SIMPLE_fld_name(1) = gl%components(nrsubst)
	!SIMPLE_fld_name(2) = gl%substcasnr(nrsubst)
	sub_part = get_sub_string_lines(SIMPLE_file,SIMPLE_fld_name,errorflag)
	if (errorflag /= 0) return
	!if (size(sub_part%p) > 1) then
	!    continue
	!endif

	if (.not.allocated(gl%substshortname)) then
		allocate(gl%substshortname(gl%ncomp))
		allocate(gl%substfullname(gl%ncomp))
		allocate(gl%substchemname,gl%substsynonym,mold=gl%substfullname)
		allocate(gl%substcasnr,mold=gl%substshortname)
	end if

	if (gl%mix_type /= 12) then
		call read_to_gl(gl,sub_part,nrsubst,2,errorflag)
		gl%REQ(nrsubst) = 8.3144598D0  !J / mol K
		gl%tred(nrsubst) = gl%tc(nrsubst)
		gl%wm(nrsubst) = gl%wm(nrsubst)/gl%factor
		!New variables for Mathias Copeman              !Andreas Jäger, January 2017
		!Only use Mathias Copeman parameters if the PSRK is used (mix_type == 22)
		if (gl%mix_type /= 22) then
			gl%Cji_cubic(:,nrsubst) = 0.D0
		end if
		if (gl%VLE_needed) then
			!Andreas, November 2015: PURE SRK equation needs to be evaluated here!!
			!Set tredmix to tc(nrsubst)
			!THE NEXT LINE IS EXTREMELY IMPORTANT if the SRK is used in mixtures, because without this line the pure fluid is evaluated
			!at the wrong temperature and thus the critical density and all other calculated densities are wrong if tredmix = tc(nrsubst) is missing!
			!Andreas Jäger, September 2017
			gl%tredmix = gl%tc(nrsubst)
			if (gl%eq_type(nrsubst) == 2) then
				gl%rhored(nrsubst) = rho_SRK (gl, gl%tc(nrsubst), gl%pc(nrsubst), 0, nrsubst)!tc(nrsubst), pc(nrsubst), 1,
			elseif (gl%eq_type(nrsubst) == 3) then
				gl%rhored(nrsubst) = rho_PR (gl, gl%tc(nrsubst), gl%pc(nrsubst), 0, nrsubst)!tc(nrsubst), pc(nrsubst), 1,
			elseif (gl%eq_type(nrsubst) == 4) then
				gl%rhored(nrsubst) = 1.d0/((0.2905d0 - 0.085d0 * gl%accen(nrsubst)) * gl%tc(nrsubst) * gl%REQ(nrsubst) / (gl%pc(nrsubst)*1.d6)) !rho_LKP (tc(nrsubst), pc(nrsubst), 0, 1)!tc(nrsubst), pc(nrsubst), 1,
			endif
			gl%rhoc(nrsubst) = gl%rhored(nrsubst)
			if ((gl%ttp(nrsubst) > 1.d-14) .and. (gl%ptp(nrsubst) > 1.d-14)) then !Sebastian: Theresa neu funktioniert nicht wenn ttp und ptp = 0
				if (gl%eq_type(nrsubst) == 2) then
					gl%rhotp(nrsubst) = rho_SRK (gl, gl%ttp(nrsubst), gl%ptp(nrsubst), 0, 1) !Theresa neu
				elseif (gl%eq_type(nrsubst) == 3) then
					gl%rhotp(nrsubst) = rho_PR (gl, gl%ttp(nrsubst), gl%ptp(nrsubst), 0, 1) !Theresa neu
				elseif (gl%eq_type(nrsubst) == 4) then
					continue
				end if
			end if
			!Set back tredmix, Andreas September 2017
			gl%tredmix = 1.D0
		end if
		gl%pmaxfluid(nrsubst) = 1000.d0    !Needed for PhaseDet_TD (mixtures) (Dummy value)
		gl%tminfluid(nrsubst) = gl%ttp(nrsubst)       !Needed for Flash calculations, arbitrary value!! Andreas, May 2016
		!Changed from 1.d0 to ttp because for PH and PS Flash calculations ran into problems. Sebastian, October 2016
		gl%tmaxfluid(nrsubst) = 2000.d0    !Needed for Flash calculations, arbitrary value!! Andreas, May 2016

		!TEST_AJ - For paper Bell & Jäger 2016 values for SRK similar to Stradi et al. have been chosen August 2016
		!REQ(nrsubst) = 8.31451D0
		if (((dabs(gl%B_cv0(nrsubst)) > 1.D-8) .and. (dabs(gl%C_cv0(nrsubst)) > 1.D-8))) gl%cpmodel(nrsubst) = .true.
		gl%refstate(nrsubst) = "OT0"
		gl%tref(nrsubst) = 298.15d0
		gl%pref(nrsubst) = 0.101325d0  !"-1" shows that pref must be calculated later in another subroutine
		gl%href(nrsubst) = 0.d0 ![J/kg]     -later changed with wm into [J/mol]
		gl%sref(nrsubst) = 0.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
		gl%tcp0red(nrsubst) = gl%tc(nrsubst)
		gl%cp0red(nrsubst) = 8.3144598D0
		if (gl%eq_type(nrsubst) == 2) then
			gl%bi_SRK(nrsubst) = 0.08664D0 * gl%REQ(nrsubst) * gl%tc(nrsubst) / (gl%pc(nrsubst)*1.D6)
		elseif (gl%eq_type(nrsubst) == 3) then
			gl%bi_PR(nrsubst) = 0.0778D0*gl%R_PR * gl%tc(nrsubst) / (gl%pc(nrsubst)*1.D6)
		elseif (gl%eq_type(nrsubst) == 4) then
			gl%vc_LKP(nrsubst) = 1.d0/gl%rhoc(nrsubst)
			if (gl%ncomp  == 1) Then
				gl%zcLKP=(gl%pc(1)*1.d6)/(gl%tc(1)*gl%REQ(1)*gl%rhoc(1))
			end if
		end if
		gl%rhomaxfluid(nrsubst) = 1.D0/gl%bi_SRK(nrsubst)
		if (.not. gl%hold_limits) gl%ttp(nrsubst) = gl%ttp(nrsubst)*0.1d0    !small triple point temperature for the extrapolation of the VLE curves
		!     (if "&" occurs in the input variable)
		!New variables for group contribution methods   !Andreas Jäger, January 2017
		gl%nr_of_groups_mix = gl%nr_of_groups_mix + gl%nr_of_groups_i(nrsubst)
	elseif (gl%mix_type == 12) then
		read(sub_part%p(1),*,iostat= errorflag) gl%substfullname(nrsubst), gl%substshortname(nrsubst), gl%substcasnr(nrsubst), gl%nr_of_groups_i(nrsubst), gl%subgroups_ik_list(nrsubst,1:30)
		gl%nr_of_groups_mix = gl%nr_of_groups_mix + gl%nr_of_groups_i(nrsubst)
	endif



	if (gl%eq_type(nrsubst) == 2 .and. gl%components(nrsubst) == 'oil') then
		!Thu 2014/05
		!cp0-model of Joback (1986)
		!derivation 1...2%, increase with more complexe molecule structures
		!source: Proofs Wärmeatlas, Teil D1, p.21

		gl%cpmodel = .true.

		!set numbers of moleculegroup of substance
		N_m(1)=4         !CH3 (1)
		N_m(2)=26        !CH2 (2)
		N_m(3)=1         !CH  (3)
		N_m(4)=2         !C   (4)

		dA(1) = 19.5d0
		dB(1) = -8.08d-3
		dC(1) = 1.53d-4
		dD(1) = -9.67d-8
		dA(2) = -0.909d0
		dB(2) = 9.5d-2
		dC(2) = -5.44d-5
		dD(2) = 1.19d-8
		dA(3) = -23.d0
		dB(3) = 2.04d-1
		dC(3) = -2.65d-4
		dD(3) = 1.2d-7
		dA(4) = -66.2d0
		dB(4) = 4.27d-01
		dC(4) = -6.41d-04
		dD(4) = 3.01d-7

		do i=1,4
			gl%cp0_A=gl%cp0_A+N_m(i)*dA(i)
			gl%cp0_B=gl%cp0_B+N_m(i)*dB(i)
			gl%cp0_C=gl%cp0_C+N_m(i)*dC(i)
			gl%cp0_D=gl%cp0_D+N_m(i)*dD(i)
		end do

		! setting up common CPP module variables
		gl%ncp0poly(nrsubst) = 4
		gl%tcp0red(nrsubst) = 1.d0
		gl%cp0red(nrsubst) = 8.3144621d0
		gl%cp0exp(1,nrsubst) = 0.d0
		gl%cp0exp(2,nrsubst) = 1.d0
		gl%cp0exp(3,nrsubst) = 2.d0
		gl%cp0exp(4,nrsubst) = 3.d0
		gl%CP0COEFF(1,nrsubst) = (gl%cp0_A-37.93d0)/gl%cp0red(nrsubst)
		gl%CP0COEFF(2,nrsubst) = (gl%cp0_B+0.21d0)/gl%cp0red(nrsubst)
		gl%CP0COEFF(3,nrsubst) = (gl%cp0_C-3.91d-4)/gl%cp0red(nrsubst)
		gl%CP0COEFF(4,nrsubst) = (gl%cp0_D+2.06d-7)/gl%cp0red(nrsubst)

	end if


	end subroutine read_simples

	subroutine read_simples_mix(gl,nrsubst,j,path,folder,is_hc,errorflag)
	character(*) :: path,folder
	integer:: nrsubst,j,k,errorflag
	type(type_gl) :: gl
	type(file_t) :: SIMPLE_file_mix
	type(return_pos_t):: sub_part
	character(255) :: CASNR_1, CASNR_2  !In TREND the CASNR of the substances is used for the identification of mixing rules
	double precision :: kij_read, lij_read
	logical :: is_hc

	SIMPLE_file_mix = file_t(is_hc)

	! Body of read_simples_mix
	kij_read = 0.d0
	lij_read = 0.d0
	call get_content(SIMPLE_file_mix,path,gl%components,folder)
	if (SIMPLE_file_mix%errorflag /= 0) then
		errorflag = SIMPLE_file_mix%errorflag
		return
	endif
	sub_part = get_sub_string_lines(SIMPLE_file_mix,(/gl%substcasnr(nrsubst)/),errorflag)


	! cas number was not found in mix file
	if(errorflag .eq. 0) then

		SIMPLE_file_mix%content = ""
		SIMPLE_file_mix%content(1:size(sub_part%POS_ID,1)) = sub_part%p
		sub_part = get_sub_string_lines(SIMPLE_file_mix,(/gl%substcasnr(j)/),errorflag)
		if (size(sub_part%p,1) == 0) then
			!no binary interaction parameter for given combination of CASNRs
			errorflag = 0
			return
		endif
		read(sub_part%p(1),*,iostat = errorflag)  CASNR_1, CASNR_2, kij_read, lij_read
		if (errorflag /= 0) read(sub_part%p(1),*,iostat = errorflag)  CASNR_1, CASNR_2, kij_read
		if (errorflag /= 0) return
		if (gl%mix_type == 2) then
			gl%kij_SRK(j,nrsubst) = kij_read
			gl%kij_SRK(nrsubst,j) = kij_read
			gl%lij_SRK(j,nrsubst) = lij_read
			gl%lij_SRK(nrsubst,j) = lij_read
		elseif (gl%mix_type == 3) then
			gl%kij_PR(j,nrsubst) = kij_read
			gl%kij_PR(nrsubst,j) = kij_read
			gl%lij_PR(j,nrsubst) = lij_read
			gl%lij_PR(nrsubst,j) = lij_read
		elseif (gl%mix_type == 4) then
			gl%kij_LKP(j,nrsubst) = kij_read
			gl%kij_LKP(nrsubst,j) = kij_read
		elseif (gl%mix_type == 6) then
			gl%kij_PCSAFT(j,nrsubst) = kij_read
			gl%kij_PCSAFT(nrsubst,j) = kij_read
		end if
		k = index(sub_part%p(1),'!')
		if (k /= 0) then
			if (.not.allocated(gl%litref)) allocate(gl%litref)
			if (allocated(gl%litref)) then
				gl%litref%lit_ref_mix = sub_part%p(1)(k+1:)
			end if
		endif
	end if

	if(errorflag == -45257) then
		errorflag = 0
	end if

	end subroutine read_simples_mix

	subroutine read_gen_to_gl(gl,sub_part,nrsubst,errorflag)

	type(type_gl):: gl
	type(return_pos_t):: sub_part
	integer :: nrsubst,errorflag

	read(sub_part%p(1),*,iostat = errorflag) gl%substfullname(nrsubst), gl%substshortname(nrsubst), gl%substcasnr(nrsubst), gl%wm(nrsubst), gl%accen(nrsubst), gl%rhoc(nrsubst), &
		& gl%tc(nrsubst), gl%ptp(nrsubst), gl%ttp(nrsubst), gl%polfac_eq
	end subroutine read_gen_to_gl

	subroutine read_gen_eq(gl,path,folder,nrsubst,is_hc,errorflag)
	character(*) :: path,folder
	integer:: nrsubst,errorflag,gen_i,iter,ncomp_orig
	double precision, dimension(:,:), allocatable :: a_sun1   !factor for calculation of the generalized EOS of Sun & Ely
	double precision, dimension(:,:), allocatable :: a_sun2   !factor for calculation of the generalized EOS of Sun & Ely
	double precision :: tc_est_inp, pc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_
	logical :: is_hc ! info if file is hc or normal path file

	type(type_gl) :: gl
	type(file_t) :: GEN_EQ_file
	type(return_pos_t):: sub_part
	! Body of read_GEN_EQS

	if (.not.allocated(a_sun1)) allocate(a_sun1(100,30))
	if (.not.allocated(a_sun2)) allocate(a_sun2(100,30))

	if (.not.allocated(gl%substshortname)) then
		allocate(gl%substshortname(gl%ncomp))
		allocate(gl%substfullname(gl%ncomp))
		allocate(gl%substchemname,gl%substsynonym,mold=gl%substfullname)
		allocate(gl%substcasnr,mold=gl%substshortname)
	end if

	! Init the file variable
	GEN_EQ_file = file_t(is_hc)

	call get_content(GEN_EQ_file,path,gl%components,folder)

	if (GEN_EQ_file%errorflag /= 0) then
		errorflag = GEN_EQ_file%errorflag
		return
	endif

	sub_part = get_sub_string_lines(GEN_EQ_file,(/gl%components(nrsubst)/),errorflag)

	if (errorflag .ne. 0) then !if (size(sub_part%p,1) == 0) then
		!Component not found
		errorflag = -8878
		return
	endif

	call read_to_gl(gl,sub_part,nrsubst,5,errorflag)


	gl%rhoc(nrsubst) = gl%rhoc(nrsubst)*gl%factor

	gl%pmaxfluid(nrsubst) = 1d6!1.d12    !Needed for PhaseDet_TD (mixtures) (Dummy value)
	gl%REQ(nrsubst) = 8.31451D0
	gl%refstate(nrsubst) = "OT0"
	gl%tref(nrsubst) = 298.15d0
	gl%pref(nrsubst) = 0.101325d0
	gl%href(nrsubst) = 0.d0 ![J/kg]     -later changed with wm into [J/mol]
	gl%sref(nrsubst) = 0.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
	gl%tcp0red(nrsubst) = gl%tc(nrsubst)
	gl%cp0red(nrsubst) = 8.31451D0
	gl%rhomaxfluid(nrsubst) = 6.5D0*gl%rhoc(nrsubst)  !dummy value
	!Triple point information
	gl%rhotp(nrsubst) = 1.d0 !rho_SRK(ttp(nrsubst), ptp(nrsubst), 0, 1)  Theresa: there is no pc defined so far!

	if (gl%Eq_type(nrsubst) == 51) then    !calculation based on Alexandrov

		if (gl%accen(nrsubst) < 0.d0) then
			errorflag = -6001
			return
		end if



		call initialize_Gen_Eq_Igor(gl)
		gl%eos_coeff%nreg(nrsubst) = gl%anz_term_igor
		do gen_i = 1, gl%anz_term_igor
			gl%eos_coeff%ni(gen_i, nrsubst) = gl%coefficientsi(gen_i, 1) + (gl%coefficientsi(gen_i, 2)*gl%accen(nrsubst)) + (gl%coefficientsi(gen_i, 3)*gl%accen(nrsubst)**(gl%coefficientsi(gen_i, 4)))
			gl%eos_coeff%di(gen_i, nrsubst) = gl%exponentsi(gen_i, 1)
			gl%eos_coeff%ti(gen_i, nrsubst) = gl%exponentsi(gen_i, 2)
			gl%eos_coeff%p_i(gen_i, nrsubst) = gl%exponentsi(gen_i, 3)
			if (gl%eos_coeff%p_i(gen_i, nrsubst) /= 0) gl%eos_coeff%gama(gen_i, nrsubst) = 1.d0
			!write (71,2209) gl%eos_coeff%ni(gen_i, nrsubst), gl%eos_coeff%ti(gen_i, nrsubst), gl%eos_coeff%di(gen_i, nrsubst), gl%eos_coeff%p_i(gen_i, nrsubst)
		end do

	elseif (gl%Eq_type(nrsubst) == 52) then  !calculation based on Span/Wagner

		call initialize_Gen_Eq_Span(gl)
		gl%eos_coeff%nreg(nrsubst) = gl%anz_term_span

		do gen_i = 1, gl%anz_term_span
			gl%eos_coeff%ni(gen_i, nrsubst) = gl%coefficientss(gen_i, 1) + (gl%coefficientss(gen_i, 2)*gl%accen(nrsubst)) + (gl%coefficientss(gen_i, 3)*(gl%accen(nrsubst)**4))
			gl%eos_coeff%ti(gen_i, nrsubst) = gl%exponentss(gen_i, 1)
			gl%eos_coeff%di(gen_i, nrsubst) = gl%exponentss(gen_i, 2)
			gl%eos_coeff%gama(gen_i, nrsubst) = gl%exponentss(gen_i, 3)
			gl%eos_coeff%p_i(gen_i, nrsubst) = gl%exponentss(gen_i, 4)
			!write (71,2209) gl%eos_coeff%ni(gen_i, nrsubst), gl%eos_coeff%ti(gen_i, nrsubst), gl%eos_coeff%di(gen_i, nrsubst), gl%eos_coeff%p_i(gen_i, nrsubst)
		end do

	elseif (gl%Eq_type(nrsubst) == 53) then  !calculation based on Sun/Ely

		call initialize_Gen_Eq_Sun(gl)
		gl%eos_coeff%nreg(nrsubst) = gl%anz_term_sun

		do gen_i = 1, gl%anz_term_sun
			a_sun1(gen_i, nrsubst) = (gl%coefficientssun(gen_i,2)-gl%coefficientssun(gen_i,1))/(gl%acenfactor_2-gl%acenfactor_1)
			a_sun2(gen_i, nrsubst) = 1/(gl%polarfactor_3-gl%polarfactor_1)*((gl%coefficientssun(gen_i,3)-gl%coefficientssun(gen_i,1))-((gl%acenfactor_3-gl%acenfactor_1)/(gl%acenfactor_2-gl%acenfactor_1)*(gl%coefficientssun(gen_i,2)-gl%coefficientssun(gen_i,1))))

			gl%eos_coeff%ni(gen_i, nrsubst) = gl%coefficientssun(gen_i,1)+((gl%accen(nrsubst)-gl%acenfactor_1)*a_sun1(gen_i, nrsubst))+((gl%polfac_eq-gl%polarfactor_1)*a_sun2(gen_i, nrsubst))
			gl%eos_coeff%ti(gen_i, nrsubst) = gl%exponentssun(gen_i,2)
			gl%eos_coeff%di(gen_i, nrsubst) = gl%exponentssun(gen_i,1)
			gl%eos_coeff%p_i(gen_i, nrsubst) = gl%exponentssun(gen_i,3)
			if (gl%eos_coeff%p_i(gen_i, nrsubst) /= 0) gl%eos_coeff%gama(gen_i, nrsubst) = 1.d0
		end do

	end if

	tc_est_inp = gl%tc(nrsubst)
	rhoc_est_inp = gl%rhoc(nrsubst)
	gl%rhored(nrsubst) =  gl%rhoc(nrsubst)
	gl%Tred(nrsubst) =  gl%tc(nrsubst)
	!Warum das hier steht ??? auskommentiert
	!call find_crit_tpd (gl,tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errorflag, nrsubst)
	!if (errorflag == 0) then
	!    gl%tc(nrsubst) = tc_eos_
	!    gl%rhoc(nrsubst) = rhoc_eos_
	!    gl%pc(nrsubst) = pc_eos_
	!end if


	!if (ncomp > 1) call reduced_parameters_calc(300.D0) !Dummy temperature 300 K for the SRK
	if (gl%VLE_needed) then
		ncomp_orig = gl%ncomp
		!gl%ncomp = 1
		gl%rhoredmix = gl%rhored(nrsubst)
		gl%Tredmix = gl%Tred(nrsubst)
		gl%pc(nrsubst) = p_calc(gl,gl%tc(nrsubst),gl%rhoc(nrsubst),nrsubst)
		gl%rhotp(nrsubst) = rho_SRK(gl, gl%ttp(nrsubst), gl%ptp(nrsubst), 0, 1)
		!gl%ncomp = ncomp_orig
	end if
	!exit

	!Ideal part oil
	if (gl%components(nrsubst) == 'oil') then
		gl%cpmodel(nrsubst) = .true.
		! setting up common CPP module variables
		gl%tcp0red(nrsubst) = 1.d0
		gl%cp0red(nrsubst) = 8.3144621d0
		gl%ncp0poly(nrsubst) = 1
		gl%CP0COEFF(1,nrsubst) = 0.414793476829D+03
		gl%cp0exp(1,nrsubst) = 0.d0
		gl%ncp0pl(nrsubst) = 1
		gl%CP0COEFF(2,nrsubst) = 0.878393464574D+03
		gl%cp0exp(2,nrsubst) = 2281.40745904d0
	end if

	end subroutine read_gen_eq

	subroutine read_saft_to_gl(gl,sub_part,nrsubst,errorflag)

	type(type_gl):: gl
	type(return_pos_t):: sub_part
	integer :: nrsubst,errorflag,pos
	character(2):: model_PCSAFT
	character(50), dimension(50):: ref_tmp
	ref_tmp = ''

	!SRK: WARNING: AT THE MOMENT ONLY 30 GROUPS ARE READ IN. IF MORE GROUPS ARE NEEDED, THE FOLLOWING READ STATEMENT AS WELL AS THE FILE srk.fld NEED TO BE ADJUSTED
	read(sub_part%p(1),*,iostat= errorflag) gl%substfullname(nrsubst), gl%substshortname(nrsubst), model_PCSAFT, gl%substcasnr(nrsubst), gl%wm(nrsubst), &
		& gl%mPCSAFT(nrsubst), gl%sigPCSAFT(nrsubst), gl%epskPCSAFT(nrsubst), gl%epsabkPCSAFT(nrsubst), gl%kabPCSAFT(nrsubst), gl%QPCSAFTQ(nrsubst), gl%nPCSAFTQ(nrsubst), &
		& gl%MyPCSAFTD(nrsubst), gl%nPCSAFTD(nrsubst),  gl%accen(nrsubst), gl%rhoc(nrsubst), gl%tc(nrsubst), gl%ptp(nrsubst), gl%ttp(nrsubst), gl%A_cv0(nrsubst),&
		& gl%B_cv0(nrsubst), gl%C_cv0(nrsubst), gl%D_cv0(nrsubst),gl%E_cv0(nrsubst),gl%F_cv0(nrsubst),gl%G_cv0(nrsubst), ref_tmp
	gl%substfullname(nrsubst) = Replace_Text(gl%substfullname(nrsubst),'_',' ')

	pos = index(gl%substfullname(nrsubst),'#')
	if (pos.ne.0) then
		gl%substfullname(nrsubst)=gl%substfullname(nrsubst)(pos+1:)
	end if
	gl%litref%lit_ref_saft(nrsubst) =  Copy_a2s(ref_tmp)

	gl%SAFTmodel_assoc(nrsubst) = model_PCSAFT(2:2)
	read (model_PCSAFT(1:1),*) gl%n_sites(nrsubst)

	!errorflag -1 means end-of-file -> sub_part%p(1) is only one row in size -> end-of-file after reading one line
	if (errorflag == -1) errorflag = 0
	end subroutine read_saft_to_gl


	subroutine read_saft(gl,path,folder,nrsubst,is_hc,errorflag)
	character(*) :: path,folder
	integer:: nrsubst,errorflag,i, iter, riter
	type(type_gl) :: gl
	type(file_t) :: SAFT_file,SAFT_sub_part
	type(return_pos_t),target:: sub_part
	character(50), dimension(1):: SAFT_fld_name, substring
	double precision :: tc_est_inp, pc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_
	logical :: is_hc ! info if file is hc or normal path file

	!DEC$ IF DEFINED(WO_PCSAFT)
	errorfld = -7878
	!DEC$ ELSE

	! Body of read_saft

	! Init the file variable
	SAFT_file = file_t(is_hc)

	SAFT_fld_name = '#' // Replace_Text(trim(gl%components(nrsubst)),' ','_')
	call get_content(SAFT_file,path,SAFT_fld_name,folder)
	if (SAFT_file%errorflag /= 0) then
		errorflag = SAFT_file%errorflag
		return
	endif
	sub_part = get_sub_string_lines(SAFT_file,SAFT_fld_name,errorflag)
	SAFT_sub_part%content = sub_part%p
	SAFT_sub_part%nr_lines = size(SAFT_sub_part%content)
	sub_part = get_sub_string_lines(SAFT_sub_part,SAFT_fld_name,errorflag)

	if(.not.allocated(gl%litref)) allocate(gl%litref)

	if (.not.allocated(gl%substshortname)) then
		allocate(gl%substshortname(gl%ncomp))
		allocate(gl%substfullname(gl%ncomp))
		allocate(gl%substchemname,gl%substsynonym,mold=gl%substfullname)
		allocate(gl%substcasnr,mold=gl%substshortname)
	end if

	call read_to_gl(gl,sub_part,nrsubst,6,errorflag)


	if ((trim(gl%substfullname(nrsubst)) == trim(gl%components(nrsubst))) .OR. (trim(gl%substshortname(nrsubst)) == trim(gl%components(nrsubst)))  &
		& .OR. (gl%substcasnr(nrsubst) == gl%components(nrsubst))) then


		if(gl%unitin == 3) then
			gl%Factor = 1.D0
			gl%factorpress=1.D6
			gl%factortrans=1.D0
			gl%factorrbwr=1.D0
			gl%factor_VS5eta=1.D0
			gl%factor_VS5slj=1.D0
		else                         !for calculating real fluids
			gl%Factor = 1.D3
			gl%factorpress=1.D0
			gl%factortrans=1.D6
			gl%factorrbwr=1.D2
			gl%factor_VS5eta=1.D5
			gl%factor_VS5slj=1.D1
		end if

		gl%tred(nrsubst) = 1.d0 ! derivatives with respect to temperature can be used but need to be transformed
		gl%rhored(nrsubst) = 1.d0 ! derivatives with respect to density can be used
		gl%Eq_type(nrsubst) = 6
		gl%Mix_type = 6
		gl%pmaxfluid(nrsubst) = 1000.d0    !Needed for PhaseDet_TD (mixtures) (Dummy value)
		gl%wm(nrsubst) = gl%wm(nrsubst)/gl%factor
		if (gl%tc(nrsubst) < 5.d0) then
			gl%REQ(nrsubst) = 1.D0
		else
			gl%REQ(nrsubst) = 8.31451D0
		end if
		if (((dabs(gl%B_cv0(nrsubst)) > 1.D-8) .and. (dabs(gl%C_cv0(nrsubst)) > 1.D-8))) gl%cpmodel(nrsubst) = .true.
		gl%refstate(nrsubst) = "OT0"
		gl%tref(nrsubst) = 298.15d0
		gl%pref(nrsubst) = 0.101325d0  !"-1" shows that pref must be calculated later in another subroutine
		gl%href(nrsubst) = 0.d0 ![J/kg]     -later changed with wm into [J/mol]
		gl%sref(nrsubst) = 0.d3   ![J/(kg K)] -later changed with wm into [J/(mol K)]
		gl%tcp0red(nrsubst) = 0d0!gl%tc(nrsubst)
		gl%cp0red(nrsubst) = 8.31451D0

		gl%rhoredmix = 1.d0


		!special case elongation given
		if (gl%mPCSAFT(nrsubst) < 0.d0 .and. (gl%QPCSAFTQ(nrsubst) > 0.d0 .or. (gl%QPCSAFTQ(nrsubst) == 0.d0 .and. gl%MyPCSAFTD(nrsubst) == 0.d0))) then
			gl%mPCSAFT(nrsubst) = 1.d0 + 0.2177d0*dabs(gl%mPCSAFT(nrsubst)/gl%sigPCSAFT(nrsubst)) + 3.1498d0*dabs(gl%mPCSAFT(nrsubst)/gl%sigPCSAFT(nrsubst))**2 &
				& - 3.6738d0*dabs(gl%mPCSAFT(nrsubst)/gl%sigPCSAFT(nrsubst))**3 + 1.3063d0*dabs(gl%mPCSAFT(nrsubst)/gl%sigPCSAFT(nrsubst))**4
		else if (gl%mPCSAFT(nrsubst) < 0.d0 .and. gl%MyPCSAFTD(nrsubst) > 0.d0) then
			gl%mPCSAFT(nrsubst) = 1.d0 + 0.1795d0*dabs(gl%mPCSAFT(nrsubst)/gl%sigPCSAFT(nrsubst)) + 3.3283d0*(gl%mPCSAFT(nrsubst)/gl%sigPCSAFT(nrsubst))**2 &
				& - 3.8855d0*(dabs(gl%mPCSAFT(nrsubst)/gl%sigPCSAFT(nrsubst)))**3 + 1.3777d0*(gl%mPCSAFT(nrsubst)/gl%sigPCSAFT(nrsubst))**4
		end if


		!determine true critical point of EOS
		tc_est_inp = gl%tc(nrsubst)
		rhoc_est_inp = gl%rhoc(nrsubst)*gl%factor
		!if (trim(gl%components(nrsubst)) == 'pec5') rhoc_est_inp = 550.d0
		!if (trim(gl%components(nrsubst)) == 'pec7') rhoc_est_inp = 410.d0

		call find_crit_tpd (gl,tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errorflag, nrsubst)


		riter = 1
		rhoc_est_inp = gl%rhoc(nrsubst)*gl%factor
		do while (errorflag .ne. 0)
			errorflag = 0
			rhoc_est_inp = rhoc_est_inp*1.02d0
			call find_crit_tpd (gl,tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errorflag, nrsubst)
			riter = riter + 1
			if (errorflag == 0 .or. riter > 105) exit
		end do

		riter = 1
		rhoc_est_inp = gl%rhoc(nrsubst)*gl%factor
		do while (errorflag  .ne. 0)
			errorflag = 0
			rhoc_est_inp = rhoc_est_inp*0.98d0
			call find_crit_tpd (gl,tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errorflag, nrsubst)
			riter = riter + 1
			if (errorflag == 0 .or. riter > 105) exit
		end do

		if (errorflag .ne. 0) then  !calculate the critical parameters from the Lennard-Jones fluid if the given ones are too far off
			errorflag = 0
			!pc_est_inp = 0.13d0 * 1.3806504d-23 * gl%epskPCSAFT(nrsubst)/(1.d-30*gl%sigPCSAFT(nrsubst)**3)*1.d-6
			tc_est_inp = 1.32d0 * gl%epskPCSAFT(nrsubst)
			rhoc_est_inp = 0.31d0 / (6.02214d-7*gl%sigPCSAFT(nrsubst)**3) /2.d0    !TE /2.d0
			call find_crit_tpd (gl,tc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_, iter, errorflag, nrsubst)
			if (errorflag .ne. 0) then
				return
			end if
		end if

		gl%tc(nrsubst) = tc_eos_
		gl%rhoc(nrsubst) = rhoc_eos_
		gl%pc(nrsubst) = pc_eos_

		if (gl%ttp(nrsubst) < 1.d-13) then   !calculate the triple point temperaturefrom the Lennard-Jones fluid if there is nothing given
			gl%ttp(nrsubst) = 0.66d0 * gl%epskPCSAFT(nrsubst)
			gl%rhotp(nrsubst) = 0.865d0 / (6.02214d-7*gl%sigPCSAFT(nrsubst)**3)
			gl%ptp(nrsubst) = P_CALC(gl, gl%ttp(nrsubst), gl%rhotp(nrsubst), nrsubst)
			gl%rhomaxfluid(nrsubst) = gl%rhotp(nrsubst) * 1.5d0
		else
			gl%rhomaxfluid(nrsubst) = 4.d0 * gl%rhoc(nrsubst)
			gl%rhotp(nrsubst) = rhomix_calc (gl, gl%ttp(nrsubst), gl%ptp(nrsubst), 0.d0, 1, nrsubst)
			if (gl%rhotp(nrsubst) < gl%rhoc(nrsubst)) then
				gl%rhotp(nrsubst) = 3.d0 * gl%rhoc(nrsubst)
				gl%rhomaxfluid(nrsubst) = 4.d0 * gl%rhoc(nrsubst)
			else
				gl%rhomaxfluid(nrsubst) = 1.5d0 * gl%rhotp(nrsubst)
			end if
		end if

		!rhomaxfluid(nrsubst) = 3.9d0*rhoc(nrsubst)

	end if

	!DEC$ END IF
	end subroutine read_saft

	!---
	subroutine read_rkm(gl,path,folder,nrsubst,is_hc,errorflag)
	character(*) :: path,folder
	integer:: nrsubst,errorflag,i, j
	type(type_gl) :: gl
	type(file_t) :: RKM_file
	type(return_pos_t):: sub_part
	character(30)::dummy
	character(10), dimension(1):: RKM_text
	logical :: is_hc ! info if file is hc or normal path file

	! Init the file variable
	RKM_file = file_t(is_hc)

	!path_file = path
	if(is_hc) then !If its HC we need to pass the names of the files as name and not the name of the component
		call get_content(RKM_file,trim(path) // '/rkm_correction_factors.fld',(/'RKM_correction_factors'/),folder)
	else
		call get_content(RKM_file,trim(path) // '/rkm_correction_factors.fld',(/gl%components(nrsubst)/),folder)
	endif


	if (RKM_file%errorflag /= 0) then
		errorflag = RKM_file%errorflag
		return
	endif
	read(RKM_file%content(3:13),*)gl%k1(0,0:10), gl%k1(1,0:10), gl%k1(2,0:10), gl%k1(3,0:10), gl%k1(4,0:10), gl%k1(5,0:10), gl%k1(6,0:10), gl%k1(7,0:10), gl%k1(8,0:10), gl%k1(9,0:10), gl%k1(10,0:10)
	read(RKM_file%content(16:26),*) gl%k2(0,0:10), gl%k2(1,0:10), gl%k2(2,0:10), gl%k2(3,0:10), gl%k2(4,0:10), gl%k2(5,0:10), gl%k2(6,0:10), gl%k2(7,0:10), gl%k2(8,0:10), gl%k2(9,0:10), gl%k2(10,0:10)

	gl%k1(1:10,1:10) = gl%k1(1:10,1:10) * 1.d-3
	gl%k2(1:10,1:10) = gl%k2(1:10,1:10) * 1.d-3

	!path_file = path
	if(is_hc) then !If its HC we need to pass the names of the files as name and not the name of the component
		call get_content(RKM_file,trim(path) // '/rkm_correction_factors.fld',(/'RKM_fluid_information'/),folder)
	else
		call get_content(RKM_file,trim(path) // '/rkm_fluid_information.fld',(/gl%components(nrsubst)/),folder)
	end if
	
	do i = 1,3
		read(RKM_file%content(i+1),*) gl%RKM_fluids(1,i),gl%RKM_fluids(2,i),gl%RKM_fluids(3,i),gl%RKM_fluids(4,i),gl%RKM_fluids(5,i),gl%RKM_fluids(6,i),gl%RKM_fluids(7,i),gl%RKM_fluids(8,i)
	enddo

	read(RKM_file%content(6),*) gl%RKM_Tc(1),gl%RKM_Tc(2),gl%RKM_Tc(3),gl%RKM_Tc(4),gl%RKM_Tc(5),gl%RKM_Tc(6),gl%RKM_Tc(7),gl%RKM_Tc(8)

	read(RKM_file%content(8),*) gl%RKM_M(1),gl%RKM_M(2),gl%RKM_M(3),gl%RKM_M(4),gl%RKM_M(5),gl%RKM_M(6),gl%RKM_M(7),gl%RKM_M(8)
	RKM_text = '140'
	sub_part = get_sub_string_lines(RKM_file,RKM_text,errorflag)

	Do i = 0,sub_part%POS_ID(1)-10
		read(RKM_file%content(i+10),*) gl%RKM_ps_me(i,0), gl%RKM_ps_me(i,1)
	enddo

	Do i = 0,sub_part%POS_ID(2)-(sub_part%POS_ID(1)+2)
		read(RKM_file%content(i+sub_part%POS_ID(1)+2),*) gl%RKM_V(i,0), gl%RKM_V(i,1), gl%RKM_V(i,2), gl%RKM_V(i,3), gl%RKM_V(i,4), gl%RKM_V(i,5), gl%RKM_V(i,6), gl%RKM_V(i,7), gl%RKM_V(i,8)
	enddo
	!##############################################################################################
	!Sort Mole Fractions
	!##############################################################################################

	do i = 1, 8
		do j = 1, 8
			if ((gl%components(j) == gl%RKM_fluids(i,1)) .or. (gl%components(j) == gl%RKM_fluids(i,2)) .or. (gl%components(j) == gl%RKM_fluids(i,3))) then
				gl%wm(j) = gl%RKM_M(i)
				exit
			end if
		end do
	end do

	gl%wm(:) = gl%wm(:) /gl%factor           ![kg/mol]

	end subroutine read_rkm

	!---
	subroutine read_costald(gl,path,folder,nrsubst,is_hc,errorflag)
	character(*) :: path,folder
	integer:: nrsubst,errorflag
	type(type_gl) :: gl
	type(file_t) :: COSTALD_file
	type(return_pos_t):: sub_part
	character(30)::dummy, dummy1, dummy2
	character(10), dimension(1):: COSTALD_CAS_NR
	logical :: is_hc ! info if file is hc or normal path file

	! Init the file variable
	COSTALD_file = file_t(is_hc)

	COSTALD_CAS_NR = gl%substcasnr(nrsubst)

	! Get the content of the current file
	call get_content(COSTALD_file,path,(/gl%components(nrsubst)/),folder)

	! Check if an error occured
	if (COSTALD_file%errorflag /= 0) then
		errorflag = COSTALD_file%errorflag
		return
	endif

	! Subpart where CASNR is
	sub_part = get_sub_string_lines(COSTALD_file,COSTALD_CAS_NR,errorflag)

	read(sub_part%p(1),*,iostat=errorflag) dummy, dummy1, gl%omega(nrsubst), gl%vstar(nrsubst)

	!fluid name maybe consists of two words
	if (errorflag .ne. 0) then
		errorflag = 0
		read(sub_part%p(1),*,iostat=errorflag) dummy, dummy1, dummy2, gl%omega(nrsubst), gl%vstar(nrsubst)
	end if

	if (errorflag .ne. 0) errorflag = -8001

	end subroutine read_costald

	end module other_reader