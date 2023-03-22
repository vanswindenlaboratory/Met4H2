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
    
    
    module ge_reader
    use file
    use other_reader
    use module_all_types
    contains

    subroutine read_unifac_mix(gl,nrsubst,nrsubst_compo2,path,folder,is_hc,errorflag)
    implicit none
    character(*) :: path,folder
    character(255) :: filename
    integer:: nrsubst,nrsubst_,i,j,k,m,n,errorflag, nrsubst_compo2
    type(type_gl):: gl
    type(file_t),allocatable :: GE_file_mix
    type(return_pos_t):: sub_part
    integer, allocatable, dimension(:,:) :: group_matrix

    logical :: exists,is_hc

    integer :: ge_unit

    !New variables for gE-models
    integer :: nr_main_group, nr_sub_group, subgroup_found
    double precision:: vdW_vol_Rk, vdW_ar_Qk
    integer:: nr_of_groups_i_read, errval
    integer, dimension(100) :: subgroups_ik_list_read
    double precision, dimension(3):: Cji_cubic_read
    integer:: line_nr, n_gr, m_gr, n_gr_read, m_gr_read, leng
    double precision :: a_nm_read, b_nm_read, c_nm_read, a_mn_read, b_mn_read, c_mn_read
    double precision :: tc_est_inp, pc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_



    !For Helm+gE model, read subgroups of fluids from file
    if (gl%mix_type == 12) then
        if (nrsubst == 1 .and. nrsubst_compo2 == 2) then
            do nrsubst_ = 1, gl%ncomp   !Loop over all components in the mixture
                !Check if groups file exist at given path
                filename = trim(path) // 'ge_models/unifac/helm_ge/groups.par'
                if (is_hc) filename = 'groups'
                !2020-12, Andy: kann eventuell gelöscht werden
                call read_simples(gl,filename,folder,nrsubst_,is_hc,errorflag)
                if (errorflag /= 0) return

            end do
        endif


        !Andreas Jäger, April 2018. Use linear mixing rules as default
        gl%LIN_or_LB = 2
        if (gl%LIN_or_LB == 1) then
            !LB-mixing rules
            gl%rfbetat(nrsubst,nrsubst_compo2) = 1.d0
            gl%rfbetarho(nrsubst,nrsubst_compo2) = 1.d0
            gl%rfgammat(nrsubst,nrsubst_compo2) = 1.d0
            gl%rfgammarho(nrsubst,nrsubst_compo2) = 1.d0
            gl%rfbetat(nrsubst_compo2,nrsubst) = 1.d0
            gl%rfbetarho(nrsubst_compo2,nrsubst) = 1.d0
            gl%rfgammat(nrsubst_compo2,nrsubst) = 1.d0
            gl%rfgammarho(nrsubst_compo2,nrsubst) = 1.d0
        elseif (gl%LIN_or_LB == 2) then
            !Linear mixing rules
            gl%rfbetat(nrsubst,nrsubst_compo2) = 1.d0
            gl%rfbetarho(nrsubst,nrsubst_compo2) = 1.d0
            gl%rfbetat(nrsubst_compo2,nrsubst) = 1.d0
            gl%rfbetarho(nrsubst_compo2,nrsubst) = 1.d0
            gl%rfgammat(nrsubst,nrsubst_compo2) = 0.5d0*(gl%tc(nrsubst)+gl%tc(nrsubst_compo2))/dsqrt(gl%tc(nrsubst)*gl%tc(nrsubst_compo2))
            gl%rfgammat(nrsubst_compo2,nrsubst) = gl%rfgammat(nrsubst,nrsubst_compo2)
            gl%rfgammarho(nrsubst,nrsubst_compo2) = 4.d0*(1.d0/gl%rhoc(nrsubst)+1.d0/gl%rhoc(nrsubst_compo2))/(gl%rhoc(nrsubst)**(-1.d0/3.d0)+gl%rhoc(nrsubst_compo2)**(-1.d0/3.d0))**3
            gl%rfgammarho(nrsubst_compo2,nrsubst) = gl%rfgammarho(nrsubst,nrsubst_compo2)
        else
            errorflag = -11223344   !Dummy, Sebastian und Sven gefällt das .-D
            return
        End if
    end if

    !Determine the start indices for the groups of each pure fluid when they are written in one common vector
    gl%group_start_index(1) = 1
    Do i = 2, gl%ncomp
        gl%group_start_index(i) = gl%group_start_index(i-1) + gl%nr_of_groups_i(i-1)
    End do
    m = 0

    !Read the van der Waals volume RK and the van der Waals area Qk of the groups in the mixture
    if (gl%mix_type == 12) then
        if (is_hc) then
            folder = 'helm_ge'
            filename = 'rk_qk'
        else
            filename = trim(path) // 'ge_models/unifac/helm_ge/rk_qk.par'
        endif
    elseif (gl%mix_type == 22) then
        if (is_hc) then
            folder = 'psrk'
            filename = 'rk_qk'
        else
            filename = trim(path) // 'ge_models/unifac/psrk/rk_qk.par'
        endif
    end if
    allocate(GE_file_mix)
    call get_content(GE_file_mix,filename,gl%components,folder)
    if (GE_file_mix%errorflag /= 0) then
        errorflag = GE_file_mix%errorflag
        return
    endif
    if (.not. allocated(group_matrix))allocate(group_matrix(GE_file_mix%nr_lines,2))
    do i = 2,GE_file_mix%nr_lines-2
        read(GE_file_mix%content(i),*)group_matrix(i,1),group_matrix(i,2)
    enddo
    if(.not.allocated(gl%R_ik)) then
        allocate(gl%R_ik(gl%ncomp,maxval(gl%nr_of_groups_i)))
        allocate(gl%Q_ik,gl%v_ik,mold=gl%R_ik)
    end if
    do i=1,gl%ncomp
        do k = 1, gl%nr_of_groups_i(i)
            subgroup_found = findloc(group_matrix(:,2),value = gl%subgroups_ik_list(i,k),dim=1)
            read(GE_file_mix%content(subgroup_found),*,iostat= errorflag) group_matrix(subgroup_found,1),group_matrix(subgroup_found,2),gl%R_ik(i,k),gl%Q_ik(i,k)
            !gl%R_ik(i,k) = vdW_vol_Rk
            !gl%Q_ik(i,k) = vdW_ar_Qk
            gl%v_ik(i,k) = 1   !ALL GROUPS ARE TREATED SEPERATELY AT THE MOMENT, EVEN IF THEY EXIST MULTIPLE TIMES IN A SUBSTANCE. THIS COULD BE CHANGED IN THE FUTURE BY COUNTING HOW OFTEN CERTAIN GROUPS OCCUR, BUT IT IS NOT CLEAR AT THE MOMENT IF THIS IS ADVANTAGEOUS OR NOT
            !sort the main group into the correct position of the vector maingroups_mix_list
            gl%maingroups_mix_list(gl%group_start_index(i) + k - 1) = group_matrix(subgroup_found,1)
            !Count how many groups have already been assigned. If all groups have been found, quit reading
            m = m + 1
            if (m == gl%nr_of_groups_mix) exit
        end do
        if (m == gl%nr_of_groups_mix) exit
    end do

    if (m < gl%nr_of_groups_mix) then
        errorflag = -8879
        return
    end if
    deallocate(GE_file_mix)

    !Read the binary interactions parameters "anm", "bnm", and "cnm" for all groups in the mixture
    if (gl%mix_type == 12) then
        if (is_hc) then
            folder = 'helm_ge'
            filename = 'Interac'
        else
            filename = trim(path) // 'ge_models/unifac/helm_ge/Interac.par'
        endif
    elseif (gl%mix_type == 22) then
        if (is_hc) then
            folder = 'psrk'
            filename = 'Interac'
        else
            filename = trim(path) // 'ge_models/unifac/psrk/Interac.par'
        endif
    end if
    allocate(GE_file_mix)
    call get_content(GE_file_mix,filename,gl%components,folder)
    if (GE_file_mix%errorflag /= 0) then
        errorflag = GE_file_mix%errorflag
        return
    endif
    !The binary parameters file has a well defined structure. In the first column the main groups are listed in numerical order (1...85, Index "n"), in the second column
    !the main groups interacting with the main groups of the first column are given, starting at (maingroup_column1 + 1), index m. After that, the up to six interaction parameters
    !are given: anm, bnm, cnm, amn, bmn, cmn
    !Example: maingroup n   maingroup m     anm     bnm     cnm     amn     bmn     cmn
    !                   1           2       ....
    !                   1           3       ....
    !                   1           4       ....
    !                   ...
    !                   1           85      ....
    !                   2           3       ....
    !                   2           4       ....
    !                   ...
    !Thus, if the groups n and m are given for which the interaction parameters are needed, the line where this information stands can be calculated with the formula
    ! line = (n-1)*85-n*(n-1)/2+(m-n)
    !allocate a_nm var
    if(.not.allocated(gl%a_nm)) then
        allocate(gl%a_nm(gl%nr_of_groups_mix,gl%nr_of_groups_mix))
        allocate(gl%b_nm,gl%c_nm,mold=gl%a_nm)
    end if
    Do i = 1, gl%nr_of_groups_mix
        n = gl%maingroups_mix_list(i)
        Do j = i, gl%nr_of_groups_mix
            m = gl%maingroups_mix_list(j)
            if (n < m) then     !This is the way that the groups are stored in the file, no switching necessary
                n_gr = n
                m_gr = m
                line_nr = (n_gr-1)*85-n_gr*(n_gr-1)/2+(m_gr-n_gr)
                if (errorflag == 0) then
                    read(GE_file_mix%content(line_nr),*,iostat = errorflag) n_gr_read, m_gr_read, gl%a_nm(i,j), gl%b_nm(i,j), gl%c_nm(i,j), gl%a_nm(j,i), gl%b_nm(j,i), gl%c_nm(j,i)
                    if ((n_gr /= n_gr_read) .or. (m_gr /= m_gr_read) .or. (errorflag /= 0)) then
                        errorflag = -7897
                        return
                    end if

                else
                    errorflag = -7897
                    return
                end if
            elseif (n > m) then              !Switch groups, such that n_gr < m_gr
                n_gr = m
                m_gr = n
                line_nr = (n_gr-1)*85-n_gr*(n_gr-1)/2+(m_gr-n_gr)
                if (errorflag == 0) then
                    read(GE_file_mix%content(line_nr),*,iostat = errorflag) n_gr_read, m_gr_read, gl%a_nm(j,i), gl%b_nm(j,i), gl%c_nm(j,i), gl%a_nm(i,j), gl%b_nm(i,j), gl%c_nm(i,j)
                    if ((n_gr /= n_gr_read) .or. (m_gr /= m_gr_read) .or. (errorflag /= 0)) then
                        errorflag = -7897
                        return
                    end if
                else
                    errorflag = -7897
                    return
                end if
            elseif (n == m) then       !Interaction of same main group -> all interaction parameters 0
                gl%a_nm(j,i) = 0.D0
                gl%b_nm(j,i) = 0.D0
                gl%c_nm(j,i) = 0.D0
                gl%a_nm(i,j) = 0.D0
                gl%b_nm(i,j) = 0.D0
                gl%c_nm(i,j) = 0.D0
            end if
        end do
    end do
    end subroutine read_unifac_mix


    subroutine read_cosmo_mix(gl,nrsubst,nrsubst_compo2,path,folder,is_hc,errorflag)
    implicit none
    character(*) :: path,folder
    character(255) :: filename
    integer:: nrsubst,nrsubst_,i,j,k,m,n,errorflag, nrsubst_compo2
    type(type_gl):: gl
    !type(file_t),allocatable :: GE_file_mix
    type(file_t),allocatable :: GE_file_mix
    type(return_pos_t):: sub_part
    integer, allocatable, dimension(:,:) :: group_matrix

    logical :: exists,is_hc

    integer :: ge_unit

    !New variables for gE-models
    integer :: nr_main_group, nr_sub_group, subgroup_found
    double precision:: vdW_vol_Rk, vdW_ar_Qk
    integer:: nr_of_groups_i_read, errval
    integer, dimension(100) :: subgroups_ik_list_read
    double precision, dimension(3):: Cji_cubic_read
    integer:: line_nr, n_gr, m_gr, n_gr_read, m_gr_read, leng
    double precision :: a_nm_read, b_nm_read, c_nm_read, a_mn_read, b_mn_read, c_mn_read
    double precision :: tc_est_inp, pc_est_inp, rhoc_est_inp, tc_eos_, rhoc_eos_, pc_eos_

    !New variables for COSMO_SAC model
    integer, dimension(30) :: index_nr
    character(len=255) :: CASNR_COSMO
    character(4), dimension(30) :: index_char
    character(len=255) :: dummy_cosmo   !used to jump over the not needed lines of the inputfile
    character(1) :: spt     !sigma-profile type
    integer :: interval, profile_count, err_read
    integer, dimension(30,3) :: counter_min, counter_max

    !Electrostatic interaction
    double precision :: sigmaacc, sigmadon              !Acceptor and Donator (hydrogen bonding interactions) in calculation of electrostatic interaction
    double precision, dimension(:,:),allocatable :: delta_w                         !Electrostatic interactions
    double precision :: eo, eps, sigma_hb, aeff
    double precision :: coord, c_hb, fpol, alpha, alpha_prime
    double precision :: c_es, c_OH_OH, c_OT_OT, c_OH_OT
    integer :: l
    !genralized equations
    integer:: gen_i
    double precision :: p_i_double



    !------------------------------------------------
    !Set COSMO-SAC Version
    gl%cosmo%COSMO_ver =  2
    !analytical derivations for COSMO-SAC Version 1?
    gl%cosmo%analytical = .false.
    !successive substitution or Newton-Raphson methode for solving segment activity coefficients; Newton-Raphson = true; successive substitution = false
    gl%cosmo%solv_method = .false.
    !------------------------------------------------

    interval=51
    filename = trim(path) // 'ge_models/cosmo-sac/VT-2005_Sigma_Profile_Database_Index_v2_orig.txt'


    ! Init the file variable
    allocate(GE_file_mix)
    GE_file_mix = file_t(is_hc)

    call get_content(GE_file_mix,filename,gl%components,folder)
    if (GE_file_mix%errorflag /= 0) then
        errorflag = GE_file_mix%errorflag
        return
    endif

    do l = 1, gl%ncomp

        index_nr(l) = 0
        index_char(l) = ""
        sub_part = get_sub_string_lines(GE_file_mix,(/gl%substcasnr(l)/),errorflag)
        do k = 1,size(sub_part%p,1)
            !if a casnr contains a substring that is another casnr then you have to find the correct entry
            Read(sub_part%p(k),*,iostat = errorflag) index_nr(l), CASNR_COSMO, dummy_cosmo, gl%cosmo%Vcosmo(l)
            if (CASNR_COSMO == gl%substcasnr(l)) exit
        enddo

        write(index_char(l),'(i4)') index_nr(l)
        if (index_nr(l) < 10) then
            index_char(l) = "000" // index_char(l)(4:4)
        end if
        if ((index_nr(l) < 100) .and. (index_nr(l) > 9)) then
            index_char(l) = "00" // index_char(l)(3:4)
        end if
        if ((index_nr(l) < 1000) .and. (index_nr(l) > 99)) then
            index_char(l) = "0" // index_char(l)(2:4)
        end if

    end do

    deallocate(GE_file_mix)

    !---------------------------------------------------------------------------------------------------------------------------------
    !Read Sigma-Profiles
    !---------------------------------------------------------------------------------------------------------------------------------
    if (gl%cosmo%COSMO_ver == 1) then
        !Read Sigma-Profile Files for components and calculate COSMO surface area Acosmo (COSMO-SAC Version 1)
        do l = 1, gl%ncomp

            if(allocated(GE_file_mix)) then
                deallocate(GE_file_mix)
            end if

            filename = trim(path) // 'ge_models/cosmo-sac/sigma-profiles/VT2005-' // index_char(l) // '.txt'
            allocate(GE_file_mix)
            GE_file_mix = file_t(is_hc)

            call get_content(GE_file_mix,filename,gl%components,folder)

            if (GE_file_mix%errorflag /= 0) then
                errorflag = -7898
                return
            end if
            gl%cosmo%Acosmo(l) = 0
            Do k = 1, interval
                Read(GE_file_mix%content(k),*,iostat = errorflag) gl%cosmo%counter(k), gl%cosmo%sigma(k, l)
                gl%cosmo%Acosmo(l) = gl%cosmo%Acosmo(l) + gl%cosmo%sigma(k, l)
                if (errorflag /= 0) then
                    errorflag = -7898
                    return
                end if
            end do
            !November 2018, Erik
            !read min and max values for occupied intervals
            Read(GE_file_mix%content(k),*,iostat = errorflag) counter_min(l,1)
            Read(GE_file_mix%content(k+1),*,iostat = errorflag) counter_max(l,1)
            profile_count = 1

            !Read total Volume of cavity
            read(GE_file_mix%content(k+2),*,iostat = errorflag) gl%cosmo%Vcosmo(l)
            deallocate(GE_file_mix)
        end do
    elseif ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3) .or. (gl%cosmo%COSMO_ver == 4)) then
        !Read Sigma-Profile Files for components and calculate COSMO surface area Acosmo (COSMO-SAC Version 2 and 3)
        !update version 4 (Xiong), Feb. 7th 2018
        do l = 1, gl%ncomp

            if(allocated(GE_file_mix)) then
                deallocate(GE_file_mix)
            end if

            gl%cosmo%Acosmo(l) = 0
            if (gl%cosmo%COSMO_ver == 4) then
                filename = trim(path) // 'ge_models/cosmo-sac/sigma-profiles_ver_4/VT2005-' // index_char(l) // '.txt'
            else
                filename = trim(path) // 'ge_models/cosmo-sac/sigma-profiles3/VT2005-' // index_char(l) // '.txt'
            end if
            allocate(GE_file_mix)
            GE_file_mix = file_t(is_hc)
            call get_content(GE_file_mix,filename,gl%components,folder)

            if (GE_file_mix%errorflag /= 0) then
                errorflag = -7898
                return
            end if
            do i = 1, 3
                do k = 1, interval
                    Read(GE_file_mix%content(interval*(i-1)+k),*,iostat = errorflag) gl%cosmo%counter_v23(k, i), gl%cosmo%sigma_v23(k, l, i)
                    gl%cosmo%Acosmo(l) = gl%cosmo%Acosmo(l) + gl%cosmo%sigma_v23(k, l, i)
                    if (errorflag /= 0) then
                        errorflag = -7898
                        return
                    end if
                end do
            end do
            !Read dispersive interaction parameters
            if ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3)) then
                read(GE_file_mix%content(154),*,iostat = errorflag) gl%cosmo%molecule_type(l)
                read(GE_file_mix%content(155),*,iostat = errorflag) gl%cosmo%eps_molecule(l)
                if (errorflag /= 0) then
                    errorflag = -7898
                    return
                end if
            elseif (gl%cosmo%COSMO_ver == 4) then
                read(GE_file_mix%content(1),*,iostat = errorflag) gl%cosmo%molecule_type(l)
                do j = 1, 17
                    read(GE_file_mix%content(j+1),*,iostat = errorflag) dummy_cosmo, gl%cosmo%m_tau(l,j)
                    if (errorflag /= 0) then
                        errorflag = -7898
                        return
                    end if
                end do
            end if
            !November 2018, Erik
            !read min and max values for occupied intervals
            do i = 1, 3
                Read(GE_file_mix%content(156+(i-1)*2),*,iostat = errorflag) counter_min(l,i)
                Read(GE_file_mix%content(157+(i-1)*2),*,iostat = errorflag) counter_max(l,i)
            end do
            profile_count = 3

            !Read total Volume of cavity
            read(GE_file_mix%content(162),*,iostat = errorflag) gl%cosmo%Vcosmo(l)
        end do
    end if

    !---------------------------------------------------------------------------------------------------------------------------------
    !November 2018, Erik
    !min and max value for occupied intervals for all compounds
    !Min
    do i = 1, profile_count
        gl%cosmo%interval_min(i) = 51
        do l = 1, gl.ncomp
            if ((counter_min(l,i) < gl%cosmo%interval_min(i)) .and. (counter_min(l,i) > 0.D0))  then
                gl%cosmo%interval_min(i) = counter_min(l,i)
            end if
        end do
        !set value for selected profile to 0, if not occupied at all
        if (gl%cosmo%interval_min(i) == 51) then
            gl%cosmo%interval_min(i) = 0
        end if
    end do
    !Max
    do i = 1, profile_count
        gl%cosmo%interval_max(i) = 0
        do l = 1, gl.ncomp
            if (counter_max(l,i) > gl%cosmo%interval_max(i)) then
                gl%cosmo%interval_max(i) = counter_max(l,i)
            end if
        end do
    end do

    !Erik, December 2018
    !to increase calculation speed allocate large matrices with the lowest size
    !very extensive to change, because everything would have to be rearranged
    !gl%cosmo%int = gl%cosmo%interval_max(1) - gl%cosmo%interval_min(1) + 1

    !---------------------------------------------------------------------------------------------------------------------------------
    !Calculate electrostatic interaction
    !---------------------------------------------------------------------------------------------------------------------------------

    if ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3)) then    !Hsieh
        aeff = 7.25    !surface area of a standard segment in Angstrom^2
    elseif (gl%cosmo%COSMO_ver == 1) then    !Mullins
        aeff = 7.5    !surface area of a standard segment in Angstrom^2
    elseif (gl%cosmo%COSMO_ver == 4) then    !Xiong
        aeff = 7.9    !surface area of a standard segment in Angstrom^2
    end if

    !Version 1

    if (gl%cosmo%COSMO_ver == 1) then
        !Calculate electrostatic interaction COSMO-SAC Version 1
        eo=2.395D-4
        eps=3.667D0   !(LIN AND SANDLER USE A CONSTANT FPOL WHICH YEILDS EPS=3.68)
        sigma_hb=0.0084D0
        c_hb=85580.D0
        fpol = (eps-1.D0)/(eps+0.5D0)
        alpha = (0.3D0*aeff**(1.5D0))/(eo)
        alpha_prime = fpol*alpha
        do j = gl%cosmo%interval_min(1), gl%cosmo%interval_max(1)
            do k = gl%cosmo%interval_min(1), gl%cosmo%interval_max(1)
                if (gl%cosmo%counter(j)>=gl%cosmo%counter(k)) then
                    sigmaacc= gl%cosmo%counter(j)
                    sigmadon = gl%cosmo%counter(k)
                end if
                if (gl%cosmo%counter(j)<gl%cosmo%counter(k)) then
                    sigmadon = gl%cosmo%counter(j)
                    sigmaacc = gl%cosmo%counter(k)
                end if
                gl%cosmo%delta_w_gl(j,k) = (alpha_prime/2.0)*(gl%cosmo%counter(j)+gl%cosmo%counter(k))**2.0 + c_hb *   &
                    max(0.D0,(sigmaacc - sigma_hb)) * min(0.D0,(sigmadon + sigma_hb))
            end do
        end do
    end if

    !Version 2 to 4
    if ((gl%cosmo%COSMO_ver == 2) .or. (gl%cosmo%COSMO_ver == 3) .or. (gl%cosmo%COSMO_ver == 4)) then
        c_OH_OH = 4013.78D0 * 4184.D0     !hydrogen bonding interaction coefficient
        c_OT_OT = 932.31D0 * 4184.D0      !hydrogen bonding interaction coefficient
        c_OH_OT = 3016.43D0 * 4184.D0     !hydrogen bonding interaction coefficient
        !in Versions 2 and 3 of COSMO-SAC electrostatic interactions are temperature dependend. This is formulated in c_ES and has to be calculated
        !in COSMO_SAC_CALC. The rest can be calculated here
        do k = 1, 3
            if (gl%cosmo%interval_min(k) == 0) then
                cycle
            end if
            do l = 1, 3
                if (gl%cosmo%interval_min(l) == 0) then
                    cycle
                end if
                do i = gl%cosmo%interval_min(k), gl%cosmo%interval_max(k)
                    do j = gl%cosmo%interval_min(l), gl%cosmo%interval_max(l)
                        if ((k == 2) .and. (l == 2) .and. (gl%cosmo%counter_v23(i,k) * gl%cosmo%counter_v23(j,l) < 0)) then
                            c_hb = c_OH_OH
                        elseif ((((k == 2) .and. (l == 3)) .or. ((k == 3) .and. (l == 2))) .and. (gl%cosmo%counter_v23(i,k) * gl%cosmo%counter_v23(j,l) < 0)) then
                            c_hb = c_OH_OT
                        elseif ((k == 3) .and. (l == 3) .and. (gl%cosmo%counter_v23(i,k) * gl%cosmo%counter_v23(j,l) < 0)) then
                            c_hb = c_OT_OT
                        else
                            c_hb = 0.D0
                        end if
                        gl%cosmo%sum_square_counter(i,j,k,l) = (gl%cosmo%counter_v23(i,k) + gl%cosmo%counter_v23(j,l)) ** 2.D0
                        gl%cosmo%hb_inter(i,j,k,l) = c_hb * (gl%cosmo%counter_v23(i,k) - gl%cosmo%counter_v23(j,l)) ** 2.D0
                    end do
                end do
            end do
        end do

    end if

    !---------------------------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------------------------------------------------------
    !Read NBT experimental Data for pure fluids (Version 4 (Xiong)
    if (gl%cosmo%COSMO_ver == 4) then
        filename = trim(path) // 'ge_models/cosmo-sac/NBT_pure.txt'
        do l = 1, gl%ncomp
            open(unit=12, file=trim(filename), status='old', action='read', iostat=errorflag)
            if (errorflag /= 0) then
                errorflag = -7898
                return
            end if
            read(12,*,iostat = errorflag) dummy_cosmo
            read(12,*,iostat = errorflag) dummy_cosmo
            do i = 1, 2000
                read(12,*,iostat = errorflag) dummy_cosmo, gl%cosmo%NBT_cosmo(l), dummy_cosmo, dummy_cosmo, dummy_cosmo, dummy_cosmo, CASNR_COSMO
                if (trim(CASNR_COSMO) == trim(gl%substcasnr(l))) then
                    exit
                end if
                if (errorflag /= 0) then
                    errorflag = -7898
                    close(12)
                    return
                end if
            end do
            close(12)
        end do

    end if

    !if (LIN_or_LB == 1) then
    !    !LB-mixing rules
    !    rfbetat(1,2) = 1.d0
    !    rfbetarho(1,2) = 1.d0
    !    rfgammat(1,2) = 1.d0
    !    rfgammarho(1,2) = 1.d0
    !    rfbetat(2,1) = 1.d0
    !    rfbetarho(2,1) = 1.d0
    !    rfgammat(2,1) = 1.d0
    !    rfgammarho(2,1) = 1.d0
    !elseif (LIN_or_LB == 2) then
    !Linear mixing rules
            gl%rfbetat(nrsubst,nrsubst_compo2) = 1.d0
            gl%rfbetarho(nrsubst,nrsubst_compo2) = 1.d0
            gl%rfbetat(nrsubst_compo2,nrsubst) = 1.d0
            gl%rfbetarho(nrsubst_compo2,nrsubst) = 1.d0
            gl%rfgammat(nrsubst,nrsubst_compo2) = 0.5d0*(gl%tc(nrsubst)+gl%tc(nrsubst_compo2))/dsqrt(gl%tc(nrsubst)*gl%tc(nrsubst_compo2))
            gl%rfgammat(nrsubst_compo2,nrsubst) = gl%rfgammat(nrsubst,nrsubst_compo2)
            gl%rfgammarho(nrsubst,nrsubst_compo2) = 4.d0*(1.d0/gl%rhoc(nrsubst)+1.d0/gl%rhoc(nrsubst_compo2))/(gl%rhoc(nrsubst)**(-1.d0/3.d0)+gl%rhoc(nrsubst_compo2)**(-1.d0/3.d0))**3
            gl%rfgammarho(nrsubst_compo2,nrsubst) = gl%rfgammarho(nrsubst,nrsubst_compo2)
    !else
    !    errorflag = -11223344   !Dummy
    !    return
    !End if
    end subroutine read_cosmo_mix
    end module ge_reader