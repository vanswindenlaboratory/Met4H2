! This file was generated -  Dont edit!
module hard_coded_files
use binary_mix_files
use costald_files
use eos_cg_2016_binary_mix_files
use eos_cg_2016_pure_fld_files
use eos_cg_2019_binary_mix_files
use eos_cg_2019_pure_fld_files
use eos_lng_mix_files
use eos_lng_pure_fld_files
use pure_fld_files
use gen_eq_file
use gerg_2008_mix_files
use pure_gerg_fld_files
use unifac_helm_ge_files
use unifac_psrk_files
use lkp_fld_file
use lkp_mix_file
use pr_fld_file
use pr_mix_file
use rkm_files
use pc_saft_pure_files
use pc_saft_mix_files
use srk_fld_file
use srk_mix_file
implicit none
type ptr
character(500),POINTER,dimension(:)::p
end type
type hc_t
    type(ptr) :: act_file_ptr
end type
!**********************************************************************************************************************
!interface for standard constructor
interface hc_t
module procedure init
end interface
!**********************************************************************************************************************
contains
type(hc_t) function init(name,folder)
character(*),dimension(:)::name
character(*) :: folder
if(folder .eq. 'binary_mix_files') then
call get_b_mix_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'costald') then
call get_costald_file(init%act_file_ptr%p,name)
elseif(folder .eq. 'eos-cg-2016_mix') then
call get_eos_cg_2016_b_mix_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'eos-cg-2016_pure') then
call get_eos_cg_2016_pure_fld_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'eos-cg-2019_mix') then
call get_eos_cg_2019_b_mix_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'eos-cg-2019_pure') then
call get_eos_cg_2019_pure_fld_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'eos-lng_mix') then
call get_eos_lng_b_mix_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'eos-lng_pure') then
call get_eos_lng_pure_fld_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'fluids_pure') then
call get_pure_fld_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'gen_eq') then
call get_gen_eq_file(init%act_file_ptr%p,name)
elseif(folder .eq. 'gerg-2008_mix') then
call get_gerg_2008_b_mix_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'gerg-2008_pure') then
call get_pure_gerg_fld_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'helm_ge') then
call get_unifac_helm_ge_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'psrk') then
call get_unifac_psrk_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'lkp_pure') then
call get_pure_lkp_file(init%act_file_ptr%p,name)
elseif(folder .eq. 'lkp_mix') then
call get_mix_lkp_file(init%act_file_ptr%p,name)
elseif(folder .eq. 'pr_pure') then
call get_pure_pr_file(init%act_file_ptr%p,name)
elseif(folder .eq. 'pr_mix') then
call get_mix_pr_file(init%act_file_ptr%p,name)
elseif(folder .eq. 'rkm') then
call get_rkm_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'saft_pure') then
call get_pc_saft_pure_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'binary_saft') then
call get_pc_saft_mix_files(init%act_file_ptr%p,name)
elseif(folder .eq. 'srk_pure') then
call get_pure_srk_file(init%act_file_ptr%p,name)
elseif(folder .eq. 'srk_mix') then
call get_mix_srk_file(init%act_file_ptr%p,name)
else
init%act_file_ptr%p=> NULL()
end if
end function init
end module hard_coded_files
