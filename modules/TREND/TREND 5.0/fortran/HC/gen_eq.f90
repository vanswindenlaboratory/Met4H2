! This file was generated -  Dont edit!
module gen_eq_file
! #######################################################
character(256), TARGET :: gen_eq(6)
DATA gen_eq(1) / 'FLUID           ALTERNATIVE NAME  	CAS-NR       	M / (g/mol)	Acentr. Fac.    rhocrit / mol/l Tcrit / K	Ptr / MPa	Ttr / K		polfac_eq         ' /
DATA gen_eq(2) / 'oil		oil			111-11-1	1773.0		0.0120177537	0.14028156	1413.697920257	0.000001D-07	0.1		-0.0020000000' /
DATA gen_eq(3) / 'acetone		acetone			67-64-1		58.07914	0.307		4.7		508.1		0.2326D-05	178.500		0' /
DATA gen_eq(4) / 'co2		co2			124-38-9	44.0098		0.224		7.3773		304.1282	0.51795D+00	216.592		0' /
DATA gen_eq(5) / 'methane         methane	                74-82-8         16.0428	        0.011	        4.5992          190.564         0.11696D-01     90.6941         0' /
DATA gen_eq(6) / '@END' /
contains
!----------------------------------------------------------------------------------
subroutine get_gen_eq_file(ptr,name)
implicit none
character(256),pointer,dimension(:)::ptr
character(*),dimension(:) :: name
ptr=>gen_eq
end subroutine get_gen_eq_file
end module gen_eq_file
