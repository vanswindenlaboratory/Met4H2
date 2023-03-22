! This file was generated -  Dont edit!
module pc_saft_mix_files
! #######################################################
character(256), TARGET :: saft(25)
DATA saft(1) / '!In this file the binary mixing parameters for PC-SAFT fluid comibations are stored' /
DATA saft(2) / '!The order is as follows: ' /
DATA saft(3) / '! CAS-NR Fluid1		CAS-NR Fluid 2		kij		!Comments' /
DATA saft(4) / '64-17-5			7732-18-5		-0.0616		!ethanol-water 				Ref: Christoph Held, TU-Dortmund, personal communication (2018)' /
DATA saft(5) / '124-38-9		15834-04-5		0.0881		!CO2-PEC5 					Ref: O. Fandi?o, E.R. L?pez, L. Lugo, J. Garc?a, J. Fern?ndez, J. Supercrit. Fluids (2010) ' /
DATA saft(6) / '67-64-1			124-18-5		0.033		!acetone-decane 			Ref: Gross & Vrabec (2006), AIChE J. 52:1194-1204  Acetone mit PCP! ' /
DATA saft(7) / '67-64-1			106-97-8		0.037		!acetone-n-butane 			Ref: Gross & Vrabec (2006), AIChE J. 52:1194-1204  Acetone mit PCP! ' /
DATA saft(8) / '67-64-1			74-84-0			0.037		!acetone-ethane 			Ref: Gross & Vrabec (2006), AIChE J. 52:1194-1204  Acetone mit PCP! ' /
DATA saft(9) / '67-64-1			109-66-0		0.024		!acetone-pentane 			Ref: Gross & Vrabec (2006), AIChE J. 52:1194-1204  Acetone mit PCP! ' /
DATA saft(10) / '67-68-5			108-88-3		-0.012		!dmso-toluene	 			Ref: Gross & Vrabec (2006), AIChE J. 52:1194-1204  Acetone mit PCP! ' /
DATA saft(11) / '124-38-9		106-97-8		0.036		!CO2-n-butane 				Ref: Gross (2005), AIChE J. 51:2556-2568  CO2 mit PCP! ' /
DATA saft(12) / '124-38-9		142-82-5		0.039		!CO2-n-heptane 				Ref: Gross (2005), AIChE J. 51:2556-2568  CO2 mit PCP!' /
DATA saft(13) / '124-38-9		108-88-3		0.03		!CO2-toluene 				Ref: Gross (2005), AIChE J. 51:2556-2568  CO2 mit PCP!' /
DATA saft(14) / '67-56-1			78-83-1			0.05		!methanol-isobutanol		Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(15) / '67-56-1			110-82-7		0.051		!methanol-cyclohexan		Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(16) / '67-56-1			111-87-5		0.020		!methanol-1-octanol			Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(17) / '67-56-1			75-28-5			0.05		!methanol-isobutane			Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(18) / '64-17-5			106-97-8 		0.028		!ethanol-n-butane			Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(19) / '71-23-8			71-43-2			0.020		!1-propanol-benzene			Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(20) / '71-23-8			100-41-4		0.023		!1-propanol-ethylbenzene	Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(21) / '67-63-0			71-43-2			0.021		!2-propanol-benzene			Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(22) / '71-36-3			106-97-8		0.015		!1-butanol-n-butane			Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(23) / '71-41-0			71-43-2			0.0135		!1-pentanol-benzene			Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(24) / '7732-18-5		71-41-0l		0.016		!water-1-pentanol			Ref: Gross & Sadowski (2002), Ind. Eng. Chem. Res. 41:5510-5515' /
DATA saft(25) / '@END' /
contains
!----------------------------------------------------------------------------------
subroutine get_pc_saft_mix_files(ptr,name)
implicit none
character(256),pointer,dimension(:)::ptr
character(*),dimension(:) :: name
ptr=>saft
end subroutine get_pc_saft_mix_files
end module pc_saft_mix_files
