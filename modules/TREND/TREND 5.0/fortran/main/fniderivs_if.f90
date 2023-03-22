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

    ! module for file fniderivs.f90
    module fniderivs_module
    !global use inclusion
    use module_all_types


    interface	
    module SUBROUTINE FNIDERIVS(gl,T, D, GETDERIVS, SETDERIVS, nrsubst)
    implicit none
    type(type_gl) :: gl
     integer, dimension(nderivsi)::GETDERIVS 
     DOUBLE PRECISION, dimension(nderivsi)::SETDERIVS
     integer::nrsubst

    double precision ::T, D
    END SUBROUTINE FNIDERIVS
 


    module subroutine ref_calc(gl,nrsubst, error)
    implicit none
    type(type_gl) :: gl
    integer:: nrsubst
    integer:: error
    end subroutine

 
 
    module SUBROUTINE IDEALCONSTS(gl,nrsubst)
    implicit none
    type(type_gl) :: gl
    integer::nrsubst
    END SUBROUTINE
  
 
    module SUBROUTINE MIXDERIVSFNI(gl,T, D, GETDERIVS, SETMIXDERIVS)
    implicit none
    type(type_gl) :: gl
    integer, dimension(nderivsi)::GETDERIVS 
    DOUBLE PRECISION, dimension(nderivsi)::SETMIXDERIVS
    double precision ::T, D
    END SUBROUTINE MIXDERIVSFNI

    end interface

    contains

    end module fniderivs_module
