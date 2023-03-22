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

    module module_regula_falsi_support
    type type_additional_parameters
        Double Precision :: a_p(65)
        integer :: i = 0 !initialization is needed to set all values of a_p = 0.d0
    end type type_additional_parameters

    !interface type_additional_parameters
    !module procedure init
    !end interface

    contains

    !type (type_additional_parameters) function init()
    !init%a_p = 0.d0
    !end function init

    end module module_regula_falsi_support

    ! module for file subroutine_regula_falsi.f90
    module module_regula_falsi
    !global use inclusion
    use module_all_types
    use module_regula_falsi_support


    contains




    Subroutine Regula_Falsi(gl,Zerofunction, Xroot, Xstart_min, Xstart_max, Delta_allowed, Xmin_allowed, Xmax_allowed, &
        &                       Max_iterations, Iterations, Errorflag, t_params)
    !
    ! Modified Regula Falsi algorithm optimized for finding roots in thermodynamic property functions.
    ! The subroutine first searchs for an interval including exactly one root. Then it applies a Regula Falsi like algorithm
    ! to find the root. The algorithm is written in a way that the found interval cannot be left again to make sure that it
    ! does not run into physically unreasonable solutions. The starting intervall Xstart_min / Xstart_max is essential for this
    ! purpose. Xmin_allowed / Xmax_allowed define an intervall that is considered physically reasonable. Roots outside of this
    ! range are avoided in any case.
    ! R. Span, Denmark, 09-2009


    implicit none

    type(type_gl) :: gl

    ABSTRACT INTERFACE
    DOUBLE PRECISION FUNCTION Zerofunction_PROTOTYPE(gl,X,Parameters)
    use module_all_types
    use module_regula_falsi_support
    type(type_gl) :: gl
    DOUBLE PRECISION :: X
    type(type_additional_parameters) :: parameters
    END FUNCTION Zerofunction_PROTOTYPE
    END INTERFACE
    PROCEDURE(Zerofunction_PROTOTYPE), POINTER, intent(in) :: Zerofunction


    !                               the value Xroot with Zerofunction = 0 has to be found
    Double Precision :: Xmin_allowed, Xmax_allowed !Minimum and maximum value allowed for the independant variable X
    Double Precision :: Xstart_min, Xstart_max !Starting interval for the variable X, within which the root is searched for
    Double Precision :: Xakt ! Current value of the independent variable
    Double Precision :: Xmin, Xmax ! Minimum and maximum value of the interval currently investigated
    Double Precision :: F_Xmin, F_Xmax, F_Xakt ! Value of Zerofunktion at the corresponding values of X
    Double Precision :: Xroot ! Result returned by the root finding algorithm
    Double Precision :: Delta_allowed ! Allowed uncertainty of the iteration in X
    Double Precision :: Deriv ! Numerical derivative of Zerofunction calculated using the functional values at Xmin and Xmax
    Double Precision :: Deviation_X ! Checks how far the current value Xakt is from the root
    Double Precision :: Xtest, F_Xtest ! Values at a point calculated to improve the speed of convergence
    Double Precision :: F_Xmin_alt, Xmin_alt, F_Xmax_alt, Xmax_alt  ! Parameters used to avoid the enlargement on the side of the interval that makes things worse
    type(type_additional_parameters) :: t_params ! Field with up to 65 additional parameters forwarded to the respective
    ! zerofunction. These parameters are not changed during the iteration process but allow for forwarding routine specific
    ! constant values (e.g. T, compositions, ...). The vector is not used in Regula_Falsi, except for the call of zerofunction.

    integer :: Iterations, Max_Iterations ! Actually needed number of iterations (returned value) and allowed number of iterations
    integer :: Errorflag ! Errorflag returned by the root finding algorithm
    integer :: Stop_min, Stop_max ! Flags that stop the enlargement of the interval in the respective direction
    ! Errorflag = 0 = >  Iteration succeeded
    ! Errorflag = 1 = >  Number of allowed iterations exceeded
    ! Errorflag = 2 = >  No root found without conflicts with Xmin_allowed or Xmax_allowed
    ! Errorflag = 3 = >  The initial interval for X had to be enlarged or reduced to find a root

    integer :: Int_Errorflag ! Internal errorflag
    integer :: Count_enlarge, Count_reduce ! Internal counters for the case, that the initial interval needs to be enlarged or reduced

    Logical :: Root_In, Root_Found, additional_checks ! Define whether the root is within the intervall Xmin, Xmax or whether it has been found

    ! Reset some values

    Iterations = 0
    Errorflag = 0
    Int_Errorflag = 0
    Count_enlarge = 0
    Count_reduce = 0
    F_Xmin = 0.D0
    F_Xmax = 0.D0
    Xroot = 0.D0
    Root_In = .False.
    Root_Found = .False.
    Stop_min = 0
    Stop_max = 0
    Xakt = 0.D0
    F_Xakt = 0.D0
    Xmin = 0.d0
    Xmax = 0.d0
    Deriv = 0.d0
    Deviation_X = 0.d0
    Xtest = 0.d0
    F_Xtest = 0.d0
    F_Xmin_alt = 0.d0
    Xmin_alt = 0.d0
    F_Xmax_alt = 0.d0
    Xmax_alt = 0.d0
    additional_checks = .false.

    ! Transfer starting values for Xmin and Xmax

    Xmin = Xstart_min
    Xmax = Xstart_max

    ! In case both the limits for Xmin, Xmax exceed the allowed range (unreasonable call of the routine) the allowed values are assigned
    ! and the internal errorflag is set to try again with smaller intervals

    If (Xmin < Xmin_allowed) Int_Errorflag = 1
    If (Xmin < Xmin_allowed) Xmin = Xmin_allowed
    If (Xmax > Xmax_allowed) Int_Errorflag = Int_Errorflag + 1
    If (Xmax > Xmax_allowed) Xmax = Xmax_allowed

    If (Int_Errorflag == 1) Int_errorflag = 0   !Just one wrong starting value could be adapted, continue in a regular way
    If (Int_Errorflag == 2) Int_Errorflag = 1   !In case both values were wrong attempt to use sub intervals

    ! Check whether there is exactly one root between Xmin and Xmax

    IF (Int_Errorflag == 0) F_Xmin = Zerofunction(gl,Xmin, t_params)
    !Check whether the zerofunction returns an errorvalue (parameters(65) != 0?)
    if (t_params%a_p(65) < -1000.D0) then
        errorflag = int(t_params%a_p(65))
        return
    end if
    IF (Int_Errorflag == 0) F_Xmax = Zerofunction(gl,Xmax, t_params)
    !Check whether the zerofunction returns an errorvalue (parameters(65) != 0?)
    if (t_params%a_p(65) < -1000.D0) then
        errorflag = int(t_params%a_p(65))
        return
    end if

    If ((F_Xmin * F_Xmax) < 0.D0) Root_In = .True.

    ! If no root is found between Xstart_min and Xstart_max, it is assumed that the interval was chosen too small
    ! The interval is enlarged in steps of (Xstart_min + XStart_max)/2 to each side until a root is included
    ! The process is repeated up to ten times or until a root is included or the allowable limits of X are exceeded

    Do While ((.Not.Root_In).And.(Int_Errorflag == 0).And.(Count_enlarge <= 10))

        Xmin_alt = Xmin
        Xmax_alt = Xmax
        F_Xmin_alt = F_Xmin
        F_Xmax_alt = F_Xmax

        If (Stop_min == 0) Xmin = Xmin - (Xstart_max - Xstart_min) * 0.5d0
        If (Stop_max == 0) Xmax = Xmax + (Xstart_max - Xstart_min) * 0.5d0
        Count_enlarge = Count_enlarge + 1

        ! Check, whether Xmin and Xmax get out of the allowed range

        !If (Xmin < Xmin_allowed) Int_Errorflag = 1
        !If (Xmax > Xmax_allowed) Int_Errorflag = 1
        If ((Xmin - Xmin_allowed) < -1e-11) Int_Errorflag = 1
        If ((Xmax - Xmax_allowed) > 1e-11) Int_Errorflag = 1
        !
        If ((Int_Errorflag == 0).AND.(Stop_min == 0))F_Xmin = Zerofunction(gl,Xmin, t_params)
        !Check whether the zerofunction returns an errorvalue (parameters(65) != 0?)
        if (t_params%a_p(65) < -1000.D0) then
            errorflag = int(t_params%a_p(65))
            return
        end if
        If ((Int_Errorflag == 0).AND.(Stop_max == 0))F_Xmax = Zerofunction(gl,Xmax, t_params)
        !Check whether the zerofunction returns an errorvalue (parameters(65) != 0?)
        if (t_params%a_p(65) < -1000.D0) then
            errorflag = int(t_params%a_p(65))
            return
        end if

        ! Reset enlargments of the intervall on the side which makes the results worse

        If ((DABS(F_Xmin_alt) < DABS(F_Xmin)).AND.((F_Xmin_alt*F_Xmin) > 0.d0)) Stop_min = 1
        If (Stop_min == 1) Xmin = Xmin_alt
        If (Stop_min == 1) F_Xmin = F_Xmin_alt

        If ((DABS(F_Xmax_alt) < DABS(F_Xmax)).AND.((F_Xmax_alt*F_Xmax) > 0.d0)) Stop_max = 1
        If (Stop_max == 1) Xmax = Xmax_alt
        If (Stop_max == 1) F_Xmax = F_Xmax_alt

        ! Check, whether a root is within the new limits
        If ((F_Xmin * F_Xmax) < 0.D0) Root_In = .True.

    End Do ! End of the enlargement loop

    ! Reset the internal Errorflag

    Int_Errorflag = 0

    ! If no root was found by enlarging the interval, there may be two roots within the initial interval.
    ! To check out this posibility, the original interval is subdivided starting in the midle of the original interval.
    ! The size of the new interval starts with 1/10 of the original interval and is enlarged then.

    Do While ((.Not.Root_In).And.(Int_Errorflag == 0).And.(Count_reduce < 9))
        if(iterations == (Max_Iterations - 2))then
            continue
        end if
        Count_reduce = Count_reduce + 1
        Xmin = (Xstart_max + Xstart_min)/2.D0 - (Xstart_max - Xstart_min)/20.D0 * Count_reduce
        Xmax = (Xstart_max + Xstart_min)/2.D0 + (Xstart_max - Xstart_min)/20.D0 * Count_reduce

        ! Check, whether Xmin and Xmax get out of the allowed range

        If (Xmin < Xmin_allowed) Int_Errorflag = 1
        If (Xmax > Xmax_allowed) Int_Errorflag = 1

        If (Int_Errorflag == 0) F_Xmin = Zerofunction(gl,Xmin, t_params)
        If (Int_Errorflag == 0) F_Xmax = Zerofunction(gl,Xmax, t_params)

        ! Check, whether a root is within the new limits

        If ((F_Xmin * F_Xmax) < 0.D0) Root_In = .True.

    End Do  ! End of the reduction loop

    ! If no root could be found in the original interval, by enlarging or by reducing of the interval the
    ! routine returns with Errorflag = 2

    If (.NOT.Root_In) Then
        Errorflag = 2
        Return
    End If

    ! A single root has been identified between the current values of Xmin and Xmax.
    ! A modified Regula Falsi algorithm is started to find the root.

    Do While ((.NOT.Root_Found).and.(Iterations <= Max_Iterations))

        Iterations = Iterations + 1   ! Count iteration steps

        Deriv = (F_Xmax - F_Xmin) / (Xmax - Xmin)  ! Calculate the "derivative" for the current interval
        if (abs(deriv) > 1.d-12) then
            Deviation_X = F_Xakt / Deriv !Deviation_X is the distance to the linear approximated root from the akt position on the x-axis
            !SH 10/2018: If derivative is large make additional checks (see line 300) if root found or not
            if (abs(deriv) > 1.d10) then
                additional_checks = .true.
            else
                additional_checks = .false.
            endif
        else
            Deriv = sign(1.d-12,Deriv)
        end if

        If (DABS(F_Xmin) < DABS(F_Xmax)) then   ! Calculate Xakt starting from the X value with the smaller function value

            If ((DABS(DERIV) > 0).And.(MOD(Iterations,5) /= 0)) Then          ! Calculate Xakt starting from the lower limit of the interval
                Xakt = Xmin + Dabs(F_Xmin / Deriv)                                 ! Use this algorithm only three out of four times to avoid problems in regions
            ELSE                                                                ! with strong curvature
                Xakt = (Xmin + Xmax)/2.D0
            End If

        ELSE

            If ((DABS(DERIV) > 0).And.(MOD(Iterations,5) /= 0)) Then          ! Calculate Xakt starting from the upper limit of the interval
                Xakt = Xmax - Dabs(F_Xmax / Deriv)                                 ! Use this algorithm only three out of four times to avoid problems in regions
            ELSE                                                                ! with strong curvature
                Xakt = (Xmin + Xmax)/2.D0
            End If

        END IF

        If ((Xakt >= Xmax).OR.(Xakt <= Xmin)) Xakt = (Xmin + Xmax)/ 2.D0    ! In case Xakt gets out of the interval where
        !                                                                       the root has been found the interval is just split

        F_Xakt = Zerofunction(gl,Xakt, t_params)
        !Check whether the zerofunction returns an errorvalue (parameters(65) != 0?)
        if (t_params%a_p(65) < -1000.D0) then
            errorflag = int(t_params%a_p(65))
            return
        end if

        ! Check whether the root has been found
        ! The criterion Delta_allowed is applied to the deviation in X

        if (abs(deriv) > 1.d-12) then
            Deviation_X = F_Xakt / Deriv
            !SH 10/2018: If derivative is large make additional checks (see line 300) if root found or not
            if (abs(deriv) > 1.d10) then
                additional_checks = .true.
            else
                additional_checks = .false.
            endif
        else
            Deriv = sign(1.d-12,Deriv)
        end if

        ! In the Maxwell iteration the difference in F_Xakt has to become smaller than the defined limit.
        ! At the appropriate point in VLEpure the Parameter(65) is set to -999.D0 to indicate this difference.

        If (Dabs(t_params%a_p(65)+ 999.D0) < 1.D-6) Deviation_X = F_Xakt

        !SH & TN: old code:
        !If (DABS(Deviation_X) < Delta_Allowed) Then
        !new code:
        !for very large oscillations Deviation_X can gain very small values although xmin and xmax deviate from the real root,
        !thus the relative deviations are checked if additional_parameters
        If ( ((DABS(Deviation_X) < Delta_Allowed) .and. (.not. additional_checks)) .or. &
            & ((DABS(Deviation_X) < Delta_Allowed) .and. (additional_checks) .and. ((dabs((xakt-xmin)/xakt) < Delta_Allowed) .or. (dabs((xakt-xmax)/xakt) < Delta_Allowed))) ) Then
            Root_found = .true.
            Xroot = Xakt + Deviation_X  ! Xroot usually is a better guess for the root then Xakt

        Else  ! Root has not yet been found, new interval has to be defined

            ! If Deviation_X is already smaller than (Xmax-Xmin)/10 it may be worth to calculate new boundaries

            If (DABS(Deviation_X) < ((Xmax-Xmin)/10.D0)) Then

                ! Use the numerical "derivative" between the pair of points that is closer together

                If (((Xakt-Xmin) < (Xmax-Xakt)).And.((Xakt-Xmin) > 0.D0)) Deriv = (F_Xakt-F_Xmin)/(Xakt-Xmin)
                If (((Xakt-Xmin) > (Xmax-Xakt)).And.((Xmax-Xakt) > 0.D0)) Deriv = (F_Xmax-F_Xakt)/(Xmax-Xakt)
                ! check if the numerical derivative becomes zero
                if (abs(deriv) > 1.d-12) then
                    Deviation_X = F_Xakt / Deriv
                else
                    Deriv = sign(1.d-12,Deriv)
                end if
                IF((F_Xmin*F_Xakt) < 0.D0) Then ! Root is located between Xmin and Xakt
                    Xtest = Xakt - 2.D0*DABS(Deviation_X) ! If Zerofunction is linear, the root is in the middle of the intervall now

                    If (Xtest < ((Xmin+Xakt)/2.D0)) Xtest = (Xmin+Xakt)/2.D0  ! Avoid disadvantageous redefinition of the interval
                    ! that is possible in cases with strong curvature
                    F_Xtest = Zerofunction(gl,Xtest, t_params)
                    !Check whether the zerofunction returns an errorvalue (parameters(65) != 0?)
                    if (t_params%a_p(65) < -1000.D0) then
                        errorflag = int(t_params%a_p(65))
                        return
                    end if
                    If ((F_Xtest*F_Xakt) < 0.D0) Xmin = Xtest ! Use Xtest as lower limit, if the root is between Xtest and Xakt
                    If ((F_Xtest*F_Xakt) < 0.D0) F_Xmin = F_Xtest

                Else IF((F_Xmax*F_Xakt) < 0.D0) Then  ! Root is located between Xakt and Xmax
                    Xtest = Xakt + 2.D0*DABS(Deviation_X)

                    If (Xtest > ((Xmax+Xakt)/2.D0)) Xtest = (Xmax+Xakt)/2.D0  ! Avoid disadvantageous redefinition of the interval
                    ! that is possible in cases with strong curvature
                    F_Xtest = Zerofunction(gl,Xtest, t_params)
                    !Check whether the zerofunction returns an errorvalue (parameters(65) != 0?)
                    if (t_params%a_p(65) < -1000.D0) then
                        errorflag = int(t_params%a_p(65))
                        return
                    end if
                    If ((F_Xtest*F_Xakt) < 0.D0) Xmax = Xtest
                    If ((F_Xtest*F_Xakt) < 0.D0) F_Xmax = F_Xtest

                End IF

            End IF  ! End of the loop possibly minimizing the size of the currently used interval

            ! Defining the new interval

            IF((F_Xmin*F_Xakt) < 0.D0) Then  ! Root is located between Xmin and Xakt
                Xmax = Xakt
                F_Xmax = F_Xakt
            Else IF((F_Xmax*F_Xakt) < 0.D0) Then   ! Root is located between Xakt and Xmax
                Xmin = Xakt
                F_Xmin = F_Xakt
            End If

        End If ! Redefinition of the interval completed

    End Do ! Iteration loop completed

    ! Interval needed to be enlarged or reduced

    IF ((Count_enlarge + Count_reduce) > 0) Errorflag = 3

    ! No root was fond within Max_Iterations

    IF(.NOT.Root_found) Errorflag = 1

    Return
    End




    end module module_regula_falsi
