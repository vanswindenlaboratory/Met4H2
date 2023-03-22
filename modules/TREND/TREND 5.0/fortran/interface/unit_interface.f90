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
    !!          Neumann, T.; Pohl S. M.; Semrau, B.; Thol, M. (2020):
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

    ! --------------------------------------------------------
    !Module for getting the correct unit for a given proptype
    ! Sven Pohl, 29.09.2020
    ! --------------------------------------------------------
    module prop_unit_m

    use utility_module
    ! use unit_convertion

    implicit none

    integer, parameter :: size_arr = 125

    type char
        character(:) , allocatable :: n
        integer :: id
    end type

    ! Type declaration
    type list_t
        type(char) , dimension(size_arr,3) :: cnt
    contains
    procedure ,pass :: create_list
    end type

    type conv_t
        type(char) , dimension(size_arr) :: cnt
    contains
    procedure ,pass :: create_conv
    end type


    !Interface for constructor
    interface list_t
    procedure :: create_list
    end interface list_t
    interface conv_t
    procedure :: create_conv
    end interface conv_t

    contains

    ! Function for transfer unit_char to unit_integer
    integer function get_unitin_m(unit_char)

    character(len=*) :: unit_char

    if(trim(unit_char) == "specific") then
        get_unitin_m=2
    elseif(trim(unit_char) == "molar") then
        get_unitin_m=1
    elseif(trim(unit_char) == "reduced") then
        get_unitin_m=3
    else
        get_unitin_m = 0
    end if
    end function

    ! Create list with units
    logical function create_list(this)
    class(list_t) :: this

    this%cnt(1,1)%n='K';					this%cnt(1,2)%n='K';					this%cnt(1,3)%n='-'
    this%cnt(2,1)%n='mol/m3';				this%cnt(2,2)%n='kg/m3';				this%cnt(2,3)%n='-'
    this%cnt(3,1)%n='MPa';					this%cnt(3,2)%n='MPa';					this%cnt(3,3)%n='-'
    this%cnt(4,1)%n='J/mol';				this%cnt(4,2)%n='J/kg';					this%cnt(4,3)%n='-'
    this%cnt(5,1)%n='J/(mol K)';			this%cnt(5,2)%n='J/(kg K)';				this%cnt(5,3)%n='-'
    this%cnt(6,1)%n='m/s';					this%cnt(6,2)%n='m/s';					this%cnt(6,3)%n='-'
    this%cnt(7,1)%n='m6/mol2';				this%cnt(7,2)%n='m6/kg2';				this%cnt(7,3)%n='-'
    this%cnt(8,1)%n='-';					this%cnt(8,2)%n='-';					this%cnt(8,3)%n='-'
    this%cnt(9,1)%n='mol/mol';				this%cnt(9,2)%n='kg/kg';				this%cnt(9,3)%n='-'
    this%cnt(10,1)%n='mN/m';				this%cnt(10,2)%n='mN/m';				this%cnt(10,3)%n='-'
    this%cnt(11,1)%n='µPa s';				this%cnt(11,2)%n='µPa s';				this%cnt(11,3)%n='-'
    this%cnt(12,1)%n='W/(m K)';				this%cnt(12,2)%n='W/(m K)';				this%cnt(12,3)%n='-'
    this%cnt(13,1)%n='m9/mol3';				this%cnt(13,2)%n='m9/kg3';				this%cnt(13,3)%n='-'
    this%cnt(14,1)%n='1/MPa';				this%cnt(14,2)%n='1/MPa';				this%cnt(14,3)%n='-'
    this%cnt(15,1)%n='m3/mol';				this%cnt(15,2)%n='m3/kg';				this%cnt(15,3)%n='-'
    this%cnt(16,1)%n='MPa/K';				this%cnt(16,2)%n='MPa/K';				this%cnt(16,3)%n='-'
    this%cnt(17,1)%n='m3/(mol K)';			this%cnt(17,2)%n='m3/(kg K)';			this%cnt(17,3)%n='-'
    this%cnt(18,1)%n='m6/(mol2 K)';			this%cnt(18,2)%n='m6/(kg2 K)';			this%cnt(18,3)%n='-'
    this%cnt(19,1)%n='(MPa m3)/mol';		this%cnt(19,2)%n='(MPa m3)/kg';			this%cnt(19,3)%n='-'
    this%cnt(20,1)%n='MPa/(mol/m3)2';		this%cnt(20,2)%n='MPa/(kg/m3)2';		this%cnt(20,3)%n='-'
    this%cnt(21,1)%n='(MPa m3)/(mol K)';	this%cnt(21,2)%n='(MPa m3)/(kg K)';		this%cnt(21,3)%n='-'
    this%cnt(22,1)%n='MPa/K2';				this%cnt(22,2)%n='MPa/K2';				this%cnt(22,3)%n='-'
    this%cnt(23,1)%n='1/K';					this%cnt(23,2)%n='1/K';					this%cnt(23,3)%n='-'
    this%cnt(24,1)%n='J/(mol m3)';			this%cnt(24,2)%n='J/(kg m3)';			this%cnt(24,3)%n='-'
    this%cnt(25,1)%n='mol/(m3 K)';			this%cnt(25,2)%n='kg/(m3 K)';			this%cnt(25,3)%n='-'
    this%cnt(26,1)%n='K/MPa';				this%cnt(26,2)%n='K/MPa';				this%cnt(26,3)%n='-'
    this%cnt(27,1)%n='kg/mol';				this%cnt(27,2)%n='kg/mol';				this%cnt(27,3)%n='-'

    create_list = .true.

    end

    ! Create list with conversions
    logical function create_conv(this)
    class(conv_t) :: this

    this%cnt(1)%n = "t";          this%cnt(21)%n="a02";        this%cnt(41)%n="dvir";         this%cnt(61)%n="d2pddt";      this%cnt(81)%n="tjtinv";           this%cnt(101)%n="pvap";     this%cnt(121)%n= 'x20';
    this%cnt(1)%id =1;            this%cnt(21)%id=8;           this%cnt(41)%id=13;            this%cnt(61)%id=21;           this%cnt(81)%id=1;                 this%cnt(101)%id=3;         this%cnt(121)%id= 9;
    this%cnt(2)%n="d";            this%cnt(22)%n="a03";        this%cnt(42)%n="riem";         this%cnt(62)%n="d2pdt2";      this%cnt(82)%n="ttrip";            this%cnt(102)%n=  'x1';     this%cnt(122)%n= 'nh';
    this%cnt(2)%id=2;             this%cnt(22)%id=8;           this%cnt(42)%id=8;             this%cnt(62)%id=22;           this%cnt(82)%id=1;                 this%cnt(102)%id= 9;        this%cnt(122)%id= 8;
    this%cnt(3)%n="p";            this%cnt(23)%n="a10";        this%cnt(43)%n="compt";        this%cnt(63)%n="dddt";        this%cnt(83)%n="ptrip";            this%cnt(103)%n=  'x2';     this%cnt(123)%n= 'theta_s^o';
    this%cnt(3)%id=3;             this%cnt(23)%id=8;           this%cnt(43)%id=14;            this%cnt(63)%id=25;           this%cnt(83)%id=3;                 this%cnt(103)%id= 9;        this%cnt(123)%id= 8;
    this%cnt(4)%n="u";            this%cnt(24)%n="a11";        this%cnt(44)%n="comps";        this%cnt(64)%n="dspin";       this%cnt(84)%n="tcrit";            this%cnt(104)%n= 'x3';      this%cnt(124)%n= 'theta_l^o';
    this%cnt(4)%id=4;             this%cnt(24)%id=8;           this%cnt(44)%id=14;            this%cnt(64)%id=2;            this%cnt(84)%id=1;                 this%cnt(104)%id= 9;        this%cnt(124)%id= 8;
    this%cnt(5)%n="h";            this%cnt(25)%n="a12";        this%cnt(45)%n="throt";        this%cnt(65)%n="phase";       this%cnt(85)%n="pcrit";            this%cnt(105)%n= 'x4';      this%cnt(125)%n= 'beta';
    this%cnt(5)%id=4;             this%cnt(25)%id=8;           this%cnt(45)%id=15;            this%cnt(65)%id=8;            this%cnt(85)%id=3;                 this%cnt(105)%id= 9;        this%cnt(125)%id= 9;
    this%cnt(6)%n="s";            this%cnt(26)%n="a20";        this%cnt(46)%n="exps";         this%cnt(66)%n="mw";          this%cnt(86)%n="dcrit";            this%cnt(106)%n= 'x5';
    this%cnt(6)%id=5;             this%cnt(26)%id=8;           this%cnt(46)%id=8;             this%cnt(66)%id=27;           this%cnt(86)%id=2;                 this%cnt(106)%id= 9;
    this%cnt(7)%n="g";            this%cnt(27)%n="a21";        this%cnt(47)%n="expt";         this%cnt(67)%n="ttp";         this%cnt(87)%n="acentricfac";      this%cnt(107)%n= 'x6';
    this%cnt(7)%id=4;             this%cnt(27)%id=8;           this%cnt(47)%id=8;             this%cnt(67)%id=1;            this%cnt(87)%id=8;                 this%cnt(107)%id= 9;
    this%cnt(8)%n="a";            this%cnt(28)%n="a30";        this%cnt(48)%n="dudv";         this%cnt(68)%n="ptp";         this%cnt(88)%n="af";               this%cnt(108)%n= 'x7';
    this%cnt(8)%id=4;             this%cnt(28)%id=8;           this%cnt(48)%id=24;            this%cnt(68)%id=3;            this%cnt(88)%id=8;                 this%cnt(108)%id= 9;
    this%cnt(9)%n="cp";           this%cnt(29)%n="ur";         this%cnt(49)%n="jtco";         this%cnt(69)%n="tceos";       this%cnt(89)%n="tmin";             this%cnt(109)%n= 'x8';
    this%cnt(9)%id=5;             this%cnt(29)%id=4;           this%cnt(49)%id=26;            this%cnt(69)%id=1;            this%cnt(89)%id=1;                 this%cnt(109)%id= 9;
    this%cnt(10)%n="cv";          this%cnt(30)%n="hr";         this%cnt(50)%n="dpdt";         this%cnt(70)%n="pceos";       this%cnt(90)%n="tmax";             this%cnt(110)%n= 'x9';
    this%cnt(10)%id=5;            this%cnt(30)%id=4;           this%cnt(50)%id=16;            this%cnt(70)%id=3;            this%cnt(90)%id=1;                 this%cnt(110)%id= 9;
    this%cnt(11)%n="ws";          this%cnt(31)%n="sr";         this%cnt(51)%n="ve";           this%cnt(71)%n="dceos";       this%cnt(91)%n="pmin";             this%cnt(111)%n= 'x10';
    this%cnt(11)%id=6;            this%cnt(31)%id=5;           this%cnt(51)%id=15;            this%cnt(71)%id=2;            this%cnt(91)%id=3;                 this%cnt(111)%id= 9;
    this%cnt(12)%n="bvir";        this%cnt(32)%n="cvr";        this%cnt(52)%n="he";           this%cnt(72)%n="afeos";       this%cnt(92)%n="pmax";             this%cnt(112)%n= 'x11';
    this%cnt(12)%id=15;           this%cnt(32)%id=5;           this%cnt(52)%id=4;             this%cnt(72)%id=8;            this%cnt(92)%id=3;                 this%cnt(112)%id= 9;
    this%cnt(13)%n="cvir";        this%cnt(33)%n="cpotr";      this%cnt(53)%n="ge";           this%cnt(73)%n="tmineos";     this%cnt(93)%n="dmax";             this%cnt(113)%n= 'x12';
    this%cnt(13)%id=7;            this%cnt(33)%id=4;           this%cnt(53)%id=4;             this%cnt(73)%id=1;            this%cnt(93)%id=2;                 this%cnt(113)%id= 9;
    this%cnt(14)%n="cp0";         this%cnt(34)%n="gruen";      this%cnt(54)%n="b12";          this%cnt(74)%n="tmaxeos";     this%cnt(94)%n="casnr";            this%cnt(114)%n= 'x13';
    this%cnt(14)%id=5;            this%cnt(34)%id=8;           this%cnt(54)%id=15;            this%cnt(74)%id=1;            this%cnt(94)%id=8;                 this%cnt(114)%id= 9;
    this%cnt(15)%n="q";           this%cnt(35)%n="pip";        this%cnt(55)%n="gammagd";      this%cnt(75)%n="pmineos";     this%cnt(95)%n="tmelt";            this%cnt(115)%n= 'x14';
    this%cnt(15)%id=9;            this%cnt(35)%id=8;           this%cnt(55)%id=8;             this%cnt(75)%id=3;            this%cnt(95)%id=1;                 this%cnt(115)%id= 9;
    this%cnt(16)%n="z";           this%cnt(36)%n="de";         this%cnt(56)%n="pr";           this%cnt(76)%n="pmaxeos";     this%cnt(96)%n="pmelt";            this%cnt(116)%n= 'x15';
    this%cnt(16)%id=8;            this%cnt(36)%id=8;           this%cnt(56)%id=3;             this%cnt(76)%id=3;            this%cnt(96)%id=3;                 this%cnt(116)%id= 9;
    this%cnt(17)%n="pnum";        this%cnt(37)%n="st";         this%cnt(57)%n="dbdt";         this%cnt(77)%n="dmaxeos";     this%cnt(97)%n="tsub";             this%cnt(117)%n= 'x16';
    this%cnt(17)%id=3;            this%cnt(37)%id=10;          this%cnt(57)%id=17;            this%cnt(77)%id=2;            this%cnt(97)%id=1;                 this%cnt(117)%id= 9;
    this%cnt(18)%n="wsnum";       this%cnt(38)%n="eta";        this%cnt(58)%n="dcdt";         this%cnt(78)%n="casfromfluid";this%cnt(98)%n="psub";             this%cnt(118)%n= 'x17';
    this%cnt(18)%id=6;            this%cnt(38)%id=11;          this%cnt(58)%id=18;            this%cnt(78)%id=8;            this%cnt(98)%id=3;                 this%cnt(118)%id= 9;
    this%cnt(19)%n="a00";         this%cnt(39)%n="tcx";        this%cnt(59)%n="dpdd";         this%cnt(79)%n="tboyle";      this%cnt(99)%n="dliq";             this%cnt(119)%n= 'x18';
    this%cnt(19)%id=8;            this%cnt(39)%id=12;          this%cnt(59)%id=19;            this%cnt(79)%id=1;            this%cnt(99)%id=2;                 this%cnt(119)%id= 9;
    this%cnt(20)%n="a01";         this%cnt(40)%n="vexp";       this%cnt(60)%n="d2pdd2";       this%cnt(80)%n="tjinv";       this%cnt(100)%n="dvap";            this%cnt(120)%n= 'x19';
    this%cnt(20)%id=8;            this%cnt(40)%id=23;          this%cnt(60)%id=20;            this%cnt(80)%id=1;            this%cnt(100)%id=2;                this%cnt(120)%id= 9;

    create_conv  = .true.
    end



    ! Main function here
    function PROP_UNIT_INTERFACE(unitdefinition_arg, unittype_arg_in, size_char, errorflag)

    type(list_t)  :: unit_list
    type(conv_t)  :: conver
    character(*) :: unitdefinition_arg,unittype_arg_in
    character(:) , allocatable :: unittype_arg
    integer :: errorflag,id_conv,id_unit,loop_unit
    integer :: size_char, unitinput
    character(size_char)  ::  PROP_UNIT_INTERFACE
    logical :: a,b

    unittype_arg = unittype_arg_in



    if(len_trim(unitdefinition_arg) > 0) then

        call uppertolower_char(unittype_arg,len(unittype_arg))

        call uppertolower_char(unitdefinition_arg,len(unitdefinition_arg))

        unitinput = check_valid_input(unitdefinition_arg)
        if(unitinput < 0) then
            errorflag = unitinput
            return
        end if

        a =  unit_list%create_list()

        b =  conver%create_conv()

        id_unit = get_unitin_m(unitdefinition_arg)

        id_conv = 0
        if(len_trim(unittype_arg) > 0) then
            do loop_unit=1,size(conver%cnt)
                if(trim(unittype_arg) .eq. conver%cnt(loop_unit)%n) then
                    id_conv = conver%cnt(loop_unit)%id
                    exit
                end if
            end do
        end if

        if(id_conv == 0 .or. id_unit == 0) then
            PROP_UNIT_INTERFACE = 'TREND Unit error'
            errorflag = -123457
        else
            PROP_UNIT_INTERFACE = unit_list%cnt(id_conv,id_unit)%n
            errorflag = 0
        end if

    else ! unitdefinition_arg is empty
        PROP_UNIT_INTERFACE = '-9935'
        errorflag  = -9935
    end if


    end function PROP_UNIT_INTERFACE

    !*******************************************************************************************
    integer function check_valid_input(unit_char)
    character(*), intent(in) :: unit_char

    if(trim(unit_char) == "specific") then
        check_valid_input = 1
    elseif(trim(unit_char) == "molar") then
        check_valid_input = 1
    elseif(trim(unit_char) == "reduced") then
        check_valid_input = 1
    else
        check_valid_input = -123456
        return
    end if

    end function check_valid_input
    !*********************************************************************************************


    end module prop_unit_m



    module calctype_setter_m

    use module_all_types
    use utility_module

    implicit none

    integer, parameter :: nr_props = 85

    type calctype_t

        integer, dimension(nr_props)  :: calctype_internal =  (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, &
            27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65, &
            66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85 /)

        character(10) , dimension(nr_props)  :: calctypes = (/"t","d","p","u","h","s","g","a","cp","cv","ws","bvir","cvir","cp0","q","z","pnum","wsnum", &
            "a00","a01","a02","a03","a10","a11","a12","a20","a21","a30","ur","hr","sr","cvr","cpotr","gruen","pip","de","st","eta","tcx","vexp","dvir", &
            "riem","compt","comps","throt","exps","expt","dudv","jtco","dpdt","ve","he","ge","b12","gammagd","pr","dbdt","dcdt","dpdd","d2pdd2","d2pddt", &
            "d2pdt2","dddt","dspin","phase","af","tdensmax","tboyle","tjinv","tjtinv","dpipdd","d2pipdd2","d2bdt2","d2cdt2","d2ddt2","d3bdt3","d3cdt3","d3ddt3",  &
            "d4bdt4","d4cdt4","d4ddt4","b_calc_num","c_calc_num","d_calc_num","neff"/)

        integer , dimension(nr_props)  :: converttype = (/0,2,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,5,0,0,0,3,0,0,1,0,0,3,1,1,3,0, &
            0,3,4,3,4,3,0,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

        logical , dimension(nr_props)  :: virial = (/.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.true. ,.true. ,.false.,.false., &
            .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
            .false.,.false.,.false.,.false.,.false.,.false.,.true. ,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
            .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
            .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)

        integer , dimension(nr_props) :: transport = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)

        logical , dimension(nr_props) :: ref_state = (/.false.,.false.,.false.,.true. ,.true. ,.true. ,.true. ,.true. ,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
            .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
            .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.true. ,.true. ,.false.,.false.,.false.,.false.,.false., &
            .false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
            .false.,.false.,.false.,.false.,.false.,.false./)


    end type

    contains

    subroutine set_calctype_to_gl(gl_in,calctype_char,calctype_internal)

    type(type_gl) :: gl_in
    integer :: index,calctype_internal
    character(*) :: calctype_char
    type(calctype_t) :: c

    call uppertolower_char(calctype_char,len(calctype_char))

    index = findloc(c%calctypes(:),trim(calctype_char),1)

    if(index == 0) then

    end if

    gl_in%calc_ref =    c%ref_state(index)
    gl_in%converttype = c%converttype(index)
    gl_in%vir = c%virial(index)
    calctype_internal = c%calctype_internal(index)
    gl_in%transport =c%transport(index)

    end subroutine



    end module

