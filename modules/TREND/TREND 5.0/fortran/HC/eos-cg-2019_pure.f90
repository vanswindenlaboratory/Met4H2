! This file was generated -  Dont edit!
module eos_cg_2019_pure_fld_files
! #######################################################
character(256), TARGET :: argon(884)
DATA argon(1) / 'argon              !short name' /
DATA argon(2) / '7440-37-1          !CAS number' /
DATA argon(3) / 'argon              !full name' /
DATA argon(4) / 'Ar                 !chemical formula' /
DATA argon(5) / 'R-740              !synonym' /
DATA argon(6) / '39.948             !molecular weight [g/mol]' /
DATA argon(7) / '83.8058            !triple point temperature [K]' /
DATA argon(8) / '87.302             !normal boiling point [K]' /
DATA argon(9) / '150.687            !critical temperature [K]' /
DATA argon(10) / '4863.0             !critical pressure [kPa]' /
DATA argon(11) / '13.4074            !critical density [mol/L]' /
DATA argon(12) / '-0.00219           !acentric factor' /
DATA argon(13) / '0.0                !dipole moment [Debye]' /
DATA argon(14) / 'OT0                !default reference state' /
DATA argon(15) / '298.15  101.325  6197.0  154.737  !tref, Pref, Href, Sref' /
DATA argon(16) / '8.0                !version number' /
DATA argon(17) / '1951               !UN Number' /
DATA argon(18) / 'other              !family' /
DATA argon(19) / '0.0                !heating value (gross or superior) [kJ/mol]' /
DATA argon(20) / 'A1                 !Safety Group (ASHRAE Standard 34, 2010)' /
DATA argon(21) / '' /
DATA argon(22) / '' /
DATA argon(23) / '' /
DATA argon(24) / '#EOS               !equation of state specification' /
DATA argon(25) / 'FEQ  Helmholtz equation of state for argon of Tegeler et al. (1999).' /
DATA argon(26) / '?LITERATURE REFERENCE \' /
DATA argon(27) / '?Tegeler, Ch., Span, R., and Wagner, W.,' /
DATA argon(28) / '? "A New Equation of State for Argon Covering the Fluid Region for' /
DATA argon(29) / '? Temperatures from the Melting Line to 700 K at Pressures up to 1000 MPa,"' /
DATA argon(30) / '? J. Phys. Chem. Ref. Data, 28(3):779-850, 1999.' /
DATA argon(31) / '?\' /
DATA argon(32) / '?The estimated uncertainty in density is less than 0.02% for pressures up' /
DATA argon(33) / '?to 12 MPa and temperatures up to 340 K with the exception of the' /
DATA argon(34) / '?critical region and less than 0.03% for pressures up to 30 MPa and' /
DATA argon(35) / '?temperatures between 235 and 520 K. Elsewhere, the uncertainty in' /
DATA argon(36) / '?density is generally within 0.2%.  In the region with densities up to' /
DATA argon(37) / '?half the critical density and for temperatures between 90 and 450 K, the' /
DATA argon(38) / '?estimated uncertainty of calculated speeds of sound is in general less' /
DATA argon(39) / '?than 0.02%.  In the liquid and supercritical regions, the uncertainty is' /
DATA argon(40) / '?less than 1%.  The uncertainty in heat capacities is within 0.3% for the' /
DATA argon(41) / '?vapor and 2% for the liquid.  The formulation gives reasonable' /
DATA argon(42) / '?extrapolation behavior up to very high pressures (50 GPa) and' /
DATA argon(43) / '?temperatures (17000 K).' /
DATA argon(44) / '?\' /
DATA argon(45) / '!end of info section' /
DATA argon(46) / '83.8058            !lower temperature limit [K]' /
DATA argon(47) / '2000.0             !upper temperature limit [K]' /
DATA argon(48) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(49) / '50.65              !maximum density [mol/L]' /
DATA argon(50) / 'CPP                                    !pointer to Cp0 model' /
DATA argon(51) / '39.948                                 !molecular weight [g/mol]' /
DATA argon(52) / '83.8058                                !triple point temperature [K]' /
DATA argon(53) / '68.891                                 !pressure at triple point [kPa]' /
DATA argon(54) / '35.465                                 !density at triple point [mol/L]' /
DATA argon(55) / '87.302                                 !normal boiling point temperature [K]' /
DATA argon(56) / '-0.00219                               !acentric factor' /
DATA argon(57) / '150.687      4863.0       13.40742965  !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA argon(58) / '150.687                   13.40742965  !reducing parameters [K, mol/L]' /
DATA argon(59) / '8.31451                                !gas constant [J/mol-K]' /
DATA argon(60) / '      37  4      4  12      0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA argon(61) / ' 0.887223049900d-01  0.000   1.00    0 !a(i),t(i),d(i),l(i)' /
DATA argon(62) / ' 0.705148051673d+00  0.250   1.00    0' /
DATA argon(63) / '-0.168201156541d+01  1.000   1.00    0' /
DATA argon(64) / '-0.149090144315d+00  2.750   1.00    0' /
DATA argon(65) / '-0.120248046009d+00  4.000   1.00    0' /
DATA argon(66) / '-0.121649787986d+00  0.000   2.00    0' /
DATA argon(67) / ' 0.400359336268d+00  0.250   2.00    0' /
DATA argon(68) / '-0.271360626991d+00  0.750   2.00    0' /
DATA argon(69) / ' 0.242119245796d+00  2.750   2.00    0' /
DATA argon(70) / ' 0.578895831856d-02  0.000   3.00    0' /
DATA argon(71) / '-0.410973356153d-01  2.000   3.00    0' /
DATA argon(72) / ' 0.247107615416d-01  0.750   4.00    0' /
DATA argon(73) / '-0.321813917507d+00  3.000   1.00    1' /
DATA argon(74) / ' 0.332300176958d+00  3.500   1.00    1' /
DATA argon(75) / ' 0.310199862873d-01  1.000   3.00    1' /
DATA argon(76) / '-0.307770860024d-01  2.000   4.00    1' /
DATA argon(77) / ' 0.938911374196d-01  4.000   4.00    1' /
DATA argon(78) / '-0.906432106820d-01  3.000   5.00    1' /
DATA argon(79) / '-0.457783492767d-03  0.000   7.00    1' /
DATA argon(80) / '-0.826597290252d-04  0.500  10.00    1' /
DATA argon(81) / ' 0.130134156031d-03  1.000  10.00    1' /
DATA argon(82) / '-0.113978400020d-01  1.000   2.00    2' /
DATA argon(83) / '-0.244551699605d-01  7.000   2.00    2' /
DATA argon(84) / '-0.643240671760d-01  5.000   4.00    2' /
DATA argon(85) / ' 0.588894710937d-01  6.000   4.00    2' /
DATA argon(86) / '-0.649335521130d-03  6.000   8.00    2' /
DATA argon(87) / '-0.138898621584d-01 10.000   3.00    3' /
DATA argon(88) / ' 0.404898392969d+00 13.000   5.00    3' /
DATA argon(89) / '-0.386125195947d+00 14.000   5.00    3' /
DATA argon(90) / '-0.188171423322d+00 11.000   6.00    3' /
DATA argon(91) / ' 0.159776475965d+00 14.000   6.00    3' /
DATA argon(92) / ' 0.539855185139d-01  8.000   7.00    3' /
DATA argon(93) / '-0.289534179580d-01 14.000   7.00    3' /
DATA argon(94) / '-0.130254133814d-01  6.000   8.00    3' /
DATA argon(95) / ' 0.289486967758d-02  7.000   9.00    3' /
DATA argon(96) / '-0.226471343048d-02 24.000   5.00    4' /
DATA argon(97) / ' 0.176164561964d-02 22.000   6.00    4' /
DATA argon(98) / ' 0.585524544828d-02  3.000   2.00    2 2 -20.0  -250.0  1.11  1.  0.  0.  0.' /
DATA argon(99) / '-0.692519082700d+00  1.000   1.00    2 2 -20.0  -375.0  1.14  1.  0.  0.  0.' /
DATA argon(100) / ' 0.153154900305d+01  0.000   2.00    2 2 -20.0  -300.0  1.17  1.  0.  0.  0.' /
DATA argon(101) / '-0.273804474498d-02  0.000   3.00    2 2 -20.0  -225.0  1.11  1.  0.  0.  0.' /
DATA argon(102) / '' /
DATA argon(103) / '' /
DATA argon(104) / '#AUX               !auxiliary model specification' /
DATA argon(105) / 'CPP  ideal gas heat capacity function' /
DATA argon(106) / '?LITERATURE REFERENCE \' /
DATA argon(107) / '?Tegeler, Ch., Span, R., and Wagner, W.,' /
DATA argon(108) / '? "A New Equation of State for Argon Covering the Fluid Region for' /
DATA argon(109) / '? Temperatures from the Melting Line to 700 K at Pressures up to 1000 MPa,"' /
DATA argon(110) / '? J. Phys. Chem. Ref. Data, 28(3):779-850, 1999.' /
DATA argon(111) / '?\' /
DATA argon(112) / '!end of info section' /
DATA argon(113) / '83.8058            !lower temperature limit [K]' /
DATA argon(114) / '700.0              !upper temperature limit [K]' /
DATA argon(115) / '0.0                !upper pressure limit [kPa]' /
DATA argon(116) / '0.0                !maximum density [mol/L]' /
DATA argon(117) / '1.0          8.31451                   !reducing parameters for T, Cp0' /
DATA argon(118) / '  1  0    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA argon(119) / ' 0.25d1            0.00' /
DATA argon(120) / '' /
DATA argon(121) / '' /
DATA argon(122) / '@EOS               !equation of state specification' /
DATA argon(123) / 'FEK  Helmholtz equation of state for argon of Kunz and Wagner (2004).' /
DATA argon(124) / '?LITERATURE REFERENCE \' /
DATA argon(125) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA argon(126) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA argon(127) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA argon(128) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA argon(129) / '?\' /
DATA argon(130) / '!end of info section' /
DATA argon(131) / '83.8058            !lower temperature limit [K]' /
DATA argon(132) / '700.0              !upper temperature limit [K]' /
DATA argon(133) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(134) / '50.65              !maximum density [mol/L]' /
DATA argon(135) / 'PHK                                    !pointer to Cp0 model' /
DATA argon(136) / '39.948                                 !molecular weight [g/mol]' /
DATA argon(137) / '83.8058                                !triple point temperature [K]' /
DATA argon(138) / '69.03                                  !pressure at triple point [kPa]' /
DATA argon(139) / '35.5                                   !density at triple point [mol/L]' /
DATA argon(140) / '87.29                                  !normal boiling point temperature [K]' /
DATA argon(141) / '-0.0006                                !acentric factor' /
DATA argon(142) / '150.687      4879.8      13.407429659  !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA argon(143) / '150.687                  13.407429659  !reducing parameters [K, mol/L]' /
DATA argon(144) / '8.314472                               !gas constant [J/mol-K]' /
DATA argon(145) / '  12  4      0  0      0  0            !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA argon(146) / ' 0.85095714803969       0.250  1.  0' /
DATA argon(147) / '-0.24003222943480d1     1.125  1.  0' /
DATA argon(148) / ' 0.54127841476466       1.500  1.  0' /
DATA argon(149) / ' 0.16919770692538d-1    1.375  2.  0' /
DATA argon(150) / ' 0.68825965019035d-1    0.250  3.  0' /
DATA argon(151) / ' 0.21428032815338d-3    0.875  7.  0' /
DATA argon(152) / ' 0.17429895321992       0.625  2.  1' /
DATA argon(153) / '-0.33654495604194d-1    1.750  5.  1' /
DATA argon(154) / '-0.13526799857691       3.625  1.  2' /
DATA argon(155) / '-0.16387350791552d-1    3.625  4.  2' /
DATA argon(156) / '-0.24987666851475d-1    14.5   3.  3' /
DATA argon(157) / ' 0.88769204815709d-2    12.0   4.  3' /
DATA argon(158) / '' /
DATA argon(159) / '' /
DATA argon(160) / '#AUX               !auxiliary model specification' /
DATA argon(161) / 'PHK  Helmholtz form for the ideal-gas state for argon of Kunz and Wagner (2004).' /
DATA argon(162) / '?LITERATURE REFERENCE \' /
DATA argon(163) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA argon(164) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA argon(165) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA argon(166) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA argon(167) / '?\' /
DATA argon(168) / '!end of info section' /
DATA argon(169) / '0.                 !lower temperature limit [K]' /
DATA argon(170) / '1000.0             !upper temperature limit [K]' /
DATA argon(171) / '0.0                !upper pressure limit [kPa]' /
DATA argon(172) / '0.0                !maximum density [mol/L]' /
DATA argon(173) / '1 2  0  0 0  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA argon(174) / '    1.5          1.             !ai, ti for [ai*log(tau**ti)] terms' /
DATA argon(175) / '    8.3166315    0.             !aj, ti for [ai*tau**ti] terms' /
DATA argon(176) / '   -4.9465026    1.' /
DATA argon(177) / '' /
DATA argon(178) / '' /
DATA argon(179) / '@EOS               !equation of state specification' /
DATA argon(180) / 'FE1  Helmholtz equation of state for argon of Stewart and Jacobsen (1989).' /
DATA argon(181) / '?LITERATURE REFERENCE \' /
DATA argon(182) / '?Stewart, R.B. and Jacobsen, R.T,' /
DATA argon(183) / '? "Thermodynamic Properties of Argon from the Triple Point to 1200 K at' /
DATA argon(184) / '? Pressures to 1000 MPa,"' /
DATA argon(185) / '? J. Phys. Chem. Ref. Data, 18(2):639-798, 1989.' /
DATA argon(186) / '?\' /
DATA argon(187) / '!end of info section' /
DATA argon(188) / '83.804             !lower temperature limit [K]' /
DATA argon(189) / '1200.0             !upper temperature limit [K]' /
DATA argon(190) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(191) / '45.814             !maximum density [mol/L]' /
DATA argon(192) / 'CP1                                    !pointer to Cp0 model' /
DATA argon(193) / '39.948                                 !molecular weight [g/mol]' /
DATA argon(194) / '83.804                                 !triple point temperature [K]' /
DATA argon(195) / '68.961                                 !pressure at triple point [kPa]' /
DATA argon(196) / '35.475                                 !density at triple point [mol/L]' /
DATA argon(197) / '87.293                                 !normal boiling point temperature [K]' /
DATA argon(198) / '-0.004                                 !acentric factor' /
DATA argon(199) / '150.6633     4860.0       13.29        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA argon(200) / '150.6633                  13.29        !reducing parameters [K, mol/L]' /
DATA argon(201) / '8.31434                                !gas constant [J/mol-K]' /
DATA argon(202) / '      28  4      0  0       0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA argon(203) / ' 0.791867571500d+00  0.250   1.00    0 !a(i),t(i),d(i),l(i)' /
DATA argon(204) / '-0.163334615100d+01  1.000   1.00    0' /
DATA argon(205) / '-0.439530293000d+00  3.000   1.00    0' /
DATA argon(206) / ' 0.103389999900d+00  4.000   1.00    0' /
DATA argon(207) / ' 0.206180166400d+00  0.250   2.00    0' /
DATA argon(208) / '-0.288868177600d+00  1.000   2.00    0' /
DATA argon(209) / ' 0.439801055000d+00  2.500   2.00    0' /
DATA argon(210) / '-0.842955039100d-01  3.500   2.00    0' /
DATA argon(211) / '-0.215565865400d+00  0.750   3.00    0' /
DATA argon(212) / ' 0.478650909900d+00  1.000   3.00    0' /
DATA argon(213) / '-0.352588459300d+00  1.500   3.00    0' /
DATA argon(214) / ' 0.301507369200d-01  2.500   3.00    0' /
DATA argon(215) / ' 0.298767905900d-01  1.000   4.00    0' /
DATA argon(216) / '-0.152256858300d-01  2.000   4.00    0' /
DATA argon(217) / ' 0.743578578600d-03  2.000   6.00    0' /
DATA argon(218) / ' 0.709954162400d-01  5.000   1.00    3' /
DATA argon(219) / '-0.290423718500d-01  7.000   1.00    3' /
DATA argon(220) / '-0.622307852500d-01  5.000   2.00    2' /
DATA argon(221) / ' 0.141089518700d-03 22.000   2.00    4' /
DATA argon(222) / '-0.148124178300d-02 16.000   2.00    6' /
DATA argon(223) / ' 0.302334278400d-01 10.000   3.00    3' /
DATA argon(224) / '-0.612678468500d-01 14.000   3.00    3' /
DATA argon(225) / ' 0.270996709000d-01 16.000   3.00    3' /
DATA argon(226) / ' 0.941103440500d-01  4.000   4.00    2' /
DATA argon(227) / '-0.729164511400d-02  8.000   4.00    2' /
DATA argon(228) / '-0.158631497600d-02 10.000   4.00    4' /
DATA argon(229) / ' 0.951094881300d-03  5.000   8.00    2' /
DATA argon(230) / ' 0.778618184400d-03  6.000   8.00    2' /
DATA argon(231) / '' /
DATA argon(232) / '' /
DATA argon(233) / '#AUX               !auxiliary model specification' /
DATA argon(234) / 'CP1  ideal gas heat capacity function' /
DATA argon(235) / '?LITERATURE REFERENCE \' /
DATA argon(236) / '?Stewart, R.B. and Jacobsen, R.T,' /
DATA argon(237) / '? "Thermodynamic Properties of Argon from the Triple Point to 1200 K at' /
DATA argon(238) / '? Pressures to 1000 MPa,"' /
DATA argon(239) / '? J. Phys. Chem. Ref. Data, 18(2):639-798, 1989.' /
DATA argon(240) / '?\' /
DATA argon(241) / '!end of info section' /
DATA argon(242) / '83.804             !lower temperature limit [K]' /
DATA argon(243) / '1200.0             !upper temperature limit [K]' /
DATA argon(244) / '0.0                !upper pressure limit [kPa]' /
DATA argon(245) / '0.0                !maximum density [mol/L]' /
DATA argon(246) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA argon(247) / '  1  0    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA argon(248) / ' 2.5               0.00' /
DATA argon(249) / '' /
DATA argon(250) / '' /
DATA argon(251) / '@EOS               !equation of state specification' /
DATA argon(252) / 'FES  short Helmholtz equation of state for argon of Span and Wagner (2003).' /
DATA argon(253) / '?LITERATURE REFERENCE \' /
DATA argon(254) / '?Span, R. and Wagner, W.' /
DATA argon(255) / '? "Equations of State for Technical Applications. II. Results for Nonpolar Fluids,"' /
DATA argon(256) / '? Int. J. Thermophys., 24(1):41-109, 2003.' /
DATA argon(257) / '?\' /
DATA argon(258) / '?The uncertainties of the equation of state are approximately 0.2% (to' /
DATA argon(259) / '?0.5% at high pressures) in density, 1% (in the vapor phase) to 2% in' /
DATA argon(260) / '?heat capacity, 1% (in the vapor phase) to 2% in the speed of sound, and' /
DATA argon(261) / '?0.2% in vapor pressure, except in the critical region.' /
DATA argon(262) / '?\' /
DATA argon(263) / '!end of info section' /
DATA argon(264) / '83.8058            !lower temperature limit [K]' /
DATA argon(265) / '600.0              !upper temperature limit [K]' /
DATA argon(266) / '100000.0           !upper pressure limit [kPa]' /
DATA argon(267) / '50.65              !maximum density [mol/L]' /
DATA argon(268) / 'CPP                                    !pointer to Cp0 model' /
DATA argon(269) / '39.948                                 !molecular weight [g/mol]' /
DATA argon(270) / '83.8058                                !triple point temperature [K]' /
DATA argon(271) / '69.026                                 !pressure at triple point [kPa]' /
DATA argon(272) / '35.498                                 !density at triple point [mol/L]' /
DATA argon(273) / '87.289                                 !normal boiling point temperature [K]' /
DATA argon(274) / '-0.002                                 !acentric factor' /
DATA argon(275) / '150.687      4863.0       13.40743     !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA argon(276) / '150.687                   13.40743     !reducing parameters [K, mol/L]' /
DATA argon(277) / '8.31451                                !gas constant [J/mol-K]' /
DATA argon(278) / '      12  4      0  0       0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA argon(279) / ' 0.850957150000E+00  0.25    1.0     0 !a(i),t(i),d(i),l(i)' /
DATA argon(280) / '-0.240032230000E+01  1.125   1.0     0' /
DATA argon(281) / ' 0.541278410000E+00  1.5     1.0     0' /
DATA argon(282) / ' 0.169197710000E-01  1.375   2.0     0' /
DATA argon(283) / ' 0.688259650000E-01  0.25    3.0     0' /
DATA argon(284) / ' 0.214280330000E-03  0.875   7.0     0' /
DATA argon(285) / ' 0.174298950000E+00  0.625   2.0     1' /
DATA argon(286) / '-0.336544960000E-01  1.75    5.0     1' /
DATA argon(287) / '-0.135268000000E+00  3.625   1.0     2' /
DATA argon(288) / '-0.163873510000E-01  3.625   4.0     2' /
DATA argon(289) / '-0.249876670000E-01 14.5     3.0     3' /
DATA argon(290) / ' 0.887692050000E-02 12.0     4.0     3' /
DATA argon(291) / '' /
DATA argon(292) / '' /
DATA argon(293) / '@EOS               !equation of state specification' /
DATA argon(294) / 'BWR  MBWR equation of state for argon of Younglove (1982).' /
DATA argon(295) / '?LITERATURE REFERENCE \' /
DATA argon(296) / '?Younglove, B.A.,' /
DATA argon(297) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA argon(298) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA argon(299) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA argon(300) / '?\' /
DATA argon(301) / '!end of info section' /
DATA argon(302) / '83.80              !lower temperature limit [K]' /
DATA argon(303) / '400.0              !upper temperature limit [K]' /
DATA argon(304) / '101000.0           !upper pressure limit [kPa]' /
DATA argon(305) / '50.65              !maximum density [mol/L]' /
DATA argon(306) / 'CP2                                    !pointer to Cp0 model' /
DATA argon(307) / '39.948                                 !molecular weight [g/mol]' /
DATA argon(308) / '83.80                                  !triple point temperature [K]' /
DATA argon(309) / '68.906                                 !pressure at triple point [kPa]' /
DATA argon(310) / '35.4                                   !density at triple point [mol/L]' /
DATA argon(311) / '87.302                                 !normal boiling point temperature [K]' /
DATA argon(312) / '-0.002                                 !acentric factor' /
DATA argon(313) / '150.86       4905.8       13.41        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA argon(314) / '150.86                    13.41        !reducing parameters [K, mol/L]' /
DATA argon(315) / '13.418                                 !gamma' /
DATA argon(316) / '0.0831434                              !gas constant [L-bar/mol-K]' /
DATA argon(317) / '      32       1                       !Nterm, Ncoeff per term' /
DATA argon(318) / ' -0.6569731294d-03     0.1822957801d-00     -0.3649470141d+01' /
DATA argon(319) / '  0.1232012107d+03    -0.8613578274d+04      0.7978579691d-04' /
DATA argon(320) / ' -0.2911489110d-01     0.7581821758d+01      0.8780488169d+04' /
DATA argon(321) / '  0.1423145989d-06     0.1674146131d-02     -0.3200447909d-00' /
DATA argon(322) / '  0.2561766372d-04    -0.5475934941d-03     -0.4505032058d-00' /
DATA argon(323) / '  0.2013254653d-04    -0.1678941273d-06      0.4207329271d-03' /
DATA argon(324) / ' -0.5444212996d-05    -0.8004855011d+04     -0.1319304201d+06' /
DATA argon(325) / ' -0.4954923930d+02     0.8092132177d+05     -0.9870104061d-01' /
DATA argon(326) / '  0.2020441562d+01    -0.1637417205d-03     -0.7038944136d-00' /
DATA argon(327) / ' -0.1154324539d-06     0.1555990117d-04     -0.1492178536d-09' /
DATA argon(328) / ' -0.1001356071d-07     0.2933963216d-06' /
DATA argon(329) / '' /
DATA argon(330) / '' /
DATA argon(331) / '#AUX               !auxiliary model specification' /
DATA argon(332) / 'CP2  ideal gas heat capacity function of Younglove' /
DATA argon(333) / '?LITERATURE REFERENCE \' /
DATA argon(334) / '?Younglove, B.A.,' /
DATA argon(335) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA argon(336) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA argon(337) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA argon(338) / '?\' /
DATA argon(339) / '!end of info section' /
DATA argon(340) / '83.80              !lower temperature limit [K]' /
DATA argon(341) / '600.0              !upper temperature limit [K]' /
DATA argon(342) / '0.0                !upper pressure limit [kPa]' /
DATA argon(343) / '0.0                !maximum density [mol/L]' /
DATA argon(344) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA argon(345) / '  1  0    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA argon(346) / '  2.5000000000d+0       0.00d0' /
DATA argon(347) / '' /
DATA argon(348) / '' /
DATA argon(349) / '@EOS' /
DATA argon(350) / 'PRT  translated Peng-Robinson equation' /
DATA argon(351) / '?LITERATURE REFERENCES \' /
DATA argon(352) / '?  volume translation of Peng Robinson EOS' /
DATA argon(353) / '?  translation computed so that sat. liquid density at Tr=0.7 matches FEQ  Helmholtz equation' /
DATA argon(354) / '?  of state for Ar of Tegeler et al. (1999).' /
DATA argon(355) / '!end of info section' /
DATA argon(356) / '83.8058            !lower temperature limit [K]' /
DATA argon(357) / '2000.0             !upper temperature limit [K]' /
DATA argon(358) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(359) / '50.65              !maximum density [mol/L]' /
DATA argon(360) / 'CPP                !pointer to Cp0 model' /
DATA argon(361) / '39.948             !molecular weight [g/mol]' /
DATA argon(362) / '-0.00219           !acentric factor' /
DATA argon(363) / '150.687            !critical temperature [K]' /
DATA argon(364) / '4863.0             !critical pressure [kPa]' /
DATA argon(365) / '13.4074            !critical density [mol/L]' /
DATA argon(366) / '8.314472           !gas constant [J/mol-K]' /
DATA argon(367) / '1                  !Number of parameters' /
DATA argon(368) / '-0.0034' /
DATA argon(369) / '' /
DATA argon(370) / '' /
DATA argon(371) / '#TCX               !thermal conductivity model specification' /
DATA argon(372) / 'TC1  pure fluid thermal conductivity model of Lemmon and Jacobsen (2004).' /
DATA argon(373) / '?LITERATURE REFERENCE \' /
DATA argon(374) / '?Lemmon, E.W. and Jacobsen, R.T,' /
DATA argon(375) / '? "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon,' /
DATA argon(376) / '? and Air,"' /
DATA argon(377) / '? Int. J. Thermophys., 25:21-69, 2004.' /
DATA argon(378) / '?\' /
DATA argon(379) / '?The uncertainty for the dilute gas is 2% with increasing uncertainties' /
DATA argon(380) / '?near the triple point.  For the non-dilute gas, the uncertainty is 2%' /
DATA argon(381) / '?for temperatures greater than 170 K. The uncertainty is 3% at' /
DATA argon(382) / '?temperatures less than the critical point and 5% in the critical region,' /
DATA argon(383) / '?except for states very near the critical point.' /
DATA argon(384) / '?\' /
DATA argon(385) / '!end of info section' /
DATA argon(386) / '83.8058            !lower temperature limit [K]' /
DATA argon(387) / '2000.0             !upper temperature limit [K]' /
DATA argon(388) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(389) / '50.65              !maximum density [mol/L]' /
DATA argon(390) / '2   0              !# terms for dilute gas function:  numerator, denominator' /
DATA argon(391) / '150.687   1.0d-3   !reducing parameters for T, tcx' /
DATA argon(392) / ' 0.8158   -97.0    !coeff, power in T' /
DATA argon(393) / '-0.4320     0.77' /
DATA argon(394) / '7   0              !# terms for background gas function:  numerator, denominator' /
DATA argon(395) / '150.687   13.40742965 1.0d-3    !reducing parameters for T, rho, tcx' /
DATA argon(396) / '13.73           0.0  1.0  0.0 !coeff, powers of T, rho, exp(rho)' /
DATA argon(397) / ' 10.07        0.0    2.   0.' /
DATA argon(398) / ' 0.7375       0.0    4.   0.' /
DATA argon(399) / '-33.96       -0.8    5.   2.' /
DATA argon(400) / ' 20.47       -1.2    6.   2.' /
DATA argon(401) / '-2.274       -0.8    9.   2.' /
DATA argon(402) / '-3.973       -0.5    1.   4.' /
DATA argon(403) / 'TK3                !pointer to critical enhancement auxiliary function' /
DATA argon(404) / '' /
DATA argon(405) / '' /
DATA argon(406) / '#AUX               !thermal conductivity critical enhancement model' /
DATA argon(407) / 'TK3  thermal conductivity critical enhancement of Lemmon and Jacobsen (2004).' /
DATA argon(408) / '?LITERATURE REFERENCE \' /
DATA argon(409) / '?\' /
DATA argon(410) / '!end of info section' /
DATA argon(411) / '83.8058            !lower temperature limit [K]' /
DATA argon(412) / '2000.0             !upper temperature limit [K]' /
DATA argon(413) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(414) / '50.65              !maximum density [mol/L]' /
DATA argon(415) / '9  0  0  0         !# terms:  terms, spare, spare, spare' /
DATA argon(416) / '1.0    1.0  1.0    !reducing par for T, rho, tcx (mW/m-K)' /
DATA argon(417) / '0.630d0            !gnu (universal exponent)' /
DATA argon(418) / '1.2415d0           !gamma (universal exponent)' /
DATA argon(419) / '1.01d0             !R0 (universal amplitude)' /
DATA argon(420) / ' 0.065d0           !z (universal exponent--not used for t.c., only viscosity)' /
DATA argon(421) / ' 1.00d0            !c (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA argon(422) / ' 0.13E-09          !xi0 (amplitude) [m]' /
DATA argon(423) / ' 0.55E-01          !gam0 (amplitude) [-]' /
DATA argon(424) / ' 0.32E-09          !qd_inverse (modified effective cutoff parameter) [m]' /
DATA argon(425) / '301.374            !tref (reference temperature) [K]' /
DATA argon(426) / '' /
DATA argon(427) / '' /
DATA argon(428) / '#ETA               !viscosity model specification' /
DATA argon(429) / 'VS1  pure fluid viscosity model of Lemmon and Jacobsen (2004).' /
DATA argon(430) / '?LITERATURE REFERENCE \' /
DATA argon(431) / '?Lemmon, E.W. and Jacobsen, R.T,' /
DATA argon(432) / '? "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon,' /
DATA argon(433) / '? and Air,"' /
DATA argon(434) / '? Int. J. Thermophys., 25:21-69, 2004.' /
DATA argon(435) / '?\' /
DATA argon(436) / '?The uncertainty is 0.5% in the dilute gas.  Away from the dilute gas' /
DATA argon(437) / '?(pressures greater than 1 MPa and in the liquid), the uncertainties are' /
DATA argon(438) / '?as low as 1% between 270 and 300 K at pressures less than 100 MPa, and' /
DATA argon(439) / '?increase outside that range.  The uncertainties are around 2% at' /
DATA argon(440) / '?temperatures of 180 K and higher.  Below this and away from the critical' /
DATA argon(441) / '?region, the uncertainties steadily increase to around 5% at the triple' /
DATA argon(442) / '?points of the fluids.  The uncertainties in the critical region are' /
DATA argon(443) / '?higher.' /
DATA argon(444) / '?\' /
DATA argon(445) / '!end of info section' /
DATA argon(446) / '83.8058            !lower temperature limit [K]' /
DATA argon(447) / '2000.0             !upper temperature limit [K]' /
DATA argon(448) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(449) / '50.65              !maximum density [mol/L]' /
DATA argon(450) / '1                  !number of terms associated with dilute-gas function' /
DATA argon(451) / 'CI1                !pointer to reduced effective collision cross-section model' /
DATA argon(452) / '0.335              !Lennard-Jones coefficient sigma [nm]' /
DATA argon(453) / '143.2              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA argon(454) / '1.0    1.0         !reducing parameters for T, eta' /
DATA argon(455) / '0.168729283  0.5   !Chapman-Enskog term' /
DATA argon(456) / '0                  !number of terms for initial density dependence' /
DATA argon(457) / '0 6 0 0 0 0        !# resid terms:  close-packed density;  simple poly; numerator of rational poly; denominator of rat. poly; numerator of exponential; denominator of exponential' /
DATA argon(458) / '150.687   13.40742965 1.0           !reducing parameters for T, rho, eta' /
DATA argon(459) / ' 12.19       -0.42   1.   0.   0    !simple polynomial terms' /
DATA argon(460) / ' 13.99        0.0    2.   0.   0' /
DATA argon(461) / ' 0.005027    -0.95  10.   0.   0' /
DATA argon(462) / '-18.93       -0.5    5.   0.   2' /
DATA argon(463) / '-6.698       -0.9    1.   0.   4' /
DATA argon(464) / '-3.827       -0.8    2.   0.   4' /
DATA argon(465) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA argon(466) / '' /
DATA argon(467) / '' /
DATA argon(468) / '#AUX               !collision integral specification' /
DATA argon(469) / 'CI1  collision integral model of Lemmon and Jacobsen (2004).' /
DATA argon(470) / '?LITERATURE REFERENCE \' /
DATA argon(471) / '?\' /
DATA argon(472) / '!end of info section' /
DATA argon(473) / '1.0                !lower temperature limit [K]' /
DATA argon(474) / '10000.0            !upper temperature limit [K]' /
DATA argon(475) / '0.0                !(dummy) upper pressure limit' /
DATA argon(476) / '0.0                !(dummy) maximum density' /
DATA argon(477) / '5                  !number of terms' /
DATA argon(478) / '  0.431      0     !coeff, power of Tstar' /
DATA argon(479) / ' -0.4623     1' /
DATA argon(480) / '  0.08406    2' /
DATA argon(481) / '  0.005341   3' /
DATA argon(482) / ' -0.00331    4' /
DATA argon(483) / '' /
DATA argon(484) / '' /
DATA argon(485) / '@ETA               !viscosity model specification' /
DATA argon(486) / 'VS3  pure fluid viscosity model of Younglove and Hanley (1986).' /
DATA argon(487) / '?REFERENCE \' /
DATA argon(488) / '?Younglove, B.A. and Hanley, H.J.M.,' /
DATA argon(489) / '? "The Viscosity and Thermal Conductivity Coefficients of Gaseous and' /
DATA argon(490) / '? Liquid Argon,"' /
DATA argon(491) / '? J. Phys. Chem. Ref. Data, 15(4):1323-1337, 1986.' /
DATA argon(492) / '?\' /
DATA argon(493) / '?The uncertainty in viscosity is 2% below 100 MPa and 3% for higher pressures.' /
DATA argon(494) / '?\' /
DATA argon(495) / '!end of info section' /
DATA argon(496) / '55.0               !lower temperature limit [K]' /
DATA argon(497) / '2500.0             !upper temperature limit [K]' /
DATA argon(498) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(499) / '50.65              !maximum density [mol/L]' /
DATA argon(500) / '9   0              !# terms for dilute gas function:  numerator, denominator' /
DATA argon(501) / '1.0     1.0        !reducing parameters for T, eta' /
DATA argon(502) / '-0.8973188257d5   -1.00d0   !coeff, power in T' /
DATA argon(503) / ' 0.8259113473d5   -0.666666666667d0' /
DATA argon(504) / '-0.2766475915d5   -0.333333333333d0' /
DATA argon(505) / ' 0.3068539784d4    0.00d0' /
DATA argon(506) / ' 0.4553103615d3    0.333333333333d0' /
DATA argon(507) / '-0.1793443839d3    0.666666666667d0' /
DATA argon(508) / ' 0.2272225106d2    1.00d0' /
DATA argon(509) / '-0.1350672796d1    1.333333333333d0' /
DATA argon(510) / ' 0.3183693230d-1   1.666666666667d0' /
DATA argon(511) / '11  3              !# terms for background gas function:  numerator, denominator' /
DATA argon(512) / '1.0    1.0    1.0                         !reducing par for T, rho (rho_c), eta' /
DATA argon(513) / ' 0.5927733783   0.0  1.0  0.0 !coeff, powers of T, rho, spare for future use' /
DATA argon(514) / '-0.4251221169d2   -1.00d0   1.00d0   0.' /
DATA argon(515) / '-0.2698477165d-1   0.00d0   2.00d0   0.' /
DATA argon(516) / ' 0.3727762288d2   -1.00d0   2.00d0   0.' /
DATA argon(517) / '-0.3958508720d4   -2.00d0   2.00d0   0.' /
DATA argon(518) / ' 0.3636730841d-2   0.00d0   3.00d0   0.' /
DATA argon(519) / '-0.2633471347d1   -1.00d0   3.00d0   0.' /
DATA argon(520) / ' 0.2936563322d3   -2.00d0   3.00d0   0.' /
DATA argon(521) / '-0.3811869019d-4   0.00d0   4.00d0   0.' /
DATA argon(522) / ' 0.4451947464d-1  -1.00d0   4.00d0   0.' /
DATA argon(523) / '-0.5385874487d1   -2.00d0   4.00d0   0.' /
DATA argon(524) / ' 1.0000000000      0.00d0   0.00d0   0.' /
DATA argon(525) / '-0.1115054926d-1   0.00d0   1.00d0   0.' /
DATA argon(526) / '-0.1328893444d1   -1.00d0   1.00d0   0.' /
DATA argon(527) / 'NUL                !pointer to critical enhancement auxiliary function' /
DATA argon(528) / '' /
DATA argon(529) / '' /
DATA argon(530) / '@TCX               !thermal conductivity model specification' /
DATA argon(531) / 'TC3  pure fluid thermal conductivity model of Younglove (1982).' /
DATA argon(532) / '?LITERATURE REFERENCE \' /
DATA argon(533) / '?Younglove, B.A.,' /
DATA argon(534) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA argon(535) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA argon(536) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA argon(537) / '?\' /
DATA argon(538) / '!end of info section' /
DATA argon(539) / '83.8058            !lower temperature limit [K]' /
DATA argon(540) / '600.0              !upper temperature limit [K]' /
DATA argon(541) / '100000.0           !upper pressure limit [kPa]' /
DATA argon(542) / '50.65              !maximum density [mol/L]' /
DATA argon(543) / '0.3297             !Lennard-Jones coefficient sigma [nm]' /
DATA argon(544) / '152.8              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA argon(545) / '0.16871158559818   !const in Eq 20 = 5/16*(k*MW/1000/pi/Na)**0.5*1.0d12' /
DATA argon(546) / ' 0                 !exponent in Eq 20 for T' /
DATA argon(547) / ' 2.64712871543D-02 !eta0 (eq 20): coeffs of {(e/kT)**((4-n)/3)}' /
DATA argon(548) / '-.216629583011974' /
DATA argon(549) / '0.709700888884514' /
DATA argon(550) / '-1.21908891344223' /
DATA argon(551) / ' 1.20168985706305' /
DATA argon(552) / '-.700084760049098' /
DATA argon(553) / '0.24816605762696' /
DATA argon(554) / '-4.79479287295D-02' /
DATA argon(555) / ' 3.93679190444D-03' /
DATA argon(556) / ' 9.64428741429D-04 !coeff for initial density dependence of viscosity (eq 21); Fv(1)' /
DATA argon(557) / ' 3.02391316601D-04 !Fv(2)' /
DATA argon(558) / ' 1                 !Fv(3)' /
DATA argon(559) / ' 152.8             !Fv(4)' /
DATA argon(560) / '-33.327027332      !coefficients for residual viscosity, eqs (22 - 25)' /
DATA argon(561) / '-355.59415848      !Ev(2)' /
DATA argon(562) / ' 22.2441164817987  !Ev(3)' /
DATA argon(563) / ' 1663.62775376509  !Ev(4)' /
DATA argon(564) / ' 0                 !Ev(5)' /
DATA argon(565) / ' 0                 !Ev(6)' /
DATA argon(566) / ' 0                 !Ev(7)' /
DATA argon(567) / ' 25.0325423049965  !Ev(8)' /
DATA argon(568) / ' 1.7124            !F' /
DATA argon(569) / '0.00000003669      !rm' /
DATA argon(570) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA argon(571) / '' /
DATA argon(572) / '' /
DATA argon(573) / '@ETA               !viscosity model specification' /
DATA argon(574) / 'VS2  pure fluid viscosity model of Younglove (1982).' /
DATA argon(575) / '?LITERATURE REFERENCE \' /
DATA argon(576) / '?Younglove, B.A.,' /
DATA argon(577) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA argon(578) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA argon(579) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA argon(580) / '?\' /
DATA argon(581) / '!end of info section' /
DATA argon(582) / '83.8058            !lower temperature limit [K]' /
DATA argon(583) / '600.0              !upper temperature limit [K]' /
DATA argon(584) / '100000.0           !upper pressure limit [kPa]' /
DATA argon(585) / '50.65              !maximum density [mol/L]' /
DATA argon(586) / 'CI2                !pointer to collision integral model' /
DATA argon(587) / '0.3297             !Lennard-Jones coefficient sigma [nm]' /
DATA argon(588) / '152.8              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA argon(589) / '0.16871158559818   !const in Eq 19 = 5/16*(k*MW/1000/pi/Na)**0.5*1.0d12' /
DATA argon(590) / ' 0                 !exponent in Eq 20 for T' /
DATA argon(591) / ' 5.85384107393D-03 !coeff for initial density dependence of viscosity (eq 21); Fv(1)' /
DATA argon(592) / '-3.09546765250D-03 !Fv(2)' /
DATA argon(593) / ' 1.4               !Fv(3)' /
DATA argon(594) / ' 152.8             !Fv(4)' /
DATA argon(595) / '-12.313579086      !coefficients for residual viscosity, eqs (22 - 25)' /
DATA argon(596) / ' 40.136071933      !Ev(2)' /
DATA argon(597) / ' 11.6160872385243  !Ev(3)' /
DATA argon(598) / '-413.04094973717   !Ev(4)' /
DATA argon(599) / ' 4.13624595833D-02 !Ev(5)' /
DATA argon(600) / ' 7.96883967907912  !Ev(6)' /
DATA argon(601) / ' 234.196850483958  !Ev(7)' /
DATA argon(602) / ' 13.4424752177831  !Ev(8)' /
DATA argon(603) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA argon(604) / '' /
DATA argon(605) / '' /
DATA argon(606) / '@AUX               !collision integral specification' /
DATA argon(607) / 'CI2  collision integral model of Younglove (1982).' /
DATA argon(608) / '?LITERATURE REFERENCE \' /
DATA argon(609) / '?Younglove, B.A.,' /
DATA argon(610) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA argon(611) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA argon(612) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA argon(613) / '?\' /
DATA argon(614) / '!end of info section' /
DATA argon(615) / '83.8058            !lower temperature limit [K]' /
DATA argon(616) / '625.0              !upper temperature limit [K]' /
DATA argon(617) / '0.0                !(dummy) upper pressure limit' /
DATA argon(618) / '0.0                !(dummy) maximum density' /
DATA argon(619) / '9                  !number of terms' /
DATA argon(620) / ' 25.7830291943396  !eta0 (eq 20): coeffs of {(e/kT)**((4-n)/3)}' /
DATA argon(621) / '-234.320222858983' /
DATA argon(622) / ' 814.636688705024' /
DATA argon(623) / '-1452.04353466585' /
DATA argon(624) / ' 1467.17535558104' /
DATA argon(625) / '-870.164951237067' /
DATA argon(626) / ' 313.024934147423' /
DATA argon(627) / '-61.2072628957372' /
DATA argon(628) / ' 5.07700488990665' /
DATA argon(629) / '' /
DATA argon(630) / '' /
DATA argon(631) / '#STN        !surface tension specification' /
DATA argon(632) / 'ST1  surface tension model of Mulero et al. (2012)' /
DATA argon(633) / '?LITERATURE REFERENCE \' /
DATA argon(634) / '?Mulero, A., Cachadina, I., and Parra, M.I.' /
DATA argon(635) / '? "Recommended Correlations for the Surface Tension of Common Fluids,"' /
DATA argon(636) / '? J. Phys. Chem. Ref. Data, 41, 043105, 2012.' /
DATA argon(637) / '?\' /
DATA argon(638) / '!end of info section' /
DATA argon(639) / '0.0                !lower temperature limit [K]' /
DATA argon(640) / '150.687            !upper temperature limit [K]' /
DATA argon(641) / '0.0                !(dummy) upper pressure limit' /
DATA argon(642) / '0.0                !(dummy) maximum density' /
DATA argon(643) / '1                           !number of terms in surface tension model' /
DATA argon(644) / '150.687                     !critical temperature used in fit (dummy)' /
DATA argon(645) / ' 0.037       1.25           !sigma0 and n' /
DATA argon(646) / '' /
DATA argon(647) / '' /
DATA argon(648) / '#DE         !dielectric constant specification' /
DATA argon(649) / 'DE3  dielectric constant model of Harvey and Lemmon (2005).' /
DATA argon(650) / '?LITERATURE REFERENCE \' /
DATA argon(651) / '?Harvey, A.H. and Lemmon, E.W.' /
DATA argon(652) / '? "Method for Estimating the Dielectric Constant of Natural Gas Mixtures,"' /
DATA argon(653) / '? Int. J. Thermophys., 26(1):31-46, 2005.' /
DATA argon(654) / '?\' /
DATA argon(655) / '!end of info section' /
DATA argon(656) / '0.0                !lower temperature limit [K]' /
DATA argon(657) / '2000.0             !upper temperature limit [K]' /
DATA argon(658) / '0.0                !(dummy) upper pressure limit' /
DATA argon(659) / '0.0                !(dummy) maximum density' /
DATA argon(660) / '273.16 1000.0 1.0  !reducing parameters for t and d' /
DATA argon(661) / '0 1 3 0 0 0                         !number of terms in dielectric constant model' /
DATA argon(662) / ' 4.1414           0.    1.    0.    !coef, t exp, d exp' /
DATA argon(663) / ' 1.597            0.    2.    0.' /
DATA argon(664) / ' 0.262            1.    2.    0.' /
DATA argon(665) / '-117.9            0.    3.1   0.' /
DATA argon(666) / '' /
DATA argon(667) / '' /
DATA argon(668) / '#MLT        !melting line specification' /
DATA argon(669) / 'ML1  melting line model of Tegeler et al. (1999).' /
DATA argon(670) / '?LITERATURE REFERENCE \' /
DATA argon(671) / '?Tegeler, Ch., Span, R., and Wagner, W.,' /
DATA argon(672) / '? "A New Equation of State for Argon Covering the Fluid Region for' /
DATA argon(673) / '? Temperatures from the Melting Line to 700 K at Pressures up to 1000 MPa,"' /
DATA argon(674) / '? J. Phys. Chem. Ref. Data, 28(3):779-850, 1999.' /
DATA argon(675) / '?\' /
DATA argon(676) / '!end of info section' /
DATA argon(677) / '83.8058            !lower temperature limit [K]' /
DATA argon(678) / '700.0              !upper temperature limit [K]' /
DATA argon(679) / '0.0                !(dummy) upper pressure limit' /
DATA argon(680) / '0.0                !(dummy) maximum density' /
DATA argon(681) / '83.8058  68.891    !reducing temperature and pressure' /
DATA argon(682) / '5 0 0 0 0 0                 !number of terms in melting line equation' /
DATA argon(683) / ' 1.             0.          !coefficients and exponents' /
DATA argon(684) / '-7476.26651     1.05' /
DATA argon(685) / ' 9959.06125     1.275' /
DATA argon(686) / ' 7476.26651     0.' /
DATA argon(687) / '-9959.06125     0.' /
DATA argon(688) / '' /
DATA argon(689) / '' /
DATA argon(690) / '#SBL        !sublimation line specification' /
DATA argon(691) / 'SB3  sublimation line model of Lemmon (2002).' /
DATA argon(692) / '?LITERATURE REFERENCE \' /
DATA argon(693) / '?Lemmon, E.W., 2002.' /
DATA argon(694) / '?\' /
DATA argon(695) / '!end of info section' /
DATA argon(696) / '83.8058            !lower temperature limit [K]' /
DATA argon(697) / '83.8058            !upper temperature limit [K]' /
DATA argon(698) / '0.0                !(dummy) upper pressure limit' /
DATA argon(699) / '0.0                !(dummy) maximum density' /
DATA argon(700) / '83.8058  68.891    !reducing temperature and pressure' /
DATA argon(701) / '0 1 0 0 0 0                 !number of terms in sublimation line equation' /
DATA argon(702) / '-11.1307        1.          !coefficients and exponents' /
DATA argon(703) / '' /
DATA argon(704) / '' /
DATA argon(705) / '#PS         !vapor pressure equation' /
DATA argon(706) / 'PS5  vapor pressure equation of Tegeler et al. (1999).' /
DATA argon(707) / '?LITERATURE REFERENCE \' /
DATA argon(708) / '?See EOS' /
DATA argon(709) / '?\' /
DATA argon(710) / '!end of info section' /
DATA argon(711) / '83.8058            !lower temperature limit [K]' /
DATA argon(712) / '150.687            !upper temperature limit [K]' /
DATA argon(713) / '0.0                !(dummy) upper pressure limit' /
DATA argon(714) / '0.0                !(dummy) maximum density' /
DATA argon(715) / '150.687 4863.0     !reducing parameters' /
DATA argon(716) / '4 0 0 0 0 0        !number of terms in equation' /
DATA argon(717) / '-5.9409785    1.0  !coefficients and exponents' /
DATA argon(718) / ' 1.3553888    1.5' /
DATA argon(719) / '-0.4649761    2.0' /
DATA argon(720) / '-1.5399043    4.5' /
DATA argon(721) / '' /
DATA argon(722) / '' /
DATA argon(723) / '#DL         !saturated liquid density equation' /
DATA argon(724) / 'DL3  saturated liquid density equation of Tegeler et al. (1999).' /
DATA argon(725) / '?LITERATURE REFERENCE \' /
DATA argon(726) / '?See EOS' /
DATA argon(727) / '?\' /
DATA argon(728) / '!end of info section' /
DATA argon(729) / '83.8058            !lower temperature limit [K]' /
DATA argon(730) / '150.687            !upper temperature limit [K]' /
DATA argon(731) / '0.0                !(dummy) upper pressure limit' /
DATA argon(732) / '0.0                !(dummy) maximum density' /
DATA argon(733) / '150.687 13.40742965 !reducing parameters' /
DATA argon(734) / '4 0 0 0 0 0         !number of terms in equation' /
DATA argon(735) / ' 1.5004264    0.334 !coefficients and exponents' /
DATA argon(736) / '-0.3138129    0.6666666666666' /
DATA argon(737) / ' 0.086461622  2.3333333333333' /
DATA argon(738) / '-0.041477525  4.0' /
DATA argon(739) / '' /
DATA argon(740) / '' /
DATA argon(741) / '#DV         !saturated vapor density equation' /
DATA argon(742) / 'DV5  saturated vapor density equation of Lemmon (2010).' /
DATA argon(743) / '?LITERATURE REFERENCE \' /
DATA argon(744) / '?Equation of Tegeler appears to be wrong, and new equation was fitted here.' /
DATA argon(745) / '?\' /
DATA argon(746) / '!end of info section' /
DATA argon(747) / '83.8058            !lower temperature limit [K]' /
DATA argon(748) / '150.687            !upper temperature limit [K]' /
DATA argon(749) / '0.0                !(dummy) upper pressure limit' /
DATA argon(750) / '0.0                !(dummy) maximum density' /
DATA argon(751) / '150.687 13.40742965   !reducing parameters' /
DATA argon(752) / '4 0 0 0 0 0           !number of terms in equation' /
DATA argon(753) / '-0.29182D+01   0.72' /
DATA argon(754) / ' 0.97930D-01   1.25' /
DATA argon(755) / '-0.13721D+01   0.32' /
DATA argon(756) / '-0.22898D+01   4.34' /
DATA argon(757) / '' /
DATA argon(758) / '' /
DATA argon(759) / '@END' /
DATA argon(760) / 'c        1         2         3         4         5         6         7         8' /
DATA argon(761) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
DATA argon(762) / '' /
DATA argon(763) / '' /
DATA argon(764) / '@TCX               !thermal conductivity model specification' /
DATA argon(765) / 'TC1  pure fluid thermal conductivity model of Younglove and Hanley (1986).' /
DATA argon(766) / '?REFERENCE \' /
DATA argon(767) / '?Younglove, B.A. and Hanley, H.J.M.,' /
DATA argon(768) / '? "The Viscosity and Thermal Conductivity Coefficients of Gaseous and' /
DATA argon(769) / '? Liquid Argon,"' /
DATA argon(770) / '? J. Phys. Chem. Ref. Data, 15(4):1323-1337, 1986.' /
DATA argon(771) / '?\' /
DATA argon(772) / '!end of info section' /
DATA argon(773) / '83.8058            !lower temperature limit [K]' /
DATA argon(774) / '700.0              !upper temperature limit [K]' /
DATA argon(775) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(776) / '50.65              !maximum density [mol/L]' /
DATA argon(777) / '9   0              !# terms for dilute gas function:  numerator, denominator' /
DATA argon(778) / '1.0     1.0d-3     !reducing parameters for T, tcx' /
DATA argon(779) / '-0.6700976192d5   -1.00d0   !coeff, power in T' /
DATA argon(780) / ' 0.6152255283d5   -0.666666666667d0' /
DATA argon(781) / '-0.2049218286d5   -0.333333333333d0' /
DATA argon(782) / ' 0.2216966254d4    0.00d0' /
DATA argon(783) / ' 0.3579189325d3    0.333333333333d0' /
DATA argon(784) / '-0.1364658914d3    0.666666666667d0' /
DATA argon(785) / ' 0.1718671649d2    1.00d0' /
DATA argon(786) / '-0.1018933154d1    1.333333333333d0' /
DATA argon(787) / ' 0.2397996932d-1   1.666666666667d0' /
DATA argon(788) / '11  3              !# terms for background gas function:  numerator, denominator' /
DATA argon(789) / '1.0    1.0    1.0d-3                      !reducing par for T, rho (rho_c), tcx' /
DATA argon(790) / ' 0.1536300190d1 0.0  1.0  0.0 !coeff, powers of T, rho, spare for future use' /
DATA argon(791) / '-0.2332533199d3   -1.00d0   1.00d0   0.' /
DATA argon(792) / '-0.3027085824d-1   0.00d0   2.00d0   0.' /
DATA argon(793) / ' 0.1896279196d2   -1.00d0   2.00d0   0.' /
DATA argon(794) / ' 0.1054230664d2   -2.00d0   2.00d0   0.' /
DATA argon(795) / ' 0.2588139028d-4   0.00d0   3.00d0   0.' /
DATA argon(796) / '-0.4546798772d0   -1.00d0   3.00d0   0.' /
DATA argon(797) / ' 0.4320206998d1   -2.00d0   3.00d0   0.' /
DATA argon(798) / ' 0.1593643304d-4   0.00d0   4.00d0   0.' /
DATA argon(799) / ' 0.1262253904d-3  -1.00d0   4.00d0   0.' /
DATA argon(800) / '-0.2937213042d-2  -2.00d0   4.00d0   0.' /
DATA argon(801) / ' 1.0000000000      0.00d0   0.00d0   0.' /
DATA argon(802) / '-0.2262773007d-1   0.00d0   1.00d0   0.' /
DATA argon(803) / '-0.1445619495d0   -1.00d0   1.00d0   0.' /
DATA argon(804) / 'TK4                !pointer to critical enhancement auxiliary function' /
DATA argon(805) / '' /
DATA argon(806) / '' /
DATA argon(807) / '#AUX               !thermal conductivity critical enhancement model' /
DATA argon(808) / 'TK4  thermal conductivity critical enhancement' /
DATA argon(809) / '?LITERATURE REFERENCE \' /
DATA argon(810) / '?Younglove, B.A. and Hanley, H.J.M.,' /
DATA argon(811) / '? "The Viscosity and Thermal Conductivity Coefficients of Gaseous and' /
DATA argon(812) / '? Liquid Argon,"' /
DATA argon(813) / '? J. Phys. Chem. Ref. Data, 15(4):1323-1337, 1986.' /
DATA argon(814) / '?\' /
DATA argon(815) / '!end of info section' /
DATA argon(816) / '83.8058            !lower temperature limit [K]' /
DATA argon(817) / '700.0              !upper temperature limit [K]' /
DATA argon(818) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(819) / '50.65              !maximum density [mol/L]' /
DATA argon(820) / '6 0 0 0            !# terms' /
DATA argon(821) / '150.86 4905.8 13.410 1.d-3      !reducing parameters' /
DATA argon(822) / '1.02' /
DATA argon(823) / '1.380658d-23' /
DATA argon(824) / '0.46807' /
DATA argon(825) / '39.8' /
DATA argon(826) / '5.45' /
DATA argon(827) / '6.0795d-1' /
DATA argon(828) / '' /
DATA argon(829) / '' /
DATA argon(830) / '#TCX               !thermal conductivity model specification' /
DATA argon(831) / 'TC1  pure fluid thermal conductivity model of Perkins et al. (1991).' /
DATA argon(832) / '?LITERATURE REFERENCE \' /
DATA argon(833) / '?Perkins, R.A., Friend, D.G., Roder, H.M., and Nieto de Castro, C.A.,' /
DATA argon(834) / '? "Thermal Conductivity Surface of Argon:  A Fresh Analysis,"' /
DATA argon(835) / '? Int. J. Thermophys., 12(6):965-984, 1991.' /
DATA argon(836) / '?\' /
DATA argon(837) / '?The uncertainty in thermal conductivity is 2.2%.' /
DATA argon(838) / '?\' /
DATA argon(839) / '!end of info section' /
DATA argon(840) / '55.0               !lower temperature limit [K]' /
DATA argon(841) / '2500.0             !upper temperature limit [K]' /
DATA argon(842) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(843) / '50.65              !maximum density [mol/L]' /
DATA argon(844) / '9   0              !# terms for dilute gas function:  numerator, denominator' /
DATA argon(845) / '1.0     1.0d-3     !reducing parameters for T, tcx' /
DATA argon(846) / ' .1225067272d+5 -1.00d0   !coeff, power in T' /
DATA argon(847) / '-.9096222831d+4 -0.66666666666666d0' /
DATA argon(848) / ' .2744958263d+4 -0.33333333333333d0' /
DATA argon(849) / '-.4170419051d+3  0.00d0' /
DATA argon(850) / ' .2527591169d+2  0.33333333333333d0' /
DATA argon(851) / ' .1604421067d+1  0.66666666666666d0' /
DATA argon(852) / '-.2618841031d+0  1.00d0' /
DATA argon(853) / ' .1381696924d-1  1.33333333333333d0' /
DATA argon(854) / '-.2463115922d-3  1.66666666666666d0' /
DATA argon(855) / '4   0              !# terms for background gas function:  numerator, denominator' /
DATA argon(856) / '1.0      1.0      1.0           !reducing par for T, rho (rho_c), tcx' /
DATA argon(857) / ' 0.757894d-3    0.0  1.0  0.0 !coeff, powers of T, rho, spare for future use' /
DATA argon(858) / ' 0.612624d-4   0.   2.   0.' /
DATA argon(859) / '-0.205353d-5   0.   3.   0.' /
DATA argon(860) / ' 0.745621d-7   0.   4.   0.' /
DATA argon(861) / 'TK4                !pointer to critical enhancement auxiliary function' /
DATA argon(862) / '' /
DATA argon(863) / '' /
DATA argon(864) / '#AUX               !thermal conductivity critical enhancement model' /
DATA argon(865) / 'TK4  thermal conductivity critical enhancement' /
DATA argon(866) / '?LITERATURE REFERENCE \' /
DATA argon(867) / '?Younglove, B.A. and Hanley, H.J.M.,' /
DATA argon(868) / '? "The Viscosity and Thermal Conductivity Coefficients of Gaseous and' /
DATA argon(869) / '? Liquid Argon,"' /
DATA argon(870) / '? J. Phys. Chem. Ref. Data, 15(4):1323-1337, 1986.' /
DATA argon(871) / '?\' /
DATA argon(872) / '!end of info section' /
DATA argon(873) / '55.0               !lower temperature limit [K]' /
DATA argon(874) / '700.0              !upper temperature limit [K]' /
DATA argon(875) / '1000000.0          !upper pressure limit [kPa]' /
DATA argon(876) / '50.65              !maximum density [mol/L]' /
DATA argon(877) / '6 0 0 0            !# terms' /
DATA argon(878) / '150.86 4905.8 13.410 1.d-3      !reducing parameters' /
DATA argon(879) / '1.02' /
DATA argon(880) / '1.380658d-23' /
DATA argon(881) / '0.46807' /
DATA argon(882) / '39.8' /
DATA argon(883) / '5.45' /
DATA argon(884) / '6.0795d-1' /
! #######################################################
character(256), TARGET :: f106_98_9____1butene(126)
DATA f106_98_9____1butene(1) / '106-98-9    1butene' /
DATA f106_98_9____1butene(2) / '64-19-7     acacid' /
DATA f106_98_9____1butene(3) / '67-64-1     acetone' /
DATA f106_98_9____1butene(4) / '1           air' /
DATA f106_98_9____1butene(5) / '7664-41-7   ammonia' /
DATA f106_98_9____1butene(6) / '7440-37-1   argon' /
DATA f106_98_9____1butene(7) / '71-43-2     benzene' /
DATA f106_98_9____1butene(8) / '106-97-8    butane' /
DATA f106_98_9____1butene(9) / '1120-21-4   c11' /
DATA f106_98_9____1butene(10) / '112-40-3    chlorine' /
DATA f106_98_9____1butene(11) / '108-87-2    c1cc6' /
DATA f106_98_9____1butene(12) / '590-18-1    c2butene' /
DATA f106_98_9____1butene(13) / '1678-92-8   c3cc6' /
DATA f106_98_9____1butene(14) / '2314-97-8   cf3i' /
DATA f106_98_9____1butene(15) / '7782-50-5   cl2' /
DATA f106_98_9____1butene(16) / '630-08-0    co' /
DATA f106_98_9____1butene(17) / '124-38-9    co2' /
DATA f106_98_9____1butene(18) / '463-58-1    cos' /
DATA f106_98_9____1butene(19) / '110-82-7    cyclohex' /
DATA f106_98_9____1butene(20) / '287-92-3    cyclopen' /
DATA f106_98_9____1butene(21) / '75-19-4     cyclopro' /
DATA f106_98_9____1butene(22) / '7782-39-0   d2' /
DATA f106_98_9____1butene(23) / '7789-20-0   d2o' /
DATA f106_98_9____1butene(24) / '556-67-2    d4' /
DATA f106_98_9____1butene(25) / '541-02-6    d5' /
DATA f106_98_9____1butene(26) / '540-97-6    d6' /
DATA f106_98_9____1butene(27) / '111-42-2    dea' /
DATA f106_98_9____1butene(28) / '124-18-5    decane' /
DATA f106_98_9____1butene(29) / '60-29-7     dee' /
DATA f106_98_9____1butene(30) / '616-38-6    dmc' /
DATA f106_98_9____1butene(31) / '115-10-6    dme' /
DATA f106_98_9____1butene(32) / '100-41-4    ebenzene' /
DATA f106_98_9____1butene(33) / '1300-21-6   edc' /
DATA f106_98_9____1butene(34) / '74-84-0     ethane' /
DATA f106_98_9____1butene(35) / '64-17-5     ethanol' /
DATA f106_98_9____1butene(36) / '74-85-1     ethylene' /
DATA f106_98_9____1butene(37) / '75-21-8     etyoxide' /
DATA f106_98_9____1butene(38) / '7782-41-4   fluorine' /
DATA f106_98_9____1butene(39) / '7783-06-4   h2s' /
DATA f106_98_9____1butene(40) / '7647-01-0   hcl' /
DATA f106_98_9____1butene(41) / '7440-59-7   helium' /
DATA f106_98_9____1butene(42) / '142-82-5    heptane' /
DATA f106_98_9____1butene(43) / '110-54-3    hexane' /
DATA f106_98_9____1butene(44) / '1333-74-0   hydrogen' /
DATA f106_98_9____1butene(45) / '115-11-7    ibutene' /
DATA f106_98_9____1butene(46) / '107-83-5    ihexane' /
DATA f106_98_9____1butene(47) / '540-84-1    ioctane' /
DATA f106_98_9____1butene(48) / '78-78-4     ipentane' /
DATA f106_98_9____1butene(49) / '75-28-5     isobutan' /
DATA f106_98_9____1butene(50) / '7439-90-9   krypton' /
DATA f106_98_9____1butene(51) / '141-43-5    mea' /
DATA f106_98_9____1butene(52) / '141-62-8    md2m' /
DATA f106_98_9____1butene(53) / '141-63-9    md3m' /
DATA f106_98_9____1butene(54) / '107-52-8    md4m' /
DATA f106_98_9____1butene(55) / '107-51-7    mdm' /
DATA f106_98_9____1butene(56) / '74-82-8     methane' /
DATA f106_98_9____1butene(57) / '67-56-1     methanol' /
DATA f106_98_9____1butene(58) / '112-63-0    mlinolea' /
DATA f106_98_9____1butene(59) / '301-00-8    mlinolen' /
DATA f106_98_9____1butene(60) / '107-46-0    mm' /
DATA f106_98_9____1butene(61) / '112-62-9    moleate' /
DATA f106_98_9____1butene(62) / '112-39-0    mpalmita' /
DATA f106_98_9____1butene(63) / '112-61-8    mstearat' /
DATA f106_98_9____1butene(64) / '108-38-3    mxylene' /
DATA f106_98_9____1butene(65) / '10024-97-2  n2o' /
DATA f106_98_9____1butene(66) / '7440-01-9   neon' /
DATA f106_98_9____1butene(67) / '463-82-1    neopentn' /
DATA f106_98_9____1butene(68) / '7783-54-2   nf3' /
DATA f106_98_9____1butene(69) / '7727-37-9   nitrogen' /
DATA f106_98_9____1butene(70) / '111-84-2    nonane' /
DATA f106_98_9____1butene(71) / '111-65-9    octane' /
DATA f106_98_9____1butene(72) / '1333-74-0o  orthohyd' /
DATA f106_98_9____1butene(73) / '7782-44-7   oxygen' /
DATA f106_98_9____1butene(74) / '95-47-6     oxylene' /
DATA f106_98_9____1butene(75) / '1333-74-0p  parahyd' /
DATA f106_98_9____1butene(76) / '109-66-0    pentane' /
DATA f106_98_9____1butene(77) / '74-98-6     propane' /
DATA f106_98_9____1butene(78) / '115-07-1    propylen' /
DATA f106_98_9____1butene(79) / '74-99-7     propyne' /
DATA f106_98_9____1butene(80) / '106-42-3    pxylene' /
DATA f106_98_9____1butene(81) / '75-69-4     r11' /
DATA f106_98_9____1butene(82) / '76-13-1     r113' /
DATA f106_98_9____1butene(83) / '76-14-2     r114' /
DATA f106_98_9____1butene(84) / '76-15-3     r115' /
DATA f106_98_9____1butene(85) / '76-16-4     r116' /
DATA f106_98_9____1butene(86) / '75-71-8     r12' /
DATA f106_98_9____1butene(87) / '116-15-4    r1216' /
DATA f106_98_9____1butene(88) / '306-83-2    r123' /
DATA f106_98_9____1butene(89) / '102687-65-0 r1233zd' /
DATA f106_98_9____1butene(90) / '754-12-1    r1234yf' /
DATA f106_98_9____1butene(91) / '29118-24-9  r1234ze' /
DATA f106_98_9____1butene(92) / '2837-89-0   r124' /
DATA f106_98_9____1butene(93) / '354-33-6    r125' /
DATA f106_98_9____1butene(94) / '75-72-9     r13' /
DATA f106_98_9____1butene(95) / '811-97-2    r134a' /
DATA f106_98_9____1butene(96) / '75-73-0     r14' /
DATA f106_98_9____1butene(97) / '1717-00-6   r141b' /
DATA f106_98_9____1butene(98) / '75-68-3     r142b' /
DATA f106_98_9____1butene(99) / '420-46-2    r143a' /
DATA f106_98_9____1butene(100) / '75-37-6     r152a' /
DATA f106_98_9____1butene(101) / '353-36-6    r161' /
DATA f106_98_9____1butene(102) / '75-43-4     r21' /
DATA f106_98_9____1butene(103) / '76-19-7     r218' /
DATA f106_98_9____1butene(104) / '75-45-6     r22' /
DATA f106_98_9____1butene(105) / '431-89-0    r227ea' /
DATA f106_98_9____1butene(106) / '75-46-7     r23' /
DATA f106_98_9____1butene(107) / '431-63-0    r236ea' /
DATA f106_98_9____1butene(108) / '690-39-1    r236fa' /
DATA f106_98_9____1butene(109) / '679-86-7    r245ca' /
DATA f106_98_9____1butene(110) / '460-73-1    r245fa' /
DATA f106_98_9____1butene(111) / '75-10-5     r32' /
DATA f106_98_9____1butene(112) / '406-58-6    r365mfc' /
DATA f106_98_9____1butene(113) / '74-87-3     r40' /
DATA f106_98_9____1butene(114) / '593-53-3    r41' /
DATA f106_98_9____1butene(115) / '115-25-3    rc318' /
DATA f106_98_9____1butene(116) / '421-14-7    re143a' /
DATA f106_98_9____1butene(117) / '22410-44-2  re245cb2' /
DATA f106_98_9____1butene(118) / '1885-48-9   re245fa2' /
DATA f106_98_9____1butene(119) / '375-03-1    re347mcc' /
DATA f106_98_9____1butene(120) / '2551-62-4   sf6' /
DATA f106_98_9____1butene(121) / '7446-09-5   so2' /
DATA f106_98_9____1butene(122) / '624-64-6    t2butene' /
DATA f106_98_9____1butene(123) / '108-88-3    toluene' /
DATA f106_98_9____1butene(124) / '75-01-4     vnychlrd' /
DATA f106_98_9____1butene(125) / '7732-18-5   water' /
DATA f106_98_9____1butene(126) / '7440-63-3   xenon' /
! #######################################################
character(256), TARGET :: Chlorine(373)
DATA Chlorine(1) / 'Chlorine             !Short name' /
DATA Chlorine(2) / '7782-50-5            !CAS number' /
DATA Chlorine(3) / 'Chlorine             !Full name' /
DATA Chlorine(4) / 'Cl2                  !Chemical formula {Cl2}' /
DATA Chlorine(5) / 'Chlorine             !Synonym' /
DATA Chlorine(6) / '70.906               !Molar mass [g/mol]' /
DATA Chlorine(7) / '172.17               !Triple point temperature [K]' /
DATA Chlorine(8) / '239.198              !Normal boiling point [K]' /
DATA Chlorine(9) / '416.8654             !Critical temperature [K]' /
DATA Chlorine(10) / '7642.4               !Critical pressure [kPa]' /
DATA Chlorine(11) / '8.06                 !Critical density [mol/L]' /
DATA Chlorine(12) / '0.070                !Acentric factor' /
DATA Chlorine(13) / '0.                   !Dipole moment [Debye]; Reid, Prausnitz, & Poling, McGraw-Hill (1987)' /
DATA Chlorine(14) / 'NBP                  !Default reference state' /
DATA Chlorine(15) / '10.0                 !Version number' /
DATA Chlorine(16) / '1017                 !UN Number                                                 :UN:' /
DATA Chlorine(17) / 'other                !Family                                                    :Family:' /
DATA Chlorine(18) / '0.0                  !Heating value (upper) [kJ/mol]                            :Heat:' /
DATA Chlorine(19) / '1S/Cl2/c1-2                               !Standard InChI String                :InChi:' /
DATA Chlorine(20) / 'KZBUYRJDOAKODT-UHFFFAOYSA-N               !Standard InChI Key                   :InChiKey:' /
DATA Chlorine(21) / '????                                      !Alternative fluid for mixing rules   :AltID:' /
DATA Chlorine(22) / '828f1c80                                  !Hash number from InChI Key           :Hash:' /
DATA Chlorine(23) / '' /
DATA Chlorine(24) / '' /
DATA Chlorine(25) / '' /
DATA Chlorine(26) / '' /
DATA Chlorine(27) / '' /
DATA Chlorine(28) / '' /
DATA Chlorine(29) / '' /
DATA Chlorine(30) / '________________________________________________________________________________' /
DATA Chlorine(31) / '' /
DATA Chlorine(32) / '#EOS   !---Equation of state---' /
DATA Chlorine(33) / 'FEQ    !Helmholtz equation of state for chlorine of Herrig et al. (2018).' /
DATA Chlorine(34) / ':TRUECRITICALPOINT:  416.8654   7.949829      !True EOS critical point [K, mol/L] (where dP/dD=0 and d^2P/dD^2=0 at constant T)' /
DATA Chlorine(35) / ':DOI:' /
DATA Chlorine(36) / '?' /
DATA Chlorine(37) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(38) / '?Thol, M., Herrig, S., Span, R., and Lemmon, E.W.' /
DATA Chlorine(39) / '? to be submitted to Phys. Chem. Chem. Phys., 2020.' /
DATA Chlorine(40) / '?' /
DATA Chlorine(41) / '?In the liquid phase, homogeneous densities can be obtained from the EOS with an' /
DATA Chlorine(42) / '? uncertainty of 0.1 %, wtih slightly higher uncertainties close to the phase' /
DATA Chlorine(43) / '? boundary.  Homogeneous gas densities are represented within 0.25 %. The EOS' /
DATA Chlorine(44) / '? represents the most accurate experimental speed-of-sound data in the gas phase' /
DATA Chlorine(45) / '? with deviations within 0.05 %. Due to the limited experimental database, the' /
DATA Chlorine(46) / '? estimated uncertainty of calculated sound speeds in the liquid phase is 0.6 %.' /
DATA Chlorine(47) / '? The uncertainty of calculated vapor pressures is 0.4 % at temperatures up to' /
DATA Chlorine(48) / '? 275 K and 1 % at higher temperatures up 320 K. The uncertainties for all' /
DATA Chlorine(49) / '? properties increase at higher temperatures where there are no reliable' /
DATA Chlorine(50) / '? experimental data.' /
DATA Chlorine(51) / '?' /
DATA Chlorine(52) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(53) / '172.17             !Lower temperature limit [K]' /
DATA Chlorine(54) / '440.               !Upper temperature limit [K]' /
DATA Chlorine(55) / '20000.             !Upper pressure limit [kPa]' /
DATA Chlorine(56) / '24.61              !Maximum density [mol/L]' /
DATA Chlorine(57) / 'CPP                                    !Pointer to Cp0 model' /
DATA Chlorine(58) / '70.906                                 !Molar mass [g/mol]' /
DATA Chlorine(59) / '172.17                                 !Triple point temperature [K]' /
DATA Chlorine(60) / '1.3808                                 !Pressure at triple point [kPa]' /
DATA Chlorine(61) / '24.6                                   !Density at triple point [mol/L]' /
DATA Chlorine(62) / '239.198                                !Normal boiling point temperature [K]' /
DATA Chlorine(63) / '0.070                                  !Acentric factor' /
DATA Chlorine(64) / '416.8654      7642.4       8.06        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Chlorine(65) / '416.8654                   8.06        !Reducing parameters [K, mol/L]' /
DATA Chlorine(66) / '8.314462618                            !Gas constant [J/mol-K]' /
DATA Chlorine(67) / '  10  4   5 12   0 0    0 0  0 0  0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Chlorine(68) / '  0.0245017    1.0     4.  0.          !a(i),t(i),d(i),l(i)' /
DATA Chlorine(69) / '  0.9132904    0.196   1.  0.' /
DATA Chlorine(70) / ' -1.72309      1.      1.  0.' /
DATA Chlorine(71) / ' -0.3359344    1.08    2.  0.' /
DATA Chlorine(72) / '  0.1200495    0.39    3.  0.' /
DATA Chlorine(73) / ' -1.214889     1.64    1.  2.' /
DATA Chlorine(74) / ' -0.10167      3.2     3.  2.' /
DATA Chlorine(75) / '  0.6196819    1.32    2.  1.' /
DATA Chlorine(76) / ' -0.6578512    2.163   2.  2.' /
DATA Chlorine(77) / ' -0.009159452  0.93    7.  1.' /
DATA Chlorine(78) / '  1.909418     0.872   1.  2. 2.    -0.969    -1.22    1.142   0.88     0. 0. 0.' /
DATA Chlorine(79) / ' -0.07163412   2.08    1.  2. 2.    -1.89     -6.8     1.22    0.73     0. 0. 0.' /
DATA Chlorine(80) / ' -0.1893345    1.6     3.  2. 2.    -1.32     -3.5     1.552   0.28     0. 0. 0.' /
DATA Chlorine(81) / ' -0.5698469    1.37    2.  2. 2.    -1.012    -1.276   1.135   0.863    0. 0. 0.' /
DATA Chlorine(82) / ' -0.8964496    1.05    2.  2. 2.    -0.98     -1.6     0.754   0.554    0. 0. 0.' /
DATA Chlorine(83) / '                                      eta      beta    gamma   epsilon' /
DATA Chlorine(84) / '                                   EXP[eta*(delta-epsilon)^2+beta*(tau-gamma)^2]' /
DATA Chlorine(85) / '' /
DATA Chlorine(86) / '' /
DATA Chlorine(87) / '#AUX   !---Auxiliary function for Cp0' /
DATA Chlorine(88) / 'CPP    !Ideal gas heat capacity function for chlorine of Thol et al. (2020).' /
DATA Chlorine(89) / '?' /
DATA Chlorine(90) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(91) / '?Herrig, S., Thol, M., and Span, R., 2018.' /
DATA Chlorine(92) / '?' /
DATA Chlorine(93) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(94) / '0.                 !' /
DATA Chlorine(95) / '10000.             !' /
DATA Chlorine(96) / '0.                 !' /
DATA Chlorine(97) / '0.                 !' /
DATA Chlorine(98) / '1.0    8.314462618 !Reducing parameters for T, Cp0' /
DATA Chlorine(99) / '1 3   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Chlorine(100) / ' 3.5        0.0' /
DATA Chlorine(101) / ' 1.0256     800.0' /
DATA Chlorine(102) / ' 0.067756   3000.0' /
DATA Chlorine(103) / ' 0.14068    8200.0' /
DATA Chlorine(104) / '' /
DATA Chlorine(105) / '' /
DATA Chlorine(106) / '#AUX   !---Auxiliary function for PX0' /
DATA Chlorine(107) / 'PX0    !Helmholtz energy ideal-gas function for chlorine of Thol et al. (2020).' /
DATA Chlorine(108) / '?' /
DATA Chlorine(109) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(110) / '?Thol, M., Herrig, S., Span, R., and Lemmon, E.W. (2020)' /
DATA Chlorine(111) / '?' /
DATA Chlorine(112) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(113) / '1 2  3  0 0  0 0 0               !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau))' /
DATA Chlorine(114) / '  2.5                   1.0      !ai, ti for [ai*log(tau**ti)] terms' /
DATA Chlorine(115) / ' -3.95390162183112      0.0      !aj, ti for [ai*tau**ti] terms' /
DATA Chlorine(116) / '  3.83990498815427      1.0      !aj, ti for [ai*tau**ti] terms' /
DATA Chlorine(117) / '  1.0256     800.0               !aj, ti for [ai*log(1-exp(-ti/T)] terms' /
DATA Chlorine(118) / '  0.067756   3000.0' /
DATA Chlorine(119) / '  0.14068    8200.0' /
DATA Chlorine(120) / '' /
DATA Chlorine(121) / '' /
DATA Chlorine(122) / '--------------------------------------------------------------------------------' /
DATA Chlorine(123) / '' /
DATA Chlorine(124) / '@EOS    !---Equation of state---' /
DATA Chlorine(125) / 'FE1     !Helmholtz equation of state for chlorine of Bonsen (2002).' /
DATA Chlorine(126) / '          ?' /
DATA Chlorine(127) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(128) / '          ?Bonsen (2002)' /
DATA Chlorine(129) / '          ?' /
DATA Chlorine(130) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(131) / '          172.12             !Lower temperature limit [K]' /
DATA Chlorine(132) / '          600.               !Upper temperature limit [K]' /
DATA Chlorine(133) / '          100000.0           !Upper pressure limit [kPa]' /
DATA Chlorine(134) / '          24.45              !Maximum density [mol/L]' /
DATA Chlorine(135) / '          CP1                                    !Pointer to Cp0 model' /
DATA Chlorine(136) / '          70.906                                 !Molar mass [g/mol]' /
DATA Chlorine(137) / '          172.12                                 !Triple point temperature [K]' /
DATA Chlorine(138) / '          1.4                                    !Pressure at triple point [kPa]' /
DATA Chlorine(139) / '          24.4                                   !Density at triple point [mol/L]' /
DATA Chlorine(140) / '          239.1                                  !Normal boiling point temperature [K]' /
DATA Chlorine(141) / '          0.174                                  !Acentric factor' /
DATA Chlorine(142) / '          417.0         9700.0       8.1375342   !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Chlorine(143) / '          417.0                      8.1375342   !Reducing parameters [K, mol/L]' /
DATA Chlorine(144) / '          8.314                                  !Gas constant [J/mol-K]' /
DATA Chlorine(145) / '            12  4    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Chlorine(146) / '           0.9178215309        0.25      1.  0.  !a(i),t(i),d(i),l(i)' /
DATA Chlorine(147) / '          -2.301431792         1.125     1.  0.' /
DATA Chlorine(148) / '           0.4517970045        1.5       1.  0.' /
DATA Chlorine(149) / '          -0.04310372015       1.375     2.  0.' /
DATA Chlorine(150) / '           0.07777809415       0.25      3.  0.' /
DATA Chlorine(151) / '           0.0001900921704     0.875     7.  0.' /
DATA Chlorine(152) / '           0.006949800797      0.625     2.  1.' /
DATA Chlorine(153) / '          -0.02150536297       1.75      5.  1.' /
DATA Chlorine(154) / '          -0.2178807183        3.625     1.  2.' /
DATA Chlorine(155) / '           0.0168814077        3.625     4.  2.' /
DATA Chlorine(156) / '           0.01549021349      14.5       3.  3.' /
DATA Chlorine(157) / '           0.01722723733      12.0       4.  3.' /
DATA Chlorine(158) / '' /
DATA Chlorine(159) / '' /
DATA Chlorine(160) / '@AUX    !---Auxiliary function for Cp0' /
DATA Chlorine(161) / 'CP1     !Ideal gas heat capacity function for chlorine.' /
DATA Chlorine(162) / '          ?' /
DATA Chlorine(163) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(164) / '          ?Bonsen' /
DATA Chlorine(165) / '          ?' /
DATA Chlorine(166) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(167) / '          0.                 !' /
DATA Chlorine(168) / '          10000.             !' /
DATA Chlorine(169) / '          0.                 !' /
DATA Chlorine(170) / '          0.                 !' /
DATA Chlorine(171) / '          1.0     8.314      !Reducing parameters for T, Cp0' /
DATA Chlorine(172) / '          4 1   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Chlorine(173) / '           2.9711935         0.' /
DATA Chlorine(174) / '           0.00498353462     1.' /
DATA Chlorine(175) / '          -0.000004177662906 2.' /
DATA Chlorine(176) / '           0.000000001318031484 3.' /
DATA Chlorine(177) / '          -0.812960306       2014.02967' /
DATA Chlorine(178) / '' /
DATA Chlorine(179) / '' /
DATA Chlorine(180) / '@AUX    !---Auxiliary function for Cp0' /
DATA Chlorine(181) / 'CP2     !Ideal gas heat capacity function for chlorine.' /
DATA Chlorine(182) / '          ?' /
DATA Chlorine(183) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(184) / '          ?Martin and Longpre, JCED, 29:466-473, 1984.' /
DATA Chlorine(185) / '          ?' /
DATA Chlorine(186) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(187) / '          0.                 !' /
DATA Chlorine(188) / '          10000.             !' /
DATA Chlorine(189) / '          0.                 !' /
DATA Chlorine(190) / '          0.                 !' /
DATA Chlorine(191) / '          1.0     4.184      !Reducing parameters for T, Cp0' /
DATA Chlorine(192) / '          4 0   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Chlorine(193) / '           10.265952         0.0' /
DATA Chlorine(194) / '          -0.00078085907     1.0' /
DATA Chlorine(195) / '          -709.60655        -1.0' /
DATA Chlorine(196) / '           40821.249        -2.0' /
DATA Chlorine(197) / '' /
DATA Chlorine(198) / '' /
DATA Chlorine(199) / '' /
DATA Chlorine(200) / '' /
DATA Chlorine(201) / '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' /
DATA Chlorine(202) / '' /
DATA Chlorine(203) / '#TRN   !---ECS Transport---' /
DATA Chlorine(204) / 'ECS    !Extended Corresponding States model (Propane reference) fit to limited data for chlorine.' /
DATA Chlorine(205) / ':DOI: 10.6028/NIST.IR.8209' /
DATA Chlorine(206) / '?' /
DATA Chlorine(207) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(208) / '?Huber, M.L., "Models for the Viscosity, Thermal Conductivity, and Surface Tension' /
DATA Chlorine(209) / '? of Selected Pure Fluids as Implemented in REFPROP v10.0," NISTIR 8209, 2018.' /
DATA Chlorine(210) / '? doi: 10.6028/NIST.IR.8209' /
DATA Chlorine(211) / '?' /
DATA Chlorine(212) / '?THERMAL CONDUCTIVITY' /
DATA Chlorine(213) / '? Chaikin, A.M. and Markevich, A.M., Zh. Fiz. Khim., 32:116-20, 1958.' /
DATA Chlorine(214) / '? Ho, C.Y., Powell, R.W., and Liley, P.E., "Thermal Conductivity of the Elements," J. Phys. Chem. Ref. Data, 1(2):279, 1972.' /
DATA Chlorine(215) / '?' /
DATA Chlorine(216) / '?Estimated uncertainty in the gas phase is <10% based on comparisons with the data' /
DATA Chlorine(217) / '? of Chaikin and Markevich.  Estimated uncertainty in the liquid phase is difficult' /
DATA Chlorine(218) / '? to assess due to lack of experimental data, estimated to be <10% along saturation' /
DATA Chlorine(219) / '? based on agreement with recommended values of Ho et al.' /
DATA Chlorine(220) / '?' /
DATA Chlorine(221) / '?VISCOSITY' /
DATA Chlorine(222) / '? Steacie, E.W.R. and Johnson, F.M.G., "The Viscosities of the Liquid Halogens," J. Am. Chem. Soc., 47:754-762, 1925.' /
DATA Chlorine(223) / '?' /
DATA Chlorine(224) / '?Estimated uncertainty in the gas phase is <10%, based on comparisons with the' /
DATA Chlorine(225) / '? data of Trautz, M. and Ruf, F., Ann. Phys., 20(5):127, 1934.  Estimated uncertainty' /
DATA Chlorine(226) / '? along the liquid saturation line is <10%, based on comparisons with the data of Steacie and Johnson.' /
DATA Chlorine(227) / '?' /
DATA Chlorine(228) / '?The Lennard-Jones parameters were taken from Hirschfelder, J.O., Curtiss, C.F., and Bird, R.B., "Molecular Theory of Gases and Liquids," John Wiley and Sons, Inc., New York, 1245 pp, 1954. doi: 10.1002/pol.1955.120178311' /
DATA Chlorine(229) / '?' /
DATA Chlorine(230) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(231) / '172.17             !Lower temperature limit [K]' /
DATA Chlorine(232) / '440.0              !Upper temperature limit [K]' /
DATA Chlorine(233) / '20000.0            !Upper pressure limit [kPa]' /
DATA Chlorine(234) / '24.61              !Maximum density [mol/L]' /
DATA Chlorine(235) / 'FEQ PROPANE.FLD' /
DATA Chlorine(236) / 'VS1                !Model for reference fluid viscosity' /
DATA Chlorine(237) / 'TC1                !Model for reference fluid thermal conductivity' /
DATA Chlorine(238) / 'NUL                !Large molecule identifier' /
DATA Chlorine(239) / '1                  !Lennard-Jones flag (0 or 1) (0 => use estimates)' /
DATA Chlorine(240) / '0.44               !Lennard-Jones coefficient sigma [nm]' /
DATA Chlorine(241) / '257.0              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA Chlorine(242) / '1  0  0                  !Number of terms in f_int term in Eucken correlation, spare1, spare2' /
DATA Chlorine(243) / ' 0.0029        0. 0. 0.  !Coefficient, power of T, spare1, spare2' /
DATA Chlorine(244) / '2  0  0                  !Number of terms in psi (visc shape factor): poly,spare1,spare2' /
DATA Chlorine(245) / ' 1.269         0. 0. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Chlorine(246) / '-0.08947       0. 1. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Chlorine(247) / '2  0  0                  !Number of terms in chi (t.c. shape factor): poly,spare1,spare2' /
DATA Chlorine(248) / ' 1.24341       0. 0. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Chlorine(249) / '-0.0812555     0. 1. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Chlorine(250) / 'TK3                !Pointer to critical enhancement auxiliary function' /
DATA Chlorine(251) / '' /
DATA Chlorine(252) / '' /
DATA Chlorine(253) / '#AUX   !---Auxiliary function for the thermal conductivity critical enhancement' /
DATA Chlorine(254) / 'TK3    !Simplified thermal conductivity critical enhancement for chlorine of Perkins et al. (2013).' /
DATA Chlorine(255) / '?' /
DATA Chlorine(256) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(257) / '?Perkins, R.A., Sengers, J.V., Abdulagatov, I.M., and Huber, M.L.,' /
DATA Chlorine(258) / '? "Simplified Model for the Critical Thermal-Conductivity Enhancement in Molecular Fluids,"' /
DATA Chlorine(259) / '? Int. J. Thermophys., 34(2):191-212, 2013. doi: 10.1007/s10765-013-1409-z' /
DATA Chlorine(260) / '?' /
DATA Chlorine(261) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(262) / '0.                 !' /
DATA Chlorine(263) / '10000.             !' /
DATA Chlorine(264) / '0.                 !' /
DATA Chlorine(265) / '0.                 !' /
DATA Chlorine(266) / '9 0 0 0            !# terms:  terms, spare, spare, spare' /
DATA Chlorine(267) / '1.0  1.0  1.0      !Reducing parameters for T, rho, tcx [mW/(m-K)]' /
DATA Chlorine(268) / '0.63               !Nu (universal exponent)' /
DATA Chlorine(269) / '1.239              !Gamma (universal exponent)' /
DATA Chlorine(270) / '1.02               !R0 (universal amplitude)' /
DATA Chlorine(271) / '0.063              !Z (universal exponent--not used for t.c., only viscosity)' /
DATA Chlorine(272) / '1.0                !C (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA Chlorine(273) / '0.179e-9           !Xi0 (amplitude) [m]' /
DATA Chlorine(274) / '0.056              !Gam0 (amplitude) [-]' /
DATA Chlorine(275) / '0.486e-9           !Qd_inverse (modified effective cutoff parameter) [m]; generic number, not fitted to data' /
DATA Chlorine(276) / '625.30             !Tref (reference temperature)=1.5*Tc [K]' /
DATA Chlorine(277) / '' /
DATA Chlorine(278) / '' /
DATA Chlorine(279) / '' /
DATA Chlorine(280) / '' /
DATA Chlorine(281) / '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' /
DATA Chlorine(282) / '' /
DATA Chlorine(283) / '#STN   !---Surface tension---' /
DATA Chlorine(284) / 'ST1    !Surface tension model for chlorine of Huber (2018).' /
DATA Chlorine(285) / ':DOI: 10.6028/NIST.IR.8209' /
DATA Chlorine(286) / '?' /
DATA Chlorine(287) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(288) / '?Huber, M.L., "Models for the Viscosity, Thermal Conductivity, and Surface Tension' /
DATA Chlorine(289) / '? of Selected Pure Fluids as Implemented in REFPROP v10.0," NISTIR 8209, 2018.' /
DATA Chlorine(290) / '? doi: 10.6028/NIST.IR.8209' /
DATA Chlorine(291) / '?' /
DATA Chlorine(292) / '?Fit to experimental data of:' /
DATA Chlorine(293) / '? Johnson, F.M.G. and McIntosh, D., "Liquid Chlorine," J. Am. Chem. Soc., 31(10):1138-1144, 1909.' /
DATA Chlorine(294) / '?' /
DATA Chlorine(295) / '?Estimated uncertainty is 5%.' /
DATA Chlorine(296) / '?' /
DATA Chlorine(297) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(298) / '0.                 !' /
DATA Chlorine(299) / '10000.             !' /
DATA Chlorine(300) / '0.                 !' /
DATA Chlorine(301) / '0.                 !' /
DATA Chlorine(302) / '1                  !Number of terms in surface tension model' /
DATA Chlorine(303) / '416.8654           !Critical temperature used in fit (dummy)' /
DATA Chlorine(304) / '0.0783601 1.28083  !Sigma0 and n' /
DATA Chlorine(305) / '' /
DATA Chlorine(306) / '' /
DATA Chlorine(307) / '#PS    !---Vapor pressure---' /
DATA Chlorine(308) / 'PS5    !Vapor pressure equation for chlorine of Thol et al. (2020).' /
DATA Chlorine(309) / '?' /
DATA Chlorine(310) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(311) / '?Functional Form:  P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T] where Theta=1-T/Tc, Tc and Pc' /
DATA Chlorine(312) / '? are the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Chlorine(313) / '?' /
DATA Chlorine(314) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(315) / '0.                 !' /
DATA Chlorine(316) / '10000.             !' /
DATA Chlorine(317) / '0.                 !' /
DATA Chlorine(318) / '0.                 !' /
DATA Chlorine(319) / '416.8654  7642.4   !Reducing parameters' /
DATA Chlorine(320) / '6 0 0 0 0 0        !Number of terms in equation' /
DATA Chlorine(321) / '-6.1289    1.0' /
DATA Chlorine(322) / ' 1.5112    1.5' /
DATA Chlorine(323) / '-1.4523    2.0' /
DATA Chlorine(324) / '-5.6038    5.94' /
DATA Chlorine(325) / ' 3.9923    7.0' /
DATA Chlorine(326) / '-1.2651    14.8' /
DATA Chlorine(327) / '' /
DATA Chlorine(328) / '' /
DATA Chlorine(329) / '#DL    !---Saturated liquid density---' /
DATA Chlorine(330) / 'DL1    !Saturated liquid density equation for chlorine of Thol et al. (2020).' /
DATA Chlorine(331) / '?' /
DATA Chlorine(332) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(333) / '?Functional Form:  D=Dc*[1+SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Chlorine(334) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Chlorine(335) / '?' /
DATA Chlorine(336) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(337) / '0.                 !' /
DATA Chlorine(338) / '10000.             !' /
DATA Chlorine(339) / '0.                 !' /
DATA Chlorine(340) / '0.                 !' /
DATA Chlorine(341) / '416.8654 8.06      !Reducing parameters' /
DATA Chlorine(342) / '4 0 0 0 0 0        !Number of terms in equation' /
DATA Chlorine(343) / ' 0.9662    0.234' /
DATA Chlorine(344) / ' 1.7744    0.68' /
DATA Chlorine(345) / '-0.23081   1.3' /
DATA Chlorine(346) / ' 0.47213   3.35' /
DATA Chlorine(347) / '' /
DATA Chlorine(348) / '' /
DATA Chlorine(349) / '#DV    !---Saturated vapor density---' /
DATA Chlorine(350) / 'DV3    !Saturated vapor density equation for chlorine of Thol et al. (2020).' /
DATA Chlorine(351) / '?' /
DATA Chlorine(352) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(353) / '?Functional Form:  D=Dc*EXP[SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Chlorine(354) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Chlorine(355) / '?' /
DATA Chlorine(356) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Chlorine(357) / '0.                 !' /
DATA Chlorine(358) / '10000.             !' /
DATA Chlorine(359) / '0.                 !' /
DATA Chlorine(360) / '0.                 !' /
DATA Chlorine(361) / '416.8654  8.06     !Reducing parameters' /
DATA Chlorine(362) / '6 0 0 0 0 0        !Number of terms in equation' /
DATA Chlorine(363) / '-1.7673    0.3' /
DATA Chlorine(364) / '-5.173     0.994' /
DATA Chlorine(365) / '-12.539    2.7' /
DATA Chlorine(366) / '-37.552    6.155' /
DATA Chlorine(367) / '-64.404    12.4' /
DATA Chlorine(368) / '-151.49    24.0' /
DATA Chlorine(369) / '' /
DATA Chlorine(370) / '' /
DATA Chlorine(371) / '@END' /
DATA Chlorine(372) / 'c        1         2         3         4         5         6         7         8' /
DATA Chlorine(373) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
! #######################################################
character(256), TARGET :: carbon_monoxide(324)
DATA carbon_monoxide(1) / 'carbon monoxide    !short name' /
DATA carbon_monoxide(2) / '630-08-0           !CAS number' /
DATA carbon_monoxide(3) / 'carbon monoxide    !full name' /
DATA carbon_monoxide(4) / 'CO                 !chemical formula' /
DATA carbon_monoxide(5) / 'carbon oxide       !synonym' /
DATA carbon_monoxide(6) / '28.0101            !molecular weight [g/mol]' /
DATA carbon_monoxide(7) / '68.16              !triple point temperature [K]' /
DATA carbon_monoxide(8) / '81.64              !normal boiling point [K]' /
DATA carbon_monoxide(9) / '132.86             !critical temperature [K]' /
DATA carbon_monoxide(10) / '3494.0             !critical pressure [kPa]' /
DATA carbon_monoxide(11) / '10.85              !critical density [mol/L]' /
DATA carbon_monoxide(12) / '0.0497             !acentric factor' /
DATA carbon_monoxide(13) / '0.1                !dipole moment [Debye]; Reid, Prausnitz, & Poling, McGraw-Hill (1987)' /
DATA carbon_monoxide(14) / 'NBP                !default reference state' /
DATA carbon_monoxide(15) / '8.0                !version number' /
DATA carbon_monoxide(16) / '1016               !UN Number' /
DATA carbon_monoxide(17) / 'other              !family' /
DATA carbon_monoxide(18) / '282.98             !heating value (gross or superior) [kJ/mol]' /
DATA carbon_monoxide(19) / '' /
DATA carbon_monoxide(20) / '' /
DATA carbon_monoxide(21) / '#EOS               !equation of state specification' /
DATA carbon_monoxide(22) / 'FEQ  short Helmholtz equation of state for carbon monoxide of Lemmon and Span (2006).' /
DATA carbon_monoxide(23) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(24) / '?Lemmon, E.W. and Span, R.,' /
DATA carbon_monoxide(25) / '? "Short Fundamental Equations of State for 20 Industrial Fluids,"' /
DATA carbon_monoxide(26) / '? J. Chem. Eng. Data, 51:785-850, 2006.' /
DATA carbon_monoxide(27) / '?\' /
DATA carbon_monoxide(28) / '?The equation of state is valid from the triple point to 500 K with' /
DATA carbon_monoxide(29) / '?pressures to 100 MPa. At higher pressures, the deviations from the equation' /
DATA carbon_monoxide(30) / '?increase rapidly and it is not recommended to use the equation above 100' /
DATA carbon_monoxide(31) / '?MPa. The uncertainties in the equation are 0.3% in density (approaching 1%' /
DATA carbon_monoxide(32) / '?near the critical point), 0.2% in vapor pressure, and 2% in heat' /
DATA carbon_monoxide(33) / '?capacities.  The uncertainty in the speed of sound is unknown.' /
DATA carbon_monoxide(34) / '?\' /
DATA carbon_monoxide(35) / '!end of info section' /
DATA carbon_monoxide(36) / '68.16              !lower temperature limit [K]' /
DATA carbon_monoxide(37) / '500.0              !upper temperature limit [K]' /
DATA carbon_monoxide(38) / '100000.0           !upper pressure limit [kPa]' /
DATA carbon_monoxide(39) / '33.84              !maximum density [mol/L]' /
DATA carbon_monoxide(40) / 'CPP                                    !pointer to Cp0 model' /
DATA carbon_monoxide(41) / '28.0101                                !molecular weight [g/mol]' /
DATA carbon_monoxide(42) / '68.16                                  !triple point temperature [K]' /
DATA carbon_monoxide(43) / '15.53                                  !pressure at triple point [kPa]' /
DATA carbon_monoxide(44) / '30.33                                  !density at triple point [mol/L]' /
DATA carbon_monoxide(45) / '81.64                                  !normal boiling point temperature [K]' /
DATA carbon_monoxide(46) / '0.0497                                 !acentric factor' /
DATA carbon_monoxide(47) / '132.86        3494.0     10.85         !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA carbon_monoxide(48) / '132.86                   10.85         !reducing parameters [K, mol/L]' /
DATA carbon_monoxide(49) / '8.314472                               !gas constant [J/mol-K]' /
DATA carbon_monoxide(50) / '  12  4      0  0      0  0            !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA carbon_monoxide(51) / '  0.90554         0.25    1.0   0      !a(i),t(i),d(i),l(i)' /
DATA carbon_monoxide(52) / ' -2.4515          1.125   1.0   0' /
DATA carbon_monoxide(53) / '  0.53149         1.5     1.0   0' /
DATA carbon_monoxide(54) / '  0.024173        1.375   2.0   0' /
DATA carbon_monoxide(55) / '  0.072156        0.25    3.0   0' /
DATA carbon_monoxide(56) / '  0.00018818      0.875   7.0   0' /
DATA carbon_monoxide(57) / '  0.19405         0.625   2.0   1' /
DATA carbon_monoxide(58) / ' -0.043268        1.75    5.0   1' /
DATA carbon_monoxide(59) / ' -0.12778         3.625   1.0   2' /
DATA carbon_monoxide(60) / ' -0.027896        3.625   4.0   2' /
DATA carbon_monoxide(61) / ' -0.034154       14.5     3.0   3' /
DATA carbon_monoxide(62) / '  0.016329       12.0     4.0   3' /
DATA carbon_monoxide(63) / '' /
DATA carbon_monoxide(64) / '' /
DATA carbon_monoxide(65) / '#AUX               !auxiliary model specification' /
DATA carbon_monoxide(66) / 'CPP  ideal gas heat capacity function' /
DATA carbon_monoxide(67) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(68) / '?Lemmon, E.W. and Span, R. (see eos for reference)' /
DATA carbon_monoxide(69) / '?\' /
DATA carbon_monoxide(70) / '!end of info section' /
DATA carbon_monoxide(71) / '50.0               !lower temperature limit [K]' /
DATA carbon_monoxide(72) / '5000.0             !upper temperature limit [K]' /
DATA carbon_monoxide(73) / '0.0                !upper pressure limit [kPa]' /
DATA carbon_monoxide(74) / '0.0                !maximum density [mol/L]' /
DATA carbon_monoxide(75) / '1.0          8.314472                  !reducing parameters for T, Cp0' /
DATA carbon_monoxide(76) / '  2  1    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA carbon_monoxide(77) / ' 3.5            0.0' /
DATA carbon_monoxide(78) / ' 0.22311E-6     1.5' /
DATA carbon_monoxide(79) / ' 1.0128      3089.0' /
DATA carbon_monoxide(80) / '' /
DATA carbon_monoxide(81) / '' /
DATA carbon_monoxide(82) / '#AUX               !auxiliary model specification' /
DATA carbon_monoxide(83) / 'PH0  Helmholtz form for the ideal-gas state' /
DATA carbon_monoxide(84) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(85) / '?Lemmon, E.W. and Span, R. (see eos for reference)' /
DATA carbon_monoxide(86) / '?\' /
DATA carbon_monoxide(87) / '!end of info section' /
DATA carbon_monoxide(88) / '50.0               !lower temperature limit [K]' /
DATA carbon_monoxide(89) / '5000.0             !upper temperature limit [K]' /
DATA carbon_monoxide(90) / '0.0                !upper pressure limit [kPa]' /
DATA carbon_monoxide(91) / '0.0                !maximum density [mol/L]' /
DATA carbon_monoxide(92) / '1 3  1  0 0  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau))' /
DATA carbon_monoxide(93) / '    2.5000000000    1.0000000000   !ai, ti for [ai*log(tau**ti)] terms' /
DATA carbon_monoxide(94) / '   -3.3728318564    0.0000000000   !aj, ti for [ai*tau**ti] terms' /
DATA carbon_monoxide(95) / '    3.3683460039    1.0000000000' /
DATA carbon_monoxide(96) / '   -0.0000911127   -1.5000000000' /
DATA carbon_monoxide(97) / '    1.0128000000  -23.2500376336   !aj, ti for [ai*log(1-exp(ti*tau)] terms' /
DATA carbon_monoxide(98) / '' /
DATA carbon_monoxide(99) / '' /
DATA carbon_monoxide(100) / '@EOS               !equation of state specification' /
DATA carbon_monoxide(101) / 'FEK  short Helmholtz equation of state for carbon monoxide of Lemmon and Span (2006).' /
DATA carbon_monoxide(102) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(103) / '?Lemmon, E.W. and Span, R.,' /
DATA carbon_monoxide(104) / '? "Short Fundamental Equations of State for 20 Industrial Fluids,"' /
DATA carbon_monoxide(105) / '? J. Chem. Eng. Data, 51:785-850, 2006.' /
DATA carbon_monoxide(106) / '?\' /
DATA carbon_monoxide(107) / '?The equation of state is valid from the triple point to 500 K with' /
DATA carbon_monoxide(108) / '?pressures to 100 MPa. At higher pressures, the deviations from the equation' /
DATA carbon_monoxide(109) / '?increase rapidly and it is not recommended to use the equation above 100' /
DATA carbon_monoxide(110) / '?MPa. The uncertainties in the equation are 0.3% in density (approaching 1%' /
DATA carbon_monoxide(111) / '?near the critical point), 0.2% in vapor pressure, and 2% in heat' /
DATA carbon_monoxide(112) / '?capacities.  The uncertainty in the speed of sound is unknown.' /
DATA carbon_monoxide(113) / '?\' /
DATA carbon_monoxide(114) / '!end of info section' /
DATA carbon_monoxide(115) / '68.16              !lower temperature limit [K]' /
DATA carbon_monoxide(116) / '500.0              !upper temperature limit [K]' /
DATA carbon_monoxide(117) / '100000.0           !upper pressure limit [kPa]' /
DATA carbon_monoxide(118) / '33.84              !maximum density [mol/L]' /
DATA carbon_monoxide(119) / 'PHK                                    !pointer to Cp0 model' /
DATA carbon_monoxide(120) / '28.0101                                !molecular weight [g/mol]' /
DATA carbon_monoxide(121) / '68.16                                  !triple point temperature [K]' /
DATA carbon_monoxide(122) / '15.45                                  !pressure at triple point [kPa]' /
DATA carbon_monoxide(123) / '30.33                                  !density at triple point [mol/L]' /
DATA carbon_monoxide(124) / '81.64                                  !normal boiling point temperature [K]' /
DATA carbon_monoxide(125) / '0.0497                                 !acentric factor' /
DATA carbon_monoxide(126) / '132.86        3494.0     10.85         !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA carbon_monoxide(127) / '132.86                   10.85         !reducing parameters [K, mol/L]' /
DATA carbon_monoxide(128) / '8.314472                               !gas constant [J/mol-K]' /
DATA carbon_monoxide(129) / '  12  4      0  0      0  0            !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA carbon_monoxide(130) / '  0.90554         0.25    1.0   0      !a(i),t(i),d(i),l(i)' /
DATA carbon_monoxide(131) / ' -2.4515          1.125   1.0   0' /
DATA carbon_monoxide(132) / '  0.53149         1.5     1.0   0' /
DATA carbon_monoxide(133) / '  0.024173        1.375   2.0   0' /
DATA carbon_monoxide(134) / '  0.072156        0.25    3.0   0' /
DATA carbon_monoxide(135) / '  0.00018818      0.875   7.0   0' /
DATA carbon_monoxide(136) / '  0.19405         0.625   2.0   1' /
DATA carbon_monoxide(137) / ' -0.043268        1.75    5.0   1' /
DATA carbon_monoxide(138) / ' -0.12778         3.625   1.0   2' /
DATA carbon_monoxide(139) / ' -0.027896        3.625   4.0   2' /
DATA carbon_monoxide(140) / ' -0.034154       14.5     3.0   3' /
DATA carbon_monoxide(141) / '  0.016329       12.0     4.0   3' /
DATA carbon_monoxide(142) / '' /
DATA carbon_monoxide(143) / '' /
DATA carbon_monoxide(144) / '#AUX               !auxiliary model specification' /
DATA carbon_monoxide(145) / 'PHK  Helmholtz form for the ideal-gas state for carbon monoxide of Kunz and Wagner (2004).' /
DATA carbon_monoxide(146) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(147) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA carbon_monoxide(148) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA carbon_monoxide(149) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA carbon_monoxide(150) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA carbon_monoxide(151) / '?\' /
DATA carbon_monoxide(152) / '!end of info section' /
DATA carbon_monoxide(153) / '0.                 !lower temperature limit [K]' /
DATA carbon_monoxide(154) / '1000.0             !upper temperature limit [K]' /
DATA carbon_monoxide(155) / '0.0                !upper pressure limit [kPa]' /
DATA carbon_monoxide(156) / '0.0                !maximum density [mol/L]' /
DATA carbon_monoxide(157) / '1 2  0  1 1  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA carbon_monoxide(158) / '    2.50055      1.             !ai, ti for [ai*log(tau**ti)] terms' /
DATA carbon_monoxide(159) / '   10.813340744  0.             !aj, ti for [ai*tau**ti] terms' /
DATA carbon_monoxide(160) / '  -19.834733959  1.' /
DATA carbon_monoxide(161) / '   -0.00493      5.302762306    !aj, ti for cosh and sinh terms' /
DATA carbon_monoxide(162) / '    1.02865     11.6698028' /
DATA carbon_monoxide(163) / '' /
DATA carbon_monoxide(164) / '' /
DATA carbon_monoxide(165) / '@EOS               !equation of state specification' /
DATA carbon_monoxide(166) / 'BWR  MBWR equation of state for carbon monoxide of McCarty (1989).' /
DATA carbon_monoxide(167) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(168) / '?McCarty, R.D.,' /
DATA carbon_monoxide(169) / '? "Correlations for the Thermophysical Properties of Carbon Monoxide,"' /
DATA carbon_monoxide(170) / '? National Institute of Standards and Technology, Boulder, CO, 1989.' /
DATA carbon_monoxide(171) / '?\' /
DATA carbon_monoxide(172) / '?N.B.  all temperatures on IPTS-68' /
DATA carbon_monoxide(173) / '?\' /
DATA carbon_monoxide(174) / '!end of info section' /
DATA carbon_monoxide(175) / '68.16              !lower temperature limit [K]' /
DATA carbon_monoxide(176) / '1000.0             !upper temperature limit [K]' /
DATA carbon_monoxide(177) / '30000.0            !upper pressure limit [kPa]' /
DATA carbon_monoxide(178) / '30.250             !maximum density [mol/L]' /
DATA carbon_monoxide(179) / 'CP1                                    !pointer to Cp0 model' /
DATA carbon_monoxide(180) / '28.011                                 !molecular weight [g/mol]' /
DATA carbon_monoxide(181) / '68.16                                  !triple point temperature [K]' /
DATA carbon_monoxide(182) / '15.423                                 !pressure at triple point [kPa]' /
DATA carbon_monoxide(183) / '30.249                                 !density at triple point [mol/L]' /
DATA carbon_monoxide(184) / '81.632                                 !normal boiling point temperature [K]' /
DATA carbon_monoxide(185) / '0.051                                  !acentric factor' /
DATA carbon_monoxide(186) / '132.8        3493.5       10.85        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA carbon_monoxide(187) / '132.8                     10.85        !reducing parameters [K, mol/L]' /
DATA carbon_monoxide(188) / '10.85                                  !gamma' /
DATA carbon_monoxide(189) / '0.0831434                              !gas constant [L-bar/mol-K]' /
DATA carbon_monoxide(190) / '      32       1                       !Nterm, Ncoeff per term' /
DATA carbon_monoxide(191) / '  0.8845582109949d-02    -0.2236741566840d-00     0.1742275796442d+01' /
DATA carbon_monoxide(192) / ' -0.2169146998363d+03     0.1721504267082d+04    -0.3990514770703d-04' /
DATA carbon_monoxide(193) / '  0.1036880040451d-00    -0.3376308165071d+02     0.2061895161095d+05' /
DATA carbon_monoxide(194) / '  0.2993711656350d-05     0.1856003597097d-02    -0.2114419664527d-00' /
DATA carbon_monoxide(195) / ' -0.2436986935194d-05    -0.1858029609177d-02    -0.1734563867767d+01' /
DATA carbon_monoxide(196) / '  0.1509970839260d-03    -0.2282721433205d-05     0.2202780295674d-02' /
DATA carbon_monoxide(197) / ' -0.3313357789163d-04    -0.1473412120276d+05    -0.3141136651147d+06' /
DATA carbon_monoxide(198) / ' -0.1451168999234d+03     0.6323441221817d+05    -0.2203560539926d-00' /
DATA carbon_monoxide(199) / ' -0.2087738308480d+02    -0.1508165207553d-02     0.2740740634030d+01' /
DATA carbon_monoxide(200) / '  0.8687687989627d-06    -0.1451419251928d-03    -0.3040346241285d-08' /
DATA carbon_monoxide(201) / '  0.4712050805815d-08    -0.2639772456566d-05' /
DATA carbon_monoxide(202) / '' /
DATA carbon_monoxide(203) / '' /
DATA carbon_monoxide(204) / '#AUX               !auxiliary model specification' /
DATA carbon_monoxide(205) / 'CP1  ideal gas heat capacity function' /
DATA carbon_monoxide(206) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(207) / '?McCarty, R.D.,' /
DATA carbon_monoxide(208) / '? "Correlations for the Thermophysical Properties of Carbon Monoxide,"' /
DATA carbon_monoxide(209) / '? National Institute of Standards and Technology, Boulder, CO, 1989.' /
DATA carbon_monoxide(210) / '?\' /
DATA carbon_monoxide(211) / '!end of info section' /
DATA carbon_monoxide(212) / '68.16              !lower temperature limit [K]' /
DATA carbon_monoxide(213) / '1000.0             !upper temperature limit [K]' /
DATA carbon_monoxide(214) / '0.0                !upper pressure limit [kPa]' /
DATA carbon_monoxide(215) / '0.0                !maximum density [mol/L]' /
DATA carbon_monoxide(216) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA carbon_monoxide(217) / '  7  1    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA carbon_monoxide(218) / '-0.20871594d+5   -3.00' /
DATA carbon_monoxide(219) / ' 0.89208708d+3   -2.00' /
DATA carbon_monoxide(220) / '-0.14157993d+2   -1.00' /
DATA carbon_monoxide(221) / ' 0.36028218d+1    0.00' /
DATA carbon_monoxide(222) / '-0.34021345d-3    1.00' /
DATA carbon_monoxide(223) / ' 0.44616091d-6    2.00' /
DATA carbon_monoxide(224) / '-0.15154703d-9    3.00' /
DATA carbon_monoxide(225) / ' 0.90426143d+0  30000.00' /
DATA carbon_monoxide(226) / '' /
DATA carbon_monoxide(227) / '' /
DATA carbon_monoxide(228) / '' /
DATA carbon_monoxide(229) / '#STN        !surface tension specification' /
DATA carbon_monoxide(230) / 'ST1  surface tension model of Mulero et al. (2012)' /
DATA carbon_monoxide(231) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(232) / '?Mulero, A., Cachadina, I., and Parra, M.I.' /
DATA carbon_monoxide(233) / '? "Recommended Correlations for the Surface Tension of Common Fluids,"' /
DATA carbon_monoxide(234) / '? J. Phys. Chem. Ref. Data, 41, 043105, 2012.' /
DATA carbon_monoxide(235) / '?\' /
DATA carbon_monoxide(236) / '!end of info section' /
DATA carbon_monoxide(237) / '0.0                !lower temperature limit [K]' /
DATA carbon_monoxide(238) / '132.86             !upper temperature limit [K]' /
DATA carbon_monoxide(239) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_monoxide(240) / '0.0                !(dummy) maximum density' /
DATA carbon_monoxide(241) / '1                           !number of terms in surface tension model' /
DATA carbon_monoxide(242) / '132.86                      !critical temperature used in fit (dummy)' /
DATA carbon_monoxide(243) / ' 0.02843     1.148          !sigma0 and n' /
DATA carbon_monoxide(244) / '' /
DATA carbon_monoxide(245) / '' /
DATA carbon_monoxide(246) / '#MLT        !melting line specification' /
DATA carbon_monoxide(247) / 'ML1  melting line model of Barreiros et al. (1982)' /
DATA carbon_monoxide(248) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(249) / '?Barreiros, S.F., Calado, J.C.G., Nunes da Ponte, M.' /
DATA carbon_monoxide(250) / '? "The melting curve of carbon monoxide,"' /
DATA carbon_monoxide(251) / '? J. Chem. Thermodyn., 14:1197-8, 1982.' /
DATA carbon_monoxide(252) / '?\' /
DATA carbon_monoxide(253) / '!end of info section' /
DATA carbon_monoxide(254) / '68.16              !lower temperature limit [K]' /
DATA carbon_monoxide(255) / '1000.0             !upper temperature limit [K]' /
DATA carbon_monoxide(256) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_monoxide(257) / '0.0                !(dummy) maximum density' /
DATA carbon_monoxide(258) / '1.     1000.       !reducing temperature and pressure' /
DATA carbon_monoxide(259) / '2 0 0 0 0 0                 !number of terms in melting line equation' /
DATA carbon_monoxide(260) / ' -142.941       0.          !coefficients and exponents' /
DATA carbon_monoxide(261) / ' 0.0195608      2.10747' /
DATA carbon_monoxide(262) / '' /
DATA carbon_monoxide(263) / '' /
DATA carbon_monoxide(264) / '#PS         !vapor pressure equation' /
DATA carbon_monoxide(265) / 'PS5  vapor pressure equation of Lemmon (2010).' /
DATA carbon_monoxide(266) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(267) / '?Lemmon, C.K. and Lemmon, E.W., 2010.' /
DATA carbon_monoxide(268) / '?\' /
DATA carbon_monoxide(269) / '!end of info section' /
DATA carbon_monoxide(270) / '68.16              !lower temperature limit [K]' /
DATA carbon_monoxide(271) / '132.86             !upper temperature limit [K]' /
DATA carbon_monoxide(272) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_monoxide(273) / '0.0                !(dummy) maximum density' /
DATA carbon_monoxide(274) / '132.86  3494       !reducing parameters' /
DATA carbon_monoxide(275) / '5 0 0 0 0 0        !number of terms in equation' /
DATA carbon_monoxide(276) / '-0.61192D+01       1.0' /
DATA carbon_monoxide(277) / ' 0.10411D+01       1.5' /
DATA carbon_monoxide(278) / '-0.62162D+01       3.9' /
DATA carbon_monoxide(279) / ' 0.10437D+02       4.6' /
DATA carbon_monoxide(280) / '-0.76813D+01       5.4' /
DATA carbon_monoxide(281) / '' /
DATA carbon_monoxide(282) / '' /
DATA carbon_monoxide(283) / '#DL         !saturated liquid density equation' /
DATA carbon_monoxide(284) / 'DL1  saturated liquid density equation of Lemmon (2010).' /
DATA carbon_monoxide(285) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(286) / '?Lemmon, C.K. and Lemmon, E.W., 2010.' /
DATA carbon_monoxide(287) / '?\' /
DATA carbon_monoxide(288) / '!end of info section' /
DATA carbon_monoxide(289) / '68.16              !lower temperature limit [K]' /
DATA carbon_monoxide(290) / '132.86             !upper temperature limit [K]' /
DATA carbon_monoxide(291) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_monoxide(292) / '0.0                !(dummy) maximum density' /
DATA carbon_monoxide(293) / '132.86  10.85      !reducing parameters' /
DATA carbon_monoxide(294) / '5 0 0 0 0 0        !number of terms in equation' /
DATA carbon_monoxide(295) / ' 0.29570D+01       0.398       !coefficients and exponents' /
DATA carbon_monoxide(296) / '-0.42880D+01       0.735' /
DATA carbon_monoxide(297) / ' 0.87643D+01       1.08' /
DATA carbon_monoxide(298) / '-0.84001D+01       1.5' /
DATA carbon_monoxide(299) / ' 0.36372D+01       1.9' /
DATA carbon_monoxide(300) / '' /
DATA carbon_monoxide(301) / '' /
DATA carbon_monoxide(302) / '#DV         !saturated vapor density equation' /
DATA carbon_monoxide(303) / 'DV3  saturated vapor density equation of Lemmon (2010).' /
DATA carbon_monoxide(304) / '?LITERATURE REFERENCE \' /
DATA carbon_monoxide(305) / '?Lemmon, C.K. and Lemmon, E.W., 2010.' /
DATA carbon_monoxide(306) / '?\' /
DATA carbon_monoxide(307) / '!end of info section' /
DATA carbon_monoxide(308) / '68.16              !lower temperature limit [K]' /
DATA carbon_monoxide(309) / '132.86             !upper temperature limit [K]' /
DATA carbon_monoxide(310) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_monoxide(311) / '0.0                !(dummy) maximum density' /
DATA carbon_monoxide(312) / '132.86  10.85      !reducing parameters' /
DATA carbon_monoxide(313) / '6 0 0 0 0 0        !number of terms in equation' /
DATA carbon_monoxide(314) / '-0.25439D+01       0.395         !coefficients and exponents' /
DATA carbon_monoxide(315) / '-0.55601D+01       1.21' /
DATA carbon_monoxide(316) / '-0.85276D+01       3.0' /
DATA carbon_monoxide(317) / '-0.51163D+01       3.5' /
DATA carbon_monoxide(318) / '-0.17701D+02       6.0' /
DATA carbon_monoxide(319) / '-0.29858D+02       8.0' /
DATA carbon_monoxide(320) / '' /
DATA carbon_monoxide(321) / '' /
DATA carbon_monoxide(322) / '@END' /
DATA carbon_monoxide(323) / 'c        1         2         3         4         5         6         7         8' /
DATA carbon_monoxide(324) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
! #######################################################
character(256), TARGET :: carbon_dioxide(683)
DATA carbon_dioxide(1) / 'carbon dioxide     !short name' /
DATA carbon_dioxide(2) / '124-38-9           !CAS number' /
DATA carbon_dioxide(3) / 'carbon dioxide     !full name' /
DATA carbon_dioxide(4) / 'CO2                !chemical formula' /
DATA carbon_dioxide(5) / 'R-744              !synonym' /
DATA carbon_dioxide(6) / '44.0098            !molecular weight [g/mol]' /
DATA carbon_dioxide(7) / '216.592            !triple point temperature [K]' /
DATA carbon_dioxide(8) / '194.686            !normal boiling point [K]' /
DATA carbon_dioxide(9) / '304.1282           !critical temperature [K]' /
DATA carbon_dioxide(10) / '7377.3             !critical pressure [kPa]' /
DATA carbon_dioxide(11) / '10.6249            !critical density [mol/L]' /
DATA carbon_dioxide(12) / '0.22394            !acentric factor' /
DATA carbon_dioxide(13) / '0.0                !dipole moment [Debye]' /
DATA carbon_dioxide(14) / 'IIR                !default reference state' /
DATA carbon_dioxide(15) / '9.1                !version number' /
DATA carbon_dioxide(16) / '1013               !UN Number' /
DATA carbon_dioxide(17) / 'other              !family' /
DATA carbon_dioxide(18) / '0.0                !heating value (gross or superior) [kJ/mol]' /
DATA carbon_dioxide(19) / '1.                 !GWP (IPCC 2007)' /
DATA carbon_dioxide(20) / '40000.             !RCL (ppm v/v, ASHRAE Standard 34, 2010)' /
DATA carbon_dioxide(21) / 'A1                 !Safety Group (ASHRAE Standard 34, 2010)' /
DATA carbon_dioxide(22) / '' /
DATA carbon_dioxide(23) / '' /
DATA carbon_dioxide(24) / '#EOS               !equation of state specification' /
DATA carbon_dioxide(25) / 'FEQ  Helmholtz equation of state for carbon dioxide of Span and Wagner (1996).' /
DATA carbon_dioxide(26) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(27) / '?Span, R. and Wagner, W.,' /
DATA carbon_dioxide(28) / '? "A New Equation of State for Carbon Dioxide Covering the Fluid Region' /
DATA carbon_dioxide(29) / '? from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa,"' /
DATA carbon_dioxide(30) / '? J. Phys. Chem. Ref. Data, 25(6):1509-1596, 1996.' /
DATA carbon_dioxide(31) / '?\' /
DATA carbon_dioxide(32) / '?At pressures up to 30 MPa and temperatures up to 523 K, the estimated' /
DATA carbon_dioxide(33) / '?uncertainty ranges from 0.03% to 0.05% in density, 0.03% (in the vapor)' /
DATA carbon_dioxide(34) / '?to 1% in the speed of sound (0.5% in the liquid) and 0.15% (in the' /
DATA carbon_dioxide(35) / '?vapor) to 1.5% (in the liquid) in heat capacity.  Special interest has' /
DATA carbon_dioxide(36) / '?been focused on the description of the critical region and the' /
DATA carbon_dioxide(37) / '?extrapolation behavior of the formulation (to the limits of chemical' /
DATA carbon_dioxide(38) / '?stability).' /
DATA carbon_dioxide(39) / '?\' /
DATA carbon_dioxide(40) / '!end of info section' /
DATA carbon_dioxide(41) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(42) / '2000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(43) / '800000.0           !upper pressure limit [kPa]' /
DATA carbon_dioxide(44) / '37.24              !maximum density [mol/L]' /
DATA carbon_dioxide(45) / 'CPP                                    !pointer to Cp0 model' /
DATA carbon_dioxide(46) / '44.0098                                !molecular weight [g/mol]' /
DATA carbon_dioxide(47) / '216.592                                !triple point temperature [K]' /
DATA carbon_dioxide(48) / '517.95                                 !pressure at triple point [kPa]' /
DATA carbon_dioxide(49) / '26.777                                 !density at triple point [mol/L]' /
DATA carbon_dioxide(50) / '194.686                                !normal boiling point temperature [K]' /
DATA carbon_dioxide(51) / '0.22394                                !acentric factor' /
DATA carbon_dioxide(52) / '304.1282     7377.3       10.6249063   !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA carbon_dioxide(53) / '304.1282                  10.6249063   !reducing parameters [K, mol/L]' /
DATA carbon_dioxide(54) / '8.31451                                !gas constant [J/mol-K]' /
DATA carbon_dioxide(55) / '      34  4      8  12      0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA carbon_dioxide(56) / ' 0.388568232032d+00  0.000   1.00    0 !a(i),t(i),d(i),l(i)' /
DATA carbon_dioxide(57) / ' 0.293854759427d+01  0.750   1.00    0' /
DATA carbon_dioxide(58) / '-0.558671885349d+01  1.000   1.00    0' /
DATA carbon_dioxide(59) / '-0.767531995925d+00  2.000   1.00    0' /
DATA carbon_dioxide(60) / ' 0.317290055804d+00  0.750   2.00    0' /
DATA carbon_dioxide(61) / ' 0.548033158978d+00  2.000   2.00    0' /
DATA carbon_dioxide(62) / ' 0.122794112203d+00  0.750   3.00    0' /
DATA carbon_dioxide(63) / ' 0.216589615432d+01  1.500   1.00    1' /
DATA carbon_dioxide(64) / ' 0.158417351097d+01  1.500   2.00    1' /
DATA carbon_dioxide(65) / '-0.231327054055d+00  2.500   4.00    1' /
DATA carbon_dioxide(66) / ' 0.581169164314d-01  0.000   5.00    1' /
DATA carbon_dioxide(67) / '-0.553691372054d+00  1.500   5.00    1' /
DATA carbon_dioxide(68) / ' 0.489466159094d+00  2.000   5.00    1' /
DATA carbon_dioxide(69) / '-0.242757398435d-01  0.000   6.00    1' /
DATA carbon_dioxide(70) / ' 0.624947905017d-01  1.000   6.00    1' /
DATA carbon_dioxide(71) / '-0.121758602252d+00  2.000   6.00    1' /
DATA carbon_dioxide(72) / '-0.370556852701d+00  3.000   1.00    2' /
DATA carbon_dioxide(73) / '-0.167758797004d-01  6.000   1.00    2' /
DATA carbon_dioxide(74) / '-0.119607366380d+00  3.000   4.00    2' /
DATA carbon_dioxide(75) / '-0.456193625088d-01  6.000   4.00    2' /
DATA carbon_dioxide(76) / ' 0.356127892703d-01  8.000   4.00    2' /
DATA carbon_dioxide(77) / '-0.744277271321d-02  6.000   7.00    2' /
DATA carbon_dioxide(78) / '-0.173957049024d-02  0.000   8.00    2' /
DATA carbon_dioxide(79) / '-0.218101212895d-01  7.000   2.00    3' /
DATA carbon_dioxide(80) / ' 0.243321665592d-01 12.000   3.00    3' /
DATA carbon_dioxide(81) / '-0.374401334235d-01 16.000   3.00    3' /
DATA carbon_dioxide(82) / ' 0.143387157569d+00 22.000   5.00    4' /
DATA carbon_dioxide(83) / '-0.134919690833d+00 24.000   5.00    4' /
DATA carbon_dioxide(84) / '-0.231512250535d-01 16.000   6.00    4' /
DATA carbon_dioxide(85) / ' 0.123631254929d-01 24.000   7.00    4' /
DATA carbon_dioxide(86) / ' 0.210583219729d-02  8.000   8.00    4' /
DATA carbon_dioxide(87) / '-0.339585190264d-03  2.000  10.00    4' /
DATA carbon_dioxide(88) / ' 0.559936517716d-02 28.000   4.00    5' /
DATA carbon_dioxide(89) / '-0.303351180556d-03 14.000   8.00    6' /
DATA carbon_dioxide(90) / '-0.213654886883d+03  1.000   2.00    2 2  -25.  -325.  1.16 1.   0.   0.  0.' /
DATA carbon_dioxide(91) / ' 0.266415691493d+05  0.000   2.00    2 2  -25.  -300.  1.19 1.   0.   0.  0.' /
DATA carbon_dioxide(92) / '-0.240272122046d+05  1.000   2.00    2 2  -25.  -300.  1.19 1.   0.   0.  0.' /
DATA carbon_dioxide(93) / '-0.283416034240d+03  3.000   3.00    2 2  -15.  -275.  1.25 1.   0.   0.  0.' /
DATA carbon_dioxide(94) / ' 0.212472844002d+03  3.000   3.00    2 2  -20.  -275.  1.22 1.   0.   0.  0.' /
DATA carbon_dioxide(95) / '-0.666422765408d+00  0.000   1.00    2 2  0.875  0.300 0.70 10.0 275. 0.3 3.5' /
DATA carbon_dioxide(96) / ' 0.726086323499d+00  0.000   1.00    2 2  0.925  0.300 0.70 10.0 275. 0.3 3.5' /
DATA carbon_dioxide(97) / ' 0.550686686128d-01  0.000   1.00    2 2  0.875  0.300 0.70 12.5 275. 1.0 3.' /
DATA carbon_dioxide(98) / '' /
DATA carbon_dioxide(99) / '' /
DATA carbon_dioxide(100) / '#AUX               !auxiliary model specification' /
DATA carbon_dioxide(101) / 'CPP  ideal gas heat capacity function' /
DATA carbon_dioxide(102) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(103) / '?Span, R. and Wagner, W.,' /
DATA carbon_dioxide(104) / '? "A New Equation of State for Carbon Dioxide Covering the Fluid Region' /
DATA carbon_dioxide(105) / '? from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa,"' /
DATA carbon_dioxide(106) / '? J. Phys. Chem. Ref. Data, 25(6):1509-1596, 1996.' /
DATA carbon_dioxide(107) / '?\' /
DATA carbon_dioxide(108) / '!end of info section' /
DATA carbon_dioxide(109) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(110) / '2000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(111) / '0.0                !upper pressure limit [kPa]' /
DATA carbon_dioxide(112) / '0.0                !maximum density [mol/L]' /
DATA carbon_dioxide(113) / '1.0          8.31451                   !reducing parameters for T, Cp0' /
DATA carbon_dioxide(114) / '  1  5    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA carbon_dioxide(115) / ' 0.35000000d+01    0.00' /
DATA carbon_dioxide(116) / '1.99427042         958.49956' /
DATA carbon_dioxide(117) / '0.621052475       1858.80115' /
DATA carbon_dioxide(118) / '0.411952928       2061.10114' /
DATA carbon_dioxide(119) / '1.04028922        3443.89908' /
DATA carbon_dioxide(120) / '0.0832767753      8238.20035' /
DATA carbon_dioxide(121) / '' /
DATA carbon_dioxide(122) / '' /
DATA carbon_dioxide(123) / '@EOS               !equation of state specification' /
DATA carbon_dioxide(124) / 'FEK  Helmholtz equation of state for carbon dioxide of Kunz and Wagner (2004).' /
DATA carbon_dioxide(125) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(126) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA carbon_dioxide(127) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA carbon_dioxide(128) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA carbon_dioxide(129) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA carbon_dioxide(130) / '?\' /
DATA carbon_dioxide(131) / '!end of info section' /
DATA carbon_dioxide(132) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(133) / '1100.0             !upper temperature limit [K]' /
DATA carbon_dioxide(134) / '800000.0           !upper pressure limit [kPa]' /
DATA carbon_dioxide(135) / '37.24              !maximum density [mol/L]' /
DATA carbon_dioxide(136) / 'PHK                                    !pointer to Cp0 model' /
DATA carbon_dioxide(137) / '44.0095                                !molecular weight [g/mol]' /
DATA carbon_dioxide(138) / '216.592                                !triple point temperature [K]' /
DATA carbon_dioxide(139) / '517.94                                 !pressure at triple point [kPa]' /
DATA carbon_dioxide(140) / '26.78                                  !density at triple point [mol/L]' /
DATA carbon_dioxide(141) / '185.36                                 !normal boiling point temperature [K]' /
DATA carbon_dioxide(142) / ' 0.225                                 !acentric factor' /
DATA carbon_dioxide(143) / '304.1282     7377.3      10.624978698  !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA carbon_dioxide(144) / '304.1282                 10.624978698  !reducing parameters [K, mol/L]' /
DATA carbon_dioxide(145) / '8.314472                               !gas constant [J/mol-K]' /
DATA carbon_dioxide(146) / '  22  4      0  0      0  0            !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA carbon_dioxide(147) / ' 0.52646564804653       0.00   1.  0' /
DATA carbon_dioxide(148) / '-0.14995725042592d1     1.25   1.  0' /
DATA carbon_dioxide(149) / ' 0.27329786733782       1.625  2.  0' /
DATA carbon_dioxide(150) / ' 0.12949500022786       0.375  3.  0' /
DATA carbon_dioxide(151) / ' 0.15404088341841       0.375  3.  1' /
DATA carbon_dioxide(152) / '-0.58186950946814       1.375  3.  1' /
DATA carbon_dioxide(153) / '-0.18022494838296       1.125  4.  1' /
DATA carbon_dioxide(154) / '-0.95389904072812d-1    1.375  5.  1' /
DATA carbon_dioxide(155) / '-0.80486819317679d-2    0.125  6.  1' /
DATA carbon_dioxide(156) / '-0.35547751273090d-1    1.625  6.  1' /
DATA carbon_dioxide(157) / '-0.28079014882405       3.75   1.  2' /
DATA carbon_dioxide(158) / '-0.82435890081677d-1    3.5    4.  2' /
DATA carbon_dioxide(159) / ' 0.10832427979006d-1    7.5    1.  3' /
DATA carbon_dioxide(160) / '-0.67073993161097d-2    8.0    1.  3' /
DATA carbon_dioxide(161) / '-0.46827907600524d-2    6.0    3.  3' /
DATA carbon_dioxide(162) / '-0.28359911832177d-1   16.0    3.  3' /
DATA carbon_dioxide(163) / ' 0.19500174744098d-1   11.0    4.  3' /
DATA carbon_dioxide(164) / '-0.21609137507166      24.0    5.  5' /
DATA carbon_dioxide(165) / ' 0.43772794926972      26.0    5.  5' /
DATA carbon_dioxide(166) / '-0.22130790113593      28.0    5.  5' /
DATA carbon_dioxide(167) / ' 0.15190189957331d-1   24.0    5.  6' /
DATA carbon_dioxide(168) / '-0.15380948953300d-1   26.0    5.  6' /
DATA carbon_dioxide(169) / '' /
DATA carbon_dioxide(170) / '' /
DATA carbon_dioxide(171) / '#AUX               !auxiliary model specification' /
DATA carbon_dioxide(172) / 'PHK  Helmholtz form for the ideal-gas state for carbon dioxide of Kunz and Wagner (2004).' /
DATA carbon_dioxide(173) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(174) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA carbon_dioxide(175) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA carbon_dioxide(176) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA carbon_dioxide(177) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA carbon_dioxide(178) / '?\' /
DATA carbon_dioxide(179) / '!end of info section' /
DATA carbon_dioxide(180) / '0.                 !lower temperature limit [K]' /
DATA carbon_dioxide(181) / '1000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(182) / '0.0                !upper pressure limit [kPa]' /
DATA carbon_dioxide(183) / '0.0                !maximum density [mol/L]' /
DATA carbon_dioxide(184) / '1 2  0  2 2  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA carbon_dioxide(185) / '    2.50002      1.             !ai, ti for [ai*log(tau**ti)] terms' /
DATA carbon_dioxide(186) / '   11.925152758  0.             !aj, ti for [ai*tau**ti] terms' /
DATA carbon_dioxide(187) / '  -16.118762264  1.' /
DATA carbon_dioxide(188) / '    1.06044     -2.844425476    !aj, ti for cosh and sinh terms' /
DATA carbon_dioxide(189) / '   -0.01393      1.12159609' /
DATA carbon_dioxide(190) / '    2.04452      3.022758166' /
DATA carbon_dioxide(191) / '    2.03366      1.589964364' /
DATA carbon_dioxide(192) / '' /
DATA carbon_dioxide(193) / '' /
DATA carbon_dioxide(194) / '@EOS               !equation of state specification' /
DATA carbon_dioxide(195) / 'BWR  MBWR equation of state for carbon dioxide of Ely et al. (1987).' /
DATA carbon_dioxide(196) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(197) / '?Ely, J.F., Magee, J.W., and Haynes, W.M.,' /
DATA carbon_dioxide(198) / '? "Thermophysical properties for special high CO2 content mixtures,"' /
DATA carbon_dioxide(199) / '? Research Report RR-110, Gas Processors Association, Tulsa, OK, 1987.' /
DATA carbon_dioxide(200) / '?\' /
DATA carbon_dioxide(201) / '?Note:  This report contains both MBWR and FEQ (referred to as the Schmidt-Wagner' /
DATA carbon_dioxide(202) / '? equation of state in the report) equations.  The FEQ (Schmidt-Wagner) will' /
DATA carbon_dioxide(203) / '? give slightly better numbers very close to the critical point but for most' /
DATA carbon_dioxide(204) / '? calculations, the MBWR is the recommended equation.' /
DATA carbon_dioxide(205) / '?\' /
DATA carbon_dioxide(206) / '!end of info section' /
DATA carbon_dioxide(207) / '216.58             !lower temperature limit [K]' /
DATA carbon_dioxide(208) / '440.1              !upper temperature limit [K]' /
DATA carbon_dioxide(209) / '40000.0            !upper pressure limit [kPa]' /
DATA carbon_dioxide(210) / '27.778             !maximum density [mol/L]' /
DATA carbon_dioxide(211) / 'CP1                                    !pointer to Cp0 model' /
DATA carbon_dioxide(212) / '44.0098                                !molecular weight [g/mol]' /
DATA carbon_dioxide(213) / '216.58                                 !triple point temperature [K]' /
DATA carbon_dioxide(214) / '518.2                                  !pressure at triple point [kPa]' /
DATA carbon_dioxide(215) / '26.778                                 !density at triple point [mol/L]' /
DATA carbon_dioxide(216) / '194.75                                 !normal boiling point temperature [K]' /
DATA carbon_dioxide(217) / '0.22394                                !acentric factor' /
DATA carbon_dioxide(218) / '304.21       7384.325     10.60        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA carbon_dioxide(219) / '304.21                    10.60        !reducing parameters [K, mol/L]' /
DATA carbon_dioxide(220) / '10.60                                  !gamma' /
DATA carbon_dioxide(221) / '0.0831434                              !gas constant [L-bar/mol-K]' /
DATA carbon_dioxide(222) / '      32       1                       !Nterm, Ncoeff per term' /
DATA carbon_dioxide(223) / '  -0.981851065838d-02   0.995062267309d+00  -0.228380160313d+02' /
DATA carbon_dioxide(224) / '   0.281827634529d+04  -0.347001262699d+06   0.394706709102d-03' /
DATA carbon_dioxide(225) / '  -0.325550000110d+00   0.484320083063d+01  -0.352181542995d+06' /
DATA carbon_dioxide(226) / '  -0.324053603343d-04   0.468596684665d-01  -0.754547012075d+01' /
DATA carbon_dioxide(227) / '  -0.381894354016d-04  -0.442192933859d-01   0.516925168095d+02' /
DATA carbon_dioxide(228) / '   0.212450985237d-02  -0.261009474785d-04  -0.888533388977d-01' /
DATA carbon_dioxide(229) / '   0.155226179403d-02   0.415091004940d+06  -0.110173967489d+08' /
DATA carbon_dioxide(230) / '   0.291990583344d+04   0.143254606508d+08   0.108574207533d+02' /
DATA carbon_dioxide(231) / '  -0.247799657039d+03   0.199293590763d-01   0.102749908059d+03' /
DATA carbon_dioxide(232) / '   0.377618865158d-04  -0.332276512346d-02   0.179196707121d-07' /
DATA carbon_dioxide(233) / '   0.945076627807d-05  -0.123400943061d-02' /
DATA carbon_dioxide(234) / '' /
DATA carbon_dioxide(235) / '' /
DATA carbon_dioxide(236) / '@EOS               !equation of state specification' /
DATA carbon_dioxide(237) / 'FE1  Helmholtz equation of state for carbon dioxide of Ely et al. (1987).' /
DATA carbon_dioxide(238) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(239) / '?Ely, J.F., Magee, J.W., and Haynes, W.M.,' /
DATA carbon_dioxide(240) / '? "Thermophysical properties for special high CO2 content mixtures,"' /
DATA carbon_dioxide(241) / '? Research Report RR-110, Gas Processors Association, Tulsa, OK, 1987.' /
DATA carbon_dioxide(242) / '?\' /
DATA carbon_dioxide(243) / '?Note:  This report contains both MBWR and FEQ (referred to as the Schmidt-Wagner' /
DATA carbon_dioxide(244) / '? equation of state in the report) equations.  The FEQ (Schmidt-Wagner) will' /
DATA carbon_dioxide(245) / '? give slightly better numbers very close to the critical point but for most' /
DATA carbon_dioxide(246) / '? calculations, the MBWR is the recommended equation.' /
DATA carbon_dioxide(247) / '?\' /
DATA carbon_dioxide(248) / '!end of info section' /
DATA carbon_dioxide(249) / '216.58             !lower temperature limit [K]' /
DATA carbon_dioxide(250) / '1000.              !upper temperature limit [K]' /
DATA carbon_dioxide(251) / '100000.0           !upper pressure limit [kPa]' /
DATA carbon_dioxide(252) / '26.776             !maximum density [mol/L]' /
DATA carbon_dioxide(253) / 'CP1                                    !pointer to Cp0 model' /
DATA carbon_dioxide(254) / '44.0098                                !molecular weight [g/mol]' /
DATA carbon_dioxide(255) / '216.58                                 !triple point temperature [K]' /
DATA carbon_dioxide(256) / '518.03                                 !pressure at triple point [kPa]' /
DATA carbon_dioxide(257) / '26.776                                 !density at triple point [mol/L]' /
DATA carbon_dioxide(258) / '194.75                                 !normal boiling point temperature [K]' /
DATA carbon_dioxide(259) / '0.22394                                !acentric factor' /
DATA carbon_dioxide(260) / '304.13       7375.21      10.63        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA carbon_dioxide(261) / '304.13                    10.63        !reducing parameters [K, mol/L]' /
DATA carbon_dioxide(262) / '8.31434                                !gas constant [J/mol-K]' /
DATA carbon_dioxide(263) / '      32  4      0  0       0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA carbon_dioxide(264) / '  0.485497428986d+00   0.000   1.00    0  !a(i),t(i),d(i),l(i)' /
DATA carbon_dioxide(265) / ' -0.191900462349d+01   1.500   1.00    0' /
DATA carbon_dioxide(266) / '  0.451739876847d+00   2.500   1.00    0' /
DATA carbon_dioxide(267) / '  0.838475229022d-02  -0.500   2.00    0' /
DATA carbon_dioxide(268) / '  0.310719428397d+00   1.500   2.00    0' /
DATA carbon_dioxide(269) / ' -0.183619563850d+00   2.000   2.00    0' /
DATA carbon_dioxide(270) / '  0.448878785519d-01   0.000   3.00    0' /
DATA carbon_dioxide(271) / ' -0.362211893044d-01   1.000   3.00    0' /
DATA carbon_dioxide(272) / ' -0.169827491865d-01   2.500   3.00    0' /
DATA carbon_dioxide(273) / '  0.803504394396d-03   0.000   6.00    0' /
DATA carbon_dioxide(274) / '  0.320223641512d-03   2.000   7.00    0' /
DATA carbon_dioxide(275) / ' -0.658956249553d-05   5.000   7.00    0' /
DATA carbon_dioxide(276) / ' -0.461991678692d-04   2.000   8.00    0' /
DATA carbon_dioxide(277) / ' -0.385989029443d+00   5.000   1.00    2' /
DATA carbon_dioxide(278) / '  0.131878614095d+00   6.000   1.00    2' /
DATA carbon_dioxide(279) / '  0.109639470331d+00   3.500   2.00    2' /
DATA carbon_dioxide(280) / ' -0.310044422115d-01   5.500   2.00    2' /
DATA carbon_dioxide(281) / ' -0.989797992915d-01   3.000   3.00    2' /
DATA carbon_dioxide(282) / ' -0.222934996927d-01   7.000   3.00    2' /
DATA carbon_dioxide(283) / ' -0.225488505376d-01   6.000   5.00    2' /
DATA carbon_dioxide(284) / ' -0.595661202393d-02   8.500   6.00    2' /
DATA carbon_dioxide(285) / ' -0.219959964099d-01   4.000   7.00    2' /
DATA carbon_dioxide(286) / '  0.140330955537d-01   6.500   8.00    2' /
DATA carbon_dioxide(287) / ' -0.315424157971d-02   5.500  10.00    2' /
DATA carbon_dioxide(288) / '  0.443394060420d-03  22.000   2.00    4' /
DATA carbon_dioxide(289) / ' -0.487628903103d-02  11.000   3.00    4' /
DATA carbon_dioxide(290) / ' -0.311643343682d-01  18.000   3.00    4' /
DATA carbon_dioxide(291) / '  0.226083669848d-01  11.000   4.00    4' /
DATA carbon_dioxide(292) / '  0.186651858191d-01  23.000   4.00    4' /
DATA carbon_dioxide(293) / ' -0.399277963883d+00  17.000   5.00    4' /
DATA carbon_dioxide(294) / '  0.464945130861d+00  18.000   5.00    4' /
DATA carbon_dioxide(295) / ' -0.817090055061d-01  23.000   5.00    4' /
DATA carbon_dioxide(296) / '' /
DATA carbon_dioxide(297) / '' /
DATA carbon_dioxide(298) / '#AUX               !auxiliary model specification' /
DATA carbon_dioxide(299) / 'CP1  ideal gas heat capacity function of Ely et al.' /
DATA carbon_dioxide(300) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(301) / '?Ely, J.F., Magee, J.W., and Haynes, W.M.,' /
DATA carbon_dioxide(302) / '? "Thermophysical properties for special high CO2 content mixtures,"' /
DATA carbon_dioxide(303) / '? Research Report RR-110, Gas Processors Association, Tulsa, OK, 1987.' /
DATA carbon_dioxide(304) / '?\' /
DATA carbon_dioxide(305) / '!end of info section' /
DATA carbon_dioxide(306) / '200.0              !lower temperature limit [K]' /
DATA carbon_dioxide(307) / '1000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(308) / '0.0                !upper pressure limit [kPa]' /
DATA carbon_dioxide(309) / '0.0                !maximum density [mol/L]' /
DATA carbon_dioxide(310) / '1.0          8.31441                   !reducing parameters for T, Cp0' /
DATA carbon_dioxide(311) / '  1  3    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA carbon_dioxide(312) / '3.50d0           0.00                  !c(i), power of T' /
DATA carbon_dioxide(313) / ' 2.00d0        960.11d0                !=omega_1 (degenerate mode--taken twice)' /
DATA carbon_dioxide(314) / ' 1.00d0       1932.00d0                !=omega_2' /
DATA carbon_dioxide(315) / ' 1.00d0       3380.20d0                !=omega_3' /
DATA carbon_dioxide(316) / '' /
DATA carbon_dioxide(317) / '' /
DATA carbon_dioxide(318) / '#AUX               !auxiliary model specification' /
DATA carbon_dioxide(319) / 'CP2  ideal gas heat capacity function of McCarty' /
DATA carbon_dioxide(320) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(321) / '?McCarty, R.D.,' /
DATA carbon_dioxide(322) / '? "Correlations for the Thermophysical Properties of Carbon Dioxide,"' /
DATA carbon_dioxide(323) / '? Unpublished correlation, National Institute of Standards' /
DATA carbon_dioxide(324) / '? and Technology, Boulder, 1988.' /
DATA carbon_dioxide(325) / '?\' /
DATA carbon_dioxide(326) / '?Use of this Cp0 equation in conjunction with Elys BWR will produce numbers' /
DATA carbon_dioxide(327) / '? identical to those calculated in NIST12, Version 3.0.' /
DATA carbon_dioxide(328) / '?\' /
DATA carbon_dioxide(329) / '!end of info section' /
DATA carbon_dioxide(330) / '216.58             !lower temperature limit [K]' /
DATA carbon_dioxide(331) / '1000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(332) / '0.0                !upper pressure limit [kPa]' /
DATA carbon_dioxide(333) / '0.0                !maximum density [mol/L]' /
DATA carbon_dioxide(334) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA carbon_dioxide(335) / '  7  1    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA carbon_dioxide(336) / ' -0.9379061144997d+07   -3.00d0' /
DATA carbon_dioxide(337) / '  0.2028666045159d+06   -2.00d0' /
DATA carbon_dioxide(338) / ' -0.1595185479613d+04   -1.00d0' /
DATA carbon_dioxide(339) / '  0.8189952737742d+01    0.00d0' /
DATA carbon_dioxide(340) / ' -0.1298281615271d-02    1.00d0' /
DATA carbon_dioxide(341) / '  0.1037209193687d-05    2.00d0' /
DATA carbon_dioxide(342) / ' -0.3399620971158d-09    3.00d0' /
DATA carbon_dioxide(343) / '  0.6961565991385d+00 30000.d0' /
DATA carbon_dioxide(344) / '' /
DATA carbon_dioxide(345) / '' /
DATA carbon_dioxide(346) / '@EOS               !equation of state specification' /
DATA carbon_dioxide(347) / 'FES  short Helmholtz equation of state for carbon dioxide of Span and Wagner (2003).' /
DATA carbon_dioxide(348) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(349) / '?Span, R. and Wagner, W.' /
DATA carbon_dioxide(350) / '? "Equations of State for Technical Applications. III. Results for Polar Fluids,"' /
DATA carbon_dioxide(351) / '? Int. J. Thermophys., 24(1):111-162, 2003.' /
DATA carbon_dioxide(352) / '?\' /
DATA carbon_dioxide(353) / '?The uncertainties of the equation of state are approximately 0.2% (to' /
DATA carbon_dioxide(354) / '?0.5% at high pressures) in density, 1% (in the vapor phase) to 2% in' /
DATA carbon_dioxide(355) / '?heat capacity, 1% (in the vapor phase) to 2% in the speed of sound, and' /
DATA carbon_dioxide(356) / '?0.2% in vapor pressure, except in the critical region.' /
DATA carbon_dioxide(357) / '?\' /
DATA carbon_dioxide(358) / '!end of info section' /
DATA carbon_dioxide(359) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(360) / '600.0              !upper temperature limit [K]' /
DATA carbon_dioxide(361) / '100000.0           !upper pressure limit [kPa]' /
DATA carbon_dioxide(362) / '37.24              !maximum density [mol/L]' /
DATA carbon_dioxide(363) / 'CPP                                    !pointer to Cp0 model' /
DATA carbon_dioxide(364) / '44.01                                  !molecular weight [g/mol]' /
DATA carbon_dioxide(365) / '216.592                                !triple point temperature [K]' /
DATA carbon_dioxide(366) / '517.86                                 !pressure at triple point [kPa]' /
DATA carbon_dioxide(367) / '26.795                                 !density at triple point [mol/L]' /
DATA carbon_dioxide(368) / '185.3                                  !normal boiling point temperature [K]' /
DATA carbon_dioxide(369) / '0.225                                  !acentric factor' /
DATA carbon_dioxide(370) / '304.1282     7377.3       10.624858    !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA carbon_dioxide(371) / '304.1282                  10.624858    !reducing parameters [K, mol/L]' /
DATA carbon_dioxide(372) / '8.31451                                !gas constant [J/mol-K]' /
DATA carbon_dioxide(373) / '      12  4      0  0       0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA carbon_dioxide(374) / ' 0.898751080000E+00  0.25    1.0     0 !a(i),t(i),d(i),l(i)' /
DATA carbon_dioxide(375) / '-0.212819850000E+01  1.25    1.0     0' /
DATA carbon_dioxide(376) / '-0.681903200000E-01  1.5     1.0     0' /
DATA carbon_dioxide(377) / ' 0.763553060000E-01  0.25    3.0     0' /
DATA carbon_dioxide(378) / ' 0.220532530000E-03  0.875   7.0     0' /
DATA carbon_dioxide(379) / ' 0.415418230000E+00  2.375   1.0     1' /
DATA carbon_dioxide(380) / ' 0.713356570000E+00  2.0     2.0     1' /
DATA carbon_dioxide(381) / ' 0.303542340000E-03  2.125   5.0     1' /
DATA carbon_dioxide(382) / '-0.366431430000E+00  3.5     1.0     2' /
DATA carbon_dioxide(383) / '-0.144077810000E-02  6.5     1.0     2' /
DATA carbon_dioxide(384) / '-0.891667070000E-01  4.75    4.0     2' /
DATA carbon_dioxide(385) / '-0.236998870000E-01 12.5     2.0     3' /
DATA carbon_dioxide(386) / '' /
DATA carbon_dioxide(387) / '' /
DATA carbon_dioxide(388) / '#TCX               !thermal conductivity model specification' /
DATA carbon_dioxide(389) / 'TC1  pure fluid thermal conductivity model of Vesovic et al. (1990).' /
DATA carbon_dioxide(390) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(391) / '?Vesovic, V., Wakeham, W.A., Olchowy, G.A., Sengers, J.V., Watson, J.T.R.,' /
DATA carbon_dioxide(392) / '? and Millat, J.,' /
DATA carbon_dioxide(393) / '? "The transport properties of carbon dioxide,"' /
DATA carbon_dioxide(394) / '? J. Phys. Chem. Ref. Data, 19:763-808, 1990.' /
DATA carbon_dioxide(395) / '?\' /
DATA carbon_dioxide(396) / '?Note:  Vesovic et al. use a crossover equation of state to compute derivatives' /
DATA carbon_dioxide(397) / '? in the critical region; the default EOS is used here.  Also, their' /
DATA carbon_dioxide(398) / '? "simplified" critical enhancement for thermal conductivity is used.' /
DATA carbon_dioxide(399) / '?\' /
DATA carbon_dioxide(400) / '?The uncertainty in thermal conductivity is less than 5%.' /
DATA carbon_dioxide(401) / '?\' /
DATA carbon_dioxide(402) / '!end of info section' /
DATA carbon_dioxide(403) / '216.58             !lower temperature limit [K]' /
DATA carbon_dioxide(404) / '2000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(405) / '800000.0           !upper pressure limit [kPa]' /
DATA carbon_dioxide(406) / '37.24              !maximum density [mol/L]' /
DATA carbon_dioxide(407) / '2   6              !# terms for dilute gas function:  numerator, denominator' /
DATA carbon_dioxide(408) / '251.196     0.001        !reducing parameters for T (=eps/k), tcx (orig in mW/m-K)' /
DATA carbon_dioxide(409) / ' 7.5378307d+0   0.50d0   !coeff (=0.475598*SQRT(eps/k)), power in T (T* in this case)' /
DATA carbon_dioxide(410) / ' 4.8109652d-2 -99.00d0   !power -99 indicates: mult above numerator term by [1 + coeff*(Cp0 - 2.5*R)], where coeff = 0.4/R' /
DATA carbon_dioxide(411) / ' 0.4226159d+0   0.00d0   !denominator is Eq 30 in Vesovic' /
DATA carbon_dioxide(412) / ' 0.6280115d+0  -1.00d0' /
DATA carbon_dioxide(413) / '-0.5387661d+0  -2.00d0' /
DATA carbon_dioxide(414) / ' 0.6735941d+0  -3.00d0' /
DATA carbon_dioxide(415) / '-0.4362677d+0  -6.00d0' /
DATA carbon_dioxide(416) / ' 0.2255388d+0  -7.00d0' /
DATA carbon_dioxide(417) / '4   0              !# terms for background gas function:  numerator, denominator' /
DATA carbon_dioxide(418) / '1.0     2.272221d-2   1.0d-3              !reducing par for T, rho, tcx (orig corr in kg/m**3, mW/m-K)' /
DATA carbon_dioxide(419) / ' 2.447164d-2    0.0  1.0  0.0 !coeff, powers of T, rho, spare for future use' /
DATA carbon_dioxide(420) / ' 8.705605d-05   0.00d0   2.00d0   0.00d0' /
DATA carbon_dioxide(421) / '-6.547950d-08   0.00d0   3.00d0   0.00d0' /
DATA carbon_dioxide(422) / ' 6.594919d-11   0.00d0   4.00d0   0.00d0' /
DATA carbon_dioxide(423) / 'TK3                !pointer to critical enhancement auxiliary function' /
DATA carbon_dioxide(424) / '' /
DATA carbon_dioxide(425) / '' /
DATA carbon_dioxide(426) / '#AUX               !thermal conductivity critical enhancement model' /
DATA carbon_dioxide(427) / 'TK3  simplified thermal conductivity critical enhancement of Olchowy and Sengers' /
DATA carbon_dioxide(428) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(429) / '?Olchowy, G.A. and Sengers, J.V.,' /
DATA carbon_dioxide(430) / '? "A simplified representation for the thermal conductivity of fluids in the' /
DATA carbon_dioxide(431) / '? critical region,"' /
DATA carbon_dioxide(432) / '? Int. J. Thermophysics, 10:417-426, 1989.' /
DATA carbon_dioxide(433) / '?\' /
DATA carbon_dioxide(434) / '?as applied to CO2 by:' /
DATA carbon_dioxide(435) / '?\' /
DATA carbon_dioxide(436) / '?Vesovic, V., Wakeham, W.A., Olchowy, G.A., Sengers, J.V., Watson, J.T.R.,' /
DATA carbon_dioxide(437) / '? and Millat, J.,' /
DATA carbon_dioxide(438) / '? "The transport properties of carbon dioxide,"' /
DATA carbon_dioxide(439) / '? J. Phys. Chem. Ref. Data, 19:763-808, 1990.' /
DATA carbon_dioxide(440) / '?\' /
DATA carbon_dioxide(441) / '!end of info section' /
DATA carbon_dioxide(442) / '216.58             !lower temperature limit [K]' /
DATA carbon_dioxide(443) / '2000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(444) / '800000.0           !upper pressure limit [kPa]' /
DATA carbon_dioxide(445) / '37.24              !maximum density [mol/L]' /
DATA carbon_dioxide(446) / '9  0  0  0         !# terms:  CO2-terms, spare, spare, spare' /
DATA carbon_dioxide(447) / '1.0     1.0     1.0      !reducing par for T, rho, tcx (mW/m-K)' /
DATA carbon_dioxide(448) / ' 0.630d+00         !gnu (universal exponent)' /
DATA carbon_dioxide(449) / ' 1.2415d+00        !gamma (universal exponent)' /
DATA carbon_dioxide(450) / ' 1.01d+00          !R0 (universal amplitude)' /
DATA carbon_dioxide(451) / ' 0.065d+00         !z (universal exponent--not used for t.c., only viscosity)' /
DATA carbon_dioxide(452) / ' 1.00d+00          !c (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA carbon_dioxide(453) / ' 1.5d-10           !xi0 (amplitude) [m]' /
DATA carbon_dioxide(454) / ' 0.052d+00         !gam0 (amplitude) [-]' /
DATA carbon_dioxide(455) / ' 0.40d-09          !qd_inverse (modified effective cutoff parameter) [m]' /
DATA carbon_dioxide(456) / ' 450.0d+00         !tref (reference temperature) [K]' /
DATA carbon_dioxide(457) / '' /
DATA carbon_dioxide(458) / '' /
DATA carbon_dioxide(459) / '#ETA               !viscosity model specification' /
DATA carbon_dioxide(460) / 'VS1  pure fluid viscosity model of Fenghour et al. (1998).' /
DATA carbon_dioxide(461) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(462) / '?Fenghour, A., Wakeham, W.A., Vesovic, V.,' /
DATA carbon_dioxide(463) / '? "The Viscosity of Carbon Dioxide,"' /
DATA carbon_dioxide(464) / '? J. Phys. Chem. Ref. Data, 27:31-44, 1998.' /
DATA carbon_dioxide(465) / '?\' /
DATA carbon_dioxide(466) / '?The uncertainty in viscosity ranges from 0.3% in the dilute gas near room' /
DATA carbon_dioxide(467) / '?temperature to 5% at the highest pressures.' /
DATA carbon_dioxide(468) / '?\' /
DATA carbon_dioxide(469) / '!end of info section' /
DATA carbon_dioxide(470) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(471) / '2000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(472) / '800000.0           !upper pressure limit [kPa]' /
DATA carbon_dioxide(473) / '37.24              !maximum density [mol/L]' /
DATA carbon_dioxide(474) / '1                  !number of terms associated with dilute-gas function' /
DATA carbon_dioxide(475) / 'CI1                !pointer to reduced effective collision cross-section model' /
DATA carbon_dioxide(476) / '1.                 !Lennard-Jones coefficient sigma [nm] (Not used for CO2)' /
DATA carbon_dioxide(477) / '251.196            !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA carbon_dioxide(478) / '1.0    1.0         !reducing parameters for T, eta' /
DATA carbon_dioxide(479) / '1.00697d0   0.50d0 !Chapman-Enskog term' /
DATA carbon_dioxide(480) / '0                  !number of terms for initial density dependence' /
DATA carbon_dioxide(481) / '0 5 0 0 0 0        !# resid terms:  close-packed density;  simple poly; numerator of rational poly; denominator of rat. poly; numerator of exponential; denominator of exponential' /
DATA carbon_dioxide(482) / '251.196 0.0227222 1.0                   !reducing parameters for T (= eps/k), rho, eta' /
DATA carbon_dioxide(483) / ' 0.4071119d-2    0.00  1.00  0.00  0  !d_11; powers of tau, del, del0; power of del in exponential [0 indicated no exponential term present]' /
DATA carbon_dioxide(484) / ' 0.7198037d-4    0.00  2.00  0.00  0  !d_21' /
DATA carbon_dioxide(485) / ' 0.2411697d-16  -3.00  6.00  0.00  0  !d_64' /
DATA carbon_dioxide(486) / ' 0.2971072d-22   0.00  8.00  0.00  0  !d_81' /
DATA carbon_dioxide(487) / '-0.1627888d-22  -1.00  8.00  0.00  0  !d_82' /
DATA carbon_dioxide(488) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA carbon_dioxide(489) / '' /
DATA carbon_dioxide(490) / '' /
DATA carbon_dioxide(491) / '#AUX               !reduced effective collision cross-section model specification' /
DATA carbon_dioxide(492) / 'CI1  reduced effective collision cross-section model (empirical form in terms of log(T*))' /
DATA carbon_dioxide(493) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(494) / '?Fenghour, A., Wakeham, W.A., Vesovic, V.,' /
DATA carbon_dioxide(495) / '? "The Viscosity of Carbon Dioxide,"' /
DATA carbon_dioxide(496) / '? J. Phys. Chem. Ref. Data, 27:31-44, 1998.' /
DATA carbon_dioxide(497) / '?\' /
DATA carbon_dioxide(498) / '!end of info section' /
DATA carbon_dioxide(499) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(500) / '2000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(501) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_dioxide(502) / '0.0                !(dummy) maximum density' /
DATA carbon_dioxide(503) / '5                  !number of terms' /
DATA carbon_dioxide(504) / ' 0.235156d0     0  !coeff, power of Tstar' /
DATA carbon_dioxide(505) / '-0.491266d0     1' /
DATA carbon_dioxide(506) / ' 5.211155d-2    2' /
DATA carbon_dioxide(507) / ' 5.347906d-2    3' /
DATA carbon_dioxide(508) / '-1.537102d-2    4' /
DATA carbon_dioxide(509) / '' /
DATA carbon_dioxide(510) / '' /
DATA carbon_dioxide(511) / '@ETA               !viscosity model specification' /
DATA carbon_dioxide(512) / 'VS4  pure fluid generalized friction theory viscosity model of Quinones-Cisneros and Deiters (2006).' /
DATA carbon_dioxide(513) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(514) / '? Quinones-Cisneros, S.E. and Deiters, U.K.' /
DATA carbon_dioxide(515) / '? "Generalization of the Friction Theory for Viscosity Modeling,"' /
DATA carbon_dioxide(516) / '? J. Phys. Chem. B, 110:12820-12834, 2006.' /
DATA carbon_dioxide(517) / '?' /
DATA carbon_dioxide(518) / '?The uncertainty in viscosity ranges from 0.3% in the dilute gas near room' /
DATA carbon_dioxide(519) / '?temperature to 5% at the highest pressures.' /
DATA carbon_dioxide(520) / '!end of info section' /
DATA carbon_dioxide(521) / '216.58             !lower temperature limit [K]' /
DATA carbon_dioxide(522) / '1000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(523) / '100000.0           !upper pressure limit [kPa]' /
DATA carbon_dioxide(524) / '37.24              !maximum density [mol/L]' /
DATA carbon_dioxide(525) / '5 0 0 0 0 0        !number of terms associated with dilute-gas function' /
DATA carbon_dioxide(526) / 'NUL                !pointer to reduced effective collision cross-section model;not used' /
DATA carbon_dioxide(527) / '0.3751             !Lennard-Jones coefficient sigma [nm] (not used)' /
DATA carbon_dioxide(528) / '251.196            !Lennard-Jones coefficient epsilon/kappa [K] (not used)' /
DATA carbon_dioxide(529) / '304.1282  1.0      !reducing parameters for T, eta' /
DATA carbon_dioxide(530) / '0.00d0   0.50d0    !Chapman-Enskog term; not used here' /
DATA carbon_dioxide(531) / '69.18424d0   0.0d0 !empirical terms for eta0' /
DATA carbon_dioxide(532) / '-215.8618d0  0.25d0' /
DATA carbon_dioxide(533) / '210.94362d0  0.5d0' /
DATA carbon_dioxide(534) / '-49.0494d0   0.75d0' /
DATA carbon_dioxide(535) / '0                  !number of terms for initial density dependence; not yet used.' /
DATA carbon_dioxide(536) / ' 1.19805d-4   -1.25861d-4  5.48871d-5   0.0d0 0.0d0    !a(0),a(1),a(2)' /
DATA carbon_dioxide(537) / ' 3.15921d-5   -2.60469d-5  7.09199d-6   0.0d0 0.0d0    !b(0),b(1),b(2)' /
DATA carbon_dioxide(538) / ' 1.80689d-5   -7.41742d-6  0.0d0        0.0d0 0.0d0    !c(0),c(1),c(2)' /
DATA carbon_dioxide(539) / '-2.31066d-9    0.0d0       5.42486d-10  0.0d0 0.0d0    !A(0),A(1),A(2)' /
DATA carbon_dioxide(540) / ' 1.04558d-8   -2.20758d-9  0.0d0        0.0d0 0.0d0    !B(0),B(1),B(2)' /
DATA carbon_dioxide(541) / ' 1.03255d-6   -8.56207d-7  3.84384d-7   0.0d0 0.0d0    !C(0),C(1),C(2)' /
DATA carbon_dioxide(542) / ' 0.0d0         0.0d0       0.0d0        0.0d0 0.0d0    !D(0),D(1),D(2)' /
DATA carbon_dioxide(543) / ' 0.0d0         0.0d0       0.0d0        0.0d0 0.0d0    !E(0),E(1),E(2)' /
DATA carbon_dioxide(544) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA carbon_dioxide(545) / '' /
DATA carbon_dioxide(546) / '' /
DATA carbon_dioxide(547) / '#STN        !surface tension specification' /
DATA carbon_dioxide(548) / 'ST1  surface tension model of Mulero et al. (2012)' /
DATA carbon_dioxide(549) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(550) / '?Mulero, A., Cachadina, I., and Parra, M.I.' /
DATA carbon_dioxide(551) / '? "Recommended Correlations for the Surface Tension of Common Fluids,"' /
DATA carbon_dioxide(552) / '? J. Phys. Chem. Ref. Data, 41, 043105, 2012.' /
DATA carbon_dioxide(553) / '?\' /
DATA carbon_dioxide(554) / '!end of info section' /
DATA carbon_dioxide(555) / '0.0                !lower temperature limit [K]' /
DATA carbon_dioxide(556) / '304.128            !upper temperature limit [K]' /
DATA carbon_dioxide(557) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_dioxide(558) / '0.0                !(dummy) maximum density' /
DATA carbon_dioxide(559) / '1                           !number of terms in surface tension model' /
DATA carbon_dioxide(560) / '304.128                     !critical temperature used in fit (dummy)' /
DATA carbon_dioxide(561) / ' 0.07863     1.254          !sigma0 and n' /
DATA carbon_dioxide(562) / '' /
DATA carbon_dioxide(563) / '' /
DATA carbon_dioxide(564) / '#DE         !dielectric constant specification' /
DATA carbon_dioxide(565) / 'DE3  dielectric constant model of Harvey and Lemmon (2005).' /
DATA carbon_dioxide(566) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(567) / '?Harvey, A.H. and Lemmon, E.W.' /
DATA carbon_dioxide(568) / '? "Method for Estimating the Dielectric Constant of Natural Gas Mixtures,"' /
DATA carbon_dioxide(569) / '? Int. J. Thermophys., 26(1):31-46, 2005.' /
DATA carbon_dioxide(570) / '?\' /
DATA carbon_dioxide(571) / '!end of info section' /
DATA carbon_dioxide(572) / '0.0                !lower temperature limit [K]' /
DATA carbon_dioxide(573) / '2000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(574) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_dioxide(575) / '0.0                !(dummy) maximum density' /
DATA carbon_dioxide(576) / '273.16 1000.0 1.0  !reducing parameters for t and d' /
DATA carbon_dioxide(577) / '0 2 4 0 0 0                         !number of terms in dielectric constant model' /
DATA carbon_dioxide(578) / ' 7.3455           0.    1.    0.    !coef, t exp, d exp' /
DATA carbon_dioxide(579) / ' 0.00335          1.    1.    0.' /
DATA carbon_dioxide(580) / ' 83.93            0.    2.    0.' /
DATA carbon_dioxide(581) / ' 145.1            1.    2.    0.' /
DATA carbon_dioxide(582) / '-578.8            0.    2.55  0.' /
DATA carbon_dioxide(583) / '-1012.0           1.    2.55  0.' /
DATA carbon_dioxide(584) / '' /
DATA carbon_dioxide(585) / '' /
DATA carbon_dioxide(586) / '#MLT        !melting line specification' /
DATA carbon_dioxide(587) / 'ML1  melting line model of Span and Wagner (1996).' /
DATA carbon_dioxide(588) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(589) / '?Span, R. and Wagner, W.,' /
DATA carbon_dioxide(590) / '? "A New Equation of State for Carbon Dioxide Covering the Fluid Region' /
DATA carbon_dioxide(591) / '? from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa,"' /
DATA carbon_dioxide(592) / '? J. Phys. Chem. Ref. Data, 25(6):1509-1596, 1996.' /
DATA carbon_dioxide(593) / '?\' /
DATA carbon_dioxide(594) / '!end of info section' /
DATA carbon_dioxide(595) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(596) / '1100.0             !upper temperature limit [K]' /
DATA carbon_dioxide(597) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_dioxide(598) / '0.0                !(dummy) maximum density' /
DATA carbon_dioxide(599) / '216.592  517.95    !reducing temperature and pressure' /
DATA carbon_dioxide(600) / '1 2 0 0 0 0                 !number of terms in melting line equation' /
DATA carbon_dioxide(601) / ' 1.             0.          !coefficients and exponents' /
DATA carbon_dioxide(602) / '1955.539        1.' /
DATA carbon_dioxide(603) / '2055.4593       2.' /
DATA carbon_dioxide(604) / '' /
DATA carbon_dioxide(605) / '' /
DATA carbon_dioxide(606) / '#SBL        !sublimation line specification' /
DATA carbon_dioxide(607) / 'SB3  sublimation line model of Span and Wagner (1996).' /
DATA carbon_dioxide(608) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(609) / '?Span, R. and Wagner, W.,' /
DATA carbon_dioxide(610) / '? "A New Equation of State for Carbon Dioxide Covering the Fluid Region' /
DATA carbon_dioxide(611) / '? from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa,"' /
DATA carbon_dioxide(612) / '? J. Phys. Chem. Ref. Data, 25(6):1509-1596, 1996.' /
DATA carbon_dioxide(613) / '?\' /
DATA carbon_dioxide(614) / '!end of info section' /
DATA carbon_dioxide(615) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(616) / '216.592            !upper temperature limit [K]' /
DATA carbon_dioxide(617) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_dioxide(618) / '0.0                !(dummy) maximum density' /
DATA carbon_dioxide(619) / '216.592 517.95     !reducing temperature and pressure' /
DATA carbon_dioxide(620) / '0 3 0 0 0 0                 !number of terms in sublimation line equation' /
DATA carbon_dioxide(621) / '-14.740846      1.          !coefficients and exponents' /
DATA carbon_dioxide(622) / '  2.4327015     1.9' /
DATA carbon_dioxide(623) / ' -5.3061778     2.9' /
DATA carbon_dioxide(624) / '' /
DATA carbon_dioxide(625) / '' /
DATA carbon_dioxide(626) / '#PS         !vapor pressure equation' /
DATA carbon_dioxide(627) / 'PS5  vapor pressure equation of Span and Wagner (1996)' /
DATA carbon_dioxide(628) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(629) / '?See EOS' /
DATA carbon_dioxide(630) / '?\' /
DATA carbon_dioxide(631) / '!end of info section' /
DATA carbon_dioxide(632) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(633) / '2000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(634) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_dioxide(635) / '0.0                !(dummy) maximum density' /
DATA carbon_dioxide(636) / '304.1282  7377.3   !reducing parameters' /
DATA carbon_dioxide(637) / '4 0 0 0 0 0        !number of terms in equation' /
DATA carbon_dioxide(638) / '-7.0602087    1.0  !coefficients and exponents' /
DATA carbon_dioxide(639) / ' 1.9391218    1.5' /
DATA carbon_dioxide(640) / '-1.6463597    2.0' /
DATA carbon_dioxide(641) / '-3.2995634    4.0' /
DATA carbon_dioxide(642) / '' /
DATA carbon_dioxide(643) / '' /
DATA carbon_dioxide(644) / '#DL         !saturated liquid density equation' /
DATA carbon_dioxide(645) / 'DL4  saturated liquid density equation of Span and Wagner (1996)' /
DATA carbon_dioxide(646) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(647) / '?See EOS' /
DATA carbon_dioxide(648) / '?\' /
DATA carbon_dioxide(649) / '!end of info section' /
DATA carbon_dioxide(650) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(651) / '2000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(652) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_dioxide(653) / '0.0                !(dummy) maximum density' /
DATA carbon_dioxide(654) / '304.1282  10.6249  !reducing parameters' /
DATA carbon_dioxide(655) / '4 0 0 0 0 0        !number of terms in equation' /
DATA carbon_dioxide(656) / ' 1.92451080 1.02   !coefficients and exponents' /
DATA carbon_dioxide(657) / '-0.62385555 1.5' /
DATA carbon_dioxide(658) / '-0.32731127 5.0' /
DATA carbon_dioxide(659) / ' 0.39245142 5.5' /
DATA carbon_dioxide(660) / '' /
DATA carbon_dioxide(661) / '' /
DATA carbon_dioxide(662) / '#DV         !saturated vapor density equation' /
DATA carbon_dioxide(663) / 'DV4  saturated vapor density equation of Span and Wagner (1996)' /
DATA carbon_dioxide(664) / '?LITERATURE REFERENCE \' /
DATA carbon_dioxide(665) / '?See EOS' /
DATA carbon_dioxide(666) / '?\' /
DATA carbon_dioxide(667) / '!end of info section' /
DATA carbon_dioxide(668) / '216.592            !lower temperature limit [K]' /
DATA carbon_dioxide(669) / '2000.0             !upper temperature limit [K]' /
DATA carbon_dioxide(670) / '0.0                !(dummy) upper pressure limit' /
DATA carbon_dioxide(671) / '0.0                !(dummy) maximum density' /
DATA carbon_dioxide(672) / '304.1282  10.6249  !reducing parameters' /
DATA carbon_dioxide(673) / '5 0 0 0 0 0        !number of terms in equation' /
DATA carbon_dioxide(674) / '-1.7074879  1.02   !coefficients and exponents' /
DATA carbon_dioxide(675) / '-0.8227467  1.5' /
DATA carbon_dioxide(676) / '-4.6008549  3.0' /
DATA carbon_dioxide(677) / '-10.111178  7.0' /
DATA carbon_dioxide(678) / '-29.742252 14.0' /
DATA carbon_dioxide(679) / '' /
DATA carbon_dioxide(680) / '' /
DATA carbon_dioxide(681) / '@END' /
DATA carbon_dioxide(682) / 'c        1         2         3         4         5         6         7         8' /
DATA carbon_dioxide(683) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
! #######################################################
character(256), TARGET :: Diethanolamine(322)
DATA Diethanolamine(1) / 'Diethanolamine       !Short name' /
DATA Diethanolamine(2) / '111-42-2             !CAS number' /
DATA Diethanolamine(3) / '2,2-Iminodiethanol  !Full name' /
DATA Diethanolamine(4) / 'HN(CH2CH2OH)2        !Chemical formula {C4H11NO2}' /
DATA Diethanolamine(5) / 'bis(2-hydroxyethyl)Amine  !Synonym' /
DATA Diethanolamine(6) / '105.1356             !Molar mass [g/mol]' /
DATA Diethanolamine(7) / '301.1                !Triple point temperature [K]' /
DATA Diethanolamine(8) / '541.234              !Normal boiling point [K]' /
DATA Diethanolamine(9) / '736.5                !Critical temperature [K]' /
DATA Diethanolamine(10) / '4950.75              !Critical pressure [kPa]' /
DATA Diethanolamine(11) / '3.3                  !Critical density [mol/L]' /
DATA Diethanolamine(12) / '1.013                !Acentric factor' /
DATA Diethanolamine(13) / '2.78994              !Dipole moment [Debye] Ikada, E., Hida, Y., Okamoto, H., Hagino, J., Koizumi, N., "Dielectric Properties of Ethanolamines, " Bull. Inst. Chem. Res., Kyoto Univ., 46, 5, 239-247 (1969).' /
DATA Diethanolamine(14) / 'NBP                  !Default reference state' /
DATA Diethanolamine(15) / '10.0                 !Version number' /
DATA Diethanolamine(16) / '3082                 !UN Number                                                 :UN:' /
DATA Diethanolamine(17) / 'other                !Family                                                    :Family:' /
DATA Diethanolamine(18) / '????                 !Heating value (upper) [kJ/mol]                            :Heat:' /
DATA Diethanolamine(19) / '1S/C4H11NO2/c6-3-1-5-2-4-7/h5-7H,1-4H2    !Standard InChI String                :InChi:' /
DATA Diethanolamine(20) / 'ZBCBWPMODOFKDW-UHFFFAOYSA-N               !Standard InChI Key                   :InChiKey:' /
DATA Diethanolamine(21) / '????                                      !Alternative fluid for mixing rules   :AltID:' /
DATA Diethanolamine(22) / '393cd060                                  !Hash number from InChI Key           :Hash:' /
DATA Diethanolamine(23) / '' /
DATA Diethanolamine(24) / '' /
DATA Diethanolamine(25) / '' /
DATA Diethanolamine(26) / '' /
DATA Diethanolamine(27) / '!The fluid files contain general information about the fluid in the first 15 to 20 lines, followed by sections for the' /
DATA Diethanolamine(28) / '! equations of state, transport equations, and auxiliary equations.  Equations of state are listed first.  The NIST recommended' /
DATA Diethanolamine(29) / '! equations begin with a hash mark (#).  The secondary equations begin with the @ symbol.  These symbols can be swapped to' /
DATA Diethanolamine(30) / '! select a secondary equation as primary and the primary as secondary.  The equation of state section also contains auxiliary' /
DATA Diethanolamine(31) / '! equations for the ideal gas heat capacity or ideal gas Helmholtz energy.  Below the equations of state (both primary and' /
DATA Diethanolamine(32) / '! secondary) are the transport equations, first viscosity and then thermal conductivity.  These are then followed by the' /
DATA Diethanolamine(33) / '! secondary equations if available.  The transport section also contains auxiliary equations required to calculate either the' /
DATA Diethanolamine(34) / '! dilute gas state or the critical enhancement.  At the end of the file are additional but not necessary auxiliary equations,' /
DATA Diethanolamine(35) / '! including simple equations for the vapor pressure, saturated liquid and vapor densities, melting line (for some fluids), and' /
DATA Diethanolamine(36) / '! sublimation line (for even fewer fluids).  This section also contains the equations for dielectric constant and surface' /
DATA Diethanolamine(37) / '! tension if available.  The sections are divided by different symbols (these being _-+=^*~) to aid the eye in locating a' /
DATA Diethanolamine(38) / '! particular section.  Secondary equations are indented 10 spaces to avoid confusion with the NIST recommended equations.  The' /
DATA Diethanolamine(39) / '! end of the fluid file is marked with @END.  Anything below that is ignored.' /
DATA Diethanolamine(40) / '' /
DATA Diethanolamine(41) / '' /
DATA Diethanolamine(42) / '! compiled by S. Herrig, Thermodynamics, Ruhr-Universitaet Bochum, Germany' /
DATA Diethanolamine(43) / '! 11-11-15  MK, Original version.' /
DATA Diethanolamine(44) / '! 03-10-18  SH, Add equation of state of Herrig et al. (2018)' /
DATA Diethanolamine(45) / '! 03-13-18 MLH, Add dipole moment, preliminary transport.' /
DATA Diethanolamine(46) / '' /
DATA Diethanolamine(47) / '' /
DATA Diethanolamine(48) / '' /
DATA Diethanolamine(49) / '' /
DATA Diethanolamine(50) / '________________________________________________________________________________' /
DATA Diethanolamine(51) / '' /
DATA Diethanolamine(52) / '#EOS   !---Equation of state---' /
DATA Diethanolamine(53) / 'FEQ    !Helmholtz equation of state for diethanolamine of Herrig et al. (2018).' /
DATA Diethanolamine(54) / ':TRUECRITICALPOINT:  736.5      3.3           !True EOS critical point [K, mol/L] (where dP/dD=0 and d^2P/dD^2=0 at constant T)' /
DATA Diethanolamine(55) / ':DOI:' /
DATA Diethanolamine(56) / '?' /
DATA Diethanolamine(57) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(58) / '?Herrig, S., Thol, M., Kortmann, M., Lemmon, E.W., and Span, R.,' /
DATA Diethanolamine(59) / '? unpublished equation, 2018.' /
DATA Diethanolamine(60) / '?' /
DATA Diethanolamine(61) / '?The experimental database that was available to fit the EOS was limited to' /
DATA Diethanolamine(62) / '? measurements in the liquid phase at atmospheric pressure.  At these conditions' /
DATA Diethanolamine(63) / '? and at temperatures up 435 K, the estimated uncertainty of calculated' /
DATA Diethanolamine(64) / '? homogeneous densities is 0.2 %. The uncertainty of calculated speed-of-sound' /
DATA Diethanolamine(65) / '? data is 1 % at temperatures between 295 K and 325 K. Calculated vapor pressures' /
DATA Diethanolamine(66) / '? are accurate to within 4 % for temperatures between 400 K and 540 K. The' /
DATA Diethanolamine(67) / '? uncertainties in all properties increase for lower and higher temperatures' /
DATA Diethanolamine(68) / '? where there are no reliable data sets to validate the equation.  Since its' /
DATA Diethanolamine(69) / '? extrapolation behavior was carefully constrained, the EOS will also give' /
DATA Diethanolamine(70) / '? qualitatively reasonable results beyond the experimentally covered regions.' /
DATA Diethanolamine(71) / '?' /
DATA Diethanolamine(72) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(73) / '301.1              !Lower temperature limit [K]' /
DATA Diethanolamine(74) / '740.0              !Upper temperature limit [K]' /
DATA Diethanolamine(75) / '5000.0             !Upper pressure limit [kPa]' /
DATA Diethanolamine(76) / '10.4               !Maximum density [mol/L]' /
DATA Diethanolamine(77) / 'CPP                                    !Pointer to Cp0 model' /
DATA Diethanolamine(78) / '105.1356                               !Molar mass [g/mol]' /
DATA Diethanolamine(79) / '301.1                                  !Triple point temperature [K]' /
DATA Diethanolamine(80) / '0.00013396                             !Pressure at triple point [kPa]' /
DATA Diethanolamine(81) / '10.39                                  !Density at triple point [mol/L]' /
DATA Diethanolamine(82) / '541.234                                !Normal boiling point temperature [K]' /
DATA Diethanolamine(83) / '1.013                                  !Acentric factor' /
DATA Diethanolamine(84) / '736.5         4950.75        3.3       !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Diethanolamine(85) / '736.5                        3.3       !Reducing parameters [K, mol/L]' /
DATA Diethanolamine(86) / '8.3144598                              !Gas constant [J/mol-K]' /
DATA Diethanolamine(87) / '  10  4   5 12   0 0    0 0  0 0  0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Diethanolamine(88) / '  0.066088158  1.0     4.  0.          !a(i),t(i),d(i),l(i)' /
DATA Diethanolamine(89) / '  6.1059245    0.507   1.  0.' /
DATA Diethanolamine(90) / ' -7.0526968    0.907   1.  0.' /
DATA Diethanolamine(91) / ' -0.29739545   1.22    2.  0.' /
DATA Diethanolamine(92) / '  0.11592105   0.649   3.  0.' /
DATA Diethanolamine(93) / ' -1.8616953    2.14    1.  2.' /
DATA Diethanolamine(94) / ' -0.97392153   2.89    3.  2.' /
DATA Diethanolamine(95) / '  0.14690655   1.54    2.  1.' /
DATA Diethanolamine(96) / ' -0.63284478   3.34    2.  2.' /
DATA Diethanolamine(97) / ' -0.037820123  0.998   7.  1.' /
DATA Diethanolamine(98) / '  2.774726     0.92    1.  2. 2.   -0.8971  -0.9691  1.216   0.6694   0. 0. 0.' /
DATA Diethanolamine(99) / ' -1.0230468    1.16    1.  2. 2.   -1.499   -1.518   0.6775  0.6466   0. 0. 0.' /
DATA Diethanolamine(100) / ' -0.19552536   1.43    3.  2. 2.   -1.681   -1.328   0.7815  0.6669   0. 0. 0.' /
DATA Diethanolamine(101) / ' -0.14997584   1.24    2.  2. 2.   -1.661   -1.23    0.8796  1.189    0. 0. 0.' /
DATA Diethanolamine(102) / ' -0.23421918   0.801   2.  2. 2.   -1.245   -1.112   1.357   1.248    0. 0. 0.' /
DATA Diethanolamine(103) / '                                      eta      beta    gamma   epsilon' /
DATA Diethanolamine(104) / '                                   EXP[eta*(delta-epsilon)^2+beta*(tau-gamma)^2]' /
DATA Diethanolamine(105) / '' /
DATA Diethanolamine(106) / '' /
DATA Diethanolamine(107) / '#AUX   !---Auxiliary function for Cp0' /
DATA Diethanolamine(108) / 'CPP    !Ideal gas heat capacity function for diethanolamine of Herrig et al. (2018).' /
DATA Diethanolamine(109) / '?' /
DATA Diethanolamine(110) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(111) / '?Herrig, S., Thol, M., Kortmann, M., Lemmon, E.W., and Span, R., 2018.' /
DATA Diethanolamine(112) / '?' /
DATA Diethanolamine(113) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(114) / '0.                 !' /
DATA Diethanolamine(115) / '10000.             !' /
DATA Diethanolamine(116) / '0.                 !' /
DATA Diethanolamine(117) / '0.                 !' /
DATA Diethanolamine(118) / '1.0     8.3144598  !Reducing parameters for T, Cp0' /
DATA Diethanolamine(119) / '1 2   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Diethanolamine(120) / '4.0       0.0' /
DATA Diethanolamine(121) / '4.25     96.0' /
DATA Diethanolamine(122) / '37.7   1165.0' /
DATA Diethanolamine(123) / '' /
DATA Diethanolamine(124) / '' /
DATA Diethanolamine(125) / '#AUX   !---Auxiliary function for PX0' /
DATA Diethanolamine(126) / 'PX0    !Helmholtz energy ideal-gas function for diethanolamine of Herrig et al. (2018).' /
DATA Diethanolamine(127) / '?' /
DATA Diethanolamine(128) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(129) / '?Herrig, S., Thol, M., Kortmann, M., Lemmon, E.W., and Span, R., 2018.' /
DATA Diethanolamine(130) / '?' /
DATA Diethanolamine(131) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(132) / '1 2  2  0 0  0 0 0               !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau))' /
DATA Diethanolamine(133) / '  3.0                   1.0      !ai, ti for [ai*log(tau**ti)] terms' /
DATA Diethanolamine(134) / ' 19.5095527798094963    0.0      !aj, ti for [ai*tau**ti] terms' /
DATA Diethanolamine(135) / ' -2.9703090858805892    1.0      !aj, ti for [ai*tau**ti] terms' /
DATA Diethanolamine(136) / ' 4.25     96.0                   !aj, ti for [ai*log(1-exp(-ti/T)] terms' /
DATA Diethanolamine(137) / ' 37.7   1165.0' /
DATA Diethanolamine(138) / '' /
DATA Diethanolamine(139) / '' /
DATA Diethanolamine(140) / '' /
DATA Diethanolamine(141) / '' /
DATA Diethanolamine(142) / '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' /
DATA Diethanolamine(143) / '' /
DATA Diethanolamine(144) / '#TRN   !---ECS Transport---' /
DATA Diethanolamine(145) / 'ECS    !Extended Corresponding States model (Propane reference)' /
DATA Diethanolamine(146) / ':DOI: 10.6028/NIST.IR.8209' /
DATA Diethanolamine(147) / '?' /
DATA Diethanolamine(148) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(149) / '?Huber, M.L., "Models for the Viscosity, Thermal Conductivity, and Surface Tension' /
DATA Diethanolamine(150) / '? of Selected Pure Fluids as Implemented in REFPROP v10.0," NISTIR 8209, 2018.' /
DATA Diethanolamine(151) / '? doi: 10.6028/NIST.IR.8209' /
DATA Diethanolamine(152) / '?' /
DATA Diethanolamine(153) / '?VISCOSITY' /
DATA Diethanolamine(154) / '? ECS parameters based on fitting the data of:' /
DATA Diethanolamine(155) / '? DiGuilio, R.M., Lee, R.J., Schaeffer, S.T., Brasher, L.L., and Teja, A.S., "Densities and Viscosities of the Ethanolamines," J. Chem. Eng. Data, 37:239-242, 1992. doi: 10.1021/je00006a028' /
DATA Diethanolamine(156) / '? Teng, T.T., Maham, Y., Hepler, L.G., and Mather, A.E., "Viscosity of Aqueous Solutions of N-Methyldiethanolamine and of Diethanolamine," J. Chem. Eng. Data, 39:290-293, 1994. doi: 10.1021/je00014a021' /
DATA Diethanolamine(157) / '? Aguila-Hernandez, J., Trejo, A., Garcia-Flores, B.E., Molnar, R., "Viscometric and Volumetric Behaviour of Binary Mixtures of Sulfolane and N-Methylpyrrolidone with Monoethanolamine and Diethanolamine in the Range 303-373 K," Fluid Phase Equilib., 267, 172-180, 2008. doi: 10.1016/j.fluid.2008.02.023' /
DATA Diethanolamine(158) / '? Blanco, A., Garcia-Abuin, A., Gomez-Diaz, D., Navaza, J., and Villaverde, O., J. Chem. Eng. Data, 58:653-659, 2013.' /
DATA Diethanolamine(159) / '? Haghtalab, A., Shojaeian, A., "Volumetric and Viscometric Behaviour of the Binary Systems of N-Methyldiethanolamine and Diethanolamine with 1-Butyl-3-Methylimidazolium Acetate at Various Temperatures," J. Chem. Thermodyn., 68:128-137, 2014. doi: 10.1016/j.jct.2013.09.001' /
DATA Diethanolamine(160) / '? The estimated uncertainty of the viscosity of the liquid phase at atmospheric pressure over the temperature range from 305 K to 423 K is 40%,' /
DATA Diethanolamine(161) / '? and the estimated uncertainty of the gas phase is 20%.' /
DATA Diethanolamine(162) / '?' /
DATA Diethanolamine(163) / '?THERMAL CONDUCTIVITY' /
DATA Diethanolamine(164) / '? ECS parameters based on fitting data of:' /
DATA Diethanolamine(165) / '? DiGuilio, R.M., McGregor, W.L., and Teja, A.S., "Thermal Conductivities of the Ethanolamines," J. Chem. Eng. Data, 37:242-245, 1992. doi: 10.1021/je00006a029' /
DATA Diethanolamine(166) / '? The estimated uncertainty of the thermal conductivity of the liquid phase at saturation over 295 K - 442 K is 2%; for the gas phase 20%, larger near critical.' /
DATA Diethanolamine(167) / '?' /
DATA Diethanolamine(168) / '?The Lennard-Jones parameters were estimated with the method of Chung, T.H., Ajlan, M., Lee, L.L., and Starling, K.E., "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties," Ind. Eng. Chem. Res., 27:671-679, 1988.' /
DATA Diethanolamine(169) / '?' /
DATA Diethanolamine(170) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(171) / '305.0              !Lower temperature limit [K] innaccurate viscosity results at lower temperatures' /
DATA Diethanolamine(172) / '740.0              !Upper temperature limit [K]' /
DATA Diethanolamine(173) / '5000.0             !Upper pressure limit [kPa]' /
DATA Diethanolamine(174) / '10.4               !Maximum density [mol/L]' /
DATA Diethanolamine(175) / 'FEQ PROPANE.FLD' /
DATA Diethanolamine(176) / 'VS1                !Model for reference fluid viscosity' /
DATA Diethanolamine(177) / 'TC1                !Model for reference fluid thermal conductivity' /
DATA Diethanolamine(178) / 'BIG                !Large molecule identifier' /
DATA Diethanolamine(179) / '0.74 0. 0. 0.      !Large molecule parameters' /
DATA Diethanolamine(180) / '1                  !Lennard-Jones flag (0 or 1) (0 => use estimates)' /
DATA Diethanolamine(181) / '0.5434             !Lennard-Jones coefficient sigma [nm] for ECS method' /
DATA Diethanolamine(182) / '584.85             !Lennard-Jones coefficient epsilon/kappa [K] for ECS method' /
DATA Diethanolamine(183) / '1  0  0            !Number of terms in f_int term in Eucken correlation, spare1, spare2' /
DATA Diethanolamine(184) / '0.00132   0. 0. 0. !Coefficient, power of T, spare1, spare2' /
DATA Diethanolamine(185) / '2  0  0            !Number of terms in psi (visc shape factor): poly,spare1,spare2' /
DATA Diethanolamine(186) / '0.996593  0. 0. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Diethanolamine(187) / '0.0399708 0. 1. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Diethanolamine(188) / '2  0  0            !Number of terms in chi (t.c. shape factor): poly,spare1,spare2' /
DATA Diethanolamine(189) / '1.47408   0. 0. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Diethanolamine(190) / '-0.123082 0. 1. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Diethanolamine(191) / 'TK3                !Pointer to critical enhancement auxiliary function' /
DATA Diethanolamine(192) / '' /
DATA Diethanolamine(193) / '' /
DATA Diethanolamine(194) / '#AUX   !---Auxiliary function for the thermal conductivity critical enhancement' /
DATA Diethanolamine(195) / 'TK3    !Simplified thermal conductivity critical enhancement for diethanolamine of Perkins et al. (2013).' /
DATA Diethanolamine(196) / '?' /
DATA Diethanolamine(197) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(198) / '?Perkins, R.A., Sengers, J.V., Abdulagatov, I.M., and Huber, M.L.,' /
DATA Diethanolamine(199) / '? "Simplified Model for the Critical Thermal-Conductivity Enhancement in Molecular Fluids,"' /
DATA Diethanolamine(200) / '? Int. J. Thermophys., 34(2):191-212, 2013. doi: 10.1007/s10765-013-1409-z' /
DATA Diethanolamine(201) / '?' /
DATA Diethanolamine(202) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(203) / '0.                 !' /
DATA Diethanolamine(204) / '10000.             !' /
DATA Diethanolamine(205) / '0.                 !' /
DATA Diethanolamine(206) / '0.                 !' /
DATA Diethanolamine(207) / '9 0 0 0            !# terms:  CO2-terms, spare, spare, spare' /
DATA Diethanolamine(208) / '1.0  1.0  1.0      !Reducing parameters for T, rho, tcx [mW/(m-K)]' /
DATA Diethanolamine(209) / '0.63               !Nu (universal exponent)' /
DATA Diethanolamine(210) / '1.239              !Gamma (universal exponent)' /
DATA Diethanolamine(211) / '1.02               !R0 (universal amplitude)' /
DATA Diethanolamine(212) / '0.063              !Z (universal exponent--not used for t.c., only viscosity)' /
DATA Diethanolamine(213) / '1.0                !C (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA Diethanolamine(214) / '0.185e-9           !Xi0 (amplitude) [m]' /
DATA Diethanolamine(215) / '0.068              !Gam0 (amplitude) [-]' /
DATA Diethanolamine(216) / '0.662e-9           !Qd_inverse (modified effective cutoff parameter) [m]; arbitrary guess' /
DATA Diethanolamine(217) / '1104.75            !Tref (reference temperature)=1.5*Tc [K]' /
DATA Diethanolamine(218) / '' /
DATA Diethanolamine(219) / '' /
DATA Diethanolamine(220) / '' /
DATA Diethanolamine(221) / '' /
DATA Diethanolamine(222) / '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' /
DATA Diethanolamine(223) / '' /
DATA Diethanolamine(224) / '#STN   !---Surface tension---' /
DATA Diethanolamine(225) / 'ST1    !Surface tension model for dioethanolamine of Huber (2018).' /
DATA Diethanolamine(226) / ':DOI: 10.6028/NIST.IR.8209' /
DATA Diethanolamine(227) / '?' /
DATA Diethanolamine(228) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(229) / '?Huber, M.L., "Models for the Viscosity, Thermal Conductivity, and Surface Tension' /
DATA Diethanolamine(230) / '? of Selected Pure Fluids as Implemented in REFPROP v10.0," NISTIR 8209, 2018.' /
DATA Diethanolamine(231) / '? doi: 10.6028/NIST.IR.8209' /
DATA Diethanolamine(232) / '?' /
DATA Diethanolamine(233) / '?Fit at NIST to data including:' /
DATA Diethanolamine(234) / '? Fu, D. and Zhong, Z.-k., "Experimental Study on the Surface Tension of Diethanolamine- N-Methyldiethanolamine-Water Mixtures," Acta Chim. Sinica, 68, 1241-1246, 2010.' /
DATA Diethanolamine(235) / '? Lopez, A. B., Garcia-Abuin, A., Gomez-Diaz, D., La Rubia, M. D., and Navaza, J. M., "Density, Speed of Sound, Viscosity, Refractive Index and Surface Tension of N-Methyl-2-Pyrrolidone + Diethanolamine (or Triethanolamine) from T = (293.15 to 323.15) K.," J. Chem. Thermodyn., 61:1-6, 2013. doi: 10.1016/j.jct.2013.01.020' /
DATA Diethanolamine(236) / '? Blanco, A., Garcia-Abuin, A., Gomez-Diaz, D., Navaza, J., Villaverde, O., "Density, Speed of Sound, Viscosity, Surface Tension, and Excess Volume of N-Ethyl-2-Pyrrolidone + Ethanolamine (or Diethanolamine or Triethanolamine) from T = (293.15 to 323.15) K," J. Chem. Eng. Data, 58:653-659, 2013. doi: 10.1021/je301123j' /
DATA Diethanolamine(237) / '?' /
DATA Diethanolamine(238) / '?Estimated uncertainty over 293 - 333 K is 2%.' /
DATA Diethanolamine(239) / '?' /
DATA Diethanolamine(240) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(241) / '0.                 !' /
DATA Diethanolamine(242) / '10000.             !' /
DATA Diethanolamine(243) / '0.                 !' /
DATA Diethanolamine(244) / '0.                 !' /
DATA Diethanolamine(245) / '1                  !Number of terms in surface tension model' /
DATA Diethanolamine(246) / '736.5              !Critical temperature (dummy)' /
DATA Diethanolamine(247) / '0.0859443  1.15945 !Sigma0 and n' /
DATA Diethanolamine(248) / '' /
DATA Diethanolamine(249) / '' /
DATA Diethanolamine(250) / '#PS    !---Vapor pressure---' /
DATA Diethanolamine(251) / 'PS5    !Vapor pressure equation for diethanolamine of Herrig et al. (2018).' /
DATA Diethanolamine(252) / '?' /
DATA Diethanolamine(253) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(254) / '?Herrig, S., Thol, M., Kortmann, M., Lemmon, E.W., and Span, R., 2018.' /
DATA Diethanolamine(255) / '?' /
DATA Diethanolamine(256) / '?Functional Form:  P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T] where Theta=1-T/Tc, Tc and Pc' /
DATA Diethanolamine(257) / '? are the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Diethanolamine(258) / '?' /
DATA Diethanolamine(259) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(260) / '0.                 !' /
DATA Diethanolamine(261) / '10000.             !' /
DATA Diethanolamine(262) / '0.                 !' /
DATA Diethanolamine(263) / '0.                 !' /
DATA Diethanolamine(264) / '736.5   4950.75    !Reducing parameters' /
DATA Diethanolamine(265) / '5 0 0 0 0 0        !Number of terms in equation' /
DATA Diethanolamine(266) / '-14.957   1.0      !Coefficients and exponents' /
DATA Diethanolamine(267) / ' 10.220   1.33' /
DATA Diethanolamine(268) / '-5.9784   1.8' /
DATA Diethanolamine(269) / '-4.9680   3.0' /
DATA Diethanolamine(270) / '-3.1300  10.3' /
DATA Diethanolamine(271) / '' /
DATA Diethanolamine(272) / '' /
DATA Diethanolamine(273) / '#DL    !---Saturated liquid density---' /
DATA Diethanolamine(274) / 'DL1    !Saturated liquid density equation for diethanolamine of Herrig et al. (2018).' /
DATA Diethanolamine(275) / '?' /
DATA Diethanolamine(276) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(277) / '?Herrig, S., Thol, M., Kortmann, M., Lemmon, E.W., and Span, R., 2018.' /
DATA Diethanolamine(278) / '?' /
DATA Diethanolamine(279) / '?Functional Form:  D=Dc*[1+SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Diethanolamine(280) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Diethanolamine(281) / '?' /
DATA Diethanolamine(282) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(283) / '0.                 !' /
DATA Diethanolamine(284) / '10000.             !' /
DATA Diethanolamine(285) / '0.                 !' /
DATA Diethanolamine(286) / '0.                 !' /
DATA Diethanolamine(287) / '736.5   3.3        !Reducing parameters' /
DATA Diethanolamine(288) / '5 0 0 0 0 0        !Number of terms in equation' /
DATA Diethanolamine(289) / '-0.1814   0.02     !Coefficients and exponents' /
DATA Diethanolamine(290) / ' 2.3567   0.26' /
DATA Diethanolamine(291) / ' 1.8100   1.43' /
DATA Diethanolamine(292) / '-3.0300   2.0' /
DATA Diethanolamine(293) / ' 1.8740   2.6' /
DATA Diethanolamine(294) / '' /
DATA Diethanolamine(295) / '' /
DATA Diethanolamine(296) / '#DV    !---Saturated vapor density---' /
DATA Diethanolamine(297) / 'DV3    !Saturated vapor density equation for diethanolamine of Herrig et al. (2018).' /
DATA Diethanolamine(298) / '?' /
DATA Diethanolamine(299) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(300) / '?Herrig, S., Thol, M., Kortmann, M., Lemmon, E.W., and Span, R., 2018.' /
DATA Diethanolamine(301) / '?' /
DATA Diethanolamine(302) / '?Functional Form:  D=Dc*EXP[SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Diethanolamine(303) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Diethanolamine(304) / '?' /
DATA Diethanolamine(305) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Diethanolamine(306) / '0.                 !' /
DATA Diethanolamine(307) / '10000.             !' /
DATA Diethanolamine(308) / '0.                 !' /
DATA Diethanolamine(309) / '0.                 !' /
DATA Diethanolamine(310) / '736.5   3.3        !Reducing parameters' /
DATA Diethanolamine(311) / '6 0 0 0 0 0        !Number of terms in equation' /
DATA Diethanolamine(312) / ' 0.108    0.09     !Coefficients and exponents' /
DATA Diethanolamine(313) / '-4.888    0.395' /
DATA Diethanolamine(314) / '-12.70    1.467' /
DATA Diethanolamine(315) / '-47.505   3.7' /
DATA Diethanolamine(316) / '-93.42    8.05' /
DATA Diethanolamine(317) / '-221.0   16.2' /
DATA Diethanolamine(318) / '' /
DATA Diethanolamine(319) / '' /
DATA Diethanolamine(320) / '@END' /
DATA Diethanolamine(321) / 'c        1         2         3         4         5         6         7         8' /
DATA Diethanolamine(322) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
! #######################################################
character(256), TARGET :: Hydrogen_sulfide(756)
DATA Hydrogen_sulfide(1) / 'Hydrogen sulfide     !Short name' /
DATA Hydrogen_sulfide(2) / '7783-06-4            !CAS number' /
DATA Hydrogen_sulfide(3) / 'Hydrogen sulfide     !Full name' /
DATA Hydrogen_sulfide(4) / 'H2S                  !Chemical formula {H2S}' /
DATA Hydrogen_sulfide(5) / 'Dihydrogen monosulfide !Synonym' /
DATA Hydrogen_sulfide(6) / '34.08088             !Molar mass [g/mol]' /
DATA Hydrogen_sulfide(7) / '187.7                !Triple point temperature [K]' /
DATA Hydrogen_sulfide(8) / '212.85               !Normal boiling point [K]' /
DATA Hydrogen_sulfide(9) / '373.1                !Critical temperature [K]' /
DATA Hydrogen_sulfide(10) / '9000.0               !Critical pressure [kPa]' /
DATA Hydrogen_sulfide(11) / '10.19                !Critical density [mol/L]' /
DATA Hydrogen_sulfide(12) / '0.1005               !Acentric factor' /
DATA Hydrogen_sulfide(13) / '0.97                 !Dipole moment [Debye]; R.D. Nelson, D.R. Lide, and A.A. Maryott, "Selected Values of Electric Dipole Moments for Molecules in the Gas Phase," NSRDS-NBS 10, National Reference Data Series, US Government Printing Office, Washington, 1967.' /
DATA Hydrogen_sulfide(14) / 'NBP                  !Default reference state' /
DATA Hydrogen_sulfide(15) / '10.0                 !Version number' /
DATA Hydrogen_sulfide(16) / '1053                 !UN Number                                                 :UN:' /
DATA Hydrogen_sulfide(17) / 'other                !Family                                                    :Family:' /
DATA Hydrogen_sulfide(18) / '562.01               !Heating value (upper) [kJ/mol]                            :Heat:' /
DATA Hydrogen_sulfide(19) / '1S/H2S/h1H2                               !Standard InChI String                :InChi:' /
DATA Hydrogen_sulfide(20) / 'RWSOTUBLDIXVET-UHFFFAOYSA-N               !Standard InChI Key                   :InChiKey:' /
DATA Hydrogen_sulfide(21) / '????                                      !Alternative fluid for mixing rules   :AltID:' /
DATA Hydrogen_sulfide(22) / 'c6c03020                                  !Hash number from InChI Key           :Hash:' /
DATA Hydrogen_sulfide(23) / '' /
DATA Hydrogen_sulfide(24) / '' /
DATA Hydrogen_sulfide(25) / '' /
DATA Hydrogen_sulfide(26) / '' /
DATA Hydrogen_sulfide(27) / '!The fluid files contain general information about the fluid in the first 15 to 20 lines, followed by sections for the' /
DATA Hydrogen_sulfide(28) / '! equations of state, transport equations, and auxiliary equations.  Equations of state are listed first.  The NIST recommended' /
DATA Hydrogen_sulfide(29) / '! equations begin with a hash mark (#).  The secondary equations begin with the @ symbol.  These symbols can be swapped to' /
DATA Hydrogen_sulfide(30) / '! select a secondary equation as primary and the primary as secondary.  The equation of state section also contains auxiliary' /
DATA Hydrogen_sulfide(31) / '! equations for the ideal gas heat capacity or ideal gas Helmholtz energy.  Below the equations of state (both primary and' /
DATA Hydrogen_sulfide(32) / '! secondary) are the transport equations, first viscosity and then thermal conductivity.  These are then followed by the' /
DATA Hydrogen_sulfide(33) / '! secondary equations if available.  The transport section also contains auxiliary equations required to calculate either the' /
DATA Hydrogen_sulfide(34) / '! dilute gas state or the critical enhancement.  At the end of the file are additional but not necessary auxiliary equations,' /
DATA Hydrogen_sulfide(35) / '! including simple equations for the vapor pressure, saturated liquid and vapor densities, melting line (for some fluids), and' /
DATA Hydrogen_sulfide(36) / '! sublimation line (for even fewer fluids).  This section also contains the equations for dielectric constant and surface' /
DATA Hydrogen_sulfide(37) / '! tension if available.  The sections are divided by different symbols (these being _-+=^*~) to aid the eye in locating a' /
DATA Hydrogen_sulfide(38) / '! particular section.  Secondary equations are indented 10 spaces to avoid confusion with the NIST recommended equations.  The' /
DATA Hydrogen_sulfide(39) / '! end of the fluid file is marked with @END.  Anything below that is ignored.' /
DATA Hydrogen_sulfide(40) / '' /
DATA Hydrogen_sulfide(41) / '' /
DATA Hydrogen_sulfide(42) / '! compiled by E.W. Lemmon, NIST Physical and Chemical Properties Division, Boulder, Colorado' /
DATA Hydrogen_sulfide(43) / '! 04-02-98 EWL, Original version.' /
DATA Hydrogen_sulfide(44) / '! 11-18-98 EWL, Add equation of state of Polt et al. (1992).' /
DATA Hydrogen_sulfide(45) / '! 03-07-00 EWL, Add DDMIX transport properties.' /
DATA Hydrogen_sulfide(46) / '! 06-25-01 EWL, Add Lemmon and Span short EOS.' /
DATA Hydrogen_sulfide(47) / '! 03-13-03 EWL, Replace cp0 equation.' /
DATA Hydrogen_sulfide(48) / '! 03-12-04 EWL, Update EOS.' /
DATA Hydrogen_sulfide(49) / '! 04-19-04 AHH, Change dipole moment.' /
DATA Hydrogen_sulfide(50) / '! 05-28-04 MLH, Add TK3.' /
DATA Hydrogen_sulfide(51) / '! 08-26-04 EWL, Add Sakoda equation of state.' /
DATA Hydrogen_sulfide(52) / '! 12-02-06 MLH, Update LJ for ECS.' /
DATA Hydrogen_sulfide(53) / '! 03-05-07 MLH, Add FT for viscosity.' /
DATA Hydrogen_sulfide(54) / '! 06-28-09 MLH, Refit dilute gas thermal conductivity with DIPPR numbers after E. Vogel demonstrated it was off by 20%.' /
DATA Hydrogen_sulfide(55) / '! 11-14-09 EWL, Duplicate FEQ as FEK and use PHK so as to work with GERG-2008.' /
DATA Hydrogen_sulfide(56) / '! 12-06-12 EWL, Add surface tension coefficients of Mulero et al. (2012).' /
DATA Hydrogen_sulfide(57) / '! 02-05-17 EWL, Revise ancillaries to fit EOS calculations better.' /
DATA Hydrogen_sulfide(58) / '! 11-12-17 MLH, Replace old DDMIX thermal conductivity with preliminary ecs' /
DATA Hydrogen_sulfide(59) / '! 02-28-18 IHB, Add sublimation line model.' /
DATA Hydrogen_sulfide(60) / '' /
DATA Hydrogen_sulfide(61) / '' /
DATA Hydrogen_sulfide(62) / '' /
DATA Hydrogen_sulfide(63) / '' /
DATA Hydrogen_sulfide(64) / '________________________________________________________________________________' /
DATA Hydrogen_sulfide(65) / '' /
DATA Hydrogen_sulfide(66) / '#EOS   !---Equation of state---' /
DATA Hydrogen_sulfide(67) / 'FEQ    !Helmholtz equation of state for hydrogen sulfide of Lemmon and Span (2006).' /
DATA Hydrogen_sulfide(68) / ':TRUECRITICALPOINT:  373.1     10.19          !True EOS critical point [K, mol/L] (where dP/dD=0 and d^2P/dD^2=0 at constant T)' /
DATA Hydrogen_sulfide(69) / ':DOI: 10.1021/je050186n' /
DATA Hydrogen_sulfide(70) / '?' /
DATA Hydrogen_sulfide(71) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(72) / '?Lemmon, E.W. and Span, R.,' /
DATA Hydrogen_sulfide(73) / '? "Short Fundamental Equations of State for 20 Industrial Fluids,"' /
DATA Hydrogen_sulfide(74) / '? J. Chem. Eng. Data, 51(3):785-850, 2006. doi: 10.1021/je050186n' /
DATA Hydrogen_sulfide(75) / '?' /
DATA Hydrogen_sulfide(76) / '?The uncertainties in density are 0.1% in the liquid phase below the' /
DATA Hydrogen_sulfide(77) / '? critical temperature, 0.4% in the vapor phase, 1% at supercritical' /
DATA Hydrogen_sulfide(78) / '? temperatures up to 500 K, and 2.5% at higher temperatures.  Uncertainties' /
DATA Hydrogen_sulfide(79) / '? will be higher near the critical point. The uncertainty in vapor pressure' /
DATA Hydrogen_sulfide(80) / '? is 0.25%, and the uncertainty in heat capacities is estimated to be 1%.' /
DATA Hydrogen_sulfide(81) / '?' /
DATA Hydrogen_sulfide(82) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(83) / '187.7              !Lower temperature limit [K]' /
DATA Hydrogen_sulfide(84) / '760.0              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(85) / '170000.0           !Upper pressure limit [kPa]' /
DATA Hydrogen_sulfide(86) / '29.12              !Maximum density [mol/L]' /
DATA Hydrogen_sulfide(87) / 'CPP                                    !Pointer to Cp0 model' /
DATA Hydrogen_sulfide(88) / '34.08088                               !Molar mass [g/mol]' /
DATA Hydrogen_sulfide(89) / '187.7                                  !Triple point temperature [K]' /
DATA Hydrogen_sulfide(90) / '23.25                                  !Pressure at triple point [kPa]' /
DATA Hydrogen_sulfide(91) / '29.12                                  !Density at triple point [mol/L]' /
DATA Hydrogen_sulfide(92) / '212.85                                 !Normal boiling point temperature [K]' /
DATA Hydrogen_sulfide(93) / '0.1005                                 !Acentric factor' /
DATA Hydrogen_sulfide(94) / '373.1         9000.0      10.19        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Hydrogen_sulfide(95) / '373.1                     10.19        !Reducing parameters [K, mol/L]' /
DATA Hydrogen_sulfide(96) / '8.314472                               !Gas constant [J/mol-K]' /
DATA Hydrogen_sulfide(97) / '  12  4   0 0    0 0    0 0  0 0  0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Hydrogen_sulfide(98) / ' 0.87641     0.25    1.  0.            !a(i),t(i),d(i),l(i)' /
DATA Hydrogen_sulfide(99) / '-2.0367      1.125   1.  0.' /
DATA Hydrogen_sulfide(100) / ' 0.21634     1.5     1.  0.' /
DATA Hydrogen_sulfide(101) / '-0.050199    1.375   2.  0.' /
DATA Hydrogen_sulfide(102) / ' 0.066994    0.25    3.  0.' /
DATA Hydrogen_sulfide(103) / ' 0.00019076  0.875   7.  0.' /
DATA Hydrogen_sulfide(104) / ' 0.20227     0.625   2.  1.' /
DATA Hydrogen_sulfide(105) / '-0.0045348   1.75    5.  1.' /
DATA Hydrogen_sulfide(106) / '-0.22230     3.625   1.  2.' /
DATA Hydrogen_sulfide(107) / '-0.034714    3.625   4.  2.' /
DATA Hydrogen_sulfide(108) / '-0.014885   14.5     3.  3.' /
DATA Hydrogen_sulfide(109) / ' 0.0074154  12.0     4.  3.' /
DATA Hydrogen_sulfide(110) / '' /
DATA Hydrogen_sulfide(111) / '' /
DATA Hydrogen_sulfide(112) / '#AUX   !---Auxiliary function for Cp0' /
DATA Hydrogen_sulfide(113) / 'CPP    !Ideal gas heat capacity function for hydrogen sulfide of Lemmon and Span (2006).' /
DATA Hydrogen_sulfide(114) / '?' /
DATA Hydrogen_sulfide(115) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(116) / '?Lemmon, E.W. and Span, R., 2006.' /
DATA Hydrogen_sulfide(117) / '?' /
DATA Hydrogen_sulfide(118) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(119) / '0.                 !' /
DATA Hydrogen_sulfide(120) / '10000.             !' /
DATA Hydrogen_sulfide(121) / '0.                 !' /
DATA Hydrogen_sulfide(122) / '0.                 !' /
DATA Hydrogen_sulfide(123) / '1.0     8.314472   !Reducing parameters for T, Cp0' /
DATA Hydrogen_sulfide(124) / '2 2   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Hydrogen_sulfide(125) / ' 4.0               0.0' /
DATA Hydrogen_sulfide(126) / ' 0.0000014327       1.5' /
DATA Hydrogen_sulfide(127) / ' 1.1364            1823.0' /
DATA Hydrogen_sulfide(128) / ' 1.9721            3965.0' /
DATA Hydrogen_sulfide(129) / '' /
DATA Hydrogen_sulfide(130) / '' /
DATA Hydrogen_sulfide(131) / '#AUX   !---Auxiliary function for PX0' /
DATA Hydrogen_sulfide(132) / 'PX0    !Helmholtz energy ideal-gas function for hydrogen sulfide of Lemmon and Span (2006).' /
DATA Hydrogen_sulfide(133) / '?' /
DATA Hydrogen_sulfide(134) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(135) / '?Lemmon, E.W. and Span, R., 2006.' /
DATA Hydrogen_sulfide(136) / '?' /
DATA Hydrogen_sulfide(137) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(138) / '1 3  2  0 0  0 0 0               !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau))' /
DATA Hydrogen_sulfide(139) / '  3.0                   1.0      !ai, ti for [ai*log(tau**ti)] terms' /
DATA Hydrogen_sulfide(140) / ' -4.0740759852320352    0.0      !aj, ti for [ai*tau**ti] terms' /
DATA Hydrogen_sulfide(141) / '  3.7632130988924168    1.0      !aj, ti for [ai*tau**ti] terms' /
DATA Hydrogen_sulfide(142) / '  0.0000014327   -1.5' /
DATA Hydrogen_sulfide(143) / '  1.1364        1823.0           !aj, ti for [ai*log(1-exp(-ti/T)] terms' /
DATA Hydrogen_sulfide(144) / '  1.9721        3965.0' /
DATA Hydrogen_sulfide(145) / '' /
DATA Hydrogen_sulfide(146) / '' /
DATA Hydrogen_sulfide(147) / '#AUX   !---Auxiliary function for PH0' /
DATA Hydrogen_sulfide(148) / 'PH0    !Ideal gas Helmholtz form for hydrogen sulfide.' /
DATA Hydrogen_sulfide(149) / '?' /
DATA Hydrogen_sulfide(150) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(151) / '?Lemmon, E.W. and Span, R., 2006.' /
DATA Hydrogen_sulfide(152) / '?' /
DATA Hydrogen_sulfide(153) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(154) / '0.                 !' /
DATA Hydrogen_sulfide(155) / '10000.             !' /
DATA Hydrogen_sulfide(156) / '0.                 !' /
DATA Hydrogen_sulfide(157) / '0.                 !' /
DATA Hydrogen_sulfide(158) / '1 3  2  0 0  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA Hydrogen_sulfide(159) / ' 3.0               1.0                 !ai, ti for [ai*log(tau**ti)] terms' /
DATA Hydrogen_sulfide(160) / '-4.0740770957      0.0                 !aj, ti for [ai*tau**ti] terms' /
DATA Hydrogen_sulfide(161) / ' 3.7632137341      1.0' /
DATA Hydrogen_sulfide(162) / '-0.0027533528     -1.5' /
DATA Hydrogen_sulfide(163) / ' 1.1364           -4.8860895202        !aj, ti for [ai*log(1-exp(ti*tau)] terms' /
DATA Hydrogen_sulfide(164) / ' 1.9721           -10.6271777003' /
DATA Hydrogen_sulfide(165) / '' /
DATA Hydrogen_sulfide(166) / '' /
DATA Hydrogen_sulfide(167) / '' /
DATA Hydrogen_sulfide(168) / '' /
DATA Hydrogen_sulfide(169) / '--------------------------------------------------------------------------------' /
DATA Hydrogen_sulfide(170) / '' /
DATA Hydrogen_sulfide(171) / '@EOS    !---Equation of state---' /
DATA Hydrogen_sulfide(172) / 'FEK     !Helmholtz equation of state for hydrogen sulfide of Lemmon and Span (2006).' /
DATA Hydrogen_sulfide(173) / '          ?' /
DATA Hydrogen_sulfide(174) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(175) / '          ?Lemmon, E.W. and Span, R.,' /
DATA Hydrogen_sulfide(176) / '          ? "Short Fundamental Equations of State for 20 Industrial Fluids,"' /
DATA Hydrogen_sulfide(177) / '          ? J. Chem. Eng. Data, 51(3):785-850, 2006. doi: 10.1021/je050186n' /
DATA Hydrogen_sulfide(178) / '          ?' /
DATA Hydrogen_sulfide(179) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(180) / '          187.7              !Lower temperature limit [K]' /
DATA Hydrogen_sulfide(181) / '          760.0              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(182) / '          170000.0           !Upper pressure limit [kPa]' /
DATA Hydrogen_sulfide(183) / '          29.12              !Maximum density [mol/L]' /
DATA Hydrogen_sulfide(184) / '          PHK                                    !Pointer to Cp0 model' /
DATA Hydrogen_sulfide(185) / '          34.08088                               !Molar mass [g/mol]' /
DATA Hydrogen_sulfide(186) / '          187.7                                  !Triple point temperature [K]' /
DATA Hydrogen_sulfide(187) / '          23.3                                   !Pressure at triple point [kPa]' /
DATA Hydrogen_sulfide(188) / '          29.12                                  !Density at triple point [mol/L]' /
DATA Hydrogen_sulfide(189) / '          212.85                                 !Normal boiling point temperature [K]' /
DATA Hydrogen_sulfide(190) / '          0.1005                                 !Acentric factor' /
DATA Hydrogen_sulfide(191) / '          373.1         9000.0      10.19        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Hydrogen_sulfide(192) / '          373.1                     10.19        !Reducing parameters [K, mol/L]' /
DATA Hydrogen_sulfide(193) / '          8.314472                               !Gas constant [J/mol-K]' /
DATA Hydrogen_sulfide(194) / '            12  4    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Hydrogen_sulfide(195) / '           0.87641     0.25    1.  0.            !a(i),t(i),d(i),l(i)' /
DATA Hydrogen_sulfide(196) / '          -2.0367      1.125   1.  0.' /
DATA Hydrogen_sulfide(197) / '           0.21634     1.5     1.  0.' /
DATA Hydrogen_sulfide(198) / '          -0.050199    1.375   2.  0.' /
DATA Hydrogen_sulfide(199) / '           0.066994    0.25    3.  0.' /
DATA Hydrogen_sulfide(200) / '           0.00019076  0.875   7.  0.' /
DATA Hydrogen_sulfide(201) / '           0.20227     0.625   2.  1.' /
DATA Hydrogen_sulfide(202) / '          -0.0045348   1.75    5.  1.' /
DATA Hydrogen_sulfide(203) / '          -0.22230     3.625   1.  2.' /
DATA Hydrogen_sulfide(204) / '          -0.034714    3.625   4.  2.' /
DATA Hydrogen_sulfide(205) / '          -0.014885   14.5     3.  3.' /
DATA Hydrogen_sulfide(206) / '           0.0074154  12.0     4.  3.' /
DATA Hydrogen_sulfide(207) / '' /
DATA Hydrogen_sulfide(208) / '' /
DATA Hydrogen_sulfide(209) / '@AUX    !---Auxiliary function for PH0' /
DATA Hydrogen_sulfide(210) / 'PHK     !Ideal gas Helmholtz form for hydrogen sulfide of Kunz and Wagner (2004).' /
DATA Hydrogen_sulfide(211) / '          ?' /
DATA Hydrogen_sulfide(212) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(213) / '          ?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA Hydrogen_sulfide(214) / '          ? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA Hydrogen_sulfide(215) / '          ? and Other Mixtures," GERG Technical Monograph 15,' /
DATA Hydrogen_sulfide(216) / '          ? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA Hydrogen_sulfide(217) / '          ?' /
DATA Hydrogen_sulfide(218) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(219) / '          0.                 !' /
DATA Hydrogen_sulfide(220) / '          10000.             !' /
DATA Hydrogen_sulfide(221) / '          0.                 !' /
DATA Hydrogen_sulfide(222) / '          0.                 !' /
DATA Hydrogen_sulfide(223) / '          1 2  0 1 1  0 0 0  !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA Hydrogen_sulfide(224) / '           3.0               1.0                 !ai, ti for [ai*log(tau**ti)] terms' /
DATA Hydrogen_sulfide(225) / '           9.336197742       0.0                 !aj, ti for [ai*tau**ti] terms' /
DATA Hydrogen_sulfide(226) / '          -16.266508995      1.0' /
DATA Hydrogen_sulfide(227) / '          -1.00243           2.27065398          !aj, ti for cosh and sinh terms' /
DATA Hydrogen_sulfide(228) / '           3.11942           4.914580541' /
DATA Hydrogen_sulfide(229) / '' /
DATA Hydrogen_sulfide(230) / '' /
DATA Hydrogen_sulfide(231) / '@EOS    !---Equation of state---' /
DATA Hydrogen_sulfide(232) / 'FE1     !Helmholtz equation of state for hydrogen sulfide of Sakoda and Uematsu (2004).' /
DATA Hydrogen_sulfide(233) / '          ?' /
DATA Hydrogen_sulfide(234) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(235) / '          ?Sakoda, N. and Uematsu, M.' /
DATA Hydrogen_sulfide(236) / '          ? "A Thermodynamic Property Model for Fluid Phase Hydrogen Sulfide,"' /
DATA Hydrogen_sulfide(237) / '          ? Int. J. Thermophys., 25(3):709-737, 2004.' /
DATA Hydrogen_sulfide(238) / '          ?' /
DATA Hydrogen_sulfide(239) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(240) / '          187.67             !Lower temperature limit [K]' /
DATA Hydrogen_sulfide(241) / '          760.0              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(242) / '          170000.0           !Upper pressure limit [kPa]' /
DATA Hydrogen_sulfide(243) / '          29.13              !Maximum density [mol/L]' /
DATA Hydrogen_sulfide(244) / '          PH1                                    !Pointer to Cp0 model' /
DATA Hydrogen_sulfide(245) / '          34.08088                               !Molar mass [g/mol]' /
DATA Hydrogen_sulfide(246) / '          187.67                                 !Triple point temperature [K]' /
DATA Hydrogen_sulfide(247) / '          23.3                                   !Pressure at triple point [kPa]' /
DATA Hydrogen_sulfide(248) / '          29.12                                  !Density at triple point [mol/L]' /
DATA Hydrogen_sulfide(249) / '          212.88                                 !Normal boiling point temperature [K]' /
DATA Hydrogen_sulfide(250) / '          0.1039                                 !Acentric factor' /
DATA Hydrogen_sulfide(251) / '          373.37        8962.91     10.2         !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Hydrogen_sulfide(252) / '          373.37                    10.2         !Reducing parameters [K, mol/L]' /
DATA Hydrogen_sulfide(253) / '          8.314472                               !Gas constant [J/mol-K]' /
DATA Hydrogen_sulfide(254) / '            23  4    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Hydrogen_sulfide(255) / '           0.1545780       0.241     1.  0.      !a(i),t(i),d(i),l(i)' /
DATA Hydrogen_sulfide(256) / '          -1.717693        0.705     1.  0.' /
DATA Hydrogen_sulfide(257) / '          -1.595211        1.0       1.  0.' /
DATA Hydrogen_sulfide(258) / '           2.046589        0.626     2.  0.' /
DATA Hydrogen_sulfide(259) / '          -1.690358        1.12      2.  0.' /
DATA Hydrogen_sulfide(260) / '           0.9483623       1.630     2.  0.' /
DATA Hydrogen_sulfide(261) / '          -0.06800772      0.21      3.  0.' /
DATA Hydrogen_sulfide(262) / '           0.004372273     3.08      4.  0.' /
DATA Hydrogen_sulfide(263) / '           0.3788552e-4    0.827     8.  0.' /
DATA Hydrogen_sulfide(264) / '          -0.368098e-4     3.05      9.  0.' /
DATA Hydrogen_sulfide(265) / '           0.8710726e-5    3.05     10.  0.' /
DATA Hydrogen_sulfide(266) / '           0.6886876       0.110     1.  1.' /
DATA Hydrogen_sulfide(267) / '           2.751922        1.07      1.  1.' /
DATA Hydrogen_sulfide(268) / '          -1.492558        1.95      1.  1.' /
DATA Hydrogen_sulfide(269) / '           0.9202832       0.142     2.  1.' /
DATA Hydrogen_sulfide(270) / '          -0.2103469       2.130     5.  1.' /
DATA Hydrogen_sulfide(271) / '           0.001084359     4.92      1.  2.' /
DATA Hydrogen_sulfide(272) / '           0.03754723      1.75      4.  2.' /
DATA Hydrogen_sulfide(273) / '          -0.05885793      3.97      4.  2.' /
DATA Hydrogen_sulfide(274) / '          -0.02329265     11.8       3.  3.' /
DATA Hydrogen_sulfide(275) / '          -0.00012726     10.0       8.  3.' /
DATA Hydrogen_sulfide(276) / '          -0.01336824      9.83      2.  4.' /
DATA Hydrogen_sulfide(277) / '           0.01053057     14.2       3.  4.' /
DATA Hydrogen_sulfide(278) / '' /
DATA Hydrogen_sulfide(279) / '' /
DATA Hydrogen_sulfide(280) / '@AUX    !---Auxiliary function for PH0' /
DATA Hydrogen_sulfide(281) / 'PH1     !Ideal gas Helmholtz form for hydrogen sulfide.' /
DATA Hydrogen_sulfide(282) / '          ?' /
DATA Hydrogen_sulfide(283) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(284) / '          ?Sakoda, N., Uematsu, M.' /
DATA Hydrogen_sulfide(285) / '          ?' /
DATA Hydrogen_sulfide(286) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(287) / '          0.                 !' /
DATA Hydrogen_sulfide(288) / '          10000.             !' /
DATA Hydrogen_sulfide(289) / '          0.                 !' /
DATA Hydrogen_sulfide(290) / '          0.                 !' /
DATA Hydrogen_sulfide(291) / '          1 2  2 0 0  0 0 0  !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA Hydrogen_sulfide(292) / '           3.0               1.0                 !ai, ti for [ai*log(tau**ti)] terms' /
DATA Hydrogen_sulfide(293) / '           7.881037          0.0                 !aj, ti for [ai*tau**ti] terms' /
DATA Hydrogen_sulfide(294) / '          -3.20986           1.0' /
DATA Hydrogen_sulfide(295) / '           0.9767422        -4.506266            !aj, ti for [ai*log(1-exp(ti*tau)] terms' /
DATA Hydrogen_sulfide(296) / '           2.151898         -10.15526' /
DATA Hydrogen_sulfide(297) / '' /
DATA Hydrogen_sulfide(298) / '' /
DATA Hydrogen_sulfide(299) / '@EOS    !---Equation of state---' /
DATA Hydrogen_sulfide(300) / 'FE2     !Helmholtz equation of state for hydrogen sulfide of Polt et al. (1992).' /
DATA Hydrogen_sulfide(301) / '          ?' /
DATA Hydrogen_sulfide(302) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(303) / '          ?Polt, A., Platzer, B., and Maurer, G.,' /
DATA Hydrogen_sulfide(304) / '          ? "Parameter der thermischen Zustandsgleichung von Bender fuer 14' /
DATA Hydrogen_sulfide(305) / '          ? mehratomige reine Stoffe,"' /
DATA Hydrogen_sulfide(306) / '          ? Chem. Tech. (Leipzig), 44(6):216-224, 1992.' /
DATA Hydrogen_sulfide(307) / '          ?' /
DATA Hydrogen_sulfide(308) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(309) / '          187.7              !Lower temperature limit [K]' /
DATA Hydrogen_sulfide(310) / '          760.0              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(311) / '          142000.0           !Upper pressure limit [kPa]' /
DATA Hydrogen_sulfide(312) / '          29.1               !Maximum density [mol/L]' /
DATA Hydrogen_sulfide(313) / '          CP2                                    !Pointer to Cp0 model' /
DATA Hydrogen_sulfide(314) / '          34.076                                 !Molar mass [g/mol]' /
DATA Hydrogen_sulfide(315) / '          187.7                                  !Triple point temperature [K]' /
DATA Hydrogen_sulfide(316) / '          23.85                                  !Pressure at triple point [kPa]' /
DATA Hydrogen_sulfide(317) / '          29.07                                  !Density at triple point [mol/L]' /
DATA Hydrogen_sulfide(318) / '          212.84                                 !Normal boiling point temperature [K]' /
DATA Hydrogen_sulfide(319) / '          0.0956                                 !Acentric factor' /
DATA Hydrogen_sulfide(320) / '          373.6         9008.0      10.18312     !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Hydrogen_sulfide(321) / '          373.6                     10.18312     !Reducing parameters [K, mol/L]' /
DATA Hydrogen_sulfide(322) / '          8.3143                                 !Gas constant [J/mol-K]' /
DATA Hydrogen_sulfide(323) / '            22  5    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Hydrogen_sulfide(324) / '           1.35782366339      3.  0.  0.  0.     !a(i),t(i),d(i),l(i)' /
DATA Hydrogen_sulfide(325) / '          -1.53224981014      4.  0.  0.  0.' /
DATA Hydrogen_sulfide(326) / '           0.329107661253     5.  0.  0.  0.' /
DATA Hydrogen_sulfide(327) / '           1.95802782279      0.  1.  0.  0.' /
DATA Hydrogen_sulfide(328) / '          -3.01125182071      1.  1.  0.  0.' /
DATA Hydrogen_sulfide(329) / '          -1.26614059078      2.  1.  0.  0.' /
DATA Hydrogen_sulfide(330) / '           1.29960331548      3.  1.  0.  0.' /
DATA Hydrogen_sulfide(331) / '          -0.185645977138     4.  1.  0.  0.' /
DATA Hydrogen_sulfide(332) / '          -1.60919744092      0.  2.  0.  0.' /
DATA Hydrogen_sulfide(333) / '           2.34395817019      1.  2.  0.  0.' /
DATA Hydrogen_sulfide(334) / '          -0.378573094883     2.  2.  0.  0.' /
DATA Hydrogen_sulfide(335) / '           0.758423219040     0.  3.  0.  0.' /
DATA Hydrogen_sulfide(336) / '          -0.973372615169     1.  3.  0.  0.' /
DATA Hydrogen_sulfide(337) / '          -0.120786235447     0.  4.  0.  0.' /
DATA Hydrogen_sulfide(338) / '           0.209004959689     1.  4.  0.  0.' /
DATA Hydrogen_sulfide(339) / '          -0.00919656385346   1.  5.  0.  0.' /
DATA Hydrogen_sulfide(340) / '          -1.35782366339      3.  0.  2.  0.9873538' /
DATA Hydrogen_sulfide(341) / '           1.53224981014      4.  0.  2.  0.9873538' /
DATA Hydrogen_sulfide(342) / '          -0.329107661253     5.  0.  2.  0.9873538' /
DATA Hydrogen_sulfide(343) / '           0.891427552242     3.  2.  2.  0.9873538' /
DATA Hydrogen_sulfide(344) / '          -2.04776100441      4.  2.  2.  0.9873538' /
DATA Hydrogen_sulfide(345) / '           1.01366381241      5.  2.  2.  0.9873538' /
DATA Hydrogen_sulfide(346) / '' /
DATA Hydrogen_sulfide(347) / '' /
DATA Hydrogen_sulfide(348) / '@AUX    !---Auxiliary function for Cp0' /
DATA Hydrogen_sulfide(349) / 'CP2     !Ideal gas heat capacity function for hydrogen sulfide.' /
DATA Hydrogen_sulfide(350) / '          ?' /
DATA Hydrogen_sulfide(351) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(352) / '          ?Polt, A., Platzer, B., and Maurer, G.,' /
DATA Hydrogen_sulfide(353) / '          ?' /
DATA Hydrogen_sulfide(354) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(355) / '          0.                 !' /
DATA Hydrogen_sulfide(356) / '          10000.             !' /
DATA Hydrogen_sulfide(357) / '          0.                 !' /
DATA Hydrogen_sulfide(358) / '          0.                 !' /
DATA Hydrogen_sulfide(359) / '          1.0     8.3143     !Reducing parameters for T, Cp0' /
DATA Hydrogen_sulfide(360) / '          5 0   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Hydrogen_sulfide(361) / '           4.1012105         0.0' /
DATA Hydrogen_sulfide(362) / '          -0.0016720073      1.0' /
DATA Hydrogen_sulfide(363) / '           0.0000075303152   2.0' /
DATA Hydrogen_sulfide(364) / '          -0.62421053e-8     3.0' /
DATA Hydrogen_sulfide(365) / '           0.18098453e-11    4.0' /
DATA Hydrogen_sulfide(366) / '' /
DATA Hydrogen_sulfide(367) / '' /
DATA Hydrogen_sulfide(368) / '@EOS    !---Equation of state---' /
DATA Hydrogen_sulfide(369) / 'FE3     !Helmholtz equation of state for hydrogen sulfide of Starling (1973).' /
DATA Hydrogen_sulfide(370) / '          ?' /
DATA Hydrogen_sulfide(371) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(372) / '          ?Starling, K.E.,' /
DATA Hydrogen_sulfide(373) / '          ? "Fluid Thermodynamic Properties for Light Petroleum Systems,"' /
DATA Hydrogen_sulfide(374) / '          ? Gulf Publishing Company, 1973.' /
DATA Hydrogen_sulfide(375) / '          ?' /
DATA Hydrogen_sulfide(376) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(377) / '          187.7              !Lower temperature limit [K]' /
DATA Hydrogen_sulfide(378) / '          589.0              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(379) / '          55000.0            !Upper pressure limit [kPa]' /
DATA Hydrogen_sulfide(380) / '          29.578             !Maximum density [mol/L]' /
DATA Hydrogen_sulfide(381) / '          CP3                                    !Pointer to Cp0 model' /
DATA Hydrogen_sulfide(382) / '          34.08                                  !Molar mass [g/mol]' /
DATA Hydrogen_sulfide(383) / '          187.7                                  !Triple point temperature [K]' /
DATA Hydrogen_sulfide(384) / '          23.85                                  !Pressure at triple point [kPa]' /
DATA Hydrogen_sulfide(385) / '          29.07                                  !Density at triple point [mol/L]' /
DATA Hydrogen_sulfide(386) / '          213.142                                !Normal boiling point temperature [K]' /
DATA Hydrogen_sulfide(387) / '          0.0956                                 !Acentric factor' /
DATA Hydrogen_sulfide(388) / '          373.6         9008.0      10.16725352  !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Hydrogen_sulfide(389) / '          373.6                     10.16725352  !Reducing parameters [K, mol/L]' /
DATA Hydrogen_sulfide(390) / '          8.3159524                              !Gas constant [J/mol-K]' /
DATA Hydrogen_sulfide(391) / '            13  5    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Hydrogen_sulfide(392) / '           1.10928333109      3.  0.  0.  0.     !a(i),t(i),d(i),l(i)' /
DATA Hydrogen_sulfide(393) / '           0.188834546108     0.  1.  0.  0.' /
DATA Hydrogen_sulfide(394) / '          -0.930906931583     1.  1.  0.  0.' /
DATA Hydrogen_sulfide(395) / '          -0.411249591635     3.  1.  0.  0.' /
DATA Hydrogen_sulfide(396) / '           0.0140676923412    4.  1.  0.  0.' /
DATA Hydrogen_sulfide(397) / '          -0.169077883177e-4  5.  1.  0.  0.' /
DATA Hydrogen_sulfide(398) / '           0.510265859853     0.  2.  0.  0.' /
DATA Hydrogen_sulfide(399) / '          -0.572402742986     1.  2.  0.  0.' /
DATA Hydrogen_sulfide(400) / '          -0.000828859606622  2.  2.  0.  0.' /
DATA Hydrogen_sulfide(401) / '           0.00971664064871   1.  5.  0.  0.' /
DATA Hydrogen_sulfide(402) / '           0.140700425434e-4  2.  5.  0.  0.' /
DATA Hydrogen_sulfide(403) / '          -1.10928333109      3.  0.  2.  0.48524558' /
DATA Hydrogen_sulfide(404) / '          -0.26913741657      3.  2.  2.  0.48524558' /
DATA Hydrogen_sulfide(405) / '' /
DATA Hydrogen_sulfide(406) / '' /
DATA Hydrogen_sulfide(407) / '@AUX    !---Auxiliary function for Cp0' /
DATA Hydrogen_sulfide(408) / 'CP3     !Ideal gas heat capacity function for hydrogen sulfide.' /
DATA Hydrogen_sulfide(409) / '          ?' /
DATA Hydrogen_sulfide(410) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(411) / '          ?Starling, K.E.,' /
DATA Hydrogen_sulfide(412) / '          ?' /
DATA Hydrogen_sulfide(413) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(414) / '          0.                 !' /
DATA Hydrogen_sulfide(415) / '          10000.             !' /
DATA Hydrogen_sulfide(416) / '          0.                 !' /
DATA Hydrogen_sulfide(417) / '          0.                 !' /
DATA Hydrogen_sulfide(418) / '          1.0     4.184      !Reducing parameters for T, Cp0' /
DATA Hydrogen_sulfide(419) / '          1 0   2 2   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Hydrogen_sulfide(420) / '           7.9468     0.0' /
DATA Hydrogen_sulfide(421) / '           2032994.7     -2.0    843.792   -1.0  -2.0' /
DATA Hydrogen_sulfide(422) / '          -3504495.7     -2.0    1102.23   -1.0  -2.0' /
DATA Hydrogen_sulfide(423) / '          -15769.761     -2.0    433.801   -1.0  -2.0' /
DATA Hydrogen_sulfide(424) / '           13861204.0    -2.0    1481.43   -1.0  -2.0' /
DATA Hydrogen_sulfide(425) / '' /
DATA Hydrogen_sulfide(426) / '' /
DATA Hydrogen_sulfide(427) / '' /
DATA Hydrogen_sulfide(428) / '' /
DATA Hydrogen_sulfide(429) / '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' /
DATA Hydrogen_sulfide(430) / '' /
DATA Hydrogen_sulfide(431) / '#ETA   !---Viscosity---' /
DATA Hydrogen_sulfide(432) / 'VS4    !Friction theory viscosity model for hydrogen sulfide of Schmidt (2008).' /
DATA Hydrogen_sulfide(433) / ':DOI: 10.1021/ef700701h' /
DATA Hydrogen_sulfide(434) / '?' /
DATA Hydrogen_sulfide(435) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(436) / '?Schmidt, K.A.G., Carroll, J.J., Quinones-Cisneros, S.E., and Kvamme, B.,' /
DATA Hydrogen_sulfide(437) / '? "Hydrogen Sulphide Viscosity Model," proceedings of the' /
DATA Hydrogen_sulfide(438) / '? 86th Annual GPA Convention, March 11-14, San Antonio, TX, 2007.' /
DATA Hydrogen_sulfide(439) / '? See also: Schmidt, K.A.G., Quinones-Cisneros, S.E., Carroll, J.J., and Kvamme, B.,' /
DATA Hydrogen_sulfide(440) / '? "Hydrogen Sulfide Viscosity Modeling," Energy & Fuels, 22, 3424-3434.' /
DATA Hydrogen_sulfide(441) / '?' /
DATA Hydrogen_sulfide(442) / '?The correlation agrees with available experimental data with' /
DATA Hydrogen_sulfide(443) / '? an average absolute percent deviation of 1% over the temperature' /
DATA Hydrogen_sulfide(444) / '? range 190-600 K at atmospheric pressure, and along the liquid saturation' /
DATA Hydrogen_sulfide(445) / '? boundary. At pressures of 100 MPa the uncertainty is estimated' /
DATA Hydrogen_sulfide(446) / '? to be on the order of 10%.' /
DATA Hydrogen_sulfide(447) / '?' /
DATA Hydrogen_sulfide(448) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(449) / '187.7              !Lower temperature limit [K]' /
DATA Hydrogen_sulfide(450) / '760.0              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(451) / '170000.0           !Upper pressure limit [kPa]' /
DATA Hydrogen_sulfide(452) / '29.12              !Maximum density [mol/L]' /
DATA Hydrogen_sulfide(453) / '4 0 0 0 0 0        !Number of terms associated with dilute-gas function' /
DATA Hydrogen_sulfide(454) / 'NUL                !Pointer to reduced effective collision cross-section model; not used' /
DATA Hydrogen_sulfide(455) / '0.36237            !Lennard-Jones coefficient sigma [nm] (not used)' /
DATA Hydrogen_sulfide(456) / '301.1              !Lennard-Jones coefficient epsilon/kappa [K] (not used)' /
DATA Hydrogen_sulfide(457) / ' 373.1     1.0     !Reducing parameters for T, eta' /
DATA Hydrogen_sulfide(458) / ' 0.0       0.5     !Chapman-Enskog term; not used here' /
DATA Hydrogen_sulfide(459) / ' 43.6694   0.0     !Empirical terms for eta0' /
DATA Hydrogen_sulfide(460) / '-121.530   0.25' /
DATA Hydrogen_sulfide(461) / ' 93.5279   0.5' /
DATA Hydrogen_sulfide(462) / '0                  !Number of terms for initial density dependence' /
DATA Hydrogen_sulfide(463) / ' 5.46919e-5  -7.32295e-6  -7.35622e-6  0.  0.      !  a(0),a(1),a(2)' /
DATA Hydrogen_sulfide(464) / ' 4.56159e-5  -1.82572e-5  -6.59654e-6  0.  0.      !  b(0),b(1),b(2)' /
DATA Hydrogen_sulfide(465) / '-4.33882e-6   6.13716e-6   0.0         0.  0.      !  c(0),c(1),c(2)' /
DATA Hydrogen_sulfide(466) / ' 6.67324e-9  -2.16365e-9   0.0         0.  0.      !  A(0),A(1),A(2)' /
DATA Hydrogen_sulfide(467) / '-1.53973e-9   2.17652e-9   0.0         0.  0.      !  B(0),B(1),B(2)' /
DATA Hydrogen_sulfide(468) / ' 3.54228e-7  -4.76258e-8   0.0         0.  0.      !  C(0),C(1),C(2)' /
DATA Hydrogen_sulfide(469) / ' 0.0          0.0          0.0         0.  0.      !  D(0),D(1),D(2)' /
DATA Hydrogen_sulfide(470) / ' 0.0          0.0          0.0         0.  0.      !  E(0),E(1),E(2)' /
DATA Hydrogen_sulfide(471) / 'NUL                !Pointer to the viscosity critical enhancement auxiliary function (none used)' /
DATA Hydrogen_sulfide(472) / '' /
DATA Hydrogen_sulfide(473) / '' /
DATA Hydrogen_sulfide(474) / '' /
DATA Hydrogen_sulfide(475) / '' /
DATA Hydrogen_sulfide(476) / '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' /
DATA Hydrogen_sulfide(477) / '@TRN    !---ECS Transport---' /
DATA Hydrogen_sulfide(478) / 'ECS     !Extended Corresponding States model (Propane reference) PREDICTIVE MODEL' /
DATA Hydrogen_sulfide(479) / ':DOI: 10.6028/NIST.IR.8209' /
DATA Hydrogen_sulfide(480) / '          ?' /
DATA Hydrogen_sulfide(481) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(482) / '          ?Huber, M.L., "Models for the Viscosity, Thermal Conductivity, and Surface Tension' /
DATA Hydrogen_sulfide(483) / '          ? of Selected Pure Fluids as Implemented in REFPROP v10.0," NISTIR 8209, 2018.' /
DATA Hydrogen_sulfide(484) / '          ? doi: 10.6028/NIST.IR.8209' /
DATA Hydrogen_sulfide(485) / '          ?' /
DATA Hydrogen_sulfide(486) / '          ?VISCOSITY' /
DATA Hydrogen_sulfide(487) / '          ? Model not fit.' /
DATA Hydrogen_sulfide(488) / '          ? The estimated uncertainty of the thermal conductivity of the liquid phase and gas phases is 20%, larger near critical.' /
DATA Hydrogen_sulfide(489) / '          ?' /
DATA Hydrogen_sulfide(490) / '          ?THERMAL CONDUCTIVITY' /
DATA Hydrogen_sulfide(491) / '          ? Predictive values only- NO EXPERIMENTAL DATA FOUND FOR LIQUID PHASE' /
DATA Hydrogen_sulfide(492) / '          ? The estimated uncertainty of the thermal conductivity of the liquid phase is 50%.' /
DATA Hydrogen_sulfide(493) / '          ? Estimated uncertainty of the gas phase is 2% for T< 1000 K based on comparisons with ab-initio values of Hellman et al., J. Chem. Eng. Data, 1312-1317, 2012.' /
DATA Hydrogen_sulfide(494) / '          ?' /
DATA Hydrogen_sulfide(495) / '          ?The Lennard-Jones parameters were estimated with the method of Chung, T.H., Ajlan, M., Lee, L.L., and Starling, K.E., "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties," Ind. Eng. Chem. Res., 27:671-679, 1988.' /
DATA Hydrogen_sulfide(496) / '          ?' /
DATA Hydrogen_sulfide(497) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(498) / '          187.7              !Lower temperature limit [K]' /
DATA Hydrogen_sulfide(499) / '          760.0              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(500) / '          170000.0           !Upper pressure limit [kPa]' /
DATA Hydrogen_sulfide(501) / '          29.12              !Maximum density [mol/L]' /
DATA Hydrogen_sulfide(502) / '          FEQ PROPANE.FLD' /
DATA Hydrogen_sulfide(503) / '          VS1                !Model for reference fluid viscosity' /
DATA Hydrogen_sulfide(504) / '          TC1                !Model for reference fluid thermal conductivity' /
DATA Hydrogen_sulfide(505) / '          BIG                !Large molecule identifier' /
DATA Hydrogen_sulfide(506) / '          1.00 0. 0. 0.      !Large molecule parameters' /
DATA Hydrogen_sulfide(507) / '          1                  !Lennard-Jones flag (0 or 1) (0 => use estimates)' /
DATA Hydrogen_sulfide(508) / '          0.3732             !Lennard-Jones coefficient sigma [nm] for ECS method' /
DATA Hydrogen_sulfide(509) / '          296.28             !Lennard-Jones coefficient epsilon/kappa [K] for ECS method' /
DATA Hydrogen_sulfide(510) / '          3  0  0                  !Number of terms in f_int term in Eucken correlation, spare1, spare2' /
DATA Hydrogen_sulfide(511) / '           1.5603e-04       0. 0. 0.  !Coefficient, power of T, spare1, spare2' /
DATA Hydrogen_sulfide(512) / '           1.78874e-6       1. 0. 0.  !Coefficient, power of T, spare1, spare2' /
DATA Hydrogen_sulfide(513) / '           -6.75136e-10     2. 0. 0.  !Coefficient, power of T, spare1, spare2' /
DATA Hydrogen_sulfide(514) / '          2  0  0                  !Number of terms in psi (visc shape factor): poly,spare1,spare2' /
DATA Hydrogen_sulfide(515) / '           1.0           0. 0. 0.  !Coefficient, power of T, spare1, spare2' /
DATA Hydrogen_sulfide(516) / '           0.0           0. 1. 0.  !Coefficient, power of T, spare1, spare2' /
DATA Hydrogen_sulfide(517) / '          2  0  0                  !Number of terms in chi (t.c. shape factor): poly,spare1,spare2' /
DATA Hydrogen_sulfide(518) / '            1.0          0. 0. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Hydrogen_sulfide(519) / '            0.0          0. 1. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Hydrogen_sulfide(520) / '          TK3                !Pointer to critical enhancement auxiliary function' /
DATA Hydrogen_sulfide(521) / '' /
DATA Hydrogen_sulfide(522) / '' /
DATA Hydrogen_sulfide(523) / '#AUX   !---Auxiliary function for the thermal conductivity critical enhancement' /
DATA Hydrogen_sulfide(524) / 'TK3    !Simplified thermal conductivity critical enhancement for hydrogen sulfide of Perkins et al. (2013).' /
DATA Hydrogen_sulfide(525) / '?' /
DATA Hydrogen_sulfide(526) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(527) / '?Perkins, R.A., Sengers, J.V., Abdulagatov, I.M., and Huber, M.L.,' /
DATA Hydrogen_sulfide(528) / '? "Simplified Model for the Critical Thermal-Conductivity Enhancement in Molecular Fluids,"' /
DATA Hydrogen_sulfide(529) / '? Int. J. Thermophys., 34(2):191-212, 2013. doi: 10.1007/s10765-013-1409-z' /
DATA Hydrogen_sulfide(530) / '?' /
DATA Hydrogen_sulfide(531) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(532) / '0.                 !' /
DATA Hydrogen_sulfide(533) / '10000.             !' /
DATA Hydrogen_sulfide(534) / '0.                 !' /
DATA Hydrogen_sulfide(535) / '0.                 !' /
DATA Hydrogen_sulfide(536) / '9 0 0 0            !# terms:  terms, spare, spare, spare' /
DATA Hydrogen_sulfide(537) / '1.0  1.0  1.0      !Reducing parameters for T, rho, tcx [mW/(m-K)]' /
DATA Hydrogen_sulfide(538) / '0.63               !Nu (universal exponent)' /
DATA Hydrogen_sulfide(539) / '1.239              !Gamma (universal exponent)' /
DATA Hydrogen_sulfide(540) / '1.02               !R0 (universal amplitude)' /
DATA Hydrogen_sulfide(541) / '0.063              !Z (universal exponent--not used for t.c., only viscosity)' /
DATA Hydrogen_sulfide(542) / '1.0                !C (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA Hydrogen_sulfide(543) / '0.164e-9           !Xi0 (amplitude) [m]' /
DATA Hydrogen_sulfide(544) / '0.058              !Gam0 (amplitude) [-]' /
DATA Hydrogen_sulfide(545) / '0.447e-9           !Qd_inverse (modified effective cutoff parameter) [m]' /
DATA Hydrogen_sulfide(546) / '559.65             !Tref (reference temperature) [K]' /
DATA Hydrogen_sulfide(547) / '' /
DATA Hydrogen_sulfide(548) / '' /
DATA Hydrogen_sulfide(549) / '' /
DATA Hydrogen_sulfide(550) / '' /
DATA Hydrogen_sulfide(551) / '================================================================================' /
DATA Hydrogen_sulfide(552) / '' /
DATA Hydrogen_sulfide(553) / '@TCX    !---Thermal conductivity---' /
DATA Hydrogen_sulfide(554) / 'TC1     !Pure fluid thermal conductivity model from NIST14 for hydrogen sulfide.' /
DATA Hydrogen_sulfide(555) / '          ?' /
DATA Hydrogen_sulfide(556) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(557) / '          ?Dense fluid coefficients are taken from NIST14, Version 9.08' /
DATA Hydrogen_sulfide(558) / '          ? Dilute gas refit with data from DIPPR diadem Aug 2008' /
DATA Hydrogen_sulfide(559) / '          ?' /
DATA Hydrogen_sulfide(560) / '          ?Critical enhancement model of Olchowy and Sengers added. Estimated uncertainty,' /
DATA Hydrogen_sulfide(561) / '          ? except near the critical region, is 4-6%' /
DATA Hydrogen_sulfide(562) / '          ?' /
DATA Hydrogen_sulfide(563) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(564) / '          187.7              !Lower temperature limit [K]' /
DATA Hydrogen_sulfide(565) / '          760.0              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(566) / '          170000.0           !Upper pressure limit [kPa]' /
DATA Hydrogen_sulfide(567) / '          29.13              !Maximum density [mol/L]' /
DATA Hydrogen_sulfide(568) / '          3   0              !# terms for dilute gas function:  numerator, denominator' /
DATA Hydrogen_sulfide(569) / '           1.0         1.    !Reducing parameters for T, tcx' /
DATA Hydrogen_sulfide(570) / '          -0.008415471 0.    !Coefficient, power in T' /
DATA Hydrogen_sulfide(571) / '           8.264077e-5 1.    !Coefficient, power in T' /
DATA Hydrogen_sulfide(572) / '          -2.526101e-8 2.    !Coefficient, power in T' /
DATA Hydrogen_sulfide(573) / '          6   0              !# terms for background gas function:  numerator, denominator' /
DATA Hydrogen_sulfide(574) / '           373.4       10.2      0.001           !Reducing parameters for T, rho, tcx' /
DATA Hydrogen_sulfide(575) / '           21.7827447865  0. 1. 0. !Coefficient, powers of T, rho, exp(rho)' /
DATA Hydrogen_sulfide(576) / '           10.8880930411  0. 3. 0.' /
DATA Hydrogen_sulfide(577) / '          -7.45794247629  0. 4. 0.' /
DATA Hydrogen_sulfide(578) / '           3.65609005216 -1. 4. 0.' /
DATA Hydrogen_sulfide(579) / '           1.89362258187  0. 5. 0.' /
DATA Hydrogen_sulfide(580) / '          -1.10975687736 -1. 5. 0.' /
DATA Hydrogen_sulfide(581) / '          TK3                !Pointer to critical enhancement auxiliary function' /
DATA Hydrogen_sulfide(582) / '' /
DATA Hydrogen_sulfide(583) / '' /
DATA Hydrogen_sulfide(584) / '' /
DATA Hydrogen_sulfide(585) / '' /
DATA Hydrogen_sulfide(586) / '********************************************************************************' /
DATA Hydrogen_sulfide(587) / '' /
DATA Hydrogen_sulfide(588) / '@ETA    !---Viscosity---' /
DATA Hydrogen_sulfide(589) / 'VS2     !Pure fluid viscosity model from NIST14 for hydrogen sulfide.' /
DATA Hydrogen_sulfide(590) / '          ?' /
DATA Hydrogen_sulfide(591) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(592) / '          ?Coefficients are taken from NIST14, Version 9.08' /
DATA Hydrogen_sulfide(593) / '          ?' /
DATA Hydrogen_sulfide(594) / '          ?Estimated uncertainty is 2 %.' /
DATA Hydrogen_sulfide(595) / '          ?' /
DATA Hydrogen_sulfide(596) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(597) / '          187.7              !Lower temperature limit [K]' /
DATA Hydrogen_sulfide(598) / '          760.0              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(599) / '          170000.0           !Upper pressure limit [kPa]' /
DATA Hydrogen_sulfide(600) / '          29.13              !Maximum density [mol/L]' /
DATA Hydrogen_sulfide(601) / '          CI0                !Pointer to collision integral model' /
DATA Hydrogen_sulfide(602) / '          0.36237            !Lennard-Jones coefficient sigma [nm]' /
DATA Hydrogen_sulfide(603) / '          301.1              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA Hydrogen_sulfide(604) / '           0.1558117         !Const' /
DATA Hydrogen_sulfide(605) / '           0.5               !Exponent for T' /
DATA Hydrogen_sulfide(606) / '           0.0               !Coefficient for initial density dependence of viscosity' /
DATA Hydrogen_sulfide(607) / '           0.0' /
DATA Hydrogen_sulfide(608) / '           0.0' /
DATA Hydrogen_sulfide(609) / '           100.0' /
DATA Hydrogen_sulfide(610) / '          -12.3286304189940  !Coefficients for residual viscosity' /
DATA Hydrogen_sulfide(611) / '           782.29421491' /
DATA Hydrogen_sulfide(612) / '           11.840322553' /
DATA Hydrogen_sulfide(613) / '          -10401.582791' /
DATA Hydrogen_sulfide(614) / '          -0.0482407464' /
DATA Hydrogen_sulfide(615) / '           69.709031672' /
DATA Hydrogen_sulfide(616) / '           256.31792390' /
DATA Hydrogen_sulfide(617) / '           10.2' /
DATA Hydrogen_sulfide(618) / '          NUL                !Pointer to the viscosity critical enhancement auxiliary function (none used)' /
DATA Hydrogen_sulfide(619) / '' /
DATA Hydrogen_sulfide(620) / '' /
DATA Hydrogen_sulfide(621) / '' /
DATA Hydrogen_sulfide(622) / '' /
DATA Hydrogen_sulfide(623) / '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' /
DATA Hydrogen_sulfide(624) / '' /
DATA Hydrogen_sulfide(625) / '#DE    !---Dielectric constant---' /
DATA Hydrogen_sulfide(626) / 'DE5    !Dielectric constant model for H2S of Harvey and Mountain (2017).' /
DATA Hydrogen_sulfide(627) / ':DOI: 10.1007/s10765-017-2279-6' /
DATA Hydrogen_sulfide(628) / '?' /
DATA Hydrogen_sulfide(629) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(630) / '?Harvey, A.H. and Mountain, R.D.,' /
DATA Hydrogen_sulfide(631) / '? "Correlations for the Dielectric Constants of H2S, SO2, and SF6,"' /
DATA Hydrogen_sulfide(632) / '? Int. J. Thermophys., 38:147, 2017.' /
DATA Hydrogen_sulfide(633) / '?' /
DATA Hydrogen_sulfide(634) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(635) / '187.7              !Lower temperature limit [K]' /
DATA Hydrogen_sulfide(636) / '760.0              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(637) / '0.                 !' /
DATA Hydrogen_sulfide(638) / '0.                 !' /
DATA Hydrogen_sulfide(639) / '380.0 29.0 1.0     !Reducing parameters for T and D' /
DATA Hydrogen_sulfide(640) / '4 1 1 1 0 0        !Number of terms in dielectric constant model' /
DATA Hydrogen_sulfide(641) / ' 3.66e-24  0.0   0.0    0.    !  alpha (cm^3)' /
DATA Hydrogen_sulfide(642) / ' 0.978325  0.0   0.0    0.    !  mu (debye)' /
DATA Hydrogen_sulfide(643) / ' 0.241     0.0   0.0    0.    !  cu' /
DATA Hydrogen_sulfide(644) / ' 1.18      0.0   0.0    0.    !  cg' /
DATA Hydrogen_sulfide(645) / ' 1.        0.0   2.83   0.5   !  f' /
DATA Hydrogen_sulfide(646) / ' 1.        0.0   2.546  0.9   !  g1' /
DATA Hydrogen_sulfide(647) / ' 1.        0.0   1.883  3.5   !  g2' /
DATA Hydrogen_sulfide(648) / '' /
DATA Hydrogen_sulfide(649) / '' /
DATA Hydrogen_sulfide(650) / '#STN   !---Surface tension---' /
DATA Hydrogen_sulfide(651) / 'ST1    !Surface tension model for hydrogen sulfide of Mulero et al. (2012).' /
DATA Hydrogen_sulfide(652) / ':DOI: 10.1063/1.4768782' /
DATA Hydrogen_sulfide(653) / '?' /
DATA Hydrogen_sulfide(654) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(655) / '?Mulero, A., Cachadina, I., and Parra, M.I.,' /
DATA Hydrogen_sulfide(656) / '? "Recommended Correlations for the Surface Tension of Common Fluids,"' /
DATA Hydrogen_sulfide(657) / '? J. Phys. Chem. Ref. Data, 41(4), 043105, 2012. doi: 10.1063/1.4768782' /
DATA Hydrogen_sulfide(658) / '?' /
DATA Hydrogen_sulfide(659) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(660) / '0.                 !' /
DATA Hydrogen_sulfide(661) / '10000.             !' /
DATA Hydrogen_sulfide(662) / '0.                 !' /
DATA Hydrogen_sulfide(663) / '0.                 !' /
DATA Hydrogen_sulfide(664) / '1                  !Number of terms in surface tension model' /
DATA Hydrogen_sulfide(665) / '373.1              !Critical temperature used in fit (dummy)' /
DATA Hydrogen_sulfide(666) / '0.078557  1.2074   !Sigma0 and n' /
DATA Hydrogen_sulfide(667) / '' /
DATA Hydrogen_sulfide(668) / '' /
DATA Hydrogen_sulfide(669) / '#SBL   !---Sublimation line---' /
DATA Hydrogen_sulfide(670) / 'SB2    !Sublimation line model for hydrogen sulfide of Fray and Schmitt (2009).' /
DATA Hydrogen_sulfide(671) / ':DOI: 10.1016/j.pss.2009.09.011' /
DATA Hydrogen_sulfide(672) / '?' /
DATA Hydrogen_sulfide(673) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(674) / '? Based on N. Fray and B. Schmitt, Planet. Space Sci. 57, 2053-2080, 2009.' /
DATA Hydrogen_sulfide(675) / '? Modified to match the triple point of the equation of state.' /
DATA Hydrogen_sulfide(676) / '?' /
DATA Hydrogen_sulfide(677) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(678) / '0.                 !' /
DATA Hydrogen_sulfide(679) / '187.7              !Upper temperature limit [K]' /
DATA Hydrogen_sulfide(680) / '0.                 !' /
DATA Hydrogen_sulfide(681) / '0.                 !' /
DATA Hydrogen_sulfide(682) / '1.0  1000.0        !Reducing temperature and pressure' /
DATA Hydrogen_sulfide(683) / '5 0 0 0 0 0        !Number of terms in sublimation line equation' /
DATA Hydrogen_sulfide(684) / ' 6.6247   0.0      !Coefficients and exponents' /
DATA Hydrogen_sulfide(685) / '-7.260e2 -1.0' /
DATA Hydrogen_sulfide(686) / '-3.504e5 -2.0' /
DATA Hydrogen_sulfide(687) / ' 2.724e7 -3.0' /
DATA Hydrogen_sulfide(688) / '-8.582e8 -4.0' /
DATA Hydrogen_sulfide(689) / '' /
DATA Hydrogen_sulfide(690) / '' /
DATA Hydrogen_sulfide(691) / '#PS    !---Vapor pressure---' /
DATA Hydrogen_sulfide(692) / 'PS5    !Vapor pressure equation for hydrogen sulfide of Lemmon (2006).' /
DATA Hydrogen_sulfide(693) / '?' /
DATA Hydrogen_sulfide(694) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(695) / '?Functional Form:  P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T] where Theta=1-T/Tc, Tc and Pc' /
DATA Hydrogen_sulfide(696) / '? are the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Hydrogen_sulfide(697) / '?' /
DATA Hydrogen_sulfide(698) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(699) / '0.                 !' /
DATA Hydrogen_sulfide(700) / '10000.             !' /
DATA Hydrogen_sulfide(701) / '0.                 !' /
DATA Hydrogen_sulfide(702) / '0.                 !' /
DATA Hydrogen_sulfide(703) / '373.1     9000.0   !Reducing parameters' /
DATA Hydrogen_sulfide(704) / '4 0 0 0 0 0        !Number of terms in equation' /
DATA Hydrogen_sulfide(705) / '-6.5884    1.0' /
DATA Hydrogen_sulfide(706) / ' 2.1582    1.5' /
DATA Hydrogen_sulfide(707) / '-1.6054    2.0' /
DATA Hydrogen_sulfide(708) / '-2.3870    4.8' /
DATA Hydrogen_sulfide(709) / '' /
DATA Hydrogen_sulfide(710) / '' /
DATA Hydrogen_sulfide(711) / '#DL    !---Saturated liquid density---' /
DATA Hydrogen_sulfide(712) / 'DL1    !Saturated liquid density equation for hydrogen sulfide of Lemmon (2006).' /
DATA Hydrogen_sulfide(713) / '?' /
DATA Hydrogen_sulfide(714) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(715) / '?Functional Form:  D=Dc*[1+SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Hydrogen_sulfide(716) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Hydrogen_sulfide(717) / '?' /
DATA Hydrogen_sulfide(718) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(719) / '0.                 !' /
DATA Hydrogen_sulfide(720) / '10000.             !' /
DATA Hydrogen_sulfide(721) / '0.                 !' /
DATA Hydrogen_sulfide(722) / '0.                 !' /
DATA Hydrogen_sulfide(723) / '373.1       10.19  !Reducing parameters' /
DATA Hydrogen_sulfide(724) / '5 0 0 0 0 0        !Number of terms in equation' /
DATA Hydrogen_sulfide(725) / ' 9.28051   0.578' /
DATA Hydrogen_sulfide(726) / '-21.7006   0.85' /
DATA Hydrogen_sulfide(727) / ' 30.1375   1.15' /
DATA Hydrogen_sulfide(728) / '-22.1012   1.5' /
DATA Hydrogen_sulfide(729) / ' 7.1598    1.9' /
DATA Hydrogen_sulfide(730) / '' /
DATA Hydrogen_sulfide(731) / '' /
DATA Hydrogen_sulfide(732) / '#DV    !---Saturated vapor density---' /
DATA Hydrogen_sulfide(733) / 'DV3    !Saturated vapor density equation for hydrogen sulfide of Lemmon (2006).' /
DATA Hydrogen_sulfide(734) / '?' /
DATA Hydrogen_sulfide(735) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(736) / '?Functional Form:  D=Dc*EXP[SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Hydrogen_sulfide(737) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Hydrogen_sulfide(738) / '?' /
DATA Hydrogen_sulfide(739) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_sulfide(740) / '0.                 !' /
DATA Hydrogen_sulfide(741) / '10000.             !' /
DATA Hydrogen_sulfide(742) / '0.                 !' /
DATA Hydrogen_sulfide(743) / '0.                 !' /
DATA Hydrogen_sulfide(744) / '373.1       10.19  !Reducing parameters' /
DATA Hydrogen_sulfide(745) / '6 0 0 0 0 0        !Number of terms in equation' /
DATA Hydrogen_sulfide(746) / '-6.49882   0.551' /
DATA Hydrogen_sulfide(747) / ' 7.98446   0.88' /
DATA Hydrogen_sulfide(748) / '-10.9661   1.23' /
DATA Hydrogen_sulfide(749) / '-16.6445   3.46' /
DATA Hydrogen_sulfide(750) / '-44.924    7.2' /
DATA Hydrogen_sulfide(751) / '-208.126   17.0' /
DATA Hydrogen_sulfide(752) / '' /
DATA Hydrogen_sulfide(753) / '' /
DATA Hydrogen_sulfide(754) / '@END' /
DATA Hydrogen_sulfide(755) / 'c        1         2         3         4         5         6         7         8' /
DATA Hydrogen_sulfide(756) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
! #######################################################
character(256), TARGET :: Hydrogen_chloride(392)
DATA Hydrogen_chloride(1) / 'Hydrogen chloride    !Short name' /
DATA Hydrogen_chloride(2) / '7647-01-0            !CAS number' /
DATA Hydrogen_chloride(3) / 'Hydrogen chloride    !Full name' /
DATA Hydrogen_chloride(4) / 'HCl                  !Chemical formula' /
DATA Hydrogen_chloride(5) / 'Hydrogen chloride    !Synonym' /
DATA Hydrogen_chloride(6) / '36.46094             !Molar mass [g/mol]' /
DATA Hydrogen_chloride(7) / '159.07               !Triple point temperature [K]' /
DATA Hydrogen_chloride(8) / '188.173              !Normal boiling point [K]' /
DATA Hydrogen_chloride(9) / '324.68               !Critical temperature [K]' /
DATA Hydrogen_chloride(10) / '8313.5               !Critical pressure [kPa]' /
DATA Hydrogen_chloride(11) / '11.87                !Critical density [mol/L]' /
DATA Hydrogen_chloride(12) / '0.129                !Acentric factor' /
DATA Hydrogen_chloride(13) / '1.079                !Dipole moment [Debye]; DIPPR DIADEM 2012' /
DATA Hydrogen_chloride(14) / 'NBP                  !Default reference state' /
DATA Hydrogen_chloride(15) / '10.0                 !Version number' /
DATA Hydrogen_chloride(16) / '1789                 !UN Number                                                 :UN:' /
DATA Hydrogen_chloride(17) / 'other                !Family                                                    :Family:' /
DATA Hydrogen_chloride(18) / '????                 !Heating value (upper) [kJ/mol]                            :Heat:' /
DATA Hydrogen_chloride(19) / '1S/ClH/h1H                                !Standard InChI String                :InChi:' /
DATA Hydrogen_chloride(20) / 'VEXZGXHMUGYJMC-UHFFFAOYSA-N               !Standard InChI Key                   :InChiKey:' /
DATA Hydrogen_chloride(21) / '????                                      !Alternative fluid for mixing rules   :AltID:' /
DATA Hydrogen_chloride(22) / '74b17450                                  !Hash number from InChI Key           :Hash:' /
DATA Hydrogen_chloride(23) / '' /
DATA Hydrogen_chloride(24) / '' /
DATA Hydrogen_chloride(25) / '' /
DATA Hydrogen_chloride(26) / '' /
DATA Hydrogen_chloride(27) / '!The fluid files contain general information about the fluid in the first 15 to 20 lines, followed by sections for the' /
DATA Hydrogen_chloride(28) / '! equations of state, transport equations, and auxiliary equations.  Equations of state are listed first.  The NIST recommended' /
DATA Hydrogen_chloride(29) / '! equations begin with a hash mark (#).  The secondary equations begin with the @ symbol.  These symbols can be swapped to' /
DATA Hydrogen_chloride(30) / '! select a secondary equation as primary and the primary as secondary.  The equation of state section also contains auxiliary' /
DATA Hydrogen_chloride(31) / '! equations for the ideal gas heat capacity or ideal gas Helmholtz energy.  Below the equations of state (both primary and' /
DATA Hydrogen_chloride(32) / '! secondary) are the transport equations, first viscosity and then thermal conductivity.  These are then followed by the' /
DATA Hydrogen_chloride(33) / '! secondary equations if available.  The transport section also contains auxiliary equations required to calculate either the' /
DATA Hydrogen_chloride(34) / '! dilute gas state or the critical enhancement.  At the end of the file are additional but not necessary auxiliary equations,' /
DATA Hydrogen_chloride(35) / '! including simple equations for the vapor pressure, saturated liquid and vapor densities, melting line (for some fluids), and' /
DATA Hydrogen_chloride(36) / '! sublimation line (for even fewer fluids).  This section also contains the equations for dielectric constant and surface' /
DATA Hydrogen_chloride(37) / '! tension if available.  The sections are divided by different symbols (these being _-+=^*~) to aid the eye in locating a' /
DATA Hydrogen_chloride(38) / '! particular section.  Secondary equations are indented 10 spaces to avoid confusion with the NIST recommended equations.  The' /
DATA Hydrogen_chloride(39) / '! end of the fluid file is marked with @END.  Anything below that is ignored.' /
DATA Hydrogen_chloride(40) / '' /
DATA Hydrogen_chloride(41) / '' /
DATA Hydrogen_chloride(42) / '! compiled by M. Thol, Thermodynamics, Ruhr-Universitaet Bochum, Germany' /
DATA Hydrogen_chloride(43) / '! 05-02-11  MT, Original version.' /
DATA Hydrogen_chloride(44) / '! 05-03-11  MT, Add ancillary equations.' /
DATA Hydrogen_chloride(45) / '! 04-06-13 EWL, Add dipole moment.' /
DATA Hydrogen_chloride(46) / '! 12-24-13 EWL, Add truncated coefficients from publication.' /
DATA Hydrogen_chloride(47) / '! 03-27-14 MLH, Add preliminary transport.' /
DATA Hydrogen_chloride(48) / '! 04-17-14 EWL, Add surface tension coefficients of Mulero et al. (2014).' /
DATA Hydrogen_chloride(49) / '! 02-19-15  MT, Add final equation of state.' /
DATA Hydrogen_chloride(50) / '! 02-16-17  KG, Add ancillary equations.' /
DATA Hydrogen_chloride(51) / '! 04-03-17 MLH, Revise transport.' /
DATA Hydrogen_chloride(52) / '! 07-31-17  MT, Add second final equation of state.' /
DATA Hydrogen_chloride(53) / '! 08-04-17 MLH, revise transport' /
DATA Hydrogen_chloride(54) / '! 12-30-17 MLH, tweak in enhancement and LJ parameters' /
DATA Hydrogen_chloride(55) / '' /
DATA Hydrogen_chloride(56) / '' /
DATA Hydrogen_chloride(57) / '' /
DATA Hydrogen_chloride(58) / '' /
DATA Hydrogen_chloride(59) / '________________________________________________________________________________' /
DATA Hydrogen_chloride(60) / '' /
DATA Hydrogen_chloride(61) / '#EOS   !---Equation of state---' /
DATA Hydrogen_chloride(62) / 'FEQ    !Helmholtz equation of state for hydrogen chloride of Thol et al. (2018).' /
DATA Hydrogen_chloride(63) / ':TRUECRITICALPOINT:  324.68    11.87          !True EOS critical point [K, mol/L] (where dP/dD=0 and d^2P/dD^2=0 at constant T)' /
DATA Hydrogen_chloride(64) / ':DOI: 10.1021/acs.jced.7b0103' /
DATA Hydrogen_chloride(65) / '?' /
DATA Hydrogen_chloride(66) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(67) / '?Thol, M., Dubberke, F.H., Baumh?gger, E., Span, R., and Vrabec, J.,' /
DATA Hydrogen_chloride(68) / '? to be submitted to J. Chem. Eng. Data, 2018.' /
DATA Hydrogen_chloride(69) / '?' /
DATA Hydrogen_chloride(70) / '?Based on the available experimental data, the equation is valid from the triple' /
DATA Hydrogen_chloride(71) / '? point temperature of 159.07 K to 480 K and up to a maximum pressure of 40 MPa.' /
DATA Hydrogen_chloride(72) / '? This range can be extended to 670 K and 200 MPa with higher uncertainties based' /
DATA Hydrogen_chloride(73) / '? on the data of Franck et al.  The uncertainties in density are 1.5% in the' /
DATA Hydrogen_chloride(74) / '? gaseous region, 0.5% in the liquid region, and 1% in the supercritical region.' /
DATA Hydrogen_chloride(75) / '? At higher temperatures, pressures, and densities, where only the data of Franck' /
DATA Hydrogen_chloride(76) / '? et al. are available, the uncertainty increases to at least 6%.  The uncertainty' /
DATA Hydrogen_chloride(77) / '? of the second virial coefficient calculated with the present equation of state' /
DATA Hydrogen_chloride(78) / '? is estimated to be 25 cm^3/mol for T < 300 K and 15 cm^3/mol for higher' /
DATA Hydrogen_chloride(79) / '? temperatures.  The uncertainties in vapor pressure are 0.5% for T < 250 K and 1%' /
DATA Hydrogen_chloride(80) / '? for higher temperatures.  Based on only limited information, the uncertainty in' /
DATA Hydrogen_chloride(81) / '? the saturated liquid density is 0.5% for T < 240 K and 1% for higher' /
DATA Hydrogen_chloride(82) / '? temperatures.  For the saturated vapor density, no reliable measurements are' /
DATA Hydrogen_chloride(83) / '? available.  The uncertainty is 0.3% for the speed of sound except for the' /
DATA Hydrogen_chloride(84) / '? low-density region where the uncertainty increases up to 1%.' /
DATA Hydrogen_chloride(85) / '?' /
DATA Hydrogen_chloride(86) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(87) / '159.07             !Lower temperature limit [K]' /
DATA Hydrogen_chloride(88) / '670.0              !Upper temperature limit [K]' /
DATA Hydrogen_chloride(89) / '200000.0           !Upper pressure limit [kPa]' /
DATA Hydrogen_chloride(90) / '34.5               !Maximum density [mol/L]' /
DATA Hydrogen_chloride(91) / 'CPP                                    !Pointer to Cp0 model' /
DATA Hydrogen_chloride(92) / '36.46094                               !Molar mass [g/mol]' /
DATA Hydrogen_chloride(93) / '159.07                                 !Triple point temperature [K]' /
DATA Hydrogen_chloride(94) / '13.8284                                !Pressure at triple point [kPa]' /
DATA Hydrogen_chloride(95) / '34.401                                 !Density at triple point [mol/L]' /
DATA Hydrogen_chloride(96) / '188.173                                !Normal boiling point temperature [K]' /
DATA Hydrogen_chloride(97) / '0.129                                  !Acentric factor' /
DATA Hydrogen_chloride(98) / '324.68        8313.5      11.87        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Hydrogen_chloride(99) / '324.68                    11.87        !Reducing parameters [K, mol/L]' /
DATA Hydrogen_chloride(100) / '8.3144598                              !Gas constant [J/mol-K]' /
DATA Hydrogen_chloride(101) / '  10  4   5 12   0 0    0 0  0 0  0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Hydrogen_chloride(102) / '  0.01952802   1.0     4.  0.          !a(i),t(i),d(i),l(i)' /
DATA Hydrogen_chloride(103) / '  1.926809     0.553   1.  0.' /
DATA Hydrogen_chloride(104) / ' -2.835744     1.037   1.  0.' /
DATA Hydrogen_chloride(105) / ' -0.2276121    0.817   2.  0.' /
DATA Hydrogen_chloride(106) / '  0.08843713   0.378   3.  0.' /
DATA Hydrogen_chloride(107) / ' -2.433471     1.523   1.  2.' /
DATA Hydrogen_chloride(108) / ' -0.2636625    2.656   3.  2.' /
DATA Hydrogen_chloride(109) / '  0.6307008    1.338   2.  1.' /
DATA Hydrogen_chloride(110) / ' -0.6382638    2.828   2.  2.' /
DATA Hydrogen_chloride(111) / ' -0.006851438  0.75    7.  1.' /
DATA Hydrogen_chloride(112) / '  7.363661     0.644   1.  2. 2.   -1.141   -0.95    1.56    0.855    0. 0. 0.' /
DATA Hydrogen_chloride(113) / ' -1.262993     2.892   1.  2. 2.   -1.162   -0.92    1.14    0.91     0. 0. 0.' /
DATA Hydrogen_chloride(114) / ' -0.006539739  0.76    3.  2. 2.   -34.6    -1550.   1.06    0.942    0. 0. 0.' /
DATA Hydrogen_chloride(115) / ' -0.8752692    1.323   2.  2. 2.   -1.175   -1.2     0.94    0.702    0. 0. 0.' /
DATA Hydrogen_chloride(116) / ' -3.224835     0.693   2.  2. 2.   -0.99    -0.89    1.25    0.487    0. 0. 0.' /
DATA Hydrogen_chloride(117) / '                                      eta      beta    gamma   epsilon' /
DATA Hydrogen_chloride(118) / '                                   EXP[eta*(delta-epsilon)^2+beta*(tau-gamma)^2]' /
DATA Hydrogen_chloride(119) / '' /
DATA Hydrogen_chloride(120) / '' /
DATA Hydrogen_chloride(121) / '#AUX   !---Auxiliary function for Cp0' /
DATA Hydrogen_chloride(122) / 'CPP    !Ideal gas heat capacity function for hydrogen chloride of Thol et al. (2018).' /
DATA Hydrogen_chloride(123) / '?' /
DATA Hydrogen_chloride(124) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(125) / '?Thol, M., Dubberke, F.H., Baumh?gger, E., Span, R., and Vrabec, J., 2018.' /
DATA Hydrogen_chloride(126) / '?' /
DATA Hydrogen_chloride(127) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(128) / '0.                 !' /
DATA Hydrogen_chloride(129) / '10000.             !' /
DATA Hydrogen_chloride(130) / '0.                 !' /
DATA Hydrogen_chloride(131) / '0.                 !' /
DATA Hydrogen_chloride(132) / '1.0     8.3144598  !Reducing parameters for T, Cp0' /
DATA Hydrogen_chloride(133) / '1 3   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Hydrogen_chloride(134) / ' 3.5        0.0' /
DATA Hydrogen_chloride(135) / ' 0.0033327  300.0' /
DATA Hydrogen_chloride(136) / ' 0.935243   4000.0' /
DATA Hydrogen_chloride(137) / ' 0.209996   6300.0' /
DATA Hydrogen_chloride(138) / '' /
DATA Hydrogen_chloride(139) / '' /
DATA Hydrogen_chloride(140) / '#AUX   !---Auxiliary function for PX0' /
DATA Hydrogen_chloride(141) / 'PX0    !Helmholtz energy ideal-gas function for hydrogen chloride of Thol et al. (2018).' /
DATA Hydrogen_chloride(142) / '?' /
DATA Hydrogen_chloride(143) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(144) / '?Thol, M., Dubberke, F.H., Baumh?gger, E., Span, R., and Vrabec, J., 2018.' /
DATA Hydrogen_chloride(145) / '?' /
DATA Hydrogen_chloride(146) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(147) / '1 2  3  0 0  0 0 0               !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau))' /
DATA Hydrogen_chloride(148) / '  2.5                   1.0      !ai, ti for [ai*log(tau**ti)] terms' /
DATA Hydrogen_chloride(149) / ' -4.0690445266807433    0.0      !aj, ti for [ai*tau**ti] terms' /
DATA Hydrogen_chloride(150) / '  4.0257768311312594    1.0      !aj, ti for [ai*tau**ti] terms' /
DATA Hydrogen_chloride(151) / '  0.0033327  300.0               !aj, ti for [ai*log(1-exp(-ti/T)] terms' /
DATA Hydrogen_chloride(152) / '  0.935243   4000.0' /
DATA Hydrogen_chloride(153) / '  0.209996   6300.0' /
DATA Hydrogen_chloride(154) / '' /
DATA Hydrogen_chloride(155) / '' /
DATA Hydrogen_chloride(156) / '' /
DATA Hydrogen_chloride(157) / '' /
DATA Hydrogen_chloride(158) / '--------------------------------------------------------------------------------' /
DATA Hydrogen_chloride(159) / '' /
DATA Hydrogen_chloride(160) / '@EOS    !---Equation of state---' /
DATA Hydrogen_chloride(161) / 'FE1     !Helmholtz equation of state for hydrogen chloride of Thol et al. (2014).' /
DATA Hydrogen_chloride(162) / '          ?' /
DATA Hydrogen_chloride(163) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(164) / '          ?Thol, M., Piazza, L., and Span, R.' /
DATA Hydrogen_chloride(165) / '          ? "A New Functional Form for Equations of State for Some Polar and Weakly Associating Fluids,"' /
DATA Hydrogen_chloride(166) / '          ? Int. J. Thermophys., 35:783-811, 2014.' /
DATA Hydrogen_chloride(167) / '          ?' /
DATA Hydrogen_chloride(168) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(169) / '          159.01             !Lower temperature limit [K]' /
DATA Hydrogen_chloride(170) / '          330.0              !Upper temperature limit [K]' /
DATA Hydrogen_chloride(171) / '          20000.0            !Upper pressure limit [kPa]' /
DATA Hydrogen_chloride(172) / '          34.4               !Maximum density [mol/L]' /
DATA Hydrogen_chloride(173) / '          CP1                                    !Pointer to Cp0 model' /
DATA Hydrogen_chloride(174) / '          36.460939                              !Molar mass [g/mol]' /
DATA Hydrogen_chloride(175) / '          131.1                                  !Triple point temperature [K]' /
DATA Hydrogen_chloride(176) / '          14.033                                 !Pressure at triple point [kPa]' /
DATA Hydrogen_chloride(177) / '          34.3                                   !Density at triple point [mol/L]' /
DATA Hydrogen_chloride(178) / '          188.199                                !Normal boiling point temperature [K]' /
DATA Hydrogen_chloride(179) / '          0.128                                  !Acentric factor' /
DATA Hydrogen_chloride(180) / '          324.55        8274.9      11.271514    !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Hydrogen_chloride(181) / '          324.55                    11.271514    !Reducing parameters [K, mol/L]' /
DATA Hydrogen_chloride(182) / '          8.314472                               !Gas constant [J/mol-K]' /
DATA Hydrogen_chloride(183) / '            16  4    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Hydrogen_chloride(184) / '          -0.40937325         -0.75      1.  0.  !a(i),t(i),d(i),l(i)' /
DATA Hydrogen_chloride(185) / '           0.943994574        -0.25      1.  0.' /
DATA Hydrogen_chloride(186) / '          -1.78830477          1.25      1.  0.' /
DATA Hydrogen_chloride(187) / '           0.128619044         0.75      2.  0.' /
DATA Hydrogen_chloride(188) / '           0.00439018427      -1.0       3.  0.' /
DATA Hydrogen_chloride(189) / '           0.0130480908       -0.375     3.  0.' /
DATA Hydrogen_chloride(190) / '           0.00169387782       1.25      5.  0.' /
DATA Hydrogen_chloride(191) / '           0.75155906          2.375     1.  1.' /
DATA Hydrogen_chloride(192) / '          -0.800007427         3.0       1.  1.' /
DATA Hydrogen_chloride(193) / '           0.430935939         2.625     2.  1.' /
DATA Hydrogen_chloride(194) / '           0.00454319457       1.875     5.  1.' /
DATA Hydrogen_chloride(195) / '          -0.152172259         4.5       1.  2.' /
DATA Hydrogen_chloride(196) / '          -0.0436174059        5.75      3.  2.' /
DATA Hydrogen_chloride(197) / '          -0.00970625964       5.375     4.  2.' /
DATA Hydrogen_chloride(198) / '           0.0101144098        2.75      5.  2.' /
DATA Hydrogen_chloride(199) / '           0.00376991644      14.5       2.  3.' /
DATA Hydrogen_chloride(200) / '' /
DATA Hydrogen_chloride(201) / '' /
DATA Hydrogen_chloride(202) / '@AUX    !---Auxiliary function for Cp0' /
DATA Hydrogen_chloride(203) / 'CP1     !Ideal gas heat capacity function for hydrogen chloride of Thol et al. (2014).' /
DATA Hydrogen_chloride(204) / '          ?' /
DATA Hydrogen_chloride(205) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(206) / '          ?Thol, M., Piazza, L., and Span, R.' /
DATA Hydrogen_chloride(207) / '          ?' /
DATA Hydrogen_chloride(208) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(209) / '          0.                 !' /
DATA Hydrogen_chloride(210) / '          10000.             !' /
DATA Hydrogen_chloride(211) / '          0.                 !' /
DATA Hydrogen_chloride(212) / '          0.                 !' /
DATA Hydrogen_chloride(213) / '          1.0     8.314472   !Reducing parameters for T, Cp0' /
DATA Hydrogen_chloride(214) / '          3 1   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Hydrogen_chloride(215) / '           3.5               0.0' /
DATA Hydrogen_chloride(216) / '           0.00002557348     1.0' /
DATA Hydrogen_chloride(217) / '          -4.567927e-8       2.0' /
DATA Hydrogen_chloride(218) / '           1.054392          4028.112' /
DATA Hydrogen_chloride(219) / '' /
DATA Hydrogen_chloride(220) / '' /
DATA Hydrogen_chloride(221) / '' /
DATA Hydrogen_chloride(222) / '' /
DATA Hydrogen_chloride(223) / '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' /
DATA Hydrogen_chloride(224) / '' /
DATA Hydrogen_chloride(225) / '#TRN   !---ECS Transport---' /
DATA Hydrogen_chloride(226) / 'ECS    !Extended Corresponding States model (Propane reference) for hydrogen chloride.' /
DATA Hydrogen_chloride(227) / ':DOI: 10.6028/NIST.IR.8209' /
DATA Hydrogen_chloride(228) / '?' /
DATA Hydrogen_chloride(229) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(230) / '?*** ESTIMATION METHOD *** NOT STANDARD REFERENCE QUALITY ***' /
DATA Hydrogen_chloride(231) / '?Huber, M.L., "Models for the Viscosity, Thermal Conductivity, and Surface Tension' /
DATA Hydrogen_chloride(232) / '? of Selected Pure Fluids as Implemented in REFPROP v10.0," NISTIR 8209, 2018.' /
DATA Hydrogen_chloride(233) / '? doi: 10.6028/NIST.IR.8209' /
DATA Hydrogen_chloride(234) / '?' /
DATA Hydrogen_chloride(235) / '?VISCOSITY' /
DATA Hydrogen_chloride(236) / '? Comparisons with the data of Krynicki, K., Hennel, J.W., "Viscosity of Liquid Ammonia and Hydrogen Chloride," Acta Phys. Pol., 24(8):269, 1963,' /
DATA Hydrogen_chloride(237) / '? suggest an estimated uncertainty of 10% for the saturated liquid phase above 240 K.' /
DATA Hydrogen_chloride(238) / '?' /
DATA Hydrogen_chloride(239) / '?THERMAL CONDUCTIVITY' /
DATA Hydrogen_chloride(240) / '? Predictive model. Limited experimental data. Values based on method of extended corresponding states; estimated uncertainty approximately 10-20%.' /
DATA Hydrogen_chloride(241) / '?' /
DATA Hydrogen_chloride(242) / '?The Lennard-Jones parameters were estimated with the method of Chung.' /
DATA Hydrogen_chloride(243) / '?' /
DATA Hydrogen_chloride(244) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(245) / '159.07             !Lower temperature limit [K]' /
DATA Hydrogen_chloride(246) / '700.0              !Upper temperature limit [K]' /
DATA Hydrogen_chloride(247) / '50000.0            !Upper pressure limit [kPa]' /
DATA Hydrogen_chloride(248) / '34.39              !Maximum density [mol/L]' /
DATA Hydrogen_chloride(249) / 'FEQ PROPANE.FLD' /
DATA Hydrogen_chloride(250) / 'VS1                !Model for reference fluid viscosity' /
DATA Hydrogen_chloride(251) / 'TC1                !Model for reference fluid thermal conductivity' /
DATA Hydrogen_chloride(252) / 'NUL                !Large molecule identifier' /
DATA Hydrogen_chloride(253) / '1                  !Lennard-Jones flag (0 or 1) (0 => use estimates)' /
DATA Hydrogen_chloride(254) / '0.355              !Lennard-Jones coefficient sigma [nm]' /
DATA Hydrogen_chloride(255) / '257.8              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA Hydrogen_chloride(256) / '1  0  0                  !Number of terms in f_int term in Eucken correlation, spare1, spare2' /
DATA Hydrogen_chloride(257) / ' 0.0006        0. 0. 0.  !Coefficient, power of T, spare1, spare2 1.32' /
DATA Hydrogen_chloride(258) / '4  0  0                  !Number of terms in psi (visc shape factor): poly,spare1,spare2' /
DATA Hydrogen_chloride(259) / ' 0.615877      0. 0. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Hydrogen_chloride(260) / ' 0.55609       0. 1. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Hydrogen_chloride(261) / '-0.337867      0. 2. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Hydrogen_chloride(262) / ' 0.0681029     0. 3. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Hydrogen_chloride(263) / '2  0  0                  !Number of terms in chi (t.c. shape factor): poly,spare1,spare2' /
DATA Hydrogen_chloride(264) / ' 1.57373       0. 0. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Hydrogen_chloride(265) / '-0.17681       0. 1. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Hydrogen_chloride(266) / 'TK3                !Pointer to critical enhancement auxiliary function' /
DATA Hydrogen_chloride(267) / '' /
DATA Hydrogen_chloride(268) / '' /
DATA Hydrogen_chloride(269) / '#AUX   !---Auxiliary function for the thermal conductivity critical enhancement' /
DATA Hydrogen_chloride(270) / 'TK3    !Simplified thermal conductivity critical enhancement for hydrogen chloride of Perkins et al. (2013).' /
DATA Hydrogen_chloride(271) / '?' /
DATA Hydrogen_chloride(272) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(273) / '?Perkins, R.A., Sengers, J.V., Abdulagatov, I.M., and Huber, M.L.,' /
DATA Hydrogen_chloride(274) / '? "Simplified Model for the Critical Thermal-Conductivity Enhancement in Molecular Fluids,"' /
DATA Hydrogen_chloride(275) / '? Int. J. Thermophys., 34(2):191-212, 2013. doi: 10.1007/s10765-013-1409-z' /
DATA Hydrogen_chloride(276) / '?' /
DATA Hydrogen_chloride(277) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(278) / '0.                 !' /
DATA Hydrogen_chloride(279) / '10000.             !' /
DATA Hydrogen_chloride(280) / '0.                 !' /
DATA Hydrogen_chloride(281) / '0.                 !' /
DATA Hydrogen_chloride(282) / '9 0 0 0            !# terms:  CO2-terms, spare, spare, spare' /
DATA Hydrogen_chloride(283) / '1.0  1.0  1.0      !Reducing parameters for T, rho, tcx [mW/(m-K)]' /
DATA Hydrogen_chloride(284) / '0.63               !Nu (universal exponent)' /
DATA Hydrogen_chloride(285) / '1.239              !Gamma (universal exponent)' /
DATA Hydrogen_chloride(286) / '1.02               !R0 (universal amplitude)' /
DATA Hydrogen_chloride(287) / '0.063              !Z (universal exponent--not used for t.c., only viscosity)' /
DATA Hydrogen_chloride(288) / '1.0                !C (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA Hydrogen_chloride(289) / '0.154e-9           !Xi0 (amplitude) [m]' /
DATA Hydrogen_chloride(290) / '0.054              !Gam0 (amplitude) [-]' /
DATA Hydrogen_chloride(291) / '0.424e-9           !Qd_inverse (modified effective cutoff parameter) [m]; estimated-not fitted to data' /
DATA Hydrogen_chloride(292) / '487.0              !Tref (reference temperature)=1.5*Tc [K]' /
DATA Hydrogen_chloride(293) / '' /
DATA Hydrogen_chloride(294) / '' /
DATA Hydrogen_chloride(295) / '' /
DATA Hydrogen_chloride(296) / '' /
DATA Hydrogen_chloride(297) / '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' /
DATA Hydrogen_chloride(298) / '' /
DATA Hydrogen_chloride(299) / '#STN   !---Surface tension---' /
DATA Hydrogen_chloride(300) / 'ST1    !Surface tension model for hydrogen chloride of Mulero et al. (2014).' /
DATA Hydrogen_chloride(301) / ':DOI: 10.1063/1.4878755' /
DATA Hydrogen_chloride(302) / '?' /
DATA Hydrogen_chloride(303) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(304) / '?Mulero, A. and Cachadina, I.,' /
DATA Hydrogen_chloride(305) / '? "Recommended Correlations for the Surface Tension of Several Fluids' /
DATA Hydrogen_chloride(306) / '? Included in the REFPROP Program,"' /
DATA Hydrogen_chloride(307) / '? J. Phys. Chem. Ref. Data, 43, 023104, 2014.' /
DATA Hydrogen_chloride(308) / '? doi: 10.1063/1.4878755' /
DATA Hydrogen_chloride(309) / '?' /
DATA Hydrogen_chloride(310) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(311) / '0.                 !' /
DATA Hydrogen_chloride(312) / '10000.             !' /
DATA Hydrogen_chloride(313) / '0.                 !' /
DATA Hydrogen_chloride(314) / '0.                 !' /
DATA Hydrogen_chloride(315) / '1                  !Number of terms in surface tension model' /
DATA Hydrogen_chloride(316) / '324.55             !Critical temperature used in fit (dummy)' /
DATA Hydrogen_chloride(317) / '0.05994   1.0953   !Sigma0 and n' /
DATA Hydrogen_chloride(318) / '' /
DATA Hydrogen_chloride(319) / '' /
DATA Hydrogen_chloride(320) / '#PS    !---Vapor pressure---' /
DATA Hydrogen_chloride(321) / 'PS5    !Vapor pressure equation for hydrogen chloride of Thol et al. (2018)' /
DATA Hydrogen_chloride(322) / '?' /
DATA Hydrogen_chloride(323) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(324) / '?Thol, M., Dubberke, F.H., Baumh?gger, E., Span, R., and Vrabec, J., 2018.' /
DATA Hydrogen_chloride(325) / '?' /
DATA Hydrogen_chloride(326) / '?Functional Form:  P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T] where Theta=1-T/Tc, Tc and Pc' /
DATA Hydrogen_chloride(327) / '? are the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Hydrogen_chloride(328) / '?' /
DATA Hydrogen_chloride(329) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(330) / '0.                 !' /
DATA Hydrogen_chloride(331) / '10000.             !' /
DATA Hydrogen_chloride(332) / '0.                 !' /
DATA Hydrogen_chloride(333) / '0.                 !' /
DATA Hydrogen_chloride(334) / '324.68  8313.5     !Reducing parameters' /
DATA Hydrogen_chloride(335) / '5 0 0 0 0 0        !Number of terms in equation' /
DATA Hydrogen_chloride(336) / '-6.730   1.0       !Coefficients and exponents' /
DATA Hydrogen_chloride(337) / ' 1.464   1.5' /
DATA Hydrogen_chloride(338) / '-1.994   3.12' /
DATA Hydrogen_chloride(339) / ' 1.283   3.95' /
DATA Hydrogen_chloride(340) / '-2.062   4.8' /
DATA Hydrogen_chloride(341) / '' /
DATA Hydrogen_chloride(342) / '' /
DATA Hydrogen_chloride(343) / '#DL    !---Saturated liquid density---' /
DATA Hydrogen_chloride(344) / 'DL1    !Saturated liquid density equation for hydrogen chloride of Thol et al. (2018)' /
DATA Hydrogen_chloride(345) / '?' /
DATA Hydrogen_chloride(346) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(347) / '?Thol, M., Dubberke, F.H., Baumh?gger, E., Span, R., and Vrabec, J., 2018.' /
DATA Hydrogen_chloride(348) / '?' /
DATA Hydrogen_chloride(349) / '?Functional Form:  D=Dc*[1+SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Hydrogen_chloride(350) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Hydrogen_chloride(351) / '?' /
DATA Hydrogen_chloride(352) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(353) / '0.                 !' /
DATA Hydrogen_chloride(354) / '10000.             !' /
DATA Hydrogen_chloride(355) / '0.                 !' /
DATA Hydrogen_chloride(356) / '0.                 !' /
DATA Hydrogen_chloride(357) / '324.68  11.87      !Reducing parameters' /
DATA Hydrogen_chloride(358) / '5 0 0 0 0 0        !Number of terms in equation' /
DATA Hydrogen_chloride(359) / ' 2.547   0.418     !Coefficients and exponents' /
DATA Hydrogen_chloride(360) / '-0.631   1.12' /
DATA Hydrogen_chloride(361) / ' 1.750   1.86' /
DATA Hydrogen_chloride(362) / '-1.922   2.66' /
DATA Hydrogen_chloride(363) / ' 1.030   3.57' /
DATA Hydrogen_chloride(364) / '' /
DATA Hydrogen_chloride(365) / '' /
DATA Hydrogen_chloride(366) / '#DV    !---Saturated vapor density---' /
DATA Hydrogen_chloride(367) / 'DV3    !Saturated vapor density equation for hydrogen chloride of Thol et al. (2018)' /
DATA Hydrogen_chloride(368) / '?' /
DATA Hydrogen_chloride(369) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(370) / '?Thol, M., Dubberke, F.H., Baumh?gger, E., Span, R., and Vrabec, J., 2018.' /
DATA Hydrogen_chloride(371) / '?' /
DATA Hydrogen_chloride(372) / '?Functional Form:  D=Dc*EXP[SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Hydrogen_chloride(373) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Hydrogen_chloride(374) / '?' /
DATA Hydrogen_chloride(375) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Hydrogen_chloride(376) / '0.                 !' /
DATA Hydrogen_chloride(377) / '10000.             !' /
DATA Hydrogen_chloride(378) / '0.                 !' /
DATA Hydrogen_chloride(379) / '0.                 !' /
DATA Hydrogen_chloride(380) / '324.68  11.87      !Reducing parameters' /
DATA Hydrogen_chloride(381) / '6 0 0 0 0 0        !Number of terms in equation' /
DATA Hydrogen_chloride(382) / '-2.5676   0.417    !Coefficients and exponents' /
DATA Hydrogen_chloride(383) / '-4.1055   0.923' /
DATA Hydrogen_chloride(384) / '-12.068   2.57' /
DATA Hydrogen_chloride(385) / '-29.03    5.54' /
DATA Hydrogen_chloride(386) / '-54.93   10.5' /
DATA Hydrogen_chloride(387) / '-222.7   23.3' /
DATA Hydrogen_chloride(388) / '' /
DATA Hydrogen_chloride(389) / '' /
DATA Hydrogen_chloride(390) / '@END' /
DATA Hydrogen_chloride(391) / 'c        1         2         3         4         5         6         7         8' /
DATA Hydrogen_chloride(392) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
! #######################################################
character(256), TARGET :: hydrogen_normal(582)
DATA hydrogen_normal(1) / 'hydrogen (normal)  !short name' /
DATA hydrogen_normal(2) / '1333-74-0          !CAS number' /
DATA hydrogen_normal(3) / 'hydrogen (normal)  !full name' /
DATA hydrogen_normal(4) / 'H2                 !chemical formula' /
DATA hydrogen_normal(5) / 'R-702              !synonym' /
DATA hydrogen_normal(6) / '2.01588            !molecular weight [g/mol]' /
DATA hydrogen_normal(7) / '13.957             !triple point temperature [K]' /
DATA hydrogen_normal(8) / '20.369             !normal boiling point [K]' /
DATA hydrogen_normal(9) / '33.145             !critical temperature [K]' /
DATA hydrogen_normal(10) / '1296.4             !critical pressure [kPa]' /
DATA hydrogen_normal(11) / '15.508             !critical density [mol/L]' /
DATA hydrogen_normal(12) / '-0.219             !acentric factor' /
DATA hydrogen_normal(13) / '0.0                !dipole moment [Debye]' /
DATA hydrogen_normal(14) / 'NBP                !default reference state' /
DATA hydrogen_normal(15) / '9.1                !version number' /
DATA hydrogen_normal(16) / '1049               !UN Number' /
DATA hydrogen_normal(17) / 'other              !family' /
DATA hydrogen_normal(18) / '285.83             !heating value (gross or superior) [kJ/mol]' /
DATA hydrogen_normal(19) / 'A3                 !Safety Group (ASHRAE Standard 34, 2010)' /
DATA hydrogen_normal(20) / '' /
DATA hydrogen_normal(21) / '' /
DATA hydrogen_normal(22) / '#EOS               !equation of state specification' /
DATA hydrogen_normal(23) / 'FEQ  Helmholtz equation of state for normal hydrogen of Leachman et al. (2009).' /
DATA hydrogen_normal(24) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(25) / '?Leachman, J.W., Jacobsen, R.T, Penoncello, S.G., Lemmon, E.W.' /
DATA hydrogen_normal(26) / '?"Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and' /
DATA hydrogen_normal(27) / '?Orthohydrogen,"' /
DATA hydrogen_normal(28) / '?J. Phys. Chem. Ref. Data, 38(3):721-748, 2009.' /
DATA hydrogen_normal(29) / '?\' /
DATA hydrogen_normal(30) / '?The uncertainty in density is 0.1% at temperatures from the triple point' /
DATA hydrogen_normal(31) / '?to 250 K and at pressures up to 40 MPa, except in the critical region,' /
DATA hydrogen_normal(32) / '?where an uncertainty of 0.2% in pressure is generally attained.  In the' /
DATA hydrogen_normal(33) / '?region between 250 and 450 K and at pressures to 300 MPa, the' /
DATA hydrogen_normal(34) / '?uncertainty in density is 0.04%.  At temperatures between 450 and 1000' /
DATA hydrogen_normal(35) / '?K, the uncertainty in density increases to 1%.  At pressures between 300' /
DATA hydrogen_normal(36) / '?and 2000 MPa, the uncertainty in density is 8%.  Speed of sound data are' /
DATA hydrogen_normal(37) / '?represented within 0.5% below 100 MPa. The estimated uncertainty for' /
DATA hydrogen_normal(38) / '?heat capacities is 1.0%.  The estimated uncertainties of vapor pressures' /
DATA hydrogen_normal(39) / '?and saturated liquid densities calculated using the Maxwell criterion' /
DATA hydrogen_normal(40) / '?are 0.2% for each property.' /
DATA hydrogen_normal(41) / '?\' /
DATA hydrogen_normal(42) / '!end of info section' /
DATA hydrogen_normal(43) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(44) / '1000.0             !upper temperature limit [K]' /
DATA hydrogen_normal(45) / '2000000.0          !upper pressure limit [kPa]' /
DATA hydrogen_normal(46) / '102.0              !maximum density [mol/L]' /
DATA hydrogen_normal(47) / 'CPP                                    !pointer to Cp0 model' /
DATA hydrogen_normal(48) / '2.01588                                !molecular weight [g/mol]' /
DATA hydrogen_normal(49) / '13.957                                 !triple point temperature [K]' /
DATA hydrogen_normal(50) / '7.357                                  !pressure at triple point [kPa]' /
DATA hydrogen_normal(51) / '38.2                                   !density at triple point [mol/L]' /
DATA hydrogen_normal(52) / '20.369                                 !normal boiling point temperature [K]' /
DATA hydrogen_normal(53) / '-0.219                                 !acentric factor' /
DATA hydrogen_normal(54) / '33.145        1296.4      15.508       !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA hydrogen_normal(55) / '33.145                    15.508       !reducing parameters [K, mol/L]' /
DATA hydrogen_normal(56) / '8.314472                               !gas constant [J/mol-K]' /
DATA hydrogen_normal(57) / '   9  4      5 12      0  0            !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA hydrogen_normal(58) / ' -0.693643D+01  0.6844   1.0  0.0      !a(i),t(i),d(i),l(i)' /
DATA hydrogen_normal(59) / '  0.100000D-01  1.000    4.0  0.0' /
DATA hydrogen_normal(60) / '  0.211010D+01  0.989    1.0  0.0' /
DATA hydrogen_normal(61) / '  0.452059D+01  0.489    1.0  0.0' /
DATA hydrogen_normal(62) / '  0.732564D+00  0.803    2.0  0.0' /
DATA hydrogen_normal(63) / ' -0.134086D+01  1.1444   2.0  0.0' /
DATA hydrogen_normal(64) / '  0.130985D+00  1.409    3.0  0.0' /
DATA hydrogen_normal(65) / ' -0.777414D+00  1.754    1.0  1.0' /
DATA hydrogen_normal(66) / '  0.351944D+00  1.311    3.0  1.0' /
DATA hydrogen_normal(67) / ' -0.211716D-01  4.187    2.0  2.0  2.0  -1.685   -0.1710  0.7164  1.506   0. 0. 0.' /
DATA hydrogen_normal(68) / '  0.226312D-01  5.646    1.0  2.0  2.0  -0.489   -0.2245  1.3444  0.156   0. 0. 0.' /
DATA hydrogen_normal(69) / '  0.321870D-01  0.791    3.0  2.0  2.0  -0.103   -0.1304  1.4517  1.736   0. 0. 0.' /
DATA hydrogen_normal(70) / ' -0.231752D-01  7.249    1.0  2.0  2.0  -2.506   -0.2785  0.7204  0.670   0. 0. 0.' /
DATA hydrogen_normal(71) / '  0.557346D-01  2.986    1.0  2.0  2.0  -1.607   -0.3967  1.5445  1.662   0. 0. 0.' /
DATA hydrogen_normal(72) / '' /
DATA hydrogen_normal(73) / '' /
DATA hydrogen_normal(74) / '#AUX               !auxiliary model specification' /
DATA hydrogen_normal(75) / 'CPP  ideal gas heat capacity function' /
DATA hydrogen_normal(76) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(77) / '?Leachman, J.W., Jacobsen, R.T, Penoncello, S.G., Lemmon, E.W.' /
DATA hydrogen_normal(78) / '?"Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and' /
DATA hydrogen_normal(79) / '?Orthohydrogen,"' /
DATA hydrogen_normal(80) / '?J. Phys. Chem. Ref. Data, 38(3):721-748, 2009.' /
DATA hydrogen_normal(81) / '?\' /
DATA hydrogen_normal(82) / '!end of info section' /
DATA hydrogen_normal(83) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(84) / '1000.0             !upper temperature limit [K]' /
DATA hydrogen_normal(85) / '0.0                !upper pressure limit [kPa]' /
DATA hydrogen_normal(86) / '0.0                !maximum density [mol/L]' /
DATA hydrogen_normal(87) / '1.0          8.314472                  !reducing parameters for T, Cp0' /
DATA hydrogen_normal(88) / '  1  5    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA hydrogen_normal(89) / '  2.500          0.0' /
DATA hydrogen_normal(90) / '  1.616        531.0' /
DATA hydrogen_normal(91) / '-0.4117        751.0' /
DATA hydrogen_normal(92) / '-0.7920       1989.0' /
DATA hydrogen_normal(93) / ' 0.7580       2484.0' /
DATA hydrogen_normal(94) / '  1.217       6859.0' /
DATA hydrogen_normal(95) / '' /
DATA hydrogen_normal(96) / '' /
DATA hydrogen_normal(97) / '#AUX               !auxiliary model specification' /
DATA hydrogen_normal(98) / 'PH0  Helmholtz form for the ideal-gas state' /
DATA hydrogen_normal(99) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(100) / '?Leachman, J.W., Jacobsen, R.T, Penoncello, S.G., Lemmon, E.W.' /
DATA hydrogen_normal(101) / '?"Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and' /
DATA hydrogen_normal(102) / '?Orthohydrogen,"' /
DATA hydrogen_normal(103) / '?J. Phys. Chem. Ref. Data, 38(3):721-748, 2009.' /
DATA hydrogen_normal(104) / '?\' /
DATA hydrogen_normal(105) / '!end of info section' /
DATA hydrogen_normal(106) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(107) / '1000.0             !upper temperature limit [K]' /
DATA hydrogen_normal(108) / '0.0                !upper pressure limit [kPa]' /
DATA hydrogen_normal(109) / '0.0                !maximum density [mol/L]' /
DATA hydrogen_normal(110) / '1 2  5  0 0  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau))' /
DATA hydrogen_normal(111) / '    1.5000000000    1.0000000000   !ai, ti for [ai*log(tau**ti)] terms' /
DATA hydrogen_normal(112) / '   -1.4579856475    0.0000000000   !aj, ti for [ai*tau**ti] terms' /
DATA hydrogen_normal(113) / '    1.8880767820    1.0000000000' /
DATA hydrogen_normal(114) / '    1.6160000000  -16.0205159149   !aj, ti for [ai*log(1-exp(ti*tau)] terms' /
DATA hydrogen_normal(115) / '   -0.4117000000  -22.6580178006' /
DATA hydrogen_normal(116) / '   -0.7920000000  -60.0090511389' /
DATA hydrogen_normal(117) / '    0.7580000000  -74.9434303817' /
DATA hydrogen_normal(118) / '    1.2170000000 -206.9392065168' /
DATA hydrogen_normal(119) / '' /
DATA hydrogen_normal(120) / '' /
DATA hydrogen_normal(121) / '@EOS               !equation of state specification' /
DATA hydrogen_normal(122) / 'BWR  MBWR equation of state for hydrogen of Younglove (1982).' /
DATA hydrogen_normal(123) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(124) / '?Younglove, B.A.,' /
DATA hydrogen_normal(125) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA hydrogen_normal(126) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA hydrogen_normal(127) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA hydrogen_normal(128) / '?\' /
DATA hydrogen_normal(129) / '?The uncertainties in density are 0.1% in the liquid phase, 0.25% in the' /
DATA hydrogen_normal(130) / '?vapor phase, and 0.2% in the supercritical region.  The uncertainty in' /
DATA hydrogen_normal(131) / '?heat capacity is 3% and the uncertainty in speed of sound is 2% in the' /
DATA hydrogen_normal(132) / '?liquid phase and 1% elsewhere.' /
DATA hydrogen_normal(133) / '?\' /
DATA hydrogen_normal(134) / '!end of info section' /
DATA hydrogen_normal(135) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(136) / '400.0              !upper temperature limit [K]' /
DATA hydrogen_normal(137) / '121000.0           !upper pressure limit [kPa]' /
DATA hydrogen_normal(138) / '38.148             !maximum density [mol/L]' /
DATA hydrogen_normal(139) / 'CP1                                    !pointer to Cp0 model' /
DATA hydrogen_normal(140) / '2.01594                                !molecular weight [g/mol]' /
DATA hydrogen_normal(141) / '13.957                                 !triple point temperature [K]' /
DATA hydrogen_normal(142) / '7.70                                   !pressure at triple point [kPa]' /
DATA hydrogen_normal(143) / '38.3                                   !density at triple point [mol/L]' /
DATA hydrogen_normal(144) / '20.39                                  !normal boiling point temperature [K]' /
DATA hydrogen_normal(145) / '-0.214                                 !acentric factor' /
DATA hydrogen_normal(146) / '33.19        1315.0       14.94        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA hydrogen_normal(147) / '33.19                     14.94        !reducing parameters [K, mol/L]' /
DATA hydrogen_normal(148) / '15.6173762                             !gamma' /
DATA hydrogen_normal(149) / '0.0831434                              !gas constant [L-bar/mol-K]' /
DATA hydrogen_normal(150) / '      32       1                       !Nterm, Ncoeff per term' /
DATA hydrogen_normal(151) / '  0.4675528393416d-03     0.4289274251454d-01     -0.5164085596504d-00' /
DATA hydrogen_normal(152) / '  0.2961790279801d+01    -0.3027194968412d+02      0.1908100320379d-04' /
DATA hydrogen_normal(153) / ' -0.1339776859288d-02     0.3056473115421d-00      0.5161197159532d+02' /
DATA hydrogen_normal(154) / '  0.1999981550224d-06     0.2896367059356d-03     -0.2257803939041d-01' /
DATA hydrogen_normal(155) / ' -0.2287392761826d-05     0.2446261478645d-04     -0.1718181601119d-02' /
DATA hydrogen_normal(156) / ' -0.5465142603459d-06     0.4051941401315d-08      0.1157595123961d-05' /
DATA hydrogen_normal(157) / ' -0.1269162728389d-07    -0.4983023605519d+02     -0.1606676092098d+03' /
DATA hydrogen_normal(158) / ' -0.1926799185310d-00     0.9319894638928d+01     -0.3222596554434d-03' /
DATA hydrogen_normal(159) / '  0.1206839307669d-02    -0.3841588197470d-06     -0.4036157453608d-04' /
DATA hydrogen_normal(160) / ' -0.1250868123513d-09     0.1976107321888d-08     -0.2411883474011d-12' /
DATA hydrogen_normal(161) / ' -0.4127551498251d-12     0.8917972883610d-11' /
DATA hydrogen_normal(162) / '' /
DATA hydrogen_normal(163) / '' /
DATA hydrogen_normal(164) / '#AUX               !auxiliary model specification' /
DATA hydrogen_normal(165) / 'CP1  ideal gas heat capacity function' /
DATA hydrogen_normal(166) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(167) / '?McCarty, R.D., Hord, J., and Roder, H.M.,' /
DATA hydrogen_normal(168) / '? "Selected Properties of Hydrogen (Engineering Design Data),"' /
DATA hydrogen_normal(169) / '? NBS Monograph 168, National Bureau of Standards, Boulder, 1981.' /
DATA hydrogen_normal(170) / '?\' /
DATA hydrogen_normal(171) / '?\' /
DATA hydrogen_normal(172) / '!end of info section' /
DATA hydrogen_normal(173) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(174) / '500.0              !upper temperature limit [K]' /
DATA hydrogen_normal(175) / '0.0                !upper pressure limit [kPa]' /
DATA hydrogen_normal(176) / '0.0                !maximum density [mol/L]' /
DATA hydrogen_normal(177) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA hydrogen_normal(178) / '  17 0    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA hydrogen_normal(179) / ' 0.12155215d+11   -7.00' /
DATA hydrogen_normal(180) / '-0.36396763d+10   -6.00' /
DATA hydrogen_normal(181) / ' 0.43375265d+09   -5.00' /
DATA hydrogen_normal(182) / '-0.23085817d+08   -4.00' /
DATA hydrogen_normal(183) / '-0.38680927d+04   -3.00' /
DATA hydrogen_normal(184) / ' 0.88240136d+05   -2.00' /
DATA hydrogen_normal(185) / '-0.78587085d+04   -1.00' /
DATA hydrogen_normal(186) / ' 0.72480209d+03    0.00' /
DATA hydrogen_normal(187) / '-0.18426806d+03    0.50' /
DATA hydrogen_normal(188) / ' 0.21801550d+02    1.00' /
DATA hydrogen_normal(189) / '-0.13051820d+01    1.50' /
DATA hydrogen_normal(190) / ' 0.21003175d-01    2.00' /
DATA hydrogen_normal(191) / ' 0.23911604d-02    2.50' /
DATA hydrogen_normal(192) / '-0.18240547d-03    3.00' /
DATA hydrogen_normal(193) / ' 0.56149561d-05    3.50' /
DATA hydrogen_normal(194) / '-0.73803310d-07    4.00' /
DATA hydrogen_normal(195) / ' 0.66357755d-11    5.00' /
DATA hydrogen_normal(196) / '' /
DATA hydrogen_normal(197) / '' /
DATA hydrogen_normal(198) / '@EOS               !equation of state specification' /
DATA hydrogen_normal(199) / 'FEK  Helmholtz equation of state for hydrogen of Kunz and Wagner (2004).' /
DATA hydrogen_normal(200) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(201) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA hydrogen_normal(202) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA hydrogen_normal(203) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA hydrogen_normal(204) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA hydrogen_normal(205) / '?\' /
DATA hydrogen_normal(206) / '!end of info section' /
DATA hydrogen_normal(207) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(208) / '400.0              !upper temperature limit [K]' /
DATA hydrogen_normal(209) / '121000.0           !upper pressure limit [kPa]' /
DATA hydrogen_normal(210) / '38.148             !maximum density [mol/L]' /
DATA hydrogen_normal(211) / 'PHK                                    !pointer to Cp0 model' /
DATA hydrogen_normal(212) / '2.01588                                !molecular weight [g/mol]' /
DATA hydrogen_normal(213) / '13.957                                 !triple point temperature [K]' /
DATA hydrogen_normal(214) / '6.669                                  !pressure at triple point [kPa]' /
DATA hydrogen_normal(215) / '38.33                                  !density at triple point [mol/L]' /
DATA hydrogen_normal(216) / ' 20.38                                 !normal boiling point temperature [K]' /
DATA hydrogen_normal(217) / '-0.2187                                !acentric factor' /
DATA hydrogen_normal(218) / '33.19        1315.0      14.94         !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA hydrogen_normal(219) / '33.19                    14.94         !reducing parameters [K, mol/L]' /
DATA hydrogen_normal(220) / '8.314472                               !gas constant [J/mol-K]' /
DATA hydrogen_normal(221) / '  14  4      0  0      0  0            !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA hydrogen_normal(222) / ' 0.53579928451252d1     0.50   1.  0' /
DATA hydrogen_normal(223) / '-0.62050252530595d1     0.625  1.  0' /
DATA hydrogen_normal(224) / ' 0.13830241327086       0.375  2.  0' /
DATA hydrogen_normal(225) / '-0.71397954896129d-1    0.625  2.  0' /
DATA hydrogen_normal(226) / ' 0.15474053959733d-1    1.125  4.  0' /
DATA hydrogen_normal(227) / '-0.14976806405771       2.625  1.  1' /
DATA hydrogen_normal(228) / '-0.26368723988451d-1    0.0    5.  1' /
DATA hydrogen_normal(229) / ' 0.56681303156066d-1    0.25   5.  1' /
DATA hydrogen_normal(230) / '-0.60063958030436d-1    1.375  5.  1' /
DATA hydrogen_normal(231) / '-0.45043942027132       4.0    1.  2' /
DATA hydrogen_normal(232) / ' 0.42478840244500       4.25   1.  2' /
DATA hydrogen_normal(233) / '-0.21997640827139d-1    5.0    2.  3' /
DATA hydrogen_normal(234) / '-0.10499521374530d-1    8.0    5.  3' /
DATA hydrogen_normal(235) / '-0.28955902866816d-2    8.0    1.  5' /
DATA hydrogen_normal(236) / '' /
DATA hydrogen_normal(237) / '' /
DATA hydrogen_normal(238) / '#AUX               !auxiliary model specification' /
DATA hydrogen_normal(239) / 'PHK  Helmholtz form for the ideal-gas state for hydrogen of Kunz and Wagner (2004).' /
DATA hydrogen_normal(240) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(241) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA hydrogen_normal(242) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA hydrogen_normal(243) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA hydrogen_normal(244) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA hydrogen_normal(245) / '?\' /
DATA hydrogen_normal(246) / '!end of info section' /
DATA hydrogen_normal(247) / '0.                 !lower temperature limit [K]' /
DATA hydrogen_normal(248) / '1000.0             !upper temperature limit [K]' /
DATA hydrogen_normal(249) / '0.0                !upper pressure limit [kPa]' /
DATA hydrogen_normal(250) / '0.0                !maximum density [mol/L]' /
DATA hydrogen_normal(251) / '1 2  0  2 2  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA hydrogen_normal(252) / '    1.47906      1.             !ai, ti for [ai*log(tau**ti)] terms' /
DATA hydrogen_normal(253) / '   13.796443393  0.             !aj, ti for [ai*tau**ti] terms' /
DATA hydrogen_normal(254) / ' -175.864487294  1.' /
DATA hydrogen_normal(255) / '   -0.45444      9.84763483     !aj, ti for cosh and sinh terms' /
DATA hydrogen_normal(256) / '    1.3756      50.367279301' /
DATA hydrogen_normal(257) / '    0.95806      6.891654113' /
DATA hydrogen_normal(258) / '    1.56039     49.76529075' /
DATA hydrogen_normal(259) / '' /
DATA hydrogen_normal(260) / '' /
DATA hydrogen_normal(261) / '@EOS               !equation of state specification' /
DATA hydrogen_normal(262) / 'FE1  Helmholtz equation of state for hydrogen of Bender (1982).' /
DATA hydrogen_normal(263) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(264) / '?Bender, E.' /
DATA hydrogen_normal(265) / '?"Equation of state of normal hydrogen in the range 18 to 700 K and 1 to 500 bar,"' /
DATA hydrogen_normal(266) / '?VDI Forschungsheft N 609, pp. 15-20, 1982.' /
DATA hydrogen_normal(267) / '?\' /
DATA hydrogen_normal(268) / '?Cp0 equation was taken from McCarty et al. (1981) since Bender equation is' /
DATA hydrogen_normal(269) / '?split in two pieces from 10 to 250 K and from 250 to 600 K.' /
DATA hydrogen_normal(270) / '?\' /
DATA hydrogen_normal(271) / '!end of info section' /
DATA hydrogen_normal(272) / '18.0               !lower temperature limit [K]' /
DATA hydrogen_normal(273) / '700.0              !upper temperature limit [K]' /
DATA hydrogen_normal(274) / '50000.0            !upper pressure limit [kPa]' /
DATA hydrogen_normal(275) / '38.74              !maximum density [mol/L]  (change to 41 when melting line is added)' /
DATA hydrogen_normal(276) / 'CP1                                    !pointer to Cp0 model' /
DATA hydrogen_normal(277) / '2.01594                                !molecular weight [g/mol]' /
DATA hydrogen_normal(278) / '13.957                                 !triple point temperature [K]' /
DATA hydrogen_normal(279) / '8.736                                  !pressure at triple point [kPa]' /
DATA hydrogen_normal(280) / '38.7                                   !density at triple point [mol/L]' /
DATA hydrogen_normal(281) / '20.39                                  !normal boiling point temperature [K]' /
DATA hydrogen_normal(282) / '-0.218                                 !acentric factor' /
DATA hydrogen_normal(283) / '33.24           1303.0    15.37744     !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA hydrogen_normal(284) / '33.24                     15.37744     !reducing parameters [K, mol/L]' /
DATA hydrogen_normal(285) / '8.3143                                 !gas constant [J/mol-K]' /
DATA hydrogen_normal(286) / '  22  5      0  0      0  0            !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA hydrogen_normal(287) / '  0.133442326203D+01   3.    0.    0.    0.0' /
DATA hydrogen_normal(288) / ' -0.104116843433D+01   4.    0.    0.    0.0' /
DATA hydrogen_normal(289) / '  0.227202245707D+00   5.    0.    0.    0.0' /
DATA hydrogen_normal(290) / '  0.300374270906D+00   0.    1.    0.    0.0' /
DATA hydrogen_normal(291) / ' -0.463984214813D+00   1.    1.    0.    0.0' /
DATA hydrogen_normal(292) / ' -0.178010492282D+01   2.    1.    0.    0.0' /
DATA hydrogen_normal(293) / '  0.100460103605D+01   3.    1.    0.    0.0' /
DATA hydrogen_normal(294) / ' -0.187200622541D+00   4.    1.    0.    0.0' /
DATA hydrogen_normal(295) / '  0.980276957749D-02   0.    2.    0.    0.0' /
DATA hydrogen_normal(296) / '  0.543224866339D-01   1.    2.    0.    0.0' /
DATA hydrogen_normal(297) / ' -0.263496312610D-01   2.    2.    0.    0.0' /
DATA hydrogen_normal(298) / '  0.315432315759D-01   0.    3.    0.    0.0' /
DATA hydrogen_normal(299) / ' -0.525788294155D-01   1.    3.    0.    0.0' /
DATA hydrogen_normal(300) / ' -0.685380627808D-02   0.    4.    0.    0.0' /
DATA hydrogen_normal(301) / '  0.344540276656D-01   1.    4.    0.    0.0' /
DATA hydrogen_normal(302) / ' -0.555747275982D-03   1.    5.    0.    0.0' /
DATA hydrogen_normal(303) / ' -0.133442326203D+01   3.    0.    2.    0.711139834571' /
DATA hydrogen_normal(304) / '  0.104116843433D+01   4.    0.    2.    0.711139834571' /
DATA hydrogen_normal(305) / ' -0.227202245707D+00   5.    0.    2.    0.711139834571' /
DATA hydrogen_normal(306) / ' -0.378598758038D+00   3.    2.    2.    0.711139834571' /
DATA hydrogen_normal(307) / '  0.249888797892D+00   4.    2.    2.    0.711139834571' /
DATA hydrogen_normal(308) / ' -0.498847982876D-01   5.    2.    2.    0.711139834571' /
DATA hydrogen_normal(309) / '' /
DATA hydrogen_normal(310) / '' /
DATA hydrogen_normal(311) / '' /
DATA hydrogen_normal(312) / '#ETA               !viscosity model specification' /
DATA hydrogen_normal(313) / 'VS0  pure fluid viscosity model from symbolic regression (Muzny, Huber, Kazakov) (2013).' /
DATA hydrogen_normal(314) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(315) / '?Muzny, C.D., Huber, M.L., and Kazakov, A.F.,' /
DATA hydrogen_normal(316) / '? "Correlation for the Viscosity of normal hydrogen obtained from symbolic regression"' /
DATA hydrogen_normal(317) / '? submitted to J. Chem. Eng. Data, 2013' /
DATA hydrogen_normal(318) / '?\' /
DATA hydrogen_normal(319) / '? The estimated uncertainty is 4 % for the saturated liquid from the triple point to 31 K, with larger deviations' /
DATA hydrogen_normal(320) / '? as the critical region is approached. The estimated uncertainty is 4 % for the supercritical fluid phase at pressures to 200 MPa.' /
DATA hydrogen_normal(321) / '? For the limited range of 200 K to 400 K at pressures up to 0.1 MPa, the uncertainty is 0.1 %.' /
DATA hydrogen_normal(322) / '?\' /
DATA hydrogen_normal(323) / '!end of info section' /
DATA hydrogen_normal(324) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(325) / '2000.0             !upper temperature limit [K]' /
DATA hydrogen_normal(326) / '2000000.0          !upper pressure limit [kPa]' /
DATA hydrogen_normal(327) / '102.0              !maximum density [mol/L]' /
DATA hydrogen_normal(328) / 'H2A                !pointer to hardcoded model' /
DATA hydrogen_normal(329) / '0 0 0 0 0 0 0 0    !number of terms for various pieces' /
DATA hydrogen_normal(330) / '1.0      1.0      1.0               !reducing parameters for T, rho, eta' /
DATA hydrogen_normal(331) / 'NUL                !pointer to critical enhancement auxiliary function' /
DATA hydrogen_normal(332) / '' /
DATA hydrogen_normal(333) / '' /
DATA hydrogen_normal(334) / '#TCX               !thermal conductivity model specification' /
DATA hydrogen_normal(335) / 'TC1  pure fluid thermal conductivity model of Assael et al. (2011).' /
DATA hydrogen_normal(336) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(337) / '? Assael, M.J., Assael. J.-A.M., Huber, M.L., Perkins, R.A. and Takata, Y.' /
DATA hydrogen_normal(338) / '? "Correlation of the Thermal Conductivity of Normal and Parahydrogen' /
DATA hydrogen_normal(339) / '? from the Triple Point to 1000 K and up to 100 MPa",' /
DATA hydrogen_normal(340) / '? J. Phys. Chem. Ref. Data, Vol.40, No. 3(2011) pp.1-13.' /
DATA hydrogen_normal(341) / '?' /
DATA hydrogen_normal(342) / '? The estimated uncertainty is less than 4% from 100 K to 1000 K at pressures to 100 MPa.' /
DATA hydrogen_normal(343) / '? For temperatures from the triple point to 100 K, at pressures' /
DATA hydrogen_normal(344) / '? to 12 MPa, we estimate the uncertainty to be 7%, except near the critical point.' /
DATA hydrogen_normal(345) / '? The model behaves in a physically reasonable manner for extrapolations to pressures above' /
DATA hydrogen_normal(346) / '? 12 MPa at temperatures below 100 K, but will be subject to larger uncertainties.' /
DATA hydrogen_normal(347) / '?\' /
DATA hydrogen_normal(348) / '!end of info section' /
DATA hydrogen_normal(349) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(350) / '1000.0             !upper temperature limit [K]' /
DATA hydrogen_normal(351) / '2000000.0          !upper pressure limit [kPa]' /
DATA hydrogen_normal(352) / '102.0              !maximum density [mol/L]' /
DATA hydrogen_normal(353) / '7   4              !# terms for dilute gas function:  numerator, denominator' /
DATA hydrogen_normal(354) / ' 1.0    1.0d-03    !reducing parameters for T, tcx' /
DATA hydrogen_normal(355) / '-1.24159d+007     0.00E+00' /
DATA hydrogen_normal(356) / ' 5.04056d+006     0.10E+01' /
DATA hydrogen_normal(357) / '-4.80868d+004     0.20E+01' /
DATA hydrogen_normal(358) / ' 3.26394d+002     0.30E+01' /
DATA hydrogen_normal(359) / ' 9.56218d-002     0.40E+01' /
DATA hydrogen_normal(360) / ' 1.73488d-004     0.50e+01' /
DATA hydrogen_normal(361) / '-3.12802d-008     0.60E+01' /
DATA hydrogen_normal(362) / ' 5.04305d+006     0.00E+00' /
DATA hydrogen_normal(363) / '-2.43753d+004     0.10E+01' /
DATA hydrogen_normal(364) / ' 1.51523d+002     0.20E+01' /
DATA hydrogen_normal(365) / ' 1.00000d+000     0.30E+01' /
DATA hydrogen_normal(366) / '10  0               !# terms for background gas function:  numerator, denominator' /
DATA hydrogen_normal(367) / '33.145  15.508  1.0                          !reducing par for T, rho, tcx' /
DATA hydrogen_normal(368) / '  .36308100E-01     .00E+00     .10E+01     .00E+00' /
DATA hydrogen_normal(369) / ' -.20762900E-01     .00E+00     .20E+01     .00E+00' /
DATA hydrogen_normal(370) / '  .31481000E-01     .00E+00     .30E+01     .00E+00' /
DATA hydrogen_normal(371) / ' -.14309700E-01     .00E+00     .40E+01     .00E+00' /
DATA hydrogen_normal(372) / '  .17498000E-02     .00E+00     .50E+01     .00E+00' /
DATA hydrogen_normal(373) / '  .18337000E-02     .10E+01     .10E+01     .00E+00' /
DATA hydrogen_normal(374) / ' -.88671600E-02     .10E+01     .20E+01     .00E+00' /
DATA hydrogen_normal(375) / '  .15826000E-01     .10E+01     .30E+01     .00E+00' /
DATA hydrogen_normal(376) / ' -.10628300E-01     .10E+01     .40E+01     .00E+00' /
DATA hydrogen_normal(377) / '  .28067300E-02     .10E+01     .50E+01     .00E+00' /
DATA hydrogen_normal(378) / 'TK3                !pointer to critical enhancement auxiliary function' /
DATA hydrogen_normal(379) / '' /
DATA hydrogen_normal(380) / '' /
DATA hydrogen_normal(381) / '#AUX               !thermal conductivity critical enhancement model' /
DATA hydrogen_normal(382) / 'TK3  thermal conductivity critical enhancement of Assael et al. (2011).' /
DATA hydrogen_normal(383) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(384) / '?\' /
DATA hydrogen_normal(385) / '!end of info section' /
DATA hydrogen_normal(386) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(387) / '1000.0             !upper temperature limit [K]' /
DATA hydrogen_normal(388) / '2000000.0          !upper pressure limit [kPa]' /
DATA hydrogen_normal(389) / '102.0              !maximum density [mol/L]' /
DATA hydrogen_normal(390) / '9  0  0  0         !# terms:  terms, spare, spare, spare' /
DATA hydrogen_normal(391) / '1.0    1.0  1.0    !reducing par for T, rho, tcx (mW/m-K)' /
DATA hydrogen_normal(392) / ' 0.630d+00         !gnu (universal exponent)' /
DATA hydrogen_normal(393) / ' 1.2415d+00        !gamma (universal exponent)' /
DATA hydrogen_normal(394) / ' 1.01d+00          !R0 (universal amplitude)' /
DATA hydrogen_normal(395) / ' 0.065d+00         !z (universal exponent--not used for t.c., only viscosity)' /
DATA hydrogen_normal(396) / ' 1.00d+00          !c (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA hydrogen_normal(397) / ' 1.5d-10           !xi0 (amplitude) [m]' /
DATA hydrogen_normal(398) / ' 0.052d+00         !gam0 (amplitude) [-]' /
DATA hydrogen_normal(399) / ' 0.40d-09          !qd_inverse (modified effective cutoff parameter) [m]' /
DATA hydrogen_normal(400) / '49.7175d0          !tref (reference temperature) [K]' /
DATA hydrogen_normal(401) / '' /
DATA hydrogen_normal(402) / '' /
DATA hydrogen_normal(403) / '' /
DATA hydrogen_normal(404) / '@ETA               !viscosity model specification' /
DATA hydrogen_normal(405) / 'VS1  pure fluid viscosity model of Vargaftik et al. (1996).' /
DATA hydrogen_normal(406) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(407) / '?Vargaftik, N.B., Vinogradov, Y.K. and Yargin, V.S.' /
DATA hydrogen_normal(408) / '? "Handbook of Physical Properties of Liquids and Gases", Begell House, NY 1996' /
DATA hydrogen_normal(409) / '?' /
DATA hydrogen_normal(410) / '?\' /
DATA hydrogen_normal(411) / '!end of info section' /
DATA hydrogen_normal(412) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(413) / '2000.0             !upper temperature limit [K]' /
DATA hydrogen_normal(414) / '2000000.0          !upper pressure limit [kPa]' /
DATA hydrogen_normal(415) / '102.0              !maximum density [mol/L]' /
DATA hydrogen_normal(416) / '9                  !number of terms associated with dilute-gas function' /
DATA hydrogen_normal(417) / 'NUL                !pointer to reduced effective collision cross-section model' /
DATA hydrogen_normal(418) / '0.2827             !Lennard-Jones coefficient sigma [nm];not used' /
DATA hydrogen_normal(419) / '59.7               !Lennard-Jones coefficient epsilon/kappa [K];not used' /
DATA hydrogen_normal(420) / '32.938d0     1.0d0 !reducing parameters for T, eta' /
DATA hydrogen_normal(421) / '0.0d0        0.5d0 !Chapman-Enskog term; not used here' /
DATA hydrogen_normal(422) / '-2.1505d-1   -1.5d0' /
DATA hydrogen_normal(423) / '10.727d-1    -1.0d0' /
DATA hydrogen_normal(424) / '-16.935d-1   -0.5d0' /
DATA hydrogen_normal(425) / '0.0d0         0.0d0' /
DATA hydrogen_normal(426) / '22.702d-1     0.5d0' /
DATA hydrogen_normal(427) / '2.2123d-1     1.0d0' /
DATA hydrogen_normal(428) / '0.34163d-1    1.5d0' /
DATA hydrogen_normal(429) / '-0.043206d-1  2.0d0' /
DATA hydrogen_normal(430) / '0                  !number of terms for initial density dependence' /
DATA hydrogen_normal(431) / '0 12 0 0 0 0       !# resid terms:  close-packed density;  simple poly; numerator of rational poly; denominator of rat. poly; numerator of exponential; denominator of exponential' /
DATA hydrogen_normal(432) / '32.938d0 15.556d0 1.0d0                   !reducing parameters for T (= eps/k), rho, eta' /
DATA hydrogen_normal(433) / '-9.22703d-1   0.00  1.00  0.00  0  ! powers of tau, del, del0; power of del in exponential [0 indicated no exponential term present]' /
DATA hydrogen_normal(434) / ' 6.41602d0   -1.00  1.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(435) / '-5.98018d0   -2.00  1.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(436) / ' 2.89715d-1  -3.00  1.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(437) / ' 2.36429d0    0.00  2.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(438) / '-2.78870d-1   0.00  3.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(439) / '-1.10595d1   -1.00  3.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(440) / ' 1.11582d1   -2.00  3.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(441) / ' 7.18928d0   -1.00  4.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(442) / '-7.76971d0   -2.00  4.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(443) / '-1.21827d0   -1.00  5.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(444) / ' 1.47193d0   -2.00  5.00  0.00  0  ! powers of tau, del, del0' /
DATA hydrogen_normal(445) / 'NUL                !pointer to critical enhancement auxiliary' /
DATA hydrogen_normal(446) / '' /
DATA hydrogen_normal(447) / '' /
DATA hydrogen_normal(448) / '#STN        !surface tension specification' /
DATA hydrogen_normal(449) / 'ST1  surface tension model of Mulero et al. (2012)' /
DATA hydrogen_normal(450) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(451) / '?Mulero, A., Cachadina, I., and Parra, M.I.' /
DATA hydrogen_normal(452) / '? "Recommended Correlations for the Surface Tension of Common Fluids,"' /
DATA hydrogen_normal(453) / '? J. Phys. Chem. Ref. Data, 41, 043105, 2012.' /
DATA hydrogen_normal(454) / '?\' /
DATA hydrogen_normal(455) / '!end of info section' /
DATA hydrogen_normal(456) / '0.0                !lower temperature limit [K]' /
DATA hydrogen_normal(457) / '33.145             !upper temperature limit [K]' /
DATA hydrogen_normal(458) / '0.0                !(dummy) upper pressure limit' /
DATA hydrogen_normal(459) / '0.0                !(dummy) maximum density' /
DATA hydrogen_normal(460) / '3                           !number of terms in surface tension model' /
DATA hydrogen_normal(461) / '33.145                      !critical temperature used in fit (dummy)' /
DATA hydrogen_normal(462) / '-1.4165      0.63882        !sigma0 and n' /
DATA hydrogen_normal(463) / ' 0.746383    0.659804' /
DATA hydrogen_normal(464) / ' 0.675625    0.619149' /
DATA hydrogen_normal(465) / '' /
DATA hydrogen_normal(466) / '' /
DATA hydrogen_normal(467) / '#DE         !dielectric constant specification' /
DATA hydrogen_normal(468) / 'DE3  dielectric constant model of Harvey and Lemmon (2005).' /
DATA hydrogen_normal(469) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(470) / '?Harvey, A.H. and Lemmon, E.W.' /
DATA hydrogen_normal(471) / '? "Method for Estimating the Dielectric Constant of Natural Gas Mixtures,"' /
DATA hydrogen_normal(472) / '? Int. J. Thermophys., 26(1):31-46, 2005.' /
DATA hydrogen_normal(473) / '?\' /
DATA hydrogen_normal(474) / '!end of info section' /
DATA hydrogen_normal(475) / '0.0                !lower temperature limit [K]' /
DATA hydrogen_normal(476) / '2000.0             !upper temperature limit [K]' /
DATA hydrogen_normal(477) / '0.0                !(dummy) upper pressure limit' /
DATA hydrogen_normal(478) / '0.0                !(dummy) maximum density' /
DATA hydrogen_normal(479) / '273.16 1000.0 1.0  !reducing parameters for t and d' /
DATA hydrogen_normal(480) / '0 2 3 0 0 0                         !number of terms in dielectric constant model' /
DATA hydrogen_normal(481) / ' 2.0306           0.    1.    0.    !coef, t exp, d exp' /
DATA hydrogen_normal(482) / ' 0.0056           1.    1.    0.    !coef, t exp, d exp' /
DATA hydrogen_normal(483) / ' 0.181            0.    2.    0.' /
DATA hydrogen_normal(484) / ' 0.021            1.    2.    0.' /
DATA hydrogen_normal(485) / '-7.4              0.    3.    0.' /
DATA hydrogen_normal(486) / '' /
DATA hydrogen_normal(487) / '' /
DATA hydrogen_normal(488) / '#MLT        !melting line specification' /
DATA hydrogen_normal(489) / 'ML1  melting line model' /
DATA hydrogen_normal(490) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(491) / '?preliminary equation, 2007.' /
DATA hydrogen_normal(492) / '?\' /
DATA hydrogen_normal(493) / '!end of info section' /
DATA hydrogen_normal(494) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(495) / '400.0              !upper temperature limit [K]' /
DATA hydrogen_normal(496) / '0.0                !(dummy) upper pressure limit' /
DATA hydrogen_normal(497) / '0.0                !(dummy) maximum density' /
DATA hydrogen_normal(498) / '13.957   7.3578    !reducing temperature and pressure' /
DATA hydrogen_normal(499) / '1 2 0 0 0 0                 !number of terms in melting line equation' /
DATA hydrogen_normal(500) / ' 1.       0.0               !coefficients and exponents' /
DATA hydrogen_normal(501) / ' 5626.3   1.0' /
DATA hydrogen_normal(502) / ' 2717.2   1.83' /
DATA hydrogen_normal(503) / '' /
DATA hydrogen_normal(504) / '' /
DATA hydrogen_normal(505) / '#SBL        !sublimation line specification' /
DATA hydrogen_normal(506) / 'SB3  sublimation line model of Lemmon (2003).' /
DATA hydrogen_normal(507) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(508) / '?Lemmon, E.W., 2003.' /
DATA hydrogen_normal(509) / '?\' /
DATA hydrogen_normal(510) / '!end of info section' /
DATA hydrogen_normal(511) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(512) / '13.957             !upper temperature limit [K]' /
DATA hydrogen_normal(513) / '0.0                !(dummy) upper pressure limit' /
DATA hydrogen_normal(514) / '0.0                !(dummy) maximum density' /
DATA hydrogen_normal(515) / '13.957  7.70       !reducing temperature and pressure' /
DATA hydrogen_normal(516) / '0 1 0 0 0 0                 !number of terms in sublimation line equation' /
DATA hydrogen_normal(517) / '-8.065      0.93            !coefficients and exponents' /
DATA hydrogen_normal(518) / '' /
DATA hydrogen_normal(519) / '' /
DATA hydrogen_normal(520) / '#PS         !vapor pressure equation' /
DATA hydrogen_normal(521) / 'PS5  vapor pressure equation' /
DATA hydrogen_normal(522) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(523) / '?Leachman, J.W., Jacobsen, R.T, Penoncello, S.G., Lemmon, E.W.' /
DATA hydrogen_normal(524) / '?"Fundamental Equations of State for Parahydrogen, Normal Hydrogen, and' /
DATA hydrogen_normal(525) / '?Orthohydrogen,"' /
DATA hydrogen_normal(526) / '?J. Phys. Chem. Ref. Data, 38(3):721-748, 2009.' /
DATA hydrogen_normal(527) / '?\' /
DATA hydrogen_normal(528) / '!end of info section' /
DATA hydrogen_normal(529) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(530) / '33.145             !upper temperature limit [K]' /
DATA hydrogen_normal(531) / '0.0                !(dummy) upper pressure limit' /
DATA hydrogen_normal(532) / '0.0                !(dummy) maximum density' /
DATA hydrogen_normal(533) / '33.145    1296.4   !reducing parameters' /
DATA hydrogen_normal(534) / ' 4 0 0 0 0 0       !number of terms in equation' /
DATA hydrogen_normal(535) / '-0.489789D+01      1.0' /
DATA hydrogen_normal(536) / ' 0.988558D+00      1.5' /
DATA hydrogen_normal(537) / ' 0.349689D+00      2.0' /
DATA hydrogen_normal(538) / ' 0.499356D+00      2.85' /
DATA hydrogen_normal(539) / '' /
DATA hydrogen_normal(540) / '' /
DATA hydrogen_normal(541) / '#DL         !saturated liquid density equation' /
DATA hydrogen_normal(542) / 'DL1  saturated liquid density equation of Lemmon (2010).' /
DATA hydrogen_normal(543) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(544) / '?Lemmon, C.K. and Lemmon, E.W., 2010.' /
DATA hydrogen_normal(545) / '?\' /
DATA hydrogen_normal(546) / '!end of info section' /
DATA hydrogen_normal(547) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(548) / '33.145             !upper temperature limit [K]' /
DATA hydrogen_normal(549) / '0.0                !(dummy) upper pressure limit' /
DATA hydrogen_normal(550) / '0.0                !(dummy) maximum density' /
DATA hydrogen_normal(551) / '33.145  15.508     !reducing parameters' /
DATA hydrogen_normal(552) / '5 0 0 0 0 0        !number of terms in equation' /
DATA hydrogen_normal(553) / ' 0.15456D+02   0.62          !coefficients and exponents' /
DATA hydrogen_normal(554) / '-0.41720D+02   0.83' /
DATA hydrogen_normal(555) / ' 0.50276D+02   1.05' /
DATA hydrogen_normal(556) / '-0.27947D+02   1.3' /
DATA hydrogen_normal(557) / ' 0.56718D+01   1.6' /
DATA hydrogen_normal(558) / '' /
DATA hydrogen_normal(559) / '' /
DATA hydrogen_normal(560) / '#DV         !saturated vapor density equation' /
DATA hydrogen_normal(561) / 'DV3  saturated vapor density equation of Lemmon (2010).' /
DATA hydrogen_normal(562) / '?LITERATURE REFERENCE \' /
DATA hydrogen_normal(563) / '?Lemmon, C.K. and Lemmon, E.W., 2010.' /
DATA hydrogen_normal(564) / '?\' /
DATA hydrogen_normal(565) / '!end of info section' /
DATA hydrogen_normal(566) / '13.957             !lower temperature limit [K]' /
DATA hydrogen_normal(567) / '33.145             !upper temperature limit [K]' /
DATA hydrogen_normal(568) / '0.0                !(dummy) upper pressure limit' /
DATA hydrogen_normal(569) / '0.0                !(dummy) maximum density' /
DATA hydrogen_normal(570) / '33.145  15.508     !reducing parameters' /
DATA hydrogen_normal(571) / '6 0 0 0 0 0        !number of terms in equation' /
DATA hydrogen_normal(572) / '-0.29962D+01   0.466         !coefficients and exponents' /
DATA hydrogen_normal(573) / '-0.16724D+02   2.' /
DATA hydrogen_normal(574) / ' 0.15819D+02   2.4' /
DATA hydrogen_normal(575) / '-0.16852D+02   4.' /
DATA hydrogen_normal(576) / ' 0.34586D+02   7.' /
DATA hydrogen_normal(577) / '-0.53754D+02   8.' /
DATA hydrogen_normal(578) / '' /
DATA hydrogen_normal(579) / '' /
DATA hydrogen_normal(580) / '@END' /
DATA hydrogen_normal(581) / 'c        1         2         3         4         5         6         7         8' /
DATA hydrogen_normal(582) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
! #######################################################
character(256), TARGET :: Monoethanolamine(321)
DATA Monoethanolamine(1) / 'Monoethanolamine     !Short name' /
DATA Monoethanolamine(2) / '141-43-5             !CAS number' /
DATA Monoethanolamine(3) / 'Ethanolamine         !Full name' /
DATA Monoethanolamine(4) / 'HOCH2CH2NH2          !Chemical formula {C2H7NO}' /
DATA Monoethanolamine(5) / '2-Aminoethanol       !Synonym' /
DATA Monoethanolamine(6) / '61.0831              !Molar mass [g/mol]' /
DATA Monoethanolamine(7) / '283.7                !Triple point temperature [K]' /
DATA Monoethanolamine(8) / '443.564              !Normal boiling point [K]' /
DATA Monoethanolamine(9) / '671.4                !Critical temperature [K]' /
DATA Monoethanolamine(10) / '8125.0               !Critical pressure [kPa]' /
DATA Monoethanolamine(11) / '5.39                 !Critical density [mol/L]' /
DATA Monoethanolamine(12) / '0.573                !Acentric factor' /
DATA Monoethanolamine(13) / '2.36992              !Dipole moment [Debye]  Ikada, E., Hida, Y., Okamoto, H., Hagino, J., Koizumi, N., "Dielectric Properties of Ethanolamines, " Bull. Inst. Chem. Res., Kyoto Univ., 46, 5, 239-247 (1969).' /
DATA Monoethanolamine(14) / 'NBP                  !Default reference state' /
DATA Monoethanolamine(15) / '10.0                 !Version number' /
DATA Monoethanolamine(16) / '2491                 !UN Number                                                 :UN:' /
DATA Monoethanolamine(17) / 'other                !Family                                                    :Family:' /
DATA Monoethanolamine(18) / '????                 !Heating value (upper) [kJ/mol]                            :Heat:' /
DATA Monoethanolamine(19) / '1S/C2H7NO/c3-1-2-4/h4H,1-3H2              !Standard InChI String                :InChi:' /
DATA Monoethanolamine(20) / 'HZAXFHJVJLSVMW-UHFFFAOYSA-N               !Standard InChI Key                   :InChiKey:' /
DATA Monoethanolamine(21) / '????                                      !Alternative fluid for mixing rules   :AltID:' /
DATA Monoethanolamine(22) / '27b11730                                  !Hash number from InChI Key           :Hash:' /
DATA Monoethanolamine(23) / '' /
DATA Monoethanolamine(24) / '' /
DATA Monoethanolamine(25) / '!The fluid files contain general information about the fluid in the first 15 to 20 lines, followed by sections for the' /
DATA Monoethanolamine(26) / '! equations of state, transport equations, and auxiliary equations.  Equations of state are listed first.  The NIST recommended' /
DATA Monoethanolamine(27) / '! equations begin with a hash mark (#).  The secondary equations begin with the @ symbol.  These symbols can be swapped to' /
DATA Monoethanolamine(28) / '! select a secondary equation as primary and the primary as secondary.  The equation of state section also contains auxiliary' /
DATA Monoethanolamine(29) / '! equations for the ideal gas heat capacity or ideal gas Helmholtz energy.  Below the equations of state (both primary and' /
DATA Monoethanolamine(30) / '! secondary) are the transport equations, first viscosity and then thermal conductivity.  These are then followed by the' /
DATA Monoethanolamine(31) / '! secondary equations if available.  The transport section also contains auxiliary equations required to calculate either the' /
DATA Monoethanolamine(32) / '! dilute gas state or the critical enhancement.  At the end of the file are additional but not necessary auxiliary equations,' /
DATA Monoethanolamine(33) / '! including simple equations for the vapor pressure, saturated liquid and vapor densities, melting line (for some fluids), and' /
DATA Monoethanolamine(34) / '! sublimation line (for even fewer fluids).  This section also contains the equations for dielectric constant and surface' /
DATA Monoethanolamine(35) / '! tension if available.  The sections are divided by different symbols (these being _-+=^*~) to aid the eye in locating a' /
DATA Monoethanolamine(36) / '! particular section.  Secondary equations are indented 10 spaces to avoid confusion with the NIST recommended equations.  The' /
DATA Monoethanolamine(37) / '! end of the fluid file is marked with @END.  Anything below that is ignored.' /
DATA Monoethanolamine(38) / '' /
DATA Monoethanolamine(39) / '' /
DATA Monoethanolamine(40) / '! compiled by S. Herrig, Thermodynamics, Ruhr-Universitaet Bochum, Germany' /
DATA Monoethanolamine(41) / '! 03-07-18  SH, Original version.' /
DATA Monoethanolamine(42) / '! 03-09-18 MLH, Add dipole moment, preliminary transport.' /
DATA Monoethanolamine(43) / '' /
DATA Monoethanolamine(44) / '' /
DATA Monoethanolamine(45) / '' /
DATA Monoethanolamine(46) / '' /
DATA Monoethanolamine(47) / '________________________________________________________________________________' /
DATA Monoethanolamine(48) / '' /
DATA Monoethanolamine(49) / '#EOS   !---Equation of state---' /
DATA Monoethanolamine(50) / 'FEQ    !Helmholtz equation of state for monoethanolamine of Herrig et al. (2018).' /
DATA Monoethanolamine(51) / ':TRUECRITICALPOINT:  671.4      5.39          !True EOS critical point [K, mol/L] (where dP/dD=0 and d^2P/dD^2=0 at constant T)' /
DATA Monoethanolamine(52) / ':DOI:' /
DATA Monoethanolamine(53) / '?' /
DATA Monoethanolamine(54) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(55) / '?Herrig, S., Thol, M., and Span, R.,' /
DATA Monoethanolamine(56) / '? unpublished equation, 2018.' /
DATA Monoethanolamine(57) / '?' /
DATA Monoethanolamine(58) / '?The experimental database that was available to fit the EOS was limited to' /
DATA Monoethanolamine(59) / '? measurements in the liquid phase at atmospheric pressure.  At these conditions' /
DATA Monoethanolamine(60) / '? and at temperatures up 430 K, the estimated uncertainty of calculated' /
DATA Monoethanolamine(61) / '? homogeneous densities is 0.3 %. The uncertainty of calculated speed-of-sound' /
DATA Monoethanolamine(62) / '? data is within 0.25 % at temperatures between 290 K and 325 K. Calculated vapor' /
DATA Monoethanolamine(63) / '? pressures are found to be accurate within 5 % for temperatures below 355 K and' /
DATA Monoethanolamine(64) / '? within 1.5 % between 355 K and 445 K. The uncertainties are higher for' /
DATA Monoethanolamine(65) / '? increasing temperatures - there are no reliable data sets to validate the' /
DATA Monoethanolamine(66) / '? equation.  The limited experimental data for isobaric heat capacities between' /
DATA Monoethanolamine(67) / '? 300 K and 355 K deviate by about 5 % from the EOS.  Since its extrapolation' /
DATA Monoethanolamine(68) / '? behavior was carefully constrained, the EOS will also give qualitatively' /
DATA Monoethanolamine(69) / '? reasonable results beyond the experimentally covered regions.' /
DATA Monoethanolamine(70) / '?' /
DATA Monoethanolamine(71) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(72) / '283.7              !Lower temperature limit [K]' /
DATA Monoethanolamine(73) / '675.0              !Upper temperature limit [K]' /
DATA Monoethanolamine(74) / '9000.0             !Upper pressure limit [kPa]' /
DATA Monoethanolamine(75) / '16.76              !Maximum density [mol/L]' /
DATA Monoethanolamine(76) / 'CPP                                    !Pointer to Cp0 model' /
DATA Monoethanolamine(77) / '61.0831                                !Molar mass [g/mol]' /
DATA Monoethanolamine(78) / '283.7                                  !Triple point temperature [K]' /
DATA Monoethanolamine(79) / '0.015907                               !Pressure at triple point [kPa]' /
DATA Monoethanolamine(80) / '16.76                                  !Density at triple point [mol/L]' /
DATA Monoethanolamine(81) / '443.564                                !Normal boiling point temperature [K]' /
DATA Monoethanolamine(82) / '0.573                                  !Acentric factor' /
DATA Monoethanolamine(83) / '671.4         8125.0       5.39        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Monoethanolamine(84) / '671.4                      5.39        !Reducing parameters [K, mol/L]' /
DATA Monoethanolamine(85) / '8.3144598                              !Gas constant [J/mol-K]' /
DATA Monoethanolamine(86) / '  10  4   4 12   0 0    0 0  0 0  0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Monoethanolamine(87) / '  0.034371657  1.0     4.  0.          !a(i),t(i),d(i),l(i)' /
DATA Monoethanolamine(88) / '  2.804815     0.53    1.  0.' /
DATA Monoethanolamine(89) / ' -3.5328022    1.146   1.  0.' /
DATA Monoethanolamine(90) / ' -0.26052106   0.95    2.  0.' /
DATA Monoethanolamine(91) / '  0.073728099  0.35    3.  0.' /
DATA Monoethanolamine(92) / ' -0.9232864    1.47    1.  2.' /
DATA Monoethanolamine(93) / ' -0.15243636   2.8     3.  2.' /
DATA Monoethanolamine(94) / '  0.44837938   0.9     2.  1.' /
DATA Monoethanolamine(95) / ' -0.17517565   3.0     2.  2.' /
DATA Monoethanolamine(96) / ' -0.012936362  0.83    7.  1.' /
DATA Monoethanolamine(97) / '  1.0823719    1.03    1.  2. 2.   -0.71     -1.82     1.04     0.84      0. 0. 0.' /
DATA Monoethanolamine(98) / ' -0.56755523   0.76    1.  2. 2.   -1.16     -1.5      1.04     0.77      0. 0. 0.' /
DATA Monoethanolamine(99) / ' -0.38808402   0.7     3.  2. 2.   -0.733    -1.74     1.04     0.6       0. 0. 0.' /
DATA Monoethanolamine(100) / ' -6.7388446    1.04    3.  2. 2.   -4.08    -57.       1.37     0.59      0. 0. 0.' /
DATA Monoethanolamine(101) / '                                      eta      beta    gamma   epsilon' /
DATA Monoethanolamine(102) / '                                   EXP[eta*(delta-epsilon)^2+beta*(tau-gamma)^2]' /
DATA Monoethanolamine(103) / '' /
DATA Monoethanolamine(104) / '' /
DATA Monoethanolamine(105) / '#AUX   !---Auxiliary function for Cp0' /
DATA Monoethanolamine(106) / 'CPP    !Ideal gas heat capacity function for monoethanolamine of Herrig et al. (2018).' /
DATA Monoethanolamine(107) / '?' /
DATA Monoethanolamine(108) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(109) / '?Herrig, S., Thol, M., and Span, R., 2018.' /
DATA Monoethanolamine(110) / '?' /
DATA Monoethanolamine(111) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(112) / '0.                 !' /
DATA Monoethanolamine(113) / '10000.             !' /
DATA Monoethanolamine(114) / '0.                 !' /
DATA Monoethanolamine(115) / '0.                 !' /
DATA Monoethanolamine(116) / '1.0     8.3144598  !Reducing parameters for T, Cp0' /
DATA Monoethanolamine(117) / '1 2   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Monoethanolamine(118) / '4.0       0.0' /
DATA Monoethanolamine(119) / '13.7    970.0' /
DATA Monoethanolamine(120) / '11.1   3380.0' /
DATA Monoethanolamine(121) / '' /
DATA Monoethanolamine(122) / '' /
DATA Monoethanolamine(123) / '#AUX   !---Auxiliary function for PX0' /
DATA Monoethanolamine(124) / 'PX0    !Helmholtz energy ideal-gas function for monoethanolamine of Herrig et al. (2018).' /
DATA Monoethanolamine(125) / '?' /
DATA Monoethanolamine(126) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(127) / '?Herrig, S., Thol, M., and Span, R., 2018.' /
DATA Monoethanolamine(128) / '?' /
DATA Monoethanolamine(129) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(130) / '1 2  2  0 0  0 0 0               !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau))' /
DATA Monoethanolamine(131) / '  3.0                   1.0      !ai, ti for [ai*log(tau**ti)] terms' /
DATA Monoethanolamine(132) / ' -1.0371130462264091    0.0      !aj, ti for [ai*tau**ti] terms' /
DATA Monoethanolamine(133) / '  3.7839413217629914    1.0      !aj, ti for [ai*tau**ti] terms' /
DATA Monoethanolamine(134) / ' 13.7    970.0                   !aj, ti for [ai*log(1-exp(-ti/T)] terms' /
DATA Monoethanolamine(135) / ' 11.1   3380.0' /
DATA Monoethanolamine(136) / '' /
DATA Monoethanolamine(137) / '' /
DATA Monoethanolamine(138) / '' /
DATA Monoethanolamine(139) / '' /
DATA Monoethanolamine(140) / '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' /
DATA Monoethanolamine(141) / '' /
DATA Monoethanolamine(142) / '#TRN   !---ECS Transport---' /
DATA Monoethanolamine(143) / 'ECS    !Extended Corresponding States model (Propane reference)' /
DATA Monoethanolamine(144) / ':DOI: 10.6028/NIST.IR.8209' /
DATA Monoethanolamine(145) / '?' /
DATA Monoethanolamine(146) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(147) / '?Huber, M.L., "Models for the Viscosity, Thermal Conductivity, and Surface Tension' /
DATA Monoethanolamine(148) / '? of Selected Pure Fluids as Implemented in REFPROP v10.0," NISTIR 8209, 2018.' /
DATA Monoethanolamine(149) / '? doi: 10.6028/NIST.IR.8209' /
DATA Monoethanolamine(150) / '?' /
DATA Monoethanolamine(151) / '?VISCOSITY' /
DATA Monoethanolamine(152) / '? ECS parameters were based on fitting the data of:' /
DATA Monoethanolamine(153) / '? DiGuilio, R.M., Lee, R.J., Schaeffer, S.T., Brasher, L.L., and Teja, A.S., J. Chem. Eng. Data, 37:239-242, 1992.' /
DATA Monoethanolamine(154) / '? Song, J.-H., Park, S.-B., Yoon, J.-H., Lee, H., and Lee, K.-H., "Densities and Viscosities of Monoethanolamine + Ethylene Glycol + Water," J. Chem. Eng. Data, 41:1152-1154, 1996. doi: 10.1021/je9601366' /
DATA Monoethanolamine(155) / '? Blanco, A., Garcia-Abuin, A., Gomez-Diaz, D., Navaza, J., and Villaverde, O., J. Chem. Eng. Data, 58:653-659, 2013.' /
DATA Monoethanolamine(156) / '? Yin, Y., Zhu, C., and Ma, Y., "Volumetric and Viscometric Properties of Binary and Ternary Mixtures of 1-Butyl-3-Methylimidazolium Tetrafluoroborate, Monoethanolamine and Water," J. Chem. Thermodyn., 102:413-428, 2016. doi: 10.1016/j.jct.2016.07.041' /
DATA Monoethanolamine(157) / '? The estimated uncertainty of the viscosity of the liquid phase at atmospheric pressure over the temperature range from 293 K to 424 K is 3%,' /
DATA Monoethanolamine(158) / '? and the estimated uncertainty of the gas phases is 20%.' /
DATA Monoethanolamine(159) / '?' /
DATA Monoethanolamine(160) / '?THERMAL CONDUCTIVITY' /
DATA Monoethanolamine(161) / '? ECS parameters based on fitting data of' /
DATA Monoethanolamine(162) / '? DiGuilio, R.M., McGregor, W.L., and Teja, A.S., J. Chem. Eng. Data, 37:242-245, 1992.' /
DATA Monoethanolamine(163) / '? The estimated uncertainty of the thermal conductivity of the liquid phase at saturation over 298 K - 447 K is 2%; for the gas phase 20%, larger near critical.' /
DATA Monoethanolamine(164) / '?' /
DATA Monoethanolamine(165) / '?The Lennard-Jones parameters were estimated with the method of Chung, T.H., Ajlan, M., Lee, L.L., and Starling, K.E., "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties," Ind. Eng. Chem. Res., 27:671-679, 1988.' /
DATA Monoethanolamine(166) / '?' /
DATA Monoethanolamine(167) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(168) / '283.7              !Lower temperature limit [K]' /
DATA Monoethanolamine(169) / '675.0              !Upper temperature limit [K]' /
DATA Monoethanolamine(170) / '9000.0             !Upper pressure limit [kPa]' /
DATA Monoethanolamine(171) / '16.76              !Maximum density [mol/L]' /
DATA Monoethanolamine(172) / 'FEQ PROPANE.FLD' /
DATA Monoethanolamine(173) / 'VS1                !Model for reference fluid viscosity' /
DATA Monoethanolamine(174) / 'TC1                !Model for reference fluid thermal conductivity' /
DATA Monoethanolamine(175) / 'BIG                !Large molecule identifier' /
DATA Monoethanolamine(176) / '0.88 0. 0. 0.      !Large molecule parameters' /
DATA Monoethanolamine(177) / '1                  !Lennard-Jones flag (0 or 1) (0 => use estimates)' /
DATA Monoethanolamine(178) / '0.4614             !Lennard-Jones coefficient sigma [nm] for ECS method' /
DATA Monoethanolamine(179) / '533.15             !Lennard-Jones coefficient epsilon/kappa [K] for ECS method' /
DATA Monoethanolamine(180) / '1  0  0              !Number of terms in f_int term in Eucken correlation, spare1, spare2' /
DATA Monoethanolamine(181) / ' 0.00132    0. 0. 0. !Coefficient, power of T, spare1, spare2' /
DATA Monoethanolamine(182) / '3  0  0              !Number of terms in psi (visc shape factor): poly,spare1,spare2' /
DATA Monoethanolamine(183) / ' 1.18676    0. 0. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Monoethanolamine(184) / '-0.260695   0. 1. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Monoethanolamine(185) / ' 7.89293e-2 0. 2. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Monoethanolamine(186) / '2  0  0              !Number of terms in chi (t.c. shape factor): poly,spare1,spare2' /
DATA Monoethanolamine(187) / ' 1.61924    0. 0. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Monoethanolamine(188) / '-0.210496   0. 1. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Monoethanolamine(189) / 'TK3                !Pointer to critical enhancement auxiliary function' /
DATA Monoethanolamine(190) / '' /
DATA Monoethanolamine(191) / '' /
DATA Monoethanolamine(192) / '#AUX   !---Auxiliary function for the thermal conductivity critical enhancement' /
DATA Monoethanolamine(193) / 'TK3    !Simplified thermal conductivity critical enhancement for monoethanolamine of Perkins et al. (2013).' /
DATA Monoethanolamine(194) / '?' /
DATA Monoethanolamine(195) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(196) / '?Perkins, R.A., Sengers, J.V., Abdulagatov, I.M., and Huber, M.L.,' /
DATA Monoethanolamine(197) / '? "Simplified Model for the Critical Thermal-Conductivity Enhancement in Molecular Fluids,"' /
DATA Monoethanolamine(198) / '? Int. J. Thermophys., 34(2):191-212, 2013. doi: 10.1007/s10765-013-1409-z' /
DATA Monoethanolamine(199) / '?' /
DATA Monoethanolamine(200) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(201) / '0.                 !' /
DATA Monoethanolamine(202) / '10000.             !' /
DATA Monoethanolamine(203) / '0.                 !' /
DATA Monoethanolamine(204) / '0.                 !' /
DATA Monoethanolamine(205) / '9 0 0 0            !# terms:  CO2-terms, spare, spare, spare' /
DATA Monoethanolamine(206) / '1.0  1.0  1.0      !Reducing parameters for T, rho, tcx [mW/(m-K)]' /
DATA Monoethanolamine(207) / '0.63               !Nu (universal exponent)' /
DATA Monoethanolamine(208) / '1.239              !Gamma (universal exponent)' /
DATA Monoethanolamine(209) / '1.02               !R0 (universal amplitude)' /
DATA Monoethanolamine(210) / '0.063              !Z (universal exponent--not used for t.c., only viscosity)' /
DATA Monoethanolamine(211) / '1.0                !C (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA Monoethanolamine(212) / '0.173e-9           !Xi0 (amplitude) [m]' /
DATA Monoethanolamine(213) / '0.065              !Gam0 (amplitude) [-]' /
DATA Monoethanolamine(214) / '0.559e-9           !Qd_inverse (modified effective cutoff parameter) [m]; arbitrary guess' /
DATA Monoethanolamine(215) / '1007.1             !Tref (reference temperature)=1.5*Tc [K]' /
DATA Monoethanolamine(216) / '' /
DATA Monoethanolamine(217) / '' /
DATA Monoethanolamine(218) / '' /
DATA Monoethanolamine(219) / '' /
DATA Monoethanolamine(220) / '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' /
DATA Monoethanolamine(221) / '' /
DATA Monoethanolamine(222) / '#STN   !---Surface tension---' /
DATA Monoethanolamine(223) / 'ST1    !Surface tension model for monoethanolamine of Huber (2018).' /
DATA Monoethanolamine(224) / ':DOI: 10.6028/NIST.IR.8209' /
DATA Monoethanolamine(225) / '?' /
DATA Monoethanolamine(226) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(227) / '?Huber, M.L., "Models for the Viscosity, Thermal Conductivity, and Surface Tension' /
DATA Monoethanolamine(228) / '? of Selected Pure Fluids as Implemented in REFPROP v10.0," NISTIR 8209, 2018.' /
DATA Monoethanolamine(229) / '? doi: 10.6028/NIST.IR.8209' /
DATA Monoethanolamine(230) / '?' /
DATA Monoethanolamine(231) / '?Fit at NIST to data including:' /
DATA Monoethanolamine(232) / '? Han, J., Jin, J., Eimer, D.A., and Melaaen, M.C., "Density of Water (1) + Monoethanolamine (2) + CO2 (3) from (298.15 to 413.15) K and Surface Tension of Water (1) + Monoethanolamine (2) from (303.15 to 333.15) K," J. Chem. Eng. Data, 57:1095-1103, 2012. doi: 10.1021/je2010038' /
DATA Monoethanolamine(233) / '? Blanco, A., Garcia-Abuin, A., Gomez-Diaz, D., Navaza, J., and Villaverde, O., J. Chem. Eng. Data, 58:653-659, 2013.' /
DATA Monoethanolamine(234) / '? Jayarathna, S.A., Weerasooriya, A., Dayarathna, S., Eimer, D.A., and Melaaen, M.C., "Densities and Surface Tensions of CO2 Loaded Aqueous Monoethanolamine Solutions with r = (0.2 to 0.7) at T = (303.15 to 333.15) K," J. Chem. Eng. Data, 58:986-992, 2013. doi: 10.1021/je301279x' /
DATA Monoethanolamine(235) / '?' /
DATA Monoethanolamine(236) / '?Estimated uncertainty over 293 - 393 K is 2%.' /
DATA Monoethanolamine(237) / '?' /
DATA Monoethanolamine(238) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(239) / '0.                 !' /
DATA Monoethanolamine(240) / '10000.             !' /
DATA Monoethanolamine(241) / '0.                 !' /
DATA Monoethanolamine(242) / '0.                 !' /
DATA Monoethanolamine(243) / '1                  !Number of terms in surface tension model' /
DATA Monoethanolamine(244) / '671.4              !Critical temperature (dummy)' /
DATA Monoethanolamine(245) / '0.0776613 0.801525 !Sigma0 and n' /
DATA Monoethanolamine(246) / '' /
DATA Monoethanolamine(247) / '' /
DATA Monoethanolamine(248) / '#PS    !---Vapor pressure---' /
DATA Monoethanolamine(249) / 'PS5    !Vapor pressure equation for monoethanolamine of Herrig et al. (2018).' /
DATA Monoethanolamine(250) / '?' /
DATA Monoethanolamine(251) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(252) / '?Herrig, S., Thol, M., and Span, R., 2018.' /
DATA Monoethanolamine(253) / '?' /
DATA Monoethanolamine(254) / '?Functional Form:  P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T] where Theta=1-T/Tc, Tc and Pc' /
DATA Monoethanolamine(255) / '? are the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Monoethanolamine(256) / '?' /
DATA Monoethanolamine(257) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(258) / '0.                 !' /
DATA Monoethanolamine(259) / '10000.             !' /
DATA Monoethanolamine(260) / '0.                 !' /
DATA Monoethanolamine(261) / '0.                 !' /
DATA Monoethanolamine(262) / '671.4   8125.0     !Reducing parameters' /
DATA Monoethanolamine(263) / '5 0 0 0 0 0        !Number of terms in equation' /
DATA Monoethanolamine(264) / '-9.5739   1.0      !Coefficients and exponents' /
DATA Monoethanolamine(265) / ' 5.593    1.5' /
DATA Monoethanolamine(266) / '-5.15     1.9' /
DATA Monoethanolamine(267) / '-5.047    3.7' /
DATA Monoethanolamine(268) / '-4.690   13.0' /
DATA Monoethanolamine(269) / '' /
DATA Monoethanolamine(270) / '' /
DATA Monoethanolamine(271) / '#DL    !---Saturated liquid density---' /
DATA Monoethanolamine(272) / 'DL1    !Saturated liquid density equation for monoethanolamine of Herrig et al. (2018).' /
DATA Monoethanolamine(273) / '?' /
DATA Monoethanolamine(274) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(275) / '?Herrig, S., Thol, M., and Span, R., 2018.' /
DATA Monoethanolamine(276) / '?' /
DATA Monoethanolamine(277) / '?Functional Form:  D=Dc*[1+SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Monoethanolamine(278) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Monoethanolamine(279) / '?' /
DATA Monoethanolamine(280) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(281) / '0.                 !' /
DATA Monoethanolamine(282) / '10000.             !' /
DATA Monoethanolamine(283) / '0.                 !' /
DATA Monoethanolamine(284) / '0.                 !' /
DATA Monoethanolamine(285) / '671.4  5.39        !Reducing parameters' /
DATA Monoethanolamine(286) / '6 0 0 0 0 0        !Number of terms in equation' /
DATA Monoethanolamine(287) / '-0.103   0.06      !Coefficients and exponents' /
DATA Monoethanolamine(288) / ' 2.116   0.3' /
DATA Monoethanolamine(289) / ' 3.495   1.5' /
DATA Monoethanolamine(290) / '-4.452   2.0' /
DATA Monoethanolamine(291) / ' 1.795   2.9' /
DATA Monoethanolamine(292) / ' 0.741   19.0' /
DATA Monoethanolamine(293) / '' /
DATA Monoethanolamine(294) / '' /
DATA Monoethanolamine(295) / '#DV    !---Saturated vapor density---' /
DATA Monoethanolamine(296) / 'DV3    !Saturated vapor density equation for monoethanolamine of Herrig et al. (2018).' /
DATA Monoethanolamine(297) / '?' /
DATA Monoethanolamine(298) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(299) / '?Herrig, S., Thol, M., and Span, R., 2018.' /
DATA Monoethanolamine(300) / '?' /
DATA Monoethanolamine(301) / '?Functional Form:  D=Dc*EXP[SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Monoethanolamine(302) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Monoethanolamine(303) / '?' /
DATA Monoethanolamine(304) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Monoethanolamine(305) / '0.                 !' /
DATA Monoethanolamine(306) / '10000.             !' /
DATA Monoethanolamine(307) / '0.                 !' /
DATA Monoethanolamine(308) / '0.                 !' /
DATA Monoethanolamine(309) / '671.4   5.39       !Reducing parameters' /
DATA Monoethanolamine(310) / '6 0 0 0 0 0        !Number of terms in equation' /
DATA Monoethanolamine(311) / ' 0.1     0.07      !Coefficients and exponents' /
DATA Monoethanolamine(312) / '-2.907   0.34' /
DATA Monoethanolamine(313) / '-7.405   1.1' /
DATA Monoethanolamine(314) / '-25.86   3.0' /
DATA Monoethanolamine(315) / '-70.87   6.34' /
DATA Monoethanolamine(316) / '-210.8  14.6' /
DATA Monoethanolamine(317) / '' /
DATA Monoethanolamine(318) / '' /
DATA Monoethanolamine(319) / '@END' /
DATA Monoethanolamine(320) / 'c        1         2         3         4         5         6         7         8' /
DATA Monoethanolamine(321) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
! #######################################################
character(256), TARGET :: Methane(1032)
DATA Methane(1) / 'Methane              !Short name' /
DATA Methane(2) / '74-82-8              !CAS number' /
DATA Methane(3) / 'Methane              !Full name' /
DATA Methane(4) / 'CH4                  !Chemical formula {CH4}' /
DATA Methane(5) / 'R-50                 !Synonym' /
DATA Methane(6) / '16.0428              !Molar mass [g/mol]' /
DATA Methane(7) / '90.6941              !Triple point temperature [K]' /
DATA Methane(8) / '111.667              !Normal boiling point [K]' /
DATA Methane(9) / '190.564              !Critical temperature [K]' /
DATA Methane(10) / '4599.2               !Critical pressure [kPa]' /
DATA Methane(11) / '10.139               !Critical density [mol/L]' /
DATA Methane(12) / '0.01142              !Acentric factor' /
DATA Methane(13) / '0.0                  !Dipole moment [Debye]; (exactly zero due to symmetry)' /
DATA Methane(14) / 'NBP                  !Default reference state' /
DATA Methane(15) / '10.0                 !Version number' /
DATA Methane(16) / '1971, 1972           !UN Number                                                 :UN:' /
DATA Methane(17) / 'n-alkane             !Family                                                    :Family:' /
DATA Methane(18) / '890.58               !Heating value (upper) at 25 C [kJ/mol] (ISO 6976:2016)    :Heat:' /
DATA Methane(19) / '25.                  !GWP (IPCC 2007)                                           :GWP:' /
DATA Methane(20) / 'A3                   !Safety Group (ASHRAE Standard 34, 2010)                   :Safety:' /
DATA Methane(21) / '1S/CH4/h1H4                               !Standard InChI String                :InChi:' /
DATA Methane(22) / 'VNWKTOKETHGBQD-UHFFFAOYSA-N               !Standard InChI Key                   :InChiKey:' /
DATA Methane(23) / '????                                      !Alternative fluid for mixing rules   :AltID:' /
DATA Methane(24) / '8ae7a700                                  !Hash number from InChI Key           :Hash:' /
DATA Methane(25) / '' /
DATA Methane(26) / '' /
DATA Methane(27) / '' /
DATA Methane(28) / '' /
DATA Methane(29) / '!The fluid files contain general information about the fluid in the first 15 to 20 lines, followed by sections for the' /
DATA Methane(30) / '! equations of state, transport equations, and auxiliary equations.  Equations of state are listed first.  The NIST recommended' /
DATA Methane(31) / '! equations begin with a hash mark (#).  The secondary equations begin with the @ symbol.  These symbols can be swapped to' /
DATA Methane(32) / '! select a secondary equation as primary and the primary as secondary.  The equation of state section also contains auxiliary' /
DATA Methane(33) / '! equations for the ideal gas heat capacity or ideal gas Helmholtz energy.  Below the equations of state (both primary and' /
DATA Methane(34) / '! secondary) are the transport equations, first viscosity and then thermal conductivity.  These are then followed by the' /
DATA Methane(35) / '! secondary equations if available.  The transport section also contains auxiliary equations required to calculate either the' /
DATA Methane(36) / '! dilute gas state or the critical enhancement.  At the end of the file are additional but not necessary auxiliary equations,' /
DATA Methane(37) / '! including simple equations for the vapor pressure, saturated liquid and vapor densities, melting line (for some fluids), and' /
DATA Methane(38) / '! sublimation line (for even fewer fluids).  This section also contains the equations for dielectric constant and surface' /
DATA Methane(39) / '! tension if available.  The sections are divided by different symbols (these being _-+=^*~) to aid the eye in locating a' /
DATA Methane(40) / '! particular section.  Secondary equations are indented 10 spaces to avoid confusion with the NIST recommended equations.  The' /
DATA Methane(41) / '! end of the fluid file is marked with @END.  Anything below that is ignored.' /
DATA Methane(42) / '' /
DATA Methane(43) / '' /
DATA Methane(44) / '! compiled by E.W. Lemmon, NIST Physical and Chemical Properties Division, Boulder, Colorado' /
DATA Methane(45) / '! 01-22-97 EWL, Original version.' /
DATA Methane(46) / '! 06-24-98 EWL, Add Younglove and Ely BWR equation of state.' /
DATA Methane(47) / '! 11-18-98 EWL, Add equation of state of Friend et al. (1989).' /
DATA Methane(48) / '! 11-01-99 EWL, Add Span 12 term short equation of state.' /
DATA Methane(49) / '! 01-26-00 EWL, Add Friend transport equations, but keep Younglove viscosity eq. as default since Friend eq. has an anomaly above 100 MPa.' /
DATA Methane(50) / '! 07-23-02 EWL, Add sublimation line.' /
DATA Methane(51) / '! 08-05-04 EWL, Add Harvey and Lemmon dielectric correlation.' /
DATA Methane(52) / '! 10-13-04 MLH, Add family.' /
DATA Methane(53) / '! 07-14-05 MLH, Add Vogel(2000) viscosity correlation.' /
DATA Methane(54) / '! 12-02-06 MLH, Update LJ for ECS.' /
DATA Methane(55) / '! 01-05-07 MLH, Add VS4 model, new VS1 model of Vogel, moved Friend VS1 model to EOF.' /
DATA Methane(56) / '! 03-05-07 EWL, Add ancillary equations.' /
DATA Methane(57) / '! 03-09-07 MLH, Add final FT model coefficients.' /
DATA Methane(58) / '! 02-14-08 MLH, Add TK6 block for ECS for mixture calculations.' /
DATA Methane(59) / '! 09-02-10 MLH, Add new VS4 model for viscosity feb2010 model.' /
DATA Methane(60) / '! 04-11-12 MLH, Add extra blank FT coeff for consistent formatting.' /
DATA Methane(61) / '! 12-06-12 EWL, Add surface tension coefficients of Mulero et al. (2012).' /
DATA Methane(62) / '! 05-15-17 EWL, Change the hard coded CH4 model to the TK7 reverse Polish notation.' /
DATA Methane(63) / '' /
DATA Methane(64) / '' /
DATA Methane(65) / '' /
DATA Methane(66) / '' /
DATA Methane(67) / '________________________________________________________________________________' /
DATA Methane(68) / '' /
DATA Methane(69) / '#EOS   !---Equation of state---' /
DATA Methane(70) / 'FEQ    !Helmholtz equation of state for methane of Setzmann and Wagner (1991).' /
DATA Methane(71) / ':TRUECRITICALPOINT:  190.564   10.139128      !True EOS critical point [K, mol/L] (where dP/dD=0 and d^2P/dD^2=0 at constant T)' /
DATA Methane(72) / ':DOI: 10.1063/1.555898' /
DATA Methane(73) / '?' /
DATA Methane(74) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(75) / '?Setzmann, U. and Wagner, W.,' /
DATA Methane(76) / '? "A New Equation of State and Tables of Thermodynamic Properties for Methane' /
DATA Methane(77) / '? Covering the Range from the Melting Line to 625 K at Pressures up to 1000 MPa,"' /
DATA Methane(78) / '? J. Phys. Chem. Ref. Data, 20(6):1061-1151, 1991.' /
DATA Methane(79) / '?' /
DATA Methane(80) / '?The uncertainties in density are 0.03% for pressures below 12 MPa and' /
DATA Methane(81) / '? temperatures below 350 K and up to 0.07% for pressures less than 50 MPa.' /
DATA Methane(82) / '? For the speed of sound, the uncertainty ranges from 0.03% (in the vapor' /
DATA Methane(83) / '? phase) to 0.3% depending on temperature and pressure.  Heat capacities' /
DATA Methane(84) / '? may be generally calculated within an uncertainty of 1%.' /
DATA Methane(85) / '?' /
DATA Methane(86) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(87) / '90.6941            !Lower temperature limit [K]' /
DATA Methane(88) / '625.0              !Upper temperature limit [K]' /
DATA Methane(89) / '1000000.0          !Upper pressure limit [kPa]' /
DATA Methane(90) / '40.072             !Maximum density [mol/L]' /
DATA Methane(91) / 'CPP                                    !Pointer to Cp0 model' /
DATA Methane(92) / '16.0428                                !Molar mass [g/mol]' /
DATA Methane(93) / '90.6941                                !Triple point temperature [K]' /
DATA Methane(94) / '11.696                                 !Pressure at triple point [kPa]' /
DATA Methane(95) / '28.142                                 !Density at triple point [mol/L]' /
DATA Methane(96) / '111.667                                !Normal boiling point temperature [K]' /
DATA Methane(97) / '0.01142                                !Acentric factor' /
DATA Methane(98) / '190.564       4599.2      10.139128    !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Methane(99) / '190.564                   10.139128    !Reducing parameters [K, mol/L]' /
DATA Methane(100) / '8.31451                                !Gas constant [J/mol-K]' /
DATA Methane(101) / '  36  4   4 12   0 0    0 0  0 0  0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Methane(102) / ' 0.04367901028        -0.5      1.  0. !a(i),t(i),d(i),l(i)' /
DATA Methane(103) / ' 0.6709236199          0.5      1.  0.' /
DATA Methane(104) / '-1.765577859           1.0      1.  0.' /
DATA Methane(105) / ' 0.8582330241          0.5      2.  0.' /
DATA Methane(106) / '-1.206513052           1.0      2.  0.' /
DATA Methane(107) / ' 0.512046722           1.5      2.  0.' /
DATA Methane(108) / '-0.0004000010791       4.5      2.  0.' /
DATA Methane(109) / '-0.01247842423         0.0      3.  0.' /
DATA Methane(110) / ' 0.03100269701         1.0      4.  0.' /
DATA Methane(111) / ' 0.001754748522        3.0      4.  0.' /
DATA Methane(112) / '-0.3171921605e-5       1.0      8.  0.' /
DATA Methane(113) / '-0.224034684e-5        3.0      9.  0.' /
DATA Methane(114) / ' 0.2947056156e-6       3.0     10.  0.' /
DATA Methane(115) / ' 0.1830487909          0.0      1.  1.' /
DATA Methane(116) / ' 0.1511883679          1.0      1.  1.' /
DATA Methane(117) / '-0.4289363877          2.0      1.  1.' /
DATA Methane(118) / ' 0.06894002446         0.0      2.  1.' /
DATA Methane(119) / '-0.01408313996         0.0      4.  1.' /
DATA Methane(120) / '-0.0306305483          2.0      5.  1.' /
DATA Methane(121) / '-0.02969906708         2.0      6.  1.' /
DATA Methane(122) / '-0.01932040831         5.0      1.  2.' /
DATA Methane(123) / '-0.1105739959          5.0      2.  2.' /
DATA Methane(124) / ' 0.09952548995         5.0      3.  2.' /
DATA Methane(125) / ' 0.008548437825        2.0      4.  2.' /
DATA Methane(126) / '-0.06150555662         4.0      4.  2.' /
DATA Methane(127) / '-0.04291792423         12.0     3.  3.' /
DATA Methane(128) / '-0.0181320729          8.0      5.  3.' /
DATA Methane(129) / ' 0.0344590476          10.0     5.  3.' /
DATA Methane(130) / '-0.00238591945         10.0     8.  3.' /
DATA Methane(131) / '-0.01159094939         10.0     2.  4.' /
DATA Methane(132) / ' 0.06641693602         14.0     3.  4.' /
DATA Methane(133) / '-0.0237154959          12.0     4.  4.' /
DATA Methane(134) / '-0.03961624905         18.0     4.  4.' /
DATA Methane(135) / '-0.01387292044         22.0     4.  4.' /
DATA Methane(136) / ' 0.03389489599         18.0     5.  4.' /
DATA Methane(137) / '-0.002927378753        14.0     6.  4.' /
DATA Methane(138) / ' 0.9324799946e-4       2.0      2.  2. 2.     -20.0    -200.0   1.07    1.0      0. 0. 0.' /
DATA Methane(139) / '-6.287171518           0.0      0.  2. 2.     -40.0    -250.0   1.11    1.0      0. 0. 0.' /
DATA Methane(140) / ' 12.71069467           1.0      0.  2. 2.     -40.0    -250.0   1.11    1.0      0. 0. 0.' /
DATA Methane(141) / '-6.423953466           2.0      0.  2. 2.     -40.0    -250.0   1.11    1.0      0. 0. 0.' /
DATA Methane(142) / '                                               eta      beta    gamma   epsilon' /
DATA Methane(143) / '                                          EXP[eta*(delta-epsilon)^2+beta*(tau-gamma)^2]' /
DATA Methane(144) / '' /
DATA Methane(145) / '' /
DATA Methane(146) / '#AUX   !---Auxiliary function for Cp0' /
DATA Methane(147) / 'CPP    !Ideal gas heat capacity function for methane of Setzmann and Wagner (1991).' /
DATA Methane(148) / '?' /
DATA Methane(149) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(150) / '?Setzmann, U. and Wagner, W., 1991.' /
DATA Methane(151) / '?' /
DATA Methane(152) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(153) / '0.                 !' /
DATA Methane(154) / '10000.             !' /
DATA Methane(155) / '0.                 !' /
DATA Methane(156) / '0.                 !' /
DATA Methane(157) / '1.0     8.31451    !Reducing parameters for T, Cp0' /
DATA Methane(158) / '1 5   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Methane(159) / ' 4.0016     0.0' /
DATA Methane(160) / ' 0.008449   648.0' /
DATA Methane(161) / ' 4.6942     1957.0' /
DATA Methane(162) / ' 3.4865     3895.0' /
DATA Methane(163) / ' 1.6572     5705.0' /
DATA Methane(164) / ' 1.4115     15080.0' /
DATA Methane(165) / '' /
DATA Methane(166) / '' /
DATA Methane(167) / '#AUX   !---Auxiliary function for PX0' /
DATA Methane(168) / 'PX0    !Helmholtz energy ideal-gas function for methane of Setzmann and Wagner (1991).' /
DATA Methane(169) / '?' /
DATA Methane(170) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(171) / '?Setzmann, U. and Wagner, W., 1991.' /
DATA Methane(172) / '?' /
DATA Methane(173) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(174) / '1 2  5  0 0  0 0 0               !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau))' /
DATA Methane(175) / '  3.0016                1.0      !ai, ti for [ai*log(tau**ti)] terms' /
DATA Methane(176) / ' -2.9705496667947529    0.0      !aj, ti for [ai*tau**ti] terms' /
DATA Methane(177) / '  2.8907453831087553    1.0      !aj, ti for [ai*tau**ti] terms' /
DATA Methane(178) / '  0.008449   648.0               !aj, ti for [ai*log(1-exp(-ti/T)] terms' /
DATA Methane(179) / '  4.6942     1957.0' /
DATA Methane(180) / '  3.4865     3895.0' /
DATA Methane(181) / '  1.6572     5705.0' /
DATA Methane(182) / '  1.4115     15080.0' /
DATA Methane(183) / '' /
DATA Methane(184) / '' /
DATA Methane(185) / '#AUX   !---Auxiliary function for PH0' /
DATA Methane(186) / 'PH0    !Ideal gas Helmholtz form for methane.' /
DATA Methane(187) / '?' /
DATA Methane(188) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(189) / '?Setzmann, U. and Wagner, W., 1991.' /
DATA Methane(190) / '?' /
DATA Methane(191) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(192) / '0.                 !' /
DATA Methane(193) / '10000.             !' /
DATA Methane(194) / '0.                 !' /
DATA Methane(195) / '0.                 !' /
DATA Methane(196) / '1 2  5  0 0  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA Methane(197) / ' 3.0016            1.0                 !ai, ti for [ai*log(tau**ti)] terms' /
DATA Methane(198) / '-2.9705496668      0.0                 !aj, ti for [ai*tau**ti] terms' /
DATA Methane(199) / ' 2.8907453831      1.0' /
DATA Methane(200) / ' 0.008449         -3.4004324007        !aj, ti for [ai*log(1-exp(ti*tau)] terms' /
DATA Methane(201) / ' 4.6942           -10.2695157532' /
DATA Methane(202) / ' 3.4865           -20.43932747' /
DATA Methane(203) / ' 1.6572           -29.9374488361' /
DATA Methane(204) / ' 1.4115           -79.1335194475' /
DATA Methane(205) / '' /
DATA Methane(206) / '' /
DATA Methane(207) / '' /
DATA Methane(208) / '' /
DATA Methane(209) / '--------------------------------------------------------------------------------' /
DATA Methane(210) / '' /
DATA Methane(211) / '@EOS    !---Equation of state---' /
DATA Methane(212) / 'FEK     !Helmholtz equation of state for methane of Kunz and Wagner (2004).' /
DATA Methane(213) / '          ?' /
DATA Methane(214) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(215) / '          ?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA Methane(216) / '          ? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA Methane(217) / '          ? and Other Mixtures," GERG Technical Monograph 15,' /
DATA Methane(218) / '          ? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA Methane(219) / '          ?' /
DATA Methane(220) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(221) / '          90.6941            !Lower temperature limit [K]' /
DATA Methane(222) / '          625.0              !Upper temperature limit [K]' /
DATA Methane(223) / '          1000000.0          !Upper pressure limit [kPa]' /
DATA Methane(224) / '          40.072             !Maximum density [mol/L]' /
DATA Methane(225) / '          PHK                                    !Pointer to Cp0 model' /
DATA Methane(226) / '          16.04246                               !Molar mass [g/mol]' /
DATA Methane(227) / '          90.6941                                !Triple point temperature [K]' /
DATA Methane(228) / '          11.698                                 !Pressure at triple point [kPa]' /
DATA Methane(229) / '          28.146                                 !Density at triple point [mol/L]' /
DATA Methane(230) / '          111.66                                 !Normal boiling point temperature [K]' /
DATA Methane(231) / '          0.0114                                 !Acentric factor' /
DATA Methane(232) / '          190.564       4599.2      10.139342719 !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Methane(233) / '          190.564                   10.139342719 !Reducing parameters [K, mol/L]' /
DATA Methane(234) / '          8.314472                               !Gas constant [J/mol-K]' /
DATA Methane(235) / '            24  4    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Methane(236) / '           0.57335704239162      0.125     1.   0.' /
DATA Methane(237) / '          -1.6760687523730       1.125     1.   0.' /
DATA Methane(238) / '           0.23405291834916      0.375     2.   0.' /
DATA Methane(239) / '          -0.21947376343441      1.125     2.   0.' /
DATA Methane(240) / '           0.016369201404128     0.625     4.   0.' /
DATA Methane(241) / '           0.01500440638928      1.5       4.   0.' /
DATA Methane(242) / '           0.098990489492918     0.625     1.   1.' /
DATA Methane(243) / '           0.58382770929055      2.625     1.   1.' /
DATA Methane(244) / '          -0.74786867560390      2.75      1.   1.' /
DATA Methane(245) / '           0.30033302857974      2.125     2.   1.' /
DATA Methane(246) / '           0.20985543806568      2.0       3.   1.' /
DATA Methane(247) / '          -0.018590151133061     1.75      6.   1.' /
DATA Methane(248) / '          -0.15782558339049      4.50      2.   2.' /
DATA Methane(249) / '           0.12716735220791      4.75      3.   2.' /
DATA Methane(250) / '          -0.032019743894346     5.0       3.   2.' /
DATA Methane(251) / '          -0.068049729364536     4.0       4.   2.' /
DATA Methane(252) / '           0.024291412853736     4.5       4.   2.' /
DATA Methane(253) / '           0.0051440451639444    7.5       2.   3.' /
DATA Methane(254) / '          -0.019084949733532     14.0      3.   3.' /
DATA Methane(255) / '           0.0055229677241291    11.5      4.   3.' /
DATA Methane(256) / '          -0.0044197392976085    26.0      5.   6.' /
DATA Methane(257) / '           0.040061416708429     28.0      6.   6.' /
DATA Methane(258) / '          -0.033752085907575     30.0      6.   6.' /
DATA Methane(259) / '          -0.0025127658213357    16.0      7.   6.' /
DATA Methane(260) / '' /
DATA Methane(261) / '' /
DATA Methane(262) / '@AUX    !---Auxiliary function for PH0' /
DATA Methane(263) / 'PHK     !Ideal gas Helmholtz form for methane of Kunz and Wagner (2004).' /
DATA Methane(264) / '          ?' /
DATA Methane(265) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(266) / '          ?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA Methane(267) / '          ?' /
DATA Methane(268) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(269) / '          0.                 !' /
DATA Methane(270) / '          10000.             !' /
DATA Methane(271) / '          0.                 !' /
DATA Methane(272) / '          0.                 !' /
DATA Methane(273) / '          1 2  0 2 2  0 0 0  !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA Methane(274) / '           3.00088           1.0                 !ai, ti for [ai*log(tau**ti)] terms' /
DATA Methane(275) / '           19.597508817      0.0                 !aj, ti for [ai*tau**ti] terms' /
DATA Methane(276) / '          -83.959667892      1.0' /
DATA Methane(277) / '          -0.0046            0.936220902         !aj, ti for cosh and sinh terms' /
DATA Methane(278) / '           4.46921           5.722644361' /
DATA Methane(279) / '           0.76315           4.306474465' /
DATA Methane(280) / '           8.74432           5.577233895' /
DATA Methane(281) / '' /
DATA Methane(282) / '' /
DATA Methane(283) / '@EOS    !---Equation of state---' /
DATA Methane(284) / 'FE1     !Helmholtz equation of state for methane of Friend et al. (1989).' /
DATA Methane(285) / '          ?' /
DATA Methane(286) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(287) / '          ?Friend, D.G., Ely, J.F., and Ingham, H.,' /
DATA Methane(288) / '          ? "Thermophysical Properties of Methane,"' /
DATA Methane(289) / '          ? J. Phys. Chem. Ref. Data, 18(2):583-638, 1989.' /
DATA Methane(290) / '          ?' /
DATA Methane(291) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(292) / '          90.6854            !Lower temperature limit [K]' /
DATA Methane(293) / '          620.0              !Upper temperature limit [K]' /
DATA Methane(294) / '          100000.0           !Upper pressure limit [kPa]' /
DATA Methane(295) / '          29.714             !Maximum density [mol/L]' /
DATA Methane(296) / '          CP1                                    !Pointer to Cp0 model' /
DATA Methane(297) / '          16.043                                 !Molar mass [g/mol]' /
DATA Methane(298) / '          90.6854                                !Triple point temperature [K]' /
DATA Methane(299) / '          11.694                                 !Pressure at triple point [kPa]' /
DATA Methane(300) / '          28.145                                 !Density at triple point [mol/L]' /
DATA Methane(301) / '          111.66                                 !Normal boiling point temperature [K]' /
DATA Methane(302) / '          0.0086                                 !Acentric factor' /
DATA Methane(303) / '          190.551       4599.2      10.139       !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Methane(304) / '          190.551                   10.139       !Reducing parameters [K, mol/L]' /
DATA Methane(305) / '          8.31451                                !Gas constant [J/mol-K]' /
DATA Methane(306) / '            32  4    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Methane(307) / '           0.384436099659      0.0       1.  0.  !a(i),t(i),d(i),l(i)' /
DATA Methane(308) / '          -1.796925988         1.5       1.  0.' /
DATA Methane(309) / '           0.329444947369      2.5       1.  0.' /
DATA Methane(310) / '           0.0226312728442    -0.5       2.  0.' /
DATA Methane(311) / '           0.0759236768798     1.5       2.  0.' /
DATA Methane(312) / '           0.0693758447259     2.0       2.  0.' /
DATA Methane(313) / '           0.0241163263947     0.0       3.  0.' /
DATA Methane(314) / '           0.0107009920854     1.0       3.  0.' /
DATA Methane(315) / '          -0.0380933275164     2.5       3.  0.' /
DATA Methane(316) / '           0.000471537561143   0.0       6.  0.' /
DATA Methane(317) / '           0.000556607678805   2.0       7.  0.' /
DATA Methane(318) / '           0.548759346533e-6   5.0       7.  0.' /
DATA Methane(319) / '          -0.999632699967e-4   2.0       8.  0.' /
DATA Methane(320) / '          -0.128087979280      5.0       1.  2.' /
DATA Methane(321) / '           0.0380198873377     6.0       1.  2.' /
DATA Methane(322) / '           0.139226650551      3.5       2.  2.' /
DATA Methane(323) / '          -0.0874996348859     5.5       2.  2.' /
DATA Methane(324) / '          -0.0033489416576     3.0       3.  2.' /
DATA Methane(325) / '          -0.0517576297122     7.0       3.  2.' /
DATA Methane(326) / '           0.0252835179116     6.0       5.  2.' /
DATA Methane(327) / '           0.00051870320595    8.5       6.  2.' /
DATA Methane(328) / '          -0.00166770594525    4.0       7.  2.' /
DATA Methane(329) / '          -0.000607401927389   6.5       8.  2.' /
DATA Methane(330) / '          -0.972915359991e-4   5.5      10.  2.' /
DATA Methane(331) / '          -0.298844010462e-4  22.0       2.  4.' /
DATA Methane(332) / '          -0.0130940111124    11.0       3.  4.' /
DATA Methane(333) / '           0.0198175833798    18.0       3.  4.' /
DATA Methane(334) / '           0.0208465762327    11.0       4.  4.' /
DATA Methane(335) / '          -0.0358025052631    23.0       4.  4.' /
DATA Methane(336) / '          -0.203486851741     17.0       5.  4.' /
DATA Methane(337) / '           0.215964755088     18.0       5.  4.' /
DATA Methane(338) / '          -0.00429340628249   23.0       5.  4.' /
DATA Methane(339) / '' /
DATA Methane(340) / '' /
DATA Methane(341) / '@AUX    !---Auxiliary function for Cp0' /
DATA Methane(342) / 'CP1     !Ideal gas heat capacity function for methane.' /
DATA Methane(343) / '          ?' /
DATA Methane(344) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(345) / '          ?Friend, D.G., Ely, J.F., and Ingham, H.,' /
DATA Methane(346) / '          ?' /
DATA Methane(347) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(348) / '          0.                 !' /
DATA Methane(349) / '          10000.             !' /
DATA Methane(350) / '          0.                 !' /
DATA Methane(351) / '          0.                 !' /
DATA Methane(352) / '          1.0     8.31451    !Reducing parameters for T, Cp0' /
DATA Methane(353) / '          4 1   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Methane(354) / '           3.5998324         0.0' /
DATA Methane(355) / '           0.2614717613495   0.3333333333333' /
DATA Methane(356) / '          -0.05671028952515  0.6666666666667' /
DATA Methane(357) / '           0.004105505612671 1.0' /
DATA Methane(358) / '           4.7206715         2009.15202' /
DATA Methane(359) / '' /
DATA Methane(360) / '' /
DATA Methane(361) / '@EOS    !---Equation of state---' /
DATA Methane(362) / 'FES     !Helmholtz equation of state for methane of Span and Wagner (2003).' /
DATA Methane(363) / '          ?' /
DATA Methane(364) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(365) / '          ?Span, R. and Wagner, W.' /
DATA Methane(366) / '          ? "Equations of State for Technical Applications. II. Results for Nonpolar Fluids,"' /
DATA Methane(367) / '          ? Int. J. Thermophys., 24(1):41-109, 2003. doi: 10.1023/A:1022310214958' /
DATA Methane(368) / '          ?' /
DATA Methane(369) / '          ?The uncertainties of the equation of state are approximately 0.2% (to' /
DATA Methane(370) / '          ? 0.5% at high pressures) in density, 1% (in the vapor phase) to 2% in' /
DATA Methane(371) / '          ? heat capacity, 1% (in the vapor phase) to 2% in the speed of sound, and' /
DATA Methane(372) / '          ? 0.2% in vapor pressure, except in the critical region.' /
DATA Methane(373) / '          ?' /
DATA Methane(374) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(375) / '          90.6941            !Lower temperature limit [K]' /
DATA Methane(376) / '          600.0              !Upper temperature limit [K]' /
DATA Methane(377) / '          100000.0           !Upper pressure limit [kPa]' /
DATA Methane(378) / '          40.072             !Maximum density [mol/L]' /
DATA Methane(379) / '          CPP                                    !Pointer to Cp0 model' /
DATA Methane(380) / '          16.043                                 !Molar mass [g/mol]' /
DATA Methane(381) / '          90.6941                                !Triple point temperature [K]' /
DATA Methane(382) / '          11.661                                 !Pressure at triple point [kPa]' /
DATA Methane(383) / '          28.167                                 !Density at triple point [mol/L]' /
DATA Methane(384) / '          111.66                                 !Normal boiling point temperature [K]' /
DATA Methane(385) / '          0.011                                  !Acentric factor' /
DATA Methane(386) / '          190.564       4599.0      10.139001    !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Methane(387) / '          190.564                   10.139001    !Reducing parameters [K, mol/L]' /
DATA Methane(388) / '          8.31451                                !Gas constant [J/mol-K]' /
DATA Methane(389) / '            12  4    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Methane(390) / '           0.89269676      0.25      1.  0.      !a(i),t(i),d(i),l(i)' /
DATA Methane(391) / '          -2.5438282       1.125     1.  0.' /
DATA Methane(392) / '           0.64980978      1.5       1.  0.' /
DATA Methane(393) / '           0.020793471     1.375     2.  0.' /
DATA Methane(394) / '           0.070189104     0.25      3.  0.' /
DATA Methane(395) / '           0.00023700378   0.875     7.  0.' /
DATA Methane(396) / '           0.16653334      0.625     2.  1.' /
DATA Methane(397) / '          -0.043855669     1.75      5.  1.' /
DATA Methane(398) / '          -0.1572678       3.625     1.  2.' /
DATA Methane(399) / '          -0.035311675     3.625     4.  2.' /
DATA Methane(400) / '          -0.029570024    14.5       3.  3.' /
DATA Methane(401) / '           0.014019842    12.0       4.  3.' /
DATA Methane(402) / '' /
DATA Methane(403) / '' /
DATA Methane(404) / '@EOS    !---Equation of state---' /
DATA Methane(405) / 'BWR     !MBWR equation of state for methane of Younglove and Ely (1987).' /
DATA Methane(406) / '          ?' /
DATA Methane(407) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(408) / '          ?Younglove, B.A. and Ely, J.F.,' /
DATA Methane(409) / '          ? "Thermophysical properties of fluids. II. Methane, ethane, propane,' /
DATA Methane(410) / '          ? isobutane and normal butane,"' /
DATA Methane(411) / '          ? J. Phys. Chem. Ref. Data, 16:577-798, 1987.' /
DATA Methane(412) / '          ? All temperatures on IPTS-68' /
DATA Methane(413) / '          ?' /
DATA Methane(414) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(415) / '          90.68              !Lower temperature limit [K]' /
DATA Methane(416) / '          600.0              !Upper temperature limit [K]' /
DATA Methane(417) / '          200000.0           !Upper pressure limit [kPa]' /
DATA Methane(418) / '          36.2029            !Maximum density [mol/L]' /
DATA Methane(419) / '          CP2                                    !Pointer to Cp0 model' /
DATA Methane(420) / '          16.043                                 !Molar mass [g/mol]' /
DATA Methane(421) / '          90.68                                  !Triple point temperature [K]' /
DATA Methane(422) / '          11.744                                 !Pressure at triple point [kPa]' /
DATA Methane(423) / '          28.147                                 !Density at triple point [mol/L]' /
DATA Methane(424) / '          111.667                                !Normal boiling point temperature [K]' /
DATA Methane(425) / '          0.011                                  !Acentric factor' /
DATA Methane(426) / '          190.53        4597.97     10.15        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Methane(427) / '          190.53                    10.15        !Reducing parameters [K, mol/L]' /
DATA Methane(428) / '          10.15                                  !gamma' /
DATA Methane(429) / '          0.0831434                              !Gas constant [L-bar/mol-K]' /
DATA Methane(430) / '                32       1                       !Nterm, Ncoeff per term' /
DATA Methane(431) / '           0.9898937956e-4       0.2199608275         -5.322788' /
DATA Methane(432) / '           202.1657962          -22343.98926           0.000106794028' /
DATA Methane(433) / '           0.001457922469       -9.265816666           2915.364732' /
DATA Methane(434) / '           0.2313546209e-5       0.001387214274        0.04780467451' /
DATA Methane(435) / '           0.0001176103833      -0.00198209673        -0.2512887756' /
DATA Methane(436) / '           0.9748899826e-4      -0.1202192137e-5       0.0004128353939' /
DATA Methane(437) / '          -0.7215842918e-5       5081.738255          -919890.3192' /
DATA Methane(438) / '          -27.32264677           749902.4351           0.01114060908' /
DATA Methane(439) / '           10.83955159          -0.0004490960312      -13.80337847' /
DATA Methane(440) / '          -0.2371902232e-6       0.0003761652197      -0.2375166954e-8' /
DATA Methane(441) / '          -0.123764079e-6        0.6766926453e-5' /
DATA Methane(442) / '' /
DATA Methane(443) / '' /
DATA Methane(444) / '@AUX    !---Auxiliary function for Cp0' /
DATA Methane(445) / 'CP2     !Ideal gas heat capacity function for methane of Younglove and Ely.' /
DATA Methane(446) / '          ?' /
DATA Methane(447) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(448) / '          ?Younglove, B.A. and Ely, J.F.,' /
DATA Methane(449) / '          ?' /
DATA Methane(450) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(451) / '          0.                 !' /
DATA Methane(452) / '          10000.             !' /
DATA Methane(453) / '          0.                 !' /
DATA Methane(454) / '          0.                 !' /
DATA Methane(455) / '          1.0     8.31434    !Reducing parameters for T, Cp0' /
DATA Methane(456) / '          7 1   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Methane(457) / '          -1804475.0507     -3.0' /
DATA Methane(458) / '           77426.666393     -2.0' /
DATA Methane(459) / '          -1324.1658754     -1.0' /
DATA Methane(460) / '           15.438149595      0.0' /
DATA Methane(461) / '          -0.051479005257    1.0' /
DATA Methane(462) / '           0.00010809172196  2.0' /
DATA Methane(463) / '          -0.65501783437e-7  3.0' /
DATA Methane(464) / '          -6.7490056171      3000.0' /
DATA Methane(465) / '' /
DATA Methane(466) / '' /
DATA Methane(467) / '' /
DATA Methane(468) / '' /
DATA Methane(469) / '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' /
DATA Methane(470) / '' /
DATA Methane(471) / '#ETA   !---Viscosity---' /
DATA Methane(472) / 'VS4    !Pure fluid generalized friction theory viscosity model for methane of Quinones-Cisneros et al. (2011). unpublished' /
DATA Methane(473) / ':DOI:' /
DATA Methane(474) / '?' /
DATA Methane(475) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(476) / '?Quinones-Cisneros, S.E., Huber, M.L., and Deiters, U.K.,' /
DATA Methane(477) / '? unpublished work, 2011.' /
DATA Methane(478) / '?' /
DATA Methane(479) / '?Detailed uncertainty analysis will be found in a future publication; however' /
DATA Methane(480) / '? in general the estimated uncertainty in viscosity varies from less than 0.3%' /
DATA Methane(481) / '? between 200-400 K for pressures less than 30 MPa, to less than 2% over the' /
DATA Methane(482) / '? rest of the fluid surface up to 100 MPa, increasing up to 5%' /
DATA Methane(483) / '? for 100 to 500 MPa, and 10% at 500 to 1000 MPa for temperatures to 625 K.' /
DATA Methane(484) / '? Above uncertainties are valid when used with the equation of state of' /
DATA Methane(485) / '? Setzmann, U. and Wagner, W., J. Phys. Chem. Ref. Data, 20(6):1061-1151, 1991.' /
DATA Methane(486) / '? The use of other equations of state may result in larger uncertainties.' /
DATA Methane(487) / '?' /
DATA Methane(488) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(489) / '90.68              !Lower temperature limit [K]' /
DATA Methane(490) / '1200.0             !Upper temperature limit [K]' /
DATA Methane(491) / '100000.0           !Upper pressure limit [kPa]' /
DATA Methane(492) / '40.072             !Maximum density [mol/L]' /
DATA Methane(493) / '6 0 0 0 0 0        !Number of terms associated with dilute-gas function' /
DATA Methane(494) / 'NUL                !Pointer to reduced effective collision cross-section model; not used' /
DATA Methane(495) / '0.36652            !Lennard-Jones coefficient sigma [nm];not used' /
DATA Methane(496) / '174.0              !Lennard-Jones coefficient epsilon/kappa [K];not used' /
DATA Methane(497) / '190.564     1.0    !Reducing parameters for T, eta' /
DATA Methane(498) / '0.0         0.5    !Chapman-Enskog term; not used here' /
DATA Methane(499) / ' 58.343920516258155    0.0' /
DATA Methane(500) / '-199.92388279110893    0.25' /
DATA Methane(501) / ' 240.35409195445984    0.5' /
DATA Methane(502) / '-113.08166560748158    0.75' /
DATA Methane(503) / ' 21.645948012444557    1.0' /
DATA Methane(504) / '0                  !Number of terms for initial density dependence' /
DATA Methane(505) / '-0.00002946520026265898 0.000011850361299482738       0.0                       0.0       0.0   !  a(0),a(1),a(2)' /
DATA Methane(506) / ' 0.00002700022529490106 -0.000032677520832951284      0.0                       0.0       0.0   !  b(0),b(1),b(2)' /
DATA Methane(507) / ' 0.00002904479739920783 -0.00001018049342159992      -3.095500930526404e-8      0.0       0.0   !  c(0),c(1),c(2)' /
DATA Methane(508) / ' 1.55372118714633e-8    -1.944037783173382e-9         0.0                       0.0       0.0   !  A(0),A(1),A(2)' /
DATA Methane(509) / '-2.6710447337075816e-9  3.2621373142076857e-9         0.0                       0.0       0.0   !  B(0),B(1),B(2)' /
DATA Methane(510) / ' 5.207541202169661e-9   1.5949945307134116e-7         3.687831977089463e-10     0.0       0.0   !  C(0),C(1),C(2)' /
DATA Methane(511) / ' 3.0218122078964884e-12 0.0                           0.0                       0.0       0.0   !  D(0),D(1),D(2)' /
DATA Methane(512) / ' 0.0                    0.0                           0.0                       0.0       0.0   !  E(0),E(1),E(2)' /
DATA Methane(513) / 'NUL                !Pointer to the viscosity critical enhancement auxiliary function (none used)' /
DATA Methane(514) / '' /
DATA Methane(515) / '' /
DATA Methane(516) / '' /
DATA Methane(517) / '' /
DATA Methane(518) / '================================================================================' /
DATA Methane(519) / '' /
DATA Methane(520) / '#TCX   !---Thermal conductivity---' /
DATA Methane(521) / 'TC1    !Pure fluid thermal conductivity model for methane of Friend et al. (1989).' /
DATA Methane(522) / ':DOI:' /
DATA Methane(523) / ':WEB: https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1325.pdf' /
DATA Methane(524) / '?' /
DATA Methane(525) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(526) / '?Friend, D.G., Ely, J.F., and Ingham, H.,' /
DATA Methane(527) / '? "Tables for the Thermophysical Properties of Methane,"' /
DATA Methane(528) / '? NIST Technical Note 1325, 1989.' /
DATA Methane(529) / '?' /
DATA Methane(530) / '?The uncertainty in thermal conductivity of the dilute gas between 130' /
DATA Methane(531) / '? and 625 K is 2.5%.  For temperatures below 130 K, the uncertainty is' /
DATA Methane(532) / '? less than 10%.  Excluding the dilute gas, the uncertainty is 2% between' /
DATA Methane(533) / '? 110 and 725 K at pressures up to 70 MPa, except near the critical' /
DATA Methane(534) / '? point which has an uncertainty of 5% or greater.  For the vapor at lower' /
DATA Methane(535) / '? temperatures and the dense liquid near the triple point, an uncertainty of' /
DATA Methane(536) / '? 10% is possible.' /
DATA Methane(537) / '?' /
DATA Methane(538) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(539) / '90.6941            !Lower temperature limit [K]' /
DATA Methane(540) / '625.0              !Upper temperature limit [K]' /
DATA Methane(541) / '1000000.0          !Upper pressure limit [kPa]' /
DATA Methane(542) / '40.072             !Maximum density [mol/L]' /
DATA Methane(543) / '3   0              !# terms for dilute gas function:  numerator, denominator' /
DATA Methane(544) / ' 174.0       0.001 !Reducing parameters for T, tcx' /
DATA Methane(545) / ' 1.45885     0.    !Coefficient, power in T' /
DATA Methane(546) / '-0.4377162  -1.' /
DATA Methane(547) / ' 0.         -96.   !Coefficient, power in T' /
DATA Methane(548) / '8   0              !# terms for background gas function:  numerator, denominator' /
DATA Methane(549) / ' 190.551  10.139    0.00629638         !Reducing parameters for T, rho, tcx' /
DATA Methane(550) / ' 1.5554612   0.  2.  0.  !Coefficient, powers of T, rho, exp(rho)' /
DATA Methane(551) / ' 1.0         0.  0. -99. !The order here is important' /
DATA Methane(552) / ' 2.4149207   0.  1.  0.' /
DATA Methane(553) / ' 0.55166331  0.  3.  0.' /
DATA Methane(554) / '-0.52837734  0.  4.  0.' /
DATA Methane(555) / ' 0.073809553 -1. 4.  0.' /
DATA Methane(556) / ' 0.24465507  0.  5.  0.' /
DATA Methane(557) / '-0.047613626 -1. 5.  0.' /
DATA Methane(558) / 'TK7                !Pointer to critical enhancement auxiliary function' /
DATA Methane(559) / '' /
DATA Methane(560) / '' /
DATA Methane(561) / '#AUX   !---Auxiliary function for the thermal conductivity critical enhancement' /
DATA Methane(562) / 'TK7    !Thermal conductivity critical enhancement for methane of Friend et al. (1989).' /
DATA Methane(563) / '?' /
DATA Methane(564) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(565) / '?Friend, D.G., Ely, J.F., and Ingham, H., 1989.' /
DATA Methane(566) / '?' /
DATA Methane(567) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(568) / '0.                 !' /
DATA Methane(569) / '10000.             !' /
DATA Methane(570) / '0.                 !' /
DATA Methane(571) / '0.                 !' /
DATA Methane(572) / '$CE RED 1 TR - =TAU  1 DR - =DEL' /
DATA Methane(573) / '$CE 1 DPDT RGAS / DENS / CNST DENS POP< =V1' /
DATA Methane(574) / '$CE CNST CNST DENS * DRED / TRED * DPDD / RGAS * DUP 0 POP> =V2' /
DATA Methane(575) / '$CE CNST ETA / TRED TEMP / SQR / V1 SQR * V2 CNST POWR *' /
DATA Methane(576) / '$CE TAU ABS SQRT CNST * CNST DEL SQR * + CNST DEL * - SIGN EXP *' /
DATA Methane(577) / '$CF' /
DATA Methane(578) / '0.001       190.564 10.139128 0. 0' /
DATA Methane(579) / '1.E-12      0. 0. 0. 0' /
DATA Methane(580) / '1.e5        0. 0. 0. 0' /
DATA Methane(581) / '0.28631     0. 0. 0. 0' /
DATA Methane(582) / '91.855      0. 0. 0. 0' /
DATA Methane(583) / '0.4681      0. 0. 0. 0' /
DATA Methane(584) / '2.646       0. 0. 0. 0' /
DATA Methane(585) / '2.678       0. 0. 0. 0' /
DATA Methane(586) / '0.637       0. 0. 0. 0' /
DATA Methane(587) / '' /
DATA Methane(588) / '' /
DATA Methane(589) / '' /
DATA Methane(590) / '' /
DATA Methane(591) / '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' /
DATA Methane(592) / '' /
DATA Methane(593) / '@TRN    !---ECS Transport---' /
DATA Methane(594) / 'ECS     !Extended Corresponding States model (Nitrogen reference); predictive mode for methane.' /
DATA Methane(595) / '          ?' /
DATA Methane(596) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(597) / '          ?Klein, S.A., McLinden, M.O., and Laesecke, A., "An Improved Extended Corresponding States Method for Estimation of Viscosity of Pure Refrigerants and Mixtures," Int. J. Refrigeration, 20(3):208-217, 1997. doi: 10.1016/S0140-7007(96)00073-4.' /
DATA Methane(598) / '          ?McLinden, M.O., Klein, S.A., and Perkins, R.A., "An Extended Corresponding States Model for the Thermal Conductivity of Refrigerants and Refrigerant Mixtures," Int. J. Refrigeration, 23(1):43-63, 2000. doi: 10.1016/S0140-7007(99)00024-9' /
DATA Methane(599) / '          ?' /
DATA Methane(600) / '          ?The Lennard-Jones parameters were taken from Friend, D.G., Ely, J.F., and Ingham, H., "Tables for the Thermophysical Properties of Methane," NIST Technical Note 1325, 1989.' /
DATA Methane(601) / '          ?' /
DATA Methane(602) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(603) / '          90.6941            !Lower temperature limit [K]' /
DATA Methane(604) / '          625.0              !Upper temperature limit [K]' /
DATA Methane(605) / '          1000000.0          !Upper pressure limit [kPa]' /
DATA Methane(606) / '          40.072             !Maximum density [mol/L]' /
DATA Methane(607) / '          FEQ NITROGEN.FLD' /
DATA Methane(608) / '          VS1                !Model for reference fluid viscosity' /
DATA Methane(609) / '          TC1                !Model for reference fluid thermal conductivity' /
DATA Methane(610) / '          NUL                !Large molecule identifier' /
DATA Methane(611) / '          1                  !Lennard-Jones flag (0 or 1) (0 => use estimates)' /
DATA Methane(612) / '          0.36652            !Lennard-Jones coefficient sigma [nm]' /
DATA Methane(613) / '          174.0              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA Methane(614) / '          1  0  0            !Number of terms in f_int term in Eucken correlation, spare1, spare2' /
DATA Methane(615) / '           0.00132  0. 0. 0. !Coefficient, power of T, spare1, spare2' /
DATA Methane(616) / '          1  0  0            !Number of terms in psi (visc shape factor): poly,spare1,spare2' /
DATA Methane(617) / '           1.0      0. 0. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Methane(618) / '          1  0  0            !Number of terms in chi (t.c. shape factor): poly,spare1,spare2' /
DATA Methane(619) / '           1.0      0. 0. 0. !Coefficient, power of Tr, power of Dr, spare' /
DATA Methane(620) / '          TK3                !Pointer to critical enhancement auxiliary function' /
DATA Methane(621) / '' /
DATA Methane(622) / '' /
DATA Methane(623) / '' /
DATA Methane(624) / '' /
DATA Methane(625) / '********************************************************************************' /
DATA Methane(626) / '' /
DATA Methane(627) / '@TCX    !---Thermal conductivity---' /
DATA Methane(628) / 'TC2     !Pure fluid thermal conductivity model for methane of Younglove and Ely (1987).' /
DATA Methane(629) / '          ?' /
DATA Methane(630) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(631) / '          ?Younglove, B.A. and Ely, J.F.,' /
DATA Methane(632) / '          ? "Thermophysical properties of fluids. II. Methane, ethane, propane,' /
DATA Methane(633) / '          ? isobutane and normal butane,"' /
DATA Methane(634) / '          ? J. Phys. Chem. Ref. Data, 16:577-798, 1987.' /
DATA Methane(635) / '          ?' /
DATA Methane(636) / '          ?The uncertainty in thermal conductivity is 5% in the liquid, 4% in the vapor,' /
DATA Methane(637) / '          ? 3% at T>Tc, and 8% in the critical region.' /
DATA Methane(638) / '          ?' /
DATA Methane(639) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(640) / '          90.6941            !Lower temperature limit [K]' /
DATA Methane(641) / '          625.0              !Upper temperature limit [K]' /
DATA Methane(642) / '          1000000.0          !Upper pressure limit [kPa]' /
DATA Methane(643) / '          40.072             !Maximum density [mol/L]' /
DATA Methane(644) / '          CI2                !Pointer to collision integral model' /
DATA Methane(645) / '          0.368              !Lennard-Jones coefficient sigma [nm]' /
DATA Methane(646) / '          168.0              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA Methane(647) / '          0.1069188          !Const in Eq 19 = 5/16*(k*MW/1000/pi/Na)**0.5*1.0d12' /
DATA Methane(648) / '           1.346953698       !Dilute gas terms (Eq 27):  Gt(1)' /
DATA Methane(649) / '          -0.3254677753      !                           Gt(2)' /
DATA Methane(650) / '           0.002325800819    !Residual terms (Eqs 26, 28-30): Et(1)' /
DATA Methane(651) / '          -0.2477927999' /
DATA Methane(652) / '           38.80593713' /
DATA Methane(653) / '          -0.1579519146e-6' /
DATA Methane(654) / '           0.003717991328' /
DATA Methane(655) / '          -0.9616989434' /
DATA Methane(656) / '          -0.03017352774' /
DATA Methane(657) / '           0.4298153386      !Et(8)' /
DATA Methane(658) / '          TK2                !Pointer to critical enhancement model (follows immediately)' /
DATA Methane(659) / '           37.42368          !Critical enhancement terms (Eqs D1-D4):  X1' /
DATA Methane(660) / '           3.16714' /
DATA Methane(661) / '           0.78035' /
DATA Methane(662) / '           0.60103           !X4' /
DATA Methane(663) / '           6.512707e-10      !Z' /
DATA Methane(664) / '           1.38054e-23       !Boltzmanns constant, k' /
DATA Methane(665) / '           0.16969859271     !Coefficient for initial density dependence of viscosity (eq 21); Fv(1)' /
DATA Methane(666) / '          -0.013337234608    !Fv(2)' /
DATA Methane(667) / '           1.4               !Fv(3)' /
DATA Methane(668) / '           168.              !Fv(4)' /
DATA Methane(669) / '          -16.20427429       !Coefficients for residual viscosity, eqs (22 - 25)' /
DATA Methane(670) / '           427.0589027       !Ev(2)  (the viscosity is also used in conductivity correlation)' /
DATA Methane(671) / '           14.02596278       !Ev(3)' /
DATA Methane(672) / '          -3916.837745       !Ev(4)' /
DATA Methane(673) / '          -0.0347709909      !Ev(5)' /
DATA Methane(674) / '           21.36542674       !Ev(6)' /
DATA Methane(675) / '           1436.802482       !Ev(7)' /
DATA Methane(676) / '' /
DATA Methane(677) / '' /
DATA Methane(678) / '@AUX    !---Auxiliary function for the thermal conductivity critical enhancement' /
DATA Methane(679) / 'TK3     !Simplified thermal conductivity critical enhancement for methane of Olchowy and Sengers (1989).' /
DATA Methane(680) / '          ?' /
DATA Methane(681) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(682) / '          ?Olchowy, G.A. and Sengers, J.V.,' /
DATA Methane(683) / '          ? "A simplified representation for the thermal conductivity of fluids in the critical region,"' /
DATA Methane(684) / '          ? Int. J. Thermophysics, 10:417-426, 1989. doi: 10.1007/BF01133538' /
DATA Methane(685) / '          ?' /
DATA Methane(686) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(687) / '          0.                 !' /
DATA Methane(688) / '          10000.             !' /
DATA Methane(689) / '          0.                 !' /
DATA Methane(690) / '          0.                 !' /
DATA Methane(691) / '          9 0 0 0            !# terms:  CO2-terms, spare, spare, spare' /
DATA Methane(692) / '          1.0  1.0  1.0      !Reducing parameters for T, rho, tcx [mW/(m-K)]' /
DATA Methane(693) / '          0.63               !Nu (universal exponent)' /
DATA Methane(694) / '          1.239              !Gamma (universal exponent)' /
DATA Methane(695) / '          1.03               !R0 (universal amplitude)' /
DATA Methane(696) / '          0.063              !Z (universal exponent--not used for t.c., only viscosity)' /
DATA Methane(697) / '          1.0                !C (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA Methane(698) / '          0.194e-9           !Xi0 (amplitude) [m]' /
DATA Methane(699) / '          0.0496             !Gam0 (amplitude) [-]' /
DATA Methane(700) / '          0.4e-9             !Qd_inverse (modified effective cutoff parameter) [m]; estmated value from matching Friend at 50 bar' /
DATA Methane(701) / '          285.846            !Tref (reference temperature)=1.5*Tc [K]' /
DATA Methane(702) / '' /
DATA Methane(703) / '' /
DATA Methane(704) / '@ETA    !---Viscosity---' /
DATA Methane(705) / 'VS1     !Pure fluid viscosity model for methane of Vogel et al. (2000).' /
DATA Methane(706) / '          ?' /
DATA Methane(707) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(708) / '          ?Vogel, E., Wilhelm, J., Kuechenmeister, C., and Jaesche, M.,' /
DATA Methane(709) / '          ? "High-precision viscosity measurements on methane,"' /
DATA Methane(710) / '          ? High Temp. - High Pressures, 32(1):73-81, 2000.' /
DATA Methane(711) / '          ?' /
DATA Methane(712) / '          ?The uncertainty in viscosity varies from 0.3% in the dilute gas between' /
DATA Methane(713) / '          ? 260-360 K, to 3.0% over the rest of the fluid surface, increasing up to 5 %' /
DATA Methane(714) / '          ? from 620 K and 100 MPa.' /
DATA Methane(715) / '          ?' /
DATA Methane(716) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(717) / '          90.6941            !Lower temperature limit [K]' /
DATA Methane(718) / '          625.0              !Upper temperature limit [K]' /
DATA Methane(719) / '          100000.0           !Upper pressure limit [kPa]' /
DATA Methane(720) / '          40.072             !Maximum density [mol/L]' /
DATA Methane(721) / '          1                  !Number of terms associated with dilute-gas function' /
DATA Methane(722) / '          CI1                !Pointer to reduced effective collision cross-section model' /
DATA Methane(723) / '          0.37333            !Lennard-Jones coefficient sigma [nm]' /
DATA Methane(724) / '          160.78             !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA Methane(725) / '          1.0      1.0       !Reducing parameters for T, eta' /
DATA Methane(726) / '          0.0855422   0.5    !Chapman-Enskog term sqrt[MW]*0.021357' /
DATA Methane(727) / '          9                  !Number of terms for initial density dependence' /
DATA Methane(728) / '          159.7   0.0306525  !Reducing parameters for T (= eps/k), etaB2 (= 0.6022137*sigma**3)' /
DATA Methane(729) / '          -19.572881  0.0    !Coefficient, power in T* = T/(eps/k)' /
DATA Methane(730) / '           219.73999   -0.25' /
DATA Methane(731) / '          -1015.3226   -0.5' /
DATA Methane(732) / '           2471.01251  -0.75' /
DATA Methane(733) / '          -3375.1717   -1.0' /
DATA Methane(734) / '           2491.6597   -1.25' /
DATA Methane(735) / '          -787.26086   -1.5' /
DATA Methane(736) / '           14.085455   -2.5' /
DATA Methane(737) / '          -0.34664158  -5.5' /
DATA Methane(738) / '          1 9 1 2 0 0        !# resid terms:  close-packed density;  simple poly; numerator of rational poly; denominator of rat. poly; numerator of exponential; denominator of exponential' /
DATA Methane(739) / '          190.564   10.139        1.0            !Reducing parameters for T, rho, eta' /
DATA Methane(740) / '           3.10860501398   0.  0.  0.  0' /
DATA Methane(741) / '          -3.02256904347   0.  2.  0.  0' /
DATA Methane(742) / '           17.6965130175  -1.  2.  0.  0' /
DATA Methane(743) / '           3.11150846518   0.  3.  0.  0' /
DATA Methane(744) / '          -21.5685107769  -1.  3.  0.  0' /
DATA Methane(745) / '           0.672852409238  0.  4.  0.  0' /
DATA Methane(746) / '           10.2387524315  -1.  4.  0.  0' /
DATA Methane(747) / '          -1.09330775541   0.  5.  0.  0' /
DATA Methane(748) / '          -1.20030749419  -1.  5.  0.  0' /
DATA Methane(749) / '          -21.1009923406   0.  1. -1.  0' /
DATA Methane(750) / '           21.1009923406   0.  1.  0.  0' /
DATA Methane(751) / '           1.0             0.  0.  1.  0' /
DATA Methane(752) / '          -1.0             0.  1.  0.  0' /
DATA Methane(753) / '          NUL                !Pointer to the viscosity critical enhancement auxiliary function (none used)' /
DATA Methane(754) / '' /
DATA Methane(755) / '' /
DATA Methane(756) / '@AUX    !---Auxiliary function for the collision integral' /
DATA Methane(757) / 'CI1     !Collision integral model for methane of Vogel et al. (2000).' /
DATA Methane(758) / '          ?' /
DATA Methane(759) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(760) / '          ?Vogel, E., Wilhelm, J., Kuechenmeister, C., and Jaesche, M.,' /
DATA Methane(761) / '          ?' /
DATA Methane(762) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(763) / '          0.                 !' /
DATA Methane(764) / '          10000.             !' /
DATA Methane(765) / '          0.                 !' /
DATA Methane(766) / '          0.                 !' /
DATA Methane(767) / '          5                  !Number of terms' /
DATA Methane(768) / '           0.215309028    0  !Coefficient, power of Tstar' /
DATA Methane(769) / '          -0.46256942     1' /
DATA Methane(770) / '           0.051313823    2' /
DATA Methane(771) / '           0.030320660    3' /
DATA Methane(772) / '          -0.0070047029   4' /
DATA Methane(773) / '' /
DATA Methane(774) / '' /
DATA Methane(775) / '@ETA    !---Viscosity---' /
DATA Methane(776) / 'VS2     !Pure fluid viscosity model for methane of Younglove and Ely (1987).' /
DATA Methane(777) / '          ?' /
DATA Methane(778) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(779) / '          ?Younglove, B.A. and Ely, J.F.,' /
DATA Methane(780) / '          ? "Thermophysical properties of fluids. II. Methane, ethane, propane,' /
DATA Methane(781) / '          ? isobutane and normal butane,"' /
DATA Methane(782) / '          ? J. Phys. Chem. Ref. Data, 16:577-798, 1987.' /
DATA Methane(783) / '          ? All temperatures on IPTS-68' /
DATA Methane(784) / '          ?' /
DATA Methane(785) / '          ?The uncertainty in viscosity is 2%, except in the critical region which is 5%.' /
DATA Methane(786) / '          ?' /
DATA Methane(787) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(788) / '          90.6941            !Lower temperature limit [K]' /
DATA Methane(789) / '          625.0              !Upper temperature limit [K]' /
DATA Methane(790) / '          1000000.0          !Upper pressure limit [kPa]' /
DATA Methane(791) / '          40.072             !Maximum density [mol/L]' /
DATA Methane(792) / '          CI2                !Pointer to collision integral model' /
DATA Methane(793) / '          0.368              !Lennard-Jones coefficient sigma [nm]' /
DATA Methane(794) / '          168.0              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA Methane(795) / '           0.1069188         !Const in Eq 19 = 5/16*(k*MW/1000/pi/Na)**0.5*1.0d12' /
DATA Methane(796) / '           0.5               !Exponent in Eq 19 for T' /
DATA Methane(797) / '           0.16969859271     !Coefficient for initial density dependence of viscosity (eq 21); Fv(1)' /
DATA Methane(798) / '          -0.013337234608    !Fv(2)' /
DATA Methane(799) / '           1.4               !Fv(3)' /
DATA Methane(800) / '           168.0             !Fv(4)' /
DATA Methane(801) / '          -16.20427429       !Coefficients for residual viscosity, eqs (22 - 25)' /
DATA Methane(802) / '           427.0589027       !Ev(2)' /
DATA Methane(803) / '           14.02596278       !Ev(3)' /
DATA Methane(804) / '          -3916.837745       !Ev(4)' /
DATA Methane(805) / '          -0.0347709909      !Ev(5)' /
DATA Methane(806) / '           21.36542674       !Ev(6)' /
DATA Methane(807) / '           1436.802482       !Ev(7)' /
DATA Methane(808) / '           10.15             !Ev(8)' /
DATA Methane(809) / '          NUL                !Pointer to the viscosity critical enhancement auxiliary function (none used)' /
DATA Methane(810) / '' /
DATA Methane(811) / '' /
DATA Methane(812) / '@AUX    !---Auxiliary function for the collision integral' /
DATA Methane(813) / 'CI2     !Collision integral model for methane of Younglove and Ely (1987).' /
DATA Methane(814) / '          ?' /
DATA Methane(815) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(816) / '          ?Friend, D.G., Ely, J.F., and Ingham, H.,' /
DATA Methane(817) / '          ? "Tables for the Thermophysical Properties of Methane,"' /
DATA Methane(818) / '          ? NIST Technical Note 1325, 1989.' /
DATA Methane(819) / '          ?' /
DATA Methane(820) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(821) / '          0.                 !' /
DATA Methane(822) / '          10000.             !' /
DATA Methane(823) / '          0.                 !' /
DATA Methane(824) / '          0.                 !' /
DATA Methane(825) / '          9                  !Number of terms' /
DATA Methane(826) / '          -3.0328138281   0' /
DATA Methane(827) / '           16.918880086   0' /
DATA Methane(828) / '          -37.189364917   0' /
DATA Methane(829) / '           41.288861858   0' /
DATA Methane(830) / '          -24.615921140   0' /
DATA Methane(831) / '           8.9488430959   0' /
DATA Methane(832) / '          -1.8739245042   0' /
DATA Methane(833) / '           0.20966101390  0' /
DATA Methane(834) / '          -0.009657043707 0' /
DATA Methane(835) / '' /
DATA Methane(836) / '' /
DATA Methane(837) / '' /
DATA Methane(838) / '' /
DATA Methane(839) / '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' /
DATA Methane(840) / '' /
DATA Methane(841) / '#STN   !---Surface tension---' /
DATA Methane(842) / 'ST1    !Surface tension model for methane of Mulero et al. (2012).' /
DATA Methane(843) / ':DOI: 10.1063/1.4768782' /
DATA Methane(844) / '?' /
DATA Methane(845) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(846) / '?Mulero, A., Cachadina, I., and Parra, M.I.,' /
DATA Methane(847) / '? "Recommended Correlations for the Surface Tension of Common Fluids,"' /
DATA Methane(848) / '? J. Phys. Chem. Ref. Data, 41(4), 043105, 2012. doi: 10.1063/1.4768782' /
DATA Methane(849) / '?' /
DATA Methane(850) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(851) / '0.                 !' /
DATA Methane(852) / '10000.             !' /
DATA Methane(853) / '0.                 !' /
DATA Methane(854) / '0.                 !' /
DATA Methane(855) / '3                  !Number of terms in surface tension model' /
DATA Methane(856) / '190.564            !Critical temperature used in fit (dummy)' /
DATA Methane(857) / ' 0.03825   1.191   !Sigma0 and n' /
DATA Methane(858) / '-0.006024  5.422' /
DATA Methane(859) / '-0.0007065 0.6161' /
DATA Methane(860) / '' /
DATA Methane(861) / '' /
DATA Methane(862) / '#DE    !---Dielectric constant---' /
DATA Methane(863) / 'DE3    !Dielectric constant model for methane of Harvey and Lemmon (2005).' /
DATA Methane(864) / ':DOI: 10.1007/s10765-005-2351-5' /
DATA Methane(865) / '?' /
DATA Methane(866) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(867) / '?Harvey, A.H. and Lemmon, E.W.,' /
DATA Methane(868) / '? "Method for Estimating the Dielectric Constant of Natural Gas Mixtures,"' /
DATA Methane(869) / '? Int. J. Thermophys., 26(1):31-46, 2005. doi: 10.1007/s10765-005-2351-5' /
DATA Methane(870) / '?' /
DATA Methane(871) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(872) / '0.                 !' /
DATA Methane(873) / '10000.             !' /
DATA Methane(874) / '0.                 !' /
DATA Methane(875) / '0.                 !' /
DATA Methane(876) / '273.16 1000.0 1.0  !Reducing parameters for T and D' /
DATA Methane(877) / '0 2 4 0 0 0        !Number of terms in dielectric constant model' /
DATA Methane(878) / ' 6.5443   0. 1. 0. !Coefficient, T exp, D exp' /
DATA Methane(879) / ' 0.0133   1. 1. 0.' /
DATA Methane(880) / ' 8.4578   0. 2. 0.' /
DATA Methane(881) / ' 3.7196   1. 2. 0.' /
DATA Methane(882) / '-352.97   0. 3. 0.' /
DATA Methane(883) / '-100.65   1. 3. 0.' /
DATA Methane(884) / '' /
DATA Methane(885) / '' /
DATA Methane(886) / '#MLT   !---Melting line---' /
DATA Methane(887) / 'ML1    !Melting line model for methane of Setzmann and Wagner (1991).' /
DATA Methane(888) / ':DOI: 10.1063/1.555898' /
DATA Methane(889) / '?' /
DATA Methane(890) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(891) / '?Setzmann, U. and Wagner, W.,' /
DATA Methane(892) / '? "A New Equation of State and Tables of Thermodynamic Properties for Methane' /
DATA Methane(893) / '? Covering the Range from the Melting Line to 625 K at Pressures up to 1000 MPa,"' /
DATA Methane(894) / '? J. Phys. Chem. Ref. Data, 20(6):1061-1151, 1991.' /
DATA Methane(895) / '?' /
DATA Methane(896) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(897) / '90.6941            !Lower temperature limit [K]' /
DATA Methane(898) / '625.0              !Upper temperature limit [K]' /
DATA Methane(899) / '0.                 !' /
DATA Methane(900) / '0.                 !' /
DATA Methane(901) / '90.6941  11.696    !Reducing temperature and pressure' /
DATA Methane(902) / '5 0 0 0 0 0        !Number of terms in melting line equation' /
DATA Methane(903) / ' 1.0       0.0     !Coefficients and exponents' /
DATA Methane(904) / ' 24756.8   1.85' /
DATA Methane(905) / '-7366.02   2.1' /
DATA Methane(906) / '-24756.8   0.0' /
DATA Methane(907) / ' 7366.02   0.0' /
DATA Methane(908) / '' /
DATA Methane(909) / '' /
DATA Methane(910) / '#SBL   !---Sublimation line---' /
DATA Methane(911) / 'SB3    !Sublimation line model for methane of Lemmon (2002).' /
DATA Methane(912) / ':DOI:' /
DATA Methane(913) / '?' /
DATA Methane(914) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(915) / '?Lemmon, E.W., 2002.' /
DATA Methane(916) / '?' /
DATA Methane(917) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(918) / '0.                 !' /
DATA Methane(919) / '90.6941            !Upper temperature limit [K]' /
DATA Methane(920) / '0.                 !' /
DATA Methane(921) / '0.                 !' /
DATA Methane(922) / '90.6941  11.696    !Reducing temperature and pressure' /
DATA Methane(923) / '0 1 0 0 0 0        !Number of terms in sublimation line equation' /
DATA Methane(924) / '-12.84  1.         !Coefficients and exponents' /
DATA Methane(925) / '' /
DATA Methane(926) / '' /
DATA Methane(927) / '#PS    !---Vapor pressure---' /
DATA Methane(928) / 'PS5    !Vapor pressure equation for methane of Setzmann and Wagner (1991).' /
DATA Methane(929) / '?' /
DATA Methane(930) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(931) / '?Functional Form:  P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T] where Theta=1-T/Tc, Tc and Pc' /
DATA Methane(932) / '? are the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Methane(933) / '?' /
DATA Methane(934) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(935) / '0.                 !' /
DATA Methane(936) / '10000.             !' /
DATA Methane(937) / '0.                 !' /
DATA Methane(938) / '0.                 !' /
DATA Methane(939) / '190.564 4599.2     !Reducing parameters' /
DATA Methane(940) / '4 0 0 0 0 0        !Number of terms in equation' /
DATA Methane(941) / '-6.036219  1.0' /
DATA Methane(942) / ' 1.409353  1.5' /
DATA Methane(943) / '-0.4945199 2.0' /
DATA Methane(944) / '-1.443048  4.5' /
DATA Methane(945) / '' /
DATA Methane(946) / '' /
DATA Methane(947) / '#DL    !---Saturated liquid density---' /
DATA Methane(948) / 'DL4    !Saturated liquid density equation for methane of Setzmann and Wagner (1991).' /
DATA Methane(949) / '?' /
DATA Methane(950) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(951) / '?Functional Form:  D=Dc*EXP[SUM(Ni*Theta^(ti/3))] where Theta=1-T/Tc, Tc and Dc are' /
DATA Methane(952) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Methane(953) / '?' /
DATA Methane(954) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(955) / '0.                 !' /
DATA Methane(956) / '10000.             !' /
DATA Methane(957) / '0.                 !' /
DATA Methane(958) / '0.                 !' /
DATA Methane(959) / '190.564  10.139128 !Reducing parameters' /
DATA Methane(960) / '3 0 0 0 0 0        !Number of terms in equation' /
DATA Methane(961) / ' 1.9906389   1.062' /
DATA Methane(962) / '-0.78756197  1.5' /
DATA Methane(963) / ' 0.036976723 7.5' /
DATA Methane(964) / '' /
DATA Methane(965) / '' /
DATA Methane(966) / '#DV    !---Saturated vapor density---' /
DATA Methane(967) / 'DV4    !Saturated vapor density equation for methane of Setzmann and Wagner (1991).' /
DATA Methane(968) / '?' /
DATA Methane(969) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(970) / '?Functional Form:  D=Dc*EXP[SUM(Ni*Theta^(ti/3))] where Theta=1-T/Tc, Tc and Dc are' /
DATA Methane(971) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Methane(972) / '?' /
DATA Methane(973) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(974) / '0.                 !' /
DATA Methane(975) / '10000.             !' /
DATA Methane(976) / '0.                 !' /
DATA Methane(977) / '0.                 !' /
DATA Methane(978) / '190.564  10.139128 !Reducing parameters' /
DATA Methane(979) / '6 0 0 0 0 0        !Number of terms in equation' /
DATA Methane(980) / '-1.8802840 1.062' /
DATA Methane(981) / '-2.8526531 2.5' /
DATA Methane(982) / '-3.0006480 4.5' /
DATA Methane(983) / '-5.2511690 7.5' /
DATA Methane(984) / '-13.191859 12.5' /
DATA Methane(985) / '-37.553961 23.5' /
DATA Methane(986) / '' /
DATA Methane(987) / '' /
DATA Methane(988) / '@END' /
DATA Methane(989) / 'c        1         2         3         4         5         6         7         8' /
DATA Methane(990) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
DATA Methane(991) / '' /
DATA Methane(992) / '                    !Cant have two vs1 models in at the same time so put Friend here' /
DATA Methane(993) / '                    !It is limited to low p' /
DATA Methane(994) / '' /
DATA Methane(995) / '                    @ETA               !Viscosity model specification' /
DATA Methane(996) / '                    VS1  pure fluid viscosity model of Friend et al. (1989).' /
DATA Methane(997) / '                    ?' /
DATA Methane(998) / '                    ?```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(999) / '                    ?Friend, D.G., Ely, J.F., and Ingham, H.,' /
DATA Methane(1000) / '                    ? "Tables for the Thermophysical Properties of Methane,"' /
DATA Methane(1001) / '                    ? NIST Technical Note 1325, 1989.' /
DATA Methane(1002) / '                    ?' /
DATA Methane(1003) / '                    ?The uncertainty in viscosity is 0.5% between 270 and 600 K, and 1% above' /
DATA Methane(1004) / '                    ? 600 K.  Below 270 K, the uncertainty is 2%.' /
DATA Methane(1005) / '                    ?' /
DATA Methane(1006) / '                    !```````````````````````````````````````````````````````````````````````````````' /
DATA Methane(1007) / '                    90.6941            !Lower temperature limit [K]' /
DATA Methane(1008) / '                    625.0              !Upper temperature limit [K]' /
DATA Methane(1009) / '                    1000000.0          !Upper pressure limit [kPa]' /
DATA Methane(1010) / '                    40.072             !Maximum density [mol/L]' /
DATA Methane(1011) / '                    1                  !Number of terms associated with dilute-gas function' /
DATA Methane(1012) / '                    CI2                !Pointer to reduced effective collision cross-section model' /
DATA Methane(1013) / '                    0.36652            !Lennard-Jones coefficient sigma [nm]' /
DATA Methane(1014) / '                    174.0              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA Methane(1015) / '                    174.     10.0      !Reducing parameters for T, eta' /
DATA Methane(1016) / '                    0.14105376  0.5    !Chapman-Enskog term' /
DATA Methane(1017) / '                    0                  !Number of terms for initial density dependence' /
DATA Methane(1018) / '                    0 0 9 3 0 0        !# resid terms:  close-packed density;  simple poly; numerator of rational poly; denominator of rat. poly; numerator of exponential; denominator of exponential' /
DATA Methane(1019) / '                    190.551   10.139        12.149         !Reducing parameters for T, rho, eta' /
DATA Methane(1020) / '                     0.41250137      0.  1.  0.  0' /
DATA Methane(1021) / '                    -0.14390912     -1.  1.  0.  0' /
DATA Methane(1022) / '                     0.10366993      0.  2.  0.  0' /
DATA Methane(1023) / '                     0.40287464     -1.  2.  0.  0' /
DATA Methane(1024) / '                    -0.24903524     -1.5 2.  0.  0' /
DATA Methane(1025) / '                    -0.12953131      0.  3.  0.  0' /
DATA Methane(1026) / '                     0.06575776     -2.  3.  0.  0' /
DATA Methane(1027) / '                     0.02566628      0.  4.  0.  0' /
DATA Methane(1028) / '                    -0.03716526     -1.  4.  0.  0' /
DATA Methane(1029) / '                     1.0             0.  0.  0.  0' /
DATA Methane(1030) / '                    -0.38798341      0.  1.  0.  0' /
DATA Methane(1031) / '                     0.03533815     -1.  1.  0.  0' /
DATA Methane(1032) / '                    NUL                !Pointer to the viscosity critical enhancement auxiliary function (none used)' /
! #######################################################
character(256), TARGET :: nitrogen(733)
DATA nitrogen(1) / 'nitrogen           !short name' /
DATA nitrogen(2) / '7727-37-9          !CAS number' /
DATA nitrogen(3) / 'nitrogen           !full name' /
DATA nitrogen(4) / 'N2                 !chemical formula' /
DATA nitrogen(5) / 'R-728              !synonym' /
DATA nitrogen(6) / '28.01348           !molecular weight [g/mol]' /
DATA nitrogen(7) / '63.151             !triple point temperature [K]' /
DATA nitrogen(8) / '77.355             !normal boiling point [K]' /
DATA nitrogen(9) / '126.192            !critical temperature [K]' /
DATA nitrogen(10) / '3395.8             !critical pressure [kPa]' /
DATA nitrogen(11) / '11.1839            !critical density [mol/L]' /
DATA nitrogen(12) / '0.0372             !acentric factor' /
DATA nitrogen(13) / '0.0                !dipole moment [Debye]' /
DATA nitrogen(14) / 'OT0                !default reference state' /
DATA nitrogen(15) / '298.15  101.325  8670.0  191.5  !tref, Pref, Href, Sref' /
DATA nitrogen(16) / '8.0                !version number' /
DATA nitrogen(17) / '1066, 1977         !UN Number' /
DATA nitrogen(18) / 'other              !family' /
DATA nitrogen(19) / '0.0                !heating value (gross or superior) [kJ/mol]' /
DATA nitrogen(20) / 'A1                 !Safety Group (ASHRAE Standard 34, 2010)' /
DATA nitrogen(21) / '' /
DATA nitrogen(22) / '' /
DATA nitrogen(23) / '#EOS               !equation of state specification' /
DATA nitrogen(24) / 'FEQ  Helmholtz equation of state for nitrogen of Span et al. (2000).' /
DATA nitrogen(25) / '?LITERATURE REFERENCE \' /
DATA nitrogen(26) / '?Span, R., Lemmon, E.W., Jacobsen, R.T, Wagner, W., and Yokozeki, A.' /
DATA nitrogen(27) / '?"A Reference Equation of State for the Thermodynamic Properties of' /
DATA nitrogen(28) / '? Nitrogen for Temperatures from 63.151 to 1000 K and Pressures to 2200 MPa,"' /
DATA nitrogen(29) / '? J. Phys. Chem. Ref. Data, 29(6):1361-1433, 2000.' /
DATA nitrogen(30) / '?\' /
DATA nitrogen(31) / '? see also Int. J. Thermophys., 19(4):1121-1132, 1998.' /
DATA nitrogen(32) / '?\' /
DATA nitrogen(33) / '?The uncertainty in density of the equation of state is 0.02% from the' /
DATA nitrogen(34) / '?triple point up to temperatures of 523 K and pressures up to 12 MPa and' /
DATA nitrogen(35) / '?from temperatures of 240 to 523 K at pressures less than 30 MPa. In the' /
DATA nitrogen(36) / '?range from 270 to 350 K at pressures less than 12 MPa, the uncertainty' /
DATA nitrogen(37) / '?in density is 0.01%.  The uncertainty at very high pressures (>1 GPa) is' /
DATA nitrogen(38) / '?0.6% in density.  The uncertainty in pressure in the critical region is' /
DATA nitrogen(39) / '?estimated to be 0.02%.  In the gaseous and supercritical region, the' /
DATA nitrogen(40) / '?speed of sound can be calculated with a typical uncertainty of 0.005% to' /
DATA nitrogen(41) / '?0.1%.  At liquid states and at high pressures, the uncertainty increases' /
DATA nitrogen(42) / '?to 0.5% - 1.5%.  For pressures up to 30 MPa, the estimated uncertainty' /
DATA nitrogen(43) / '?for heat capacities ranges from 0.3% at gaseous and gas like supercritical' /
DATA nitrogen(44) / '?states up to 0.8% at liquid states and at certain gaseous and supercritical' /
DATA nitrogen(45) / '?states at low temperatures.  The uncertainty is 2% for pressures up to' /
DATA nitrogen(46) / '?200 MPa and larger at higher pressures.  The estimated uncertainties of' /
DATA nitrogen(47) / '?vapor pressure, saturated liquid density, and saturated vapor density' /
DATA nitrogen(48) / '?are in general 0.02% for each property.  The formulation yields a' /
DATA nitrogen(49) / '?reasonable extrapolation behavior up to the limits of chemical stability' /
DATA nitrogen(50) / '?of nitrogen.' /
DATA nitrogen(51) / '?\' /
DATA nitrogen(52) / '!end of info section' /
DATA nitrogen(53) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(54) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(55) / '2200000.0          !upper pressure limit [kPa]' /
DATA nitrogen(56) / '53.15              !maximum density [mol/L]' /
DATA nitrogen(57) / 'CPP                                    !pointer to Cp0 model' /
DATA nitrogen(58) / '28.01348                               !molecular weight [g/mol]' /
DATA nitrogen(59) / '63.151                                 !triple point temperature [K]' /
DATA nitrogen(60) / '12.5198                                !pressure at triple point [kPa]' /
DATA nitrogen(61) / '30.957                                 !density at triple point [mol/L]' /
DATA nitrogen(62) / '77.3550                                !normal boiling point temperature [K]' /
DATA nitrogen(63) / '0.0372                                 !acentric factor' /
DATA nitrogen(64) / '126.192      3395.8       11.1839      !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA nitrogen(65) / '126.192                   11.1839      !reducing parameters [K, mol/L]' /
DATA nitrogen(66) / '8.31451                                !gas constant [J/mol-K]' /
DATA nitrogen(67) / '      32  4      4  12      0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA nitrogen(68) / ' 0.924803575275d+00  0.250   1.00    0 !a(i),t(i),d(i),l(i)' /
DATA nitrogen(69) / '-0.492448489428d+00  0.875   1.00    0' /
DATA nitrogen(70) / ' 0.661883336938d+00  0.500   2.00    0' /
DATA nitrogen(71) / '-0.192902649201d+01  0.875   2.00    0' /
DATA nitrogen(72) / '-0.622469309629d-01  0.375   3.00    0' /
DATA nitrogen(73) / ' 0.349943957581d+00  0.750   3.00    0' /
DATA nitrogen(74) / ' 0.564857472498d+00  0.500   1.00    1' /
DATA nitrogen(75) / '-0.161720005987d+01  0.750   1.00    1' /
DATA nitrogen(76) / '-0.481395031883d+00  2.000   1.00    1' /
DATA nitrogen(77) / ' 0.421150636384d+00  1.250   3.00    1' /
DATA nitrogen(78) / '-0.161962230825d-01  3.500   3.00    1' /
DATA nitrogen(79) / ' 0.172100994165d+00  1.000   4.00    1' /
DATA nitrogen(80) / ' 0.735448924933d-02  0.500   6.00    1' /
DATA nitrogen(81) / ' 0.168077305479d-01  3.000   6.00    1' /
DATA nitrogen(82) / '-0.107626664179d-02  0.000   7.00    1' /
DATA nitrogen(83) / '-0.137318088513d-01  2.750   7.00    1' /
DATA nitrogen(84) / ' 0.635466899859d-03  0.750   8.00    1' /
DATA nitrogen(85) / ' 0.304432279419d-02  2.500   8.00    1' /
DATA nitrogen(86) / '-0.435762336045d-01  4.000   1.00    2' /
DATA nitrogen(87) / '-0.723174889316d-01  6.000   2.00    2' /
DATA nitrogen(88) / ' 0.389644315272d-01  6.000   3.00    2' /
DATA nitrogen(89) / '-0.212201363910d-01  3.000   4.00    2' /
DATA nitrogen(90) / ' 0.408822981509d-02  3.000   5.00    2' /
DATA nitrogen(91) / '-0.551990017984d-04  6.000   8.00    2' /
DATA nitrogen(92) / '-0.462016716479d-01 16.000   4.00    3' /
DATA nitrogen(93) / '-0.300311716011d-02 11.000   5.00    3' /
DATA nitrogen(94) / ' 0.368825891208d-01 15.000   5.00    3' /
DATA nitrogen(95) / '-0.255856846220d-02 12.000   8.00    3' /
DATA nitrogen(96) / ' 0.896915264558d-02 12.000   3.00    4' /
DATA nitrogen(97) / '-0.441513370350d-02  7.000   5.00    4' /
DATA nitrogen(98) / ' 0.133722924858d-02  4.000   6.00    4' /
DATA nitrogen(99) / ' 0.264832491957d-03 16.000   9.00    4' /
DATA nitrogen(100) / ' 0.196688194015d+02  0.000   1.00    2 2 -20.0  -325.0  1.16  1.  0.  0.  0.' /
DATA nitrogen(101) / '-0.209115600730d+02  1.000   1.00    2 2 -20.0  -325.0  1.16  1.  0.  0.  0.' /
DATA nitrogen(102) / ' 0.167788306989d-01  2.000   3.00    2 2 -15.0  -300.0  1.13  1.  0.  0.  0.' /
DATA nitrogen(103) / ' 0.262767566274d+04  3.000   2.00    2 2 -25.0  -275.0  1.25  1.  0.  0.  0.' /
DATA nitrogen(104) / '' /
DATA nitrogen(105) / '' /
DATA nitrogen(106) / '#AUX               !auxiliary model specification' /
DATA nitrogen(107) / 'CPP  ideal gas heat capacity function' /
DATA nitrogen(108) / '?LITERATURE REFERENCE \' /
DATA nitrogen(109) / '?Span, R., Lemmon, E.W., Jacobsen, R.T, Wagner, W., and Yokozeki, A.' /
DATA nitrogen(110) / '?"A Reference Equation of State for the Thermodynamic Properties of' /
DATA nitrogen(111) / '? Nitrogen for Temperatures from 63.151 to 1000 K and Pressures to 2200 MPa,"' /
DATA nitrogen(112) / '? J. Phys. Chem. Ref. Data, 29(6):1361-1433, 2000.' /
DATA nitrogen(113) / '?\' /
DATA nitrogen(114) / '? see also Int. J. Thermophys., 19(4):1121-1132, 1998.' /
DATA nitrogen(115) / '?\' /
DATA nitrogen(116) / '!end of info section' /
DATA nitrogen(117) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(118) / '7000.0             !upper temperature limit [K]' /
DATA nitrogen(119) / '0.0                !upper pressure limit [kPa]' /
DATA nitrogen(120) / '0.0                !maximum density [mol/L]' /
DATA nitrogen(121) / '1.0          8.31451                   !reducing parameters for T, Cp0' /
DATA nitrogen(122) / '  4  1    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA nitrogen(123) / ' 3.5               0.00' /
DATA nitrogen(124) / ' 3.066469d-6       1.00' /
DATA nitrogen(125) / ' 4.70124d-9        2.00' /
DATA nitrogen(126) / '-3.987984d-13      3.00' /
DATA nitrogen(127) / ' 0.1012941d1      3364.011' /
DATA nitrogen(128) / '' /
DATA nitrogen(129) / '' /
DATA nitrogen(130) / '@EOS               !equation of state specification' /
DATA nitrogen(131) / 'FEK  Helmholtz equation of state for nitrogen of Kunz and Wagner (2004).' /
DATA nitrogen(132) / '?LITERATURE REFERENCE \' /
DATA nitrogen(133) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA nitrogen(134) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA nitrogen(135) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA nitrogen(136) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA nitrogen(137) / '?\' /
DATA nitrogen(138) / '!end of info section' /
DATA nitrogen(139) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(140) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(141) / '2200000.0          !upper pressure limit [kPa]' /
DATA nitrogen(142) / '53.15              !maximum density [mol/L]' /
DATA nitrogen(143) / 'PHK                                    !pointer to Cp0 model' /
DATA nitrogen(144) / '28.0134                                !molecular weight [g/mol]' /
DATA nitrogen(145) / '63.151                                 !triple point temperature [K]' /
DATA nitrogen(146) / '12.523                                 !pressure at triple point [kPa]' /
DATA nitrogen(147) / '30.954                                 !density at triple point [mol/L]' /
DATA nitrogen(148) / ' 77.36                                 !normal boiling point temperature [K]' /
DATA nitrogen(149) / ' 0.0373                                !acentric factor' /
DATA nitrogen(150) / '126.192      3395.8      11.1839       !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA nitrogen(151) / '126.192                  11.1839       !reducing parameters [K, mol/L]' /
DATA nitrogen(152) / '8.314472                               !gas constant [J/mol-K]' /
DATA nitrogen(153) / '  24  4      0  0      0  0            !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA nitrogen(154) / ' 0.59889711801201       0.125  1.  0' /
DATA nitrogen(155) / '-0.16941557480731d1     1.125  1.  0' /
DATA nitrogen(156) / ' 0.24579736191718       0.375  2.  0' /
DATA nitrogen(157) / '-0.23722456755175       1.125  2.  0' /
DATA nitrogen(158) / ' 0.17954918715141d-1    0.625  4.  0' /
DATA nitrogen(159) / ' 0.14592875720215d-1    1.5    4.  0' /
DATA nitrogen(160) / ' 0.10008065936206       0.625  1.  1' /
DATA nitrogen(161) / ' 0.73157115385532       2.625  1.  1' /
DATA nitrogen(162) / '-0.88372272336366       2.75   1.  1' /
DATA nitrogen(163) / ' 0.31887660246708       2.125  2.  1' /
DATA nitrogen(164) / ' 0.20766491728799       2.0    3.  1' /
DATA nitrogen(165) / '-0.19379315454158d-1    1.75   6.  1' /
DATA nitrogen(166) / '-0.16936641554983       4.50   2.  2' /
DATA nitrogen(167) / ' 0.13546846041701       4.75   3.  2' /
DATA nitrogen(168) / '-0.33066712095307d-1    5.0    3.  2' /
DATA nitrogen(169) / '-0.60690817018557d-1    4.0    4.  2' /
DATA nitrogen(170) / ' 0.12797548292871d-1    4.5    4.  2' /
DATA nitrogen(171) / ' 0.58743664107299d-2    7.5    2.  3' /
DATA nitrogen(172) / '-0.18451951971969d-1    14.0   3.  3' /
DATA nitrogen(173) / ' 0.47226622042472d-2    11.5   4.  3' /
DATA nitrogen(174) / '-0.52024079680599d-2    26.0   5.  6' /
DATA nitrogen(175) / ' 0.43563505956635d-1    28.0   6.  6' /
DATA nitrogen(176) / '-0.36251690750939d-1    30.0   6.  6' /
DATA nitrogen(177) / '-0.28974026866543d-2    16.0   7.  6' /
DATA nitrogen(178) / '' /
DATA nitrogen(179) / '' /
DATA nitrogen(180) / '#AUX               !auxiliary model specification' /
DATA nitrogen(181) / 'PHK  Helmholtz form for the ideal-gas state for nitrogen of Kunz and Wagner (2004).' /
DATA nitrogen(182) / '?LITERATURE REFERENCE \' /
DATA nitrogen(183) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA nitrogen(184) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA nitrogen(185) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA nitrogen(186) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA nitrogen(187) / '?\' /
DATA nitrogen(188) / '!end of info section' /
DATA nitrogen(189) / '0.                 !lower temperature limit [K]' /
DATA nitrogen(190) / '1000.0             !upper temperature limit [K]' /
DATA nitrogen(191) / '0.0                !upper pressure limit [kPa]' /
DATA nitrogen(192) / '0.0                !maximum density [mol/L]' /
DATA nitrogen(193) / '1 2  0  1 2  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA nitrogen(194) / '    2.50031      1.             !ai, ti for [ai*log(tau**ti)] terms' /
DATA nitrogen(195) / '   11.083407489  0.             !aj, ti for [ai*tau**ti] terms' /
DATA nitrogen(196) / '  -22.202102428  1.' /
DATA nitrogen(197) / '    0.1466      -5.393067706    !aj, ti for cosh and sinh terms' /
DATA nitrogen(198) / '    0.13732      5.25182262' /
DATA nitrogen(199) / '    0.90066     13.788988208' /
DATA nitrogen(200) / '' /
DATA nitrogen(201) / '' /
DATA nitrogen(202) / '@EOS               !equation of state specification' /
DATA nitrogen(203) / 'FE1  Helmholtz equation of state for nitrogen of Jacobsen et al. (1986).' /
DATA nitrogen(204) / '?LITERATURE REFERENCE \' /
DATA nitrogen(205) / '?Jacobsen, R.T, Stewart, R.B., and Jahangiri, M.,' /
DATA nitrogen(206) / '? "Thermodynamic properties of nitrogen from the freezing line to 2000 K at' /
DATA nitrogen(207) / '? pressures to 1000 MPa,"' /
DATA nitrogen(208) / '? J. Phys. Chem. Ref. Data, 15(2):735-909, 1986.' /
DATA nitrogen(209) / '?\' /
DATA nitrogen(210) / '!end of info section' /
DATA nitrogen(211) / '63.148             !lower temperature limit [K]' /
DATA nitrogen(212) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(213) / '1000000.0          !upper pressure limit [kPa]' /
DATA nitrogen(214) / '30.96              !maximum density [mol/L]' /
DATA nitrogen(215) / 'CP1                                    !pointer to Cp0 model' /
DATA nitrogen(216) / '28.0134                                !molecular weight [g/mol]' /
DATA nitrogen(217) / '63.148                                 !triple point temperature [K]' /
DATA nitrogen(218) / '12.52                                  !pressure at triple point [kPa]' /
DATA nitrogen(219) / '31.046                                 !density at triple point [mol/L]' /
DATA nitrogen(220) / ' 77.348                                !normal boiling point temperature [K]' /
DATA nitrogen(221) / '0.03701                                !acentric factor' /
DATA nitrogen(222) / '126.193      3397.8       11.177       !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA nitrogen(223) / '126.193                   11.177       !reducing parameters [K, mol/L]' /
DATA nitrogen(224) / '8.31434                                !gas constant [J/mol-K]' /
DATA nitrogen(225) / '      28  4      0  0       0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA nitrogen(226) / '  0.949954182700d+00   0.250    1.0  0 !a(i),t(i),d(i),l(i)' /
DATA nitrogen(227) / '  0.248171851300d+00   0.250    2.0  0' /
DATA nitrogen(228) / ' -0.204628712200d+00   0.250    3.0  0' /
DATA nitrogen(229) / ' -0.174842900800d+00   0.500    2.0  0' /
DATA nitrogen(230) / '  0.638701714800d+00   0.500    3.0  0' /
DATA nitrogen(231) / ' -0.527298616800d+00   0.750    3.0  0' /
DATA nitrogen(232) / ' -0.204974150400d+01   1.000    1.0  0' /
DATA nitrogen(233) / '  0.555138355300d-01   1.000    4.0  0' /
DATA nitrogen(234) / ' -0.819110639600d-03   1.000    6.0  0' /
DATA nitrogen(235) / ' -0.503251969900d-01   1.000    2.0  2' /
DATA nitrogen(236) / '  0.265011079800d+00   1.500    1.0  0' /
DATA nitrogen(237) / '  0.731145937200d-01   2.000    2.0  0' /
DATA nitrogen(238) / ' -0.281308071800d-01   2.000    4.0  0' /
DATA nitrogen(239) / '  0.165982356900d-02   2.000    6.0  0' /
DATA nitrogen(240) / '  0.601281781200d-01   2.000    2.0  2' /
DATA nitrogen(241) / ' -0.378544519400d+00   3.000    1.0  0' /
DATA nitrogen(242) / '  0.189529043300d+00   3.000    2.0  0' /
DATA nitrogen(243) / ' -0.700189509300d-02   3.000    4.0  0' /
DATA nitrogen(244) / ' -0.492771092700d-01   3.000    1.0  3' /
DATA nitrogen(245) / '  0.651201367900d-01   4.000    4.0  2' /
DATA nitrogen(246) / '  0.113812194200d+00   4.000    1.0  3' /
DATA nitrogen(247) / ' -0.955140963197d-01   5.000    2.0  2' /
DATA nitrogen(248) / '  0.211835414000d-01   6.000    4.0  2' /
DATA nitrogen(249) / ' -0.110072177100d-01   8.000    2.0  4' /
DATA nitrogen(250) / '  0.128443221000d-01  14.000    4.0  4' /
DATA nitrogen(251) / ' -0.105447491000d-01  18.000    4.0  4' /
DATA nitrogen(252) / ' -0.148460053800d-03  20.000    2.0  4' /
DATA nitrogen(253) / ' -0.580648346700d-02  22.000    3.0  3' /
DATA nitrogen(254) / '' /
DATA nitrogen(255) / '' /
DATA nitrogen(256) / '#AUX               !auxiliary model specification' /
DATA nitrogen(257) / 'CP1  ideal gas heat capacity function of Jacobsen et al.' /
DATA nitrogen(258) / '?LITERATURE REFERENCE \' /
DATA nitrogen(259) / '?Jacobsen, R.T, Stewart, R.B., and Jahangiri, M.,' /
DATA nitrogen(260) / '? "Thermodynamic properties of nitrogen from the freezing line to 2000 K at' /
DATA nitrogen(261) / '? pressures to 1000 MPa,"' /
DATA nitrogen(262) / '? J. Phys. Chem. Ref. Data, 15(2):735-909, 1986.' /
DATA nitrogen(263) / '?\' /
DATA nitrogen(264) / '!end of info section' /
DATA nitrogen(265) / '63.148             !lower temperature limit [K]' /
DATA nitrogen(266) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(267) / '0.0                !upper pressure limit [kPa]' /
DATA nitrogen(268) / '0.0                !maximum density [mol/L]' /
DATA nitrogen(269) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA nitrogen(270) / '  7  1    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA nitrogen(271) / '-0.837079888737d+03      -3.00         !Ni, power in T' /
DATA nitrogen(272) / ' 0.379147114487d+02      -2.00' /
DATA nitrogen(273) / '-0.601737844275d+00      -1.00' /
DATA nitrogen(274) / ' 0.350418363823d+01       0.00' /
DATA nitrogen(275) / '-0.874955653028d-05       1.00' /
DATA nitrogen(276) / ' 0.148958507239d-07       2.00' /
DATA nitrogen(277) / '-0.256370354277d-11       3.00' /
DATA nitrogen(278) / ' 0.100773735767d+01       0.33534061d4  !exponential term' /
DATA nitrogen(279) / '' /
DATA nitrogen(280) / '' /
DATA nitrogen(281) / '@EOS               !equation of state specification' /
DATA nitrogen(282) / 'FES  short Helmholtz equation of state for nitrogen of Span and Wagner (2003).' /
DATA nitrogen(283) / '?LITERATURE REFERENCE \' /
DATA nitrogen(284) / '?Span, R. and Wagner, W.' /
DATA nitrogen(285) / '? "Equations of State for Technical Applications. II. Results for Nonpolar Fluids,"' /
DATA nitrogen(286) / '? Int. J. Thermophys., 24(1):41-109, 2003.' /
DATA nitrogen(287) / '?\' /
DATA nitrogen(288) / '?The uncertainties of the equation of state are approximately 0.2% (to' /
DATA nitrogen(289) / '?0.5% at high pressures) in density, 1% (in the vapor phase) to 2% in' /
DATA nitrogen(290) / '?heat capacity, 1% (in the vapor phase) to 2% in the speed of sound, and' /
DATA nitrogen(291) / '?0.2% in vapor pressure, except in the critical region.' /
DATA nitrogen(292) / '?\' /
DATA nitrogen(293) / '!end of info section' /
DATA nitrogen(294) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(295) / '600.0              !upper temperature limit [K]' /
DATA nitrogen(296) / '100000.0           !upper pressure limit [kPa]' /
DATA nitrogen(297) / '53.15              !maximum density [mol/L]' /
DATA nitrogen(298) / 'CPP                                    !pointer to Cp0 model' /
DATA nitrogen(299) / '28.013                                 !molecular weight [g/mol]' /
DATA nitrogen(300) / '63.151                                 !triple point temperature [K]' /
DATA nitrogen(301) / '12.566                                 !pressure at triple point [kPa]' /
DATA nitrogen(302) / '30.935                                 !density at triple point [mol/L]' /
DATA nitrogen(303) / '77.356                                 !normal boiling point temperature [K]' /
DATA nitrogen(304) / '0.037                                  !acentric factor' /
DATA nitrogen(305) / '126.192      3396.0       11.1839      !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA nitrogen(306) / '126.192                   11.1839      !reducing parameters [K, mol/L]' /
DATA nitrogen(307) / '8.31451                                !gas constant [J/mol-K]' /
DATA nitrogen(308) / '      12  4      0  0       0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA nitrogen(309) / ' 0.922965670000E+00  0.25    1.0     0 !a(i),t(i),d(i),l(i)' /
DATA nitrogen(310) / '-0.255750120000E+01  1.125   1.0     0' /
DATA nitrogen(311) / ' 0.644824630000E+00  1.5     1.0     0' /
DATA nitrogen(312) / ' 0.108310200000E-01  1.375   2.0     0' /
DATA nitrogen(313) / ' 0.739241670000E-01  0.25    3.0     0' /
DATA nitrogen(314) / ' 0.235329620000E-03  0.875   7.0     0' /
DATA nitrogen(315) / ' 0.180248540000E+00  0.625   2.0     1' /
DATA nitrogen(316) / '-0.456602990000E-01  1.75    5.0     1' /
DATA nitrogen(317) / '-0.155210600000E+00  3.625   1.0     2' /
DATA nitrogen(318) / '-0.381114900000E-01  3.625   4.0     2' /
DATA nitrogen(319) / '-0.319624220000E-01 14.5     3.0     3' /
DATA nitrogen(320) / ' 0.155135320000E-01 12.0     4.0     3' /
DATA nitrogen(321) / '' /
DATA nitrogen(322) / '' /
DATA nitrogen(323) / '@EOS               !equation of state specification' /
DATA nitrogen(324) / 'BWR  MBWR equation of state for nitrogen of Younglove (1982).' /
DATA nitrogen(325) / '?LITERATURE REFERENCE \' /
DATA nitrogen(326) / '?Younglove, B.A.,' /
DATA nitrogen(327) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA nitrogen(328) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA nitrogen(329) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA nitrogen(330) / '?\' /
DATA nitrogen(331) / '!end of info section' /
DATA nitrogen(332) / '63.15              !lower temperature limit [K]' /
DATA nitrogen(333) / '1900.0             !upper temperature limit [K]' /
DATA nitrogen(334) / '1013000.0          !upper pressure limit [kPa]' /
DATA nitrogen(335) / '30.977             !maximum density [mol/L]' /
DATA nitrogen(336) / 'CP2                                    !pointer to Cp0 model' /
DATA nitrogen(337) / '28.013                                 !molecular weight [g/mol]' /
DATA nitrogen(338) / '63.15                                  !triple point temperature [K]' /
DATA nitrogen(339) / '12.463                                 !pressure at triple point [kPa]' /
DATA nitrogen(340) / '30.977                                 !density at triple point [mol/L]' /
DATA nitrogen(341) / '77.348                                 !normal boiling point temperature [K]' /
DATA nitrogen(342) / '0.03701                                !acentric factor' /
DATA nitrogen(343) / '126.26       3399.08      11.21        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA nitrogen(344) / '126.26                    11.21        !reducing parameters [K, mol/L]' /
DATA nitrogen(345) / '13.3630620956                          !gamma' /
DATA nitrogen(346) / '0.0831411                              !gas constant [L-bar/mol-K]' /
DATA nitrogen(347) / '      32       1                       !Nterm, Ncoeff per term' /
DATA nitrogen(348) / '  0.1380297474657d-02      0.1084506501349d-00     -0.2471324064362d+01' /
DATA nitrogen(349) / '  0.3455257980807d+02     -0.4279707690666d+04      0.1064911566998d-03' /
DATA nitrogen(350) / ' -0.1140867079735d-01      0.1444902497287d-03      0.1871457567553d+05' /
DATA nitrogen(351) / '  0.8218876886831d-07      0.2360990493348d-02     -0.5144803081201d-00' /
DATA nitrogen(352) / '  0.4914545013668d-04     -0.1151627162399d-02     -0.7168037246650d-00' /
DATA nitrogen(353) / '  0.7616667619500d-04     -0.1130930066213d-05      0.3736831166831d-03' /
DATA nitrogen(354) / ' -0.2039851507581d-05     -0.1719662008990d+05     -0.1213055199748d+06' /
DATA nitrogen(355) / ' -0.9881399141428d+02      0.5619886893511d+05     -0.1823043964118d-00' /
DATA nitrogen(356) / ' -0.2599826498477d+01     -0.4191893423157d-03     -0.2596406670530d-00' /
DATA nitrogen(357) / ' -0.1258683201921d-06      0.1049286599400d-04     -0.5458369305152d-09' /
DATA nitrogen(358) / ' -0.7674511670597d-08      0.5931232870994d-07' /
DATA nitrogen(359) / '' /
DATA nitrogen(360) / '' /
DATA nitrogen(361) / '#AUX               !auxiliary model specification' /
DATA nitrogen(362) / 'CP2  ideal gas heat capacity function of Younglove' /
DATA nitrogen(363) / '?LITERATURE REFERENCE \' /
DATA nitrogen(364) / '?Younglove, B.A.,' /
DATA nitrogen(365) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA nitrogen(366) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA nitrogen(367) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA nitrogen(368) / '?\' /
DATA nitrogen(369) / '!end of info section' /
DATA nitrogen(370) / '63.15              !lower temperature limit [K]' /
DATA nitrogen(371) / '1900.0             !upper temperature limit [K]' /
DATA nitrogen(372) / '0.0                !upper pressure limit [kPa]' /
DATA nitrogen(373) / '0.0                !maximum density [mol/L]' /
DATA nitrogen(374) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA nitrogen(375) / '  7  1    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA nitrogen(376) / ' -0.7352104011573d+03    -3.00d0' /
DATA nitrogen(377) / '  0.3422399804120d+02    -2.00d0' /
DATA nitrogen(378) / ' -0.5576482845676d+00    -1.00d0' /
DATA nitrogen(379) / '  0.3504042283088d+01     0.00d0' /
DATA nitrogen(380) / ' -0.1733901850810d-04     1.00d0' /
DATA nitrogen(381) / '  0.1746508497665d-07     2.00d0' /
DATA nitrogen(382) / ' -0.3568920335443d-11     3.00d0' /
DATA nitrogen(383) / '  0.1005387228088d+01  3353.4061' /
DATA nitrogen(384) / '' /
DATA nitrogen(385) / '' /
DATA nitrogen(386) / '' /
DATA nitrogen(387) / '#TCX               !thermal conductivity model specification' /
DATA nitrogen(388) / 'TC1  pure fluid thermal conductivity model of Lemmon and Jacobsen (2004).' /
DATA nitrogen(389) / '?LITERATURE REFERENCE \' /
DATA nitrogen(390) / '?Lemmon, E.W. and Jacobsen, R.T,' /
DATA nitrogen(391) / '? "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon,' /
DATA nitrogen(392) / '? and Air,"' /
DATA nitrogen(393) / '? Int. J. Thermophys., 25:21-69, 2004.' /
DATA nitrogen(394) / '?\' /
DATA nitrogen(395) / '?The uncertainty for the dilute gas is 2% with increasing uncertainties' /
DATA nitrogen(396) / '?near the triple point.  For the non-dilute gas, the uncertainty is 2%' /
DATA nitrogen(397) / '?for temperatures greater than 150 K. The uncertainty is 3% at' /
DATA nitrogen(398) / '?temperatures less than the critical point and 5% in the critical region,' /
DATA nitrogen(399) / '?except for states very near the critical point.' /
DATA nitrogen(400) / '?\' /
DATA nitrogen(401) / '!end of info section' /
DATA nitrogen(402) / '50.0               !lower temperature limit [K]' /
DATA nitrogen(403) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(404) / '2200000.0          !upper pressure limit [kPa]' /
DATA nitrogen(405) / '53.15              !maximum density [mol/L]' /
DATA nitrogen(406) / '3   0              !# terms for dilute gas function:  numerator, denominator' /
DATA nitrogen(407) / '126.192   1.0d-3   !reducing parameters for T, tcx' /
DATA nitrogen(408) / ' 1.511   -97.0     !coeff, power in T' /
DATA nitrogen(409) / ' 2.117     1.0' /
DATA nitrogen(410) / '-3.332     0.7' /
DATA nitrogen(411) / '6   0              !# terms for background gas function:  numerator, denominator' /
DATA nitrogen(412) / '126.192   11.1839     1.0d-3    !reducing parameters for T, rho, tcx' /
DATA nitrogen(413) / ' 8.862          0.0  1.0  0.0 !coeff, powers of T, rho, exp(rho)' /
DATA nitrogen(414) / ' 31.11       -0.03   2.   0.' /
DATA nitrogen(415) / '-73.13       -0.2    3.   1.' /
DATA nitrogen(416) / ' 20.03       -0.8    4.   2.' /
DATA nitrogen(417) / '-0.7096      -0.6    8.   2.' /
DATA nitrogen(418) / ' 0.2672      -1.9   10.   2.' /
DATA nitrogen(419) / 'TK3                !pointer to critical enhancement auxiliary function' /
DATA nitrogen(420) / '' /
DATA nitrogen(421) / '' /
DATA nitrogen(422) / '#AUX               !thermal conductivity critical enhancement model' /
DATA nitrogen(423) / 'TK3  thermal conductivity critical enhancement of Lemmon and Jacobsen (2004).' /
DATA nitrogen(424) / '?LITERATURE REFERENCE \' /
DATA nitrogen(425) / '?\' /
DATA nitrogen(426) / '!end of info section' /
DATA nitrogen(427) / '50.0               !lower temperature limit [K]' /
DATA nitrogen(428) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(429) / '2200000.0          !upper pressure limit [kPa]' /
DATA nitrogen(430) / '53.15              !maximum density [mol/L]' /
DATA nitrogen(431) / '9  0  0  0         !# terms:  terms, spare, spare, spare' /
DATA nitrogen(432) / '1.0    1.0  1.0    !reducing par for T, rho, tcx (mW/m-K)' /
DATA nitrogen(433) / '0.630d0            !gnu (universal exponent)' /
DATA nitrogen(434) / '1.2415d0           !gamma (universal exponent)' /
DATA nitrogen(435) / '1.01d0             !R0 (universal amplitude)' /
DATA nitrogen(436) / ' 0.065d0           !z (universal exponent--not used for t.c., only viscosity)' /
DATA nitrogen(437) / ' 1.00d0            !c (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA nitrogen(438) / ' 0.17E-09          !xi0 (amplitude) [m]' /
DATA nitrogen(439) / ' 0.55E-01          !gam0 (amplitude) [-]' /
DATA nitrogen(440) / ' 0.40E-09          !qd_inverse (modified effective cutoff parameter) [m]' /
DATA nitrogen(441) / '252.384            !tref (reference temperature) [K]' /
DATA nitrogen(442) / '' /
DATA nitrogen(443) / '' /
DATA nitrogen(444) / '#ETA               !viscosity model specification' /
DATA nitrogen(445) / 'VS1  pure fluid viscosity model of Lemmon and Jacobsen (2004).' /
DATA nitrogen(446) / '?LITERATURE REFERENCE \' /
DATA nitrogen(447) / '?Lemmon, E.W. and Jacobsen, R.T,' /
DATA nitrogen(448) / '? "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon,' /
DATA nitrogen(449) / '? and Air,"' /
DATA nitrogen(450) / '? Int. J. Thermophys., 25:21-69, 2004.' /
DATA nitrogen(451) / '?\' /
DATA nitrogen(452) / '?The uncertainty is 0.5% in the dilute gas.  Away from the dilute gas' /
DATA nitrogen(453) / '?(pressures greater than 1 MPa and in the liquid), the uncertainties are' /
DATA nitrogen(454) / '?as low as 1% between 270 and 300 K at pressures less than 100 MPa, and' /
DATA nitrogen(455) / '?increase outside that range.  The uncertainties are around 2% at' /
DATA nitrogen(456) / '?temperatures of 180 K and higher.  Below this and away from the critical' /
DATA nitrogen(457) / '?region, the uncertainties steadily increase to around 5% at the triple' /
DATA nitrogen(458) / '?points of the fluids.  The uncertainties in the critical region are' /
DATA nitrogen(459) / '?higher.' /
DATA nitrogen(460) / '!end of info section' /
DATA nitrogen(461) / '50.0               !lower temperature limit [K]' /
DATA nitrogen(462) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(463) / '2200000.0          !upper pressure limit [kPa]' /
DATA nitrogen(464) / '53.15              !maximum density [mol/L]' /
DATA nitrogen(465) / '1                  !number of terms associated with dilute-gas function' /
DATA nitrogen(466) / 'CI1                !pointer to reduced effective collision cross-section model' /
DATA nitrogen(467) / '0.3656             !Lennard-Jones coefficient sigma [nm]' /
DATA nitrogen(468) / '98.94              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA nitrogen(469) / '1.0    1.0         !reducing parameters for T, eta' /
DATA nitrogen(470) / '0.141294895  0.5   !Chapman-Enskog term' /
DATA nitrogen(471) / '0                  !number of terms for initial density dependence' /
DATA nitrogen(472) / '0 5 0 0 0 0        !# resid terms:  close-packed density;  simple poly; numerator of rational poly; denominator of rat. poly; numerator of exponential; denominator of exponential' /
DATA nitrogen(473) / '126.192   11.1839     1.0           !reducing parameters for T, rho, eta' /
DATA nitrogen(474) / ' 10.72       -0.1    2.   0.   0    !simple polynomial terms' /
DATA nitrogen(475) / ' 0.03989     -0.25  10.   0.   1' /
DATA nitrogen(476) / ' 0.001208    -3.2   12.   0.   1' /
DATA nitrogen(477) / '-7.402       -0.9    2.   0.   2' /
DATA nitrogen(478) / ' 4.620       -0.3    1.   0.   3' /
DATA nitrogen(479) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA nitrogen(480) / '' /
DATA nitrogen(481) / '' /
DATA nitrogen(482) / '#AUX               !collision integral specification' /
DATA nitrogen(483) / 'CI1  collision integral model of Lemmon and Jacobsen (2004).' /
DATA nitrogen(484) / '?LITERATURE REFERENCE \' /
DATA nitrogen(485) / '?\' /
DATA nitrogen(486) / '!end of info section' /
DATA nitrogen(487) / '1.0                !lower temperature limit [K]' /
DATA nitrogen(488) / '10000.0            !upper temperature limit [K]' /
DATA nitrogen(489) / '0.0                !(dummy) upper pressure limit' /
DATA nitrogen(490) / '0.0                !(dummy) maximum density' /
DATA nitrogen(491) / '5                  !number of terms' /
DATA nitrogen(492) / '  0.431      0     !coeff, power of Tstar' /
DATA nitrogen(493) / ' -0.4623     1' /
DATA nitrogen(494) / '  0.08406    2' /
DATA nitrogen(495) / '  0.005341   3' /
DATA nitrogen(496) / ' -0.00331    4' /
DATA nitrogen(497) / '' /
DATA nitrogen(498) / '' /
DATA nitrogen(499) / '@TCX               !thermal conductivity model specification' /
DATA nitrogen(500) / 'TC3  pure fluid thermal conductivity model of Younglove (1982).' /
DATA nitrogen(501) / '?LITERATURE REFERENCE \' /
DATA nitrogen(502) / '?Younglove, B.A.,' /
DATA nitrogen(503) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA nitrogen(504) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA nitrogen(505) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA nitrogen(506) / '?\' /
DATA nitrogen(507) / '!end of info section' /
DATA nitrogen(508) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(509) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(510) / '2200000.0          !upper pressure limit [kPa]' /
DATA nitrogen(511) / '53.15              !maximum density [mol/L]' /
DATA nitrogen(512) / '0.354              !Lennard-Jones coefficient sigma [nm]' /
DATA nitrogen(513) / '118                !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA nitrogen(514) / '0.141286429751707  !const in Eq 20 = 5/16*(k*MW/1000/pi/Na)**0.5*1.0d12' /
DATA nitrogen(515) / ' 0                 !exponent in Eq 20 for T' /
DATA nitrogen(516) / '-.15055520615565   !eta0 (eq 20): coeffs of {(e/kT)**((4-n)/3)}' /
DATA nitrogen(517) / '0.183477124982509' /
DATA nitrogen(518) / ' 1.45008451566007' /
DATA nitrogen(519) / '-4.88031780663869' /
DATA nitrogen(520) / ' 6.68390592664363' /
DATA nitrogen(521) / '-4.90242883649539' /
DATA nitrogen(522) / ' 2.02630917877999' /
DATA nitrogen(523) / '-.439826733340102' /
DATA nitrogen(524) / ' 3.91906706514D-02' /
DATA nitrogen(525) / ' 1.50938067650D-03 !coeff for initial density dependence of viscosity (eq 21); Fv(1)' /
DATA nitrogen(526) / ' 1.70975795748D-04 !Fv(2)' /
DATA nitrogen(527) / ' 1.2               !Fv(3)' /
DATA nitrogen(528) / ' 118               !Fv(4)' /
DATA nitrogen(529) / '-38.613291627      !coefficients for residual viscosity, eqs (22 - 25)' /
DATA nitrogen(530) / '-31.826109485      !Ev(2)' /
DATA nitrogen(531) / ' 26.0197970589236  !Ev(3)' /
DATA nitrogen(532) / '-27.2869897441495  !Ev(4)' /
DATA nitrogen(533) / ' 0                 !Ev(5)' /
DATA nitrogen(534) / ' 0                 !Ev(6)' /
DATA nitrogen(535) / ' 0                 !Ev(7)' /
DATA nitrogen(536) / ' 35.6938892061679  !Ev(8)' /
DATA nitrogen(537) / ' 1.67108           !F' /
DATA nitrogen(538) / '0.00000003933      !rm' /
DATA nitrogen(539) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA nitrogen(540) / '' /
DATA nitrogen(541) / '' /
DATA nitrogen(542) / '@ETA               !viscosity model specification' /
DATA nitrogen(543) / 'VS2  pure fluid viscosity model of Younglove (1982).' /
DATA nitrogen(544) / '?LITERATURE REFERENCE \' /
DATA nitrogen(545) / '?Younglove, B.A.,' /
DATA nitrogen(546) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA nitrogen(547) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA nitrogen(548) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA nitrogen(549) / '?\' /
DATA nitrogen(550) / '!end of info section' /
DATA nitrogen(551) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(552) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(553) / '2200000.0          !upper pressure limit [kPa]' /
DATA nitrogen(554) / '53.15              !maximum density [mol/L]' /
DATA nitrogen(555) / 'CI2                !pointer to collision integral model' /
DATA nitrogen(556) / '0.354              !Lennard-Jones coefficient sigma [nm]' /
DATA nitrogen(557) / '118                !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA nitrogen(558) / '0.141286429751707  !const in Eq 19 = 5/16*(k*MW/1000/pi/Na)**0.5*1.0d12' /
DATA nitrogen(559) / ' 0                 !exponent in Eq 20 for T' /
DATA nitrogen(560) / '-3.14276193277D-03 !coeff for initial density dependence of viscosity (eq 21); Fv(1)' /
DATA nitrogen(561) / ' 9.22071479907D-04 !Fv(2)' /
DATA nitrogen(562) / ' 1.4               !Fv(3)' /
DATA nitrogen(563) / ' 118               !Fv(4)' /
DATA nitrogen(564) / '-12.128154129      !coefficients for residual viscosity, eqs (22 - 25)' /
DATA nitrogen(565) / ' 68.46443564       !Ev(2)' /
DATA nitrogen(566) / ' 11.2569594404402  !Ev(3)' /
DATA nitrogen(567) / '-565.76279020055   !Ev(4)' /
DATA nitrogen(568) / ' 9.56677570672D-02 !Ev(5)' /
DATA nitrogen(569) / '-.355533724265011  !Ev(6)' /
DATA nitrogen(570) / ' 618.536783201947  !Ev(7)' /
DATA nitrogen(571) / ' 11.2435750999429  !Ev(8)' /
DATA nitrogen(572) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA nitrogen(573) / '' /
DATA nitrogen(574) / '' /
DATA nitrogen(575) / '@AUX               !collision integral specification' /
DATA nitrogen(576) / 'CI2  collision integral model of Younglove (1982).' /
DATA nitrogen(577) / '?LITERATURE REFERENCE \' /
DATA nitrogen(578) / '?Younglove, B.A.,' /
DATA nitrogen(579) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA nitrogen(580) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA nitrogen(581) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA nitrogen(582) / '?\' /
DATA nitrogen(583) / '!end of info section' /
DATA nitrogen(584) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(585) / '625.0              !upper temperature limit [K]' /
DATA nitrogen(586) / '0.0                !(dummy) upper pressure limit' /
DATA nitrogen(587) / '0.0                !(dummy) maximum density' /
DATA nitrogen(588) / '9                  !number of terms' /
DATA nitrogen(589) / '-136.985150760851  !eta0 (eq 20): coeffs of {(e/kT)**((4-n)/3)}' /
DATA nitrogen(590) / ' 734.241371453542' /
DATA nitrogen(591) / '-1655.39131952744' /
DATA nitrogen(592) / ' 2062.67809686969' /
DATA nitrogen(593) / '-1579.52439123889' /
DATA nitrogen(594) / ' 777.942880032361' /
DATA nitrogen(595) / '-232.996787901831' /
DATA nitrogen(596) / ' 40.0691427576552' /
DATA nitrogen(597) / '-2.99482706239363' /
DATA nitrogen(598) / '' /
DATA nitrogen(599) / '' /
DATA nitrogen(600) / '#STN        !surface tension specification' /
DATA nitrogen(601) / 'ST1  surface tension model of Mulero et al. (2012)' /
DATA nitrogen(602) / '?LITERATURE REFERENCE \' /
DATA nitrogen(603) / '?Mulero, A., Cachadina, I., and Parra, M.I.' /
DATA nitrogen(604) / '? "Recommended Correlations for the Surface Tension of Common Fluids,"' /
DATA nitrogen(605) / '? J. Phys. Chem. Ref. Data, 41, 043105, 2012.' /
DATA nitrogen(606) / '?\' /
DATA nitrogen(607) / '!end of info section' /
DATA nitrogen(608) / '0.0                !lower temperature limit [K]' /
DATA nitrogen(609) / '126.192            !upper temperature limit [K]' /
DATA nitrogen(610) / '0.0                !(dummy) upper pressure limit' /
DATA nitrogen(611) / '0.0                !(dummy) maximum density' /
DATA nitrogen(612) / '1                           !number of terms in surface tension model' /
DATA nitrogen(613) / '126.192                     !critical temperature used in fit (dummy)' /
DATA nitrogen(614) / ' 0.02898     1.246          !sigma0 and n' /
DATA nitrogen(615) / '' /
DATA nitrogen(616) / '' /
DATA nitrogen(617) / '#DE         !dielectric constant specification' /
DATA nitrogen(618) / 'DE3  dielectric constant model of Harvey and Lemmon (2005).' /
DATA nitrogen(619) / '?LITERATURE REFERENCE \' /
DATA nitrogen(620) / '?Harvey, A.H. and Lemmon, E.W.' /
DATA nitrogen(621) / '? "Method for Estimating the Dielectric Constant of Natural Gas Mixtures,"' /
DATA nitrogen(622) / '? Int. J. Thermophys., 26(1):31-46, 2005.' /
DATA nitrogen(623) / '?\' /
DATA nitrogen(624) / '!end of info section' /
DATA nitrogen(625) / '0.0                !lower temperature limit [K]' /
DATA nitrogen(626) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(627) / '0.0                !(dummy) upper pressure limit' /
DATA nitrogen(628) / '0.0                !(dummy) maximum density' /
DATA nitrogen(629) / '273.16 1000.0 1.0  !reducing parameters for t and d' /
DATA nitrogen(630) / '0 2 4 0 0 0                         !number of terms in dielectric constant model' /
DATA nitrogen(631) / ' 4.3872           0.    1.    0.    !coef, t exp, d exp' /
DATA nitrogen(632) / ' 0.00226          1.    1.    0.' /
DATA nitrogen(633) / ' 2.206            0.    2.    0.' /
DATA nitrogen(634) / ' 1.135            1.    2.    0.' /
DATA nitrogen(635) / '-169.0            0.    3.1   0.' /
DATA nitrogen(636) / '-35.83            1.    3.1   0.' /
DATA nitrogen(637) / '' /
DATA nitrogen(638) / '' /
DATA nitrogen(639) / '#MLT        !melting line specification' /
DATA nitrogen(640) / 'ML1  melting line model of Span et al. (2000).' /
DATA nitrogen(641) / '?LITERATURE REFERENCE \' /
DATA nitrogen(642) / '?Span, R., Lemmon, E.W., Jacobsen, R.T, Wagner, W., and Yokozeki, A.' /
DATA nitrogen(643) / '?"A Reference Equation of State for the Thermodynamic Properties of' /
DATA nitrogen(644) / '? Nitrogen for Temperatures from 63.151 to 1000 K and Pressures to 2200 MPa,"' /
DATA nitrogen(645) / '? J. Phys. Chem. Ref. Data, 29(6):1361-1433, 2000.' /
DATA nitrogen(646) / '?\' /
DATA nitrogen(647) / '? see also Int. J. Thermophys., 19(4):1121-1132, 1998.' /
DATA nitrogen(648) / '?\' /
DATA nitrogen(649) / '!end of info section' /
DATA nitrogen(650) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(651) / '2000.0             !upper temperature limit [K]' /
DATA nitrogen(652) / '0.0                !(dummy) upper pressure limit' /
DATA nitrogen(653) / '0.0                !(dummy) maximum density' /
DATA nitrogen(654) / '63.151   12.523    !reducing temperature and pressure' /
DATA nitrogen(655) / '3 0 0 0 0 0                 !number of terms in melting line equation' /
DATA nitrogen(656) / ' 1.             0.          !coefficients and exponents' /
DATA nitrogen(657) / ' 12798.61       1.78963' /
DATA nitrogen(658) / '-12798.61       0.' /
DATA nitrogen(659) / '' /
DATA nitrogen(660) / '' /
DATA nitrogen(661) / '#SBL        !sublimation line specification' /
DATA nitrogen(662) / 'SB3  sublimation line model of Lemmon (1999).' /
DATA nitrogen(663) / '?LITERATURE REFERENCE \' /
DATA nitrogen(664) / '?Lemmon, E.W., 1999.' /
DATA nitrogen(665) / '?\' /
DATA nitrogen(666) / '!end of info section' /
DATA nitrogen(667) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(668) / '63.151             !upper temperature limit [K]' /
DATA nitrogen(669) / '0.0                !(dummy) upper pressure limit' /
DATA nitrogen(670) / '0.0                !(dummy) maximum density' /
DATA nitrogen(671) / '63.151 12.523      !reducing temperature and pressure' /
DATA nitrogen(672) / '0 1 0 0 0 0                 !number of terms in sublimation line equation' /
DATA nitrogen(673) / '-13.088692      1.          !coefficients and exponents' /
DATA nitrogen(674) / '' /
DATA nitrogen(675) / '' /
DATA nitrogen(676) / '#PS         !vapor pressure equation' /
DATA nitrogen(677) / 'PS6  vapor pressure equation of Span et al. (2000).' /
DATA nitrogen(678) / '?LITERATURE REFERENCE \' /
DATA nitrogen(679) / '?See EOS' /
DATA nitrogen(680) / '?\' /
DATA nitrogen(681) / '!end of info section' /
DATA nitrogen(682) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(683) / '126.192            !upper temperature limit [K]' /
DATA nitrogen(684) / '0.0                !(dummy) upper pressure limit' /
DATA nitrogen(685) / '0.0                !(dummy) maximum density' /
DATA nitrogen(686) / '126.192 3395.8     !reducing parameters' /
DATA nitrogen(687) / '4 0 0 0 0 0                 !number of terms in equation' /
DATA nitrogen(688) / ' -0.612445284E+01    2.     !coefficients and exponents' /
DATA nitrogen(689) / '  0.126327220E+01    3.' /
DATA nitrogen(690) / ' -0.765910082E+00    5.' /
DATA nitrogen(691) / ' -0.177570564E+01    10.' /
DATA nitrogen(692) / '' /
DATA nitrogen(693) / '' /
DATA nitrogen(694) / '#DL         !saturated liquid density equation' /
DATA nitrogen(695) / 'DL4  saturated liquid density equation of Span et al. (2000).' /
DATA nitrogen(696) / '?LITERATURE REFERENCE \' /
DATA nitrogen(697) / '?See EOS' /
DATA nitrogen(698) / '?\' /
DATA nitrogen(699) / '!end of info section' /
DATA nitrogen(700) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(701) / '126.192            !upper temperature limit [K]' /
DATA nitrogen(702) / '0.0                !(dummy) upper pressure limit' /
DATA nitrogen(703) / '0.0                !(dummy) maximum density' /
DATA nitrogen(704) / '126.192 11.1839    !reducing parameters' /
DATA nitrogen(705) / '4 0 0 0 0 0                 !number of terms in equation' /
DATA nitrogen(706) / '  0.148654237E+1     0.9882 !coefficients and exponents' /
DATA nitrogen(707) / ' -0.280476066E+0     2.' /
DATA nitrogen(708) / '  0.894143085E-1     8.' /
DATA nitrogen(709) / ' -0.119879866E+0     17.5' /
DATA nitrogen(710) / '' /
DATA nitrogen(711) / '' /
DATA nitrogen(712) / '#DV         !saturated vapor density equation' /
DATA nitrogen(713) / 'DV6  saturated vapor density equation of Span et al. (2000).' /
DATA nitrogen(714) / '?LITERATURE REFERENCE \' /
DATA nitrogen(715) / '?See EOS' /
DATA nitrogen(716) / '?\' /
DATA nitrogen(717) / '!end of info section' /
DATA nitrogen(718) / '63.151             !lower temperature limit [K]' /
DATA nitrogen(719) / '126.192            !upper temperature limit [K]' /
DATA nitrogen(720) / '0.0                !(dummy) upper pressure limit' /
DATA nitrogen(721) / '0.0                !(dummy) maximum density' /
DATA nitrogen(722) / '126.192 11.1839    !reducing parameters' /
DATA nitrogen(723) / '5 0 0 0 0 0                 !number of terms in equation' /
DATA nitrogen(724) / ' -0.170127164E+1     1.02   !coefficients and exponents' /
DATA nitrogen(725) / ' -0.370402649E+1     2.5' /
DATA nitrogen(726) / '  0.129859383E+1     3.5' /
DATA nitrogen(727) / ' -0.561424977E+0     6.5' /
DATA nitrogen(728) / ' -0.268505381E+1     14.0' /
DATA nitrogen(729) / '' /
DATA nitrogen(730) / '' /
DATA nitrogen(731) / '@END' /
DATA nitrogen(732) / 'c        1         2         3         4         5         6         7         8' /
DATA nitrogen(733) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
! #######################################################
character(256), TARGET :: oxygen(649)
DATA oxygen(1) / 'oxygen             !short name' /
DATA oxygen(2) / '7782-44-7          !CAS number' /
DATA oxygen(3) / 'oxygen             !full name' /
DATA oxygen(4) / 'O2                 !chemical formula' /
DATA oxygen(5) / 'R-732              !synonym' /
DATA oxygen(6) / '31.9988            !molecular weight [g/mol]' /
DATA oxygen(7) / '54.361             !triple point temperature [K]' /
DATA oxygen(8) / '90.1878            !normal boiling point [K]' /
DATA oxygen(9) / '154.581            !critical temperature [K]' /
DATA oxygen(10) / '5043.0             !critical pressure [kPa]' /
DATA oxygen(11) / '13.63              !critical density [mol/L]' /
DATA oxygen(12) / '0.0222             !acentric factor' /
DATA oxygen(13) / '0.0                !dipole moment [Debye]' /
DATA oxygen(14) / 'OT0                !default reference state' /
DATA oxygen(15) / '298.15  101.325  8680.0  205.043  !tref, Pref, Href, Sref' /
DATA oxygen(16) / '8.0                !version number' /
DATA oxygen(17) / '1072, 1073         !UN Number' /
DATA oxygen(18) / 'other              !family' /
DATA oxygen(19) / '0.0                !heating value (gross or superior) [kJ/mol]' /
DATA oxygen(20) / '' /
DATA oxygen(21) / '' /
DATA oxygen(22) / '#EOS               !equation of state specification' /
DATA oxygen(23) / 'FEQ  Helmholtz equation of state for oxygen of Schmidt and Wagner (1985).' /
DATA oxygen(24) / '?LITERATURE REFERENCE \' /
DATA oxygen(25) / '?Schmidt, R. and Wagner, W.,' /
DATA oxygen(26) / '? "A New Form of the Equation of State for Pure Substances and its' /
DATA oxygen(27) / '? Application to Oxygen,"' /
DATA oxygen(28) / '? Fluid Phase Equilibria, 19:175-200, 1985.' /
DATA oxygen(29) / '?\' /
DATA oxygen(30) / '?also published in:' /
DATA oxygen(31) / '?\' /
DATA oxygen(32) / '?Stewart, R.B., Jacobsen, R.T, and Wagner, W.,' /
DATA oxygen(33) / '? "Thermodynamic Properties of Oxygen from the Triple Point to 300 K' /
DATA oxygen(34) / '? with Pressures to 80 MPa,"' /
DATA oxygen(35) / '? J. Phys. Chem. Ref. Data, 20(5):917-1021, 1991.' /
DATA oxygen(36) / '?\' /
DATA oxygen(37) / '?The uncertainties of the equation of state are 0.1% in density, 2% in heat' /
DATA oxygen(38) / '? capacity, and 1% in the speed of sound, except in the critical region.' /
DATA oxygen(39) / '?\' /
DATA oxygen(40) / '!end of info section' /
DATA oxygen(41) / '54.361             !lower temperature limit [K]' /
DATA oxygen(42) / '2000.0             !upper temperature limit [K]' /
DATA oxygen(43) / '82000.0            !upper pressure limit [kPa]' /
DATA oxygen(44) / '43.348             !maximum density [mol/L]' /
DATA oxygen(45) / 'CPP                                    !pointer to Cp0 model' /
DATA oxygen(46) / '31.9988                                !molecular weight [g/mol]' /
DATA oxygen(47) / '54.361                                 !triple point temperature [K]' /
DATA oxygen(48) / '0.14628                                !pressure at triple point [kPa]' /
DATA oxygen(49) / '40.816                                 !density at triple point [mol/L]' /
DATA oxygen(50) / '90.1878                                !normal boiling point temperature [K]' /
DATA oxygen(51) / '0.0222                                 !acentric factor' /
DATA oxygen(52) / '154.581      5043.0       13.63        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA oxygen(53) / '154.581                   13.63        !reducing parameters [K, mol/L]' /
DATA oxygen(54) / '8.31434                                !gas constant [J/mol-K]' /
DATA oxygen(55) / '      32  4      0  0       0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA oxygen(56) / ' 0.39837687490d+00  0.000   1.00    0  !a(i),t(i),d(i),l(i)' /
DATA oxygen(57) / '-0.18461574540d+01  1.500   1.00    0' /
DATA oxygen(58) / ' 0.41834731970d+00  2.500   1.00    0' /
DATA oxygen(59) / ' 0.23706207110d-01 -0.500   2.00    0' /
DATA oxygen(60) / ' 0.97717305730d-01  1.500   2.00    0' /
DATA oxygen(61) / ' 0.30178912940d-01  2.000   2.00    0' /
DATA oxygen(62) / ' 0.22733532120d-01  0.000   3.00    0' /
DATA oxygen(63) / ' 0.13572540860d-01  1.000   3.00    0' /
DATA oxygen(64) / '-0.40526989430d-01  2.500   3.00    0' /
DATA oxygen(65) / ' 0.54546285150d-03  0.000   6.00    0' /
DATA oxygen(66) / ' 0.51131822770d-03  2.000   7.00    0' /
DATA oxygen(67) / ' 0.29534668830d-06  5.000   7.00    0' /
DATA oxygen(68) / '-0.86876450720d-04  2.000   8.00    0' /
DATA oxygen(69) / '-0.21270825890d+00  5.000   1.00    2' /
DATA oxygen(70) / ' 0.87359419580d-01  6.000   1.00    2' /
DATA oxygen(71) / ' 0.12755091900d+00  3.500   2.00    2' /
DATA oxygen(72) / '-0.90677010640d-01  5.500   2.00    2' /
DATA oxygen(73) / '-0.35400842060d-01  3.000   3.00    2' /
DATA oxygen(74) / '-0.36232780590d-01  7.000   3.00    2' /
DATA oxygen(75) / ' 0.13276992900d-01  6.000   5.00    2' /
DATA oxygen(76) / '-0.32541118650d-03  8.500   6.00    2' /
DATA oxygen(77) / '-0.83135829320d-02  4.000   7.00    2' /
DATA oxygen(78) / ' 0.21245705590d-02  6.500   8.00    2' /
DATA oxygen(79) / '-0.83252062320d-03  5.500  10.00    2' /
DATA oxygen(80) / '-0.26261732760d-04 22.000   2.00    4' /
DATA oxygen(81) / ' 0.25995814820d-02 11.000   3.00    4' /
DATA oxygen(82) / ' 0.99846496630d-02 18.000   3.00    4' /
DATA oxygen(83) / ' 0.21999231530d-02 11.000   4.00    4' /
DATA oxygen(84) / '-0.25913504860d-01 23.000   4.00    4' /
DATA oxygen(85) / '-0.12596308480d+00 17.000   5.00    4' /
DATA oxygen(86) / ' 0.14783556370d+00 18.000   5.00    4' /
DATA oxygen(87) / '-0.10112510780d-01 23.000   5.00    4' /
DATA oxygen(88) / '' /
DATA oxygen(89) / '' /
DATA oxygen(90) / '#AUX               !auxiliary model specification' /
DATA oxygen(91) / 'CPP  ideal gas heat capacity function' /
DATA oxygen(92) / '?LITERATURE REFERENCE \' /
DATA oxygen(93) / '?Refit by Roland Span of the Schmidt and Wagner equation listed below' /
DATA oxygen(94) / '?to account for the electronic contribution up to 2000 K by using' /
DATA oxygen(95) / '?Planck-Einstein terms only.' /
DATA oxygen(96) / '?\' /
DATA oxygen(97) / '?Schmidt, R. and Wagner, W.,' /
DATA oxygen(98) / '? "A New Form of the Equation of State for Pure Substances and its' /
DATA oxygen(99) / '? Application to Oxygen,"' /
DATA oxygen(100) / '? Fluid Phase Equilibria, 19:175-200, 1985.' /
DATA oxygen(101) / '?\' /
DATA oxygen(102) / '!end of info section' /
DATA oxygen(103) / '54.361             !lower temperature limit [K]' /
DATA oxygen(104) / '2000.0             !upper temperature limit [K]' /
DATA oxygen(105) / '0.0                !upper pressure limit [kPa]' /
DATA oxygen(106) / '0.0                !maximum density [mol/L]' /
DATA oxygen(107) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA oxygen(108) / '  1  5    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA oxygen(109) / '3.51808732          0.0' /
DATA oxygen(110) / '0.102323928D+01  0.224632440D+04' /
DATA oxygen(111) / '0.784357918D+00  0.112599763D+05' /
DATA oxygen(112) / '0.337183363D-02  0.120126209D+04' /
DATA oxygen(113) / '-.170864084D-01  0.690089445D+02' /
DATA oxygen(114) / '0.463751562D-01  0.532805445D+04' /
DATA oxygen(115) / '' /
DATA oxygen(116) / '' /
DATA oxygen(117) / '#AUX               !auxiliary model specification' /
DATA oxygen(118) / 'CPx  ideal gas heat capacity function' /
DATA oxygen(119) / '?LITERATURE REFERENCE \' /
DATA oxygen(120) / '?Schmidt, R. and Wagner, W.,' /
DATA oxygen(121) / '? "A New Form of the Equation of State for Pure Substances and its' /
DATA oxygen(122) / '? Application to Oxygen,"' /
DATA oxygen(123) / '? Fluid Phase Equilibria, 19:175-200, 1985.' /
DATA oxygen(124) / '?\' /
DATA oxygen(125) / '?The electronic part of the equation of Schmidt and Wagner is not included here.' /
DATA oxygen(126) / '?\' /
DATA oxygen(127) / '!end of info section' /
DATA oxygen(128) / '54.361             !lower temperature limit [K]' /
DATA oxygen(129) / '300.0              !upper temperature limit [K]' /
DATA oxygen(130) / '0.0                !upper pressure limit [kPa]' /
DATA oxygen(131) / '0.0                !maximum density [mol/L]' /
DATA oxygen(132) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA oxygen(133) / '  3  1    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA oxygen(134) / '0.10677800d+01     -1.50' /
DATA oxygen(135) / '0.35004200d+01      0.00' /
DATA oxygen(136) / '0.16696100d-07      2.00' /
DATA oxygen(137) / '0.10125800d+01   2242.45' /
DATA oxygen(138) / '' /
DATA oxygen(139) / '' /
DATA oxygen(140) / '@EOS               !equation of state specification' /
DATA oxygen(141) / 'FEK  Helmholtz equation of state for oxygen of Kunz and Wagner (2004).' /
DATA oxygen(142) / '?LITERATURE REFERENCE \' /
DATA oxygen(143) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA oxygen(144) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA oxygen(145) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA oxygen(146) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA oxygen(147) / '?\' /
DATA oxygen(148) / '!end of info section' /
DATA oxygen(149) / '54.361             !lower temperature limit [K]' /
DATA oxygen(150) / '1000.0             !upper temperature limit [K]' /
DATA oxygen(151) / '82000.0            !upper pressure limit [kPa]' /
DATA oxygen(152) / '43.348             !maximum density [mol/L]' /
DATA oxygen(153) / 'PHK                                    !pointer to Cp0 model' /
DATA oxygen(154) / '31.9988                                !molecular weight [g/mol]' /
DATA oxygen(155) / '54.361                                 !triple point temperature [K]' /
DATA oxygen(156) / '0.1460                                 !pressure at triple point [kPa]' /
DATA oxygen(157) / '40.89                                  !density at triple point [mol/L]' /
DATA oxygen(158) / ' 90.18                                 !normal boiling point temperature [K]' /
DATA oxygen(159) / ' 0.0236                                !acentric factor' /
DATA oxygen(160) / '154.595      5061.6      13.63         !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA oxygen(161) / '154.595                  13.63         !reducing parameters [K, mol/L]' /
DATA oxygen(162) / '8.314472                               !gas constant [J/mol-K]' /
DATA oxygen(163) / '  12  4      0  0      0  0            !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA oxygen(164) / ' 0.88878286369701       0.250  1.  0' /
DATA oxygen(165) / '-0.24879433312148d1     1.125  1.  0' /
DATA oxygen(166) / ' 0.59750190775886       1.500  1.  0' /
DATA oxygen(167) / ' 0.96501817061881d-2    1.375  2.  0' /
DATA oxygen(168) / ' 0.71970428712770d-1    0.250  3.  0' /
DATA oxygen(169) / ' 0.22337443000195d-3    0.875  7.  0' /
DATA oxygen(170) / ' 0.18558686391474       0.625  2.  1' /
DATA oxygen(171) / '-0.38129368035760d-1    1.750  5.  1' /
DATA oxygen(172) / '-0.15352245383006       3.625  1.  2' /
DATA oxygen(173) / '-0.26726814910919d-1    3.625  4.  2' /
DATA oxygen(174) / '-0.25675298677127d-1    14.5   3.  3' /
DATA oxygen(175) / ' 0.95714302123668d-2    12.0   4.  3' /
DATA oxygen(176) / '' /
DATA oxygen(177) / '' /
DATA oxygen(178) / '#AUX               !auxiliary model specification' /
DATA oxygen(179) / 'PHK  Helmholtz form for the ideal-gas state for oxygen of Kunz and Wagner (2004).' /
DATA oxygen(180) / '?LITERATURE REFERENCE \' /
DATA oxygen(181) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA oxygen(182) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA oxygen(183) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA oxygen(184) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA oxygen(185) / '?\' /
DATA oxygen(186) / '!end of info section' /
DATA oxygen(187) / '0.                 !lower temperature limit [K]' /
DATA oxygen(188) / '1000.0             !upper temperature limit [K]' /
DATA oxygen(189) / '0.0                !upper pressure limit [kPa]' /
DATA oxygen(190) / '0.0                !maximum density [mol/L]' /
DATA oxygen(191) / '1 2  0  1 1  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA oxygen(192) / '    2.50146      1.             !ai, ti for [ai*log(tau**ti)] terms' /
DATA oxygen(193) / '   10.001843586  0.             !aj, ti for [ai*tau**ti] terms' /
DATA oxygen(194) / '  -14.996095135  1.' /
DATA oxygen(195) / '   -1.01334      7.223325463    !aj, ti for cosh and sinh terms' /
DATA oxygen(196) / '    1.07558     14.461722565' /
DATA oxygen(197) / '' /
DATA oxygen(198) / '' /
DATA oxygen(199) / '@EOS               !equation of state specification' /
DATA oxygen(200) / 'FES  short Helmholtz equation of state for oxygen of Span and Wagner (2003).' /
DATA oxygen(201) / '?LITERATURE REFERENCE \' /
DATA oxygen(202) / '?Span, R. and Wagner, W.' /
DATA oxygen(203) / '? "Equations of State for Technical Applications. II. Results for Nonpolar Fluids,"' /
DATA oxygen(204) / '? Int. J. Thermophys., 24(1):41-109, 2003.' /
DATA oxygen(205) / '?\' /
DATA oxygen(206) / '?The uncertainties of the equation of state are approximately 0.2% (to' /
DATA oxygen(207) / '?0.5% at high pressures) in density, 1% (in the vapor phase) to 2% in' /
DATA oxygen(208) / '?heat capacity, 1% (in the vapor phase) to 2% in the speed of sound, and' /
DATA oxygen(209) / '?0.2% in vapor pressure, except in the critical region.' /
DATA oxygen(210) / '?\' /
DATA oxygen(211) / '!end of info section' /
DATA oxygen(212) / '54.361             !lower temperature limit [K]' /
DATA oxygen(213) / '600.0              !upper temperature limit [K]' /
DATA oxygen(214) / '100000.0           !upper pressure limit [kPa]' /
DATA oxygen(215) / '43.348             !maximum density [mol/L]' /
DATA oxygen(216) / 'CPP                                    !pointer to Cp0 model' /
DATA oxygen(217) / '31.999                                 !molecular weight [g/mol]' /
DATA oxygen(218) / '54.361                                 !triple point temperature [K]' /
DATA oxygen(219) / '0.14603                                !pressure at triple point [kPa]' /
DATA oxygen(220) / '40.885                                 !density at triple point [mol/L]' /
DATA oxygen(221) / '90.182                                 !normal boiling point temperature [K]' /
DATA oxygen(222) / '0.0222                                 !acentric factor' /
DATA oxygen(223) / '154.595      5043.0       13.63        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA oxygen(224) / '154.595                   13.63        !reducing parameters [K, mol/L]' /
DATA oxygen(225) / '8.31451                                !gas constant [J/mol-K]' /
DATA oxygen(226) / '      12  4      0  0       0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA oxygen(227) / ' 0.888782860000E+00  0.25    1.0     0 !a(i),t(i),d(i),l(i)' /
DATA oxygen(228) / '-0.248794330000E+01  1.125   1.0     0' /
DATA oxygen(229) / ' 0.597501910000E+00  1.5     1.0     0' /
DATA oxygen(230) / ' 0.965018170000E-02  1.375   2.0     0' /
DATA oxygen(231) / ' 0.719704290000E-01  0.25    3.0     0' /
DATA oxygen(232) / ' 0.223374430000E-03  0.875   7.0     0' /
DATA oxygen(233) / ' 0.185586860000E+00  0.625   2.0     1' /
DATA oxygen(234) / '-0.381293680000E-01  1.75    5.0     1' /
DATA oxygen(235) / '-0.153522450000E+00  3.625   1.0     2' /
DATA oxygen(236) / '-0.267268150000E-01  3.625   4.0     2' /
DATA oxygen(237) / '-0.256752990000E-01 14.5     3.0     3' /
DATA oxygen(238) / ' 0.957143020000E-02 12.0     4.0     3' /
DATA oxygen(239) / '' /
DATA oxygen(240) / '' /
DATA oxygen(241) / '@EOS               !equation of state specification' /
DATA oxygen(242) / 'BWR  MBWR equation of state for oxygen of Younglove (1982).' /
DATA oxygen(243) / '?LITERATURE REFERENCE \' /
DATA oxygen(244) / '?Younglove, B.A.,' /
DATA oxygen(245) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA oxygen(246) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA oxygen(247) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA oxygen(248) / '?\' /
DATA oxygen(249) / '!end of info section' /
DATA oxygen(250) / '54.359             !lower temperature limit [K]' /
DATA oxygen(251) / '400.0              !upper temperature limit [K]' /
DATA oxygen(252) / '121000.0           !upper pressure limit [kPa]' /
DATA oxygen(253) / '40.820             !maximum density [mol/L]' /
DATA oxygen(254) / 'CP1                                    !pointer to Cp0 model' /
DATA oxygen(255) / '31.9988                                !molecular weight [g/mol]' /
DATA oxygen(256) / '54.359                                 !triple point temperature [K]' /
DATA oxygen(257) / '0.148                                  !pressure at triple point [kPa]' /
DATA oxygen(258) / '40.820                                 !density at triple point [mol/L]' /
DATA oxygen(259) / '90.1878                                !normal boiling point temperature [K]' /
DATA oxygen(260) / '0.0222                                 !acentric factor' /
DATA oxygen(261) / '154.581      5043.        13.63        !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA oxygen(262) / '154.581                   13.63        !reducing parameters [K, mol/L]' /
DATA oxygen(263) / '13.3630620956                          !gamma' /
DATA oxygen(264) / '0.0831411                              !gas constant [L-bar/mol-K]' /
DATA oxygen(265) / '      32       1                       !Nterm, Ncoeff per term' /
DATA oxygen(266) / ' -0.4365859650d-03      0.2005820677d-00      -0.4197909916d+01' /
DATA oxygen(267) / '  0.1878215317d+03     -0.1287473398d+05       0.1556745888d-04' /
DATA oxygen(268) / '  0.1343639359d-02     -0.2228415518d+01       0.4767792275d+04' /
DATA oxygen(269) / '  0.4790846641d-06      0.2462611107d-02      -0.1921891680d-00' /
DATA oxygen(270) / ' -0.6978320847d-05     -0.6214145909d-03      -0.1860852567d-00' /
DATA oxygen(271) / '  0.2609791417d-04     -0.2447611408d-06       0.1457743352d-03' /
DATA oxygen(272) / ' -0.1726492873d-05     -0.2384892520d+04      -0.2301807796d+06' /
DATA oxygen(273) / ' -0.2790303526d+02      0.9400577575d+05      -0.4169449637d-01' /
DATA oxygen(274) / '  0.2008497853d+01     -0.1256076520d-03      -0.6406362964d-00' /
DATA oxygen(275) / ' -0.2475580168d-07      0.1346309703d-04      -0.1161502470d-09' /
DATA oxygen(276) / ' -0.1034699798d-07      0.2365936964d-06' /
DATA oxygen(277) / '' /
DATA oxygen(278) / '' /
DATA oxygen(279) / '#AUX               !auxiliary model specification' /
DATA oxygen(280) / 'CP1  ideal gas heat capacity function of Younglove' /
DATA oxygen(281) / '?LITERATURE REFERENCE \' /
DATA oxygen(282) / '?Younglove, B.A.,' /
DATA oxygen(283) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA oxygen(284) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA oxygen(285) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA oxygen(286) / '?\' /
DATA oxygen(287) / '!end of info section' /
DATA oxygen(288) / '54.359             !lower temperature limit [K]' /
DATA oxygen(289) / '400.0              !upper temperature limit [K]' /
DATA oxygen(290) / '0.0                !upper pressure limit [kPa]' /
DATA oxygen(291) / '0.0                !maximum density [mol/L]' /
DATA oxygen(292) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA oxygen(293) / '  7  1    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA oxygen(294) / ' -0.4981998537119d+04    -3.00d0' /
DATA oxygen(295) / '  0.2302477799952d+03    -2.00d0' /
DATA oxygen(296) / ' -0.3455653235107d+01    -1.00d0' /
DATA oxygen(297) / '  0.3521876773671d+01     0.00d0' /
DATA oxygen(298) / ' -0.4354202160244d-04     1.00d0' /
DATA oxygen(299) / '  0.1346353450132d-07     2.00d0' /
DATA oxygen(300) / '  0.1620598259591d-10     3.00d0' /
DATA oxygen(301) / '  0.1031468515726d+01  2239.18105' /
DATA oxygen(302) / '' /
DATA oxygen(303) / '' /
DATA oxygen(304) / '' /
DATA oxygen(305) / '#TCX               !thermal conductivity model specification' /
DATA oxygen(306) / 'TC1  pure fluid thermal conductivity model of Lemmon and Jacobsen (2004).' /
DATA oxygen(307) / '?LITERATURE REFERENCE \' /
DATA oxygen(308) / '?Lemmon, E.W. and Jacobsen, R.T,' /
DATA oxygen(309) / '? "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon,' /
DATA oxygen(310) / '? and Air,"' /
DATA oxygen(311) / '? Int. J. Thermophys., 25:21-69, 2004.' /
DATA oxygen(312) / '?\' /
DATA oxygen(313) / '?The uncertainty for the dilute gas is 2% with increasing uncertainties' /
DATA oxygen(314) / '?near the triple point.  The uncertainties range from 3% between 270 and' /
DATA oxygen(315) / '?300 K to 5% elsewhere.  The uncertainties above 100 MPa are not known due' /
DATA oxygen(316) / '?to a lack of experimental data.' /
DATA oxygen(317) / '?\' /
DATA oxygen(318) / '!end of info section' /
DATA oxygen(319) / '54.361             !lower temperature limit [K]' /
DATA oxygen(320) / '2000.0             !upper temperature limit [K]' /
DATA oxygen(321) / '82000.0            !upper pressure limit [kPa]' /
DATA oxygen(322) / '43.348             !maximum density [mol/L]' /
DATA oxygen(323) / '3   0              !# terms for dilute gas function:  numerator, denominator' /
DATA oxygen(324) / '154.581   1.0d-3   !reducing parameters for T, tcx' /
DATA oxygen(325) / ' 1.036   -97.0     !coeff, power in T' /
DATA oxygen(326) / ' 6.283     0.9' /
DATA oxygen(327) / '-4.262     0.6' /
DATA oxygen(328) / '6   0              !# terms for background gas function:  numerator, denominator' /
DATA oxygen(329) / '154.581   13.63       1.0d-3    !reducing parameters for T, rho, tcx' /
DATA oxygen(330) / '15.31           0.0  1.0  0.0 !coeff, powers of T, rho, exp(rho)' /
DATA oxygen(331) / ' 8.898        0.0    3.   0.' /
DATA oxygen(332) / '-0.7336      -0.3    4.   0.' /
DATA oxygen(333) / ' 6.728       -4.3    5.   2.' /
DATA oxygen(334) / '-4.374       -0.5    7.   2.' /
DATA oxygen(335) / '-0.4747      -1.8   10.   2.' /
DATA oxygen(336) / 'TK3                !pointer to critical enhancement auxiliary function' /
DATA oxygen(337) / '' /
DATA oxygen(338) / '' /
DATA oxygen(339) / '#AUX               !thermal conductivity critical enhancement model' /
DATA oxygen(340) / 'TK3  thermal conductivity critical enhancement of Lemmon and Jacobsen (2004).' /
DATA oxygen(341) / '?LITERATURE REFERENCE \' /
DATA oxygen(342) / '?\' /
DATA oxygen(343) / '!end of info section' /
DATA oxygen(344) / '54.361             !lower temperature limit [K]' /
DATA oxygen(345) / '2000.0             !upper temperature limit [K]' /
DATA oxygen(346) / '82000.0            !upper pressure limit [kPa]' /
DATA oxygen(347) / '43.348             !maximum density [mol/L]' /
DATA oxygen(348) / '9  0  0  0         !# terms:  terms, spare, spare, spare' /
DATA oxygen(349) / '1.0    1.0  1.0    !reducing par for T, rho, tcx (mW/m-K)' /
DATA oxygen(350) / '0.630d0            !gnu (universal exponent)' /
DATA oxygen(351) / '1.2415d0           !gamma (universal exponent)' /
DATA oxygen(352) / '1.01d0             !R0 (universal amplitude)' /
DATA oxygen(353) / ' 0.065d0           !z (universal exponent--not used for t.c., only viscosity)' /
DATA oxygen(354) / ' 1.00d0            !c (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA oxygen(355) / ' 0.24E-09          !xi0 (amplitude) [m]' /
DATA oxygen(356) / ' 0.55E-01          !gam0 (amplitude) [-]' /
DATA oxygen(357) / ' 0.51E-09          !qd_inverse (modified effective cutoff parameter) [m]' /
DATA oxygen(358) / '309.162            !tref (reference temperature) [K]' /
DATA oxygen(359) / '' /
DATA oxygen(360) / '' /
DATA oxygen(361) / '#ETA               !viscosity model specification' /
DATA oxygen(362) / 'VS1  pure fluid viscosity model of Lemmon and Jacobsen (2004).' /
DATA oxygen(363) / '?LITERATURE REFERENCE \' /
DATA oxygen(364) / '?Lemmon, E.W. and Jacobsen, R.T,' /
DATA oxygen(365) / '? "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon,' /
DATA oxygen(366) / '? and Air,"' /
DATA oxygen(367) / '? Int. J. Thermophys., 25:21-69, 2004.' /
DATA oxygen(368) / '?\' /
DATA oxygen(369) / '?The uncertainty is 1% in the dilute gas at temperatures above 200 K, and' /
DATA oxygen(370) / '?5% in the dilute gas at lower temperatures.  The uncertainty is around' /
DATA oxygen(371) / '?2% between 270 and 300 K, and increases to 5% outside of this region.' /
DATA oxygen(372) / '?The uncertainty may be higher in the liquid near the triple point.' /
DATA oxygen(373) / '?\' /
DATA oxygen(374) / '!end of info section' /
DATA oxygen(375) / '54.361             !lower temperature limit [K]' /
DATA oxygen(376) / '2000.0             !upper temperature limit [K]' /
DATA oxygen(377) / '82000.0            !upper pressure limit [kPa]' /
DATA oxygen(378) / '43.348             !maximum density [mol/L]' /
DATA oxygen(379) / '1                  !number of terms associated with dilute-gas function' /
DATA oxygen(380) / 'CI1                !pointer to reduced effective collision cross-section model' /
DATA oxygen(381) / '0.3428             !Lennard-Jones coefficient sigma [nm]' /
DATA oxygen(382) / '118.5              !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA oxygen(383) / '1.0    1.0         !reducing parameters for T, eta' /
DATA oxygen(384) / '0.151011418  0.5   !Chapman-Enskog term' /
DATA oxygen(385) / '0                  !number of terms for initial density dependence' /
DATA oxygen(386) / '0 5 0 0 0 0        !# resid terms:  close-packed density;  simple poly; numerator of rational poly; denominator of rat. poly; numerator of exponential; denominator of exponential' /
DATA oxygen(387) / '154.581   13.63       1.0           !reducing parameters for T, rho, eta' /
DATA oxygen(388) / ' 17.67       -0.05   1.   0.   0    !simple polynomial terms' /
DATA oxygen(389) / ' 0.4042       0.0    5.   0.   0' /
DATA oxygen(390) / ' 0.0001077   -2.10  12.   0.   0' /
DATA oxygen(391) / ' 0.3510       0.0    8.   0.   1' /
DATA oxygen(392) / '-13.67       -0.5    1.   0.   2' /
DATA oxygen(393) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA oxygen(394) / '' /
DATA oxygen(395) / '' /
DATA oxygen(396) / '#AUX               !collision integral specification' /
DATA oxygen(397) / 'CI1  collision integral model of Lemmon and Jacobsen (2004).' /
DATA oxygen(398) / '?LITERATURE REFERENCE \' /
DATA oxygen(399) / '?\' /
DATA oxygen(400) / '!end of info section' /
DATA oxygen(401) / '1.0                !lower temperature limit [K]' /
DATA oxygen(402) / '10000.0            !upper temperature limit [K]' /
DATA oxygen(403) / '0.0                !(dummy) upper pressure limit' /
DATA oxygen(404) / '0.0                !(dummy) maximum density' /
DATA oxygen(405) / '5                  !number of terms' /
DATA oxygen(406) / '  0.431      0     !coeff, power of Tstar' /
DATA oxygen(407) / ' -0.4623     1' /
DATA oxygen(408) / '  0.08406    2' /
DATA oxygen(409) / '  0.005341   3' /
DATA oxygen(410) / ' -0.00331    4' /
DATA oxygen(411) / '' /
DATA oxygen(412) / '' /
DATA oxygen(413) / '@TCX               !thermal conductivity model specification' /
DATA oxygen(414) / 'TC3  pure fluid thermal conductivity model of Younglove (1982).' /
DATA oxygen(415) / '?LITERATURE REFERENCE \' /
DATA oxygen(416) / '?Younglove, B.A.,' /
DATA oxygen(417) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA oxygen(418) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA oxygen(419) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA oxygen(420) / '?\' /
DATA oxygen(421) / '!end of info section' /
DATA oxygen(422) / '54.361             !lower temperature limit [K]' /
DATA oxygen(423) / '600.0              !upper temperature limit [K]' /
DATA oxygen(424) / '100000.0           !upper pressure limit [kPa]' /
DATA oxygen(425) / '43.348             !maximum density [mol/L]' /
DATA oxygen(426) / '0.3437             !Lennard-Jones coefficient sigma [nm]' /
DATA oxygen(427) / '113                !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA oxygen(428) / '0.15099557923496   !const in Eq 20 = 5/16*(k*MW/1000/pi/Na)**0.5*1.0d12' /
DATA oxygen(429) / ' 0                 !exponent in Eq 20 for T' /
DATA oxygen(430) / '-1.41202117453516  !eta0 (eq 20): coeffs of {(e/kT)**((4-n)/3)}' /
DATA oxygen(431) / ' 8.06267523869911' /
DATA oxygen(432) / '-19.44147946395' /
DATA oxygen(433) / ' 25.78193316324' /
DATA oxygen(434) / '-20.5167203343277' /
DATA oxygen(435) / ' 10.0087040966906' /
DATA oxygen(436) / '-2.90450673487991' /
DATA oxygen(437) / '0.459605807669332' /
DATA oxygen(438) / '-3.01906029521D-02' /
DATA oxygen(439) / '0.00097916328      !coeff for initial density dependence of viscosity (eq 21); Fv(1)' /
DATA oxygen(440) / '0.00089116658      !Fv(2)' /
DATA oxygen(441) / ' 1.12              !Fv(3)' /
DATA oxygen(442) / ' 100               !Fv(4)' /
DATA oxygen(443) / '-21.520741137      !coefficients for residual viscosity, eqs (22 - 25)' /
DATA oxygen(444) / ' 473.50508788      !Ev(2)' /
DATA oxygen(445) / ' 11.9072051301147  !Ev(3)' /
DATA oxygen(446) / '-2122.44247203833  !Ev(4)' /
DATA oxygen(447) / ' 0                 !Ev(5)' /
DATA oxygen(448) / ' 0                 !Ev(6)' /
DATA oxygen(449) / ' 0                 !Ev(7)' /
DATA oxygen(450) / ' 31.251171918947   !Ev(8)' /
DATA oxygen(451) / ' 2.21064           !F' /
DATA oxygen(452) / '0.000000038896     !rm' /
DATA oxygen(453) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA oxygen(454) / '' /
DATA oxygen(455) / '' /
DATA oxygen(456) / '@ETA               !viscosity model specification' /
DATA oxygen(457) / 'VS2  pure fluid viscosity model of Younglove (1982).' /
DATA oxygen(458) / '?LITERATURE REFERENCE \' /
DATA oxygen(459) / '?Younglove, B.A.,' /
DATA oxygen(460) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA oxygen(461) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA oxygen(462) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA oxygen(463) / '?\' /
DATA oxygen(464) / '!end of info section' /
DATA oxygen(465) / '54.361             !lower temperature limit [K]' /
DATA oxygen(466) / '600.0              !upper temperature limit [K]' /
DATA oxygen(467) / '100000.0           !upper pressure limit [kPa]' /
DATA oxygen(468) / '43.348             !maximum density [mol/L]' /
DATA oxygen(469) / 'CI2                !pointer to collision integral model' /
DATA oxygen(470) / '0.3437             !Lennard-Jones coefficient sigma [nm]' /
DATA oxygen(471) / '113                !Lennard-Jones coefficient epsilon/kappa [K]' /
DATA oxygen(472) / '0.15099557923496   !const in Eq 19 = 5/16*(k*MW/1000/pi/Na)**0.5*1.0d12' /
DATA oxygen(473) / ' 0                 !exponent in Eq 20 for T' /
DATA oxygen(474) / ' 1.39279625307D-02 !coeff for initial density dependence of viscosity (eq 21); Fv(1)' /
DATA oxygen(475) / '-6.51536010579D-03 !Fv(2)' /
DATA oxygen(476) / ' 1.4               !Fv(3)' /
DATA oxygen(477) / ' 100               !Fv(4)' /
DATA oxygen(478) / '-14.45497211       !coefficients for residual viscosity, eqs (22 - 25)' /
DATA oxygen(479) / ' 243.40689667      !Ev(2)' /
DATA oxygen(480) / ' 12.9006761056004  !Ev(3)' /
DATA oxygen(481) / '-1949.07966423848  !Ev(4)' /
DATA oxygen(482) / '-5.62078436742D-02 !Ev(5)' /
DATA oxygen(483) / ' 21.3075467849104  !Ev(6)' /
DATA oxygen(484) / ' 48.9965711691056  !Ev(7)' /
DATA oxygen(485) / ' 13.5942597847419  !Ev(8)' /
DATA oxygen(486) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA oxygen(487) / '' /
DATA oxygen(488) / '' /
DATA oxygen(489) / '@AUX               !collision integral specification' /
DATA oxygen(490) / 'CI2  collision integral model of Younglove (1982).' /
DATA oxygen(491) / '?LITERATURE REFERENCE \' /
DATA oxygen(492) / '?Younglove, B.A.,' /
DATA oxygen(493) / '? "Thermophysical Properties of Fluids.  I. Argon, Ethylene,' /
DATA oxygen(494) / '? Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen,"' /
DATA oxygen(495) / '? J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.' /
DATA oxygen(496) / '?\' /
DATA oxygen(497) / '!end of info section' /
DATA oxygen(498) / '54.361             !lower temperature limit [K]' /
DATA oxygen(499) / '625.0              !upper temperature limit [K]' /
DATA oxygen(500) / '0.0                !(dummy) upper pressure limit' /
DATA oxygen(501) / '0.0                !(dummy) maximum density' /
DATA oxygen(502) / '9                  !number of terms' /
DATA oxygen(503) / '-67.2093902106092  !eta0 (eq 20): coeffs of {(e/kT)**((4-n)/3)}' /
DATA oxygen(504) / ' 277.148660965491' /
DATA oxygen(505) / '-399.192753863192' /
DATA oxygen(506) / ' 166.828729537446' /
DATA oxygen(507) / ' 143.163477478684' /
DATA oxygen(508) / '-191.767060368781' /
DATA oxygen(509) / ' 98.4332230147836' /
DATA oxygen(510) / '-22.9410694301649' /
DATA oxygen(511) / ' 2.12402264924749' /
DATA oxygen(512) / '' /
DATA oxygen(513) / '' /
DATA oxygen(514) / '#STN        !surface tension specification' /
DATA oxygen(515) / 'ST1  surface tension model of Mulero et al. (2012)' /
DATA oxygen(516) / '?LITERATURE REFERENCE \' /
DATA oxygen(517) / '?Mulero, A., Cachadina, I., and Parra, M.I.' /
DATA oxygen(518) / '? "Recommended Correlations for the Surface Tension of Common Fluids,"' /
DATA oxygen(519) / '? J. Phys. Chem. Ref. Data, 41, 043105, 2012.' /
DATA oxygen(520) / '?\' /
DATA oxygen(521) / '!end of info section' /
DATA oxygen(522) / '0.0                !lower temperature limit [K]' /
DATA oxygen(523) / '154.581            !upper temperature limit [K]' /
DATA oxygen(524) / '0.0                !(dummy) upper pressure limit' /
DATA oxygen(525) / '0.0                !(dummy) maximum density' /
DATA oxygen(526) / '1                           !number of terms in surface tension model' /
DATA oxygen(527) / '154.581                     !critical temperature used in fit (dummy)' /
DATA oxygen(528) / ' 0.03843     1.225          !sigma0 and n' /
DATA oxygen(529) / '' /
DATA oxygen(530) / '' /
DATA oxygen(531) / '#DE         !dielectric constant specification' /
DATA oxygen(532) / 'DE3  dielectric constant model of Harvey and Lemmon (2005).' /
DATA oxygen(533) / '?LITERATURE REFERENCE \' /
DATA oxygen(534) / '?Harvey, A.H. and Lemmon, E.W.' /
DATA oxygen(535) / '? "Method for Estimating the Dielectric Constant of Natural Gas Mixtures,"' /
DATA oxygen(536) / '? Int. J. Thermophys., 26(1):31-46, 2005.' /
DATA oxygen(537) / '?\' /
DATA oxygen(538) / '!end of info section' /
DATA oxygen(539) / '0.0                !lower temperature limit [K]' /
DATA oxygen(540) / '2000.0             !upper temperature limit [K]' /
DATA oxygen(541) / '0.0                !(dummy) upper pressure limit' /
DATA oxygen(542) / '0.0                !(dummy) maximum density' /
DATA oxygen(543) / '273.16 1000.0 1.0  !reducing parameters for t and d' /
DATA oxygen(544) / '0 2 4 0 0 0                         !number of terms in dielectric constant model' /
DATA oxygen(545) / ' 3.9578           0.    1.    0.    !coef, t exp, d exp' /
DATA oxygen(546) / ' 0.0065           1.    1.    0.' /
DATA oxygen(547) / ' 0.575            0.    2.    0.' /
DATA oxygen(548) / ' 1.028            1.    2.    0.' /
DATA oxygen(549) / '-8.96             0.    2.5   0.' /
DATA oxygen(550) / '-5.15             1.    2.5   0.' /
DATA oxygen(551) / '' /
DATA oxygen(552) / '' /
DATA oxygen(553) / '#MLT        !melting line specification' /
DATA oxygen(554) / 'ML2  melting line model of Schmidt and Wagner (1985).' /
DATA oxygen(555) / '?LITERATURE REFERENCE \' /
DATA oxygen(556) / '?Schmidt, R. and Wagner, W.,' /
DATA oxygen(557) / '? "A New Form of the Equation of State for Pure Substances and its' /
DATA oxygen(558) / '? Application to Oxygen,"' /
DATA oxygen(559) / '? Fluid Phase Equilibria, 19:175-200, 1985.' /
DATA oxygen(560) / '?\' /
DATA oxygen(561) / '!end of info section' /
DATA oxygen(562) / '54.361             !lower temperature limit [K]' /
DATA oxygen(563) / '300.0              !upper temperature limit [K]' /
DATA oxygen(564) / '0.0                !(dummy) upper pressure limit' /
DATA oxygen(565) / '0.0                !(dummy) maximum density' /
DATA oxygen(566) / '54.361 0.14633     !reducing temperature and pressure' /
DATA oxygen(567) / '0 4 0 0 0 0                 !number of terms in melting line equation' /
DATA oxygen(568) / '-0.32463539d+2  0.0625      !coefficients and exponents' /
DATA oxygen(569) / ' 0.14278011d+3  0.1250' /
DATA oxygen(570) / '-0.14702341d+3  0.1875' /
DATA oxygen(571) / ' 0.52001200d+2  0.2500' /
DATA oxygen(572) / '' /
DATA oxygen(573) / '' /
DATA oxygen(574) / '#SBL        !sublimation line specification' /
DATA oxygen(575) / 'SB3  sublimation line model of Lemmon (2003).' /
DATA oxygen(576) / '?LITERATURE REFERENCE \' /
DATA oxygen(577) / '?Lemmon, E.W., 2003.' /
DATA oxygen(578) / '?\' /
DATA oxygen(579) / '!end of info section' /
DATA oxygen(580) / '54.361             !lower temperature limit [K]' /
DATA oxygen(581) / '54.361             !upper temperature limit [K]' /
DATA oxygen(582) / '0.0                !(dummy) upper pressure limit' /
DATA oxygen(583) / '0.0                !(dummy) maximum density' /
DATA oxygen(584) / '54.361 0.14628     !reducing temperature and pressure' /
DATA oxygen(585) / '0 1 0 0 0 0                 !number of terms in sublimation line equation' /
DATA oxygen(586) / '-20.714   1.06              !coefficients and exponents' /
DATA oxygen(587) / '' /
DATA oxygen(588) / '' /
DATA oxygen(589) / '#PS         !vapor pressure equation' /
DATA oxygen(590) / 'PS5  vapor pressure equation of Cullimore (2010).' /
DATA oxygen(591) / '?LITERATURE REFERENCE \' /
DATA oxygen(592) / '?Cullimore, I.D., 2010.' /
DATA oxygen(593) / '?\' /
DATA oxygen(594) / '!end of info section' /
DATA oxygen(595) / '54.361             !lower temperature limit [K]' /
DATA oxygen(596) / '154.581            !upper temperature limit [K]' /
DATA oxygen(597) / '0.0                !(dummy) upper pressure limit' /
DATA oxygen(598) / '0.0                !(dummy) maximum density' /
DATA oxygen(599) / '154.581 5043.0     !reducing parameters' /
DATA oxygen(600) / '5 0 0 0 0 0        !number of terms in equation' /
DATA oxygen(601) / '-0.60595D+01   1.0' /
DATA oxygen(602) / ' 0.13050D+01   1.5' /
DATA oxygen(603) / '-0.54178D+00   2.2' /
DATA oxygen(604) / '-0.19410D+01   4.8' /
DATA oxygen(605) / ' 0.35514D+00   6.2' /
DATA oxygen(606) / '' /
DATA oxygen(607) / '' /
DATA oxygen(608) / '#DL         !saturated liquid density equation' /
DATA oxygen(609) / 'DL1  saturated liquid density equation of Cullimore (2010).' /
DATA oxygen(610) / '?LITERATURE REFERENCE \' /
DATA oxygen(611) / '?Cullimore, I.D., 2010.' /
DATA oxygen(612) / '?\' /
DATA oxygen(613) / '!end of info section' /
DATA oxygen(614) / '54.361             !lower temperature limit [K]' /
DATA oxygen(615) / '154.581            !upper temperature limit [K]' /
DATA oxygen(616) / '0.0                !(dummy) upper pressure limit' /
DATA oxygen(617) / '0.0                !(dummy) maximum density' /
DATA oxygen(618) / '154.581 13.63      !reducing parameters' /
DATA oxygen(619) / '5 0 0 0 0 0        !number of terms in equation' /
DATA oxygen(620) / ' 0.16622D+01   0.345     !coefficients and exponents' /
DATA oxygen(621) / ' 0.76846D+00   0.74' /
DATA oxygen(622) / '-0.10041D+00   1.2' /
DATA oxygen(623) / ' 0.20480D+00   2.6' /
DATA oxygen(624) / ' 0.11551D-01   7.2' /
DATA oxygen(625) / '' /
DATA oxygen(626) / '' /
DATA oxygen(627) / '#DV         !saturated vapor density equation' /
DATA oxygen(628) / 'DV3  saturated vapor density equation of Cullimore (2010).' /
DATA oxygen(629) / '?LITERATURE REFERENCE \' /
DATA oxygen(630) / '?Cullimore, I.D., 2010.' /
DATA oxygen(631) / '?\' /
DATA oxygen(632) / '!end of info section' /
DATA oxygen(633) / '54.361             !lower temperature limit [K]' /
DATA oxygen(634) / '154.581            !upper temperature limit [K]' /
DATA oxygen(635) / '0.0                !(dummy) upper pressure limit' /
DATA oxygen(636) / '0.0                !(dummy) maximum density' /
DATA oxygen(637) / '154.581 13.63      !reducing parameters' /
DATA oxygen(638) / '6 0 0 0 0 0        !number of terms in equation' /
DATA oxygen(639) / '-0.22695D+01   0.3785        !coefficients and exponents' /
DATA oxygen(640) / '-0.46578D+01   1.07' /
DATA oxygen(641) / '-0.99480D+01   2.7' /
DATA oxygen(642) / '-0.22845D+02   5.5' /
DATA oxygen(643) / '-0.45190D+02  10.0' /
DATA oxygen(644) / '-0.25101D+02  20.0' /
DATA oxygen(645) / '' /
DATA oxygen(646) / '' /
DATA oxygen(647) / '@END' /
DATA oxygen(648) / 'c        1         2         3         4         5         6         7         8' /
DATA oxygen(649) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
! #######################################################
character(256), TARGET :: Sulfur_dioxide(537)
DATA Sulfur_dioxide(1) / 'Sulfur dioxide       !Short name' /
DATA Sulfur_dioxide(2) / '7446-09-5            !CAS number' /
DATA Sulfur_dioxide(3) / 'Sulfur dioxide       !Full name' /
DATA Sulfur_dioxide(4) / 'SO2                  !Chemical formula {O2S}' /
DATA Sulfur_dioxide(5) / 'R-764                !Synonym' /
DATA Sulfur_dioxide(6) / '64.0638              !Molar mass [g/mol]' /
DATA Sulfur_dioxide(7) / '197.7                !Triple point temperature [K]' /
DATA Sulfur_dioxide(8) / '263.137              !Normal boiling point [K]' /
DATA Sulfur_dioxide(9) / '430.64               !Critical temperature [K]' /
DATA Sulfur_dioxide(10) / '7886.6               !Critical pressure [kPa]' /
DATA Sulfur_dioxide(11) / '8.078                !Critical density [mol/L]' /
DATA Sulfur_dioxide(12) / '0.256                !Acentric factor' /
DATA Sulfur_dioxide(13) / '1.6                  !Dipole moment [Debye]; Reid, Prausnitz, & Poling, McGraw-Hill (1987)' /
DATA Sulfur_dioxide(14) / 'NBP                  !Default reference state' /
DATA Sulfur_dioxide(15) / '10.0                 !Version number' /
DATA Sulfur_dioxide(16) / '1079                 !UN Number                                                 :UN:' /
DATA Sulfur_dioxide(17) / 'other                !Family                                                    :Family:' /
DATA Sulfur_dioxide(18) / '????                 !Heating value (upper) [kJ/mol]                            :Heat:' /
DATA Sulfur_dioxide(19) / 'B1                   !Safety Group (ASHRAE Standard 34, 2010)                   :Safety:' /
DATA Sulfur_dioxide(20) / '1S/O2S/c1-3-2                             !Standard InChI String                :InChi:' /
DATA Sulfur_dioxide(21) / 'RAHZWNYVWXNFOC-UHFFFAOYSA-N               !Standard InChI Key                   :InChiKey:' /
DATA Sulfur_dioxide(22) / 'e9847540  (ammonia)                       !Alternative fluid for mixing rules   :AltID:' /
DATA Sulfur_dioxide(23) / '7fad4b80                                  !Hash number from InChI Key           :Hash:' /
DATA Sulfur_dioxide(24) / '' /
DATA Sulfur_dioxide(25) / '' /
DATA Sulfur_dioxide(26) / '' /
DATA Sulfur_dioxide(27) / '' /
DATA Sulfur_dioxide(28) / '!The fluid files contain general information about the fluid in the first 15 to 20 lines, followed by sections for the' /
DATA Sulfur_dioxide(29) / '! equations of state, transport equations, and auxiliary equations.  Equations of state are listed first.  The NIST recommended' /
DATA Sulfur_dioxide(30) / '! equations begin with a hash mark (#).  The secondary equations begin with the @ symbol.  These symbols can be swapped to' /
DATA Sulfur_dioxide(31) / '! select a secondary equation as primary and the primary as secondary.  The equation of state section also contains auxiliary' /
DATA Sulfur_dioxide(32) / '! equations for the ideal gas heat capacity or ideal gas Helmholtz energy.  Below the equations of state (both primary and' /
DATA Sulfur_dioxide(33) / '! secondary) are the transport equations, first viscosity and then thermal conductivity.  These are then followed by the' /
DATA Sulfur_dioxide(34) / '! secondary equations if available.  The transport section also contains auxiliary equations required to calculate either the' /
DATA Sulfur_dioxide(35) / '! dilute gas state or the critical enhancement.  At the end of the file are additional but not necessary auxiliary equations,' /
DATA Sulfur_dioxide(36) / '! including simple equations for the vapor pressure, saturated liquid and vapor densities, melting line (for some fluids), and' /
DATA Sulfur_dioxide(37) / '! sublimation line (for even fewer fluids).  This section also contains the equations for dielectric constant and surface' /
DATA Sulfur_dioxide(38) / '! tension if available.  The sections are divided by different symbols (these being _-+=^*~) to aid the eye in locating a' /
DATA Sulfur_dioxide(39) / '! particular section.  Secondary equations are indented 10 spaces to avoid confusion with the NIST recommended equations.  The' /
DATA Sulfur_dioxide(40) / '! end of the fluid file is marked with @END.  Anything below that is ignored.' /
DATA Sulfur_dioxide(41) / '' /
DATA Sulfur_dioxide(42) / '' /
DATA Sulfur_dioxide(43) / '! compiled by E.W. Lemmon, NIST Physical and Chemical Properties Division, Boulder, Colorado' /
DATA Sulfur_dioxide(44) / '! 11-13-98 EWL, Original version.' /
DATA Sulfur_dioxide(45) / '! 01-31-02 EWL, Fit new equation of state based on data of Ihmels.' /
DATA Sulfur_dioxide(46) / '! 04-30-02 EWL, Update fit.' /
DATA Sulfur_dioxide(47) / '! 11-14-02 EWL, Update fit with new PVT data of Ihmels.' /
DATA Sulfur_dioxide(48) / '! 08-17-10 IDC, Add ancillary equations.' /
DATA Sulfur_dioxide(49) / '! 12-06-12 EWL, Add surface tension coefficients of Mulero et al. (2012).' /
DATA Sulfur_dioxide(50) / '! 03-11-13 MLH, Add ECS transport, estimation not std ref quality.' /
DATA Sulfur_dioxide(51) / '! 02-29-16 EWL, Add equation of state of Gao et al. (2016).' /
DATA Sulfur_dioxide(52) / '! 02-28-17 MLH, Revise transport.' /
DATA Sulfur_dioxide(53) / '' /
DATA Sulfur_dioxide(54) / '' /
DATA Sulfur_dioxide(55) / '' /
DATA Sulfur_dioxide(56) / '' /
DATA Sulfur_dioxide(57) / '________________________________________________________________________________' /
DATA Sulfur_dioxide(58) / '' /
DATA Sulfur_dioxide(59) / '#EOS   !---Equation of state---' /
DATA Sulfur_dioxide(60) / 'FEQ    !Helmholtz equation of state for sulfur dioxide of Gao et al. (2016).' /
DATA Sulfur_dioxide(61) / ':TRUECRITICALPOINT:  430.64     8.078         !True EOS critical point [K, mol/L] (where dP/dD=0 and d^2P/dD^2=0 at constant T)' /
DATA Sulfur_dioxide(62) / ':DOI: 10.1021/acs.jced.6b00195' /
DATA Sulfur_dioxide(63) / '?' /
DATA Sulfur_dioxide(64) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(65) / '?Gao, K., Wu, J., Zhang, P., and Lemmon, E.W.,' /
DATA Sulfur_dioxide(66) / '? "A Helmholtz Energy Equation of State for Sulfur Dioxide,"' /
DATA Sulfur_dioxide(67) / '? J. Chem. Eng. Data, 61:2859-2872, 2016. doi: 10.1021/acs.jced.6b00195' /
DATA Sulfur_dioxide(68) / '?' /
DATA Sulfur_dioxide(69) / '?The equation of state is valid from the triple point temperature of 197.7 K to' /
DATA Sulfur_dioxide(70) / '? 523 K, with pressures up to 35 MPa and densities up to 25.3 mol/l.  The' /
DATA Sulfur_dioxide(71) / '? uncertainties in density of the equation of state are 0.1 % in the liquid phase,' /
DATA Sulfur_dioxide(72) / '? 0.25 % in the vapor phase, and 1 % in the critical region.  The uncertainty in' /
DATA Sulfur_dioxide(73) / '? vapor pressure is 0.2 % and the uncertainty in saturated liquid density is 0.2 %' /
DATA Sulfur_dioxide(74) / '? below 410 K. The uncertainty in isobaric heat capacity is 2 % between 200 K and' /
DATA Sulfur_dioxide(75) / '? 290 K. In the critical region, the uncertainties are higher for all properties' /
DATA Sulfur_dioxide(76) / '? except for vapor pressure.' /
DATA Sulfur_dioxide(77) / '?' /
DATA Sulfur_dioxide(78) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(79) / '197.7              !Lower temperature limit [K]' /
DATA Sulfur_dioxide(80) / '525.0              !Upper temperature limit [K]' /
DATA Sulfur_dioxide(81) / '35000.0            !Upper pressure limit [kPa]' /
DATA Sulfur_dioxide(82) / '25.42              !Maximum density [mol/L]' /
DATA Sulfur_dioxide(83) / 'CPP                                    !Pointer to Cp0 model' /
DATA Sulfur_dioxide(84) / '64.0638                                !Molar mass [g/mol]' /
DATA Sulfur_dioxide(85) / '197.7                                  !Triple point temperature [K]' /
DATA Sulfur_dioxide(86) / '1.6661                                 !Pressure at triple point [kPa]' /
DATA Sulfur_dioxide(87) / '25.41                                  !Density at triple point [mol/L]' /
DATA Sulfur_dioxide(88) / '263.137                                !Normal boiling point temperature [K]' /
DATA Sulfur_dioxide(89) / '0.256                                  !Acentric factor' /
DATA Sulfur_dioxide(90) / '430.64        7886.6       8.078       !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Sulfur_dioxide(91) / '430.64                     8.078       !Reducing parameters [K, mol/L]' /
DATA Sulfur_dioxide(92) / '8.3144598                              !Gas constant [J/mol-K]' /
DATA Sulfur_dioxide(93) / '  10  4   6 12   0 0    0 0  0 0  0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Sulfur_dioxide(94) / '  0.01744413   1.0     4.  0.          !a(i),t(i),d(i),l(i)' /
DATA Sulfur_dioxide(95) / '  1.814878     0.45    1.  0.' /
DATA Sulfur_dioxide(96) / ' -2.246338     0.9994  1.  0.' /
DATA Sulfur_dioxide(97) / ' -0.4602906    1.      2.  0.' /
DATA Sulfur_dioxide(98) / '  0.1097049    0.45    3.  0.' /
DATA Sulfur_dioxide(99) / ' -0.9485769    2.907   1.  2.' /
DATA Sulfur_dioxide(100) / ' -0.8751541    2.992   3.  2.' /
DATA Sulfur_dioxide(101) / '  0.4228777    0.87    2.  1.' /
DATA Sulfur_dioxide(102) / ' -0.4174962    3.302   2.  2.' /
DATA Sulfur_dioxide(103) / ' -0.002903451  1.002   7.  1.' /
DATA Sulfur_dioxide(104) / '  1.64041      1.15    1.  2. 2.    -1.061    -0.967   1.276   0.832    0. 0. 0.' /
DATA Sulfur_dioxide(105) / ' -0.4103535    0.997   1.  2. 2.    -0.945    -2.538   0.738   0.69     0. 0. 0.' /
DATA Sulfur_dioxide(106) / ' -0.08316597   1.36    3.  2. 2.    -1.741    -2.758   0.71    0.35     0. 0. 0.' /
DATA Sulfur_dioxide(107) / ' -0.2728126    2.086   2.  2. 2.    -1.139    -1.062   0.997   0.961    0. 0. 0.' /
DATA Sulfur_dioxide(108) / ' -0.1075782    0.855   2.  2. 2.    -1.644    -1.039   1.35    0.981    0. 0. 0.' /
DATA Sulfur_dioxide(109) / ' -0.4348434    0.785   1.  2. 2.    -0.647    -0.41    0.919   0.333    0. 0. 0.' /
DATA Sulfur_dioxide(110) / '                                      eta      beta    gamma   epsilon' /
DATA Sulfur_dioxide(111) / '                                   EXP[eta*(delta-epsilon)^2+beta*(tau-gamma)^2]' /
DATA Sulfur_dioxide(112) / '' /
DATA Sulfur_dioxide(113) / '' /
DATA Sulfur_dioxide(114) / '#AUX   !---Auxiliary function for Cp0' /
DATA Sulfur_dioxide(115) / 'CPP    !Ideal gas heat capacity function for sulfur dioxide of Gao et al. (2016).' /
DATA Sulfur_dioxide(116) / '?' /
DATA Sulfur_dioxide(117) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(118) / '?Gao, K., Wu, J., Zhang, P., and Lemmon, E.W., 2016.' /
DATA Sulfur_dioxide(119) / '?' /
DATA Sulfur_dioxide(120) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(121) / '0.                 !' /
DATA Sulfur_dioxide(122) / '10000.             !' /
DATA Sulfur_dioxide(123) / '0.                 !' /
DATA Sulfur_dioxide(124) / '0.                 !' /
DATA Sulfur_dioxide(125) / '1.0     8.3144598  !Reducing parameters for T, Cp0' /
DATA Sulfur_dioxide(126) / '2 2   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Sulfur_dioxide(127) / ' 4.0        0.0' /
DATA Sulfur_dioxide(128) / ' 0.00007397 1.0' /
DATA Sulfur_dioxide(129) / ' 1.0875     783.0' /
DATA Sulfur_dioxide(130) / ' 1.916      1864.0' /
DATA Sulfur_dioxide(131) / '' /
DATA Sulfur_dioxide(132) / '' /
DATA Sulfur_dioxide(133) / '#AUX   !---Auxiliary function for PX0' /
DATA Sulfur_dioxide(134) / 'PX0    !Helmholtz energy ideal-gas function for sulfur dioxide of Gao et al. (2016).' /
DATA Sulfur_dioxide(135) / '?' /
DATA Sulfur_dioxide(136) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(137) / '?Gao, K., Wu, J., Zhang, P., and Lemmon, E.W., 2016.' /
DATA Sulfur_dioxide(138) / '?' /
DATA Sulfur_dioxide(139) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(140) / '1 3  2  0 0  0 0 0               !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau))' /
DATA Sulfur_dioxide(141) / '  3.0                   1.0      !ai, ti for [ai*log(tau**ti)] terms' /
DATA Sulfur_dioxide(142) / ' -4.541423325625578     0.0      !aj, ti for [ai*tau**ti] terms' /
DATA Sulfur_dioxide(143) / '  4.4732288061807504    1.0      !aj, ti for [ai*tau**ti] terms' /
DATA Sulfur_dioxide(144) / '  0.00007397 -1.0' /
DATA Sulfur_dioxide(145) / '  1.0875     783.0               !aj, ti for [ai*log(1-exp(-ti/T)] terms' /
DATA Sulfur_dioxide(146) / '  1.916      1864.0' /
DATA Sulfur_dioxide(147) / '' /
DATA Sulfur_dioxide(148) / '' /
DATA Sulfur_dioxide(149) / '#AUX   !---Auxiliary function for PH0' /
DATA Sulfur_dioxide(150) / 'PH0    !Ideal gas Helmholtz form for sulfur dioxide of Gao et al. (2016).' /
DATA Sulfur_dioxide(151) / '?' /
DATA Sulfur_dioxide(152) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(153) / '?Gao, K., Wu, J., Zhang, P., and Lemmon, E.W., 2016.' /
DATA Sulfur_dioxide(154) / '?' /
DATA Sulfur_dioxide(155) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(156) / '0.                 !' /
DATA Sulfur_dioxide(157) / '10000.             !' /
DATA Sulfur_dioxide(158) / '0.                 !' /
DATA Sulfur_dioxide(159) / '0.                 !' /
DATA Sulfur_dioxide(160) / '1 3  2  0 0  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA Sulfur_dioxide(161) / ' 3.0               1.0                 !ai, ti for [ai*log(tau**ti)] terms' /
DATA Sulfur_dioxide(162) / '-4.5414235721      0.0                 !aj, ti for [ai*tau**ti] terms' /
DATA Sulfur_dioxide(163) / ' 4.4732289572      1.0' /
DATA Sulfur_dioxide(164) / '-0.0159272204     -1.0' /
DATA Sulfur_dioxide(165) / ' 1.0875           -1.8182240386        !aj, ti for [ai*log(1-exp(ti*tau)] terms' /
DATA Sulfur_dioxide(166) / ' 1.916            -4.3284413896' /
DATA Sulfur_dioxide(167) / '' /
DATA Sulfur_dioxide(168) / '' /
DATA Sulfur_dioxide(169) / '' /
DATA Sulfur_dioxide(170) / '' /
DATA Sulfur_dioxide(171) / '--------------------------------------------------------------------------------' /
DATA Sulfur_dioxide(172) / '' /
DATA Sulfur_dioxide(173) / '@EOS    !---Equation of state---' /
DATA Sulfur_dioxide(174) / 'FE1     !Helmholtz equation of state for sulfur dioxide of Lemmon and Span (2006).' /
DATA Sulfur_dioxide(175) / '          ?' /
DATA Sulfur_dioxide(176) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(177) / '          ?Lemmon, E.W. and Span, R.,' /
DATA Sulfur_dioxide(178) / '          ? "Short Fundamental Equations of State for 20 Industrial Fluids,"' /
DATA Sulfur_dioxide(179) / '          ? J. Chem. Eng. Data, 51(3):785-850, 2006. doi: 10.1021/je050186n' /
DATA Sulfur_dioxide(180) / '          ?' /
DATA Sulfur_dioxide(181) / '          ?see also:' /
DATA Sulfur_dioxide(182) / '          ? Ihmels, E.C., Lemmon, E.W., Gmehling, J.,' /
DATA Sulfur_dioxide(183) / '          ? "An Equation of State and Compressed Liquid and Supercritical Densities for' /
DATA Sulfur_dioxide(184) / '          ? Sulfur Dioxide,"' /
DATA Sulfur_dioxide(185) / '          ? Fluid Phase Equilib., 207:111-130, 2003. doi: 10.1016/S0378-3812(03)00004-9' /
DATA Sulfur_dioxide(186) / '          ?' /
DATA Sulfur_dioxide(187) / '          ?The uncertainty in density of the equation of state ranges from 0.1% at low' /
DATA Sulfur_dioxide(188) / '          ? temperatures in the liquid and vapor to 0.5% at the highest temperatures.' /
DATA Sulfur_dioxide(189) / '          ? The uncertainty in heat capacities is 2%, and the uncertainty in vapor' /
DATA Sulfur_dioxide(190) / '          ? pressure is 0.4% at temperatures above 270 K. The uncertainty in vapor' /
DATA Sulfur_dioxide(191) / '          ? pressure increases at lower temperatures due to the lack of reliable' /
DATA Sulfur_dioxide(192) / '          ? experimental data.  In the critical region, the uncertainties are higher' /
DATA Sulfur_dioxide(193) / '          ? for all properties except vapor pressure.' /
DATA Sulfur_dioxide(194) / '          ?' /
DATA Sulfur_dioxide(195) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(196) / '          197.7              !Lower temperature limit [K]' /
DATA Sulfur_dioxide(197) / '          525.0              !Upper temperature limit [K]' /
DATA Sulfur_dioxide(198) / '          35000.0            !Upper pressure limit [kPa]' /
DATA Sulfur_dioxide(199) / '          25.30              !Maximum density [mol/L]' /
DATA Sulfur_dioxide(200) / '          CP1                                    !Pointer to Cp0 model' /
DATA Sulfur_dioxide(201) / '          64.0638                                !Molar mass [g/mol]' /
DATA Sulfur_dioxide(202) / '          197.7                                  !Triple point temperature [K]' /
DATA Sulfur_dioxide(203) / '          1.66                                   !Pressure at triple point [kPa]' /
DATA Sulfur_dioxide(204) / '          25.29                                  !Density at triple point [mol/L]' /
DATA Sulfur_dioxide(205) / '          263.13                                 !Normal boiling point temperature [K]' /
DATA Sulfur_dioxide(206) / '          0.2557                                 !Acentric factor' /
DATA Sulfur_dioxide(207) / '          430.64        7884.0       8.195       !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Sulfur_dioxide(208) / '          430.64                     8.195       !Reducing parameters [K, mol/L]' /
DATA Sulfur_dioxide(209) / '          8.314472                               !Gas constant [J/mol-K]' /
DATA Sulfur_dioxide(210) / '            12  4    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Sulfur_dioxide(211) / '           0.93061     0.25    1.  0.            !a(i),t(i),d(i),l(i)' /
DATA Sulfur_dioxide(212) / '          -1.9528      1.25    1.  0.' /
DATA Sulfur_dioxide(213) / '          -0.17467     1.5     1.  0.' /
DATA Sulfur_dioxide(214) / '           0.061524    0.25    3.  0.' /
DATA Sulfur_dioxide(215) / '           0.00017711  0.875   7.  0.' /
DATA Sulfur_dioxide(216) / '           0.21615     2.375   1.  1.' /
DATA Sulfur_dioxide(217) / '           0.51353     2.0     2.  1.' /
DATA Sulfur_dioxide(218) / '           0.010419    2.125   5.  1.' /
DATA Sulfur_dioxide(219) / '          -0.25286     3.5     1.  2.' /
DATA Sulfur_dioxide(220) / '          -0.054720    6.5     1.  2.' /
DATA Sulfur_dioxide(221) / '          -0.059856    4.75    4.  2.' /
DATA Sulfur_dioxide(222) / '          -0.016523   12.5     2.  3.' /
DATA Sulfur_dioxide(223) / '' /
DATA Sulfur_dioxide(224) / '' /
DATA Sulfur_dioxide(225) / '@AUX    !---Auxiliary function for Cp0' /
DATA Sulfur_dioxide(226) / 'CP1     !Ideal gas heat capacity function for sulfur dioxide.' /
DATA Sulfur_dioxide(227) / '          ?' /
DATA Sulfur_dioxide(228) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(229) / '          ?Lemmon, E.W. and Span, R.' /
DATA Sulfur_dioxide(230) / '          ?' /
DATA Sulfur_dioxide(231) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(232) / '          0.                 !' /
DATA Sulfur_dioxide(233) / '          10000.             !' /
DATA Sulfur_dioxide(234) / '          0.                 !' /
DATA Sulfur_dioxide(235) / '          0.                 !' /
DATA Sulfur_dioxide(236) / '          1.0     8.314472   !Reducing parameters for T, Cp0' /
DATA Sulfur_dioxide(237) / '          2 2   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Sulfur_dioxide(238) / '           4.0               0.0' /
DATA Sulfur_dioxide(239) / '           0.000072453       1.0' /
DATA Sulfur_dioxide(240) / '           1.062             775.0' /
DATA Sulfur_dioxide(241) / '           1.9401            1851.0' /
DATA Sulfur_dioxide(242) / '' /
DATA Sulfur_dioxide(243) / '' /
DATA Sulfur_dioxide(244) / '@AUX    !---Auxiliary function for PH0' /
DATA Sulfur_dioxide(245) / 'PH1     !Ideal gas Helmholtz form for sulfur dioxide.' /
DATA Sulfur_dioxide(246) / '          ?' /
DATA Sulfur_dioxide(247) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(248) / '          ?Lemmon, E.W. and Span, R.' /
DATA Sulfur_dioxide(249) / '          ?' /
DATA Sulfur_dioxide(250) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(251) / '          0.                 !' /
DATA Sulfur_dioxide(252) / '          10000.             !' /
DATA Sulfur_dioxide(253) / '          0.                 !' /
DATA Sulfur_dioxide(254) / '          0.                 !' /
DATA Sulfur_dioxide(255) / '          1 3  2 0 0  0 0 0  !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA Sulfur_dioxide(256) / '           3.0               1.0                 !ai, ti for [ai*log(tau**ti)] terms' /
DATA Sulfur_dioxide(257) / '          -4.5328346436      0.0                 !aj, ti for [ai*tau**ti] terms' /
DATA Sulfur_dioxide(258) / '           4.4777967379      1.0' /
DATA Sulfur_dioxide(259) / '          -0.01560058       -1.0' /
DATA Sulfur_dioxide(260) / '           1.062            -1.799647037         !aj, ti for [ai*log(1-exp(ti*tau)] terms' /
DATA Sulfur_dioxide(261) / '           1.9401           -4.2982537618' /
DATA Sulfur_dioxide(262) / '' /
DATA Sulfur_dioxide(263) / '' /
DATA Sulfur_dioxide(264) / '@EOS    !---Equation of state---' /
DATA Sulfur_dioxide(265) / 'FE2     !Helmholtz equation of state for sulfur dioxide of Polt (1987).' /
DATA Sulfur_dioxide(266) / '          ?' /
DATA Sulfur_dioxide(267) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(268) / '          ?Polt, A.,' /
DATA Sulfur_dioxide(269) / '          ? Zur Beschreibung der thermodynamischen Eigenschaften reiner Fluide' /
DATA Sulfur_dioxide(270) / '          ? mit "Erweiterten BWR-Gleichungen",' /
DATA Sulfur_dioxide(271) / '          ? Ph.D. Dissertation, Universitaet Kaiserslautern, Germany, 1987.' /
DATA Sulfur_dioxide(272) / '          ?' /
DATA Sulfur_dioxide(273) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(274) / '          273.0              !Lower temperature limit [K]' /
DATA Sulfur_dioxide(275) / '          523.0              !Upper temperature limit [K]' /
DATA Sulfur_dioxide(276) / '          32000.0            !Upper pressure limit [kPa]' /
DATA Sulfur_dioxide(277) / '          22.91              !Maximum density [mol/L]' /
DATA Sulfur_dioxide(278) / '          CP2                                    !Pointer to Cp0 model' /
DATA Sulfur_dioxide(279) / '          64.066                                 !Molar mass [g/mol]' /
DATA Sulfur_dioxide(280) / '          197.7                                  !Triple point temperature [K]' /
DATA Sulfur_dioxide(281) / '          11.82                                  !Pressure at triple point [kPa]' /
DATA Sulfur_dioxide(282) / '          23.0                                   !Density at triple point [mol/L] (unknown)' /
DATA Sulfur_dioxide(283) / '          256.61                                 !Normal boiling point temperature [K]' /
DATA Sulfur_dioxide(284) / '          0.23                                   !Acentric factor' /
DATA Sulfur_dioxide(285) / '          430.65        7880.0       8.1946742   !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA Sulfur_dioxide(286) / '          430.65                     8.1946742   !Reducing parameters [K, mol/L]' /
DATA Sulfur_dioxide(287) / '          8.3143                                 !Gas constant [J/mol-K]' /
DATA Sulfur_dioxide(288) / '            22  4    0  0    0 0    0 0 0 0 0 0  !# terms and # coefs/term for normal terms, Gaussian terms, and Gao terms' /
DATA Sulfur_dioxide(289) / '           0.789407019882      3.0       0.  0.  !a(i),t(i),d(i),l(i)' /
DATA Sulfur_dioxide(290) / '          -1.70449580056       4.0       0.  0.' /
DATA Sulfur_dioxide(291) / '           1.15984637964       5.0       0.  0.' /
DATA Sulfur_dioxide(292) / '          -0.576307837294      0.0       1.  0.' /
DATA Sulfur_dioxide(293) / '           2.49237283833       1.0       1.  0.' /
DATA Sulfur_dioxide(294) / '          -5.18115678632       2.0       1.  0.' /
DATA Sulfur_dioxide(295) / '           3.20766081899       3.0       1.  0.' /
DATA Sulfur_dioxide(296) / '          -1.23636065893       4.0       1.  0.' /
DATA Sulfur_dioxide(297) / '           0.0144419600938     0.0       2.  0.' /
DATA Sulfur_dioxide(298) / '          -0.153807055040      1.0       2.  0.' /
DATA Sulfur_dioxide(299) / '           0.386324300525      2.0       2.  0.' /
DATA Sulfur_dioxide(300) / '           0.292550313202      0.0       3.  0.' /
DATA Sulfur_dioxide(301) / '          -0.372445361392      1.0       3.  0.' /
DATA Sulfur_dioxide(302) / '          -0.063692433391      0.0       4.  0.' /
DATA Sulfur_dioxide(303) / '           0.0986166451596     1.0       4.  0.' /
DATA Sulfur_dioxide(304) / '          -0.00216993783055    1.0       5.  0.' /
DATA Sulfur_dioxide(305) / '          -0.789407019882      3.0       0.  2.' /
DATA Sulfur_dioxide(306) / '           1.70449580056       4.0       0.  2.' /
DATA Sulfur_dioxide(307) / '          -1.15984637964       5.0       0.  2.' /
DATA Sulfur_dioxide(308) / '          -0.480876182378      3.0       2.  2.' /
DATA Sulfur_dioxide(309) / '           1.64910076886       4.0       2.  2.' /
DATA Sulfur_dioxide(310) / '          -1.33861069604       5.0       2.  2.' /
DATA Sulfur_dioxide(311) / '' /
DATA Sulfur_dioxide(312) / '' /
DATA Sulfur_dioxide(313) / '@AUX    !---Auxiliary function for Cp0' /
DATA Sulfur_dioxide(314) / 'CP2     !Ideal gas heat capacity function for sulfur dioxide.' /
DATA Sulfur_dioxide(315) / '          ?' /
DATA Sulfur_dioxide(316) / '          ?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(317) / '          ?Polt, A.,' /
DATA Sulfur_dioxide(318) / '          ?' /
DATA Sulfur_dioxide(319) / '          !```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(320) / '          0.                 !' /
DATA Sulfur_dioxide(321) / '          10000.             !' /
DATA Sulfur_dioxide(322) / '          0.                 !' /
DATA Sulfur_dioxide(323) / '          0.                 !' /
DATA Sulfur_dioxide(324) / '          1.0     64.066     !Reducing parameters for T, Cp0' /
DATA Sulfur_dioxide(325) / '          5 0   0 0   0 0 0  !Nterms:  polynomial, exponential, cosh, sinh' /
DATA Sulfur_dioxide(326) / '           0.4021066         0.0' /
DATA Sulfur_dioxide(327) / '           0.0008734857      1.0' /
DATA Sulfur_dioxide(328) / '          -0.4596882e-6      2.0' /
DATA Sulfur_dioxide(329) / '          -0.133284e-11      3.0' /
DATA Sulfur_dioxide(330) / '           0.23785e-13       4.0' /
DATA Sulfur_dioxide(331) / '' /
DATA Sulfur_dioxide(332) / '' /
DATA Sulfur_dioxide(333) / '' /
DATA Sulfur_dioxide(334) / '' /
DATA Sulfur_dioxide(335) / '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' /
DATA Sulfur_dioxide(336) / '' /
DATA Sulfur_dioxide(337) / '#TRN   !---ECS Transport---' /
DATA Sulfur_dioxide(338) / 'ECS    !Extended Corresponding States model (Propane reference); fitted to very limited data for sulfur dioxide.' /
DATA Sulfur_dioxide(339) / ':DOI: 10.6028/NIST.IR.8209' /
DATA Sulfur_dioxide(340) / '?' /
DATA Sulfur_dioxide(341) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(342) / '?Huber, M.L., "Models for the Viscosity, Thermal Conductivity, and Surface Tension' /
DATA Sulfur_dioxide(343) / '? of Selected Pure Fluids as Implemented in REFPROP v10.0," NISTIR 8209, 2018.' /
DATA Sulfur_dioxide(344) / '? doi: 10.6028/NIST.IR.8209' /
DATA Sulfur_dioxide(345) / '?' /
DATA Sulfur_dioxide(346) / '?THERMAL CONDUCTIVITY' /
DATA Sulfur_dioxide(347) / '? The parameters for thermal conductivity were based in part on the data of:' /
DATA Sulfur_dioxide(348) / '? Kardos, A., "The Heat Conductivities of Various Liquids," Z. ges. Kalte-Ind., 41, 29-35, 1934.' /
DATA Sulfur_dioxide(349) / '? Baker, C.B. and de Haas, N., "Gas Thermal Conductivity Studies at High Temperature. III. Results for SO2," Phys. Fluids, 7:1400-1402, 1964. doi: 10.1063/1.1711394' /
DATA Sulfur_dioxide(350) / '?' /
DATA Sulfur_dioxide(351) / '?The estimated uncertainty of thermal conductivity in the liquid phase along the saturation boundary is 5%, rising to 10% at 35 MPa.' /
DATA Sulfur_dioxide(352) / '? The estimated uncertainty of thermal conductivity in the gas phase is 5%.' /
DATA Sulfur_dioxide(353) / '?' /
DATA Sulfur_dioxide(354) / '?VISCOSITY' /
DATA Sulfur_dioxide(355) / '? The ECS parameters for viscosity were based in part on the data of:' /
DATA Sulfur_dioxide(356) / '? Hartl, R., Neueder, R., and Gores, H.J., "Temperature Dependence of Association Constants of LiAlCl4 in Liquid Sulfur Dioxide," Acta Chim. Slov., 56:109-114, 2009.' /
DATA Sulfur_dioxide(357) / '? Awbery, J.H. and Griffiths, E., "The Viscosities of Some Liquid Refrigerants," Proc. Phys. Soc., London, 48:372-80, 1936.' /
DATA Sulfur_dioxide(358) / '? Stewart, W.W. and Maass, "The Coefficient of Viscosity of Sulphur Dioxide over a Low Temperature Range," Can. J. of Research, 6(5):453-457, 1932. doi: 10.1139/cjr32-035' /
DATA Sulfur_dioxide(359) / '?' /
DATA Sulfur_dioxide(360) / '?The estimated uncertainty of viscosity in the liquid phase along the saturation boundary is 5%, rising to 10% at 35 MPa.' /
DATA Sulfur_dioxide(361) / '?The estimated uncertainty of viscosity in the gas phase is 5%.' /
DATA Sulfur_dioxide(362) / '?' /
DATA Sulfur_dioxide(363) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(364) / '197.7              !Lower temperature limit [K]' /
DATA Sulfur_dioxide(365) / '525.0              !Upper temperature limit [K]' /
DATA Sulfur_dioxide(366) / '35000.0            !Upper pressure limit [kPa]' /
DATA Sulfur_dioxide(367) / '25.42              !Maximum density [mol/L]' /
DATA Sulfur_dioxide(368) / 'FEQ PROPANE.FLD' /
DATA Sulfur_dioxide(369) / 'VS1                !Model for reference fluid viscosity' /
DATA Sulfur_dioxide(370) / 'TC1                !Model for reference fluid thermal conductivity' /
DATA Sulfur_dioxide(371) / 'NUL                !Large molecule identifier' /
DATA Sulfur_dioxide(372) / '1                  !Lennard-Jones flag (0 or 1) (0 => use estimates) !from scaling R134a' /
DATA Sulfur_dioxide(373) / '0.4026             !Lennard-Jones coefficient sigma [nm] for ECS method' /
DATA Sulfur_dioxide(374) / '363.0              !Lennard-Jones coefficient epsilon/kappa [K] for ECS method' /
DATA Sulfur_dioxide(375) / '2  0  0                  !Number of terms in f_int term in Eucken correlation, spare1, spare2' /
DATA Sulfur_dioxide(376) / ' 6.60505e-4    0. 0. 0.  !Coefficient, power of T, spare1, spare2' /
DATA Sulfur_dioxide(377) / ' 7.47059e-7    1. 0. 0.  !Coefficient, power of T, spare1, spare2' /
DATA Sulfur_dioxide(378) / '2  0  0                  !Number of terms in psi (visc shape factor): poly,spare1,spare2' /
DATA Sulfur_dioxide(379) / ' 0.917778      0. 0. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Sulfur_dioxide(380) / ' 0.0248405     0. 1. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Sulfur_dioxide(381) / '2  0  0                  !Number of terms in chi (t.c. shape factor): poly,spare1,spare2' /
DATA Sulfur_dioxide(382) / ' 1.38755       0. 0. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Sulfur_dioxide(383) / '-0.128721      0. 1. 0.  !Coefficient, power of Tr, power of Dr, spare' /
DATA Sulfur_dioxide(384) / 'TK3                !Pointer to critical enhancement auxiliary function' /
DATA Sulfur_dioxide(385) / '' /
DATA Sulfur_dioxide(386) / '' /
DATA Sulfur_dioxide(387) / '#AUX   !---Auxiliary function for the thermal conductivity critical enhancement' /
DATA Sulfur_dioxide(388) / 'TK3    !Simplified thermal conductivity critical enhancement for sulfur dioxide of Perkins et al. (2013).' /
DATA Sulfur_dioxide(389) / '?' /
DATA Sulfur_dioxide(390) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(391) / '?Perkins, R.A., Sengers, J.V., Abdulagatov, I.M., and Huber, M.L.,' /
DATA Sulfur_dioxide(392) / '? "Simplified Model for the Critical Thermal-Conductivity Enhancement in Molecular Fluids,"' /
DATA Sulfur_dioxide(393) / '? Int. J. Thermophys., 34(2):191-212, 2013. doi: 10.1007/s10765-013-1409-z' /
DATA Sulfur_dioxide(394) / '?' /
DATA Sulfur_dioxide(395) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(396) / '0.                 !' /
DATA Sulfur_dioxide(397) / '10000.             !' /
DATA Sulfur_dioxide(398) / '0.                 !' /
DATA Sulfur_dioxide(399) / '0.                 !' /
DATA Sulfur_dioxide(400) / '9 0 0 0            !# terms:  CO2-terms, spare, spare, spare' /
DATA Sulfur_dioxide(401) / '1.0  1.0  1.0      !Reducing parameters for T, rho, tcx [mW/(m-K)]' /
DATA Sulfur_dioxide(402) / '0.63               !Nu (universal exponent)' /
DATA Sulfur_dioxide(403) / '1.239              !Gamma (universal exponent)' /
DATA Sulfur_dioxide(404) / '1.02               !R0 (universal amplitude)' /
DATA Sulfur_dioxide(405) / '0.063              !Z (universal exponent--not used for t.c., only viscosity)' /
DATA Sulfur_dioxide(406) / '1.0                !C (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA Sulfur_dioxide(407) / '0.167e-9           !Xi0 (amplitude) [m]' /
DATA Sulfur_dioxide(408) / '0.059              !Gam0 (amplitude) [-]' /
DATA Sulfur_dioxide(409) / '0.485e-9           !Qd_inverse (modified effective cutoff parameter) [m]' /
DATA Sulfur_dioxide(410) / '645.96             !Tref (reference temperature)=1.5*Tc [K]' /
DATA Sulfur_dioxide(411) / '' /
DATA Sulfur_dioxide(412) / '' /
DATA Sulfur_dioxide(413) / '' /
DATA Sulfur_dioxide(414) / '' /
DATA Sulfur_dioxide(415) / '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' /
DATA Sulfur_dioxide(416) / '' /
DATA Sulfur_dioxide(417) / '#DE    !---Dielectric constant---' /
DATA Sulfur_dioxide(418) / 'DE5    !Dielectric constant model for SO2 of Harvey and Mountain (2017).' /
DATA Sulfur_dioxide(419) / ':DOI: 10.1007/s10765-017-2279-6' /
DATA Sulfur_dioxide(420) / '?' /
DATA Sulfur_dioxide(421) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(422) / '?Harvey, A.H. and Mountain, R.D.,' /
DATA Sulfur_dioxide(423) / '? "Correlations for the Dielectric Constants of H2S, SO2, and SF6,"' /
DATA Sulfur_dioxide(424) / '? Int. J. Thermophys., 38:147, 2017.' /
DATA Sulfur_dioxide(425) / '?' /
DATA Sulfur_dioxide(426) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(427) / '197.7              !Lower temperature limit [K]' /
DATA Sulfur_dioxide(428) / '525.0              !Upper temperature limit [K]' /
DATA Sulfur_dioxide(429) / '0.                 !' /
DATA Sulfur_dioxide(430) / '0.                 !' /
DATA Sulfur_dioxide(431) / '450.0 25.0 1.0     !Reducing parameters for T and D' /
DATA Sulfur_dioxide(432) / '4 1 1 1 0 0        !Number of terms in dielectric constant model' /
DATA Sulfur_dioxide(433) / ' 4.09e-24  0.0   0.0    0.    !  alpha (cm^3)' /
DATA Sulfur_dioxide(434) / ' 1.63308   0.0   0.0    0.    !  mu (debye)' /
DATA Sulfur_dioxide(435) / ' 0.335     0.0   0.0    0.    !  cu' /
DATA Sulfur_dioxide(436) / ' 2.516     0.0   0.0    0.    !  cg' /
DATA Sulfur_dioxide(437) / ' 1.        0.0   1.018  0.75  !  f' /
DATA Sulfur_dioxide(438) / ' 1.        0.0   0.8972 0.98  !  g1' /
DATA Sulfur_dioxide(439) / ' 1.        0.0   0.5264 1.2   !  g2' /
DATA Sulfur_dioxide(440) / '' /
DATA Sulfur_dioxide(441) / '' /
DATA Sulfur_dioxide(442) / '#STN   !---Surface tension---' /
DATA Sulfur_dioxide(443) / 'ST1    !Surface tension model for sulfur dioxide of Mulero et al. (2012).' /
DATA Sulfur_dioxide(444) / ':DOI: 10.1063/1.4768782' /
DATA Sulfur_dioxide(445) / '?' /
DATA Sulfur_dioxide(446) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(447) / '?Mulero, A., Cachadina, I., and Parra, M.I.,' /
DATA Sulfur_dioxide(448) / '? "Recommended Correlations for the Surface Tension of Common Fluids,"' /
DATA Sulfur_dioxide(449) / '? J. Phys. Chem. Ref. Data, 41(4), 043105, 2012. doi: 10.1063/1.4768782' /
DATA Sulfur_dioxide(450) / '?' /
DATA Sulfur_dioxide(451) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(452) / '0.                 !' /
DATA Sulfur_dioxide(453) / '10000.             !' /
DATA Sulfur_dioxide(454) / '0.                 !' /
DATA Sulfur_dioxide(455) / '0.                 !' /
DATA Sulfur_dioxide(456) / '3                  !Number of terms in surface tension model' /
DATA Sulfur_dioxide(457) / '430.64             !Critical temperature used in fit (dummy)' /
DATA Sulfur_dioxide(458) / ' 0.0803    0.928   !Sigma0 and n' /
DATA Sulfur_dioxide(459) / ' 0.0139    1.57' /
DATA Sulfur_dioxide(460) / '-0.0114    0.364' /
DATA Sulfur_dioxide(461) / '' /
DATA Sulfur_dioxide(462) / '' /
DATA Sulfur_dioxide(463) / '#PS    !---Vapor pressure---' /
DATA Sulfur_dioxide(464) / 'PS5    !Vapor pressure equation for sulfur dioxide of Gao et al. (2016).' /
DATA Sulfur_dioxide(465) / '?' /
DATA Sulfur_dioxide(466) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(467) / '?Gao, K., Wu, J., Zhang, P., and Lemmon, E.W., 2016.' /
DATA Sulfur_dioxide(468) / '?' /
DATA Sulfur_dioxide(469) / '?Functional Form:  P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T] where Theta=1-T/Tc, Tc and Pc' /
DATA Sulfur_dioxide(470) / '? are the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Sulfur_dioxide(471) / '?' /
DATA Sulfur_dioxide(472) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(473) / '0.                 !' /
DATA Sulfur_dioxide(474) / '10000.             !' /
DATA Sulfur_dioxide(475) / '0.                 !' /
DATA Sulfur_dioxide(476) / '0.                 !' /
DATA Sulfur_dioxide(477) / '430.64  7886.6     !Reducing parameters' /
DATA Sulfur_dioxide(478) / '5 0 0 0 0 0        !Number of terms in equation' /
DATA Sulfur_dioxide(479) / '-7.303     1.0' /
DATA Sulfur_dioxide(480) / ' 1.9794    1.5' /
DATA Sulfur_dioxide(481) / '-2.078     2.2' /
DATA Sulfur_dioxide(482) / '-3.5446    4.7' /
DATA Sulfur_dioxide(483) / ' 0.51776   6.0' /
DATA Sulfur_dioxide(484) / '' /
DATA Sulfur_dioxide(485) / '' /
DATA Sulfur_dioxide(486) / '#DL    !---Saturated liquid density---' /
DATA Sulfur_dioxide(487) / 'DL1    !Saturated liquid density equation for sulfur dioxide of Gao et al. (2016).' /
DATA Sulfur_dioxide(488) / '?' /
DATA Sulfur_dioxide(489) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(490) / '?Gao, K., Wu, J., Zhang, P., and Lemmon, E.W., 2016.' /
DATA Sulfur_dioxide(491) / '?' /
DATA Sulfur_dioxide(492) / '?Functional Form:  D=Dc*[1+SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Sulfur_dioxide(493) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Sulfur_dioxide(494) / '?' /
DATA Sulfur_dioxide(495) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(496) / '0.                 !' /
DATA Sulfur_dioxide(497) / '10000.             !' /
DATA Sulfur_dioxide(498) / '0.                 !' /
DATA Sulfur_dioxide(499) / '0.                 !' /
DATA Sulfur_dioxide(500) / '430.64  8.078      !Reducing parameters' /
DATA Sulfur_dioxide(501) / '5 0 0 0 0 0        !Number of terms in equation' /
DATA Sulfur_dioxide(502) / ' 7.2296    0.54' /
DATA Sulfur_dioxide(503) / '-16.928    0.88' /
DATA Sulfur_dioxide(504) / ' 29.832    1.23' /
DATA Sulfur_dioxide(505) / '-27.901    1.6' /
DATA Sulfur_dioxide(506) / ' 11.085    2.0' /
DATA Sulfur_dioxide(507) / '' /
DATA Sulfur_dioxide(508) / '' /
DATA Sulfur_dioxide(509) / '#DV    !---Saturated vapor density---' /
DATA Sulfur_dioxide(510) / 'DV3    !Saturated vapor density equation for sulfur dioxide of Gao et al. (2016).' /
DATA Sulfur_dioxide(511) / '?' /
DATA Sulfur_dioxide(512) / '?```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(513) / '?Gao, K., Wu, J., Zhang, P., and Lemmon, E.W., 2016.' /
DATA Sulfur_dioxide(514) / '?' /
DATA Sulfur_dioxide(515) / '?Functional Form:  D=Dc*EXP[SUM(Ni*Theta^ti)] where Theta=1-T/Tc, Tc and Dc are' /
DATA Sulfur_dioxide(516) / '? the reducing parameters below, which are followed by rows containing Ni and ti.' /
DATA Sulfur_dioxide(517) / '?' /
DATA Sulfur_dioxide(518) / '!```````````````````````````````````````````````````````````````````````````````' /
DATA Sulfur_dioxide(519) / '0.                 !' /
DATA Sulfur_dioxide(520) / '10000.             !' /
DATA Sulfur_dioxide(521) / '0.                 !' /
DATA Sulfur_dioxide(522) / '0.                 !' /
DATA Sulfur_dioxide(523) / '430.64  8.078      !Reducing parameters' /
DATA Sulfur_dioxide(524) / '6 0 0 0 0 0        !Number of terms in equation' /
DATA Sulfur_dioxide(525) / '-7.487     0.545' /
DATA Sulfur_dioxide(526) / ' 10.118    0.85' /
DATA Sulfur_dioxide(527) / '-13.608    1.2' /
DATA Sulfur_dioxide(528) / '-25.408    3.7' /
DATA Sulfur_dioxide(529) / '-42.04     7.5' /
DATA Sulfur_dioxide(530) / '-38.668    10.0' /
DATA Sulfur_dioxide(531) / '' /
DATA Sulfur_dioxide(532) / '' /
DATA Sulfur_dioxide(533) / '@END' /
DATA Sulfur_dioxide(534) / 'c        1         2         3         4         5         6         7         8' /
DATA Sulfur_dioxide(535) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
DATA Sulfur_dioxide(536) / '                    0.4112             !Lennard-Jones coefficient sigma [nm] for ECS method' /
DATA Sulfur_dioxide(537) / '                    335.4              !Lennard-Jones coefficient epsilon/kappa [K] for ECS method' /
! #######################################################
character(256), TARGET :: water(654)
DATA water(1) / 'water              !short name' /
DATA water(2) / '7732-18-5          !CAS number' /
DATA water(3) / 'water              !full name' /
DATA water(4) / 'H2O                !chemical formula' /
DATA water(5) / 'R-718              !synonym' /
DATA water(6) / '18.015268          !molecular weight [g/mol]' /
DATA water(7) / '273.16             !triple point temperature [K]' /
DATA water(8) / '373.1243           !normal boiling point [K]' /
DATA water(9) / '647.096            !critical temperature [K]' /
DATA water(10) / '22064.0            !critical pressure [kPa]' /
DATA water(11) / '17.8737279956      !critical density [mol/L]' /
DATA water(12) / '0.3443             !acentric factor' /
DATA water(13) / '1.855              !dipole moment [Debye]' /
DATA water(14) / 'OTH                !default reference state' /
DATA water(15) / '300.0   1.0   45957.191490119709 164.00522417832 !tref, Pref, Href, Sref (corresponds to u,s = 0 @ Ttp)' /
DATA water(16) / '9.1                !version number' /
DATA water(17) / 'other              !family' /
DATA water(18) / '44.016             !heating value (gross or superior) [kJ/mol]' /
DATA water(19) / 'A1                 !Safety Group (ASHRAE Standard 34, 2010)' /
DATA water(20) / '' /
DATA water(21) / '' /
DATA water(22) / '#EOS               !equation of state specification' /
DATA water(23) / 'FEQ  Helmholtz equation of state for water of Wagner and Pruss (2002).' /
DATA water(24) / '?LITERATURE REFERENCE \' /
DATA water(25) / '?Wagner, W. and Pruss, A.,' /
DATA water(26) / '? "The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary' /
DATA water(27) / '? Water Substance for General and Scientific Use,"' /
DATA water(28) / '? J. Phys. Chem. Ref. Data, 31(2):387-535, 2002.' /
DATA water(29) / '?\' /
DATA water(30) / '?The uncertainty in density of the equation of state is 0.0001% at 1 atm' /
DATA water(31) / '?in the liquid phase, and 0.001% at other liquid states at pressures up' /
DATA water(32) / '?to 10 MPa and temperatures to 423 K. In the vapor phase, the uncertainty' /
DATA water(33) / '?is 0.05% or less.  The uncertainties rise at higher temperatures and/or' /
DATA water(34) / '?pressures, but are generally less than 0.1% in density except at extreme' /
DATA water(35) / '?conditions.  The uncertainty in pressure in the critical region is 0.1%.' /
DATA water(36) / '?The uncertainty of the speed of sound is 0.15% in the vapor and 0.1% or' /
DATA water(37) / '?less in the liquid, and increases near the critical region and at high' /
DATA water(38) / '?temperatures and pressures.  The uncertainty in isobaric heat capacity' /
DATA water(39) / '?is 0.2% in the vapor and 0.1% in the liquid, with increasing values in' /
DATA water(40) / '?the critical region and at high pressures.  The uncertainties of' /
DATA water(41) / '?saturation conditions are 0.025% in vapor pressure, 0.0025% in' /
DATA water(42) / '?saturated liquid density, and 0.1% in saturated vapor density.  The' /
DATA water(43) / '?uncertainties in the saturated densities increase substantially as the' /
DATA water(44) / '?critical region is approached.' /
DATA water(45) / '?\' /
DATA water(46) / '!end of info section' /
DATA water(47) / '273.16             !lower temperature limit [K]' /
DATA water(48) / '2000.0             !upper temperature limit [K]' /
DATA water(49) / '1000000.0          !upper pressure limit [kPa]' /
DATA water(50) / '73.96              !maximum density [mol/L]' /
DATA water(51) / 'CPP                                    !pointer to Cp0 model' /
DATA water(52) / '18.015268                              !molecular weight [g/mol]' /
DATA water(53) / '273.16                                 !triple point temperature [K]' /
DATA water(54) / '0.6116                                 !pressure at triple point [kPa]' /
DATA water(55) / '55.49696                               !density at triple point [mol/L]' /
DATA water(56) / '373.1243                               !normal boiling point temperature [K]' /
DATA water(57) / '0.3443                                 !acentric factor' /
DATA water(58) / '647.096      22064.0     17.8737279956 !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA water(59) / '647.096                  17.8737279956 !reducing parameters [K, mol/L]' /
DATA water(60) / '8.314371357587                         !gas constant [J/mol-K]' /
DATA water(61) / '      51  4      5  12      0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA water(62) / ' 0.12533547935523d-1   -0.5     1.   0 !a(i),t(i),d(i),l(i)' /
DATA water(63) / ' 0.78957634722828d+1    0.875   1.   0' /
DATA water(64) / '-0.87803203303561d+1    1.0     1.   0' /
DATA water(65) / ' 0.31802509345418       0.5     2.   0' /
DATA water(66) / '-0.26145533859358       0.75    2.   0' /
DATA water(67) / '-0.78199751687981d-2    0.375   3.   0' /
DATA water(68) / ' 0.88089493102134d-2    1.0     4.   0' /
DATA water(69) / '-0.66856572307965       4.0     1.   1' /
DATA water(70) / ' 0.20433810950965       6.0     1.   1' /
DATA water(71) / '-0.66212605039687d-4   12.0     1.   1' /
DATA water(72) / '-0.19232721156002       1.0     2.   1' /
DATA water(73) / '-0.25709043003438       5.0     2.   1' /
DATA water(74) / ' 0.16074868486251       4.0     3.   1' /
DATA water(75) / '-0.40092828925807d-1    2.0     4.   1' /
DATA water(76) / ' 0.39343422603254d-6   13.0     4.   1' /
DATA water(77) / '-0.75941377088144d-5    9.0     5.   1' /
DATA water(78) / ' 0.56250979351888d-3    3.0     7.   1' /
DATA water(79) / '-0.15608652257135d-4    4.0     9.   1' /
DATA water(80) / ' 0.11537996422951d-8   11.0    10.   1' /
DATA water(81) / ' 0.36582165144204d-6    4.0    11.   1' /
DATA water(82) / '-0.13251180074668d-11  13.0    13.   1' /
DATA water(83) / '-0.62639586912454d-9    1.0    15.   1' /
DATA water(84) / '-0.10793600908932       7.0     1.   2' /
DATA water(85) / ' 0.17611491008752d-1    1.0     2.   2' /
DATA water(86) / ' 0.22132295167546       9.0     2.   2' /
DATA water(87) / '-0.40247669763528      10.0     2.   2' /
DATA water(88) / ' 0.58083399985759      10.0     3.   2' /
DATA water(89) / ' 0.49969146990806d-2    3.0     4.   2' /
DATA water(90) / '-0.31358700712549d-1    7.0     4.   2' /
DATA water(91) / '-0.74315929710341      10.0     4.   2' /
DATA water(92) / ' 0.47807329915480      10.0     5.   2' /
DATA water(93) / ' 0.20527940895948d-1    6.0     6.   2' /
DATA water(94) / '-0.13636435110343      10.0     6.   2' /
DATA water(95) / ' 0.14180634400617d-1   10.0     7.   2' /
DATA water(96) / ' 0.83326504880713d-2    1.0     9.   2' /
DATA water(97) / '-0.29052336009585d-1    2.0     9.   2' /
DATA water(98) / ' 0.38615085574206d-1    3.0     9.   2' /
DATA water(99) / '-0.20393486513704d-1    4.0     9.   2' /
DATA water(100) / '-0.16554050063734d-2    8.0     9.   2' /
DATA water(101) / ' 0.19955571979541d-2    6.0    10.   2' /
DATA water(102) / ' 0.15870308324157d-3    9.0    10.   2' /
DATA water(103) / '-0.16388568342530d-4    8.0    12.   2' /
DATA water(104) / ' 0.43613615723811d-1   16.0     3.   3' /
DATA water(105) / ' 0.34994005463765d-1   22.0     4.   3' /
DATA water(106) / '-0.76788197844621d-1   23.0     4.   3' /
DATA water(107) / ' 0.22446277332006d-1   23.0     5.   3' /
DATA water(108) / '-0.62689710414685d-4   10.0    14.   4' /
DATA water(109) / '-0.55711118565645d-9   50.0     3.   6' /
DATA water(110) / '-0.19905718354408      44.0     6.   6' /
DATA water(111) / ' 0.31777497330738      46.0     6.   6' /
DATA water(112) / '-0.11841182425981      50.0     6.   6' /
DATA water(113) / '-0.31306260323435d+2    0.0     3.   2 2    -20.  -150.  1.21  1.   0.   0.   0.' /
DATA water(114) / ' 0.31546140237781d+2    1.0     3.   2 2    -20.  -150.  1.21  1.   0.   0.   0.' /
DATA water(115) / '-0.25213154341695d+4    4.0     3.   2 2    -20.  -250.  1.25  1.   0.   0.   0.' /
DATA water(116) / '-0.14874640856724       0.0     1.   2 2    0.85   0.3   0.32  28.  700. 0.2  3.5' /
DATA water(117) / ' 0.31806110878444       0.0     1.   2 2    0.95   0.3   0.32  32.  800. 0.2  3.5' /
DATA water(118) / '' /
DATA water(119) / '' /
DATA water(120) / '#AUX               !auxiliary model specification' /
DATA water(121) / 'CPP  ideal gas heat capacity function' /
DATA water(122) / '?LITERATURE REFERENCE \' /
DATA water(123) / '?Wagner, W. and Pruss, A.,' /
DATA water(124) / '? "The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary' /
DATA water(125) / '? Water Substance for General and Scientific Use,"' /
DATA water(126) / '? J. Phys. Chem. Ref. Data, 31(2):387-535, 2002.' /
DATA water(127) / '?\' /
DATA water(128) / '!end of info section' /
DATA water(129) / '273.16             !lower temperature limit [K]' /
DATA water(130) / '2000.0             !upper temperature limit [K]' /
DATA water(131) / '0.0                !upper pressure limit [kPa]' /
DATA water(132) / '0.0                !maximum density [mol/L]' /
DATA water(133) / '1.0          8.314371357587            !reducing parameters for T, Cp0' /
DATA water(134) / '  1  5    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA water(135) / '0.40063200d+01      0.00' /
DATA water(136) / '0.12436000d-01    833.00' /
DATA water(137) / '0.97315000d+00   2289.0' /
DATA water(138) / '0.12795000d+01   5009.00' /
DATA water(139) / '0.96956000d+00   5982.0' /
DATA water(140) / '0.24873000d+00  17800.0' /
DATA water(141) / '' /
DATA water(142) / '' /
DATA water(143) / '@EOS               !equation of state specification' /
DATA water(144) / 'FEK  Helmholtz equation of state for water of Kunz and Wagner (2004).' /
DATA water(145) / '?LITERATURE REFERENCE \' /
DATA water(146) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA water(147) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA water(148) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA water(149) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA water(150) / '?\' /
DATA water(151) / '!end of info section' /
DATA water(152) / '273.16             !lower temperature limit [K]' /
DATA water(153) / '1350.0             !upper temperature limit [K]' /
DATA water(154) / '1000000.0          !upper pressure limit [kPa]' /
DATA water(155) / '73.96              !maximum density [mol/L]' /
DATA water(156) / 'PHK                                    !pointer to Cp0 model' /
DATA water(157) / '18.01528                               !molecular weight [g/mol]' /
DATA water(158) / '273.16                                 !triple point temperature [K]' /
DATA water(159) / '0.61248                                !pressure at triple point [kPa]' /
DATA water(160) / '55.49696                               !density at triple point [mol/L]' /
DATA water(161) / '373.17                                 !normal boiling point temperature [K]' /
DATA water(162) / ' 0.345                                 !acentric factor' /
DATA water(163) / '647.096     22064.0      17.87371609   !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA water(164) / '647.096                  17.87371609   !reducing parameters [K, mol/L]' /
DATA water(165) / '8.314472                               !gas constant [J/mol-K]' /
DATA water(166) / '  16  4      0  0      0  0            !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA water(167) / ' 0.82728408749586       0.5    1.  0' /
DATA water(168) / '-0.18602220416584d1     1.25   1.  0' /
DATA water(169) / '-0.11199009613744d1     1.875  1.  0' /
DATA water(170) / ' 0.15635753976056       0.125  2.  0' /
DATA water(171) / ' 0.87375844859025       1.5    2.  0' /
DATA water(172) / '-0.36674403715731       1.0    3.  0' /
DATA water(173) / ' 0.53987893432436d-1    0.75   4.  0' /
DATA water(174) / ' 0.10957690214499d1     1.5    1.  1' /
DATA water(175) / ' 0.53213037828563d-1    0.625  5.  1' /
DATA water(176) / ' 0.13050533930825d-1    2.625  5.  1' /
DATA water(177) / '-0.41079520434476       5.0    1.  2' /
DATA water(178) / ' 0.14637443344120       4.0    2.  2' /
DATA water(179) / '-0.55726838623719d-1    4.5    4.  2' /
DATA water(180) / '-0.11201774143800d-1    3.0    4.  3' /
DATA water(181) / '-0.66062758068099d-2    4.0    1.  5' /
DATA water(182) / ' 0.46918522004538d-2    6.0    1.  5' /
DATA water(183) / '' /
DATA water(184) / '' /
DATA water(185) / '#AUX               !auxiliary model specification' /
DATA water(186) / 'PHK  Helmholtz form for the ideal-gas state for water of Kunz and Wagner (2004).' /
DATA water(187) / '?LITERATURE REFERENCE \' /
DATA water(188) / '?Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M.' /
DATA water(189) / '? "The GERG-2004 Wide-Range Equation of State for Natural Gases' /
DATA water(190) / '? and Other Mixtures," GERG Technical Monograph 15,' /
DATA water(191) / '? Fortschritt-Berichte VDI, VDI-Verlag, Duesseldorf, 2007.' /
DATA water(192) / '?\' /
DATA water(193) / '!end of info section' /
DATA water(194) / '0.                 !lower temperature limit [K]' /
DATA water(195) / '1000.0             !upper temperature limit [K]' /
DATA water(196) / '0.0                !upper pressure limit [kPa]' /
DATA water(197) / '0.0                !maximum density [mol/L]' /
DATA water(198) / '1 2  0  1 2  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau)); cosh; sinh' /
DATA water(199) / '    3.00392      1.             !ai, ti for [ai*log(tau**ti)] terms' /
DATA water(200) / '    8.20352069   0.             !aj, ti for [ai*tau**ti] terms' /
DATA water(201) / '  -11.996306443  1.' /
DATA water(202) / '   -0.98763      1.763895929    !aj, ti for cosh and sinh terms' /
DATA water(203) / '    0.01059      0.415386589' /
DATA water(204) / '    3.06904      3.874803739' /
DATA water(205) / '' /
DATA water(206) / '' /
DATA water(207) / '@AUX               !auxiliary model specification' /
DATA water(208) / 'PH0  Helmholtz form for the ideal-gas state' /
DATA water(209) / '?LITERATURE REFERENCE \' /
DATA water(210) / '?Wagner, W. and Pruss, A.,' /
DATA water(211) / '? "The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary' /
DATA water(212) / '? Water Substance for General and Scientific Use,"' /
DATA water(213) / '? J. Phys. Chem. Ref. Data, 31(2):387-535, 2002.' /
DATA water(214) / '?\' /
DATA water(215) / '!end of info section' /
DATA water(216) / '273.16             !lower temperature limit [K]' /
DATA water(217) / '1350.0             !upper temperature limit [K]' /
DATA water(218) / '0.0                !upper pressure limit [kPa]' /
DATA water(219) / '0.0                !maximum density [mol/L]' /
DATA water(220) / '1 2  5  0 0  0 0 0 !Nterms:  ai*log(tau**ti); ai*tau**ti; ai*log(1-exp(bi*tau))' /
DATA water(221) / ' 0.30063200d+01      1.0           !ai, ti for [ai*log(tau**ti)] terms' /
DATA water(222) / '-0.832044648201d+01  0.0           !aj, ti for [ai*tau**ti] terms' /
DATA water(223) / ' 0.66832105268d+01   1.0' /
DATA water(224) / ' 0.12436000d-01  -1.28728967       !aj, ti for [ai*log(1-exp(ti*tau)] terms' /
DATA water(225) / ' 0.97315000d+00  -3.53734222' /
DATA water(226) / ' 0.12795000d+01  -7.74073708' /
DATA water(227) / ' 0.96956000d+00  -9.24437796' /
DATA water(228) / ' 0.24873000d+00  -27.5075105' /
DATA water(229) / '' /
DATA water(230) / '' /
DATA water(231) / '@EOS               !equation of state specification' /
DATA water(232) / 'FE1  Helmholtz equation of state for water of Saul and Wagner (1989).' /
DATA water(233) / '?LITERATURE REFERENCE \' /
DATA water(234) / '?Saul, A. and Wagner, W.,' /
DATA water(235) / '? "A Fundamental Equation for Water Covering the Range From the' /
DATA water(236) / '? Melting Line to 1273 K at Pressures up to 25000 MPa,"' /
DATA water(237) / '? J. Phys. Chem. Ref. Data, 18(4):1537-1564, 1989.' /
DATA water(238) / '?\' /
DATA water(239) / '!end of info section' /
DATA water(240) / '273.16             !lower temperature limit [K]' /
DATA water(241) / '1273.0             !upper temperature limit [K]' /
DATA water(242) / '400000.0           !upper pressure limit [kPa]' /
DATA water(243) / '55.49              !maximum density [mol/L]' /
DATA water(244) / 'CP1                                    !pointer to Cp0 model' /
DATA water(245) / '18.01534                               !molecular weight [g/mol]' /
DATA water(246) / '273.16                                 !triple point temperature [K]' /
DATA water(247) / '0.61166                                !pressure at triple point [kPa]' /
DATA water(248) / '55.497                                 !density at triple point [mol/L]' /
DATA water(249) / '373.15                                 !normal boiling point temperature [K]' /
DATA water(250) / '0.341                                  !acentric factor' /
DATA water(251) / '647.14       22064.0      17.8737      !Tc [K], pc [kPa], rhoc [mol/L]' /
DATA water(252) / '647.14                    17.8737      !reducing parameters [K, mol/L]' /
DATA water(253) / '8.31434                                !gas constant [J/mol-K]' /
DATA water(254) / '      38  4      0  0       0  0       !# terms, # coeff/term for:  "normal" terms, critical, spare' /
DATA water(255) / ' 0.233000901300d+0  0.0     1.0     0  !a(i),t(i),d(i),l(i)' /
DATA water(256) / '-0.140209112800d+1  2.0     1.0     0' /
DATA water(257) / ' 0.117224804100d+0  0.0     2.0     0' /
DATA water(258) / '-0.185074949900d+0  1.0     2.0     0' /
DATA water(259) / ' 0.177011042200d+0  2.0     2.0     0' /
DATA water(260) / ' 0.552515179400d-1  3.0     2.0     0' /
DATA water(261) / '-0.341325738000d-3  5.0     3.0     0' /
DATA water(262) / ' 0.855727436700d-3  0.0     5.0     0' /
DATA water(263) / ' 0.371690068500d-3  1.0     5.0     0' /
DATA water(264) / '-0.130887123300d-3  3.0     6.0     0' /
DATA water(265) / ' 0.321689519900d-4  2.0     7.0     0' /
DATA water(266) / ' 0.278588103400d-6  5.0     8.0     0' /
DATA water(267) / '-0.352151113000d+0  5.0     1.0     2' /
DATA water(268) / ' 0.788191453600d-1  7.0     1.0     2' /
DATA water(269) / '-0.151966661000d-1  9.0     1.0     2' /
DATA water(270) / '-0.106845858600d+0  5.0     2.0     2' /
DATA water(271) / '-0.205504628800d+0  4.0     3.0     2' /
DATA water(272) / ' 0.914619801200d+0  6.0     3.0     2' /
DATA water(273) / ' 0.321334356900d-3 13.0     3.0     2' /
DATA water(274) / '-0.113359139100d+1  5.0     4.0     2' /
DATA water(275) / '-0.310752074900d+0  2.0     5.0     2' /
DATA water(276) / ' 0.121790152700d+1  3.0     5.0     2' /
DATA water(277) / '-0.448171083100d+0  2.0     6.0     2' /
DATA water(278) / ' 0.549421877200d-1  0.0     7.0     2' /
DATA water(279) / '-0.866522209600d-4 11.0     7.0     2' /
DATA water(280) / ' 0.384408408800d-1  1.0     8.0     2' /
DATA water(281) / ' 0.985304488400d-2  4.0     8.0     2' /
DATA water(282) / '-0.176759847200d-1  0.0     9.0     2' /
DATA water(283) / ' 0.148854922200d-2  0.0    11.0     2' /
DATA water(284) / '-0.307071906900d-2  3.0    11.0     2' /
DATA water(285) / ' 0.388080328000d-2  5.0    11.0     2' /
DATA water(286) / '-0.262750521500d-2  6.0    11.0     2' /
DATA water(287) / ' 0.525837138800d-3  7.0    11.0     2' /
DATA water(288) / '-0.171639690100d+0 13.0     2.0     3' /
DATA water(289) / ' 0.718882362400d-1 14.0     2.0     3' /
DATA water(290) / ' 0.588126835700d-1 15.0     3.0     3' /
DATA water(291) / '-0.145593888000d-1 24.0     3.0     3' /
DATA water(292) / '-0.121613940000d-1 15.0     5.0     3' /
DATA water(293) / '' /
DATA water(294) / '' /
DATA water(295) / '#AUX               !auxiliary model specification' /
DATA water(296) / 'CP1  ideal gas heat capacity function' /
DATA water(297) / '?LITERATURE REFERENCE \' /
DATA water(298) / '?Saul, A. and Wagner, W.,' /
DATA water(299) / '? "A Fundamental Equation for Water Covering the Range From the' /
DATA water(300) / '? Melting Line to 1273 K at Pressures up to 25000 MPa,"' /
DATA water(301) / '? J. Phys. Chem. Ref. Data, 18(4):1537-1564, 1989.' /
DATA water(302) / '?\' /
DATA water(303) / '!end of info section' /
DATA water(304) / '273.16             !lower temperature limit [K]' /
DATA water(305) / '1273.0             !upper temperature limit [K]' /
DATA water(306) / '0.0                !upper pressure limit [kPa]' /
DATA water(307) / '0.0                !maximum density [mol/L]' /
DATA water(308) / '1.0          8.31434                   !reducing parameters for T, Cp0' /
DATA water(309) / '  1  5    0  0    0  0  0              !Nterms:  polynomial, exponential, cosh, sinh' /
DATA water(310) / ' 0.400632d+1    0.00' /
DATA water(311) / ' 0.124360d-1    833.0' /
DATA water(312) / ' 0.973150d+0   2289.0' /
DATA water(313) / ' 0.127950d+1   5009.0' /
DATA water(314) / ' 0.969560d+0   5982.0' /
DATA water(315) / ' 0.248730d+0  17800.0' /
DATA water(316) / '' /
DATA water(317) / '' /
DATA water(318) / '' /
DATA water(319) / '#TCX               !thermal conductivity model specification' /
DATA water(320) / 'TC0  pure fluid thermal conductivity model of Huber et al. (2011).' /
DATA water(321) / '?LITERATURE REFERENCE \' /
DATA water(322) / '?International Association for the Properties of Water and Steam,' /
DATA water(323) / '? "Release on the IAPWS Formulation 2011 for the thermal conductivity of Ordinary Water Substance"' /
DATA water(324) / '? Sept. 2011, Plzen, Czech Republic.' /
DATA water(325) / '?\' /
DATA water(326) / '? M. L. Huber, R. A. Perkins, D. G. Friend, J. V. Sengers' /
DATA water(327) / '? M. J. Assael, I. N. Metaxa, K. Miyagawa, R. Hellmann and E. Vogel,' /
DATA water(328) / '? "New International Formulation for the Thermal Conductivity of H2O",' /
DATA water(329) / '? J. Phys. Chem. Ref. Data Vol.41, No.3 (2012) pp. 1-23. [http://dx.doi.org/10.1063/1.4738955]' /
DATA water(330) / '?' /
DATA water(331) / '?For the uncertainties, see the IAPWS Release or publication cited above.' /
DATA water(332) / '?\' /
DATA water(333) / '!end of info section' /
DATA water(334) / '251.165            !lower temperature limit [K]' /
DATA water(335) / '1350.0             !upper temperature limit [K]' /
DATA water(336) / '1000000.0          !upper pressure limit [kPa]' /
DATA water(337) / '73.96              !maximum density [mol/L]' /
DATA water(338) / 'H2O                !pointer to hardcoded thermal conductivity model' /
DATA water(339) / '5 0 30 0 0 0 0 0   !number of terms for various pieces' /
DATA water(340) / '647.096  17.8737279956   1.d-3    !reducing parameters for T, rho, tcx' /
DATA water(341) / ' 2.443221d-3     0.  0. 0. 0       !coeff, power in T' /
DATA water(342) / ' 1.323095d-2     1.  0. 0. 0' /
DATA water(343) / ' 6.770357d-3     2.  0. 0. 0' /
DATA water(344) / '-3.454586d-3     3.  0. 0. 0' /
DATA water(345) / ' 4.096266d-4     4.  0. 0. 0' /
DATA water(346) / ' 1.60397357      0.0   0.0   0.0   0' /
DATA water(347) / ' 2.33771842      1.0   0.0   0.0   0' /
DATA water(348) / ' 2.19650529      2.0   0.0   0.0   0' /
DATA water(349) / '-1.21051378      3.0   0.0   0.0   0' /
DATA water(350) / '-2.72033700      4.0   0.0   0.0   0' /
DATA water(351) / '-0.646013523     0.0   1.0   0.0   0' /
DATA water(352) / '-2.78843778      1.0   1.0   0.0   0' /
DATA water(353) / '-4.54580785      2.0   1.0   0.0   0' /
DATA water(354) / ' 1.60812989      3.0   1.0   0.0   0' /
DATA water(355) / ' 4.57586331      4.0   1.0   0.0   0' /
DATA water(356) / ' 0.111443906     0.0   2.0   0.0   0' /
DATA water(357) / ' 1.53616167      1.0   2.0   0.0   0' /
DATA water(358) / ' 3.55777244      2.0   2.0   0.0   0' /
DATA water(359) / '-0.621178141     3.0   2.0   0.0   0' /
DATA water(360) / '-3.18369245      4.0   2.0   0.0   0' /
DATA water(361) / ' 0.102997357     0.0   3.0   0.0   0' /
DATA water(362) / '-0.463045512     1.0   3.0   0.0   0' /
DATA water(363) / '-1.40944978      2.0   3.0   0.0   0' /
DATA water(364) / ' 0.0716373224    3.0   3.0   0.0   0' /
DATA water(365) / ' 1.11683480      4.0   3.0   0.0   0' /
DATA water(366) / '-0.0504123634    0.0   4.0   0.0   0' /
DATA water(367) / ' 0.0832827019    1.0   4.0   0.0   0' /
DATA water(368) / ' 0.275418278     2.0   4.0   0.0   0' /
DATA water(369) / ' 0.00000000      3.0   4.0   0.0   0' /
DATA water(370) / '-0.19268305      4.0   4.0   0.0   0' /
DATA water(371) / ' 0.00609859258   0.0   5.0   0.0   0' /
DATA water(372) / '-0.00719201245   1.0   5.0   0.0   0' /
DATA water(373) / '-0.0205938816    2.0   5.0   0.0   0' /
DATA water(374) / ' 0.00000000      3.0   5.0   0.0   0' /
DATA water(375) / ' 0.012913842     4.0   5.0   0.0   0' /
DATA water(376) / 'TK3                !pointer to critical enhancement auxiliary function' /
DATA water(377) / '' /
DATA water(378) / '' /
DATA water(379) / '#AUX               !thermal conductivity critical enhancement model' /
DATA water(380) / 'TK3  simplified thermal conductivity critical enhancement of Olchowy & Sengers' /
DATA water(381) / '?LITERATURE REFERENCE \' /
DATA water(382) / '? M. L. Huber, R. A. Perkins, D. G. Friend, J. V. Sengers' /
DATA water(383) / '? M. J. Assael, I. N. Metaxa, K. Miyagawa, R. Hellmann and E. Vogel,' /
DATA water(384) / '? "New International Formulation for the Thermal Conductivity of H2O",' /
DATA water(385) / '? J. Phys. Chem. Ref. Data Vol.41, No.3 (2012) pp. 1-23. [http://dx.doi.org/10.1063/1.4738955]' /
DATA water(386) / '?\' /
DATA water(387) / '!end of info section' /
DATA water(388) / '251.165            !lower temperature limit [K]' /
DATA water(389) / '20000.00           !upper temperature limit [K]' /
DATA water(390) / '1000000.0          !upper pressure limit [kPa]' /
DATA water(391) / '73.96              !maximum density [mol/L]' /
DATA water(392) / '9  0  0  0         !# terms:  critical-terms, spare, spare, spare' /
DATA water(393) / '1.0     1.0     1.0     !reducing par for T, rho, tcx' /
DATA water(394) / '0.630d0              !gnu (universal exponent)' /
DATA water(395) / '1.239d0              !gamma (universal exponent)' /
DATA water(396) / '1.01d0               !R0 (universal amplitude)was 1.03 and 1.05 and 1.01' /
DATA water(397) / '0.068d0              !z (universal exponent--not used for t.c., only viscosity)was0.069' /
DATA water(398) / '1.00d0               !c (constant in viscosity eqn = 1/[2 - (alpha + gamma)/(2*nu)], but often set to 1)' /
DATA water(399) / '1.3d-10              !xi0 (amplitude) [m]was 1.094' /
DATA water(400) / '0.060d0              !gam0 (amplitude) [-] was 0.0496' /
DATA water(401) / '0.40D-09             !qd_inverse (modified effective cutoff parameter) [m] may 22 2010 new value from Jan' /
DATA water(402) / '970.644d0            !tref (reference temperature) [= 1.5 * 647.096  K]' /
DATA water(403) / '' /
DATA water(404) / '' /
DATA water(405) / '#ETA               !viscosity model specification' /
DATA water(406) / 'VS0  pure fluid viscosity model of Huber et al. (2009).' /
DATA water(407) / '?LITERATURE REFERENCE \' /
DATA water(408) / '?International Association for the Properties of Water and Steam,' /
DATA water(409) / '? "Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance"' /
DATA water(410) / '? Sept. 2008, Berlin.' /
DATA water(411) / '?\' /
DATA water(412) / '? "New International Formulation for the viscosity of water"' /
DATA water(413) / '?  Huber, M.L., Perkins, R.A., Laesecke, A., Friend, D.G., Sengers, J.V., Assael, M.J.,' /
DATA water(414) / '?  Metaxa, I.M., Vogel, E., Mares, R. and Miyagawa, K.' /
DATA water(415) / '?  J. Phys. Chem. Ref. Data, Vol. 38, No. 2 (2009) pp. 101-125.' /
DATA water(416) / '?\' /
DATA water(417) / '?For the uncertainties, see the IAPWS Release or the publication cited above.' /
DATA water(418) / '? NOTE: To use in faster industrial mode, change critical model at end of this VS0 block' /
DATA water(419) / '? to NUL instead of I08' /
DATA water(420) / '?\' /
DATA water(421) / '!end of info section' /
DATA water(422) / '251.165            !lower temperature limit [K]' /
DATA water(423) / '1350.0             !upper temperature limit [K]' /
DATA water(424) / '1000000.0          !upper pressure limit [kPa]' /
DATA water(425) / '73.96              !maximum density [mol/L]' /
DATA water(426) / 'H2O                !pointer to hardcoded thermal conductivity model' /
DATA water(427) / '0 0 4 21 0 0 0 0   !number of terms for various pieces' /
DATA water(428) / '647.096d0  17.8737279956d0   100.0d0    !reducing parameters for T, rho, eta' /
DATA water(429) / ' 1.67752d0      0.0 0. 0. 0       !coeff, power in T' /
DATA water(430) / ' 2.20462d0      1.0 0. 0. 0' /
DATA water(431) / ' 0.6366564d0    2.0 0. 0. 0' /
DATA water(432) / '-0.241605d0     3.0 0. 0. 0' /
DATA water(433) / '  0.520094d+00 0.0 0.0 0.0 0' /
DATA water(434) / '  0.850895d-01 1.0 0.0 0.0 0' /
DATA water(435) / '  0.222531d+00 0.0 1.0 0.0 0' /
DATA water(436) / '  0.999115d+00 1.0 1.0 0.0 0' /
DATA water(437) / '  0.188797d+01 2.0 1.0 0.0 0' /
DATA water(438) / '  0.126613d+01 3.0 1.0 0.0 0' /
DATA water(439) / ' -0.281378d+00 0.0 2.0 0.0 0' /
DATA water(440) / ' -0.906851d+00 1.0 2.0 0.0 0' /
DATA water(441) / ' -0.772479d+00 2.0 2.0 0.0 0' /
DATA water(442) / '  0.161913d+00 0.0 3.0 0.0 0' /
DATA water(443) / '  0.257399d+00 1.0 3.0 0.0 0' /
DATA water(444) / ' -0.325372d-01 0.0 4.0 0.0 0' /
DATA water(445) / '  0.698452d-01 3.0 4.0 0.0 0' /
DATA water(446) / ' -0.435673d-02 3.0 6.0 0.0 0' /
DATA water(447) / ' -0.108374d+01 2.0 0.0 0.0 0' /
DATA water(448) / ' -0.289555d+00 3.0 0.0 0.0 0' /
DATA water(449) / '  0.120573d+00 5.0 1.0 0.0 0' /
DATA water(450) / ' -0.489837d+00 3.0 2.0 0.0 0' /
DATA water(451) / ' -0.257040d+00 4.0 2.0 0.0 0' /
DATA water(452) / '  0.872102d-02 4.0 5.0 0.0 0' /
DATA water(453) / ' -0.593264d-03 5.0 6.0 0.0 0' /
DATA water(454) / 'I08                !pointer to critical enhancement auxiliary function 2008 IAPWS formulation' /
DATA water(455) / '' /
DATA water(456) / '' /
DATA water(457) / '@ETA               !viscosity model specification' /
DATA water(458) / 'VS4  pure fluid generalized friction theory viscosity model of Quinones-Cisneros and Deiters (2006).' /
DATA water(459) / '?LITERATURE REFERENCE \' /
DATA water(460) / '? Quinones-Cisneros, S.E. and Deiters, U.K.' /
DATA water(461) / '? "Generalization of the Friction Theory for Viscosity Modeling,"' /
DATA water(462) / '? J. Phys. Chem. B, 110:12820-12834, 2006.' /
DATA water(463) / '?' /
DATA water(464) / '!end of info section' /
DATA water(465) / '273.16             !lower temperature limit [K]' /
DATA water(466) / '2000.0             !upper temperature limit [K]' /
DATA water(467) / '1000000.0          !upper pressure limit [kPa]' /
DATA water(468) / '73.96              !maximum density [mol/L]' /
DATA water(469) / '5 0 0 0 0 0        !number of terms associated with dilute-gas function' /
DATA water(470) / 'NUL                !pointer to reduced effective collision cross-section model;not used' /
DATA water(471) / '0.2641             !Lennard-Jones coefficient sigma [nm] for ECS method (not used)' /
DATA water(472) / '809.1              !Lennard-Jones coefficient epsilon/kappa [K] for ECS method (not used)' /
DATA water(473) / '647.096d0    1.0d0 !reducing parameters for T, eta' /
DATA water(474) / ' 0.0d0      0.5d0  !Chapman-Enskog term; not used here' /
DATA water(475) / ' 151.138d0  0.00d0 !empirical terms for eta0' /
DATA water(476) / '-444.318d0  0.25d0' /
DATA water(477) / ' 398.262d0  0.50d0' /
DATA water(478) / '-81.7008d0  0.75d0' /
DATA water(479) / '0                  !number of terms for initial density dependence; not yet used.' /
DATA water(480) / '-1.17407105202836d-05 -3.78854818708520d-07 3.56742875797909d-08 0.0d0 0.0d0 !a(0),a(1),a(2)' /
DATA water(481) / ' 1.62216397984014d-06 -8.36595322447571d-06 9.10862531286788d-08 0.0d0 0.0d0 !b(0),b(1),b(2)' /
DATA water(482) / ' 1.92706925578893d-05 -1.28679815491711d-05 0.00000000000000d+00 0.0d0 0.0d0 !c(0),c(1),c(2)' /
DATA water(483) / '-3.30144899918610d-10  0.00000000000000d+00 1.02931444103415d-11 0.0d0 0.0d0 !A(0),A(1),A(2)' /
DATA water(484) / ' 5.03139997945133d-10  1.82304182380560d-10 0.00000000000000d+00 0.0d0 0.0d0 !B(0),B(1),B(2)' /
DATA water(485) / ' 8.01449084635477d-10  5.65613687804585d-09 1.10163426018591d-10 0.0d0 0.0d0 !C(0),C(1),C(2)' /
DATA water(486) / ' 0.0d0                   0.0d0                  0.0d0            0.0d0 0.0d0 !D(0),D(1),D(2)' /
DATA water(487) / ' 0.0d0                   0.0d0                  0.0d0            0.0d0 0.0d0 !E(0),E(1),E(2)' /
DATA water(488) / 'NUL                !pointer to critical enhancement auxiliary function (none used)' /
DATA water(489) / '' /
DATA water(490) / '' /
DATA water(491) / '#STN        !surface tension specification' /
DATA water(492) / 'ST1  surface tension model from IAPWS.' /
DATA water(493) / '?LITERATURE REFERENCE \' /
DATA water(494) / '?International Association for the Properties of Water and Steam,' /
DATA water(495) / '? "Release on the surface tension of ordinary water substance,"' /
DATA water(496) / '? Physical Chemistry of Aqueous Systems:  Proceedings of the 12th' /
DATA water(497) / '? International Conference on the Properties of Water and Steam,' /
DATA water(498) / '? Orlando, Florida, September 11-16, A139-A142, 1994.' /
DATA water(499) / '?\' /
DATA water(500) / '?For the uncertainties in surface tension, see the IAPWS Release.' /
DATA water(501) / '?\' /
DATA water(502) / '!end of info section' /
DATA water(503) / '273.16             !lower temperature limit [K]' /
DATA water(504) / '647.096            !upper temperature limit [K]' /
DATA water(505) / '0.0                !(dummy) upper pressure limit' /
DATA water(506) / '0.0                !(dummy) maximum density' /
DATA water(507) / '2                           !number of terms in surface tension model' /
DATA water(508) / '647.096                     !critical temperature used in fit (dummy)' /
DATA water(509) / ' 0.2358      1.256          !sigma0 and n' /
DATA water(510) / '-0.147375    2.256' /
DATA water(511) / '' /
DATA water(512) / '' /
DATA water(513) / '#DE         !dielectric constant specification' /
DATA water(514) / 'DE2  dielectric constant model of Fernandez et al. (1997).' /
DATA water(515) / '?LITERATURE REFERENCE \' /
DATA water(516) / '?Fernandez, D.P., Goodwin, A.R.H., Lemmon, E.W., Levelt Sengers, J.M.H.,' /
DATA water(517) / '? and Williams, R.C.' /
DATA water(518) / '? "A Formulation for the Static Permittivity of Water and Steam at' /
DATA water(519) / '? Temperatures from 238 K to 873 K at Pressures up to 1200 MPa,' /
DATA water(520) / '? Including Derivatives and Debye-Huckel Coefficients,"' /
DATA water(521) / '? J. Phys. Chem. Ref. Data, 26(4):1125-1165, 1997.' /
DATA water(522) / '?\' /
DATA water(523) / '!end of info section' /
DATA water(524) / '273.16             !lower temperature limit [K]' /
DATA water(525) / '1350.0             !upper temperature limit [K]' /
DATA water(526) / '0.0                !(dummy) upper pressure limit' /
DATA water(527) / '0.0                !(dummy) maximum density' /
DATA water(528) / '647.096 17.8737279956 1.0    !reducing parameters for t, d, and p' /
DATA water(529) / '11 1 0 0 0 0                        !number of terms in dielectric constant model' /
DATA water(530) / ' 0.978224486826d+0  0.25d0 1.    0. !coef, t exp, d exp, p exp' /
DATA water(531) / '-0.957771379375d+0  1.     1.    0.' /
DATA water(532) / ' 0.237511794148d+0  2.5d0  1.    0.' /
DATA water(533) / ' 0.714692244396d+0  1.5d0  2.    0.' /
DATA water(534) / '-0.298217036956d+0  1.5d0  3.    0.' /
DATA water(535) / '-0.108863472196d+0  2.5d0  3.    0.' /
DATA water(536) / ' 0.949327488264d-1  2.     4.    0.' /
DATA water(537) / '-0.980469816509d-2  2.     5.    0.' /
DATA water(538) / ' 0.165167634970d-4  5.     6.    0.' /
DATA water(539) / ' 0.937359795772d-4  0.5d0  7.    0.' /
DATA water(540) / '-0.123179218720d-9  10.    10.   0.' /
DATA water(541) / ' 0.196096504426d-2  228.   1.    1.2d0' /
DATA water(542) / '' /
DATA water(543) / '' /
DATA water(544) / '#MLT        !melting line specification' /
DATA water(545) / 'MLH  melting line model of Wagner et al. (1994).' /
DATA water(546) / '?LITERATURE REFERENCE \' /
DATA water(547) / '?Wagner, W., Saul, A., and Pruss, A.,' /
DATA water(548) / '? "International Equations for the Pressure Along the Melting and Along' /
DATA water(549) / '? the Sublimation Curve of Ordinary Water Substance,"' /
DATA water(550) / '? J. Phys. Chem. Ref. Data, 23(3):515-527, 1994.' /
DATA water(551) / '?\' /
DATA water(552) / '!end of info section' /
DATA water(553) / '251.165            !lower temperature limit [K]' /
DATA water(554) / '1350.0             !upper temperature limit [K]' /
DATA water(555) / '0.0                !(dummy) upper pressure limit' /
DATA water(556) / '0.0                !(dummy) maximum density' /
DATA water(557) / '273.16      611.657       !reducing temperature and pressure' /
DATA water(558) / '3 1 1 1 3 0                 !number of terms in melting line equation' /
DATA water(559) / '0.119539337d7 0.30000d1     ! H2OIh, 251.165 - 273.16 K' /
DATA water(560) / '0.808183159d5 0.25750d2' /
DATA water(561) / '0.333826860d4 0.10375d3' /
DATA water(562) / '-0.299948 60                ! H2OIII, 251.165 - 256.164 K ' /
DATA water(563) / '-1.18721 8                  ! H2OV, 256.164 - 273.31 K' /
DATA water(564) / '-1.07476 4.6                ! H2OVI, 273.31 - 355 K' /
DATA water(565) / '0.173683d1 -1               ! H2OVII, 355 - 715 K' /
DATA water(566) / '-0.544606d-1 5' /
DATA water(567) / '0.806106d-7 22' /
DATA water(568) / '' /
DATA water(569) / '' /
DATA water(570) / '' /
DATA water(571) / '' /
DATA water(572) / '#SBL        !sublimation line specification' /
DATA water(573) / 'SB2  sublimation line model of Wagner et al. (2010).' /
DATA water(574) / '?LITERATURE REFERENCE \' /
DATA water(575) / '? Wagner, W., Riethmann, T., Feistel, R., and Harvey, A.H.,' /
DATA water(576) / '? "New Equations for the Sublimation Pressure and Melting Pressure of' /
DATA water(577) / '? H2O Ice Ih,"' /
DATA water(578) / '? J. Phys. Chem. Ref. Data, 40(4), 2011.' /
DATA water(579) / '?\' /
DATA water(580) / '!end of info section' /
DATA water(581) / '50.0               !lower temperature limit [K]' /
DATA water(582) / '273.16             !upper temperature limit [K]' /
DATA water(583) / '0.0                !(dummy) upper pressure limit' /
DATA water(584) / '0.0                !(dummy) maximum density' /
DATA water(585) / '273.16 0.611657    !reducing temperature and pressure' /
DATA water(586) / '3 0 0 0 0 0                 !number of terms in sublimation line equation' /
DATA water(587) / '-21.2144006 -0.99666666667  !coefficients and exponents' /
DATA water(588) / ' 27.3203819  0.20666667' /
DATA water(589) / '-6.10598130  0.70333333' /
DATA water(590) / '' /
DATA water(591) / '' /
DATA water(592) / '#PS         !vapor pressure equation' /
DATA water(593) / 'PS6  vapor pressure equation of Wagner and Pruss (2002).' /
DATA water(594) / '?LITERATURE REFERENCE \' /
DATA water(595) / '?See EOS' /
DATA water(596) / '?\' /
DATA water(597) / '!end of info section' /
DATA water(598) / '273.16             !lower temperature limit [K]' /
DATA water(599) / '647.096            !upper temperature limit [K]' /
DATA water(600) / '0.0                !(dummy) upper pressure limit' /
DATA water(601) / '0.0                !(dummy) maximum density' /
DATA water(602) / '647.096 22064.0    !reducing parameters' /
DATA water(603) / '6 0 0 0 0 0                 !number of terms in equation' /
DATA water(604) / ' -7.85951783         2.     !coefficients and exponents' /
DATA water(605) / '  1.84408259         3.' /
DATA water(606) / ' -11.7866497         6.' /
DATA water(607) / '  22.6807411         7.' /
DATA water(608) / ' -15.9618719         8.' /
DATA water(609) / '  1.80122502         15.' /
DATA water(610) / '' /
DATA water(611) / '' /
DATA water(612) / '#DL         !saturated liquid density equation' /
DATA water(613) / 'DL2  saturated liquid density equation of Wagner and Pruss (2002).' /
DATA water(614) / '?LITERATURE REFERENCE \' /
DATA water(615) / '?See EOS' /
DATA water(616) / '?\' /
DATA water(617) / '!end of info section' /
DATA water(618) / '273.16             !lower temperature limit [K]' /
DATA water(619) / '647.096            !upper temperature limit [K]' /
DATA water(620) / '0.0                !(dummy) upper pressure limit' /
DATA water(621) / '0.0                !(dummy) maximum density' /
DATA water(622) / '647.096 17.8737279956 !reducing parameters' /
DATA water(623) / '6 0 0 0 0 0                 !number of terms in equation' /
DATA water(624) / '  1.99274064         1.     !coefficients and exponents' /
DATA water(625) / '  1.09965342         2.' /
DATA water(626) / ' -0.510839303        5.' /
DATA water(627) / ' -1.75493479         16.' /
DATA water(628) / ' -45.5170352         43.' /
DATA water(629) / ' -6.74694450d5       110.' /
DATA water(630) / '' /
DATA water(631) / '' /
DATA water(632) / '#DV         !saturated vapor density equation' /
DATA water(633) / 'DV4  saturated vapor density equation of Wagner and Pruss (2002).' /
DATA water(634) / '?LITERATURE REFERENCE \' /
DATA water(635) / '?See EOS' /
DATA water(636) / '?\' /
DATA water(637) / '!end of info section' /
DATA water(638) / '273.16             !lower temperature limit [K]' /
DATA water(639) / '647.096            !upper temperature limit [K]' /
DATA water(640) / '0.0                !(dummy) upper pressure limit' /
DATA water(641) / '0.0                !(dummy) maximum density' /
DATA water(642) / '647.096 17.8737279956 !reducing parameters' /
DATA water(643) / '6 0 0 0 0 0                 !number of terms in equation' /
DATA water(644) / ' -2.03150240         1.0    !coefficients and exponents' /
DATA water(645) / ' -2.68302940         2.0' /
DATA water(646) / ' -5.38626492         4.0' /
DATA water(647) / ' -17.2991605         9.0' /
DATA water(648) / ' -44.7586581         18.5' /
DATA water(649) / ' -63.9201063         35.5' /
DATA water(650) / '' /
DATA water(651) / '' /
DATA water(652) / '@END' /
DATA water(653) / 'c        1         2         3         4         5         6         7         8' /
DATA water(654) / 'c2345678901234567890123456789012345678901234567890123456789012345678901234567890' /
contains
!----------------------------------------------------------------------------------
subroutine get_eos_cg_2019_pure_fld_files(ptr,name)
implicit none
character(256),pointer,dimension(:)::ptr
character(*),dimension(:) :: name
if(any(name.eq.'argon')) then
ptr=>argon
elseif(any(name.eq.'cas')) then
ptr=>f106_98_9____1butene
elseif(any(name.eq.'chlorine')) then
ptr=>Chlorine
elseif(any(name.eq.'co')) then
ptr=>carbon_monoxide
elseif(any(name.eq.'co2')) then
ptr=>carbon_dioxide
elseif(any(name.eq.'dea')) then
ptr=>Diethanolamine
elseif(any(name.eq.'h2s')) then
ptr=>Hydrogen_sulfide
elseif(any(name.eq.'hcl')) then
ptr=>Hydrogen_chloride
elseif(any(name.eq.'hydrogen')) then
ptr=>hydrogen_normal
elseif(any(name.eq.'mea')) then
ptr=>Monoethanolamine
elseif(any(name.eq.'methane')) then
ptr=>Methane
elseif(any(name.eq.'nitrogen')) then
ptr=>nitrogen
elseif(any(name.eq.'oxygen')) then
ptr=>oxygen
elseif(any(name.eq.'so2')) then
ptr=>Sulfur_dioxide
elseif(any(name.eq.'water')) then
ptr=>water
end if
end subroutine get_eos_cg_2019_pure_fld_files
end module eos_cg_2019_pure_fld_files
