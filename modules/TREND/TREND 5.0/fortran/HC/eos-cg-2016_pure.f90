! This file was generated -  Dont edit!
module eos_cg_2016_pure_fld_files
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
DATA f106_98_9____1butene(10) / '112-40-3    c12' /
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
subroutine get_eos_cg_2016_pure_fld_files(ptr,name)
implicit none
character(256),pointer,dimension(:)::ptr
character(*),dimension(:) :: name
if(any(name.eq.'argon')) then
ptr=>argon
elseif(any(name.eq.'cas')) then
ptr=>f106_98_9____1butene
elseif(any(name.eq.'co')) then
ptr=>carbon_monoxide
elseif(any(name.eq.'co2')) then
ptr=>carbon_dioxide
elseif(any(name.eq.'nitrogen')) then
ptr=>nitrogen
elseif(any(name.eq.'oxygen')) then
ptr=>oxygen
elseif(any(name.eq.'water')) then
ptr=>water
end if
end subroutine get_eos_cg_2016_pure_fld_files
end module eos_cg_2016_pure_fld_files
