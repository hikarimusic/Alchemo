data = '''
!
!V(angle) = Ktheta(Theta - Theta0)**2
!
!V(Urey-Bradley) = Kub(S - S0)**2
!
!Ktheta: kcal/mole/rad**2
!Theta0: degrees
!Kub: kcal/mole/A**2 (Urey-Bradley)
!S0: A
!
!atom types     Ktheta    Theta0   Kub     S0
!
H    NH2  CT1   50.000    111.00              ! From LSN HC-NH2-CT2
H    NH2  CT2   50.000    111.00              ! From LSN HC-NH2-CT2, Neutral Gly Nterminus
NH2  CT1  CT1   67.700    110.00              ! From LSN NH2-CT2-CT2
NH2  CT1  CT2   67.700    110.00              ! From LSN NH2-CT2-CT2
NH2  CT1  CT2A  67.700    110.00              ! From LSN NH2-CT2-CT2
NH2  CT1  CT3   67.700    110.00              ! From LSN NH2-CT2-CT2
CT1  CD   OH1   55.000    110.50              ! From ASPP CT2-CD-OH1
CT3  CT1  CD    52.000    108.00              ! Ala cter
NH2  CT1  HB1   38.000    109.50   50.00   2.1400 ! From LSN NH2-CT2-HA
NH2  CT1  C     50.000    107.00              ! From ALA Dipep. NH1-CT2-C
NH2  CT1  CC    50.000    107.00              ! From ALA Dipep. NH1-CT2-C, for nsam
NH2  CT1  CD    50.000    107.00              ! From ALA Dipep. NH1-CT2-C
NH2  CT2  C     50.000    107.00              ! From ALA Dipep. NH1-CT2-C, Neutral Gly Nterminus
HB2  CT1  HB2   36.000    115.00              ! from HB2  CT2  HB2
HB2  CT1  CD    50.000    109.50              ! from HB2  CT2  CD
NH1  CT1  HB2   48.000    108.00              ! from NH1  CT2  HB2

!
!Indole/Tryptophan
CAI  CAI  CA    40.000    120.00   35.00   2.41620 ! from CA CA CA
CAI  CA   CA    40.000    120.00   35.00   2.41620 ! from CA CA CA
CPT  CA   CA    50.000    113.20 ! atm, methylindole, 1/17/04
CPT  CPT  CA    50.000    110.00 ! atm, methylindole, 1/17/04
CPT  CAI  CA    50.000    113.20 ! atm, methylindole, 1/17/04
CPT  CPT  CAI   50.000    110.00 ! atm, methylindole, 1/17/04
CPT  CY   CA    85.000    106.40   25.00   2.26100 ! atm, methylindole, 1/17/04
CPT  NY   CA    85.000    112.00 ! atm, methylindole, 1/17/04
CT2  CY   CA    30.000    127.00 ! atm, methylindole, CT3  CY   CA
CT2  CY   CPT   30.000    126.70 ! atm, methylindole, 1/17/04
CT3  CY   CA    30.000    127.00 ! atm, methylindole, CT3  CY   CA
CT3  CY   CPT   30.000    126.70 ! atm, methylindole, 1/17/04
CY   CPT  CA   130.000    133.50 ! atm, methylindole, 1/17/04
CY   CPT  CAI  130.000    133.50 ! atm, methylindole, 1/17/04
CY   CPT  CPT   85.000    108.00 ! atm, methylindole, 1/17/04
CY   CT2  CT1   58.350    114.00 ! from TRP crystal, JWK
CY   CT2  CT3   58.350    114.00 ! from TRP crystal, JWK
H    NY   CA    28.000    126.00 ! trp, adm jr., 12/30/91
H    NY   CAI   28.000    126.00 ! trp, adm jr., 12/30/91
H    NY   CPT   28.000    126.00 ! trp, adm jr., 12/30/91
HA2  CT2  CY    55.000    109.50 ! atm, methylindole, 1/17/04
HA3  CT3  CY    55.000    109.50 ! atm, methylindole, 1/17/04
HP   CA   CAI   30.000    120.00   22.00   2.15250 ! from HP CA CA
HP   CAI  CA    30.000    120.00   22.00   2.15250 ! from HP CA CA
HP   CA   CPT   30.000    122.00   22.00   2.14600 ! trp, adm jr., 12/30/91
HP   CAI  CPT   30.000    122.00   22.00   2.14600 ! from HP CA CPT
HP   CA   CY    32.000    125.00   25.00   2.17300 ! JWK 05/14/91 new theta0 and r0UB from indole
HP   CY   CA    32.000    126.40   25.00   2.18600 ! trp, adm jr., 12/30/91
HP   CY   CPT   32.000    126.40   25.00   2.25500 ! JWK 05/14/91 new theta0 and r0UB from indole
NY   CA   CY    85.000    110.50   25.00   2.24000 ! trp, adm jr., 12/30/91
NY   CA   HP    32.000    125.00   25.00   2.17700 ! JWK 05/14/91 new theta0 and r0UB from indole
NY   CPT  CA   130.000    129.50 ! atm, methylindole, 1/17/04
NY   CPT  CAI  130.000    129.50 ! atm, methylindole, 1/17/04
NY   CPT  CPT   95.000    107.40 ! atm, methylindole, 1/17/04
CA   CA   CA    40.000    120.00   35.00   2.41620 ! ALLOW   ARO
                ! JES 8/25/89
CE1  CE1  CT2    48.00    123.50   !
                ! for 2-butene, yin/adm jr., 12/95
CE1  CE1  CT3    48.00    123.50   ! 
		! for 2-butene, yin/adm jr., 12/95
CE1  CT2  CT3    32.00    112.20   ! 
		! for 1-butene; from propene, yin/adm jr., 12/95
CE2  CE1  CT2    48.00    126.00   ! 
		! for 1-butene; from propene, yin/adm jr., 12/95
CE2  CE1  CT3    47.00    125.20   ! 
		! for propene, yin/adm jr., 12/95
CP1  N    C      60.000   117.0000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP2  CP1  C      52.000   112.3000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP2  CP1  CC     52.000   112.3000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP2  CP1  CD     50.000   112.3000 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP2  CP2  CP1    70.000   108.5000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  CP2  CP2    70.000   108.5000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  N    C      60.000   117.0000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  N    CP1   100.000   114.2000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  NP   CP1   100.000   111.0000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CPH2 NR1  CPH1  130.000   107.5000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
CPH2 NR2  CPH1  130.000   104.0000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
CPH2 NR3  CPH1  145.000   108.0000 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
CT1  CT1  C      52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! Alanine Dipeptide ab initio calc's (LK)
CT1  CT1  CC     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
CT1  CT1  CD     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 6/27/2012, for Thr with CT1 patch
CT1  CT1  CT1   53.350    111.00    8.00   2.56100 ! ALLOW ALI
                ! alkane update, adm jr., 3/2/92
CT1  CT2  CA     51.800   107.5000 ! ALLOW   ALI ARO
                ! PARALLH19 (JES)
CT1  CT2  CC     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
CT1  CT2  CD     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 5/02/91, acetic acid pure solvent
CT1  CT2  CPH1   58.350   113.0000 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, from CT2CT2CT, U-B omitted
CT1  CT2  CT1   58.350    113.50   11.16   2.56100 ! ALLOW   ALI
                ! alkane frequencies (MJF), alkane geometries (SF)
CT1  NH1  C      50.000   120.0000 ! ALLOW   ALI PEP POL ARO
                ! NMA Vib Modes (LK)
CT2  CA   CA     45.800   122.3000 ! ALLOW   ALI ARO
                ! PARALLH19 (JES)
CT2  CPH1 CPH1   45.800   130.0000 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FC=>CT2CA CA,BA=> CRYSTALS
CT2  CT1  C      52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! Alanine Dipeptide ab initio calc's (LK)
CT2  CT1  CC     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
CT2A CT1  CC     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
CT2  CT1  CD     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 5/02/91, acetic acid pure solvent
CT2  CT1  CT1   53.350    111.00    8.00   2.56100 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
CT2  CT2  C      52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! from CT2  CT1  C, for lactams, adm jr.
CT2  CT2  CC     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
CT3  CT2  CC     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
CT2  CT2  CD     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 5/02/91, acetic acid pure solvent
CT2A CT2  CD     52.000   108.0000 ! for GLUP, ZHU
CT2  CT2  CT1   58.350    113.50   11.16   2.56100 ! ALLOW   ALI
                ! alkane frequencies (MJF), alkane geometries (SF)
CT2  CT2  CT2   58.350    113.60   11.16   2.56100 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
CT2  CT3  CT1   58.350    113.50   11.16   2.56100 ! ALLOW   ALI
                ! alkane frequencies (MJF), alkane geometries (SF)
CT2  NC2  C      62.300   120.0000 ! ALLOW   ALI POL PEP ARO
                ! 107.5->120.0 to make planar Arg (KK)
CT2  NH1  C      50.000   120.0000 ! ALLOW   ALI PEP POL ARO
                ! NMA Vib Modes (LK)
CT2  OS   CD    40.000    109.60   30.00   2.26510 ! ALLOW  POL PEP
                ! adm jr. 5/02/91, acetic acid pure solvent
CT3  CA   CA     45.800   122.3000 ! ALLOW   ALI ARO
                ! toluene, adm jr., 3/7/92
CT3  CPH1 CPH1   45.800   130.0000 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FC=>CT2CA CA,BA=> CRYSTALS
CT3  CT1  C      52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! Alanine Dipeptide ab initio calc's (LK)
CT3  CT1  CC     52.000   108.0000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/09/92, for ALA cter
CT3  CT1  CT1   53.350    108.50    8.00   2.56100 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
CT3  CT1  CT2   53.350    114.00    8.00   2.56100 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
CT3  CT1  CT3   53.350    114.00    8.00   2.56100 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
CT3  CT2  CA     51.800   107.5000 ! ALLOW   ALI ARO
                ! ethylbenzene, adm jr., 3/7/92
CT3  CT2  CPH1   58.350   113.0000 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, from CT2CT2CT, U-B omitted
CT3  CT2  CT1   58.350    113.50   11.16   2.56100 ! ALLOW   ALI
                ! alkane frequencies (MJF), alkane geometries (SF)
CT3  CT2  CT2   58.000    115.00    8.00   2.56100 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
CT3  CT2  CT3   53.350    114.00    8.00   2.56100 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
CT3  NC2  C      62.300   120.0000 ! ALLOW   ALI POL PEP ARO
                ! methylguanidinium, adm jr., 3/26/92
CT3  NH1  C      50.000   120.0000 ! ALLOW   ALI PEP POL ARO
                ! NMA Vib Modes (LK)
CT3  OS   CD    40.000    109.60   30.00   2.26510 ! ALLOW  POL PEP
                ! adm jr. 5/02/91, acetic acid pure solvent
CT3  S    CT2    34.000    95.0000 ! ALLOW   ALI SUL ION
                ! expt. MeEtS,    3/26/92 (FL)
H    NH1  C      34.000   123.0000 ! ALLOW   PEP POL ARO
                ! NMA Vib Modes (LK)
H    NH1  CT1    35.000   117.0000 ! ALLOW   PEP POL ARO ALI
                ! NMA Vibrational Modes (LK)
H    NH1  CT2    35.000   117.0000 ! ALLOW   PEP POL ARO ALI
                ! NMA Vibrational Modes (LK)
H    NH1  CT3    35.000   117.0000 ! ALLOW   PEP POL ARO ALI
                ! NMA Vibrational Modes (LK)
H    NH2  CC     50.000   120.0000 ! ALLOW   POL PEP ARO
                ! his, adm jr. 8/13/90 acetamide geometry and vibrations
H    NH2  H      23.000   120.0000 ! ALLOW   POL
                ! adm jr. 8/13/90 acetamide geometry and vibrations
H    NR1  CPH1  30.000    125.50   20.00   2.15000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
H    NR1  CPH2  30.000    127.00   20.00   2.14000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
H    NR3  CPH1  25.000    126.00   15.00   2.13000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
H    NR3  CPH2  25.000    126.00   15.00   2.09000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
H    OH1  CA     65.000   108.0000 ! ALLOW   ALC ARO
                ! JES 8/25/89 phenol
H    OH1  CD     55.000   115.0000 ! ALLOW   ALC ARO PEP POL
                ! adm jr. 5/02/91, acetic acid pure solvent
H    OH1  CT1    57.500   106.0000 ! ALLOW   ALC ARO ALI
                ! methanol vib fit EMB 11/21/89
H    OH1  CT2    57.500   106.0000 ! ALLOW   ALC ARO ALI
                ! methanol vib fit EMB 11/21/89
H    OH1  CT3    57.500   106.0000 ! ALLOW   ALC ARO ALI
                ! methanol vib fit EMB 11/21/89
HA2  CP2  CP1   33.430    110.10   22.53   2.17900 ! ALLOW ALI PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HA2  CP2  CP2   26.500    110.10   22.53   2.17900 ! ALLOW ALI PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HA2  CP2  CP3   26.500    110.10   22.53   2.17900 ! ALLOW ALI PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HA2  CP2  HA2   35.500    109.00    5.40   1.80200 ! ALLOW ALI PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HA2  CP3  CP2   26.500    110.10   22.53   2.17900 ! ALLOW ALI PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HA2  CP3  HA2   35.500    109.00    5.40   1.80200 ! ALLOW ALI PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HA2  CS   CT3   34.600    110.10   22.53   2.17900 ! ALLOW SUL
                ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
HA2  CS   HA2   35.500    108.40   14.00   1.77500 ! ALLOW SUL
                ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
HA3  CS   HA3   35.500    108.40   14.00   1.77500 ! ALLOW SUL
                ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
HA1  CT1  C     33.000    109.50   30.00   2.16300 ! ALLOW   ALI PEP POL ARO
                ! alanine dipeptide, LK, replaced, adm jr., 5/09/91
HA1  CT1  CD    33.000    109.50   30.00   2.16300 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 5/02/91, acetic acid pure solvent
HA1  CT1  CT1   34.500    110.10   22.53   2.17900 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
HA1  CT1  CT2   34.500    110.10   22.53   2.17900 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
HA1  CT1  CT3   34.500    110.10   22.53   2.17900 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
HA1  CT1  HA1   35.500    109.00    5.40   1.80200 ! TEST for test cpd
                ! based on HA   CT2  HA
HA2  CT2  C     33.000    109.50   30.00   2.16300 ! ALLOW   ALI PEP POL ARO
                ! alanine dipeptide, LK, replaced, adm jr., 5/09/91
HA2  CT2  CA     49.300   107.5000 ! ALLOW   ALI ARO
                ! PARALLH19 (JES)
HA2  CT2  CC    33.000    109.50   30.00   2.16300 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
HA2  CT2  CD    33.000    109.50   30.00   2.16300 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 5/02/91, acetic acid pure solvent
HA2  CT2  CE1    45.00    111.50   ! 
		! for 1-butene; from propene, yin/adm jr., 12/95
HA2  CT2  CPH1   33.430   109.5000 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, from CT2CT2HA, U-B OMITTED
HA2  CT2  CT1   26.500    110.10   22.53   2.17900 ! ALLOW   ALI
                ! alkane frequencies (MJF), alkane geometries (SF)
HA2  CT2  CT2   26.500    110.10   22.53   2.17900 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
HA2  CT2  CT3   34.600    110.10   22.53   2.17900 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
HA2  CT2  HA2   35.500    109.00    5.40   1.80200 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
HA3  CT3  C     33.000    109.50   30.00   2.16300 ! ALLOW   ALI PEP POL ARO
                ! alanine dipeptide, LK, replaced, adm jr., 5/09/91
HA3  CT3  CA     49.300   107.5000 ! ALLOW   ALI ARO
                ! toluene, adm jr. 3/7/92
HA3  CT3  CC    33.000    109.50   30.00   2.16300 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
HA3  CT3  CD    33.000    109.50   30.00   2.16300 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 5/02/91, acetic acid pure solvent
HA3  CT3  CE1    42.00    111.50   ! 
		! for 2-butene, yin/adm jr., 12/95
HA3  CT3  CPH1   33.430   109.5000 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, from CT2CT2HA, U-B OMITTED
HA3  CT3  CS    34.600    110.10   22.53   2.17900 ! ALLOW SUL
                ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
HA3  CT3  CT1   33.430    110.10   22.53   2.17900 ! ALLOW   ALI
                ! alkane frequencies (MJF), alkane geometries (SF)
HA3  CT3  CT2   34.600    110.10   22.53   2.17900 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
HA3  CT3  CT3   37.500    110.10   22.53   2.17900 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
HA3  CT3  HA3   35.500    108.40    5.40   1.80200 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
HE1  CE1  CE1    52.00    119.50   ! 
		! for 2-butene, yin/adm jr., 12/95
HE1  CE1  CE2    42.00    118.00   ! 
		! for propene, yin/adm jr., 12/95
HE1  CE1  CT2    40.00    116.00   ! 
		! for 1-butene; from propene, yin/adm jr., 12/95
HE1  CE1  CT3    22.00    117.00   ! 
		! for propene, yin/adm jr., 12/95
HE2  CE2  CE1    45.00    120.50   ! 
		! for propene, yin/adm jr., 12/95
HE2  CE2  CE2    55.50    120.50   ! 
		! for ethene, yin/adm jr., 12/95
HE2  CE2  HE2    19.00    119.00   ! 
		! for propene, yin/adm jr., 12/95
HB1  CP1  C      50.000   112.0000 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HB1  CP1  CC     50.000   112.0000 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HB1  CP1  CD     50.000   112.0000 ! ALLOW PEP POL PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HB1  CP1  CP2    35.000   118.0000 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HB1  CT1  C      50.000   109.5000 ! ALLOW  PEP
                ! Alanine Dipeptide ab initio calc's (LK)
HB1  CT1  CC     50.000   109.5000 ! ALLOW  PEP POL
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
HB1  CT1  CD     50.000   109.5000 ! ALLOW  PEP POL
                ! adm jr. 5/02/91, acetic acid pure solvent
HB1  CT1  CT1    35.000   111.0000 ! ALLOW  PEP
                ! Alanine Dipeptide ab initio calc's (LK)
HB1  CT1  CT2    35.000   111.0000 ! ALLOW  PEP
                ! Alanine Dipeptide ab initio calc's (LK)
HB1  CT1  CT3    35.000   111.0000 ! ALLOW  PEP
                ! Alanine Dipeptide ab initio calc's (LK)
HB2  CT2  C      50.000   109.5000 ! ALLOW  PEP
                ! Alanine Dipeptide ab initio calc's (LK)
HB2  CT2  CC     50.000   109.5000 ! ALLOW  PEP POL
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
HB2  CT2  CD     50.000   109.5000 ! ALLOW  PEP POL
                ! adm jr. 5/02/91, acetic acid pure solvent
HB2  CT2  HB2    36.000   115.0000 ! ALLOW   PEP
                ! Alanine Dipeptide ab initio calc's (LK)
HC   NC2  C      49.000   120.0000 ! ALLOW   POL PEP ARO
                ! 35.3->49.0 GUANIDINIUM (KK)
HC   NC2  CT2    40.400   120.0000 ! ALLOW   POL ALI
                ! 107.5->120.0 to make planar Arg (KK)
HC   NC2  CT3    40.400   120.0000 ! ALLOW   POL ALI
                ! methylguanidinium, adm jr., 3/26/92
HC   NC2  HC     25.000   120.0000 ! ALLOW   POL
                ! 40.0->25.0 GUANIDINIUM (KK)
HC   NH2  CT2    50.000   111.0000 ! ALLOW   POL
                ! from HC NH2 CT3, neutral glycine, adm jr.
HC   NH2  CT3    50.000   111.0000 ! ALLOW   POL
                ! methylamine geom/freq, adm jr., 6/2/92
HC   NH2  HC     39.000   106.5000 ! ALLOW   POL
                ! 40.0->25.0 GUANIDINIUM (KK)
HC   NH3  CT1   30.000    109.50   20.00   2.07400 ! ALLOW   POL ALI
                ! new stretch and bend; methylammonium (KK 03/10/92)
HC   NH3  CT2   30.000    109.50   20.00   2.07400 ! ALLOW   POL ALI
                ! new stretch and bend; methylammonium (KK 03/10/92)
HC   NH3  CT3   30.000    109.50   20.00   2.07400 ! ALLOW   POL ALI
                ! new stretch and bend; methylammonium (KK 03/10/92)
HC   NH3  HC     44.000   109.5000 ! ALLOW   POL
                ! new stretch and bend; methylammonium (KK 03/10/92)
HC   NP   CP1   33.000    109.50    4.00   2.05600 ! ALLOW POL ALI PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HC   NP   CP3   33.000    109.50    4.00   2.05600 ! ALLOW POL ALI PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HC   NP   HC     51.000   107.5000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HP   CA   CA    30.000    120.00   22.00   2.15250 ! ALLOW   ARO
                ! JES 8/25/89 benzene
HR1  CPH1 CPH1  22.000    130.00   15.00   2.21500 ! ALLOW ARO
                ! adm jr., 6/27/90, his
HR3  CPH1 CPH1  25.000    130.00   20.00   2.20000 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HS   S    CT2    38.800    95.0000 ! ALLOW   SUL ION ALI
                ! methanethiol pure solvent, adm jr., 6/22/92
HS   S    CT3    43.000    95.0000 ! ALLOW   SUL ION ALI
                ! methanethiol pure solvent, adm jr., 6/22/92
N    C    CP1    20.000   112.5000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CT1    20.000   112.5000 ! ALLOW ALI PEP POL ARO PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CT2    20.000   112.5000 ! ALLOW ALI PEP POL ARO PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CT3    20.000   112.5000 ! ALLOW ALI PEP POL ARO PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    CP1  C      50.000   108.2000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    CP1  CC     50.000   108.2000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    CP1  CD     50.000   108.2000 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    CP1  CP2    70.000   110.8000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    CP1  HB1    48.000   112.0000 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    CP3  CP2    70.000   110.5000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    CP3  HA2    48.000   108.0000 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NC2  C    NC2   40.000    120.00   70.00   2.31
                ! mp2/6-311g** mgua vib data, adm jr., 1/04
                ! N-N distances: 2.29001, 2.31146, 2.33240
NC2  CT2  CT2    67.700   107.5000 ! ALLOW   ALI POL
                ! arg, (DS)
NC2  CT2  HA2    56.500   107.5000 ! ALLOW   ALI POL
                ! mp2/6-311g** mgua vib data, adm jr., 1/04
NC2  CT3  HA3    56.5000   107.5000 ! ALLOW   ALI POL
                ! mp2/6-311g** mgua vib data, adm jr., 1/04
NH1  C    CP1    80.000   116.5000 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH1  C    CT1    80.000   116.5000 ! ALLOW   ALI PEP POL ARO
                ! NMA Vib Modes (LK)
NH1  C    CT2    80.000   116.5000 ! ALLOW   ALI PEP POL ARO
                ! NMA Vib Modes (LK)
NH1  C    CT3    80.000   116.5000 ! ALLOW   ALI PEP POL ARO
                ! NMA Vib Modes (LK)
NH1  CT1  C      50.000   107.0000 ! ALLOW   PEP POL ARO ALI
                ! Alanine Dipeptide ab initio calc's (LK)
NH1  CT1  CC     50.000   107.0000 ! ALLOW   PEP POL ARO ALI
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
NH1  CT1  CD     50.000   107.0000 ! ALLOW   PEP POL ARO ALI
                ! adm jr. 5/02/91, acetic acid pure solvent
NH1  CT1  CT1    70.000   113.5000 ! ALLOW   ALI PEP POL ARO
                ! Alanine Dipeptide ab initio calc's (LK)
NH1  CT1  CT2    70.000   113.5000 ! ALLOW   ALI PEP POL ARO
                ! Alanine Dipeptide ab initio calc's (LK)
NH1  CT1  CT3    70.000   113.5000 ! ALLOW   ALI PEP POL ARO
                ! Alanine Dipeptide ab initio calc's (LK)
NH1  CT1  HB1    48.000   108.0000 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
NH1  CT2  C      50.000   107.0000 ! ALLOW   PEP POL ARO ALI
                ! Alanine Dipeptide ab initio calc's (LK)
NH1  CT2  CC     50.000   107.0000 ! ALLOW   PEP POL ARO ALI
                ! adm jr. 5/20/92, for asn,asp,gln,glu and cters
NH1  CT2  CD     50.000   107.0000 ! ALLOW   PEP POL ARO ALI
                ! adm jr. 5/02/91, acetic acid pure solvent
NH1  CT2  CT2    70.000   113.5000 ! ALLOW   ALI PEP POL ARO
                ! from NH1  CT1  CT2, for lactams, adm jr.
NH1  CT2  HA2    51.500   109.5000 ! ALLOW   ALI PEP POL ARO
                ! from NH1  CT3  HA, for lactams, adm jr.
NH1  CT2  HB2    48.000   108.0000 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
NH1  CT3  HA3    51.500   109.5000 ! ALLOW   ALI PEP POL ARO
                ! NMA crystal (JCS)
NH2  CC   CP1    80.000   112.5000 ! ALLOW ALI PEP POL ARO PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH2  CC   CT1   50.000    116.50   50.00   2.45000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 8/13/90 acetamide geometry and vibrations
NH2  CC   CT2   50.000    116.50   50.00   2.45000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 8/13/90 acetamide geometry and vibrations
NH2  CC   CT3   50.000    116.50   50.00   2.45000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 8/13/90 acetamide geometry and vibrations
NH2  CC   HA1   44.000    111.00   50.00   1.98000 ! ALLOW POL
                ! adm jr., 5/13/91, formamide geometry and vibrations
NH2  CT2  HB2   38.000    109.50   50.00   2.14000
                !from NH2  CT3  HA, neutral glycine, adm jr.
NH2  CT2  CD    52.000   108.0000
                !from CT2 CT2 CD, neutral glycine, adm jr.
NH2  CT2  CT2    67.700   110.0000 ! ALLOW   ALI POL
                !from NH3  CT2  CT2, neutral lysine
NH2  CT2  HA2   38.000    109.50   50.00   2.14000
                !from NH2  CT3  HA, neutral lysine
NH2  CT3  HA3   38.000    109.50   50.00   2.14000 ! ALLOW POL
                ! methylamine geom/freq, adm jr., 6/2/92
NH3  CT1  C      43.700   110.0000 ! ALLOW   PEP POL ARO ALI
                ! new aliphatics, adm jr., 2/3/92
NH3  CT1  CC     43.700   110.0000 ! ALLOW   PEP POL ARO ALI
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
NH3  CT1  CT1    67.700   110.0000 ! ALLOW   ALI POL
                ! new aliphatics, adm jr., 2/3/92
NH3  CT1  CT2    67.700   110.0000 ! ALLOW   ALI POL
                ! new aliphatics, adm jr., 2/3/92
NH3  CT1  CT3    67.700   110.0000 ! ALLOW   ALI POL
                ! new aliphatics, adm jr., 2/3/92
NH3  CT1  HB1    51.500   107.5000 ! ALLOW   ALI POL PEP
                ! new aliphatics, adm jr., 2/3/92
NH3  CT2  C      43.700   110.0000 ! ALLOW   PEP POL ARO ALI
                ! alanine (JCS)
NH3  CT2  CC     43.700   110.0000 ! ALLOW   PEP POL ARO ALI
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
NH3  CT2  CD     43.700   110.0000 ! ALLOW   PEP POL ARO ALI
                ! adm jr. 5/02/91, acetic acid pure solvent
NH3  CT2  CT2    67.700   110.0000 ! ALLOW   ALI POL
                ! methylammonium
NH3  CT2  CT3    67.700   110.0000 ! ALLOW   ALI POL
                ! ethylammonium
NH3  CT2  HA2   45.000    107.50   35.00   2.10100 ! ALLOW   ALI POL
                ! new stretch and bend; methylammonium (KK 03/10/92)
NH3  CT2  HB2    51.500   107.5000 ! ALLOW   ALI POL PEP
                ! for use on NTER -- from NH3 CT2HA (JCS) -- (LK)
NH3  CT3  HA3   45.000    107.50   35.00   2.10100 ! ALLOW   ALI POL
                ! new stretch and bend; methylammonium (KK 03/10/92)
NP   CP1  C      50.000   106.0000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NP   CP1  CC     50.000   106.0000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NP   CP1  CD     50.000   106.0000 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NP   CP1  CP2    70.000   108.5000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NP   CP1  HB1    51.500   107.5000 ! ALLOW ALI POL PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NP   CP3  CP2    70.000   108.5000 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NP   CP3  HA2    51.500   109.1500 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NR1  CPH1 CPH1  130.000   106.0000 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
NR1  CPH1 CT2    45.800   124.0000 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FC FROM CA CT2CT
NR1  CPH1 CT3    45.800   124.0000 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FC FROM CA CT2CT
NR1  CPH1 HR3   25.000    124.00   20.00   2.14000 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
NR1  CPH2 HR1   25.000    122.50   20.00   2.14000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
NR2  CPH1 CPH1  130.000   110.0000 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
NR2  CPH1 CT2    45.800   120.0000 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FC FROM CA CT2CT
NR2  CPH1 HR3   25.000    120.00   20.00   2.14000 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
NR2  CPH2 HR1   25.000    125.00   20.00   2.12000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
NR2  CPH2 NR1   130.000   112.5000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
NR3  CPH1 CPH1  145.000   108.0000 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
NR3  CPH1 CT2    45.800   122.0000 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FC FROM CA CT2CT
NR3  CPH1 HR1   22.000    122.00   15.00   2.18000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
NR3  CPH2 HR2   32.000    126.00   25.00   2.14000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
NR3  CPH2 NR3   145.000   108.0000 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
O    C    CP1    80.000   118.0000 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    CT1    80.000   121.0000 ! ALLOW   ALI PEP POL ARO
                ! Alanine Dipeptide ab initio calc's (LK)
O    C    CT2    80.000   121.0000 ! ALLOW   ALI PEP POL ARO
                ! Alanine Dipeptide ab initio calc's (LK)
O    C    CT3    80.000   121.0000 ! ALLOW   ALI PEP POL ARO
                ! Alanine Dipeptide ab initio calc's (LK)
O    C    H      50.000   121.7000 ! ALLOW   PEP POL ARO
                ! acetaldehyde (JCS)
O    C    N      80.000   122.5000 ! ALLOW PRO PEP POL ARO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    NH1    80.000   122.5000 ! ALLOW   PEP POL ARO
                ! NMA Vib Modes (LK)
O    CC   CP1    80.000   118.0000 ! ALLOW ALI PEP POL ARO PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    CC   CT1   15.000    121.00   50.00   2.44000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/10/91, acetamide update
O    CC   CT2   15.000    121.00   50.00   2.44000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/10/91, acetamide update
O    CC   CT3   15.000    121.00   50.00   2.44000 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 4/10/91, acetamide update
O    CC   HA1    44.000   122.0000 ! ALLOW POL
                ! adm jr., 5/13/91, formamide geometry and vibrations
O    CC   NH2   75.000    122.50   50.00   2.37000 ! ALLOW   POL PEP ARO
                ! adm jr. 4/10/91, acetamide update
OB   CD   CP1   70.000    125.00   20.00   2.44200 ! ALLOW ALI PEP POL ARO PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
OB   CD   CT1   70.000    125.00   20.00   2.44200 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 5/02/91, acetic acid pure solvent
OB   CD   CT2   70.000    125.00   20.00   2.44200 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 5/02/91, acetic acid pure solvent
OB   CD   CT3   70.000    125.00   20.00   2.44200 ! ALLOW   ALI PEP POL ARO
                ! adm jr. 5/02/91, acetic acid pure solvent
OC   CA   CA     40.000   120.0000 ! ALLOW  POL ARO
                ! adm jr. 8/27/91, phenoxide
OC   CC   CP1   40.000    118.00   50.00   2.38800 ! ALLOW ALI PEP POL ARO ION PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
OC   CC   CT1   40.000    118.00   50.00   2.38800 ! ALLOW   ALI PEP POL ARO ION
                ! adm jr. 7/23/91, correction, ACETATE (KK)
OC   CC   CT2   40.000    118.00   50.00   2.38800 ! ALLOW   ALI PEP POL ARO ION
                ! adm jr. 7/23/91, correction, ACETATE (KK)
OC   CC   CT3   40.000    118.00   50.00   2.38800 ! ALLOW   ALI PEP POL ARO ION
                ! adm jr. 7/23/91, correction, ACETATE (KK)
OC   CC   OC   100.000    124.00   70.00   2.22500 ! ALLOW   POL ION PEP ARO
                ! adm jr. 7/23/91, correction, ACETATE (KK)
OC   CT2  CT3    65.000   122.0000 ! ALLOW  ALC
                ! ethoxide 6-31+G* geom/freq, adm jr., 6/1/92
OC   CT2  HA2    65.000   118.3000 ! ALLOW  ALC
                ! ethoxide 6-31+G* geom/freq, adm jr., 6/1/92
OC   CT3  HA3    65.000   118.3000 ! ALLOW  ALC
                ! methoxide 6-31+G* geom/freq, adm jr., 6/1/92
OH1  CA   CA     45.200   120.0000 ! ALLOW   ARO ALC
                ! PARALLH19 WITH [122.3] (JES)
OH1  CD   CT2    55.000   110.5000 ! ALLOW   ALI PEP POL ARO ALC
                ! adm jr, 10/17/90, acetic acid vibrations
OH1  CD   CT3    55.000   110.5000 ! ALLOW   ALI PEP POL ARO ALC
                ! adm jr, 10/17/90, acetic acid vibrations
OH1  CD   OB    50.000    123.00  210.00   2.26200 ! ALLOW   PEP POL ARO ALC
                ! adm jr, 10/17/90, acetic acid vibrations
OH1  CT1  CT1    75.700   110.1000 ! ALLOW   ALI ALC ARO
                ! MeOH, EMB, 10/10/89
OH1  CT1  CT3    75.700   110.1000 ! ALLOW   ALI ALC ARO
                ! MeOH, EMB, 10/10/89
OH1  CT1  HA1    45.900   108.8900 ! ALLOW   ALI ALC ARO
                ! MeOH, EMB, 10/10/89
OH1  CT2  CT1    75.700   110.1000 ! ALLOW   ALI ALC ARO
                ! MeOH, EMB, 10/10/89
OH1  CT2  CT2    75.700   110.1000 ! ALLOW   ALI ALC ARO
                ! MeOH, EMB, 10/10/89
OH1  CT2  CT3    75.700   110.1000 ! ALLOW   ALI ALC ARO
                ! MeOH, EMB, 10/10/89
OH1  CT2  HA2    45.900   108.8900 ! ALLOW   ALI ALC ARO
                ! MeOH, EMB, 10/10/89
OH1  CT3  HA3    45.900   108.8900 ! ALLOW   ALI ALC ARO
                ! MeOH, EMB, 10/10/89
OS   CD   CP1   55.000    109.00   20.00   2.32600 ! ALLOW POL PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
OS   CD   CT1   55.000    109.00   20.00   2.32600 ! ALLOW POL PEP
                ! adm jr., 4/05/91, for PRES CT1 from methylacetate
OS   CD   CT2   55.000    109.00   20.00   2.32600 ! ALLOW POL PEP
                ! adm jr., 4/05/91, for PRES CT1 from methylacetate
OS   CD   CT3   55.000    109.00   20.00   2.32600 ! ALLOW POL PEP
                ! adm jr., 4/05/91, for PRES CT1 from methylacetate
OS   CD   OB    90.000    125.90  160.00   2.25760 ! ALLOW  PEP POL
                ! adm jr. 3/19/92, from lipid methyl acetate
OS   CT2  HA2    60.000   109.5000 ! ALLOW PEP POL
                ! adm jr. 4/05/91, for PRES CT1 from methyl acetate
OS   CT3  HA3    60.000   109.5000 ! ALLOW PEP POL
                ! adm jr. 4/05/91, for PRES CT1 from methyl acetate
S    CT2  CT1    58.000   112.5000 ! ALLOW   ALI SUL ION
                ! as in expt.MeEtS & DALC crystal,  5/15/92
S    CT2  CT2    58.000   114.5000 ! ALLOW   ALI SUL ION
                ! expt. MeEtS,     3/26/92 (FL)
S    CT2  CT3    58.000   114.5000 ! ALLOW   ALI SUL ION
                ! expt. MeEtS,     3/26/92 (FL)
S    CT2  HA2    46.100   111.3000 ! ALLOW   ALI SUL ION
                ! vib. freq. and HF/6-31G* geo. (DTN) 8/24/90
S    CT3  HA3    46.100   111.3000 ! ALLOW   ALI SUL ION
                ! vib. freq. and HF/6-31G* geo. (DTN) 8/24/90
SM   CT2  CT1    58.000   112.5000 ! ALLOW   ALI SUL ION
                ! as in expt.MeEtS & DALC crystal,  5/15/92
SM   CT2  CT3    58.000   112.5000 ! ALLOW   ALI SUL ION
                ! diethyldisulfide, as in expt.MeEtS & DALC crystal,  5/15/92
SM   CT2  HA2    38.000   111.0000 ! ALLOW   ALI SUL ION
                ! new S-S atom type 8/24/90
SM   CT3  HA3    38.000   111.0000 ! ALLOW   ALI SUL ION
                ! new S-S atom type 8/24/90
SM   SM   CT2    72.500   103.3000 ! ALLOW   ALI SUL ION
                ! expt. dimethyldisulfide,    3/26/92 (FL)
SM   SM   CT3    72.500   103.3000 ! ALLOW   ALI SUL ION
                ! expt. dimethyldisulfide,    3/26/92 (FL)
SS   CS   CT3    55.000   118.0000 ! ALLOW SUL
                ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
SS   CS   HA2    40.000   112.3000 ! ALLOW SUL
                ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
SS   CS   HA3    40.000   112.3000 ! ALLOW SUL
                ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
O    CD   HR1    75.000   121.0000 ! acetaldehyde, benzaldehyde, 3ALP, retinal
!For GLU/HSP, Zhu
NH1  CT1  CT2A   70.000   113.5000 ! from NH1  CT1  CT2
HB1  CT1  CT2A   35.000   111.0000 ! from HB1  CT1  CT2
CT2A CT1  C      52.000   108.0000 ! from CT2  CT1  C
CT1  CT2A HA2    26.500   110.1000  22.53   2.17900 ! from HA2  CT2  CT1
CT1  CT2A CT2    58.350   113.5000  11.16   2.56100 ! from CT2  CT2  CT1
HA2  CT2A HA2    35.500   109.0000   5.40   1.80200 ! from HA2  CT2  HA2
HA2  CT2A CT2    26.500   110.1000  22.53   2.17900 ! from HA2  CT2  CT2
CT2A CT2  HA2    26.500   110.1000  22.53   2.17900 ! from HA2  CT2  CT2
CT2A CT2  CC     52.000   108.0000 ! from CT2  CT2  CC
CT1  CT2A CPH1   58.350   113.0000 ! from CT1  CT2  CPH1 
HA2  CT2A CPH1   33.430   109.5000 ! from HA2  CT2  CPH1
CT2A CPH1 CPH1   45.800   130.0000 ! from CT2  CPH1 CPH1
CT2A CPH1 NR3    45.800   122.0000 ! from NR3  CPH1 CT2
!ASP, CT2->CT2A, jshim
CT1  CT2A CC     52.000   108.0000 ! from CT1  CT2  CC
HA2  CT2A CC     33.000   109.5000  30.00   2.16300 ! from HA2  CT2  CC
OC   CC   CT2A   40.000   118.0000  50.00   2.38800 ! from OC   CC   CT2
NH3  CT1  CT2A   67.700   110.0000 ! from NH3  CT1  CT2
CT2A CT1  CD     52.000   108.0000 ! from CT2  CT1  CD
! RESI CYSM and PRES CYSD
NH2  CT1  CS     67.700   110.0000 ! from NH2  CT1  CT2 , kevo
CS   CT1  C      52.000   108.0000 ! from CT2  CT1  C   , kevo
CS   CT1  CC     52.000   108.0000 ! from CT2  CT1  CC  , kevo
CS   CT1  CD     52.000   108.0000 ! from CT2  CT1  CD  , kevo
HB1  CT1  CS     35.000   111.0000 ! from HB1  CT1  CT2 , kevo
NH1  CT1  CS     70.000   113.5000 ! from NH1  CT1  CT2 , kevo
NH3  CT1  CS     67.700   110.0000 ! from NH3  CT1  CT2 , kevo
SS   CS   CT1    55.000   118.0000 ! from SS   CS   CT3 , kevo
HA2  CS   CT1    34.600   110.10    22.53   2.17900 ! from HA2 CS CT3 to be consistent with SS CS CT1, kevo
! PRES SERD
OC   CT2  CT1    65.000   122.0000 ! from OC   CT2  CT3 , kevo
'''

def temp1():
    
    # Split the string into lines
    lines = data.split('\n')

    # Process each line to remove comments (everything after the '!') and strip whitespace
    processed_lines = [line.split('!')[0].strip() for line in lines]

    # Remove empty lines
    filtered_lines = [line for line in processed_lines if line]

    # Convert the list back into a string
    result = '\n'.join(filtered_lines)

    # Print the result
    print(result)

def temp2():
    # Split the string into lines
    lines = data.strip().split('\n')
    
    # Initialize an empty list to store the transformed lines
    transformed_lines = []
    
    # Iterate through each line
    for line in lines:
        # Extract parts of the string based on fixed positions
        # Assuming the format and positions are consistent as per the example
        part1 = line[:7]  # Label
        part2 = line[17:31]  # First number to keep
        part3 = line[31:40]  # Second number to keep
        
        # Combine the parts with the correct spacing to maintain alignment
        new_line = part1 + part2 + part3
        transformed_lines.append(new_line)
    
    # Join the transformed lines back into a single string with newline characters
    transformed_string = '\n'.join(transformed_lines)
    
    print(transformed_string)


def temp3():
    lines = data.strip().split('\n')
    for line in lines:
        parts = line.split()
        if len(parts) > 6:  # Ensure there are enough parts to include Kub and S0
            try:
                Kub = float(parts[5])
                S0 = float(parts[6])
                # atom_types = ' '.join(parts[:3])
                print(f"{parts[0]:5}{parts[1]:5}{parts[2]:5}  {Kub:.2f}    {S0:.5f}")
            except ValueError:
                continue

temp3()