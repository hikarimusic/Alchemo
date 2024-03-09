data = '''
!V(dihedral) = Kchi(1 + cos(n(chi) - delta))
!
!Kchi: kcal/mole
!n: multiplicity
!delta: degrees
!
!atom types             Kchi    n   delta
!
!Neutral N and C termini
NH2  CT1  C    O        0.0000  1     0.00
NH2  CT2  C    O        0.0000  1     0.00   ! Neutral Gly Nterminus
NH2  CT1  C    NH1      0.0000  1     0.00
NH2  CT2  C    NH1      0.0000  1     0.00   ! Neutral Gly Nterminus
NH2  CT2  C    N        0.0000  1     0.00   ! NNEU
NH2  CT1  CT2A HA2      0.2000  3     0.00   ! From X    CT1  CT2  X
H    NH2  CT1  CT1      0.0000  1     0.00
H    NH2  CT1  C        0.0000  1     0.00
H    NH2  CT1  CC       0.0000  1     0.00   ! From neutral N-termini
H    NH2  CT1  CD       0.0000  1     0.00
H    NH2  CT2  C        0.0000  1     0.00   ! Neutral Gly Nterminus
H    NH2  CT2  CD       0.0000  1     0.00   ! Neutral Gly Nterminus
H    NH2  CT1  HB1      0.1100  3     0.00   ! From LSN HC-NH2-CT2-HA
H    NH2  CT2  HB2      0.1100  3     0.00   ! From LSN HC-NH2-CT2-HA, Neutral Gly Nterminus
H    NH2  CT1  CT2      0.1100  3     0.00   ! From LSN HC-NH2-CT2-CT2
H    NH2  CT1  CT2A     0.1100  3     0.00   ! From LSN HC-NH2-CT2-CT2
H    NH2  CT1  CT3      0.1100  3     0.00   ! From LSN HC-NH2-CT2-CT2
CC   CT2A CT1  NH2      0.6800  1   180.00   ! From CC   CT2A CT1  NH1
CC   CT2A CT1  NH2      0.1000  2   180.00
CC   CT2A CT1  NH2      0.3800  3     0.00
CD   CT1  CT2A CC       1.6100  1   180.00   ! From C    CT1  CT2A CC
CD   CT1  CT2A CC       1.2900  2   180.00
CD   CT1  CT2A CC       0.5900  3   180.00

!Indole/Tryptophan
CAI  CA   CA   CAI      3.1000  2   180.00 ! from CA CA CA CA
CA   CPT  CPT  CA       3.0000  2   180.00 ! atm, methylindole, 1/17/04	
CAI  CPT  CPT  CAI      3.0000  2   180.00 ! atm, methylindole, 1/17/04	
CA   CY   CPT  CA       3.0000  2   180.00 ! atm, methylindole, 1/17/04
CA   CY   CPT  CAI      3.0000  2   180.00 ! atm, methylindole, 1/17/04
CA   NY   CPT  CA       3.0000  2   180.00 ! atm, methylindole, 1/17/04
CPT  CA   CA   CA       3.0000  2   180.00 ! JWK 05/14/91 fit to indole
CPT  CPT  CA   CA       3.0000  2   180.00 ! JWK 05/14/91 fit to indole
CA   NY   CPT  CAI      3.0000  2   180.00 ! atm, methylindole, 1/17/04
CPT  CAI  CA   CA       3.0000  2   180.00 ! JWK 05/14/91 fit to indole
CPT  CPT  CAI  CA       3.0000  2   180.00 ! JWK 05/14/91 fit to indole
CPT  CPT  CY   CA       5.0000  2   180.00 ! atm, methylindole, 1/17/04
CPT  CPT  NY   CA       6.5000  2   180.00 ! atm, methylindole, 1/17/04
CT3  CY   CPT  CA       2.5000  2   180.00 ! atm, methylindole, r6r5
CT3  CY   CPT  CAI      2.5000  2   180.00 ! atm, methylindole, r6r5
CT3  CY   CPT  CPT      3.0000  2   180.00 ! atm, methylindole, meth
CT2  CY   CPT  CA       2.5000  2   180.00 ! atm, methylindole, r6r5
CT2  CY   CPT  CAI      2.5000  2   180.00 ! atm, methylindole, r6r5
CT2  CY   CPT  CPT      3.0000  2   180.00 ! atm, methylindole, meth
CY   CA   NY   CPT      6.0000  2   180.00 ! atm, methylindole, 1/17/04
CY   CPT  CA   CA       4.0000  2   180.00 ! atm, methylindole, 1/17/04
CY   CPT  CPT  CA       4.0000  2   180.00 ! atm, methylindole, 1/17/04
CY   CPT  CAI  CA       4.0000  2   180.00 ! atm, methylindole, 1/17/04
CY   CPT  CPT  CAI      4.0000  2   180.00 ! atm, methylindole, 1/17/04
H    NY   CA   CY       0.0500  2   180.00 ! atm, methylindole, 1/17/04
H    NY   CPT  CA       0.2000  2   180.00 ! atm, methylindole, 1/17/04
H    NY   CPT  CAI      0.2000  2   180.00 ! atm, methylindole, 1/17/04
H    NY   CPT  CPT      0.8500  2   180.00 ! atm, methylindole, 1/17/04
HP   CAI  CA   CA       4.2000  2   180.00 ! from HP CA CA CA
HP   CA   CA   CPT      3.0000  2   180.00 ! JWK 05/14/91 fit to indole
HP   CA   CPT  CPT      3.0000  2   180.00 ! JWK indole 05/14/91
HP   CA   CPT  CY       4.0000  2   180.00 ! atm, methylindole, 1/17/04
HP   CA   CA   CAI      4.2000  2   180.00 ! from HP CA CA CA
HP   CA   CAI  CPT      3.0000  2   180.00 ! from HP CA CA CPT
HP   CAI  CA   HP       2.4000  2   180.00 ! from HP CA CA HP
HP   CAI  CPT  CPT      3.0000  2   180.00 ! from HP CA CPT CPT
HP   CAI  CPT  CY       4.0000  2   180.00 ! from HP CA CPT CY, r6r5
HP   CA   CY   CPT      2.8000  2   180.00 ! adm jr., 12/30/91, for jwk
HP   CA   CY   CT3      1.2000  2   180.00 ! atm, methylindole
HP   CA   CY   CT2      1.2000  2   180.00 ! atm, methylindole
HP   CA   NY   CPT      2.6000  2   180.00 ! adm jr., 12/30/91, for jwk
HP   CA   NY   H        0.4000  2   180.00 ! JWK 05/14/91 fit to indole
HP   CY   CA   HP       1.0000  2   180.00 ! JWK 05/14/91 fit to indole
HP   CY   CPT  CA       2.8000  2   180.00 ! JWK 05/14/91 fit to indole
HP   CY   CPT  CAI      2.8000  2   180.00 ! JWK 05/14/91 fit to indole
HP   CY   CPT  CPT      2.6000  2   180.00 ! JWK 05/14/91 fit to indole
NY   CA   CY   CPT      5.0000  2   180.00 ! atm, methylindole, 1/17/04
NY   CA   CY   CT3      2.5000  2   180.00 ! atm, methylindole, from NY   CA   CY   CT3
NY   CA   CY   CT2      2.5000  2   180.00 ! atm, methylindole, from NY   CA   CY   CT3
NY   CA   CY   HP       3.5000  2   180.00 ! JWK indole 05/14/91
NY   CPT  CA   CA       3.0000  2   180.00 ! atm, methylindole, 1/17/04, r6r5 
NY   CPT  CA   HP       3.0000  2   180.00 ! JWK 05/14/91 fit to indole, r6r5
NY   CPT  CPT  CA       4.0000  2   180.00 ! atm, methylindole, 1/17/04, bfly
NY   CPT  CAI  CA       3.0000  2   180.00 ! atm, methylindole, 1/17/04
NY   CPT  CAI  HP       3.0000  2   180.00 ! JWK 05/14/91 fit to indole, r6r5
NY   CPT  CPT  CAI      4.0000  2   180.00 ! atm, methylindole, 1/17/04, bfly
NY   CPT  CPT  CY       6.5000  2   180.00 ! JWK 05/14/91 fit to indole,  r5 t1
CT3  CT2  CY   CA       0.3800  2     0.00 ! trp, from ethylbenzene, adm jr., 3/7/92
CT3  CT2  CY   CPT      0.2500  2   180.00 ! atm 1/14/04 3-ethylindole
CT3  CT2  CY   CPT      0.3000  3     0.00 ! atm 1/14/04 3-ethylindole
HA3  CT3  CY   CA       0.0100  3     0.00 ! atm, methylindole, 1/17/04
HA3  CT3  CY   CPT      0.2000  3     0.00 ! atm, methylindole, 1/17/04
HA2  CT2  CY   CA       0.0100  3     0.00 ! atm, methylindole, 1/17/04
HA2  CT2  CY   CPT      0.2000  3     0.00 ! atm, methylindole, 1/17/04
X    CS   SS   X        0.0000  3     0.20 ! guess
                !from methanethiol, HS S CT3 HA
                !adm jr., 7/01
C    CT1  NH1  C        0.2000  1   180.00 ! ALLOW PEP
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c
C    CT2  NH1  C        0.2000  1   180.00 ! ALLOW PEP
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c
C    N    CP1  C        0.8000  3     0.00 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CA   CA   CA   CA       3.1000  2   180.00 ! ALLOW   ARO
                ! JES 8/25/89
!CA   CT2  CT1  C        0.0400  3     0.00 ! ALLOW   ARO
                ! 2.7 kcal/mole CH3 rot in ethylbenzene, adm jr, 3/7/92
CC   CP1  N    C        0.8000  3     0.00 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CC   CT1  CT2  CA       0.0400  3     0.00 ! ALLOW   ARO
                ! 2.7 kcal/mole CH3 rot in ethylbenzene, adm jr, 3/7/92
CC   CT1  NH1  C        0.2000  1   180.00 ! ALLOW PEP POL
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c
!CC   CT2  NH1  C        0.2000  1   180.00 ! ALLOW PEP POL
!                ! Alanine dipeptide; NMA; acetate; etc. adm jr., 3/3/93c
CC   CT2  NH1  C        2.0000  1   180.00 ! ALLOW PEP POL
                ! Based on Gly3 data from graf et al, RB 7/1/11
CD   CP1  N    C        0.0000  1   180.00 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CD   CT1  NH1  C        0.2000  1   180.00 ! ALLOW PEP POL
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c
!CD   CT2  NH1  C        0.2000  1   180.00 ! ALLOW PEP POL
!               ! Alanine dipeptide; NMA; acetate; etc. backbon adm jr., 3/3/93c
CD   CT2  NH1  C        2.0000  1   180.00 ! ALLOW PEP POL
                ! Based on Gly3 data from graf et al, RB 7/1/11
CE1  CE1  CT3  HA3      0.0300  3     0.00 ! 
		! for butene, yin/adm jr., 12/95
CE2  CE1  CT2  CT3      0.5000  1   180.00 !
                ! 1-butene, adm jr., 2/00 update
CE2  CE1  CT2  CT3      1.3000  3   180.00 !
		! 1-butene, adm jr., 2/00 update
CE2  CE1  CT2  HA2      0.1200  3     0.00 ! 
		! for butene, yin/adm jr., 12/95
CE2  CE1  CT3  HA3      0.0500  3   180.00 ! 
		! for propene, yin/adm jr., 12/95
CP1  C    N    CP1      2.7500  2   180.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP1  C    N    CP1      0.3000  4     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP2  CP1  N    C        0.8000  3     0.00 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP2  CP3  N    C        0.0000  3   180.00 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP2  CP3  N    CP1      0.1000  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP2  CP3  NP   CP1      0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  N    C    CP1      2.7500  2   180.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  N    C    CP1      0.3000  4     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  N    CP1  C        0.1000  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  N    CP1  CC       0.1000  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  N    CP1  CP2      0.1000  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  NP   CP1  C        0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  NP   CP1  CC       0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  NP   CP1  CD       0.0800  3     0.00 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CP3  NP   CP1  CP2      0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CPH2 NR1  CPH1 CPH1    14.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
CPH2 NR2  CPH1 CPH1    14.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
CPH2 NR3  CPH1 CPH1    12.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
CT1  C    N    CP1      2.7500  2   180.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT1  C    N    CP1      0.3000  4     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT1  C    N    CP3      2.7500  2   180.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT1  C    N    CP3      0.3000  4     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT1  C    NH1  CT1      1.6000  1     0.00 !  ALLOW PEP
                ! Revised to adjust NMA cis/trans energy difference. (LK)
CT1  C    NH1  CT1      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
CT1  CT1  NH1  C        1.8000  1     0.00 ! ALLOW PEP
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c
CT1  NH1  C    CP1      1.6000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT1  NH1  C    CP1      2.5000  2   180.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT2  C    N    CP1      2.7500  2   180.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT2  C    N    CP1      0.3000  4     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT2  C    N    CP3      2.7500  2   180.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT2  C    N    CP3      0.3000  4     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT2  C    NH1  CT1      1.6000  1     0.00 !  ALLOW PEP
                ! Revised to adjust NMA cis/trans energy difference. (LK)
CT2  C    NH1  CT1      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
CT2  C    NH1  CT2      1.6000  1     0.00 !  ALLOW PEP
                ! Revised to adjust NMA cis/trans energy difference. (LK)
CT2  C    NH1  CT2      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
CT2  C    NH1  CT3      1.6000  1     0.00 !  ALLOW PEP
                ! from CT2  C    NH1  CT2, adm jr. 10/21/96
CT2  C    NH1  CT3      2.5000  2   180.00 !  ALLOW PEP
                ! from CT2  C    NH1  CT2, adm jr. 10/21/96
CT2  CA   CA   CA       3.1000  2   180.00 ! ALLOW   ARO
                ! JES 8/25/89 toluene and ethylbenzene
CT2  CPH1 NR1  CPH2     3.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FROM HA CPH1 NR1 CPH2
CT2  CPH1 NR2  CPH2     3.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FROM HA CPH1 NR2 CPH2
CT2  CPH1 NR3  CPH2     2.5000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
CT2  CT1  NH1  C        1.8000  1     0.00 ! ALLOW PEP
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c
CT2  CT2  CPH1 CPH1     0.4000  1     0.00 ! ALLOW ARO
                ! 4-methylimidazole 4-21G//6-31G* rot bar. ADM JR., 9/4/89
!aliphatic chain parameters compatible with the revised side-chain parameters, from all22_carb>>all27_lip>>all31
                ! lower butane gauche conformer
CT2  CT2  CT2  CT2      0.10    2   180.00 ! alkane, 4/98, adm jr.
CT2  CT2  CT2  CT2      0.15    4     0.00 ! alkane, 4/98, adm jr.
CT2  CT2  CT2  CT2      0.10    6   180.00 ! alkane, 4/98, adm jr.
CT2  CT2  CT2  CT3      0.10    2   180.00 ! alkane, 4/98, adm jr.
CT2  CT2  CT2  CT3      0.15    4     0.00 ! alkane, 4/98, adm jr.
CT2  CT2  CT2  CT3      0.10    6   180.00 ! alkane, 4/98, adm jr.
!
CT2  CT2  NH1  C        1.8000  1     0.00 ! ALLOW PEP
                ! from CT2  CT1  NH1  C, for lactams, adm jr.
CT2  NH1  C    CP1      1.6000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT2  NH1  C    CP1      2.5000  2   180.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT2  NH1  C    CT1      1.6000  1     0.00 !  ALLOW PEP
                ! Revised to adjust NMA cis/trans energy difference. (LK)
CT2  NH1  C    CT1      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
CT2  SM   SM   CT2      1.0000  1     0.00 ! ALLOW   ALI SUL ION
                ! improved CSSC dihedral in DMDS  5/15/92 (FL)
CT2  SM   SM   CT2      4.1000  2     0.00 ! ALLOW   ALI SUL ION
                ! mp 6-311G** dimethyldisulfide,  3/26/92 (FL)
CT2  SM   SM   CT2      0.9000  3     0.00 ! ALLOW   ALI SUL ION
                ! improved CSSC dihedral in DMDS  5/15/92 (FL)
CT3  C    N    CP1      2.7500  2   180.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT3  C    N    CP1      0.3000  4     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT3  C    N    CP3      2.7500  2   180.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT3  C    N    CP3      0.3000  4     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT3  C    NH1  CT1      1.6000  1     0.00 !  ALLOW PEP
                ! Revised to adjust NMA cis/trans energy difference. (LK)
CT3  C    NH1  CT1      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
CT3  C    NH1  CT2      1.6000  1     0.00 !  ALLOW PEP
                ! for acetylated GLY N-terminus, adm jr.
CT3  C    NH1  CT2      2.5000  2   180.00 !  ALLOW PEP
                ! for acetylated GLY N-terminus, adm jr.
CT3  C    NH1  CT3      1.6000  1     0.00 !  ALLOW PEP
                ! Revised to adjust NMA cis/trans energy difference. (LK)
CT3  C    NH1  CT3      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
CT3  CA   CA   CA       3.1000  2   180.00 ! ALLOW   ARO
                ! toluene, adm jr., 3/7/92
CT3  CE1  CE2  HE2      5.2000  2   180.00 ! 
		! for propene, yin/adm jr., 12/95
CT3  CPH1 NR1  CPH2     3.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FROM HA CPH1 NR1 CPH2
CT3  CT1  NH1  C        1.8000  1     0.00 ! ALLOW PEP
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c
CT3  CT2  CA   CA       0.2300  2   180.00 ! ALLOW   ARO ALI
                ! ethylbenzene ethyl rotation, adm jr. 3/7/92
CT3  CT2  CPH1 CPH1     0.2000  1     0.00 ! ALLOW ARO
                ! 4-ethylimidazole 4-21G rot bar, adm jr. 3/4/92
CT3  CT2  CPH1 CPH1     0.2700  2     0.00 ! ALLOW ARO
                ! 4-ethylimidazole 4-21G rot bar, adm jr. 3/4/92
CT3  CT2  CPH1 CPH1     0.0000  3     0.00 ! ALLOW ARO
                ! 4-ethylimidazole 4-21G rot bar, adm jr. 3/4/92
CT3  CT2  S    CT3      0.2400  1   180.00 ! ALOW    ALI SUL ION
                ! expt. MeEtS,      3/26/92 (FL)
CT3  CT2  S    CT3      0.3700  3     0.00 ! ALOW    ALI SUL ION
                ! DTN 8/24/90
CT3  NH1  C    CP1      1.6000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT3  NH1  C    CP1      2.5000  2   180.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
CT3  NH1  C    CT1      1.6000  1     0.00 !  ALLOW PEP
                ! Revised to adjust NMA cis/trans energy difference. (LK)
CT3  NH1  C    CT1      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
CT3  S    CT2  CT2      0.2400  1   180.00 ! ALOW    ALI SUL ION
                ! expt. MeEtS,      3/26/92 (FL)
CT3  S    CT2  CT2      0.3700  3     0.00 ! ALOW    ALI SUL ION
                ! expt. MeEtS,      3/26/92 (FL)
CT3  SM   SM   CT3      1.0000  1     0.00 ! ALLOW   ALI SUL ION
                ! improved CSSC dihedral in DMDS  5/15/92 (FL)
CT3  SM   SM   CT3      4.1000  2     0.00 ! ALLOW   ALI SUL ION
                ! mp 6-311G** dimethyldisulfide,   3/26/92 (FL)
CT3  SM   SM   CT3      0.9000  3     0.00 ! ALLOW   ALI SUL ION
                ! improved CSSC dihedral in DMDS  5/15/92 (FL)
H    NH1  C    CP1      2.5000  2   180.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
H    NH1  C    CT1      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
H    NH1  C    CT2      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
H    NH1  C    CT3      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
H    NH1  CT1  C        0.0000  1     0.00 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
H    NH1  CT1  CC       0.0000  1     0.00 ! ALLOW PEP POL
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
H    NH1  CT1  CD       0.0000  1     0.00 ! ALLOW PEP POL
                ! adm jr. 5/02/91, acetic acid pure solvent
H    NH1  CT1  CT1      0.0000  1     0.00 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
H    NH1  CT1  CT2      0.0000  1     0.00 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
H    NH1  CT1  CT3      0.0000  1     0.00 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
H    NH1  CT2  C        0.0000  1     0.00 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
H    NH1  CT2  CC       0.0000  1     0.00 ! ALLOW PEP POL
                ! Alanine dipeptide; NMA; acetate; etc. backbone param. RLD 3/22/92
H    NH1  CT2  CD       0.0000  1     0.00 ! ALLOW PEP POL
                ! adm jr. 5/02/91, acetic acid pure solvent
H    NH1  CT2  CT2      0.0000  1     0.00 ! ALLOW PEP
                ! from H    NH1  CT2  CT3, for lactams, adm jr.
H    NH1  CT2  CT3      0.0000  1     0.00 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
H    NH2  CC   CT1      1.4000  2   180.00 !  ALLOW   PEP POL ARO PRO
                ! adm jr. 4/10/91, acetamide update
H    NH2  CC   CT2      1.4000  2   180.00 !  ALLOW   PEP POL ARO PRO
                ! adm jr. 4/10/91, acetamide update
H    NH2  CC   CT3      1.4000  2   180.00 !  ALLOW   PEP POL ARO PRO
                ! adm jr. 4/10/91, acetamide update
H    NH2  CC   CP1      2.5000  2   180.00 ! ALLOW PEP POL ARO PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
H    NR1  CPH1 CPH1     1.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 7/20/89
H    NR1  CPH1 CT2      1.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 7/22/89, FROM HA CPH1 NR1 H
H    NR1  CPH1 CT3      1.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 7/22/89, FROM HA CPH1 NR1 H
H    NR3  CPH1 CPH1     1.4000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
H    NR3  CPH1 CT2      3.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 7/22/89, FROM HC NR3 CPH1 HA
H    NR3  CPH1 CT3      3.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 7/22/89, FROM HC NR3 CPH1 HA
H    OH1  CA   CA       0.9900  2   180.00 ! ALLOW   ARO ALC
                ! phenol OH rot bar, 3.37 kcal/mole, adm jr. 3/7/92
H    OH1  CT1  CT3      1.3300  1     0.00 ! ALLOW ALC
                ! 2-propanol OH hf/6-31g* torsional surface, adm jr., 3/2/93
H    OH1  CT1  CT3      0.1800  2     0.00 ! ALLOW ALC
                ! 2-propanol OH hf/6-31g* torsional surface, adm jr., 3/2/93
H    OH1  CT1  CT3      0.3200  3     0.00 ! ALLOW ALC
                ! 2-propanol OH hf/6-31g* torsional surface, adm jr., 3/2/93
H    OH1  CT2  CT2      1.3000  1     0.00 ! ALLOW ALC
                ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93
H    OH1  CT2  CT2      0.3000  2     0.00 ! ALLOW ALC
                ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93
H    OH1  CT2  CT2      0.4200  3     0.00 ! ALLOW ALC
                ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93
H    OH1  CT2  CT3      1.3000  1     0.00 ! ALLOW ALC
                ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93
H    OH1  CT2  CT3      0.3000  2     0.00 ! ALLOW ALC
                ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93
H    OH1  CT2  CT3      0.4200  3     0.00 ! ALLOW ALC
                ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93
HA1  CC   NH2  H        1.4000  2   180.00 !  ALLOW PEP POL
                ! adm jr. 4/10/91, acetamide update
HA2  CP3  N    C        0.0000  3   180.00 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HA2  CP3  N    CP1      0.1000  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HA2  CP3  NP   CP1      0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HA1  CT1  CT2  CA       0.0400  3     0.00 ! ALLOW   ARO
                ! 2.7 kcal/mole CH3 rot in ethylbenzene, adm jr, 3/7/92
HA2  CT2  CPH1 CPH1     0.0000  3     0.00 ! ALLOW ARO
                ! 4-methylimidazole 4-21G//6-31G* rot bar. adm jr., 9/4/89
HA2  CT2  NH1  C        0.0000  3     0.00 ! ALLOW PEP
                ! LK for autogenerate dihe, sp2-methyl, no dihedral potential
HA2  CT2  NH1  H        0.0000  3     0.00 ! ALLOW PEP
                ! LK for autogenerate dihe, sp2-methyl, no dihedral potential
HA2  CT2  S    CT3      0.2800  3     0.00 ! ALLOW   ALI SUL ION
                ! DTN 8/24/90
HA3  CT3  CPH1 CPH1     0.0000  3     0.00 ! ALLOW ARO
                ! 4-methylimidazole 4-21G//6-31G* rot bar. adm jr., 9/4/89
HA3  CT3  CS   HA2      0.1600  3     0.00 ! ALLOW SUL
                ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
HA3  CT3  CS   HA3      0.1600  3     0.00 ! ALLOW SUL
                ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
HA3  CT3  CT2  CA       0.0400  3     0.00 ! ALLOW   ARO
                ! 2.7 kcal/mole CH3 rot in ethylbenzene, adm jr, 3/7/92
HA3  CT3  NH1  C        0.0000  3     0.00 ! ALLOW PEP
                ! LK for autogenerate dihe, sp2-methyl, no dihedral potential
HA3  CT3  NH1  H        0.0000  3     0.00 ! ALLOW PEP
                ! LK for autogenerate dihe, sp2-methyl, no dihedral potential
HA3  CT3  S    CT2      0.2800  3     0.00 ! ALLOW   ALI SUL ION
                ! DTN 8/24/90
HE1  CE1  CE1  HE1      1.0000  2   180.00 ! 
                ! 2-butene, adm jr., 8/98 update
CT3  CE1  CE1  HE1      1.0000  2   180.00 !
                ! 2-butene, adm jr., 8/98 update
HE1  CE1  CE2  HE2      5.2000  2   180.00 ! 
		! for propene, yin/adm jr., 12/95
HE1  CE1  CT2  HA2      0.0000  3     0.00
		! butene, adm jr., 2/00 update
HE1  CE1  CT2  CT3      0.1200  3     0.00 ! 
		! for butene, yin/adm jr., 12/95
HE1  CE1  CT3  HA3      0.0000  3     0.00
		! butene, adm jr., 2/00 update
HE2  CE2  CE1  CT2      5.2000  2   180.00 ! 
		! for butene, yin/adm jr., 12/95
HB1  CP1  N    C        0.8000  3     0.00 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HB1  CP1  N    CP3      0.1000  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HB1  CP1  NP   CP3      0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HB1  CT1  NH1  C        0.0000  1     0.00 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
HB1  CT1  NH1  H        0.0000  1     0.00 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
HB2  CT2  NH1  C        0.0000  1     0.00 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
HB2  CT2  NH1  H        0.0000  1     0.00 ! ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
HC   NH2  CT2  HB2      0.1100  3     0.00
                !from X CT3 NH2 X, neutral glycine, adm jr.
HC   NH2  CT2  CD       0.1100  3     0.00
                !from X CT3 NH2 X, neutral glycine, adm jr.
HC   NH2  CT2  CT2       0.1100  3     0.00
                !from X CT3 NH2 X, neutral lysine
HC   NH2  CT2  HA2      0.1100  3     0.00
                !from X CT3 NH2 X, neutral lysine
HC   NP   CP1  C        0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HC   NP   CP1  CC       0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HC   NP   CP1  CD       0.0800  3     0.00 ! ALLOW PRO PEP
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HC   NP   CP1  CP2      0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HC   NP   CP1  HB1      0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HC   NP   CP3  CP2      0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HC   NP   CP3  HA2      0.0800  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
HP   CA   CA   CA       4.2000  2   180.00 ! ALLOW   ARO
                ! JES 8/25/89 benzene
HP   CA   CA   CT2      4.2000  2   180.00 ! ALLOW   ARO
                ! JES 8/25/89 toluene and ethylbenzene
HP   CA   CA   CT3      4.2000  2   180.00 ! ALLOW   ARO
                ! toluene, adm jr., 3/7/92
HP   CA   CA   HP       2.4000  2   180.00 ! ALLOW   ARO
                ! JES 8/25/89 benzene
HR1  CPH1 CPH1 CT2      1.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
HR1  CPH1 CPH1 CT3      1.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
HR1  CPH1 CPH1 HR1      1.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90, his
HR1  CPH1 NR3  CPH2     2.5000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
HR1  CPH1 NR3  H        3.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
HR1  CPH2 NR1  CPH1     3.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
HR1  CPH2 NR1  H        1.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
HR1  CPH2 NR2  CPH1     3.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
HR2  CPH2 NR3  CPH1     3.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
HR2  CPH2 NR3  H        0.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90, YES, 0.0
HR3  CPH1 CPH1 CT2      2.0000  2   180.00 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HR3  CPH1 CPH1 CT3      2.0000  2   180.00 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HR3  CPH1 CPH1 HR3      2.0000  2   180.00 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HR3  CPH1 NR1  CPH2     3.0000  2   180.00 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HR3  CPH1 NR1  H        1.0000  2   180.00 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HR3  CPH1 NR2  CPH2     3.0000  2   180.00 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HS   S    CT2  CT3      0.2400  1     0.00 ! ALLOW   ALI SUL ION
                ! ethanethiol C-C-S-H surface, adm jr., 4/18/93
HS   S    CT2  CT3      0.1500  2     0.00 ! ALLOW   ALI SUL ION
                ! ethanethiol C-C-S-H surface, adm jr., 4/18/93
HS   S    CT2  CT3      0.2700  3     0.00 ! ALLOW   ALI SUL ION
                ! ethanethiol C-C-S-H surface, adm jr., 4/18/93
HS   S    CT2  HA2      0.2000  3     0.00 ! ALLOW   ALI SUL ION
                ! methanethiol pure solvent, adm jr., 6/22/92
HS   S    CT3  HA3      0.2000  3     0.00 ! ALLOW   ALI SUL ION
                ! methanethiol pure solvent, adm jr., 6/22/92
N    C    CP1  CP2      0.4000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CP1  CP2      0.6000  2     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CP1  HB1      0.4000  1   180.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CP1  HB1      0.6000  2     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CP1  N        0.3000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CP1  N       -0.3000  4     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CT1  CT1      0.0000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CT1  CT2      0.0000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CT1  CT3      0.0000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CT1  HB1      0.0000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CT2  HB2      0.0000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    C    CT3  HA3      0.0000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
N    CT1  CT2  CA       0.0400  3     0.00 ! ALLOW   ARO
                ! 2.7 kcal/mole CH3 rot in ethylbenzene, adm jr, 3/7/92
NH1  C    CP1  CP2      0.4000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH1  C    CP1  CP2      0.6000  2     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH1  C    CP1  HB1      0.4000  1   180.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH1  C    CP1  HB1      0.6000  2     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH1  C    CP1  N        0.3000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH1  C    CP1  N       -0.3000  4     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH1  C    CT1  CT1      0.0000  1     0.00 !   ALLOW PEP
                ! ala dipeptide corrxn for new C VDW Rmin, 4/10/93 (LK)
NH1  C    CT1  CT2      0.0000  1     0.00 !   ALLOW PEP
                ! ala dipeptide corrxn for new C VDW Rmin, 4/10/93 (LK)
NH1  C    CT1  CT3      0.0000  1     0.00 !   ALLOW PEP
                ! ala dipeptide corrxn for new C VDW Rmin, 4/10/93 (LK)
NH1  C    CT1  HB1      0.0000  1     0.00 !   ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
NH1  C    CT1  NH1      0.6000  1     0.00 !   ALLOW PEP
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93
NH1  C    CT2  CT2      0.0000  1     0.00 !   ALLOW PEP
                ! from NH1  C    CT1  CT2, for lactams, adm jr.
NH1  C    CT2  HA2      0.0000  3     0.00 ! ALLOW PEP
                ! LK for autogenerate dihe, sp2-methyl, no dihedral potential
NH1  C    CT2  HB2      0.0000  1     0.00 !   ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
NH1  C    CT2  NH1      0.6000  1     0.00 !   ALLOW PEP
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93
NH1  C    CT3  HA3      0.0000  3     0.00 ! ALLOW PEP
                ! LK for autogenerate dihe, sp2-methyl, no dihedral potential
NH1  CT1  C    N        0.4000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH1  CT2  C    N        0.4000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH2  CC   CP1  CP2      0.4000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH2  CC   CP1  CP2      0.6000  2     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH2  CC   CP1  HB1      0.4000  1   180.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH2  CC   CP1  HB1      0.6000  2     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH2  CC   CP1  N        0.3000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH2  CC   CP1  N       -0.3000  4     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH2  CC   CT2  HA2      0.0000  3   180.00 ! ALLOW POL
                ! adm jr., 8/13/90 acetamide geometry and vibrations
NH3  CT1  C    N        0.4000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NH3  CT1  C    NH1      0.6000  1     0.00 ! ALLOW PEP PRO
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93
NH3  CT1  CC   NH2      0.4000  1     0.00 ! ALLOW PEP PRO
                ! Alanine dipeptide; NMA; acetate; etc. backbone param. RLD 3/22/92
NH3  CT2  C    N        0.4000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
!!!NH3  CT2  C    NH1      0.4000  1     0.00 ! ALLOW PEP PRO
!!!                ! adm jr. 3/24/92, for PRES GLYP
NH3  CT2  C    NH1      1.0000  1     0.00 ! ALLOW PEP PRO
                ! RB 1/07/11, based on graf et al Gly 3 N-ter J-couplings for PRES GLYP
NH3  CT2  CC   NH2      0.4000  1     0.00 ! ALLOW PEP PRO
                ! Alanine dipeptide; NMA; acetate; etc. backbone param. RLD 3/22/92
NP   CP1  C    N        0.3000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NP   CP1  C    NH1      0.3000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NP   CP1  CC   NH2      0.3000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NR1  CPH1 CPH1 CT2      3.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FROM NR1 CPH1 CPH1 HA
NR1  CPH1 CPH1 CT3      3.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FROM NR1 CPH1 CPH1 HA
NR1  CPH1 CPH1 HR3      3.0000  2   180.00 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
NR1  CPH1 CT2  CT2      0.1900  3     0.00 ! ALLOW ARO
                ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89
NR1  CPH1 CT2  CT3      0.1900  3     0.00 ! ALLOW ARO
                ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89
NR1  CPH1 CT2  HA2      0.1900  3     0.00 ! ALLOW ARO
                ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89
NR1  CPH1 CT3  HA3      0.1900  3     0.00 ! ALLOW ARO
                ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89
NR1  CPH2 NR2  CPH1    14.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
NR2  CPH1 CPH1 CT2      3.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FROM NR2 CPH1 CPH1 HA
NR2  CPH1 CPH1 CT3      3.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/22/89, FROM NR2 CPH1 CPH1 HA
NR2  CPH1 CPH1 HR3      3.0000  2   180.00 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
NR2  CPH1 CPH1 NR1     14.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
!NR2  CPH1 CT2  CT1      0.1900  3     0.00 ! ALLOW ARO
                ! HIS CB-CG TORSION,
NR2  CPH1 CT2  CT2      0.1900  3     0.00 ! ALLOW ARO
                ! HIS CB-CG TORSION,
NR2  CPH1 CT2  CT3      0.1900  3     0.00 ! ALLOW ARO
                ! HIS CB-CG TORSION,
NR2  CPH1 CT2  HA2      0.1900  3     0.00 ! ALLOW ARO
                ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89
NR2  CPH1 CT3  HA3      0.1900  3     0.00 ! ALLOW ARO
                ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89
NR2  CPH2 NR1  CPH1    14.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
NR2  CPH2 NR1  H        1.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
NR3  CPH1 CPH1 CT2      2.5000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
NR3  CPH1 CPH1 CT3      2.5000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
NR3  CPH1 CPH1 HR1      2.5000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
NR3  CPH1 CPH1 NR3     12.0000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
NR3  CPH1 CT2  CT2      0.1900  3     0.00 ! ALLOW ARO
                ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89
NR3  CPH1 CT2  CT3      0.1900  3     0.00 ! ALLOW ARO
                ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89
NR3  CPH1 CT2  HA2      0.1900  3     0.00 ! ALLOW ARO
                ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89
NR3  CPH1 CT3  HA3      0.1900  3     0.00 ! ALLOW ARO
                ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89
NR3  CPH2 NR3  CPH1    12.0000  2   180.00 ! ALLOW ARO
                ! his, ADM JR., 7/20/89
NR3  CPH2 NR3  H        1.4000  2   180.00 ! ALLOW ARO
                ! his, adm jr., 6/27/90
O    C    CP1  CP2      0.4000  1   180.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    CP1  CP2      0.6000  2     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    CP1  HB1      0.4000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    CP1  HB1      0.6000  2     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    CP1  N       -0.3000  4     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    CT1  CT1      1.4000  1     0.00 !   ALLOW PEP
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c
O    C    CT1  CT2      1.4000  1     0.00 !   ALLOW PEP
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c
O    C    CT1  CT3      1.4000  1     0.00 !   ALLOW PEP
                ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c
O    C    CT1  HB1      0.0000  1     0.00 !   ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
O    C    CT1  NH1      0.0000  1     0.00 !   ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
O    C    CT1  NH3      0.0000  1     0.00 ! ALLOW PEP PRO
                ! Backbone parameter set made complete RLD 8/8/90
O    C    CT2  CT2      1.4000  1     0.00 !   ALLOW PEP
                ! from O    C    CT1  CT2, for lactams, adm jr.
O    C    CT2  HA2      0.0000  3   180.00 ! ALLOW POL
                ! adm jr., 8/13/90 acetamide geometry and vibrations
O    C    CT2  HB2      0.0000  1     0.00 !   ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
O    C    CT2  NH1      0.0000  1     0.00 !   ALLOW PEP
                ! Alanine Dipeptide ab initio calc's (LK)
O    C    CT2  NH3      0.0000  1     0.00 ! ALLOW PEP PRO
                ! Backbone parameter set made complete RLD 8/8/90
O    C    CT3  HA3      0.0000  3   180.00 ! ALLOW POL
                ! adm jr., 8/13/90 acetamide geometry and vibrations
O    C    N    CP1      2.7500  2   180.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    N    CP1      0.3000  4     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    N    CP3      2.7500  2   180.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    N    CP3      0.3000  4     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    C    NH1  CT1      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
O    C    NH1  CT2      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
O    C    NH1  CT3      2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
O    C    NH1  H        2.5000  2   180.00 !  ALLOW PEP
                ! Gives appropriate NMA cis/trans barrier. (LK)
O    CC   CP1  CP2      0.4000  1   180.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    CC   CP1  CP2      0.6000  2     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    CC   CP1  HB1      0.4000  1     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    CC   CP1  HB1      0.6000  2     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    CC   CP1  N       -0.3000  4     0.00 ! ALLOW PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
O    CC   CT2  HA2      0.0000  3   180.00 ! ALLOW POL
                ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
O    CC   NH2  H        1.4000  2   180.00 !  ALLOW PEP POL ARO PRO
                ! adm jr. 4/10/91, acetamide update
OB   CD   OS   CT2      0.9650  1   180.00 ! ALLOW   PEP POL
                ! adm jr. 3/19/92, from lipid methyl acetate
OB   CD   OS   CT2      3.8500  2   180.00 ! ALLOW   PEP POL
                ! adm jr. 3/19/92, from lipid methyl acetate
OB   CD   OS   CT3      0.9650  1   180.00 ! ALLOW   PEP POL
                ! adm jr. 3/19/92, from lipid methyl acetate
OB   CD   OS   CT3      3.8500  2   180.00 ! ALLOW   PEP POL
                ! adm jr. 3/19/92, from lipid methyl acetate
OC   CA   CA   CA       3.1000  2   180.00 ! ALLOW   ARO
                ! adm jr. 8/27/91, phenoxide
OC   CA   CA   HP       4.2000  2   180.00 ! ALLOW   ARO
                ! adm jr. 8/27/91, phenoxide
OC   CC   CP1  CP2      0.1600  3     0.00 ! ALLOW PEP PRO POL
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
OC   CC   CP1  HB1      0.1600  3     0.00 ! ALLOW PEP PRO POL
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
OC   CC   CP1  N        0.1600  3     0.00 ! ALLOW PEP PRO POL
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
OC   CC   CP1  NP       0.1600  3     0.00 ! ALLOW PEP PRO POL
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
OC   CC   CT1  NH3      3.2000  2   180.00 ! ALLOW PEP PRO
                ! adm jr. 4/17/94, zwitterionic glycine
OC   CC   CT2  NH3      3.2000  2   180.00 ! ALLOW PEP PRO
                ! adm jr. 4/17/94, zwitterionic glycine
OH1  CA   CA   CA       3.1000  2   180.00 ! ALLOW   ARO
                ! JES 8/25/89 phenol
OH1  CA   CA   HP       4.2000  2   180.00 ! ALLOW   ARO
                ! JES 8/25/89 phenol
S    CT2  CT2  HA2      0.0100  3     0.00 ! ALLOW   ALI SUL ION
                ! DTN 8/24/90
SM   CT2  CT2  HA2      0.0100  3     0.00 ! ALLOW   ALI SUL ION
                ! DTN 8/24/90
SM   SM   CT2  CT1      0.3100  3     0.00 ! ALLOW  SUL ALI
                ! S-S for cys-cys, dummy parameter for now ... DTN  9/04/90
SM   SM   CT2  CT2      0.3100  3     0.00 ! ALLOW  SUL ALI
                ! S-S for cys-cys, dummy parameter for now ... DTN  9/04/90
SM   SM   CT2  CT3      0.3100  3     0.00 ! ALLOW  SUL ALI
                ! S-S for cys-cys, dummy parameter for now ... DTN  9/04/90
SM   SM   CT2  HA2      0.1580  3     0.00 ! ALLOW   ALI SUL ION
                ! expt. dimethyldisulfide,    3/26/92 (FL)
SM   SM   CT3  HA3      0.1580  3     0.00 ! ALLOW   ALI SUL ION
                ! expt. dimethyldisulfide,    3/26/92 (FL)
SS   CS   CT3  HA3      0.1500  3     0.00 ! ALLOW SUL
                ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
X    C    NC2  X        2.2500  2   180.00 ! ALLOW   PEP POL ARO
                ! 9.0->2.25 GUANIDINIUM (KK)
X    CD   OH1  X        2.0500  2   180.00 ! ALLOW   PEP POL ARO ALC
                ! adm jr, 10/17/90, acetic acid C-Oh rotation barrier
X    CD   OS   X        2.0500  2   180.00 ! ALLOW   PEP POL
                ! adm jr. 3/19/92, from lipid methyl acetate
X    CE1  CE1  X        0.1500  1     0.00
                ! 2-butene, adm jr., 2/00 update
X    CE1  CE1  X        8.5000  2   180.00
                ! 2-butene, adm jr., 2/00 update
X    CE2  CE2  X        4.9000  2   180.00 ! 
		! for ethene, yin/adm jr., 12/95
X    CP1  C    X        0.0000  6   180.00 ! ALLOW   POL PEP PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
X    CP1  CC   X        0.0000  6   180.00 ! ALLOW   POL PEP
                ! changed to 0.0 RLD 5/19/92
X    CP1  CD   X        0.0000  6   180.00 ! ALLOW   POL PEP
                ! Alanine dipeptide; NMA; acetate; etc. backbone param. RLD 3/22/92
X    CP1  CP2  X        0.1400  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
X    CP2  CP2  X        0.1600  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
X    CP3  CP2  X        0.1400  3     0.00 ! ALLOW PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
X    CT1  CC   X        0.0500  6   180.00 ! ALLOW   POL PEP
                ! For side chains of asp,asn,glu,gln, (n=6) from KK(LK)
X    CT1  CD   X        0.0000  6   180.00 ! ALLOW   POL PEP
                ! adm jr. 3/19/92, from lipid methyl acetate
X    CT1  CT1  X        0.2000  3     0.00 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
X    CT1  CT2  X        0.2000  3     0.00 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
X    CT1  CT3  X        0.2000  3     0.00 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
X    CT1  NH3  X        0.1000  3     0.00 ! ALLOW   ALI POL
                ! 0.715->0.10 METHYLAMMONIUM (KK)
X    CT1  OH1  X        0.1400  3     0.00 ! ALLOW   ALI ALC ARO
                ! EMB  11/21/89 methanol vib fit
X    CT1  OS   X       -0.1000  3     0.00 ! ALLOW   PEP POL
                ! adm jr. 3/19/92, from lipid methyl acetate
X    CT2  CA   X        0.0000  6     0.00 ! ALLOW   ALI ARO
                ! toluene, adm jr., 3/7/92
X    CT2  CC   X        0.0500  6   180.00 ! ALLOW   POL PEP
                ! For side chains of asp,asn,glu,gln, (n=6) from KK(LK)
X    CT2  CD   X        0.0000  6   180.00 ! ALLOW   POL PEP
                ! adm jr. 3/19/92, from lipid methyl acetate
X    CT2  CT2  X        0.1900  3     0.00 ! ALLOW   ALI
                ! alkane, 4/98, yin and mackerell
X    CT2  CT3  X        0.1600  3     0.00 ! ALLOW   ALI
                ! alkane, 4/98, yin and mackerell
X    CT2  NC2  X        0.0000  6   180.00 ! ALLOW   ALI POL
                ! methylguanidinium, adm jr., 3/26/92
X    CT2  NH3  X        0.1000  3     0.00 ! ALLOW   ALI POL
                ! 0.715->0.10 METHYLAMMONIUM (KK)
X    CT2  OH1  X        0.1400  3     0.00 ! ALLOW   ALI ALC ARO
                ! EMB  11/21/89 methanol vib fit
X    CT2  OS   X       -0.1000  3     0.00 ! ALLOW   PEP POL
                ! adm jr. 3/19/92, from lipid methyl acetate
X    CT3  CA   X        0.0000  6     0.00 ! ALLOW   ALI ARO
                ! toluene, adm jr., 3/7/92
X    CT3  CC   X        0.0500  6   180.00 ! ALLOW   POL PEP
                ! For side chains of asp,asn,glu,gln, (n=6) from KK(LK)
X    CT3  CD   X        0.0000  6   180.00 ! ALLOW   POL PEP
                ! adm jr. 3/19/92, from lipid methyl acetate
X    CT3  CT3  X        0.1525  3     0.00 ! ALLOW   ALI
                ! alkane, 4/98, yin and mackerell
X    CT3  NC2  X        0.0000  6   180.00 ! ALLOW   ALI POL
                ! methylguanidinium, adm jr., 3/26/92
X    CT3  NH2  X        0.1100  3     0.00 ! ALLOW   POL
                ! methylamine geom/freq, adm jr., 6/2/92
X    CT3  NH3  X        0.0900  3     0.00 ! ALLOW   ALI POL
                ! fine-tuned to ab initio; METHYLAMMONIUM, KK 03/10/92
X    CT3  OH1  X        0.1400  3     0.00 ! ALLOW   ALI ALC ARO
                ! EMB  11/21/89 methanol vib fit
X    CT3  OS   X       -0.1000  3     0.00 ! ALLOW   PEP POL
                ! adm jr. 3/19/92, from lipid methyl acetate

!chi1/chi2 fitting, Zhu, 2011
!directly transferred parameters
NH1  CT1  CT1  HA1      0.2000  3     0.00 ! From X    CT1  CT1  X
HB1  CT1  CT1  HA1      0.2000  3     0.00 ! From X    CT1  CT1  X
HB1  CT1  CT1  CT3      0.2000  3     0.00 ! From X    CT1  CT1  X
HA1  CT1  CT1  C        0.2000  3     0.00 ! From X    CT1  CT1  X
!
NH1  CT1  CT2  HA2      0.2000  3     0.00 ! From X    CT1  CT2  X
HB1  CT1  CT2  HA2      0.2000  3     0.00 ! From X    CT1  CT2  X
HB1  CT1  CT2  OH1      0.2000  3     0.00 ! From X    CT1  CT2  X
HB1  CT1  CT2  CT2      0.2000  3     0.00 ! From X    CT1  CT2  X
HA2  CT2  CT1  C        0.2000  3     0.00 ! From X    CT1  CT2  X
HA2  CT2  OH1  H        0.1400  3     0.00 ! From X    CT2  OH1  X      
!
CT1  CT2  CT2  HA2      0.1900  3     0.00 ! From X    CT2  CT2  X
HA2  CT2  CT2  HA2      0.1900  3     0.00 ! From X    CT2  CT2  X
HA2  CT2  CT2  CC       0.1900  3     0.00 ! From X    CT2  CT2  X
!
HB1  CT1  CT2  S        0.2000  3     0.00 ! From X    CT1  CT2  X
!Arg
CT2  CT2  CT2  HA2      0.1900  3     0.00 ! From X    CT2  CT2  X
CT2  CT2  CT2  NC2      0.1900  3     0.00 ! From X    CT2  CT2  X
CT2  CT2  NC2  HC       0.0000  6   180.00 ! From X    CT2  NC2  X
CT2  CT2  NC2  C        0.0000  6   180.00 ! From X    CT2  NC2  X
HA2  CT2  CT2  NC2      0.1900  3     0.00 ! From X    CT2  CT2  X
CT2  NC2  C    NC2      2.2500  2   180.00 ! From X    C    NC2  X
HA2  CT2  NC2  HC       0.0000  6   180.00 ! From X    CT2  NC2  X
HA2  CT2  NC2  C        0.0000  6   180.00 ! From X    CT2  NC2  X
NC2  C    NC2  HC       2.2500  2   180.00 ! From X    C    NC2  X
!Asn
HB1  CT1  CT2  CC       0.2000  3     0.00 ! From X    CT1  CT2  X
!Trp
HB1  CT1  CT2  CY       0.2000  3     0.00 ! From X    CT1  CT2  X
!Asp
HA2  CT2  CC   OC       0.0500  6   180.00 ! From X    CT2  CC   X
!Hsd/Hse
HB1  CT1  CT2  CPH1     0.2000  3     0.00 ! From X    CT1  CT2  X
!Ile,Leu,Val
CT1  CT1  CT3  HA3      0.2000  3     0.00 ! From X    CT1  CT3  X
CT1  CT1  CT2  HA2      0.2000  3     0.00 ! From X    CT1  CT2  X
HB1  CT1  CT1  CT2      0.2000  3     0.00 ! From X    CT1  CT1  X
CT1  CT2  CT3  HA3      0.1600  3     0.00 ! From X    CT2  CT3  X
HA1  CT1  CT3  HA3      0.2000  3     0.00 ! From X    CT1  CT3  X
HA1  CT1  CT2  HA2      0.2000  3     0.00 ! From X    CT1  CT2  X
HA1  CT1  CT2  CT3      0.2000  3     0.00 ! From X    CT1  CT2  X
CT3  CT1  CT2  HA2      0.2000  3     0.00 ! From X    CT1  CT2  X
CT3  CT1  CT2  CT3      0.2000  3     0.00 ! From X    CT1  CT2  X
HA3  CT3  CT1  CT2      0.2000  3     0.00 ! From X    CT1  CT3  X
HA2  CT2  CT3  HA3      0.1600  3     0.00 ! From X    CT2  CT3  X
CT1  CT2  CT1  HA1      0.2000  3     0.00 ! From X    CT1  CT2  X
HB1  CT1  CT2  CT1      0.2000  3     0.00 ! From X    CT1  CT2  X
CT3  CT1  CT3  HA3      0.2000  3     0.00 ! From X    CT1  CT3  X
!Lys
CT2  CT2  CT2  NH3      0.1900  3     0.00 ! From X    CT2  CT2  X
CT2  CT2  NH3  HC       0.1000  3     0.00 ! From X    CT2  NH3  X
HA2  CT2  CT2  NH3      0.1900  3     0.00 ! From X    CT2  CT2  X
HA2  CT2  NH3  HC       0.1000  3     0.00 ! From X    CT2  NH3  X
!Tyr/Phe
HB1  CT1  CT2  CA       0.2000  3     0.00 ! From X    CT1  CT2  X
HA2  CT2  CA   CA       0.0000  6     0.00 ! From X    CT2  CA   X
!Thr
HB1  CT1  CT1  OH1      0.2000  3     0.00 ! From X    CT1  CT1  X
HA1  CT1  OH1  H        0.1400  3     0.00 ! From X    CT1  OH1  X
OH1  CT1  CT3  HA3      0.2000  3     0.00 ! From X    CT1  CT3  X
!Gln
CT2  CT2  CC   O        0.0500  6   180.00 ! From X    CT2  CC   X
CT2  CT2  CC   NH2      0.0500  6   180.00 ! From X    CT2  CC   X
!Glu
CT2  CT2  CC   OC       0.0500  6   180.00 ! From X    CT2  CC   X
!Glu/Hsp
NH1  CT1  CT2A HA2      0.2000  3     0.00 ! From X    CT1  CT2  X
NH3  CT1  CT2A CT2      0.2000  3     0.00 ! From X    CT1  CT2  X !N terminus
CT1  CT2A CT2  HA2      0.1900  3     0.00 ! From X    CT2  CT2  X
HB1  CT1  CT2A HA2      0.2000  3     0.00 ! From X    CT1  CT2  X
HB1  CT1  CT2A CT2      0.2000  3     0.00 ! From X    CT1  CT2  X
HA2  CT2A CT1  C        0.2000  3     0.00 ! From X    CT1  CT2  X
HA2  CT2A CT1  CC       0.2000  3     0.00 ! RB: added for C-ter Glu
HA2  CT2A CT2  HA2      0.1900  3     0.00 ! From X    CT2  CT2  X
HA2  CT2A CT2  CC       0.1900  3     0.00 ! From X    CT2  CT2  X
HB1  CT1  CT2A CPH1     0.2000  3     0.00 ! From X    CT1  CT2  X
C    NH1  CT1  CT2A     1.8000  1     0.00 ! from CT2  CT1  NH1  C
H    NH1  CT1  CT2A     0.0000  1     0.00 ! from H    NH1  CT1  CT2
CT2A CT1  C    O        1.4000  1     0.00 ! from O    C    CT1  CT2
CT2A CT1  C    NH1      0.0000  1     0.00 ! NH1  C    CT1  CT2
CT2A CT1  C    N        0.0000  1     0.00 ! RB: added for GLU-PRO in UBQ
! Glup
CT1  CT2A CT2  CD       0.1900  3     0.00 ! From X    CT2  CT2  X
HA2  CT2A CT2  CD       0.1900  3     0.00 ! From X    CT2  CT2  X
CT2A CPH1 CPH1 HR1      1.0000  2   180.00 ! from HR1  CPH1 CPH1 CT2
CT2A CPH1 CPH1 NR3      2.5000  2   180.00 ! from NR3  CPH1 CPH1 CT2
CT2A CPH1 NR3  H        3.0000  2   180.00 ! from H    NR3  CPH1 CT2
CT2A CPH1 NR3  CPH2     2.5000  2   180.00 ! from CT2  CPH1 NR3  CPH2
HA2  CT2A CPH1 CPH1     0.0000  3     0.00 ! from HA2  CT2  CPH1 CPH1
HA2  CT2A CPH1 NR3      0.1900  3     0.00 ! from NR3  CPH1 CT2  HA2

! Fit dihedrals
! Variable cutoff based on QM and weighted in favor of alphaR and EXT (5:5:1)
! Shared dihedrals were fitted simultaneously

! Group-fitted for Lys/Arg/Gln/Met
C    CT1  CT2  CT2      0.3500  1   180.00 
C    CT1  CT2  CT2      0.4200  2   180.00 
C    CT1  CT2  CT2      1.9100  3   180.00 
CT2  CT2  CT1  NH1      0.8800  1   180.00 
CT2  CT2  CT1  NH1      0.0000  2   180.00 
CT2  CT2  CT1  NH1      1.9000  3     0.00 
CC   CT2  CT2  CT1      1.8400  1   180.00 
CC   CT2  CT2  CT1      0.8400  2   180.00 
CC   CT2  CT2  CT1      0.3900  3   180.00 
CT1  CT2  CT2  CT2      0.6300  1   180.00 
CT1  CT2  CT2  CT2      0.0100  2     0.00 
CT1  CT2  CT2  CT2      0.1500  3     0.00 
CT1  CT2  CT2  S        0.1400  1   180.00 
CT1  CT2  CT2  S        0.5400  2     0.00 
CT1  CT2  CT2  S        0.6900  3     0.00 
! Fitted Asn 
C    CT1  CT2  CC       1.4100  1   180.00 
C    CT1  CT2  CC       1.2900  2   180.00 
C    CT1  CT2  CC       0.5900  3   180.00 
CC   CT2  CT1  NH1      0.2800  1   180.00 
CC   CT2  CT1  NH1      0.5000  2   180.00 
CC   CT2  CT1  NH1      0.3800  3     0.00 
CT1  CT2  CC   NH2      0.6200  1   180.00 
CT1  CT2  CC   NH2      0.6600  2   180.00 
CT1  CT2  CC   NH2      0.7200  3   180.00 
CT1  CT2  CC   O        0.4200  1   180.00 
CT1  CT2  CC   O        0.1500  2   180.00 
CT1  CT2  CC   O        0.9500  3   180.00 
! Fitted Asp
C    CT1  CT2A CC       1.6100  1   180.00 
C    CT1  CT2A CC       1.2900  2   180.00 
C    CT1  CT2A CC       0.5900  3   180.00 
CC   CT2A CT1  NH1      0.6800  1   180.00 
CC   CT2A CT1  NH1      0.1000  2   180.00 
CC   CT2A CT1  NH1      0.3800  3     0.00 
CT1  CT2A CC   OC       0.8400  1     0.00 
CT1  CT2A CC   OC       0.9800  2   180.00 
CT1  CT2A CC   OC       1.4600  3     0.00 
! Fitted Cys
CT1  CT2  S    HS       0.2000  1     0.00
CT1  CT2  S    HS       0.6500  2     0.00
CT1  CT2  S    HS       0.2200  3     0.00
C    CT1  CT2  S        0.2400  1   180.00
C    CT1  CT2  S        0.7500  2   180.00
C    CT1  CT2  S        1.3500  3   180.00
NH1  CT1  CT2  S        0.3400  1     0.00
NH1  CT1  CT2  S        0.5000  2   180.00
NH1  CT1  CT2  S        1.4300  3     0.00
! Fitted Glu
CC   CT2  CT2A CT1      0.0000  1   180.00
CC   CT2  CT2A CT1      0.3800  2   180.00
CC   CT2  CT2A CT1      0.5900  3   180.00
C    CT1  CT2A CT2      0.1100  1     0.00
C    CT1  CT2A CT2      0.9800  2   180.00
C    CT1  CT2A CT2      1.6000  3   180.00
CC   CT1  CT2A CT2      1.6000  3   180.00
CT2  CT2A CT1  NH1      0.3000  1     0.00 
CT2  CT2A CT1  NH1      0.3500  2     0.00
CT2  CT2A CT1  NH1      1.7600  3     0.00
! Group-fitted for Hsd/Hse
CPH1 CPH1 CT2  CT1      1.7400  1     0.00
CPH1 CPH1 CT2  CT1      0.1500  2     0.00
CPH1 CPH1 CT2  CT1      0.7700  3   180.00
CT1  CT2  CPH1 NR1      1.4900  1     0.00
CT1  CT2  CPH1 NR1      0.0900  2   180.00
CT1  CT2  CPH1 NR1      0.7900  3   180.00
CT1  CT2  CPH1 NR2      1.0900  1     0.00
CT1  CT2  CPH1 NR2      0.0900  2     0.00
CT1  CT2  CPH1 NR2      0.6700  3   180.00
C    CT1  CT2  CPH1     0.1800  1   180.00
C    CT1  CT2  CPH1     0.6400  2   180.00
C    CT1  CT2  CPH1     0.8700  3   180.00
CPH1 CT2  CT1  NH1      0.0000  1     0.00
CPH1 CT2  CT1  NH1      0.0000  2   180.00
CPH1 CT2  CT1  NH1      0.9000  3     0.00
! Fitted Hsp
CPH1 CPH1 CT2A CT1      2.0400  1     0.00
CPH1 CPH1 CT2A CT1      0.4400  2     0.00
CPH1 CPH1 CT2A CT1      0.1300  3   180.00
CT1  CT2A CPH1 NR3      0.5300  1   180.00
CT1  CT2A CPH1 NR3      0.4200  2   180.00
CT1  CT2A CPH1 NR3      0.3000  3   180.00
C    CT1  CT2A CPH1     1.7500  1   180.00
C    CT1  CT2A CPH1     0.1300  2     0.00
C    CT1  CT2A CPH1     1.8600  3   180.00
CPH1 CT2A CT1  NH1      1.0900  1   180.00
CPH1 CT2A CT1  NH1      0.2200  2   180.00
CPH1 CT2A CT1  NH1      2.3200  3     0.00
! Group-fitted for Ile/Thr
CT1  CT1  CT2  CT3      0.3800  1   180.00
CT1  CT1  CT2  CT3      0.1300  2   180.00
CT1  CT1  CT2  CT3      0.2900  3   180.00
C    CT1  CT1  CT2      0.1000  1   180.00
C    CT1  CT1  CT2      0.5200  2   180.00
C    CT1  CT1  CT2      0.2900  3   180.00
CT2  CT1  CT1  NH1      0.1200  1   180.00
CT2  CT1  CT1  NH1      0.3600  2   180.00
CT2  CT1  CT1  NH1      0.4100  3     0.00
! Fitted Leu 
CT1  CT2  CT1  CT3      0.0500  1     0.00
CT1  CT2  CT1  CT3      0.1000  2   180.00
CT1  CT2  CT1  CT3      0.0100  3   180.00
C    CT1  CT2  CT1      0.3200  1   180.00
C    CT1  CT2  CT1      0.6100  2   180.00
C    CT1  CT2  CT1      0.7200  3   180.00
CT1  CT2  CT1  NH1      0.4800  1   180.00
CT1  CT2  CT1  NH1      0.4200  2   180.00
CT1  CT2  CT1  NH1      0.6500  3     0.00
! Group-fitted for Phe/Tyr
CA   CA   CT2  CT1      1.0700  1     0.00
CA   CA   CT2  CT1      0.2400  2   180.00
CA   CA   CT2  CT1      0.1700  3   180.00
C    CT1  CT2  CA       1.2800  1   180.00
C    CT1  CT2  CA       0.9400  2   180.00
C    CT1  CT2  CA       1.5700  3   180.00
CA   CT2  CT1  NH1      0.5200  1   180.00
CA   CT2  CT1  NH1      0.6200  2   180.00
CA   CT2  CT1  NH1      1.5800  3     0.00
! Fitted Ser
CT1  CT2  OH1  H        0.0200  1     0.00
CT1  CT2  OH1  H        0.5600  2     0.00
CT1  CT2  OH1  H        0.4900  3     0.00
C    CT1  CT2  OH1      0.6500  1   180.00
C    CT1  CT2  OH1      0.2500  2   180.00
C    CT1  CT2  OH1      1.1700  3   180.00
NH1  CT1  CT2  OH1      0.1800  1   180.00
NH1  CT1  CT2  OH1      0.1900  2   180.00
NH1  CT1  CT2  OH1      1.4600  3     0.00
! Group-fitted for Ile/Thr
CT1  CT1  OH1  H        0.1800  1     0.00
CT1  CT1  OH1  H        0.0600  2     0.00
CT1  CT1  OH1  H        0.2500  3     0.00
C    CT1  CT1  OH1      0.7900  1   180.00
C    CT1  CT1  OH1      0.3900  2   180.00
C    CT1  CT1  OH1      0.9900  3   180.00
NH1  CT1  CT1  OH1      0.0900  1     0.00
NH1  CT1  CT1  OH1      0.1900  2   180.00
NH1  CT1  CT1  OH1      0.1700  3     0.00
! Fitted Trp
CA   CY   CT2  CT1      0.0300  1     0.00
CA   CY   CT2  CT1      0.5500  2     0.00
CA   CY   CT2  CT1      0.3900  3   180.00
CPT  CY   CT2  CT1      0.3600  1   180.00
CPT  CY   CT2  CT1      0.0500  2     0.00
CPT  CY   CT2  CT1      0.1900  3   180.00
C    CT1  CT2  CY       1.0900  1   180.00
C    CT1  CT2  CY       0.5000  2   180.00
C    CT1  CT2  CY       1.1700  3   180.00
CY   CT2  CT1  NH1      0.2900  1   180.00
CY   CT2  CT1  NH1      0.6600  2   180.00
CY   CT2  CT1  NH1      1.1700  3     0.00
! Fitted Val
C    CT1  CT1  CT3      0.1400  1   180.00
C    CT1  CT1  CT3      0.2600  2   180.00
C    CT1  CT1  CT3      0.3300  3   180.00
CT3  CT1  CT1  NH1      0.1800  1     0.00
CT3  CT1  CT1  NH1      0.0600  2     0.00
CT3  CT1  CT1  NH1      0.5900  3     0.00
!ASP, CT2->CT2A, jshim
H    NH1  CT2A CC       0.0000  1     0.00
X    CT2A CC   X        0.0500  6   180.00
HB1  CT1  CT2A CC       0.2000  3     0.00
HA2  CT2A CC   OC       0.0500  6   180.00
NH3  CT1  CT2A HA2      0.2000  3     0.00
NH3  CT1  CT2A CC       0.2000  3     0.00
CC   CT2A CT1  CC       0.2000  3     0.00
!termini specific terms
CPH1 CT2A CT1  CC       0.2000  3     0.00
CPH1 CT2A CT1  NH3      0.2000  3     0.00
CPH1 CT2A CT1  CD       0.2000  3     0.00 
HA2  CT2A CT1  CD       0.2000  3     0.00     
CT2  CT2A CT1  CD       0.2000  3     0.00
! RESI CYSM and PRES CYSD
H    NH2  CT1  CS       0.1100  3     0.00 ! from H    NH2  CT1  CT2 or H    NH2  CT1  CT2 , kevo
CS   CT1  NH1  C        1.8000  1     0.00 ! from CT2  CT1  NH1  C   or CT2A CT1  NH1  C , kevo
H    NH1  CT1  CS       0.0000  1     0.00 ! from H    NH1  CT1  CT2 or H    NH1  CT1  CT2 , kevo
N    C    CT1  CS       0.0000  1     0.00 ! from N    C    CT1  CT2 or N    C    CT1  CT2 , kevo
NH1  C    CT1  CS       0.0000  1     0.00 ! from NH1  C    CT1  CT2 or NH1  C    CT1  CT2 , kevo
O    C    CT1  CS       1.4000  1     0.00 ! from O    C    CT1  CT2 or O    C    CT1  CT2 , kevo
HA2  CS   CT1  C        0.2000  3     0.00 ! from HA2  CT2  CT1  C   or HA2  CT2A CT1  C , kevo
NH1  CT1  CS   HA2      0.2000  3     0.00 ! from NH1  CT1  CT2  HA2 or NH1  CT1  CT2A HA2 , kevo
HB1  CT1  CS   HA2      0.2000  3     0.00 ! from HB1  CT1  CT2  HA2 or HB1  CT1  CT2A HA2 , kevo
HB1  CT1  CS   SS       0.2000  3     0.00 ! from HB1  CT1  CT2  S   or HB1  CT1  CT2A S , kevo
C    CT1  CS   SS       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
NH1  CT1  CS   SS       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
! Termini 
NH3  CT1  CS   HA2      0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
NH3  CT1  CS   SS       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
NH2  CT1  CS   HA2      0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
NH2  CT1  CS   SS       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
CC   CT1  CS   HA2      0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
CC   CT1  CS   SS       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
CD   CT1  CS   HA2      0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
CD   CT1  CS   SS       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
! PRES SERD
NH1  CT1  CT2  OC       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
NH2  CT1  CT2  OC       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
NH3  CT1  CT2  OC       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
C    CT1  CT2  OC       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
CC   CT1  CT2  OC       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
CD   CT1  CT2  OC       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
HB1  CT1  CT2  OC       0.2000  3     0.00 ! from X    CT1  CT2  X , kevo
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

temp1()