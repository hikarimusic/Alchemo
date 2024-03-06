def temp1():
    # Given string
    data = """
    CA   CA    305.000     1.3750 ! ALLOW   ARO
                    ! benzene, JES 8/25/89
    CE1  CE1   440.000     1.3400   ! 
            ! for butene; from propene, yin/adm jr., 12/95
    CE1  CE2   500.000     1.3420   ! 
            ! for propene, yin/adm jr., 12/95
    CE1  CT2   365.000     1.5020   ! 
            ! for butene; from propene, yin/adm jr., 12/95
    CE1  CT3   383.000     1.5040   ! 
            ! for butene, yin/adm jr., 12/95
    CE2  CE2   510.000     1.3300   ! 
            ! for ethene, yin/adm jr., 12/95
    CP1  C     250.000     1.4900 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    CP1  CC    250.000     1.4900 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    CP1  CD    200.000     1.4900 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    CP2  CP1   222.500     1.5270 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    CP2  CP2   222.500     1.5370 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    CP3  CP2   222.500     1.5370 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    CPH1 CPH1  410.000     1.3600 ! ALLOW ARO
                    ! histidine, adm jr., 6/27/90
    CT1  C     250.000     1.4900 ! ALLOW   ALI PEP POL ARO
                    ! Ala Dipeptide ab initio calc's (LK) fixed from 10/90 (5/91)
    CT1  CC    200.000     1.5220 ! ALLOW   POL
                    ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
    CT1  CD    200.000     1.5220 ! ALLOW   POL
                    ! adm jr. 5/02/91, acetic acid pure solvent
    CT1  CT1   222.500     1.5000 ! ALLOW   ALI
                    ! alkane update, adm jr., 3/2/92
    CT2  C     250.000     1.4900 ! ALLOW   ALI PEP POL ARO
                    ! Ala Dipeptide ab initio calc's (LK) fixed from 10/90 (5/91)
    CT2  CA    230.000     1.4900 ! ALLOW   ALI ARO
                    ! phe,tyr, JES 8/25/89
    CT2  CC    200.000     1.5220 ! ALLOW   POL
                    ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
    CT2  CD    200.000     1.5220 ! ALLOW   POL
                    ! adm jr. 5/02/91, acetic acid pure solvent
    CT2  CPH1  229.630     1.5000 ! ALLOW ARO
                    ! his, adm jr., 7/22/89, FC from CT2CT, BL from crystals
    CT2  CT1   222.500     1.5380 ! ALLOW   ALI
                    ! alkane update, adm jr., 3/2/92
    CT2  CT2   222.500     1.5300 ! ALLOW   ALI
                    ! alkane update, adm jr., 3/2/92
    CT3  C     250.000     1.4900 ! ALLOW   ALI PEP POL ARO
                    ! Ala Dipeptide ab initio calc's (LK) fixed from 10/90 (5/91)
    CT3  CA    230.000     1.4900 ! ALLOW   ALI ARO
                    ! toluene, adm jr. 3/7/92
    CT3  CC    200.000     1.5220 ! ALLOW   POL
                    ! adm jr. 4/05/91, for asn,asp,gln,glu and cters
    CT3  CD    200.000     1.5220 ! ALLOW   POL
                    ! adm jr. 5/02/91, acetic acid pure solvent
    CT3  CPH1  229.630     1.5000 ! ALLOW ARO
                    ! his, adm jr., 7/22/89, FC from CT2CT, BL from crystals
    CT3  CS    190.000     1.5310 ! ALLOW   SUL
                    ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
    CT3  CT1   222.500     1.5380 ! ALLOW   ALI
                    ! alkane update, adm jr., 3/2/92
    CT3  CT2   222.500     1.5280 ! ALLOW   ALI
                    ! alkane update, adm jr., 3/2/92
    CT3  CT3   222.500     1.5300 ! ALLOW   ALI
                    ! alkane update, adm jr., 3/2/92
    H    CD    330.000     1.1100 ! ALLOW   PEP POL ARO
                    ! adm jr. 5/02/91, acetic acid pure solvent
    HA1  CC    317.130     1.1000 ! ALLOW POL
                    ! adm jr., 5/13/91, formamide geometry and vibrations
    HA2  CP2   309.000     1.1110 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    HA2  CP3   309.000     1.1110 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    HA2  CS    300.000     1.1110 ! ALLOW   SUL
                    ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
    HA3  CS    300.000     1.1110 ! ALLOW   SUL
                    ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
    HA1  CT1   309.000     1.1110 ! ALLOW   ALI
                    ! alkane update, adm jr., 3/2/92
    HA2  CT2   309.000     1.1110 ! ALLOW   ALI
                    ! alkane update, adm jr., 3/2/92
    HA3  CT3   322.000     1.1110 ! ALLOW   ALI
                    ! alkane update, adm jr., 3/2/92
    !HA   CY    330.000     1.0800 ! ALLOW   ARO
                    ! JWK 05/14/91 new r0 from indole
    HE1  CE1   360.500     1.1000   ! 
            ! for propene, yin/adm jr., 12/95
    HE2  CE2   365.000     1.1000   ! 
            ! for ethene, yin/adm jr., 12/95
    HB1  CP1   330.000     1.0800 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    HB1  CT1   330.000     1.0800 ! ALLOW   PEP
                    ! Alanine Dipeptide ab initio calc's (LK)
    HB2  CT2   330.000     1.0800 ! ALLOW   PEP
                    ! Alanine Dipeptide ab initio calc's (LK)
    !HB3  CT3   330.000     1.0800 ! ALLOW   PEP
                    ! Alanine Dipeptide ab initio calc's (LK)
    HP   CA    340.000     1.0800 ! ALLOW   ARO
                    ! phe,tyr JES 8/25/89
    HR1  CPH1  375.000     1.0830 ! ALLOW ARO
                    ! his, adm jr., 6/27/90
    HR1  CPH2  340.000     1.0900 ! ALLOW ARO
                    ! his, adm jr., 6/28/29
    HR2  CPH2  333.000     1.0700 ! ALLOW ARO
                    ! his, adm jr., 6/27/90
    HR3  CPH1  365.000     1.0830 ! ALLOW ARO
                    ! adm jr., 3/24/92, maintain old aliphatic H VDW params
    N    C     260.000     1.3000 ! ALLOW PEP POL ARO PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    N    CP1   320.000     1.4340 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    N    CP3   320.000     1.4550 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    NC2  C     450.000     1.3650 ! ALLOW   PEP POL ARO
                    ! mp2/6-311g** mgua vib. data, adm jr., 1/04
    NC2  CT2   390.000     1.4900 ! ALLOW   ALI POL
                    ! mp2/6-311g** mgua vib. data, adm jr., 1/04
    NC2  CT3   390.000     1.4900 ! ALLOW   ALI POL
                    ! mp2/6-311g** mgua vib. data, adm jr., 1/04
    NC2  HC    455.000     1.0000 ! ALLOW   POL
                    ! 405.0->455.0 GUANIDINIUM (KK)
    NH1  C     370.000     1.3450 ! ALLOW   PEP POL ARO
                    ! Alanine Dipeptide ab initio calc's (LK)
    NH1  CT1   320.000     1.4300 ! ALLOW   ALI PEP POL ARO
                    ! NMA Gas & Liquid Phase IR Spectra (LK)
    NH1  CT2   320.000     1.4300 ! ALLOW   ALI PEP POL ARO
                    ! NMA Gas & Liquid Phase IR Spectra (LK)
    NH1  CT3   320.000     1.4300 ! ALLOW   ALI PEP POL ARO
                    ! NMA Gas & Liquid Phase IR Spectra (LK)
    NH1  H     440.000     0.9970 ! ALLOW   PEP POL ARO
                    ! Alanine Dipeptide ab initio calc's (LK)
    NH1  HC    405.000     0.9800 ! ALLOW   PEP POL ARO
                    ! (DS)
    NH2  CC    430.000     1.3600 ! ALLOW   PEP POL ARO
                    ! adm jr. 4/10/91, acetamide
    NH2  CT2   240.000     1.4550
                    ! from NH2  CT3, neutral glycine, adm jr.
    NH2  CT3   240.000     1.4550 ! ALLOW   POL
                    ! methylamine geom/freq, adm jr., 6/2/92
    NH2  H     480.000     1.0000 ! ALLOW   POL
                    ! adm jr. 8/13/90 acetamide geometry and vibrations
    NH2  HC    460.000     1.0000 ! ALLOW   POL
                    ! methylamine geom/freq, adm jr., 6/2/92
    NH3  CT1   200.000     1.4800 ! ALLOW   ALI POL
                    ! new stretch and bend; methylammonium (KK 03/10/92)
    NH3  CT2   200.000     1.4800 ! ALLOW   ALI POL
                    ! new stretch and bend; methylammonium (KK 03/10/92)
    NH3  CT3   200.000     1.4800 ! ALLOW   ALI POL
                    ! new stretch and bend; methylammonium (KK 03/10/92)
    NH3  HC    403.000     1.0400 ! ALLOW   POL
                    ! new stretch and bend; methylammonium (KK 03/10/92)
    NP   CP1   320.000     1.4850 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    NP   CP3   320.000     1.5020 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    NP   HC    460.000     1.0060 ! ALLOW PRO
                    ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
    NR1  CPH1  400.000     1.3800 ! ALLOW ARO
                    ! his, ADM JR., 7/20/89
    NR1  CPH2  400.000     1.3600 ! ALLOW ARO
                    ! his, ADM JR., 7/20/89
    NR1  H     466.000     1.0000 ! ALLOW ARO
                    ! his, ADM JR., 7/20/89
    NR2  CPH1  400.000     1.3800 ! ALLOW ARO
                    ! his, ADM JR., 7/20/89
    NR2  CPH2  400.000     1.3200 ! ALLOW ARO
                    ! his, ADM JR., 7/20/89
    NR3  CPH1  380.000     1.3700 ! ALLOW ARO
                    ! his, adm jr., 6/28/90
    NR3  CPH2  380.000     1.3200 ! ALLOW ARO
                    ! his, adm jr., 6/27/90
    NR3  H     453.000     1.0000 ! ALLOW ARO
                    ! his, adm jr., 6/27/90
    O    C     620.000     1.2300 ! ALLOW   PEP POL ARO
                    ! Peptide geometry, condensed phase (LK)
    O    CC    650.000     1.2300 ! ALLOW   PEP POL ARO
                    ! adm jr. 4/10/91, acetamide
    OB   CC    750.000     1.2200 ! ALLOW   PEP POL ARO
                    ! adm jr., 10/17/90, acetic acid vibrations and geom.
    OB   CD    750.000     1.2200 ! ALLOW   PEP POL ARO
                    ! adm jr. 5/02/91, acetic acid pure solvent
    OC   CA    525.000     1.2600 ! ALLOW   PEP POL ARO ION
                    ! adm jr. 8/27/91, phenoxide
    OC   CC    525.000     1.2600 ! ALLOW   PEP POL ARO ION
                    ! adm jr. 7/23/91, acetic acid
    OC   CT2   450.000     1.3300 ! ALLOW   ALC
                    ! ethoxide 6-31+G* geom/freq, adm jr., 6/1/92
    OC   CT3   450.000     1.3300 ! ALLOW   ALC
                    ! methoxide 6-31+G* geom/freq, adm jr., 6/1/92
    OH1  CA    334.300     1.4110 ! ALLOW   ARO ALC
                    ! MeOH, EMB 10/10/89,
    OH1  CD    230.000     1.4000 ! ALLOW   PEP POL ARO ALC
                    ! adm jr. 5/02/91, acetic acid pure solvent
    OH1  CT1   428.000     1.4200 ! ALLOW   ALI ALC ARO
                    ! methanol vib fit EMB 11/21/89
    OH1  CT2   428.000     1.4200 ! ALLOW   ALI ALC ARO
                    ! methanol vib fit EMB 11/21/89
    OH1  CT3   428.000     1.4200 ! ALLOW   ALI ALC ARO
                    ! methanol vib fit EMB 11/21/89
    OH1  H     545.000     0.9600 ! ALLOW   ALC ARO
                    ! EMB 11/21/89 methanol vib fit
    OS   CD    150.000     1.3340 ! ALLOW POL PEP
                    ! adm jr. 5/02/91, acetic acid pure solvent
    OS   CT3   340.000     1.4300 ! ALLOW POL PEP
                    ! adm jr., 4/05/91, for PRES CT1 from methylacetate
    S    CT2   198.000     1.8180 ! ALLOW   ALI SUL ION
                    ! fitted to C-S s   9/26/92 (FL)
    S    CT3   240.000     1.8160 ! ALLOW   ALI SUL ION
                    ! fitted to C-S s   9/26/92 (FL)
    S    HS    275.000     1.3250 ! ALLOW   SUL ION
                    ! methanethiol pure solvent, adm jr., 6/22/92
    SM   CT2   214.000     1.8160 ! ALLOW   SUL ION
                    ! improved CSSC surface in DMDS  5/15/92 (FL)
    SM   CT3   214.000     1.8160 ! ALLOW   SUL ION
                    ! improved CSSC surface in DMDS  5/15/92 (FL)
    SM   SM    173.000     2.0290 ! ALLOW   SUL ION
                    ! improved CSSC surface in DMDS  5/15/92 (FL)
    SS   CS    205.000     1.8360 ! ALLOW   SUL
                    ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92
    """

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

data = """
C      0.000000  -0.110000     2.000000
CA     0.000000  -0.070000     1.992400
CC     0.000000  -0.070000     2.000000
CD     0.000000  -0.070000     2.000000
CE1    0.000000  -0.068000     2.090000
CE2    0.000000  -0.064000     2.080000
CP1    0.000000  -0.020000     2.275000   0.000000  -0.010000     1.900000
CP2    0.000000  -0.055000     2.175000   0.000000  -0.010000     1.900000
CP3    0.000000  -0.055000     2.175000   0.000000  -0.010000     1.900000
CPH1   0.000000  -0.050000     1.800000
CPH2   0.000000  -0.050000     1.800000
CS     0.000000  -0.110000     2.200000
CPT    0.000000  -0.099000     1.860000
CY     0.000000  -0.073000     1.990000
CAI    0.000000  -0.073000     1.990000
CT     0.000000  -0.020000     2.275000 0.0 -0.01 1.9
CT1    0.000000  -0.032000     2.000000 0.0 -0.01 1.9
CT2    0.000000  -0.056000     2.010000 0.0 -0.01 1.9
CT2A   0.000000  -0.056000     2.010000 0.0 -0.01 1.9
CT3    0.000000  -0.078000     2.040000 0.0 -0.01 1.9
C3     0.000000  -0.020000     2.275000
H      0.000000  -0.046000     0.224500
HA     0.000000  -0.022000     1.320000
HB1    0.000000  -0.022000     1.320000
HB2    0.000000  -0.028000     1.340000
HE1    0.000000  -0.031000     1.250000
HE2    0.000000  -0.026000     1.260000
HC     0.000000  -0.046000     0.224500
HP     0.000000  -0.030000     1.358200   0.000000  -0.030000     1.358200
HR1    0.000000  -0.046000     0.900000
HR2    0.000000  -0.046000     0.700000
HR3    0.000000  -0.007800     1.468000
HS     0.000000  -0.100000     0.450000
HA1    0.000000  -0.045000     1.340000
HA2    0.000000  -0.034000     1.340000
HA3    0.000000  -0.024000     1.340000
N      0.000000  -0.200000     1.850000   0.000000  -0.000100     1.850000
NC2    0.000000  -0.200000     1.850000
NH1    0.000000  -0.200000     1.850000   0.000000  -0.200000     1.550000
NH2    0.000000  -0.200000     1.850000
NH3    0.000000  -0.200000     1.850000
NP     0.000000  -0.200000     1.850000
NR1    0.000000  -0.200000     1.850000
NR2    0.000000  -0.200000     1.850000
NR3    0.000000  -0.200000     1.850000
NY     0.000000  -0.200000     1.850000
O      0.000000  -0.120000     1.700000   0.000000  -0.120000     1.400000
OB     0.000000  -0.120000     1.700000   0.000000  -0.120000     1.400000
OC     0.000000  -0.120000     1.700000
OH1    0.000000  -0.152100     1.770000
OS     0.000000  -0.152100     1.770000
S      0.000000  -0.450000     2.000000
SM     0.000000  -0.380000     1.975000
SS     0.000000  -0.470000     2.200000
"""

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

temp2()