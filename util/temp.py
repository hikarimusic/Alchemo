data = '''
!
!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
!
!epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
!Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
!
!atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4
!
!carbons
C      0.000000  -0.110000     2.000000 ! ALLOW   PEP POL ARO
                ! NMA pure solvent, adm jr., 3/3/93
CA     0.000000  -0.070000     1.992400 ! ALLOW   ARO
                ! benzene (JES)
CC     0.000000  -0.070000     2.000000 ! ALLOW   PEP POL ARO
                ! adm jr. 3/3/92, acetic acid heat of solvation
CD     0.000000  -0.070000     2.000000 ! ALLOW  POL
                ! adm jr. 3/19/92, acetate a.i. and dH of solvation
CE1    0.000000  -0.068000     2.090000 ! 
		! for propene, yin/adm jr., 12/95
CE2    0.000000  -0.064000     2.080000 ! 
		! for ethene, yin/adm jr., 12/95
CP1    0.000000  -0.020000     2.275000   0.000000  -0.010000     1.900000 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
CP2    0.000000  -0.055000     2.175000   0.000000  -0.010000     1.900000 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
CP3    0.000000  -0.055000     2.175000   0.000000  -0.010000     1.900000 ! ALLOW   ALI
                ! alkane update, adm jr., 3/2/92
CPH1   0.000000  -0.050000     1.800000 ! ALLOW ARO
                ! adm jr., 10/23/91, imidazole solvation and sublimation
CPH2   0.000000  -0.050000     1.800000 ! ALLOW ARO
                ! adm jr., 10/23/91, imidazole solvation and sublimation
CS     0.000000  -0.110000     2.200000 ! ALLOW SUL
                ! methylthiolate to water and F.E. of solvation, adm jr. 6/1/92
CPT    0.000000  -0.099000     1.860000 ! atm, indole vaporization 5/05
CY     0.000000  -0.073000     1.990000 ! atm, indole vaporization 5/05
CAI    0.000000  -0.073000     1.990000 ! atm, indole vaporization 5/05
                ! TRP, JWK 08/29/89
!new alkanes atoms types for conversion to new LJ parameters for c27
CT       0.0       -0.0200    2.275 0.0 -0.01 1.9 ! 
CT1      0.0       -0.0320    2.000 0.0 -0.01 1.9 ! alkane, 4/07, viv and adm jr.
CT2      0.0       -0.0560    2.010 0.0 -0.01 1.9 ! alkane, 4/98, yin, adm jr.
CT2A     0.0       -0.0560    2.010 0.0 -0.01 1.9 ! from CT2 (GLU, HSP), 05282010, zhu
CT3      0.0       -0.0780    2.040 0.0 -0.01 1.9 ! alkane, 4/98, yin, adm jr.
!
C3     0.000000  -0.020000     2.275000 ! cyclopropane JMW  16 april 04
! hydrogens
H      0.000000  -0.046000     0.224500 ! ALLOW PEP POL SUL ARO ALC
                ! same as TIP3P hydrogen, adm jr., 7/20/89
HA     0.000000  -0.022000     1.320000 ! ALLOW PEP ALI POL SUL ARO PRO ALC
                ! methane/ethane a.i. and ethane pure solvent, adm jr, 2/3/92
HB1    0.000000  -0.022000     1.320000 ! 
                ! methane/ethane a.i. and ethane pure solvent, adm jr, 2/3/92
HB2    0.000000  -0.028000     1.340000 ! 
                ! Yin and MacKerell, adm jr., 5/30/02
HE1    0.000000  -0.031000     1.250000 ! 
		! for propene, yin/adm jr., 12/95
HE2    0.000000  -0.026000     1.260000 ! 
		! for ethene, yin/adm jr., 12/95
!HB     0.000000  -0.022000     1.320000 ! ALLOW PEP ALI POL SUL ARO PRO ALC
                ! methane/ethane a.i. and ethane pure solvent, adm jr, 2/3/92
HC     0.000000  -0.046000     0.224500 ! ALLOW POL
                ! new, small polar Hydrogen, see also adm jr. JG 8/27/89
HP     0.000000  -0.030000     1.358200   0.000000  -0.030000     1.358200 ! ALLOW ARO
                ! JES 8/25/89 values from Jorgensen fit to hydration energy
HR1    0.000000  -0.046000     0.900000 ! ALLOW ARO
                ! adm jr., 6/27/90, his
HR2    0.000000  -0.046000     0.700000 ! ALLOW ARO
                ! adm jr., 6/27/90, his
HR3    0.000000  -0.007800     1.468000 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HS     0.000000  -0.100000     0.450000 ! ALLOW SUL
                ! methanethiol pure solvent, adm jr., 6/22/92
!new alkanes atoms types for conversion to new LJ parameters for c27 (see toppar_all22_prot_aliphatic_c27.str)
HA1     0.0       -0.045     1.3400 ! alkane, viv and adm jr., 4/07
HA2     0.0       -0.034     1.3400 ! alkane, viv and adm jr., 4/07
HA3     0.0       -0.024     1.3400 ! alkane, yin and mackerell, 4/98
!nitrogens
N      0.000000  -0.200000     1.850000   0.000000  -0.000100     1.850000 ! ALLOW   PRO
                ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93
NC2    0.000000  -0.200000     1.850000 ! ALLOW   POL
                ! JG 8/27/89; note: NH1 in ARG was changed to NC2.
NH1    0.000000  -0.200000     1.850000   0.000000  -0.200000     1.550000 ! ALLOW   PEP POL ARO
                ! This 1,4 vdW allows the C5 dipeptide minimum to exist.(LK)
NH2    0.000000  -0.200000     1.850000 ! ALLOW   POL
                ! adm jr.
NH3    0.000000  -0.200000     1.850000 ! ALLOW   POL
                ! adm jr.
NP     0.000000  -0.200000     1.850000 ! ALLOW  PRO
                ! N-terminal proline; from 6-31g* +ProNH2  RLD 9/28/90
NR1    0.000000  -0.200000     1.850000 ! ALLOW ARO
                ! His, adm jr., 9/4/89
NR2    0.000000  -0.200000     1.850000 ! ALLOW ARO
                ! His, adm jr., 9/4/89
NR3    0.000000  -0.200000     1.850000 ! ALLOW ARO
                ! His, adm jr., 9/4/89
NY     0.000000  -0.200000     1.850000 ! atm, indole vaporization 5/05
! oxygens
O      0.000000  -0.120000     1.700000   0.000000  -0.120000     1.400000 ! ALLOW   PEP POL
                ! This 1,4 vdW allows the C5 dipeptide minimum to exist.(LK)
OB     0.000000  -0.120000     1.700000   0.000000  -0.120000     1.400000 ! ALLOW   PEP POL ARO
                ! adm jr., 10/17/90, acetic acid carbonyl O
OC     0.000000  -0.120000     1.700000 ! ALLOW   POL ION
                ! JG 8/27/89
OH1    0.000000  -0.152100     1.770000 ! ALLOW   ALC ARO
                ! adm jr. 8/14/90, MeOH nonbond and solvent (same as TIP3P)
OS     0.000000  -0.152100     1.770000 ! ALLOW   ALC ARO
                ! adm jr. 9/17/90, avoid O* wildcard
! sulfurs
S      0.000000  -0.450000     2.000000 ! ALLOW   SUL ION
                ! adm jr., 3/3/92, methanethiol/ethylmethylsulfide pure solvent
SM     0.000000  -0.380000     1.975000 ! ALLOW  SUL  ION
                ! adm jr., 3/3/92, dimethyldisulphide pure solvent
SS     0.000000  -0.470000     2.200000 ! ALLOW  SUL
                ! methylthiolate to water and F.E. of solvation, adm jr. 6/1/92
NBFIX
!              Emin         Rmin
!            (kcal/mol)     (A)
NC2    OC       -0.154919   3.637 !  From osmotic pressure calibration
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