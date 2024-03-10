data = '''
!
!V(improper) = Kpsi(psi - psi0)**2
!
!Kpsi: kcal/mole/rad**2
!psi0: degrees
!note that the second column of numbers (0) is ignored
!
!atom types           Kpsi                   psi0
!
HE2  HE2  CE2  CE2     3.0            0      0.00   ! 
		! for ethene, yin/adm jr., 12/95
HR1  NR1  NR2  CPH2    0.5000         0      0.0000 ! ALLOW ARO
                ! his, adm jr., 7/05/90
HR1  NR2  NR1  CPH2    0.5000         0      0.0000 ! ALLOW ARO
                ! his, adm jr., 7/05/90
HR3  CPH1 NR1  CPH1    0.5000         0      0.0000 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HR3  CPH1 NR2  CPH1    0.5000         0      0.0000 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HR3  CPH1 NR3  CPH1    1.0000         0      0.0000 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HR3  NR1  CPH1 CPH1    0.5000         0      0.0000 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
HR3  NR2  CPH1 CPH1    0.5000         0      0.0000 ! ALLOW ARO
                ! adm jr., 3/24/92, maintain old aliphatic H VDW params
N    C    CP1  CP3     0.0000         0      0.0000 ! ALLOW PRO
                ! 6-31g* AcProNH2 and ProNH2  RLD 5/19/92
NC2  X    X    C      45.0000         0      0.0000 ! ALLOW   PEP POL ARO
                ! mp2/6-311g** guan vibrational data, adm jr., 1/04
C   HC    HC   NC2      0.0          0      0.0
                ! mp2/6-311g** guan vibrational data, adm jr., 1/04
NC2  X    X    HC      -2.0          0      0.0
                ! mp2/6-311g** guan vibrational data, adm jr., 1/04
NH1  X    X    H      20.0000         0      0.0000 ! ALLOW   PEP POL ARO
                ! NMA Vibrational Modes (LK)
NH2  X    X    H       4.0000         0      0.0000 ! ALLOW   POL
                ! adm jr., 8/13/90 acetamide geometry and vibrations
NR1  CPH1 CPH2 H       0.4500         0      0.0000 ! ALLOW ARO
                ! his, adm jr., 7/05/90
NR1  CPH2 CPH1 H       0.4500         0      0.0000 ! ALLOW ARO
                ! his, adm jr., 7/05/90
NR3  CPH1 CPH2 H       1.2000         0      0.0000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
NR3  CPH2 CPH1 H       1.2000         0      0.0000 ! ALLOW ARO
                ! his, adm jr., 6/27/90
O    CP1  NH2  CC     45.0000         0      0.0000 ! ALLOW PEP POL PRO
                ! 6-31g* AcProNH2 and ProNH2  RLD 5/19/92
O    CT1  NH2  CC     45.0000         0      0.0000 ! ALLOW PEP POL
                ! adm jr., 8/13/90 acetamide geometry and vibrations
O    CT2  NH2  CC     45.0000         0      0.0000 ! ALLOW PEP POL
                ! adm jr., 8/13/90 acetamide geometry and vibrations
O    CT3  NH2  CC     45.0000         0      0.0000 ! ALLOW PEP POL
                ! adm jr., 8/13/90 acetamide geometry and vibrations
O    HA1  NH2  CC     45.0000         0      0.0000 ! ALLOW PEP POL PRO
                ! adm jr., 5/13/91, formamide geometry and vibrations
O    N    CT2  CC    120.0000         0      0.0000 ! ALLOW PEP POL PRO
                ! 6-31g* AcProNH2 and ProNH2  RLD 5/19/92
O    NH2  CP1  CC     45.0000         0      0.0000 ! ALLOW PEP POL PRO
                ! 6-31g* AcProNH2 and ProNH2  RLD 5/19/92
O    NH2  CT1  CC     45.0000         0      0.0000 ! ALLOW PEP POL
                ! adm jr., 8/13/90 acetamide geometry and vibrations
O    NH2  CT2  CC     45.0000         0      0.0000 ! ALLOW PEP POL
                ! adm jr., 8/13/90 acetamide geometry and vibrations
O    NH2  CT3  CC     45.0000         0      0.0000 ! ALLOW PEP POL
                ! adm jr., 8/13/90 acetamide geometry and vibrations
O    NH2  HA1  CC     45.0000         0      0.0000 ! ALLOW PEP POL
                ! adm jr., 5/13/91, formamide geometry and vibrations
O    X    X    C     120.0000         0      0.0000 ! ALLOW   PEP POL ARO
                ! NMA Vibrational Modes (LK)
OB   X    X    CD    100.0000         0      0.0000 ! ALLOW   ALC ARO POL
                ! adm jr., 10/17/90, acetic acid vibrations
OC   X    X    CC     96.0000         0      0.0000 ! ALLOW   PEP POL ARO ION
                ! 90.0->96.0 acetate, single impr (KK)
CC   X    X    CT1    96.0000         0      0.0000 ! ALLOW   PEP POL ARO ION
                ! 90.0->96.0 acetate, single impr (KK)
CC   X    X    CT2    96.0000         0      0.0000 ! ALLOW   PEP POL ARO ION
                ! 90.0->96.0 acetate, single impr (KK)
CC   X    X    CT3    96.0000         0      0.0000 ! ALLOW   PEP POL ARO ION
                ! 90.0->96.0 acetate, single impr (KK)
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