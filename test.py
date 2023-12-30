def parse_pdb_data(pdb_data):
    # Split the data into lines
    lines = pdb_data.split('\n')

    base_structure = []

    for line in lines:
        if line.startswith("ATOM"):
            parts = line.split()
            atom_type = parts[2]
            amino_acid = parts[3]
            chain = parts[4]
            residue_number = int(parts[5])
            x = float(parts[6])
            y = float(parts[7])
            z = float(parts[8])
            base_structure.append([atom_type, amino_acid, chain, residue_number, x, y, z])

    return base_structure

# Example usage with a PDB data string
pdb_data = """
HEADER    IMMUNE SYSTEM                           21-AUG-02   1MI0              
TITLE     CRYSTAL STRUCTURE OF THE REDESIGNED PROTEIN G VARIANT NUG2            
COMPND    MOL_ID: 1;                                                                    
SOURCE    MOL_ID: 1;                                                                                          
KEYWDS    ALPHA-BETA PROTEIN, REDESIGNED BETA-HAIRPIN, IMMUNE SYSTEM            
EXPDTA    X-RAY DIFFRACTION                                                     
AUTHOR    S.NAULI,B.KUHLMAN,I.LE TRONG,R.E.STENKAMP,D.C.TELLER,D.BAKER          
REVDAT   5   27-OCT-21 1MI0    1       SEQADV                                                                       
JRNL        AUTH   S.NAULI,B.KUHLMAN,I.LE TRONG,R.E.STENKAMP,D.C.TELLER,D.BAKER                                          
REMARK   2                                                                                                                             
DBREF  1MI0 A    5    61  PIR    A45063   A45063         328    384                        
SEQADV 1MI0 MET A   -3  PIR  A45063              EXPRESSION TAG                           
SEQRES   1 A   65  MET HIS HIS HIS HIS HIS HIS ALA MET ASP THR TYR LYS                   
FORMUL   3  HOH   *110(H2 O)                                                    
HELIX    1   1 ASP A   27  ASP A   41  1                                  15    
HELIX    2   2 ASP B   28  ALA B   40  1                                  13    
SHEET    1   A 4 THR A  18  ALA A  25  0                                        
SHEET    2   A 4 ASP A   6  VAL A  13 -1  N  TYR A   8   O  THR A  23           
SHEET    3   A 4 THR A  56  THR A  60  1  O  PHE A  57   N  VAL A  11           
SHEET    4   A 4 GLU A  47  ALA A  51 -1  N  THR A  49   O  THR A  58           
SHEET    1   B 4 THR B  19  ALA B  26  0                                        
SHEET    2   B 4 ASP B   7  LEU B  15 -1  N  TYR B   9   O  THR B  24           
SHEET    3   B 4 THR B  57  GLU B  62  1  O  VAL B  60   N  VAL B  12           
SHEET    4   B 4 GLU B  48  ALA B  52 -1  N  THR B  50   O  THR B  59           
CRYST1   47.330   73.790   39.180  90.00  96.00  90.00 C 1 2 1       8          
ORIGX1      1.000000  0.000000  0.000000        0.00000                         
ORIGX2      0.000000  1.000000  0.000000        0.00000                         
ORIGX3      0.000000  0.000000  1.000000        0.00000                         
SCALE1      0.021128  0.000000  0.002221        0.00000                         
SCALE2      0.000000  0.013552  0.000000        0.00000                         
SCALE3      0.000000  0.000000  0.025664        0.00000                         
ATOM      1  N   HIS A   1       8.499  -7.023  -3.949  1.00 73.05           N  
ATOM      2  CA  HIS A   1       9.115  -7.134  -2.593  1.00 72.10           C  
ATOM      3  C   HIS A   1       8.115  -6.777  -1.494  1.00 72.25           C  
ATOM      4  O   HIS A   1       7.801  -5.622  -1.201  1.00 72.26           O  
ATOM      5  CB  HIS A   1      10.370  -6.278  -2.496  1.00 72.09           C  
ATOM      6  N   HIS A   2       7.589  -7.812  -0.854  1.00 71.45           N  
ATOM      7  CA  HIS A   2       6.603  -7.820   0.192  1.00 69.97           C  
ATOM      8  C   HIS A   2       5.861  -6.490   0.392  1.00 69.63           C  
ATOM      9  O   HIS A   2       6.386  -5.564   1.000  1.00 69.03           O  
ATOM     10  CB  HIS A   2       7.163  -8.246   1.543  1.00 68.48           C  
ATOM     11  N   HIS A   3       4.599  -6.523  -0.023  1.00 68.52           N  
ATOM     12  CA  HIS A   3       3.681  -5.404   0.057  1.00 66.70           C  
ATOM     13  C   HIS A   3       2.638  -5.560   1.171  1.00 65.32           C  
ATOM     14  O   HIS A   3       1.525  -6.035   0.988  1.00 65.82           O  
ATOM     15  CB  HIS A   3       2.998  -5.221  -1.294  1.00 67.54           C  
ATOM     16  N   ALA A   4       3.034  -5.157   2.375  1.00 61.68           N  
ATOM     17  CA  ALA A   4       2.205  -5.125   3.560  1.00 60.72           C  
ATOM     18  C   ALA A   4       2.491  -3.821   4.322  1.00 60.32           C  
ATOM     19  O   ALA A   4       2.541  -3.776   5.553  1.00 61.16           O  
ATOM     20  CB  ALA A   4       2.608  -6.346   4.376  1.00 60.59           C  
ATOM     21  N   MET A   5       2.773  -2.731   3.597  1.00 58.55           N  
ATOM     22  CA  MET A   5       3.276  -1.524   4.217  1.00 55.89           C  
ATOM     23  C   MET A   5       2.414  -0.300   4.103  1.00 53.05           C  
ATOM     24  O   MET A   5       1.565  -0.100   3.242  1.00 55.76           O  
ATOM     25  CB  MET A   5       4.577  -1.057   3.543  1.00 56.55           C  
ATOM     26  CG  MET A   5       5.543  -2.060   3.019  1.00 56.46           C  
ATOM     27  SD  MET A   5       7.155  -2.041   3.783  1.00 58.18           S  
ATOM     28  CE  MET A   5       8.201  -2.238   2.315  1.00 56.50           C  
ATOM     29  N   ASP A   6       2.707   0.615   5.016  1.00 49.98           N  
ATOM     30  CA  ASP A   6       1.999   1.888   4.915  1.00 46.76           C  
ATOM     31  C   ASP A   6       2.970   2.888   4.290  1.00 45.15           C  
ATOM     32  O   ASP A   6       4.180   2.662   4.196  1.00 43.38           O  
ATOM     33  CB  ASP A   6       1.435   2.250   6.259  1.00 49.50           C  
ATOM     34  CG  ASP A   6       0.110   1.587   6.578  1.00 52.10           C  
ATOM     35  OD1 ASP A   6       0.081   0.424   7.004  1.00 53.10           O  
ATOM     36  OD2 ASP A   6      -0.926   2.233   6.420  1.00 54.31           O  
ATOM     37  N   THR A   7       2.423   4.015   3.837  1.00 40.51           N  
ATOM     38  CA  THR A   7       3.306   5.044   3.314  1.00 39.11           C  
ATOM     39  C   THR A   7       3.425   6.052   4.463  1.00 35.34           C  
ATOM     40  O   THR A   7       2.457   6.414   5.111  1.00 32.63           O  
ATOM     41  CB  THR A   7       2.837   5.671   2.007  1.00 42.75           C  
ATOM     42  OG1 THR A   7       3.090   4.736   0.929  1.00 43.01           O  
ATOM     43  CG2 THR A   7       3.624   6.969   1.748  1.00 40.63           C  
ATOM     44  N   TYR A   8       4.670   6.348   4.873  1.00 31.62           N  
ATOM     45  CA  TYR A   8       4.911   7.366   5.881  1.00 28.45           C  
ATOM     46  C   TYR A   8       5.567   8.638   5.297  1.00 28.09           C  
ATOM     47  O   TYR A   8       6.196   8.602   4.233  1.00 28.30           O  
ATOM     48  CB  TYR A   8       5.891   6.831   6.924  1.00 27.63           C  
ATOM     49  CG  TYR A   8       5.239   5.707   7.717  1.00 28.58           C  
ATOM     50  CD1 TYR A   8       4.594   5.911   8.897  1.00 31.31           C  
ATOM     51  CD2 TYR A   8       5.246   4.412   7.184  1.00 30.90           C  
ATOM     52  CE1 TYR A   8       3.997   4.865   9.591  1.00 32.93           C  
ATOM     53  CE2 TYR A   8       4.642   3.352   7.873  1.00 31.79           C  
ATOM     54  CZ  TYR A   8       4.051   3.589   9.061  1.00 33.37           C  
ATOM     55  OH  TYR A   8       3.434   2.552   9.766  1.00 37.13           O  
ATOM     56  N   LYS A   9       5.420   9.757   6.021  1.00 25.46           N  
ATOM     57  CA  LYS A   9       5.800  11.013   5.454  1.00 27.16           C  
ATOM     58  C   LYS A   9       6.459  11.902   6.470  1.00 24.04           C  
ATOM     59  O   LYS A   9       6.121  11.908   7.654  1.00 27.46           O  
ATOM     60  CB  LYS A   9       4.521  11.714   4.920  1.00 29.42           C  
ATOM     61  CG  LYS A   9       4.948  13.031   4.197  1.00 33.09           C  
ATOM     62  CD  LYS A   9       3.730  13.935   3.989  1.00 35.89           C  
ATOM     63  CE  LYS A   9       4.041  14.995   2.959  1.00 35.31           C  
ATOM     64  NZ  LYS A   9       2.986  15.981   2.597  1.00 40.38           N  
ATOM     65  N   LEU A  10       7.555  12.513   5.962  1.00 25.53           N  
ATOM     66  CA  LEU A  10       8.279  13.412   6.883  1.00 26.93           C  
ATOM     67  C   LEU A  10       8.085  14.863   6.408  1.00 24.40           C  
ATOM     68  O   LEU A  10       8.116  15.135   5.199  1.00 28.33           O  
ATOM     69  CB  LEU A  10       9.789  13.073   6.887  1.00 26.90           C  
ATOM     70  CG  LEU A  10      10.709  14.197   7.437  1.00 27.05           C  
ATOM     71  CD1 LEU A  10      10.457  14.256   8.951  1.00 28.22           C  
ATOM     72  CD2 LEU A  10      12.149  13.759   7.049  1.00 26.75           C  
ATOM     73  N   VAL A  11       7.712  15.777   7.313  1.00 24.67           N  
ATOM     74  CA  VAL A  11       7.502  17.177   6.762  1.00 22.53           C  
ATOM     75  C   VAL A  11       8.456  18.058   7.567  1.00 23.40           C  
ATOM     76  O   VAL A  11       8.445  18.053   8.799  1.00 28.36           O  
ATOM     77  CB  VAL A  11       6.011  17.604   6.930  1.00 23.18           C  
ATOM     78  CG1 VAL A  11       5.827  19.080   6.644  1.00 21.48           C  
ATOM     79  CG2 VAL A  11       5.156  16.728   5.983  1.00 24.68           C  
ATOM     80  N   ILE A  12       9.266  18.846   6.939  1.00 27.52           N  
ATOM     81  CA  ILE A  12      10.273  19.725   7.537  1.00 26.71           C  
ATOM     82  C   ILE A  12       9.968  21.140   7.114  1.00 24.86           C  
ATOM     83  O   ILE A  12       9.896  21.499   5.952  1.00 27.04           O  
ATOM     84  CB  ILE A  12      11.693  19.328   7.015  1.00 30.54           C  
ATOM     85  CG1 ILE A  12      12.001  17.893   7.485  1.00 30.02           C  
ATOM     86  CG2 ILE A  12      12.656  20.292   7.714  1.00 31.50           C  
ATOM     87  CD1 ILE A  12      13.315  17.277   6.949  1.00 35.28           C  
ATOM     88  N   VAL A  13       9.554  21.947   8.109  1.00 25.55           N  
ATOM     89  CA  VAL A  13       9.109  23.300   7.827  1.00 27.76           C  
ATOM     90  C   VAL A  13      10.288  24.260   8.255  1.00 29.13           C  
ATOM     91  O   VAL A  13      10.797  24.119   9.354  1.00 32.52           O  
ATOM     92  CB  VAL A  13       7.920  23.572   8.753  1.00 29.63           C  
ATOM     93  CG1 VAL A  13       7.475  25.027   8.632  1.00 32.24           C  
ATOM     94  CG2 VAL A  13       6.735  22.686   8.320  1.00 29.51           C  
TER     105      GLU A  62                                                            
END                                                                             
"""

base_structure = parse_pdb_data(pdb_data)
print(base_structure)  # This will print the extracted base structure
