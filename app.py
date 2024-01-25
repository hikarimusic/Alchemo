import torch
import numpy as np 
import pandas as pd 

from flask import Flask, jsonify, request
from flask_cors import CORS
import random

app = Flask(__name__)
CORS(app)


# pdb_data = """
# ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  
# ATOM      2  CA  ALA A   1       1.456   0.000   0.000  1.00  0.00           C  
#     ATOM      3  C   ALA A   1       1.930   0.000   1.463  1.00  0.00           C  
# ATOM      4  O   ALA A   1       1.160   0.000   2.421  1.00  0.00           O  
# ATOM      5  HN  ALA A   1      -0.495   0.091   0.883  1.00  0.00           H  
# ATOM      6  HA  ALA A   1       1.799  -0.926  -0.476  1.00  0.00           H  
# ATOM      7  CB  ALA A   1       2.010   1.208  -0.746  1.00  0.00           C  
# ATOM      8  1HB ALA A   1       1.657   1.230  -1.782  1.00  0.00           H  
# ATOM      9  2HB ALA A   1       1.695   2.143  -0.268  1.00  0.00           H  
# ATOM     10  3HB ALA A   1       3.104   1.197  -0.763  1.00  0.00           H  
# ATOM     11  HN  ALA A   1      -0.508   0.077  -0.877  1.00  0.00           H  

# ATOM     12  N   ALA A   2       3.241  -0.000   1.742  1.00  0.00           N  
# ATOM     13  CA  ALA A   2       3.690   0.000   3.127  1.00  0.00           C  
# ATOM     14  C   ALA A   2       5.228   0.000   3.127  1.00  0.00           C  
# ATOM     15  O   ALA A   2       5.902  -0.000   2.098  1.00  0.00           O  
# ATOM     16  HN  ALA A   2       3.928  -0.091   0.999  1.00  0.00           H  
# ATOM     17  HA  ALA A   2       3.343   0.926   3.600  1.00  0.00           H  
# ATOM     18  CB  ALA A   2       3.151  -1.208   3.884  1.00  0.00           C  
# ATOM     19  1HB ALA A   2       2.057  -1.230   3.869  1.00  0.00           H  
# ATOM     20  2HB ALA A   2       3.508  -2.143   3.438  1.00  0.00           H  
# ATOM     21  3HB ALA A   2       3.472  -1.197   4.931  1.00  0.00           H  
# ATOM     22  N   ALA A   3       5.898   0.000   4.287  1.00  0.00           N  
# ATOM     23  CA  ALA A   3       7.353   0.000   4.287  1.00  0.00           C  
# ATOM     24  C   ALA A   3       7.828   0.000   5.750  1.00  0.00           C  
# ATOM     25  O   ALA A   3       7.058   0.001   6.709  1.00  0.00           O  
# ATOM     26  HN  ALA A   3       5.403   0.092   5.170  1.00  0.00           H  
# ATOM     27  HA  ALA A   3       7.696  -0.926   3.811  1.00  0.00           H  
# ATOM     28  CB  ALA A   3       7.907   1.208   3.541  1.00  0.00           C  
# ATOM     29  1HB ALA A   3       7.555   1.230   2.505  1.00  0.00           H  
# ATOM     30  2HB ALA A   3       7.593   2.143   4.018  1.00  0.00           H  
# ATOM     31  3HB ALA A   3       9.002   1.197   3.524  1.00  0.00           H  
# ATOM     32  N   ALA A   4       9.138   0.000   6.029  1.00  0.00           N  
# ATOM     33  CA  ALA A   4       9.587   0.001   7.414  1.00  0.00           C  
# ATOM     34  C   ALA A   4      11.125   0.000   7.414  1.00  0.00           C  
# ATOM     35  O   ALA A   4      11.799  -0.000   6.386  1.00  0.00           O  
# ATOM     36  HN  ALA A   4       9.825  -0.092   5.286  1.00  0.00           H  
# ATOM     37  HA  ALA A   4       9.240   0.927   7.887  1.00  0.00           H  
# ATOM     38  CB  ALA A   4       9.048  -1.207   8.171  1.00  0.00           C  
# ATOM     39  1HB ALA A   4       7.954  -1.229   8.156  1.00  0.00           H  
# ATOM     40  2HB ALA A   4       9.405  -2.142   7.726  1.00  0.00           H  
# ATOM     41  3HB ALA A   4       9.370  -1.196   9.218  1.00  0.00           H  
# """ 

# pdb_data = """
# ATOM      1  N   ALA A   1       1.826   0.088   0.641  1.00  0.00           N  
# ATOM      2  CA  ALA A   1       0.429   0.558   0.495  1.00  0.00           C  
# ATOM      3  C   ALA A   1      -0.334   0.405   1.812  1.00  0.00           C  
# ATOM      4  O   ALA A   1      -0.106  -1.272  -1.087  1.00  0.00           O  
# """

AMINO_ACID = {
    "A": """
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  
ATOM      2  CA  ALA A   1       1.312   0.640   0.000  1.00  0.00           C  
ATOM      3  C   ALA A   1       2.416  -0.390   0.000  1.00  0.00           C  
ATOM      4  O   ALA A   1       2.416  -1.630   0.000  1.00  0.00           O 
ATOM      5  N+  ALA A   1       3.658   0.087   0.000  1.00  0.00           N  
"""
}


class MainApp:
    def __init__(self):
        super().__init__()
        self.protein = pd.DataFrame(columns=["ATOM", "AtomNum", "AtomType", "Residue", "Chain", "ResidueNum", 
                                             "X", "Y", "Z", "Occupancy", "TempFactor", "Element"])
        self.amino_acid_data = self.load_amino_acid_data()

    def parse_pdb_data(self, pdb_data):
        lines = pdb_data.split('\n')
        base_structure = []

        for line in lines:
            if line.startswith("ATOM"):
                record_type = line[0:6].strip()
                atom_serial_number = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                alt_loc = line[16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21].strip()
                residue_seq_number = int(line[22:26].strip())
                insertion_code = line[26].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occupancy = float(line[54:60].strip())
                temp_factor = float(line[60:66].strip())
                segment_id = line[72:76].strip()
                element_symbol = line[76:78].strip()
                charge = line[78:80].strip()
                
                base_structure.append({
                    "ATOM": record_type, 
                    "AtomNum": atom_serial_number, 
                    "AtomType": atom_name, 
                    "Residue": residue_name, 
                    "Chain": chain_id, 
                    "ResidueNum": residue_seq_number, 
                    "X": x, "Y": y, "Z": z, 
                    "Occupancy": occupancy, 
                    "TempFactor": temp_factor, 
                    "Element": element_symbol
                })

        return base_structure

    def load_amino_acid_data(self, filename="amino_acids.pdb"):
        amino_acid_data = {}
        current_aa = ''
        with open(filename, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    current_aa = line[1].strip()
                    amino_acid_data[current_aa] = ''
                else:
                    amino_acid_data[current_aa] += line
        return amino_acid_data

    def protein_from_sequence(self, protein_sequence):
        # parsed_data = self.parse_pdb_data()
        parsed_data = []
        atom_num, res_num, current_pos = 0, 0, (0., 0., 0.)
        for aa in protein_sequence:
            parsed_aa = self.parse_pdb_data(self.amino_acid_data[aa])
            for line in parsed_aa:
                line["AtomNum"] += atom_num
                atom_num += 1
                line["ResidueNum"] += res_num
                if res_num % 2 == 1:
                    line["Y"] *= -1
                    line["Z"] *= -1
                line["X"] += current_pos[0]
                line["Y"] += current_pos[1]
                line["Z"] += current_pos[2]
            res_num += 1
            current_pos = (parsed_aa[-1]["X"], parsed_aa[-1]["Y"], parsed_aa[-1]["Z"])
            parsed_aa.pop()
            parsed_data += parsed_aa
        self.protein = pd.DataFrame(parsed_data)


mainapp = MainApp()

@app.route('/generateProteinStructure', methods=['POST'])
def generate_protein_structure():
    data = request.json
    protein_sequence = data['sequence']
    mainapp.protein_from_sequence(protein_sequence)
    return jsonify({"message": "Protein structure generated"})

@app.route('/getProteinStructure', methods=['GET'])
def get_protein_structure():
    # Initialize an empty string to store PDB formatted data
    pdb_formatted_data = ""

    # Iterate over each row in the DataFrame and format it into a PDB line
    for index, row in mainapp.protein.iterrows():
        pdb_line = f"ATOM  {int(row['AtomNum']):>5} {row['AtomType']:<4} {row['Residue']:<3} {row['Chain']:>1}{int(row['ResidueNum']):>4}    {float(row['X']):>8.3f}{float(row['Y']):>8.3f}{float(row['Z']):>8.3f}{float(row['Occupancy']):>6.2f}{float(row['TempFactor']):>6.2f}          {row['Element']:>2}\n"
        pdb_formatted_data += pdb_line

    # Return the formatted PDB data
    return jsonify(pdbData=pdb_formatted_data)












# def parse_pdb_data(pdb_data):
#     lines = pdb_data.split('\n')
#     base_structure = []

#     for line in lines:
#         if line.startswith("ATOM"):
#             record_type = line[0:6].strip()
#             atom_serial_number = int(line[6:11].strip())
#             atom_name = line[12:16].strip()
#             alt_loc = line[16].strip()
#             residue_name = line[17:20].strip()
#             chain_id = line[21].strip()
#             residue_seq_number = int(line[22:26].strip())
#             insertion_code = line[26].strip()
#             x = float(line[30:38].strip())
#             y = float(line[38:46].strip())
#             z = float(line[46:54].strip())
#             occupancy = float(line[54:60].strip())
#             temp_factor = float(line[60:66].strip())
#             segment_id = line[72:76].strip()
#             element_symbol = line[76:78].strip()
#             charge = line[78:80].strip()
            
#             base_structure.append([record_type, atom_serial_number, atom_name, alt_loc, residue_name, chain_id, 
#                                    residue_seq_number, insertion_code, x, y, z, occupancy, temp_factor, 
#                                    segment_id, element_symbol, charge])

#     return base_structure

# def generate_frames(base_structure, num_frames, displacement=0.1):
#     frames = []
#     for frame in range(num_frames):
#         frame_data = "HEADER    TEST PDB\n"
#         for atom in base_structure:
#             # Displacing only the x, y, z coordinates
#             displaced_x = atom[8] + random.uniform(-displacement, displacement)
#             displaced_y = atom[9] + random.uniform(-displacement, displacement)
#             displaced_z = atom[10] + random.uniform(-displacement, displacement)

#             frame_data += "ATOM  {:>5} {:<4}{:1}{:>3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<4}{:>2}{:2}\n".format(
#                 atom[1], atom[2], atom[3], atom[4], atom[5], atom[6], atom[7], displaced_x, displaced_y, 
#                 displaced_z, atom[11], atom[12], atom[13], atom[14], atom[15]
#             )
#         frame_data += "TER\nEND\n"
#         frames.append(frame_data)
        
#     return frames

# @app.route('/generateProteinStructure', methods=['POST'])
# def generate_protein_structure():
#     mainapp.protein_from_sequence()

#     # data = request.json
#     # protein_sequence = data['sequence']
#     # print("HELLO", protein_sequence)
#     # # Code to generate PDB data from protein sequence
#     # pdb_data = "generated PDB data for protein"
#     # return jsonify(pdbData=pdb_data)

# @app.route('/generateMoleculeStructure', methods=['POST'])
# def generate_molecule_structure():
#     data = request.json
#     molecule_sequence = data['sequence']
#     # Code to generate PDB data from molecule sequence
#     pdb_data = "generated PDB data for molecule"
#     return jsonify(pdbData=pdb_data)

# @app.route('/animateProtein', methods=['POST'])
# def animate_protein():
#     content = request.json
#     pdb_data = content['pdbData']
#     base_structure = parse_pdb_data(pdb_data)
#     frames = generate_frames(base_structure, 10)  # Generate 10 frames
#     return jsonify({'frames': frames})

if __name__ == '__main__':
    mainapp = MainApp()
    app.run(debug=True)
