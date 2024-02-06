import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np 
import pandas as pd 

from flask import Flask, jsonify, request
from flask_cors import CORS
import random

app = Flask(__name__)
CORS(app)

class ForceField(nn.Module):
    def __init__(self):
        super(ForceField, self).__init__()
        self.bond_prm = None
        self.angle_prm = None
        self.dihedral_prm = None
        self.vanderwaals_prm = None
        self.electro_prm = None

        self.bond_k_b = None
        self.angle_k_t = None
        self.dihedral_k_c_d = None
        self.vanderwaals_e_r = None
        self.electro_kqq = None

    def forward(self, x):
        return x

class AtomCoordinates(nn.Module):
    def __init__(self, ):
        super(AtomCoordinates, self).__init__()
        self.coordinates = None

        self.amino_acids = None
        self.atom_names = None
        self.atom_types = None
        self.bonds = None
        self.angles = None
        self.dihedrals = None
    
    def initialize(self, coordinates):
        self.coordinates = nn.Parameter(coordinates)

    def forward(self):
        return self.coordinates

class MainApp:
    def __init__(self):
        super().__init__()
        self.amino_acid_data = self.load_amino_acid_data()
        self.protein_table = pd.DataFrame(columns=["ATOM", "AtomNum", "AtomType", "Residue", "Chain", "ResidueNum", 
                                             "X", "Y", "Z", "Occupancy", "TempFactor", "Element"])
        self.atom_coordinates = AtomCoordinates()

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

    def load_amino_acid_data(self, filename="parameter/aa_coordinates"):
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
        parsed_data = []
        atom_num, res_num, current_pos = 0, 0, (0., 0., 0.)
        for aa in protein_sequence:
            if aa not in self.amino_acid_data:
                continue
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
        self.protein_table = pd.DataFrame(parsed_data)
        self.atom_coordinates.initialize(torch.tensor(self.protein_table[['X', 'Y', 'Z']].values))


mainapp = MainApp()

@app.route('/generateProteinStructure', methods=['POST'])
def generate_protein_structure():
    data = request.json
    protein_sequence = data['sequence']
    mainapp.protein_from_sequence(protein_sequence)
    return jsonify({"message": "Protein structure generated"})

@app.route('/getProteinStructure', methods=['GET'])
def get_protein_structure():
    pdb_formatted_data = ""
    for index, row in mainapp.protein_table.iterrows():
        pdb_line = f"ATOM  {int(row['AtomNum']):>5} {row['AtomType']:<4} {row['Residue']:<3} {row['Chain']:>1}{int(row['ResidueNum']):>4}    {float(row['X']):>8.3f}{float(row['Y']):>8.3f}{float(row['Z']):>8.3f}{float(row['Occupancy']):>6.2f}{float(row['TempFactor']):>6.2f}          {row['Element']:>2}\n"
        pdb_formatted_data += pdb_line
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
