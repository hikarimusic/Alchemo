import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd 

import logging
from flask import Flask, jsonify, request
from flask_cors import CORS

from atom import AtomCoordinates
from force import ForceField

app = Flask(__name__)
CORS(app)

log = logging.getLogger('werkzeug')
log.setLevel(logging.ERROR)

class MainApp:
    def __init__(self):
        super().__init__()
        self.atom_coordinates = AtomCoordinates()
        self.force_field = ForceField()
        self.device = 'cuda:0'
        self.step = 0
        self.potential_energy = 0
        self.potential_max = 0

    def protein_from_sequence(self, protein_sequence, pdb_data=None):
        self.atom_coordinates.initialize(protein_sequence, self.device, pdb_data)
        self.force_field.initialize(self.atom_coordinates, self.device)
        self.atom_coordinates.to(self.device)
        self.force_field.to(self.device)
        self.atom_optimizer = optim.Adam(self.atom_coordinates.parameters(), lr=0.3)
        # self.atom_optimizer = optim.NAdam(self.atom_coordinates.parameters(), lr=0.1)
        # self.atom_optimizer = optim.SGD(self.atom_coordinates.parameters(), lr=0.01, momentum=0.9)
        potential_energy = torch.sum(self.force_field(self.atom_coordinates()))
        potential_energy.backward()

    def protein_render(self):
        self.atom_coordinates.render()
        return self.atom_coordinates.render()
    
    def gradient_descent(self):
        self.atom_optimizer.zero_grad()
        self.potential_energy = torch.sum(self.force_field(self.atom_coordinates()))
        self.potential_energy.backward()
        torch.nn.utils.clip_grad_norm_(self.atom_coordinates.coordinates, 100)
        self.atom_optimizer.step()
        self.atom_coordinates.addNoise(vol=0.005, damp=0.99)
    
    def molecular_dynamics(self, dt=0.1):
        self.atom_optimizer.zero_grad()
        self.potential_energy = torch.sum(self.force_field(self.atom_coordinates()))
        self.potential_energy.backward()
        torch.nn.utils.clip_grad_norm_(self.atom_coordinates.coordinates, 100)
        acceleration = -self.atom_coordinates.coordinates.grad / self.force_field.mass_m['m'].unsqueeze(1)
        self.atom_coordinates.update_positions(0.5 * dt)
        self.atom_coordinates.update_velocities(acceleration, dt)
        self.atom_coordinates.update_positions(0.5 * dt)

    def protein_simulate(self):
        self.gradient_descent()
        self.step += 1
        if self.step % 100 == 0:
            print(self.potential_energy.item())

@app.route('/loadProtein', methods=['POST'])
def loadProtein():
    data = request.json
    protein_sequence = data['sequence']
    mainapp.protein_from_sequence(protein_sequence)
    return jsonify({"message": "Protein structure generated"})

@app.route('/loadPdbId', methods=['POST'])
def loadPdbId():
    data = request.json
    pdb_data = data['pdbData']
    mainapp.protein_from_sequence("", pdb_data)
    return jsonify({"message": "Protein structure generated"})

@app.route('/showProtein', methods=['GET'])
def showProtein():
    protein_table = mainapp.protein_render()
    pdb_formatted_data = ""
    for id, row in protein_table.iterrows():
        pdb_line = f"ATOM  {int(row['AtomNum']):>5} {row['AtomType']:<4} {row['Residue']:<3} {row['Chain']:>1}{int(row['ResidueNum']):>4}    {float(row['X']):>8.3f}{float(row['Y']):>8.3f}{float(row['Z']):>8.3f}{float(row['Occupancy']):>6.2f}{float(row['TempFactor']):>6.2f}          {row['Element']:>2}\n"
        pdb_formatted_data += pdb_line
    return jsonify(pdbData=pdb_formatted_data)

@app.route('/simulateProtein', methods=['GET'])
def simulateProtein():
    speed = int(request.args.get('speed', 1))
    for i in range(speed):
        mainapp.protein_simulate()
    protein_table = mainapp.protein_render()
    pdb_formatted_data = ""
    for id, row in protein_table.iterrows():
        pdb_line = f"ATOM  {int(row['AtomNum']):>5} {row['AtomType']:<4} {row['Residue']:<3} {row['Chain']:>1}{int(row['ResidueNum']):>4}    {float(row['X']):>8.3f}{float(row['Y']):>8.3f}{float(row['Z']):>8.3f}{float(row['Occupancy']):>6.2f}{float(row['TempFactor']):>6.2f}          {row['Element']:>2}\n"
        pdb_formatted_data += pdb_line
    return jsonify(pdbData=pdb_formatted_data)


if __name__ == '__main__':
    mainapp = MainApp()
    app.run(debug=True)




























# @app.route('/simulateProtein', methods=['POST'])
# def simulateProtein():
#     content = request.json
#     pdb_data = content['pdbData']
#     base_structure = parse_pdb_data(pdb_data)
#     frames = generate_frames(base_structure, 10)  # Generate 10 frames
#     return jsonify({'frames': frames})

# @app.route('/simulateProtein', methods=['GET'])
# def simulateProtein():
#     mainapp.protein_simulate()
#     protein_table = mainapp.protein_render()
#     pdb_formatted_data = ""
#     for id, row in protein_table.iterrows():
#         pdb_line = f"ATOM  {int(row['AtomNum']):>5} {row['AtomType']:<4} {row['Residue']:<3} {row['Chain']:>1}{int(row['ResidueNum']):>4}    {float(row['X']):>8.3f}{float(row['Y']):>8.3f}{float(row['Z']):>8.3f}{float(row['Occupancy']):>6.2f}{float(row['TempFactor']):>6.2f}          {row['Element']:>2}\n"
#         pdb_formatted_data += pdb_line
#     return jsonify(pdbData=pdb_formatted_data)









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

