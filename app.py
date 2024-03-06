import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np 
import pandas as pd 

from flask import Flask, jsonify, request
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

class ForceField(nn.Module):
    def __init__(self):
        super(ForceField, self).__init__()
        self.init()
    
    def init(self):
        self.bond_prm = None
        self.angle_prm = None
        self.dihedral_prm = None
        self.vanderwaals_prm = None
        self.electro_prm = None

        self.bond_k_b = {}
        self.angle_k_t = {}
        self.dihedral_k_c_d = {}
        self.vanderwaals_e_r = {}
        self.electro_kqq = {}

    def forward(self, x):
        # bond
        pairs, kb, b0 = self.bond_k_b["bond"], self.bond_k_b["kb"], self.bond_k_b["b0"]
        vec_s = x[pairs[:, 0]]
        vec_t = x[pairs[:, 1]]
        dis = torch.sqrt(torch.sum((vec_s - vec_t) ** 2, dim=1))
        V_bond = torch.sum(kb * (dis - b0) ** 2)

        # vanderwaals
        eps_ij, rmin_ij = self.vanderwaals_e_r["eps_ij"], self.vanderwaals_e_r["rmin_ij"]
        dist = torch.cdist(x, x) + torch.eye(x.size()[0])
        rod = rmin_ij / dist
        vdw_potential = eps_ij * (torch.pow(rod, 12) - 2 * torch.pow(rod, 6))
        i_upper = torch.triu_indices(vdw_potential.size(0), vdw_potential.size(1), offset=1)
        V_vanderwaals = vdw_potential[i_upper[0], i_upper[1]].sum()

        V = V_bond + V_vanderwaals
        return V
    
    def initialize(self, atom_coordinates):
        self.init()

        # Read parameter data
        self.bond_prm = pd.read_csv("parameter/ff_bond", sep='\s+', header=None, names=["Atom1", "Atom2", "Kb", "b0"])
        self.bond_prm[['Atom1', 'Atom2']] = self.bond_prm.apply(
            lambda row: pd.Series(sorted([row['Atom1'], row['Atom2']])), axis=1
        )
        self.vanderwaals_prm = pd.read_csv("parameter/ff_vanderwaals", sep='\s+', header=None, names=["Atom", "epsilon", "Rmin/2"])

        # Build bond potential
        bonds, kb, b0 = atom_coordinates.bonds, [], []
        for x, y in atom_coordinates.bonds:
            A1 = atom_coordinates.atom_types[x]
            A2 = atom_coordinates.atom_types[y]
            A1, A2 = sorted([A1, A2])
            matching_row = self.bond_prm[(self.bond_prm['Atom1']==A1) & (self.bond_prm['Atom2']==A2)]
            kb.append(matching_row['Kb'].iloc[0])
            b0.append(matching_row['b0'].iloc[0])
        self.bond_k_b = {
            "bond": torch.tensor(bonds),
            "kb": torch.tensor(kb),
            "b0": torch.tensor(b0)
        }

        # Build vanderwaals potential
        eps, rmin = [], []
        for A in atom_coordinates.atom_types:
            matching_row = self.vanderwaals_prm[(self.vanderwaals_prm['Atom']==A)]
            eps.append(matching_row["epsilon"].iloc[0])
            rmin.append(matching_row["Rmin/2"].iloc[0])
        eps = torch.tensor(eps)
        rmin = torch.tensor(rmin)
        eps_ij = torch.sqrt(eps[:, None] * eps)
        rmin_ij = rmin[:, None] + rmin
        for bond in atom_coordinates.bonds:
            i, j = bond
            eps_ij[i, j] = 0
            eps_ij[j, i] = 0
        for angle in atom_coordinates.angles:
            i, j = angle[0], angle[2]
            eps_ij[i, j] = 0
            eps_ij[j, i] = 0
        self.vanderwaals_e_r = {
            "eps_ij": eps_ij,
            "rmin_ij": rmin_ij
        }

class AtomCoordinates(nn.Module):
    def __init__(self):
        super(AtomCoordinates, self).__init__()
        self.init()
    
    def forward(self):
        return self.coordinates
    
    def init(self):
        self.coordinates = None
        self.pdb_table = None

        self.amino_acids = []
        self.atom_names = []
        self.atom_types = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
    
    def initialize(self, protein_sequence):
        self.init()

        # Read amino acid data
        self.aa_coordinates = {}
        current_aa = ''
        with open("parameter/aa_coordinates", 'r') as file:
            for line in file:
                if line.startswith('>'):
                    current_aa = line[1].strip()
                    self.aa_coordinates[current_aa] = ''
                else:
                    self.aa_coordinates[current_aa] += line
        self.aa_connectivity = {}
        current_aa = ''
        with open("parameter/aa_connectivity", 'r') as file:
            for line in file:
                if line.startswith('>'):
                    current_aa = line[1].strip()
                    self.aa_connectivity[current_aa] = ''
                else:
                    self.aa_connectivity[current_aa] += line    

        # Build amino acid coordinates
        parsed_data = []
        atom_num, res_num, current_pos = 0, 0, (0., 0., 0.)
        for id, aa in enumerate(protein_sequence):
            if aa not in self.aa_coordinates:
                continue
            parsed_aa = self.parse_pdb(self.aa_coordinates[aa])
            if id == 0:
                for i, d in enumerate(parsed_aa):
                    if d.get("AtomType") == "H":
                        new_H = [d.copy() for i in range(3)]
                        new_a = ["H1", "H2", "H3"]
                        new_x = [ 0.000, -0.461, -0.461]
                        new_y = [-0.980,  0.327, -0.327]
                        new_z = [ 0.000,  0.800, -0.800]
                        for j in range(3):
                            new_H[j]['AtomType'] = new_a[j]
                            new_H[j]['X'] = new_x[j]
                            new_H[j]['Y'] = new_y[j]
                            new_H[j]['Z'] = new_z[j]
                        parsed_aa[i:i+1] = new_H
                        break
            res_num += 1
            for line in parsed_aa:
                atom_num += 1
                line["AtomNum"] = atom_num
                line["ResidueNum"] = res_num
                if res_num % 2 == 0:
                    line["Y"] *= -1
                    line["Z"] *= -1
                line["X"] += current_pos[0]
                line["Y"] += current_pos[1]
                line["Z"] += current_pos[2]
            if id == len(protein_sequence)-1:
                parsed_aa[-1]["AtomType"] = "OXT"
                parsed_aa[-1]["Residue"] = parsed_aa[-2]["Residue"]
                parsed_aa[-1]["Element"] = "O"
            else:
                current_pos = (parsed_aa[-1]["X"], parsed_aa[-1]["Y"], parsed_aa[-1]["Z"])
                parsed_aa.pop()
                atom_num -= 1
            parsed_data += parsed_aa
        self.pdb_table = pd.DataFrame(parsed_data)
        self.coordinates = nn.Parameter(torch.tensor(self.pdb_table[['X', 'Y', 'Z']].values))

        # Build amino acid connectivity        
        nam_idx = {"-C":-1}
        graph = []
        atom_idx = 0
        for id, aa in enumerate(protein_sequence):
            if aa not in self.aa_connectivity:
                continue
            lines = self.aa_connectivity[aa].split('\n')
            if id == 0:
                pos = lines.index("ATOM  N     NH1")
                lines[pos] = "ATOM  N     NH3"
                pos = lines.index("ATOM  H     H")
                lines[pos:pos+1] = ["ATOM  H1    HC", "ATOM  H2    HC", "ATOM  H3    HC"]
                pos = lines.index("BOND  H     N")
                lines[pos:pos+1] = ["BOND  H1    N", "BOND  H2    N", "BOND  H3    N"]
            if id == len(protein_sequence)-1:
                pos = lines.index("ATOM  C     C")
                lines[pos] = "ATOM  C     CC"
                pos = lines.index("ATOM  O     O")
                lines[pos] = "ATOM  O     OC"
                pos = lines.index("BOND  N     -C")
                lines.insert(pos, "ATOM  OXT   OC")
                lines.append("BOND  OXT   C")
            for line in lines:
                record = line.split()
                if line.startswith("ATOM"):
                    self.amino_acids.append(aa)
                    self.atom_names.append(record[1])
                    self.atom_types.append(record[2])
                    nam_idx[record[1]] = atom_idx
                    graph.append([])
                    atom_idx += 1
                    if id == len(protein_sequence)-1 and self.atom_types[-1] == "O":
                        self.atom_types[-1] = "OC"
                if line.startswith("BOND"):
                    v, u = nam_idx[record[1]], nam_idx[record[2]]
                    if u != -1:
                        graph[v].append(u)
                        graph[u].append(v)
                    self.dfs(graph, v, -1, [v])
            nam_idx["-C"] = nam_idx["C"]
    
    def addNoise(self, vol=1):
        noise = torch.randn_like(self.coordinates) * vol
        self.coordinates.data.add_(noise)
    
    def render(self):
        self.pdb_table['X'] = self.coordinates[:, 0].detach().numpy()
        self.pdb_table['Y'] = self.coordinates[:, 1].detach().numpy()
        self.pdb_table['Z'] = self.coordinates[:, 2].detach().numpy()
        return self.pdb_table

    def dfs(self, adj, v, p, arr):
        if len(arr) > 4:
            return
        if len(arr) == 2:
            self.bonds.append(arr[::-1])
        if len(arr) == 3:
            self.angles.append(arr[::-1])
        if len(arr) == 4:
            self.dihedrals.append(arr[::-1])
        for u in adj[v]:
            if u != p:
                self.dfs(adj, u, v, arr+[u])

    def parse_pdb(self, pdb_data):
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

class MainApp:
    def __init__(self):
        super().__init__()
        self.atom_coordinates = AtomCoordinates()
        self.force_field = ForceField()

    def protein_from_sequence(self, protein_sequence):
        self.atom_coordinates.initialize(protein_sequence)
        self.force_field.initialize(self.atom_coordinates)
        self.atom_optimizer = optim.Adam(self.atom_coordinates.parameters(), lr=0.1)
        potential_energy = torch.sum(self.force_field(self.atom_coordinates()))
        potential_energy.backward()
    
    def protein_render(self):
        self.atom_coordinates.render()
        return self.atom_coordinates.render()

    def protein_simulate(self):
        self.atom_optimizer.zero_grad()
        potential_energy = torch.sum(self.force_field(self.atom_coordinates()))
        potential_energy.backward()
        self.atom_optimizer.step()
        # self.atom_coordinates.addNoise(vol=0.02)

@app.route('/loadProtein', methods=['POST'])
def loadProtein():
    data = request.json
    protein_sequence = data['sequence']
    mainapp.protein_from_sequence(protein_sequence)
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

