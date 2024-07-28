import torch
import torch.nn as nn
import pandas as pd

class AtomCoordinates(nn.Module):
    def __init__(self):
        super(AtomCoordinates, self).__init__()
        self.init()
    
    def init(self):
        self.coordinates = None
        self.velocities = None
        self.pdb_table = None

        self.res_num = []
        self.amino_acids = []
        self.atom_names = []
        self.atom_types = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        
        self.protein_sequence = ""
        self.water_num = 100

    def forward(self):
        return self.coordinates
    
    def initialize(self, protein_sequence, device, pdb_data=None):
        self.init()

        # Read amino acid data
        self.aa_coordinates = {}
        current_aa = ''
        with open("parameter/aa_coordinates", 'r') as file:
            for line in file:
                if line.startswith('>'):
                    current_aa = line[1:].strip()
                    self.aa_coordinates[current_aa] = ''
                else:
                    self.aa_coordinates[current_aa] += line
        self.aa_connectivity = {}
        current_aa = ''
        with open("parameter/aa_connectivity", 'r') as file:
            for line in file:
                if line.startswith('>'):
                    current_aa = line[1:].strip()
                    self.aa_connectivity[current_aa] = ''
                else:
                    self.aa_connectivity[current_aa] += line      

        # Build amino acid coordinates
        if pdb_data is not None:
            aa_name_conv = {}
            for aa in self.aa_connectivity:
                parsed_aa = self.parse_pdb(self.aa_coordinates[aa])
                aa_name_conv[parsed_aa[0]['Residue']] = aa
            parsed_pdb = self.parse_pdb(pdb_data)
            pdb_sequence = []
            res_num = -1
            for line in parsed_pdb:
                if line['ResidueNum'] != res_num:
                    protein_sequence += aa_name_conv[line['Residue']]
                    pdb_sequence.append([])
                    res_num = line['ResidueNum']
                pdb_sequence[-1].append(line)
        parsed_data = []
        atom_num, res_num, current_pos, pos_list = 0, 0, (0., 0., 0.), []
        for id, aa in enumerate(protein_sequence):
            if aa not in self.aa_coordinates:
                continue
            parsed_aa = self.parse_pdb(self.aa_coordinates[aa])
            if id == 0:
                for i, d in enumerate(parsed_aa):
                    if aa == "P" and d.get("AtomType") == "HA":
                        new_H = [d.copy() for i in range(2)]
                        new_a = ["H1", "H2"]
                        new_x = [ 0.000, -0.461]
                        new_y = [-0.980,  0.327]
                        new_z = [ 0.000, -0.800]
                        for j in range(2):
                            new_H[j]['AtomType'] = new_a[j]
                            new_H[j]['X'] = new_x[j]
                            new_H[j]['Y'] = new_y[j]
                            new_H[j]['Z'] = new_z[j]
                        parsed_aa[i:i] = new_H
                        break
                    if d.get("AtomType") == "H":
                        new_H = [d.copy() for i in range(3)]
                        new_a = ["H1", "H2", "H3"]
                        new_x = [ 0.000, -0.461, -0.461]
                        new_y = [-0.980,  0.327,  0.327]
                        new_z = [ 0.000,  0.800, -0.800]
                        for j in range(3):
                            new_H[j]['AtomType'] = new_a[j]
                            new_H[j]['X'] = new_x[j]
                            new_H[j]['Y'] = new_y[j]
                            new_H[j]['Z'] = new_z[j]
                        parsed_aa[i:i+1] = new_H
                        break
            res_num += 1
            if pdb_data is not None:
                atom_first = pdb_sequence[id][0]
                current_pos = (atom_first['X'], atom_first['Y'], atom_first['Z'])
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
            if pdb_data is not None:
                aa_num_conv = {}
                for lid, line in enumerate(parsed_aa):
                    aa_num_conv[line['AtomType']] = lid
                parsed_pdb = pdb_sequence[id]
                for line in parsed_pdb:
                    atomtype = line['AtomType']
                    if atomtype in aa_num_conv:
                        parsed_aa[aa_num_conv[atomtype]]['X'] = line['X']
                        parsed_aa[aa_num_conv[atomtype]]['Y'] = line['Y']
                        parsed_aa[aa_num_conv[atomtype]]['Z'] = line['Z']
            for line in parsed_aa:
                pos_list.append([line["X"], line["Y"], line["Z"]])
            parsed_data += parsed_aa
        current_pos = torch.mean(torch.tensor(pos_list), dim=0)
        for i in range(self.water_num):
            parsed_aa = self.parse_pdb(self.aa_coordinates['H2O'])
            res_num += 1
            for line in parsed_aa:
                atom_num += 1
                line["AtomNum"] = atom_num
                line["ResidueNum"] = res_num
                pos_tmp = current_pos + 20 * (2 * torch.rand(3) - 1)
                line["X"] += pos_tmp[0].item()
                line["Y"] += pos_tmp[1].item()
                line["Z"] += pos_tmp[2].item()
            parsed_data += parsed_aa
        self.pdb_table = pd.DataFrame(parsed_data)
        self.coordinates = nn.Parameter(torch.tensor(self.pdb_table[['X', 'Y', 'Z']].values).to(device))
        self.velocities = torch.zeros_like(self.coordinates, device=device)
        self.protein_sequence = protein_sequence

        # Build amino acid connectivity        
        nam_idx = {"-C":-1}
        graph = []
        atom_idx = 0
        for id, aa in enumerate(protein_sequence):
            if aa not in self.aa_connectivity:
                continue
            lines = self.aa_connectivity[aa].split('\n')
            if id == 0:
                if aa == "P":
                    pos = lines.index("ATOM  N    N")
                    lines[pos] = "ATOM  N    NP"
                    pos = lines.index("ATOM  HA   HB1")
                    lines[pos:pos] = ["ATOM  H1   HC", "ATOM  H2   HC"]
                    pos = lines.index("BOND  HA   CA")
                    lines[pos:pos] = ["BOND  H1   N", "BOND  H2   N"]
                else:
                    pos = lines.index("ATOM  N    NH1")
                    lines[pos] = "ATOM  N    NH3"
                    pos = lines.index("ATOM  H    H")
                    lines[pos:pos+1] = ["ATOM  H1   HC", "ATOM  H2   HC", "ATOM  H3   HC"]
                    pos = lines.index("BOND  H    N")
                    lines[pos:pos+1] = ["BOND  H1   N", "BOND  H2   N", "BOND  H3   N"]
            if id == len(protein_sequence)-1:
                pos = lines.index("ATOM  C    C")
                lines[pos] = "ATOM  C    CC"
                pos = lines.index("ATOM  O    O")
                lines[pos] = "ATOM  O    OC"
                pos = lines.index("BOND  N    -C")
                lines.insert(pos, "ATOM  OXT  OC")
                lines.append("BOND  OXT  C")
            for line in lines:
                record = line.split()
                if line.startswith("ATOM"):
                    self.res_num.append(id)
                    self.amino_acids.append(aa)
                    self.atom_names.append(record[1])
                    self.atom_types.append(record[2])
                    nam_idx[record[1]] = atom_idx
                    graph.append([])
                    atom_idx += 1
                if line.startswith("BOND"):
                    v, u = nam_idx[record[1]], nam_idx[record[2]]
                    if u != -1:
                        graph[v].append(u)
                        graph[u].append(v)
                    self.dfs(graph, v, -1, [v])
            nam_idx["-C"] = nam_idx["C"]
        for i in range(len(graph)):
            if len(graph[i]) == 3:
                self.impropers.append([i]+graph[i])
        for id in range(self.water_num):
            aa = "H2O"
            lines = self.aa_connectivity[aa].split('\n')
            for line in lines:
                record = line.split()
                if line.startswith("ATOM"):
                    self.res_num.append(id)
                    self.amino_acids.append(aa)
                    self.atom_names.append(record[1])
                    self.atom_types.append(record[2])
                    nam_idx[record[1]] = atom_idx
                    graph.append([])
                    atom_idx += 1
                if line.startswith("BOND"):
                    v, u = nam_idx[record[1]], nam_idx[record[2]]
                    if u != -1:
                        graph[v].append(u)
                        graph[u].append(v)
                    self.dfs(graph, v, -1, [v])

    def update_positions(self, dt):
        with torch.no_grad():
            self.coordinates.add_(self.velocities * dt)

    def update_velocities(self, acceleration, dt, vol=0.005, damp=0.99):
        with torch.no_grad():
            self.velocities.add_(acceleration * dt)
            self.velocities += torch.randn_like(self.velocities) * vol
            self.velocities.mul_(damp)
            self.velocities = torch.tanh(self.velocities)

    def addNoise(self, vol=1, damp=0.99):
        self.velocities += torch.randn_like(self.coordinates)
        self.velocities *= damp
        self.coordinates.data.add_(self.velocities * vol)
    
    def render(self):
        self.pdb_table['X'] = self.coordinates[:, 0].cpu().detach().numpy()
        self.pdb_table['Y'] = self.coordinates[:, 1].cpu().detach().numpy()
        self.pdb_table['Z'] = self.coordinates[:, 2].cpu().detach().numpy()
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