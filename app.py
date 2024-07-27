import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np 
import pandas as pd 

import logging
from flask import Flask, jsonify, request
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

log = logging.getLogger('werkzeug')
log.setLevel(logging.ERROR)

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
        atom_num, res_num, current_pos = 0, 0, (0., 0., 0.)
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
            parsed_data += parsed_aa
        self.pdb_table = pd.DataFrame(parsed_data)
        self.coordinates = nn.Parameter(torch.tensor(self.pdb_table[['X', 'Y', 'Z']].values).to(device))
        self.velocities = torch.zeros_like(self.coordinates, device=device)

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
    
class ForceField(nn.Module):
    def __init__(self):
        super(ForceField, self).__init__()
        self.init()
    
    def init(self):
        self.mass_prm = None
        self.bond_prm = None
        self.angle_prm = None
        self.ureybradley_prm = None
        self.dihedral_prm = None
        self.improper_prm = None
        self.vanderwaals_prm = None
        self.electro_prm = None

        self.mass_m = {}
        self.bond_k_b = {}
        self.angle_k_t = {}
        self.ureybradley_k_s = {}
        self.dihedral_k_n_d = {}
        self.improper_k_p = {}
        self.vanderwaals_e_r = {}
        self.electro_kqq = {}

    def forward(self, x, central=False):
        # Bond
        bond, kb, b0 = self.bond_k_b["bond"], self.bond_k_b["kb"], self.bond_k_b["b0"]
        vec_s = x[bond[:, 0]]
        vec_t = x[bond[:, 1]]
        dis = torch.sqrt(torch.sum((vec_s - vec_t) ** 2, dim=1))
        V_bond = torch.sum(kb * (dis - b0) ** 2)

        # Angle
        angle, kt, t0 = self.angle_k_t["angle"], self.angle_k_t["kt"], self.angle_k_t["t0"]
        t0 = torch.deg2rad(t0)
        vec_ba = x[angle[:, 0]] - x[angle[:, 1]]
        vec_bc = x[angle[:, 2]] - x[angle[:, 1]]
        norm_ba = torch.linalg.norm(vec_ba, dim=1)
        norm_bc = torch.linalg.norm(vec_bc, dim=1)
        cos_t = (vec_ba * vec_bc).sum(dim=1) / (norm_ba * norm_bc)
        cos_t = torch.clamp(cos_t, -0.9999, 0.9999)
        theta = torch.acos(cos_t)
        V_angle = torch.sum(kt * (theta - t0) ** 2)

        # UreyBradley
        ureybradley, ku, s0 = self.ureybradley_k_s["ureybradley"], self.ureybradley_k_s["ku"], self.ureybradley_k_s["s0"]
        vec_s = x[ureybradley[:, 0]]
        vec_t = x[ureybradley[:, 1]]
        dis = torch.sqrt(torch.sum((vec_s - vec_t) ** 2, dim=1))
        V_ureybradley = torch.sum(ku * (dis - s0) ** 2)

        # Dihedral
        dihedral, kd, nd, d0 = self.dihedral_k_n_d["dihedral"], self.dihedral_k_n_d["kd"], self.dihedral_k_n_d["nd"], self.dihedral_k_n_d["d0"]
        d0 = torch.deg2rad(d0)
        vec_ab = x[dihedral[:, 1]] - x[dihedral[:, 0]]
        vec_bc = x[dihedral[:, 2]] - x[dihedral[:, 1]]
        vec_cd = x[dihedral[:, 3]] - x[dihedral[:, 2]]
        norm_abc = torch.cross(vec_ab, vec_bc, dim=1)
        norm_bcd = torch.cross(vec_bc, vec_cd, dim=1)
        cos_d = (norm_abc * norm_bcd).sum(dim=1) / (torch.linalg.norm(norm_abc, dim=1) * torch.linalg.norm(norm_bcd, dim=1))
        cos_d = torch.clamp(cos_d, -0.9999, 0.9999)
        delta = torch.acos(cos_d)
        delta = delta * torch.sign((torch.cross(norm_abc, norm_bcd) * vec_bc).sum(dim=1))
        V_dihedral = torch.sum(kd * (1 + torch.cos(nd * delta - d0)))

        # Improper
        improper, kp, p0 = self.improper_k_p["improper"], self.improper_k_p["kp"], self.improper_k_p["p0"]
        p0 = torch.deg2rad(p0)
        vec_ab = x[improper[:, 1]] - x[improper[:, 0]]
        vec_bc = x[improper[:, 2]] - x[improper[:, 1]]
        vec_cd = x[improper[:, 3]] - x[improper[:, 2]]
        norm_abc = torch.cross(vec_ab, vec_bc, dim=1)
        norm_bcd = torch.cross(vec_bc, vec_cd, dim=1)
        cos_p = (norm_abc * norm_bcd).sum(dim=1) / (torch.linalg.norm(norm_abc, dim=1) * torch.linalg.norm(norm_bcd, dim=1))
        cos_p = torch.clamp(cos_p, -0.9999, 0.9999)
        psi = torch.acos(cos_p)
        psi = psi * torch.sign((torch.cross(norm_abc, norm_bcd) * vec_bc).sum(dim=1))
        V_improper = torch.sum(kp * (psi - p0) ** 2)

        # Vanderwaals
        eps_ij, rmin_ij = self.vanderwaals_e_r["eps_ij"], self.vanderwaals_e_r["rmin_ij"]
        dist = torch.cdist(x, x) + torch.eye(x.size()[0]).to(x.device)
        rod = rmin_ij / dist
        vdw_potential = eps_ij * (torch.pow(rod, 12) - 2 * torch.pow(rod, 6))
        i_upper = torch.triu_indices(vdw_potential.size(0), vdw_potential.size(1), offset=1).to(x.device)
        V_vanderwaals = vdw_potential[i_upper[0], i_upper[1]].sum()

        # Electro
        kqq_ij = self.electro_kqq["kqq_ij"]
        dist = torch.cdist(x, x) + torch.eye(x.size()[0]).to(x.device)
        ele_potential = kqq_ij / dist
        i_upper = torch.triu_indices(ele_potential.size(0), ele_potential.size(1), offset=1).to(x.device)
        V_electro = ele_potential[i_upper[0], i_upper[1]].sum()

        # Central
        V_central = 0
        if central == True:
            dist = torch.sqrt(torch.sum((x - torch.mean(x, dim=0)) ** 2, dim=1))
            cen_potential = 0.1 * dist ** 2
            V_central = torch.sum(cen_potential)    

        V = V_bond + V_angle + V_ureybradley + V_dihedral + V_improper + V_vanderwaals + V_electro + V_central
        return V
    
    def initialize(self, atom_coordinates, device):
        self.init()

        # Build atom mass
        self.mass_prm = pd.read_csv("parameter/ff_mass", sep='\s+', header=None, names=["Atom", "m"])
        m = []
        for x in atom_coordinates.atom_types:
            matching_row = self.mass_prm[(self.mass_prm['Atom']==x)]
            m.append(matching_row['m'].iloc[0])
        self.mass_m = {
            "m": torch.tensor(m).to(device)
        }

        # Build bond potential
        self.bond_prm = pd.read_csv("parameter/ff_bond", sep='\s+', header=None, names=["Atom1", "Atom2", "Kb", "b0"])
        bond, kb, b0 = atom_coordinates.bonds, [], []

        for x, y in atom_coordinates.bonds:
            A1 = atom_coordinates.atom_types[x]
            A2 = atom_coordinates.atom_types[y]
            matching_row = self.bond_prm[(self.bond_prm['Atom1']==A1) & (self.bond_prm['Atom2']==A2)]
            if len(matching_row) == 0:
                matching_row = self.bond_prm[(self.bond_prm['Atom1']==A2) & (self.bond_prm['Atom2']==A1)]
            kb.append(matching_row['Kb'].iloc[0])
            b0.append(matching_row['b0'].iloc[0])
        self.bond_k_b = {
            "bond": torch.tensor(bond).to(device),
            "kb": torch.tensor(kb).to(device),
            "b0": torch.tensor(b0).to(device)
        }

        # Build angle potential
        self.angle_prm = pd.read_csv("parameter/ff_angle", sep='\s+', header=None, names=["Atom1", "Atom2", "Atom3", "Kt", "t0"])
        angle, kt, t0 = atom_coordinates.angles, [], []
        for x, y, z in atom_coordinates.angles:
            A1 = atom_coordinates.atom_types[x]
            A2 = atom_coordinates.atom_types[y]
            A3 = atom_coordinates.atom_types[z]
            matching_row = self.angle_prm[(self.angle_prm['Atom1']==A1) & (self.angle_prm['Atom2']==A2) & (self.angle_prm['Atom3']==A3)]
            if len(matching_row) == 0:
                matching_row = self.angle_prm[(self.angle_prm['Atom1']==A3) & (self.angle_prm['Atom2']==A2) & (self.angle_prm['Atom3']==A1)]
            kt.append(matching_row['Kt'].iloc[0])
            t0.append(matching_row['t0'].iloc[0])
        self.angle_k_t = {
            "angle": torch.tensor(angle).to(device),
            "kt": torch.tensor(kt).to(device),
            "t0": torch.tensor(t0).to(device)
        }

        # Build ureybradley potential
        self.ureybradley_prm = pd.read_csv("parameter/ff_ureybradley", sep='\s+', header=None, names=["Atom1", "Atom2", "Atom3", "Ku", "s0"])
        ureybradley, ku, s0 = [], [], []
        for x, y, z in atom_coordinates.angles:
            A1 = atom_coordinates.atom_types[x]
            A2 = atom_coordinates.atom_types[y]
            A3 = atom_coordinates.atom_types[z]
            matching_row = self.ureybradley_prm[(self.ureybradley_prm['Atom1']==A1) & (self.ureybradley_prm['Atom2']==A2) & (self.ureybradley_prm['Atom3']==A3)]
            if len(matching_row) == 0:
                matching_row = self.ureybradley_prm[(self.ureybradley_prm['Atom1']==A3) & (self.ureybradley_prm['Atom2']==A2) & (self.ureybradley_prm['Atom3']==A1)]
            if len(matching_row) > 0:
                ureybradley.append([x, z])
                ku.append(matching_row['Ku'].iloc[0])
                s0.append(matching_row['s0'].iloc[0])
        self.ureybradley_k_s = {
            "ureybradley": torch.tensor(ureybradley).to(device),
            "ku": torch.tensor(ku).to(device),
            "s0": torch.tensor(s0).to(device)
        }

        # Build dihedral potential
        self.dihedral_prm = pd.read_csv("parameter/ff_dihedral", sep='\s+', header=None, names=["Atom1", "Atom2", "Atom3", "Atom4", "Kd", "nd", "d0"])
        dihedral, kd, nd, d0 = [], [], [], []
        for x, y, z, t in atom_coordinates.dihedrals:
            A1 = atom_coordinates.atom_types[x]
            A2 = atom_coordinates.atom_types[y]
            A3 = atom_coordinates.atom_types[z]
            A4 = atom_coordinates.atom_types[t]
            matching_row = self.dihedral_prm[(self.dihedral_prm['Atom1']==A1) & (self.dihedral_prm['Atom2']==A2) & (self.dihedral_prm['Atom3']==A3) & (self.dihedral_prm['Atom4']==A4)]
            if len(matching_row) == 0:
                matching_row = self.dihedral_prm[(self.dihedral_prm['Atom1']==A4) & (self.dihedral_prm['Atom2']==A3) & (self.dihedral_prm['Atom3']==A2) & (self.dihedral_prm['Atom4']==A1)]
            if len(matching_row) == 0:
                matching_row = self.dihedral_prm[(self.dihedral_prm['Atom1']=='X') & (self.dihedral_prm['Atom2']==A2) & (self.dihedral_prm['Atom3']==A3) & (self.dihedral_prm['Atom4']=='X')]
            if len(matching_row) == 0:
                matching_row = self.dihedral_prm[(self.dihedral_prm['Atom1']=='X') & (self.dihedral_prm['Atom2']==A3) & (self.dihedral_prm['Atom3']==A2) & (self.dihedral_prm['Atom4']=='X')]              
            for id , row in matching_row.iterrows():
                dihedral.append([x, y, z, t])
                kd.append(row['Kd'])
                nd.append(row['nd'])
                d0.append(row['d0'])
        self.dihedral_k_n_d = {
            "dihedral": torch.tensor(dihedral).to(device),
            "kd": torch.tensor(kd).to(device),
            "nd": torch.tensor(nd).to(device),
            "d0": torch.tensor(d0).to(device)
        }

        # Build improper potential
        self.improper_prm = pd.read_csv("parameter/ff_improper", sep='\s+', header=None, names=["Atom1", "Atom2", "Atom3", "Atom4", "Kp", "p0"])
        improper, kp, p0 = [], [], []
        for x, y, z, t in atom_coordinates.impropers:
            A1 = atom_coordinates.atom_types[x]
            A2 = atom_coordinates.atom_types[y]
            A3 = atom_coordinates.atom_types[z]
            A4 = atom_coordinates.atom_types[t]
            matching_row = self.improper_prm[(self.improper_prm['Atom1']==A1) & (self.improper_prm['Atom2']==A2) & (self.improper_prm['Atom3']==A3) & (self.improper_prm['Atom4']==A4)]
            if len(matching_row) == 0:
                matching_row = self.improper_prm[(self.improper_prm['Atom1']==A4) & (self.improper_prm['Atom2']==A3) & (self.improper_prm['Atom3']==A2) & (self.improper_prm['Atom4']==A1)]
            if len(matching_row) == 0:
                matching_row = self.improper_prm[(self.improper_prm['Atom1']==A1) & (self.improper_prm['Atom2']=='X') & (self.improper_prm['Atom3']=='X') & (self.improper_prm['Atom4']==A4)]
            if len(matching_row) == 0:
                matching_row = self.improper_prm[(self.improper_prm['Atom1']==A4) & (self.improper_prm['Atom2']=='X') & (self.improper_prm['Atom3']=='X') & (self.improper_prm['Atom4']==A1)]
            for id , row in matching_row.iterrows():
                improper.append([x, y, z, t])
                kp.append(row['Kp'])
                p0.append(row['p0'])
        self.improper_k_p = {
            "improper": torch.tensor(improper).to(device),
            "kp": torch.tensor(kp).to(device),
            "p0": torch.tensor(p0).to(device)
        }

        # Build vanderwaals potential
        self.vanderwaals_prm = pd.read_csv("parameter/ff_vanderwaals", sep='\s+', header=None, names=["Atom", "epsilon", "Rmin/2"])
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
            "eps_ij": eps_ij.to(device),
            "rmin_ij": rmin_ij.to(device)
        }

        # Build electro potential
        self.electro_prm = {}
        current_key = None
        current_data = []
        with open("parameter/ff_electro", 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if current_key is not None and current_data:
                        self.electro_prm[current_key] = pd.DataFrame(current_data, columns=["Atom", "Q"])
                    current_key = line[1:].strip()
                    current_data = []
                else:
                    parts = line.split()
                    if len(parts) >= 2:
                        current_data.append([parts[0], float(parts[1])])
            if current_key is not None and current_data:
                self.electro_prm[current_key] = pd.DataFrame(current_data, columns=["Atom", "Q"])
        qs = []
        for id, aa, A in zip(atom_coordinates.res_num, atom_coordinates.amino_acids, atom_coordinates.atom_names):
            aa_table = self.electro_prm[aa]
            if id == 0:
                if aa == "G":
                    aa_table = pd.concat([aa_table, self.electro_prm["GLYP"]])
                if aa == "P":
                    aa_table = pd.concat([aa_table, self.electro_prm["PROP"]])
                else:
                    aa_table = pd.concat([aa_table, self.electro_prm["NTER"]])
                aa_table = aa_table.drop_duplicates(subset='Atom', keep='last').reset_index(drop=True)
            if id == atom_coordinates.res_num[-1]:
                aa_table = pd.concat([aa_table, self.electro_prm["CTER"]])
                aa_table = aa_table.drop_duplicates(subset='Atom', keep='last').reset_index(drop=True)
            matching_row = aa_table[(aa_table['Atom']==A)]
            qs.append(matching_row["Q"].iloc[0])
        qs = torch.tensor(qs)
        kqq_ij = 322 * qs[:, None] * qs
        for bond in atom_coordinates.bonds:
            i, j = bond
            kqq_ij[i, j] = 0
            kqq_ij[j, i] = 0
        for angle in atom_coordinates.angles:
            i, j = angle[0], angle[2]
            kqq_ij[i, j] = 0
            kqq_ij[j, i] = 0
        self.electro_kqq = {
            "kqq_ij": kqq_ij.to(device)
        }

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
    
    def gradient_descent(self, central=False):
        self.atom_optimizer.zero_grad()
        self.potential_energy = torch.sum(self.force_field(self.atom_coordinates(), central))
        self.potential_energy.backward()
        torch.nn.utils.clip_grad_norm_(self.atom_coordinates.coordinates, 100)
        self.atom_optimizer.step()
        self.atom_coordinates.addNoise(vol=0.003, damp=0.99)
    
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
        if self.step < 1000:
            self.gradient_descent(True)
        else:
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
    for i in range(100):
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

