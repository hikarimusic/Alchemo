import torch
import torch.nn as nn
import pandas as pd 

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

    def forward(self, x):
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

        # Box
        dist = torch.sqrt(torch.sum((x - torch.mean(x, dim=0)) ** 2, dim=1))
        box_potential = 0.1 * (dist ** 2)
        V_box = torch.sum(box_potential[dist>20])

        V = V_bond + V_angle + V_ureybradley + V_dihedral + V_improper + V_vanderwaals + V_electro + V_box
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
            if id == len(atom_coordinates.protein_sequence)-1:
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