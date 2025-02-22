from flask import Flask, request, jsonify
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem
import traceback
import torch
import torch.nn as nn
import numpy as np
from typing import List, Tuple, Dict
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class AtomTyper:
    """Handles atom typing for the Vina scoring function."""
    
    VDW_RADII = {
        'C_H': 1.9, 'C_P': 1.9, 'N_P': 1.8, 'N_D': 1.8, 'N_A': 1.8, 'N_DA': 1.8,
        'O_P': 1.7, 'O_D': 1.7, 'O_A': 1.7, 'O_DA': 1.7, 'S_P': 2.0, 'P_P': 2.1,
        'F_H': 1.5, 'Cl_H': 1.8, 'Br_H': 2.0, 'I_H': 2.2, 'Si': 2.2, 'At': 2.3,
        'Met_D': 1.2, 'C_H_CG0': 1.9, 'C_P_CG0': 1.9, 'C_H_CG1': 1.9, 'C_P_CG1': 1.9,
        'C_H_CG2': 1.9, 'C_P_CG2': 1.9, 'C_H_CG3': 1.9, 'C_P_CG3': 1.9
    }
    
    @staticmethod
    def get_atom_type(atom: Chem.Atom) -> str:
        """Determine atom type based on RDKit atom."""
        atomic_num = atom.GetAtomicNum()
        is_aromatic = atom.GetIsAromatic()
        
        if atomic_num == 6:  # Carbon
            return 'C_H' if not is_aromatic else 'C_P'
        elif atomic_num == 7:  # Nitrogen
            if atom.GetTotalNumHs() > 0:
                return 'N_D' if atom.GetTotalNumHs() == 1 else 'N_DA'
            return 'N_A'
        elif atomic_num == 8:  # Oxygen
            if atom.GetTotalNumHs() > 0:
                return 'O_D'
            return 'O_A'
        elif atomic_num == 16:  # Sulfur
            return 'S_P'
        elif atomic_num == 15:  # Phosphorus
            return 'P_P'
        elif atomic_num == 9:   # Fluorine
            return 'F_H'
        elif atomic_num == 17:  # Chlorine
            return 'Cl_H'
        elif atomic_num == 35:  # Bromine
            return 'Br_H'
        elif atomic_num == 53:  # Iodine
            return 'I_H'
        elif atomic_num == 14:  # Silicon
            return 'Si'
        elif atomic_num == 85:  # Astatine
            return 'At'
        elif atomic_num in [26, 29, 30]:  # Metals (Fe, Cu, Zn)
            return 'Met_D'
        
        return 'C_H'  # Default to carbon if unknown

class DrugConformation(nn.Module):
    def __init__(self, mol, n_conformers=100):
        super().__init__()
        self.mol = mol
        self.n_conformers = n_conformers
        
        # Find rotatable bonds
        self.rotatable_bonds = self._find_rotatable_bonds()
        n_rotatable = len(self.rotatable_bonds)
        
        # Parameters: [n_conformers, n_params]
        n_params = 6 + n_rotatable  # 3 for position, 3 for orientation, rest for rotatable bonds
        
        # Initialize parameters
        self.conf_params = nn.Parameter(torch.randn(n_conformers, n_params))
        
        # Get atom types
        self.atom_types = [AtomTyper.get_atom_type(atom) for atom in mol.GetAtoms()]
        self.vdw_radii = torch.tensor([AtomTyper.VDW_RADII[atype] for atype in self.atom_types])
        
    def _find_rotatable_bonds(self) -> List[Tuple[int, int]]:
        rotatable = []
        for bond in self.mol.GetBonds():
            if bond.IsRotatable():
                rotatable.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        return rotatable

    def _compute_atom_coordinates(self) -> torch.Tensor:
        device = self.conf_params.device
        conf = self.mol.GetConformer()
        coords = torch.tensor(conf.GetPositions(), device=device)
        coords = coords.unsqueeze(0).repeat(self.n_conformers, 1, 1)
        
        # Apply transformations in parallel
        # Translation
        translations = self.conf_params[:, :3].unsqueeze(1)
        coords = coords + translations
        
        # Rotation (using quaternions)
        angles = self.conf_params[:, 3:6]
        quats = self._euler_to_quaternion(angles)
        coords = self._rotate_points_by_quaternions(coords, quats)
        
        # Rotatable bonds
        for bond_idx, (atom1, atom2) in enumerate(self.rotatable_bonds):
            angles = self.conf_params[:, 6 + bond_idx]
            coords = self._rotate_branches_parallel(coords, atom1, atom2, angles)
        
        return coords

    def _euler_to_quaternion(self, angles: torch.Tensor) -> torch.Tensor:
        x, y, z = angles.chunk(3, dim=1)
        cx, sx = torch.cos(x/2), torch.sin(x/2)
        cy, sy = torch.cos(y/2), torch.sin(y/2)
        cz, sz = torch.cos(z/2), torch.sin(z/2)
        
        qw = cx*cy*cz + sx*sy*sz
        qx = sx*cy*cz - cx*sy*sz
        qy = cx*sy*cz + sx*cy*sz
        qz = cx*cy*sz - sx*sy*cz
        
        return torch.cat([qw, qx, qy, qz], dim=1)

    def _rotate_points_by_quaternions(self, points: torch.Tensor, quats: torch.Tensor) -> torch.Tensor:
        """Rotate points by quaternions for all conformers in parallel."""
        # points: [n_conf, n_points, 3]
        # quats: [n_conf, 4]
        
        qw, qx, qy, qz = quats.chunk(4, dim=1)  # each [n_conf, 1]
        
        # Prepare quaternion rotation matrix elements
        R00 = 1 - 2*qy*qy - 2*qz*qz
        R01 = 2*qx*qy - 2*qz*qw
        R02 = 2*qx*qz + 2*qy*qw
        R10 = 2*qx*qy + 2*qz*qw
        R11 = 1 - 2*qx*qx - 2*qz*qz
        R12 = 2*qy*qz - 2*qx*qw
        R20 = 2*qx*qz - 2*qy*qw
        R21 = 2*qy*qz + 2*qx*qw
        R22 = 1 - 2*qx*qx - 2*qy*qy
        
        # Construct rotation matrices [n_conf, 3, 3]
        R = torch.stack([
            torch.cat([R00, R01, R02], dim=1),
            torch.cat([R10, R11, R12], dim=1),
            torch.cat([R20, R21, R22], dim=1)
        ], dim=2)
        
        # Apply rotation to all points simultaneously
        return torch.bmm(points, R.transpose(1, 2))

    def _rotate_branches_parallel(self, coords: torch.Tensor, atom1: int, atom2: int, angles: torch.Tensor) -> torch.Tensor:
        """Rotate branches around bonds for all conformers in parallel."""
        # coords: [n_conf, n_atoms, 3]
        # angles: [n_conf]
        
        # Get branch atoms
        branch_atoms = self._get_branch_atoms(atom2, atom1)
        if not branch_atoms:
            return coords
            
        # Create rotation axis (bond vector)
        bond_vectors = coords[:, atom2] - coords[:, atom1]  # [n_conf, 3]
        bond_vectors = bond_vectors / torch.norm(bond_vectors, dim=1, keepdim=True)
        
        # Create rotation matrices for all conformers using Rodrigues rotation formula
        angles = angles.unsqueeze(1).unsqueeze(2)  # [n_conf, 1, 1]
        cos_angles = torch.cos(angles)
        sin_angles = torch.sin(angles)
        
        # Prepare cross product matrices
        k = bond_vectors  # [n_conf, 3]
        K = torch.zeros(self.n_conformers, 3, 3, device=coords.device)
        K[:, 0, 1] = -k[:, 2]
        K[:, 0, 2] = k[:, 1]
        K[:, 1, 0] = k[:, 2]
        K[:, 1, 2] = -k[:, 0]
        K[:, 2, 0] = -k[:, 1]
        K[:, 2, 1] = k[:, 0]
        
        # Construct rotation matrices using Rodrigues formula: R = I + sin(θ)K + (1-cos(θ))K²
        I = torch.eye(3, device=coords.device).unsqueeze(0)
        kkt = torch.bmm(k.unsqueeze(2), k.unsqueeze(1))  # [n_conf, 3, 3]
        R = (I * cos_angles + 
             K * sin_angles + 
             kkt * (1 - cos_angles))  # [n_conf, 3, 3]
        
        # Apply rotation to branch atoms only
        rel_coords = coords[:, branch_atoms] - coords[:, atom1:atom1+1]
        rotated_coords = torch.bmm(rel_coords, R.transpose(1, 2))
        coords = coords.clone()
        coords[:, branch_atoms] = rotated_coords + coords[:, atom1:atom1+1]
        
        return coords

class VinaScoringFunction:
    def __init__(self):
        self.weights = torch.tensor([
            -0.0356,  # gauss1
            -0.00516, # gauss2
            0.840,    # repulsion
            -0.0351,  # hydrophobic
            -0.587,   # hydrogen bonding
            0.0585    # Nrot weight
        ])

    def compute_score(self, ligand_coords: torch.Tensor, protein_coords: torch.Tensor, 
                     ligand_types: List[str], protein_types: List[str], n_rot: int) -> torch.Tensor:
        # Get VDW radii for all atoms
        lig_vdw = torch.tensor([AtomTyper.VDW_RADII[t] for t in ligand_types])
        prot_vdw = torch.tensor([AtomTyper.VDW_RADII[t] for t in protein_types])
        
        # Compute pairwise distances
        dists = torch.cdist(ligand_coords, protein_coords)
        
        # Compute surface distances
        surface_dists = dists - (lig_vdw.unsqueeze(1) + prot_vdw.unsqueeze(0))
        
        # Compute scoring terms
        gauss1 = torch.exp(-(surface_dists/0.5)**2)
        gauss2 = torch.exp(-((surface_dists-3.0)/2.0)**2)
        repulsion = torch.where(surface_dists < 0, surface_dists**2, torch.zeros_like(surface_dists))
        
        # Compute final score
        score = (self.weights[0] * gauss1 +
                self.weights[1] * gauss2 +
                self.weights[2] * repulsion)
        
        # Sum over all atom pairs and apply rotatable bond penalty
        total_score = score.sum(dim=(1, 2)) / (1 + self.weights[5] * n_rot)
        
        return total_score

app = Flask(__name__)
CORS(app)

# Store current state
current_state = {
    'pdb_data': None,
    'drug_data': None,
    'smile_sequence': None,
    'protein_mol': None,
    'ligand_mol': None
}

@app.route('/api/load_protein', methods=['POST'])
def load_protein():
    try:
        data = request.json
        pdb_data = data.get('pdb_data')
        
        if not pdb_data:
            return jsonify({'error': 'No PDB data provided'}), 400
        
        # Store the PDB data
        current_state['pdb_data'] = pdb_data
        
        return jsonify({
            'message': 'Protein loaded successfully',
            'pdb_data': pdb_data
        })
        
    except Exception as e:
        logger.error(f"Error in load_protein: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/api/load_drug', methods=['POST'])
def load_drug():
    try:
        data = request.json
        smile_sequence = data.get('smile_sequence')
        
        if not smile_sequence:
            return jsonify({'error': 'No SMILE sequence provided'}), 400
        
        # Convert SMILE to 3D structure
        mol = Chem.MolFromSmiles(smile_sequence)
        if mol is None:
            return jsonify({'error': 'Invalid SMILE sequence'}), 400
        
        # Generate 3D conformation
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Store molecule and convert to PDB
        current_state['ligand_mol'] = mol
        pdb_data = Chem.MolToPDBBlock(mol)
        current_state['drug_data'] = pdb_data
        current_state['smile_sequence'] = smile_sequence
        
        return jsonify({
            'message': 'Drug loaded successfully',
            'drug_data': pdb_data
        })
        
    except Exception as e:
        logger.error(f"Error in load_drug: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/api/start_docking', methods=['POST'])
def start_docking():
    try:
        if not current_state['ligand_mol'] or not current_state['pdb_data']:
            return jsonify({'error': 'Both protein and ligand must be loaded'}), 400
        
        # Initialize models
        drug_conf = DrugConformation(current_state['ligand_mol'], n_conformers=100)
        scoring_function = VinaScoringFunction()
        
        # Setup optimizer
        optimizer = torch.optim.Adam(drug_conf.parameters(), lr=0.01)
        
        logger.info("Starting optimization...")
        
        # Optimization loop
        for i in range(1000):
            optimizer.zero_grad()
            
            # Get conformer coordinates
            ligand_coords = drug_conf._compute_atom_coordinates()
            
            # Compute score
            score = scoring_function.compute_score(
                ligand_coords,
                torch.tensor(current_state['protein_coords']),
                drug_conf.atom_types,
                current_state['protein_types'],
                len(drug_conf.rotatable_bonds)
            )
            
            # Compute loss and backpropagate
            loss = score.mean()
            loss.backward()
            optimizer.step()
            
            if i % 100 == 0:
                logger.info(f"Iteration {i}, Score: {loss.item():.3f}")
        
        return jsonify({'message': 'Docking completed successfully'})
        
    except Exception as e:
        logger.error(f"Error in start_docking: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/api/get_current_state', methods=['GET'])
def get_current_state():
    try:
        return jsonify({
            'pdb_data': current_state['pdb_data'],
            'drug_data': current_state['drug_data'],
            'smile_sequence': current_state['smile_sequence']
        })
        
    except Exception as e:
        logger.error(f"Error in get_current_state: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True, port=5000)