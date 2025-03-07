from flask import Flask, request, jsonify
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem, BRICS
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
    """Enhanced atom typing based on X-score for the Vina scoring function."""
    
    VDW_RADII = {
        'C': 1.9,  # Carbon
        'N': 1.8,  # Nitrogen
        'O': 1.7,  # Oxygen
        'S': 2.0,  # Sulfur
        'P': 2.1,  # Phosphorus
        'F': 1.5,  # Fluorine
        'Cl': 1.8, # Chlorine
        'Br': 2.0, # Bromine
        'I': 2.2,  # Iodine
        'Si': 2.2, # Silicon
        'Fe': 1.2, # Iron
        'Cu': 1.2, # Copper
        'Zn': 1.2, # Zinc
        'Mg': 1.2, # Magnesium
        'Ca': 1.2, # Calcium
        'Mn': 1.2, # Manganese
        'H': 1.2   # Hydrogen
    }
    
    @staticmethod
    def get_atom_type(atom: Chem.Atom) -> str:
        """Determine atom type based on element only."""
        elements = {
            1: 'H',   # Hydrogen
            6: 'C',   # Carbon
            7: 'N',   # Nitrogen
            8: 'O',   # Oxygen
            9: 'F',   # Fluorine
            12: 'Mg', # Magnesium
            15: 'P',  # Phosphorus
            16: 'S',  # Sulfur
            17: 'Cl', # Chlorine
            20: 'Ca', # Calcium
            25: 'Mn', # Manganese
            26: 'Fe', # Iron
            29: 'Cu', # Copper
            30: 'Zn', # Zinc
            35: 'Br', # Bromine
            53: 'I'   # Iodine
        }
        
        atomic_num = atom.GetAtomicNum()
        return elements.get(atomic_num, 'C')  # Default to carbon if unknown
    
    @staticmethod
    def is_hydrophobic(atom: Chem.Atom) -> bool:
        """Determine if atom is hydrophobic following X-score rules."""
        # Carbon atoms are hydrophobic unless they're bonded to polar atoms
        if atom.GetAtomicNum() == 6:  # Carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() in [7, 8, 15, 16]:  # N, O, P, S
                    return False
            return True
        
        # Halogens (except F) are considered hydrophobic
        if atom.GetAtomicNum() in [17, 35, 53]:  # Cl, Br, I
            return True
        
        return False
    
    @staticmethod
    def is_hbond_acceptor(atom: Chem.Atom) -> bool:
        """Determine if atom is a hydrogen bond acceptor."""
        if atom.GetAtomicNum() in [7, 8]:  # N, O
            # Check if N or O has lone pair (not fully saturated with bonds)
            if atom.GetExplicitValence() < atom.GetTotalValence():
                return True
        
        return False
    
    @staticmethod
    def is_hbond_donor(atom: Chem.Atom) -> bool:
        """Determine if atom is a hydrogen bond donor."""
        # Oxygen or nitrogen with at least one hydrogen
        if atom.GetAtomicNum() in [7, 8]:  # N, O
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:  # Hydrogen
                    return True
        
        # Metals are treated as hydrogen bond donors in Vina
        if atom.GetAtomicNum() in [12, 20, 25, 26, 29, 30]:  # Mg, Ca, Mn, Fe, Cu, Zn
            return True
        
        return False
    
    @staticmethod
    def precompute_atom_properties(mol: Chem.Mol) -> tuple:
        """Precompute atom properties for a molecule to avoid repeated calculations."""
        n_atoms = mol.GetNumAtoms()
        atom_types = []
        is_hydrophobic = []
        is_hbond_donor = []
        is_hbond_acceptor = []
        
        for i in range(n_atoms):
            atom = mol.GetAtomWithIdx(i)
            atom_types.append(AtomTyper.get_atom_type(atom))
            is_hydrophobic.append(AtomTyper.is_hydrophobic(atom))
            is_hbond_donor.append(AtomTyper.is_hbond_donor(atom))
            is_hbond_acceptor.append(AtomTyper.is_hbond_acceptor(atom))
        
        return atom_types, is_hydrophobic, is_hbond_donor, is_hbond_acceptor

class VinaScoringFunction:
    def __init__(self, device, mol):
        self.device = device
        self.weights = torch.tensor([
            -0.0356,  # gauss1
            -0.00516, # gauss2
            0.840,    # repulsion
            -0.0351,  # hydrophobic
            -0.587,   # hydrogen bonding
            0.0585    # Nrot weight
        ], dtype=torch.float32, device=device)
        
        # Distance cutoff for all interactions (8Å)
        self.cutoff_distance = 8.0
        
        # Precompute 1-4 interaction exclusion mask
        self.exclusion_mask = self._compute_exclusion_mask(mol, device)
        
        # Precompute atom properties
        atom_types, is_hydrophobic, is_hbond_donor, is_hbond_acceptor = AtomTyper.precompute_atom_properties(mol)
        
        # Convert to tensors
        self.ligand_is_hydrophobic = torch.tensor(is_hydrophobic, dtype=torch.bool, device=device)
        self.ligand_is_hbond_donor = torch.tensor(is_hbond_donor, dtype=torch.bool, device=device)
        self.ligand_is_hbond_acceptor = torch.tensor(is_hbond_acceptor, dtype=torch.bool, device=device)
        
    def _compute_exclusion_mask(self, mol, device):
        """
        Precompute mask for excluding 1-4 and closer interactions.
        Returns a boolean tensor where True indicates pairs to include (not exclude).
        """
        n_atoms = mol.GetNumAtoms()
        # Initialize mask where True means we include the interaction
        mask = torch.ones((n_atoms, n_atoms), dtype=torch.bool, device=device)
        
        # Create graph representation for path finding
        graph = {}
        for bond in mol.GetBonds():
            begin = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            if begin not in graph:
                graph[begin] = set()
            if end not in graph:
                graph[end] = set()
            graph[begin].add(end)
            graph[end].add(begin)
            
        # Function to find all paths of length 1-4 between atoms
        def find_short_paths(start, max_length=4):
            excluded = set()
            visited = {start: 0}
            queue = [(start, 0)]
            
            while queue:
                current, length = queue.pop(0)
                if length >= max_length:
                    continue
                    
                for neighbor in graph.get(current, []):
                    new_length = length + 1
                    if neighbor not in visited or visited[neighbor] > new_length:
                        visited[neighbor] = new_length
                        queue.append((neighbor, new_length))
                        if new_length <= 4:  # Include 1-4 and closer interactions
                            excluded.add(neighbor)
            
            return excluded
        
        # Build exclusion mask
        for i in range(n_atoms):
            excluded = find_short_paths(i)
            for j in excluded:
                # Set False for pairs to exclude
                mask[i, j] = False
                mask[j, i] = False
                
        # Set diagonal to False (no self-interactions)
        mask.fill_diagonal_(False)
        
        return mask
    
    def _linear_interpolation(self, x, x1, y1, x2, y2):
        """Linear interpolation between two points (x1,y1) and (x2,y2)."""
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        
    def _compute_intramolecular_score(self, ligand_coords, lig_vdw):
        """
        Compute intramolecular scoring terms for ligand.
        ligand_coords: tensor of shape (n_conformers, n_atoms, 3)
        lig_vdw: tensor of shape (n_atoms,)
        """
        n_conformers = ligand_coords.shape[0]
        n_atoms = ligand_coords.shape[1]
        
        # Expand exclusion mask for batch dimension
        batch_mask = self.exclusion_mask.unsqueeze(0)
        
        # Compute pairwise distances between all ligand atoms for all conformers
        dists = torch.cdist(ligand_coords, ligand_coords)
        
        # Apply distance cutoff (8Å)
        cutoff_mask = dists < self.cutoff_distance
        cutoff_mask = cutoff_mask & batch_mask  # Combine with exclusion mask
        
        # Compute surface distances
        surface_dists = dists - (lig_vdw.unsqueeze(0).unsqueeze(2) + lig_vdw.unsqueeze(0).unsqueeze(1))
        
        # Calculate all scoring terms
        gauss1 = torch.exp(-(surface_dists/0.5)**2)
        gauss2 = torch.exp(-((surface_dists-3.0)/2.0)**2)
        repulsion = torch.where(surface_dists < 0, surface_dists**2, torch.zeros_like(surface_dists))
        
        # Hydrophobic interactions
        # Create (n_conformers, n_atoms, n_atoms) bool tensor indicating hydrophobic pairs
        hydrophobic_pairs = self.ligand_is_hydrophobic.unsqueeze(0).unsqueeze(2) & self.ligand_is_hydrophobic.unsqueeze(0).unsqueeze(1)
        
        # Compute hydrophobic score using linear interpolation
        hydrophobic_score = torch.zeros_like(surface_dists)
        mask_hydrophobic = (surface_dists < 0.5) & hydrophobic_pairs
        hydrophobic_score = torch.where(mask_hydrophobic, torch.ones_like(surface_dists), hydrophobic_score)
        
        mask_interpolate = (surface_dists >= 0.5) & (surface_dists <= 1.5) & hydrophobic_pairs
        interp_values = 1.0 - (surface_dists - 0.5) / (1.5 - 0.5)  # Linear interpolation
        hydrophobic_score = torch.where(mask_interpolate, interp_values, hydrophobic_score)
        
        # For all other hydrophobic pairs with surface distance > 1.5, score remains 0
        
        # Hydrogen bonding
        # Create tensors for donor-acceptor pairs (n_conformers, n_atoms, n_atoms)
        donor_acceptor = self.ligand_is_hbond_donor.unsqueeze(0).unsqueeze(2) & self.ligand_is_hbond_acceptor.unsqueeze(0).unsqueeze(1)
        acceptor_donor = self.ligand_is_hbond_acceptor.unsqueeze(0).unsqueeze(2) & self.ligand_is_hbond_donor.unsqueeze(0).unsqueeze(1)
        hbond_pairs = donor_acceptor | acceptor_donor
        
        # Compute hydrogen bonding score using linear interpolation
        hbond_score = torch.zeros_like(surface_dists)
        mask_hbond = (surface_dists < -0.7) & hbond_pairs
        hbond_score = torch.where(mask_hbond, torch.ones_like(surface_dists), hbond_score)
        
        mask_interpolate = (surface_dists >= -0.7) & (surface_dists <= 0) & hbond_pairs
        interp_values = (-surface_dists) / 0.7  # Linear interpolation from -0.7 to 0
        hbond_score = torch.where(mask_interpolate, interp_values, hbond_score)
        
        # For all other H-bond pairs with surface distance > 0, score remains 0
        
        # Apply weights and sum all terms
        weighted_score = (
            self.weights[0] * gauss1 +
            self.weights[1] * gauss2 +
            self.weights[2] * repulsion +
            self.weights[3] * hydrophobic_score +
            self.weights[4] * hbond_score
        )
        
        # Apply masks and sum
        weighted_score = weighted_score * cutoff_mask.float()
        return weighted_score.sum(dim=(1, 2)) / 2  # Divide by 2 to avoid double counting

    def compute_score(self, ligand_coords: torch.Tensor, protein_coords: torch.Tensor, 
                     lig_vdw: torch.Tensor, prot_vdw: torch.Tensor, 
                     prot_is_hydrophobic: torch.Tensor, 
                     prot_is_hbond_donor: torch.Tensor,
                     prot_is_hbond_acceptor: torch.Tensor,
                     n_rot: int) -> torch.Tensor:
        """
        Compute combined intermolecular and intramolecular Vina-like scoring function.
        
        Parameters:
        -----------
        ligand_coords: tensor of shape (n_conformers, n_atoms, 3)
            Coordinates of ligand atoms
        protein_coords: tensor of shape (n_protein_atoms, 3)
            Coordinates of protein atoms
        lig_vdw: tensor of shape (n_atoms,)
            Van der Waals radii of ligand atoms
        prot_vdw: tensor of shape (n_protein_atoms,)
            Van der Waals radii of protein atoms
        prot_is_hydrophobic: tensor of shape (n_protein_atoms,)
            Boolean tensor indicating which protein atoms are hydrophobic
        prot_is_hbond_donor: tensor of shape (n_protein_atoms,)
            Boolean tensor indicating which protein atoms are hydrogen bond donors
        prot_is_hbond_acceptor: tensor of shape (n_protein_atoms,)
            Boolean tensor indicating which protein atoms are hydrogen bond acceptors
        n_rot: int
            Number of rotatable bonds in ligand
        """
        # Compute distances between ligand and protein atoms
        dists = torch.cdist(ligand_coords, protein_coords)
        
        # Apply distance cutoff (8Å)
        cutoff_mask = dists < self.cutoff_distance
        
        # Compute surface distances
        surface_dists = dists - (lig_vdw.unsqueeze(1) + prot_vdw.unsqueeze(0))
        
        # Calculate steric terms
        gauss1 = torch.exp(-(surface_dists/0.5)**2)
        gauss2 = torch.exp(-((surface_dists-3.0)/2.0)**2)
        repulsion = torch.where(surface_dists < 0, surface_dists**2, torch.zeros_like(surface_dists))
        
        # Hydrophobic interactions
        # Create pairs of hydrophobic atoms
        lig_hydrophobic = self.ligand_is_hydrophobic.unsqueeze(1).expand(-1, protein_coords.shape[0])
        prot_hydrophobic = prot_is_hydrophobic.unsqueeze(0).expand(ligand_coords.shape[1], -1)
        hydrophobic_pairs = lig_hydrophobic & prot_hydrophobic
        hydrophobic_pairs = hydrophobic_pairs.unsqueeze(0).expand(ligand_coords.shape[0], -1, -1)
        
        # Compute hydrophobic score
        hydrophobic_score = torch.zeros_like(surface_dists)
        
        # If surface distance < 0.5Å, score = 1
        mask_hydrophobic = (surface_dists < 0.5) & hydrophobic_pairs
        hydrophobic_score = torch.where(mask_hydrophobic, torch.ones_like(surface_dists), hydrophobic_score)
        
        # If 0.5Å <= surface distance <= 1.5Å, linearly interpolate
        mask_interpolate = (surface_dists >= 0.5) & (surface_dists <= 1.5) & hydrophobic_pairs
        interp_values = 1.0 - (surface_dists - 0.5) / (1.5 - 0.5)  # Linear interpolation
        hydrophobic_score = torch.where(mask_interpolate, interp_values, hydrophobic_score)
        
        # Hydrogen bonding
        # Create pairs for ligand donor - protein acceptor
        lig_donor_prot_acceptor = self.ligand_is_hbond_donor.unsqueeze(1) & prot_is_hbond_acceptor.unsqueeze(0)
        lig_donor_prot_acceptor = lig_donor_prot_acceptor.unsqueeze(0).expand(ligand_coords.shape[0], -1, -1)
        
        # Create pairs for ligand acceptor - protein donor
        lig_acceptor_prot_donor = self.ligand_is_hbond_acceptor.unsqueeze(1) & prot_is_hbond_donor.unsqueeze(0)
        lig_acceptor_prot_donor = lig_acceptor_prot_donor.unsqueeze(0).expand(ligand_coords.shape[0], -1, -1)
        
        # Combine to get all H-bond pairs
        hbond_pairs = lig_donor_prot_acceptor | lig_acceptor_prot_donor
        
        # Compute hydrogen bonding score
        hbond_score = torch.zeros_like(surface_dists)
        
        # If surface distance < -0.7Å, score = 1
        mask_hbond = (surface_dists < -0.7) & hbond_pairs
        hbond_score = torch.where(mask_hbond, torch.ones_like(surface_dists), hbond_score)
        
        # If -0.7Å <= surface distance <= 0Å, linearly interpolate
        mask_interpolate = (surface_dists >= -0.7) & (surface_dists <= 0) & hbond_pairs
        interp_values = (-surface_dists) / 0.7  # Linear interpolation from -0.7 to 0
        hbond_score = torch.where(mask_interpolate, interp_values, hbond_score)
        
        # Apply weights to all terms
        weighted_score = (
            self.weights[0] * gauss1 +
            self.weights[1] * gauss2 +
            self.weights[2] * repulsion +
            self.weights[3] * hydrophobic_score +
            self.weights[4] * hbond_score
        )
        
        # Apply cutoff mask and sum over all atom pairs
        weighted_score = weighted_score * cutoff_mask.float()
        intermolecular_score = weighted_score.sum(dim=(1, 2))
        
        # Compute intramolecular score
        intramolecular_score = self._compute_intramolecular_score(ligand_coords, lig_vdw)
        
        # Add rotatable bond penalty
        rotatable_penalty = self.weights[5] * n_rot
        
        # Combine all scores
        total_score = intermolecular_score + intramolecular_score + rotatable_penalty
        
        return total_score
        
    def precompute_protein_properties(self, protein_mol):
        """
        Precompute protein atom properties for scoring.
        
        Parameters:
        -----------
        protein_mol: RDKit Mol object
            Protein molecule
            
        Returns:
        --------
        tuple of tensors: (is_hydrophobic, is_hbond_donor, is_hbond_acceptor)
            Boolean tensors indicating atom properties
        """
        _, is_hydrophobic, is_hbond_donor, is_hbond_acceptor = AtomTyper.precompute_atom_properties(protein_mol)
        
        # Convert to tensors
        prot_is_hydrophobic = torch.tensor(is_hydrophobic, dtype=torch.bool, device=self.device)
        prot_is_hbond_donor = torch.tensor(is_hbond_donor, dtype=torch.bool, device=self.device)
        prot_is_hbond_acceptor = torch.tensor(is_hbond_acceptor, dtype=torch.bool, device=self.device)
        
        return prot_is_hydrophobic, prot_is_hbond_donor, prot_is_hbond_acceptor

class DrugConformation(nn.Module):
    def __init__(self, mol, protein_coords, n_conformers=100, surface_distance=0.0):
        super().__init__()
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        logger.info(f"DrugConformation using device: {self.device}")
        
        self.mol = mol
        self.n_conformers = n_conformers
        self.protein_coords = torch.tensor(protein_coords, device=self.device, dtype=torch.float32)
        
        # Get atom types and VDW radii for ligand
        self.ligand_types = [AtomTyper.get_atom_type(atom) for atom in mol.GetAtoms()]
        self.ligand_vdw = torch.tensor(
            [AtomTyper.VDW_RADII[atype] for atype in self.ligand_types],
            device=self.device,
            dtype=torch.float32
        )
        
        # Find rotatable bonds and their associated branch atoms
        self.rotatable_bonds, self.branch_atoms = self._find_rotatable_bonds_and_branches()
        n_rotatable = len(self.rotatable_bonds)
        
        # Initialize parameters near protein surface
        n_params = 6 + n_rotatable  # 3 for position, 3 for orientation, rest for rotatable bonds
        initial_params = self._initialize_near_surface(surface_distance)
        
        # Parameters: [n_conformers, n_params]
        self.conf_params = nn.Parameter(initial_params)
        
        self.to(self.device)

    def _initialize_near_surface(self, surface_distance: float) -> torch.Tensor:
        """Initialize conformers near the protein surface."""
        n_rotatable = len(self.rotatable_bonds)
        n_params = 6 + n_rotatable
        
        # Calculate protein center and approximate radius
        protein_center = self.protein_coords.mean(dim=0)
        protein_radius = torch.max(torch.norm(self.protein_coords - protein_center, dim=1))
        
        # Generate random points on a sphere around the protein
        theta = torch.rand(self.n_conformers, device=self.device) * 2 * np.pi
        phi = torch.acos(2 * torch.rand(self.n_conformers, device=self.device) - 1)
        r = protein_radius + surface_distance
        
        x = r * torch.sin(phi) * torch.cos(theta)
        y = r * torch.sin(phi) * torch.sin(theta)
        z = r * torch.cos(phi)
        positions = torch.stack([x, y, z], dim=1) + protein_center
        
        # Random orientations (euler angles)
        orientations = torch.rand(self.n_conformers, 3, device=self.device) * 2 * np.pi
        
        # Random rotatable bond angles
        rot_angles = torch.rand(self.n_conformers, n_rotatable, device=self.device) * 2 * np.pi
        
        # Combine all parameters
        params = torch.cat([positions, orientations, rot_angles], dim=1)
        
        return params.to(torch.float32)

    def _find_rotatable_bonds_and_branches(self) -> Tuple[List[Tuple[int, int]], List[List[int]]]:
        """Find rotatable bonds and precalculate their branch atoms."""
        rotatable = []
        branches = []
        
        for bond in self.mol.GetBonds():
            if (not bond.IsInRing() and 
                bond.GetBondType() == Chem.rdchem.BondType.SINGLE):
                
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                begin_idx = begin_atom.GetIdx()
                end_idx = end_atom.GetIdx()
                
                # Skip if either atom is hydrogen
                if begin_atom.GetAtomicNum() == 1 or end_atom.GetAtomicNum() == 1:
                    continue
                    
                # Skip if both atoms have less than 2 heavy atom neighbors
                if (len([n for n in begin_atom.GetNeighbors() if n.GetAtomicNum() != 1]) < 2 or
                    len([n for n in end_atom.GetNeighbors() if n.GetAtomicNum() != 1]) < 2):
                    continue
                
                # Find branch atoms
                branch_atoms = self._get_branch_atoms(end_idx, begin_idx)
                if branch_atoms:
                    rotatable.append((begin_idx, end_idx))
                    branches.append(branch_atoms)
        
        return rotatable, branches

    def _get_branch_atoms(self, atom2: int, atom1: int) -> List[int]:
        """Find all atoms on the atom2 side of the rotatable bond."""
        visited = set()
        branch_atoms = set()
        
        def dfs(atom_idx: int) -> None:
            if atom_idx in visited:
                return
                
            visited.add(atom_idx)
            if atom_idx != atom1:
                branch_atoms.add(atom_idx)
                
            atom = self.mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx != atom1 and neighbor_idx not in visited:
                    dfs(neighbor_idx)
        
        dfs(atom2)
        return sorted(list(branch_atoms))

    def _euler_to_quaternion(self, angles: torch.Tensor) -> torch.Tensor:
        """Convert Euler angles to quaternions."""
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
        qw, qx, qy, qz = quats.chunk(4, dim=1)
        
        R00 = 1 - 2*qy*qy - 2*qz*qz
        R01 = 2*qx*qy - 2*qz*qw
        R02 = 2*qx*qz + 2*qy*qw
        R10 = 2*qx*qy + 2*qz*qw
        R11 = 1 - 2*qx*qx - 2*qz*qz
        R12 = 2*qy*qz - 2*qx*qw
        R20 = 2*qx*qz - 2*qy*qw
        R21 = 2*qy*qz + 2*qx*qw
        R22 = 1 - 2*qx*qx - 2*qy*qy
        
        R = torch.stack([
            torch.cat([R00, R01, R02], dim=1),
            torch.cat([R10, R11, R12], dim=1),
            torch.cat([R20, R21, R22], dim=1)
        ], dim=2)
        
        return torch.bmm(points, R.transpose(1, 2))

    def _rotate_branch_atoms(self, coords: torch.Tensor, atom1: int, atom2: int, 
                           angles: torch.Tensor, branch_atoms: List[int]) -> torch.Tensor:
        """Rotate branch atoms around bond using precalculated branch list."""
        if not branch_atoms:
            return coords
            
        # Create rotation axis (bond vector)
        bond_vectors = coords[:, atom2] - coords[:, atom1]
        bond_vectors = bond_vectors / torch.norm(bond_vectors, dim=1, keepdim=True)
        
        # Create rotation matrices using Rodrigues formula
        angles = angles.unsqueeze(1).unsqueeze(2)
        cos_angles = torch.cos(angles)
        sin_angles = torch.sin(angles)
        
        k = bond_vectors
        K = torch.zeros(self.n_conformers, 3, 3, device=coords.device, dtype=torch.float32)
        K[:, 0, 1] = -k[:, 2]
        K[:, 0, 2] = k[:, 1]
        K[:, 1, 0] = k[:, 2]
        K[:, 1, 2] = -k[:, 0]
        K[:, 2, 0] = -k[:, 1]
        K[:, 2, 1] = k[:, 0]
        
        I = torch.eye(3, device=coords.device, dtype=torch.float32).unsqueeze(0)
        kkt = torch.bmm(k.unsqueeze(2), k.unsqueeze(1))
        R = (I * cos_angles + 
             K * sin_angles + 
             kkt * (1 - cos_angles))
        
        rel_coords = coords[:, branch_atoms] - coords[:, atom1:atom1+1]
        rotated_coords = torch.bmm(rel_coords, R.transpose(1, 2))
        coords = coords.clone()
        coords[:, branch_atoms] = rotated_coords + coords[:, atom1:atom1+1]
        
        return coords

    def _compute_atom_coordinates(self) -> torch.Tensor:
        """Compute 3D coordinates for all atoms in all conformers."""
        conf = self.mol.GetConformer()
        coords = torch.tensor(conf.GetPositions(), device=self.device, dtype=torch.float32)
        coords = coords.unsqueeze(0).repeat(self.n_conformers, 1, 1)
                        
        # Apply rotatable bond rotations
        for bond_idx, ((atom1, atom2), branch_atoms) in enumerate(zip(self.rotatable_bonds, self.branch_atoms)):
            angles = self.conf_params[:, 6 + bond_idx]
            coords = self._rotate_branch_atoms(coords, atom1, atom2, angles, branch_atoms)

        # Apply global rotation
        angles = self.conf_params[:, 3:6]
        quats = self._euler_to_quaternion(angles)
        coords = self._rotate_points_by_quaternions(coords, quats)

        # Apply translation
        translations = self.conf_params[:, :3].unsqueeze(1)
        coords = coords + translations

        return coords

    def forward(self) -> torch.Tensor:
        """Forward pass returns the computed 3D coordinates for all conformers."""
        return self._compute_atom_coordinates()

app = Flask(__name__)
CORS(app)

# Store current state
current_state = {
    'pdb_data': None,
    'drug_data': None,
    'smile_sequence': None,
    'protein_mol': None,
    'ligand_mol': None,
    'protein_coords': None,
    'protein_types': None,
    'drug_conf': None,
    'scoring_fn': None,
    'optimizer': None,
    'protein_coords_tensor': None,
    'prot_vdw': None,
    'prot_is_hydrophobic': None,
    'prot_is_hbond_donor': None,
    'prot_is_hbond_acceptor': None,
    'best_scores': None,
    'best_coords': None,
    'step_count': 0,
    'last_improvement_step': 0,
    'stagnation_tolerance': 10,
    'global_best_scores': None,  # Track best scores across all rounds
    'global_best_coords': None,  # Track best coordinates across all rounds
    'round_best_scores': None,   # Track best scores for current round
    'total_steps': 1e10
}

def get_protein_coords_and_types(pdb_data: str) -> Tuple[np.ndarray, List[str]]:
    """Extract protein coordinates and atom types from PDB data."""
    coords = []
    atom_types = []
    
    for line in pdb_data.split('\n'):
        if line.startswith(('ATOM', 'HETATM')):
            try:
                # Extract coordinates
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])
                
                # Extract element symbol (columns 77-78)
                element = line[76:78].strip()
                if not element:  # If element field is empty, use first character of atom name
                    element = line[12:16].strip()[0]
                
                # Map some common PDB element representations
                element_map = {
                    'FE': 'Fe',
                    'CU': 'Cu',
                    'ZN': 'Zn',
                    'BR': 'Br',
                    'CL': 'Cl'
                }
                element = element_map.get(element.upper(), element.capitalize())
                
                atom_types.append(element)
                
            except (ValueError, IndexError):
                continue
                
    return np.array(coords), atom_types

@app.route('/api/load_protein', methods=['POST'])
def load_protein():
    try:
        data = request.json
        pdb_data = data.get('pdb_data')
        
        if not pdb_data:
            return jsonify({'error': 'No PDB data provided'}), 400
        
        # Extract coordinates and atom types
        coords, atom_types = get_protein_coords_and_types(pdb_data)
        
        if len(coords) == 0:
            return jsonify({'error': 'Failed to extract protein coordinates'}), 400
        
        # Store the data
        current_state['pdb_data'] = pdb_data
        current_state['protein_coords'] = coords
        current_state['protein_types'] = atom_types
        
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
        mol = Chem.AddHs(mol)  # Add hydrogens explicitly for better scoring
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

@app.route('/api/initialize_docking', methods=['POST'])
def initialize_docking():
    try:
        if not current_state['ligand_mol'] or current_state['protein_coords'] is None:
            return jsonify({'error': 'Both protein and ligand must be loaded'}), 400
            
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        logger.info(f"Using device: {device}")
        
        # Reset all state variables
        current_state.update({
            'drug_conf': None,
            'scoring_fn': None,
            'optimizer': None,
            'protein_coords_tensor': None,
            'prot_vdw': None,
            'prot_is_hydrophobic': None,
            'prot_is_hbond_donor': None,
            'prot_is_hbond_acceptor': None,
            'best_scores': None,
            'best_coords': None,
            'step_count': 0,
            'last_improvement_step': 0,
            'global_best_scores': None,
            'global_best_coords': None,
            'round_best_scores': None,
            'stagnation_tolerance': 10,
            'total_steps': 1e10
        })
        
        # Initialize drug conformation model
        # MODIFIED: Reduced n_conformers from 1000 to 10 for easier tracking
        drug_conf = DrugConformation(
            mol=current_state['ligand_mol'],
            protein_coords=current_state['protein_coords'],
            n_conformers=1000,  # Reduced significantly to make tracking easier
            surface_distance=0.0  # Initial distance from protein surface
        ).to(device)
        
        # Initialize scoring function
        scoring_fn = VinaScoringFunction(device, current_state['ligand_mol'])
        
        # Set up optimizer
        optimizer = torch.optim.SGD([drug_conf.conf_params],
                                lr=0.3,                # Aggressive learning rate
                                momentum=0.98,         # Very high momentum
                                dampening=0,           # No dampening for more aggressive updates
                                weight_decay=0.0,      # Add noise through weight decay
                                nesterov=True)         # More aggressive momentum implementation
        
        # Get protein data
        protein_coords = torch.tensor(
            current_state['protein_coords'],
            device=device,
            dtype=torch.float32
        )
        
        # Get protein VDW radii
        prot_vdw = torch.tensor(
            [AtomTyper.VDW_RADII.get(t, 1.7) for t in current_state['protein_types']],  # Default to oxygen radius if unknown
            device=device,
            dtype=torch.float32
        )
        
        # Get protein properties - simulate from protein types since we don't have a protein_mol
        prot_types = current_state['protein_types']
        is_hydrophobic = [t in ['C', 'Cl', 'Br', 'I'] for t in prot_types]
        is_hbond_donor = [t in ['N', 'O', 'Mg', 'Ca', 'Mn', 'Fe', 'Cu', 'Zn'] for t in prot_types]
        is_hbond_acceptor = [t in ['N', 'O'] for t in prot_types]
        
        prot_is_hydrophobic = torch.tensor(is_hydrophobic, device=device, dtype=torch.bool)
        prot_is_hbond_donor = torch.tensor(is_hbond_donor, device=device, dtype=torch.bool)
        prot_is_hbond_acceptor = torch.tensor(is_hbond_acceptor, device=device, dtype=torch.bool)

        # Initialize best scores and coordinates trackers
        best_scores = float('inf') * torch.ones(drug_conf.n_conformers, device=device)
        global_best_scores = float('inf') * torch.ones(drug_conf.n_conformers, device=device)
        round_best_scores = float('inf') * torch.ones(drug_conf.n_conformers, device=device)
        
        # Store all necessary state
        current_state.update({
            'drug_conf': drug_conf,
            'scoring_fn': scoring_fn,
            'optimizer': optimizer,
            'protein_coords_tensor': protein_coords,
            'prot_vdw': prot_vdw,
            'prot_is_hydrophobic': prot_is_hydrophobic,
            'prot_is_hbond_donor': prot_is_hbond_donor,
            'prot_is_hbond_acceptor': prot_is_hbond_acceptor,
            'best_scores': best_scores,
            'best_coords': None,
            'global_best_scores': global_best_scores,
            'global_best_coords': None,
            'round_best_scores': round_best_scores,
            'step_count': 0,
            'last_improvement_step': 0
        })
        
        return jsonify({'message': 'Docking initialized successfully'})
        
    except Exception as e:
        logger.error(f"Error in initialize_docking: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/api/optimize_docking', methods=['POST'])
def optimize_docking():
    try:
        # Get necessary objects from current state
        drug_conf = current_state['drug_conf']
        scoring_fn = current_state['scoring_fn']
        optimizer = current_state['optimizer']
        protein_coords = current_state['protein_coords_tensor']
        prot_vdw = current_state['prot_vdw']
        prot_is_hydrophobic = current_state['prot_is_hydrophobic']
        prot_is_hbond_donor = current_state['prot_is_hbond_donor']
        prot_is_hbond_acceptor = current_state['prot_is_hbond_acceptor']
        
        # Parse request to check if we're stopping
        data = request.json or {}
        is_stopped = data.get('is_stopped', False)
        
        # Initialize or get best scores/coords
        if current_state['global_best_scores'] is None:
            current_state['global_best_scores'] = float('inf') * torch.ones(drug_conf.n_conformers, device=drug_conf.device)
            current_state['round_best_scores'] = float('inf') * torch.ones(drug_conf.n_conformers, device=drug_conf.device)
        
        global_best_scores = current_state['global_best_scores']
        global_best_coords = current_state['global_best_coords']
        round_best_scores = current_state['round_best_scores']
        
        # Get current conformer coordinates - moved outside the if statement
        # so we have access to current ligand coordinates in all cases
        ligand_coords = drug_conf()
        
        # If not stopped, perform one optimization step
        if not is_stopped:
            optimizer.zero_grad()
            
            # Compute scores with current conformations
            scores = scoring_fn.compute_score(
                ligand_coords,
                protein_coords,
                drug_conf.ligand_vdw,
                prot_vdw,
                prot_is_hydrophobic,
                prot_is_hbond_donor,
                prot_is_hbond_acceptor,
                len(drug_conf.rotatable_bonds)
            )
            
            # Update global best scores/coordinates
            global_improved = scores < global_best_scores
            if global_improved.any():
                if global_best_coords is None:
                    global_best_coords = ligand_coords.clone()
                else:
                    global_best_coords[global_improved] = ligand_coords[global_improved].clone()
                global_best_scores[global_improved] = scores[global_improved]
                current_state['global_best_coords'] = global_best_coords
                
            # Check round improvements
            round_improved = scores < round_best_scores
            improvement_ratio = round_improved.float().mean()  # Calculate percentage of improved conformers
            significant_round_improvement = improvement_ratio > 0.05  # More than 5% improved
            
            if significant_round_improvement:
                current_state['last_improvement_step'] = current_state['step_count']
                round_best_scores[round_improved] = scores[round_improved]
            
            # Check for stagnation
            steps_since_improvement = current_state['step_count'] - current_state['last_improvement_step']
            if steps_since_improvement > current_state['stagnation_tolerance']:
                logger.info("Score stagnant - reinitializing coordinates")
                # Reinitialize coordinates
                drug_conf.conf_params.data = drug_conf._initialize_near_surface(0.0)
                # Reset round-based scores but keep global scores
                current_state['round_best_scores'] = float('inf') * torch.ones(drug_conf.n_conformers, device=drug_conf.device)
                round_best_scores = current_state['round_best_scores']
                # Reset last improvement step
                current_state['last_improvement_step'] = current_state['step_count']
            
            # Compute loss (negative because we want to minimize energy)
            loss = scores.mean()
            loss.backward()
            
            optimizer.step()
            
            # Update step count
            current_state['step_count'] += 1
        
        # Check if optimization is complete
        is_complete = current_state['step_count'] >= current_state['total_steps']
        
        # Generate conformer PDbs for visualization
        conformer_pdbs = []
        mol = current_state['ligand_mol']
        
        if is_stopped or is_complete:
            # If stopped or complete, show the single best global conformer
            sorted_indices = torch.argsort(global_best_scores)
            best_coords_np = global_best_coords[sorted_indices][:1].detach().cpu().numpy()
            scores_to_display = global_best_scores[sorted_indices][:1].detach().cpu()
        else:
            # MODIFIED: Show ALL current conformers (not just best ones)
            # Since we've reduced the number to 10, we can show all of them
            current_scores = scores  # Using scores computed earlier
            
            # Just include all current conformers instead of selecting top ones
            best_coords_np = ligand_coords.detach().cpu().numpy()
            scores_to_display = current_scores.detach().cpu()
        
        # Convert coordinates to PDB format
        for conf_coords in best_coords_np:
            mol_copy = Chem.Mol(mol)
            conf = Chem.Conformer(mol_copy.GetNumAtoms())
            
            for atom_idx in range(mol_copy.GetNumAtoms()):
                x, y, z = conf_coords[atom_idx]
                conf.SetAtomPosition(atom_idx, Chem.rdGeometry.Point3D(float(x), float(y), float(z)))
            
            mol_copy.RemoveAllConformers()
            mol_copy.AddConformer(conf)
            
            conformer_pdb = Chem.MolToPDBBlock(mol_copy)
            conformer_pdbs.append(conformer_pdb)
        
        # Clean up if complete
        if is_complete:
            current_state.update({
                'drug_conf': None,
                'scoring_fn': None,
                'optimizer': None,
                'protein_coords_tensor': None,
                'prot_vdw': None,
                'prot_is_hydrophobic': None,
                'prot_is_hbond_donor': None,
                'prot_is_hbond_acceptor': None,
                'global_best_scores': None,
                'global_best_coords': None,
                'round_best_scores': None,
                'step_count': 0,
                'last_improvement_step': 0
            })
        
        return jsonify({
            'message': 'Optimization step completed',
            'protein_pdb': current_state['pdb_data'],
            'conformer_pdbs': conformer_pdbs,
            'scores': scores_to_display.numpy().tolist(),
            'is_complete': is_complete,
            'is_stopped': is_stopped,
            'current_step': current_state['step_count'],
            'total_steps': current_state['total_steps']
        })
        
    except Exception as e:
        logger.error(f"Error in optimize_docking: {str(e)}")
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