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
        """Determine if atom is a hydrogen bond donor using implicit hydrogens."""
        # Oxygen or nitrogen with at least one hydrogen
        if atom.GetAtomicNum() in [7, 8]:  # N, O
            if atom.GetTotalNumHs() > 0:  # Check for implicit hydrogens
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

class ProteinAtomTyper:
    """Specialized atom typing for protein structures based on PDB conventions"""
    
    # Mappings of amino acid atom names to properties
    # Based on standard PDB naming conventions
    
    @staticmethod
    def get_atom_properties_from_pdb(pdb_data):
        """
        Determine atom properties directly from PDB file using amino acid knowledge
        Returns tuple of (atom_types, is_hydrophobic, is_hbond_donor, is_hbond_acceptor)
        """
        atom_types = []
        is_hydrophobic = []
        is_hbond_donor = []
        is_hbond_acceptor = []
        
        for line in pdb_data.split('\n'):
            if not line.startswith(('ATOM', 'HETATM')):
                continue
                
            try:
                # Extract element and atom name
                element = line[76:78].strip()
                if not element:  # If element field is empty, use first character of atom name
                    element = line[12:16].strip()[0]
                
                element = element.capitalize()
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                
                atom_types.append(element)
                
                # Determine properties based on atom identity and context
                # Hydrophobic atoms
                if element == 'C':
                    # Carbon is generally hydrophobic unless it's a carbonyl carbon
                    if atom_name in ['C', 'CA']:  # Backbone carbons
                        is_hydrophobic.append(False)
                    elif residue_name in ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TYR', 'TRP', 'PRO'] and atom_name.startswith('C'):
                        # Side chain carbons of hydrophobic amino acids
                        is_hydrophobic.append(True)
                    else:
                        is_hydrophobic.append(False)
                elif element in ['Cl', 'Br', 'I']:
                    is_hydrophobic.append(True)
                else:
                    is_hydrophobic.append(False)
                
                # Hydrogen bond donors
                if element == 'N':
                    # Amide/amine nitrogens are donors
                    if atom_name == 'N' or residue_name in ['LYS', 'ARG', 'HIS', 'TRP'] and atom_name.startswith('N'):
                        is_hbond_donor.append(True)
                    else:
                        is_hbond_donor.append(False)
                elif element == 'O':
                    # Hydroxyl oxygens are donors
                    if residue_name in ['SER', 'THR', 'TYR'] and atom_name.startswith('O') and not atom_name in ['O', 'OXT']:
                        is_hbond_donor.append(True)
                    else:
                        is_hbond_donor.append(False)
                elif element in ['Mg', 'Ca', 'Mn', 'Fe', 'Cu', 'Zn']:
                    # Metals can act as H-bond donors in Vina
                    is_hbond_donor.append(True)
                else:
                    is_hbond_donor.append(False)
                
                # Hydrogen bond acceptors
                if element == 'O':
                    # Almost all oxygens are acceptors
                    is_hbond_acceptor.append(True)
                elif element == 'N':
                    # Most nitrogens with free lone pairs are acceptors
                    if atom_name == 'N':  # Backbone N is usually not a good acceptor (has H)
                        is_hbond_acceptor.append(False)
                    elif residue_name in ['HIS', 'ASN', 'GLN', 'TRP']:
                        is_hbond_acceptor.append(True)
                    else:
                        is_hbond_acceptor.append(False)
                else:
                    is_hbond_acceptor.append(False)
                    
            except (ValueError, IndexError) as e:
                # Skip problematic lines
                logger.warning(f"Error processing PDB line: {e}")
                continue
        
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
        self.cutoff_distance = 1000.0
        
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
    def __init__(self, mol, protein_coords, n_conformers=100, sample_mode="grid", surface_distance=0.0, grid_spacing=2.0):
        super().__init__()
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        # self.device = 'cpu'
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
        
        # Calculate protein center and approximate radius for volume sampling
        protein_center = self.protein_coords.mean(dim=0)
        protein_radius = torch.max(torch.norm(self.protein_coords - protein_center, dim=1))
        self.protein_radius = protein_radius
        
        # Initialize parameters based on sample mode
        n_params = 6 + n_rotatable  # 3 for position, 3 for orientation, rest for rotatable bonds
        if sample_mode == "surface":
            initial_params = self._initialize_near_surface(surface_distance)
        elif sample_mode == "volume":
            initial_params = self._initialize_throughout_volume(protein_radius)
        elif sample_mode == "grid":
            initial_params = self._initialize_grid_sampling(grid_spacing)
        else:
            raise ValueError(f"Unknown sample mode: {sample_mode}")
        
        # Parameters: [n_conformers, n_params]
        self.conf_params = nn.Parameter(initial_params)
        
        self.to(self.device)
        
    def perturb_selected_conformations(self, indices_to_perturb):
        """
        Apply random perturbations to selected conformations.
        
        Parameters:
        -----------
        indices_to_perturb: torch.Tensor
            Tensor of indices of conformations to perturb
        """
        if indices_to_perturb.numel() == 0:
            return
            
        n_params = self.conf_params.shape[1]
        
        # Create perturbation: position (±2Å), orientation (±π), torsion angles (±π)
        position_perturbation = torch.randn(indices_to_perturb.size(0), 3, device=self.device) * 0.0
        orientation_perturbation = (torch.rand(indices_to_perturb.size(0), 3, device=self.device) * torch.pi) * (torch.randint(0, 2, (indices_to_perturb.size(0), 3), device=self.device) * 2 - 1)
        
        n_rotatable = n_params - 6
        if n_rotatable > 0:
            torsion_perturbation = (torch.rand(indices_to_perturb.size(0), n_rotatable, device=self.device) * torch.pi) * (torch.randint(0, 2, (indices_to_perturb.size(0), n_rotatable), device=self.device) * 2 - 1)
            perturbation = torch.cat([position_perturbation, orientation_perturbation, torsion_perturbation], dim=1)
        else:
            perturbation = torch.cat([position_perturbation, orientation_perturbation], dim=1)
        
        # Apply perturbation to all indices at once (GPU parallel)
        self.conf_params.data[indices_to_perturb] += perturbation
        
        # Ensure perturbed conformations stay within their grid cells using vectorized operations
        if hasattr(self, 'grid_cells'):
            # Create tensors of cell mins and maxs for all indices at once
            cell_mins = torch.zeros((indices_to_perturb.size(0), 3), device=self.device)
            cell_maxs = torch.zeros((indices_to_perturb.size(0), 3), device=self.device)
            
            # Extract cell boundaries (still need loop here, but just to build the tensors)
            for i, idx in enumerate(indices_to_perturb):
                cell_min, cell_max = self.grid_cells[idx]
                cell_mins[i] = cell_min
                cell_maxs[i] = cell_max
            
            # Apply clamp operation to all positions at once
            positions = self.conf_params.data[indices_to_perturb, :3]
            positions = torch.max(positions, cell_mins)
            positions = torch.min(positions, cell_maxs)
            self.conf_params.data[indices_to_perturb, :3] = positions

    def _initialize_throughout_volume(self, protein_radius):
        """Initialize conformers with positions sampled throughout the actual protein volume."""
        n_rotatable = len(self.rotatable_bonds)
        
        # Calculate bounding box of protein
        min_coords = self.protein_coords.min(dim=0)[0]
        max_coords = self.protein_coords.max(dim=0)[0]
        
        # Add a small buffer to the bounding box (20% of the range in each dimension)
        box_size = max_coords - min_coords
        buffer = box_size * 0.0
        min_coords = min_coords - buffer
        max_coords = max_coords + buffer
        
        # Generate random positions within the bounding box
        positions = min_coords + torch.rand(self.n_conformers, 3, device=self.device) * (max_coords - min_coords)
        
        # Random orientations (euler angles)
        orientations = torch.rand(self.n_conformers, 3, device=self.device) * 2 * np.pi
        
        # Use zero rotatable bond angles to maintain initial conformer rotations
        rot_angles = torch.zeros(self.n_conformers, n_rotatable, device=self.device)
        
        # Combine all parameters
        params = torch.cat([positions, orientations, rot_angles], dim=1)
        
        return params.to(torch.float32)

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

    def _initialize_grid_sampling(self, grid_spacing=2.0):
        """
        Initialize conformers on a grid throughout the protein volume.
        
        Parameters:
        -----------
        grid_spacing: float
            Distance between grid points in Angstroms
        """
        n_rotatable = len(self.rotatable_bonds)
        
        # Calculate bounding box of protein
        min_coords = self.protein_coords.min(dim=0)[0]
        max_coords = self.protein_coords.max(dim=0)[0]
        
        # Calculate number of grid points in each dimension
        box_size = max_coords - min_coords
        nx = max(1, int(box_size[0] / grid_spacing))
        ny = max(1, int(box_size[1] / grid_spacing))
        nz = max(1, int(box_size[2] / grid_spacing))
        
        logger.info(f"Creating grid with dimensions {nx}x{ny}x{nz} (total: {nx*ny*nz} points)")
        
        # Generate grid points
        self.all_grid_points = []
        for i in range(nx):
            x = min_coords[0] + (i + 0.5) * (box_size[0] / nx)
            for j in range(ny):
                y = min_coords[1] + (j + 0.5) * (box_size[1] / ny)
                for k in range(nz):
                    z = min_coords[2] + (k + 0.5) * (box_size[2] / nz)
                    point = torch.tensor([x, y, z], device=self.device)
                    
                    # Check if this grid point is near the protein
                    dists = torch.norm(self.protein_coords - point.unsqueeze(0), dim=1)
                    min_dist = torch.min(dists)
                    
                    # Only include points within a certain distance of the protein
                    if min_dist < 2 * grid_spacing:
                        self.all_grid_points.append([x.item(), y.item(), z.item()])
        
        # If no grid points are near the protein, fall back to random points
        if len(self.all_grid_points) == 0:
            logger.warning("No grid points near protein, falling back to random initialization")
            return self._initialize_throughout_volume(self.protein_radius)
        
        # Check if we need to adjust n_conformers or pad the grid points
        original_points = len(self.all_grid_points)
        logger.info(f"Generated {original_points} grid points near protein")
        
        if original_points < self.n_conformers:
            # If fewer grid points than conformers, adjust n_conformers to match grid points
            logger.info(f"Reducing n_conformers from {self.n_conformers} to {original_points} to match available grid points")
            self.n_conformers = original_points
            self.n_active_conformers = original_points
            # No padding needed
            points_needed = 0
        else:
            # Make sure the total number of grid points is a multiple of n_conformers
            points_needed = (self.n_conformers - (original_points % self.n_conformers)) % self.n_conformers
        
        if points_needed > 0:
            logger.info(f"Adding {points_needed} additional grid points to make total a multiple of {self.n_conformers}")
            # Add duplicates of first grid point with slight offsets
            for i in range(points_needed):
                base_point = self.all_grid_points[0]
                # Add small random offset (10% of grid spacing)
                offset = [(torch.rand(1).item() - 0.5) * 0.1 * grid_spacing for _ in range(3)]
                new_point = [base_point[0] + offset[0], 
                            base_point[1] + offset[1], 
                            base_point[2] + offset[2]]
                self.all_grid_points.append(new_point)
        
        self.total_points = len(self.all_grid_points)
        self.total_batches = self.total_points // self.n_conformers
        logger.info(f"Final grid has {self.total_points} points in {self.total_batches} batches of {self.n_conformers}")
        
        # Only use the first batch initially
        self.current_batch = 0
        start_idx = 0
        end_idx = min(self.n_conformers, self.total_points)
        grid_points = self.all_grid_points[start_idx:end_idx]
        self.n_active_conformers = len(grid_points)
        
        # Convert grid points to tensor
        positions = torch.tensor(grid_points, device=self.device, dtype=torch.float32)
        
        # Random orientations (euler angles)
        orientations = torch.rand(self.n_active_conformers, 3, device=self.device) * 2 * np.pi
        
        # Random rotatable bond angles
        rot_angles = torch.rand(self.n_active_conformers, n_rotatable, device=self.device) * 2 * np.pi
        
        # Store grid cell information for box constraints
        self.grid_cells = []
        self.grid_spacing = grid_spacing
        cell_factor = 1.2  # Make cells just 1.2x the grid spacing for local exploration
        for i in range(self.n_active_conformers):
            half_size = grid_spacing * cell_factor / 2.0
            center = torch.tensor(grid_points[i], device=self.device)
            cell_min = center - half_size
            cell_max = center + half_size
            self.grid_cells.append((cell_min, cell_max))
        
        # Combine all parameters
        params = torch.cat([positions, orientations, rot_angles], dim=1)
        
        logger.info(f"Initialized {self.n_active_conformers} conformers on grid")
        return params.to(torch.float32)


    def switch_to_next_batch(self):
        """
        Switch to the next batch of grid points.
        Returns True if switched to a new batch, False if no more batches.
        """
        # Move to next batch
        self.current_batch = (self.current_batch + 1) % self.total_batches
        
        # Calculate indices for this batch
        start_idx = self.current_batch * self.n_conformers
        end_idx = min(start_idx + self.n_conformers, self.total_points)
        
        # If we've processed all batches, return False
        if start_idx >= self.total_points:
            return False
        
        logger.info(f"Switching to batch {self.current_batch + 1}/{self.total_batches} with indices {start_idx}:{end_idx}")
        
        # Get grid points for this batch
        grid_points = self.all_grid_points[start_idx:end_idx]
        self.n_active_conformers = len(grid_points)
        
        # Get number of rotatable bonds
        n_rotatable = len(self.rotatable_bonds)
        
        # Create new parameter tensor
        positions = torch.tensor(grid_points, device=self.device, dtype=torch.float32)
        orientations = torch.rand(self.n_active_conformers, 3, device=self.device) * 2 * np.pi
        rot_angles = torch.rand(self.n_active_conformers, n_rotatable, device=self.device) * 2 * np.pi
        
        # Create new parameters tensor
        new_params = torch.cat([positions, orientations, rot_angles], dim=1)
        
        # Update conf_params - need to create a new Parameter since size changed
        old_params = self.conf_params
        self.conf_params = nn.Parameter(new_params)
        
        # Update grid cells
        self.grid_cells = []
        cell_factor = 1.2
        for i in range(self.n_active_conformers):
            half_size = self.grid_spacing * cell_factor / 2.0
            center = torch.tensor(grid_points[i], device=self.device)
            cell_min = center - half_size
            cell_max = center + half_size
            self.grid_cells.append((cell_min, cell_max))
        
        logger.info(f"Successfully switched to next batch with {self.n_active_conformers} conformers")
        return True

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

    def sample_random_conformations(self, protein_radius):
        """Generate new random conformations throughout the protein volume."""
        # Store the original requires_grad state
        requires_grad = self.conf_params.requires_grad
        
        # Temporarily set requires_grad to False to avoid gradient computation
        self.conf_params.requires_grad = False
        
        # Generate new random parameters
        self.conf_params.data = self._initialize_throughout_volume(protein_radius)
        
        # Compute coordinates
        coords = self()
        
        # Restore the original requires_grad state
        self.conf_params.requires_grad = requires_grad
        
        return coords

    def forward(self) -> torch.Tensor:
        """Forward pass returns the computed 3D coordinates for all conformers."""
        return self._compute_atom_coordinates()

app = Flask(__name__)
CORS(app)

# Store current state
current_state = {
    # Input data
    'pdb_data': None,
    'drug_data': None,
    'smile_sequence': None,
    'ligand_mol': None,
    'protein_coords': None,
    'protein_types': None,
    
    # Docking components
    'drug_conf': None,
    'scoring_fn': None,
    'optimizer': None,
    'protein_coords_tensor': None,
    'prot_vdw': None,
    'prot_is_hydrophobic': None,
    'prot_is_hbond_donor': None,
    'prot_is_hbond_acceptor': None,
    
    # Optimization state
    'step_count': 0,
    'total_steps': 1e10,
    'global_best_scores': None,  # Best scores across all iterations
    'global_best_coords': None,  # Best coordinates across all iterations
    
    # ILS-specific tracking variables
    'local_best_scores': None,   # Best scores since last perturbation
    'no_improvement_steps': None,  # Number of steps without significant improvement for each conformation
    'pre_perturbation_params': None, # Parameters before perturbation
    'pre_perturbation_scores': None, # Scores before perturbation
    'max_no_improvement': 100000,       # Threshold for perturbation
    'improvement_threshold': 0.01, # Require 1% improvement to be considered significant 
    'protein_radius': None
}

@app.route('/api/load_protein', methods=['POST'])
def load_protein():
    try:
        data = request.json
        pdb_data = data.get('pdb_data')
        
        if not pdb_data:
            return jsonify({'error': 'No PDB data provided'}), 400
        
        # Extract coordinates and atom types directly
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
                    
                    # Simply capitalize the element
                    atom_types.append(element.capitalize())
                    
                except (ValueError, IndexError):
                    continue
        
        if len(coords) == 0:
            return jsonify({'error': 'Failed to extract protein coordinates'}), 400
        
        # Store the data
        current_state['pdb_data'] = pdb_data
        current_state['protein_coords'] = np.array(coords)
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
        
        # Add hydrogens temporarily for better 3D geometry optimization
        mol_with_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_with_h, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol_with_h)
        mol = Chem.RemoveHs(mol_with_h)
        
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
        # device = 'cpu'
        logger.info(f"Using device: {device}")
        
        # Calculate protein center and radius for volume sampling
        protein_coords_tensor = torch.tensor(current_state['protein_coords'], device=device, dtype=torch.float32)
        protein_center = protein_coords_tensor.mean(dim=0)
        protein_radius = torch.max(torch.norm(protein_coords_tensor - protein_center, dim=1))
        
        # Get n_conformers from request if provided, otherwise use default
        n_conformers = 4000  # Set your desired number of conformers here
        grid_spacing = 3.0  # Set your desired grid spacing here
        
        # Reset all state variables
        current_state.update({
            # Docking components
            'drug_conf': None,
            'scoring_fn': None,
            'optimizer': None,
            'protein_coords_tensor': None,
            'prot_vdw': None,
            'prot_is_hydrophobic': None,
            'prot_is_hbond_donor': None,
            'prot_is_hbond_acceptor': None,
            
            # Optimization state
            'step_count': 0,
            'batch_step_count': 0,
            'steps_per_batch': 10000,  # Process each batch for 50 steps before switching
            'total_steps': 1e10,
            'global_best_scores': None,
            'global_best_coords': None,
            'global_best_params': None,
            
            # Batch processing state
            'current_batch': 0,
            'all_batches_best_scores': {},  # Store best scores for each batch
            'all_batches_best_coords': {},  # Store best coordinates for each batch
            'all_batches_best_params': {},  # Store best parameters for each batch
            
            # ILS-specific variables
            'local_best_scores': None,
            'no_improvement_steps': None,
            'pre_perturbation_params': None,
            'pre_perturbation_scores': None,
            'max_no_improvement': 10,
            'improvement_threshold': 0.05,
            'protein_radius': protein_radius.item()
        })
        
        # Initialize drug conformation model with grid sampling
        drug_conf = DrugConformation(
            mol=current_state['ligand_mol'],
            protein_coords=current_state['protein_coords'],
            n_conformers=n_conformers,
            sample_mode="grid",
            grid_spacing=grid_spacing
        ).to(device)
        
        # Get the actual number of active conformers in the current batch
        n_active_conformers = drug_conf.n_active_conformers
        logger.info(f"Active conformers in current batch: {n_active_conformers}")
        
        # Initialize scoring function
        scoring_fn = VinaScoringFunction(device, current_state['ligand_mol'])
        
        # Set up optimizer using SGD
        optimizer = torch.optim.SGD([drug_conf.conf_params],
                                lr=0.05,
                                momentum=0.9,
                                dampening=0,
                                weight_decay=0,
                                nesterov=True)
        
        # Use specialized PDB-based atom typing for proteins
        try:
            logger.info("Using specialized protein atom typing from PDB data")
            pdb_data = current_state['pdb_data']
            protein_types, is_hydrophobic, is_hbond_donor, is_hbond_acceptor = ProteinAtomTyper.get_atom_properties_from_pdb(pdb_data)
            
            prot_vdw = torch.tensor(
                [AtomTyper.VDW_RADII.get(t, 1.7) for t in protein_types],
                device=device,
                dtype=torch.float32
            )
        except Exception as e:
            logger.warning(f"Error in specialized protein typing: {str(e)}. Falling back to basic atom typing")
            prot_types = current_state['protein_types']
            is_hydrophobic = [t in ['C', 'Cl', 'Br', 'I'] for t in prot_types]
            is_hbond_donor = [t in ['N', 'O', 'Mg', 'Ca', 'Mn', 'Fe', 'Cu', 'Zn'] for t in prot_types]
            is_hbond_acceptor = [t in ['N', 'O'] for t in prot_types]
            
            prot_vdw = torch.tensor(
                [AtomTyper.VDW_RADII.get(t, 1.7) for t in prot_types],
                device=device,
                dtype=torch.float32
            )
        
        # Convert properties to tensors
        prot_is_hydrophobic = torch.tensor(is_hydrophobic, device=device, dtype=torch.bool)
        prot_is_hbond_donor = torch.tensor(is_hbond_donor, device=device, dtype=torch.bool)
        prot_is_hbond_acceptor = torch.tensor(is_hbond_acceptor, device=device, dtype=torch.bool)

        # Initialize best scores and coordinates trackers with the actual number of conformers
        global_best_scores = float('inf') * torch.ones(n_active_conformers, device=device)
        
        # Now initialize the ILS tracking variables with the actual number of conformers
        ils_tracking = {
            'local_best_scores': float('inf') * torch.ones(n_active_conformers, device=device),
            'no_improvement_steps': torch.zeros(n_active_conformers, device=device, dtype=torch.int),
            'pre_perturbation_params': None,
            'pre_perturbation_scores': None
        }
        
        # Store all necessary state
        current_state.update({
            # Docking components
            'drug_conf': drug_conf,
            'scoring_fn': scoring_fn,
            'optimizer': optimizer,
            'protein_coords_tensor': protein_coords_tensor,
            'prot_vdw': prot_vdw,
            'prot_is_hydrophobic': prot_is_hydrophobic,
            'prot_is_hbond_donor': prot_is_hbond_donor,
            'prot_is_hbond_acceptor': prot_is_hbond_acceptor,
            
            # Optimization state
            'global_best_scores': global_best_scores,
            'global_best_coords': None,
            
            # ILS tracking variables
            **ils_tracking,
            
            # Add total batches info
            'total_batches': drug_conf.total_batches,
            'total_points': drug_conf.total_points
        })
        
        return jsonify({
            'message': 'Docking initialized successfully',
            'total_grid_points': drug_conf.total_points,
            'total_batches': drug_conf.total_batches,
            'current_batch': drug_conf.current_batch,
            'active_conformers': drug_conf.n_active_conformers
        })
        
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
        
        # Initialize or get ILS tracking variables
        no_improvement_steps = current_state['no_improvement_steps']
        global_best_scores = current_state['global_best_scores']
        global_best_coords = current_state['global_best_coords']
        local_best_scores = current_state['local_best_scores']
        pre_perturbation_params = current_state['pre_perturbation_params']
        pre_perturbation_scores = current_state['pre_perturbation_scores']
        improvement_threshold = current_state['improvement_threshold']
        
        # Get batch processing state
        batch_step_count = current_state['batch_step_count']
        steps_per_batch = current_state['steps_per_batch']
        current_batch = current_state['current_batch']
        all_batches_best_scores = current_state['all_batches_best_scores']
        all_batches_best_coords = current_state['all_batches_best_coords']
        all_batches_best_params = current_state['all_batches_best_params']
        
        # Variables to track current iteration conformations
        current_ligand_coords = None
        current_scores = None
        
        # Check if it's time to switch to next batch
        should_switch_batch = False
        
        # If not stopped, perform 10 optimization steps
        if not is_stopped:
            needs_perturbation = None
            significant_improvement = None
            
            # Run 10 optimization steps in a single API call
            for step_index in range(10):
                # Check if we should switch to next batch
                if batch_step_count >= steps_per_batch and drug_conf.total_batches > 1:
                    logger.info(f"Completed {batch_step_count} steps for batch {current_batch + 1}, switching to next batch")
                    
                    # Store best results for current batch
                    all_batches_best_scores[current_batch] = global_best_scores.clone().cpu()
                    if global_best_coords is not None:
                        all_batches_best_coords[current_batch] = global_best_coords.clone().cpu()
                    
                    # Switch to next batch
                    success = drug_conf.switch_to_next_batch()
                    if success:
                        # Update current batch number
                        current_batch = drug_conf.current_batch
                        current_state['current_batch'] = current_batch
                        
                        # Create new optimizer for new parameters
                        optimizer = torch.optim.SGD([drug_conf.conf_params],
                                                lr=0.005,
                                                momentum=0.9,
                                                dampening=0,
                                                weight_decay=0,
                                                nesterov=True)
                        current_state['optimizer'] = optimizer
                        
                        # Reset ILS tracking variables for new batch
                        n_active_conformers = drug_conf.n_active_conformers
                        no_improvement_steps = torch.zeros(n_active_conformers, device=drug_conf.device, dtype=torch.int)
                        local_best_scores = float('inf') * torch.ones(n_active_conformers, device=drug_conf.device)
                        global_best_scores = float('inf') * torch.ones(n_active_conformers, device=drug_conf.device)
                        pre_perturbation_params = None
                        pre_perturbation_scores = None
                        global_best_coords = None
                        
                        # Update state
                        current_state.update({
                            'no_improvement_steps': no_improvement_steps,
                            'local_best_scores': local_best_scores,
                            'global_best_scores': global_best_scores,
                            'global_best_coords': global_best_coords,
                            'pre_perturbation_params': pre_perturbation_params,
                            'pre_perturbation_scores': pre_perturbation_scores,
                            'batch_step_count': 0
                        })
                        
                        # If this batch has been processed before, restore best scores/coords
                        if current_batch in all_batches_best_scores:
                            best_scores_cpu = all_batches_best_scores[current_batch]
                            global_best_scores = best_scores_cpu.to(drug_conf.device)
                            current_state['global_best_scores'] = global_best_scores
                            
                            if current_batch in all_batches_best_coords:
                                best_coords_cpu = all_batches_best_coords[current_batch]
                                global_best_coords = best_coords_cpu.to(drug_conf.device)
                                current_state['global_best_coords'] = global_best_coords
                        
                        should_switch_batch = True
                        # Compute new batch coordinates and scores after switching
                        ligand_coords = drug_conf()
                        current_ligand_coords = ligand_coords
                        
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
                        current_scores = scores
                        
                        # Break the loop to avoid optimization on this step
                        break
                    else:
                        logger.info("No more batches to process, continuing with current batch")
                
                # Skip optimization step if we just switched batches
                if not should_switch_batch:
                    optimizer.zero_grad()
                    
                    # Get current conformer coordinates
                    ligand_coords = drug_conf()
                    current_ligand_coords = ligand_coords  # Store current conformations
                    
                    # Compute scores
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
                    current_scores = scores  # Store current scores
                    
                    # Compute loss (negative because we want to minimize energy)
                    loss = scores.mean()
                    loss.backward()
                    
                    optimizer.step()
                    
                    # Apply box constraints to keep each conformation within its grid cell
                    if hasattr(drug_conf, 'grid_cells'):
                        # Get current positions after optimization step
                        current_positions = drug_conf.conf_params.data[:, :3]
                        
                        # Apply constraints for each conformation
                        for i in range(drug_conf.n_active_conformers):
                            # Get grid cell boundaries for this conformer
                            cell_min, cell_max = drug_conf.grid_cells[i]
                            
                            # Check if position is outside cell
                            position = current_positions[i]
                            outside_cell = (position < cell_min) | (position > cell_max)
                            
                            if outside_cell.any():
                                # Project back into cell
                                position = torch.max(position, cell_min)
                                position = torch.min(position, cell_max)
                                drug_conf.conf_params.data[i, :3] = position
                    
                    # Initialize pre_perturbation variables if they're None
                    if pre_perturbation_params is None:
                        pre_perturbation_params = drug_conf.conf_params.clone()
                        current_state['pre_perturbation_params'] = pre_perturbation_params
                    
                    if pre_perturbation_scores is None:
                        pre_perturbation_scores = scores.clone()
                        current_state['pre_perturbation_scores'] = pre_perturbation_scores
                    
                    # Update global best scores/coordinates
                    global_improved = scores < global_best_scores
                    if global_improved.any():
                        if global_best_coords is None:
                            global_best_coords = ligand_coords.clone()
                        else:
                            global_best_coords[global_improved] = ligand_coords[global_improved].clone()
                        global_best_scores[global_improved] = scores[global_improved]
                        current_state['global_best_coords'] = global_best_coords
                        
                        # Update pre_perturbation_params to match global best
                        pre_perturbation_params[global_improved] = drug_conf.conf_params[global_improved].clone()
                        pre_perturbation_scores[global_improved] = global_best_scores[global_improved].clone()
                    
                    # Check for significant improvement against local best scores
                    improvement_percentage = (scores - local_best_scores) 
                    
                    # Significant improvement if current score is better than local best by at least the threshold
                    significant_improvement = improvement_percentage < -improvement_threshold
                    
                    # Update local best scores
                    local_improved = scores < local_best_scores
                    if local_improved.any():
                        local_best_scores[local_improved] = scores[local_improved]
                    
                    # Update counters for each conformation individually
                    # Reset counters for significantly improved conformations
                    no_improvement_steps[significant_improvement] = 0
                    
                    # Increment counters for non-improving conformations
                    no_improvement_steps[~significant_improvement] += 1
                    
                    # ILS Step 1: Check if any conformations need perturbation (reached max_no_improvement)
                    needs_perturbation = no_improvement_steps >= current_state['max_no_improvement']

                    if needs_perturbation.any():
                        logger.info(f"Perturbing {needs_perturbation.sum().item()} conformations")
                        perturb_indices = torch.where(needs_perturbation)[0]
                        
                        # Create masks for different conditions using tensor operations
                        improved_after_perturbation = scores[perturb_indices] < pre_perturbation_scores[perturb_indices]
                        global_improvement = scores[perturb_indices] < global_best_scores[perturb_indices]
                        
                        # Update global best for indices that improved globally - all in parallel
                        if global_improvement.any():
                            improved_indices = perturb_indices[global_improvement]
                            global_best_scores[improved_indices] = scores[improved_indices]
                            if global_best_coords is not None:
                                global_best_coords[improved_indices] = ligand_coords[improved_indices].clone()
                            pre_perturbation_params[improved_indices] = drug_conf.conf_params[improved_indices].clone()
                            pre_perturbation_scores[improved_indices] = scores[improved_indices].clone()
                        
                        # Apply perturbations to all conformations that need it
                        drug_conf.perturb_selected_conformations(perturb_indices)
                        
                        # Reset no_improvement_steps and local_best_scores for perturbed conformations
                        no_improvement_steps[perturb_indices] = 0
                        local_best_scores[perturb_indices] = float('inf')
                    
                    # Update batch step count
                    current_state['batch_step_count'] += 1
                
                # Update total step count
                current_state['step_count'] += 1
                
                # Save the updated ILS state variables
                current_state['no_improvement_steps'] = no_improvement_steps
                current_state['local_best_scores'] = local_best_scores
                current_state['pre_perturbation_params'] = pre_perturbation_params
                current_state['pre_perturbation_scores'] = pre_perturbation_scores
                
                # Update batch_step_count for the condition at the beginning of the loop
                batch_step_count = current_state['batch_step_count']
                
                # Check if we're at the maximum number of steps
                if current_state['step_count'] >= current_state['total_steps']:
                    # Stop running more steps if we've reached the maximum
                    break
        
        # Check if optimization is complete
        is_complete = current_state['step_count'] >= current_state['total_steps']
        
        # If we just switched batches, we need current_scores and current_ligand_coords
        if should_switch_batch and current_ligand_coords is None:
            # Compute coordinates and scores for the new batch
            ligand_coords = drug_conf()
            current_ligand_coords = ligand_coords
            
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
            current_scores = scores
        
        # ===== VISUALIZATION CONFIGURATION =====
        VIEW_CURRENT_CONFORMATIONS = False
        VIEW_ALL_CONFORMERS = False
        
        if VIEW_CURRENT_CONFORMATIONS and current_ligand_coords is not None:
            conformer_coords = current_ligand_coords
            conformer_scores = current_scores
            # logger.info("Using current iteration conformations")
        else:
            conformer_coords = global_best_coords
            conformer_scores = global_best_scores
            # logger.info("Using best conformations across all iterations")
            
        if (not VIEW_ALL_CONFORMERS or is_stopped) and conformer_coords is not None:
            sorted_indices = torch.argsort(conformer_scores)
            conformer_coords = conformer_coords[sorted_indices]
            conformer_scores = conformer_scores[sorted_indices]
            # logger.info("Conformations sorted by score")
        else:
            # logger.info("Keeping original order of conformations (unsorted)")
            pass
        
        if is_stopped:
            # Return best conformer across all batches
            all_best_scores = []
            all_best_coords = []
            
            # Get best scores from current batch
            if global_best_coords is not None:
                min_idx = torch.argmin(global_best_scores)
                all_best_scores.append(global_best_scores[min_idx].item())
                all_best_coords.append(global_best_coords[min_idx].clone())
            
            # Get best scores from all previously processed batches
            for batch_idx in all_batches_best_scores:
                batch_scores = all_batches_best_scores[batch_idx]
                if batch_idx in all_batches_best_coords:
                    batch_coords = all_batches_best_coords[batch_idx]
                    min_idx = torch.argmin(batch_scores)
                    all_best_scores.append(batch_scores[min_idx].item())
                    all_best_coords.append(batch_coords[min_idx].clone())
            
            if all_best_scores:
                # Find the best score across all batches
                best_idx = np.argmin(all_best_scores)
                best_score = all_best_scores[best_idx]
                best_coords = all_best_coords[best_idx]
                
                # Use the single best conformer across all batches
                conformer_coords = best_coords.unsqueeze(0).to(drug_conf.device)
                conformer_scores = torch.tensor([best_score], device=drug_conf.device)
                n_top = 1
            else:
                n_top = 0
                
            logger.info("Returning overall best conformation across all batches")
        elif is_complete:
            if not VIEW_ALL_CONFORMERS:
                n_top = 1
            else:
                n_top = len(conformer_coords) if conformer_coords is not None else 0
        else:
            if not VIEW_ALL_CONFORMERS:
                n_top = min(10, len(conformer_coords) if conformer_coords is not None else 0)
            else:
                n_top = len(conformer_coords) if conformer_coords is not None else 0
        
        conformer_pdbs = []
        if conformer_coords is not None and n_top > 0:
            mol = current_state['ligand_mol']
            coords_np = conformer_coords[:n_top].detach().cpu().numpy()
            
            for conf_coords in coords_np:
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
                'step_count': 0,
                'batch_step_count': 0,
                'global_best_scores': None,
                'global_best_coords': None,
                'local_best_scores': None,
                'no_improvement_steps': None,
                'pre_perturbation_params': None,
                'pre_perturbation_scores': None,
                'all_batches_best_scores': {},
                'all_batches_best_coords': {}
            })
        
        # Get ILS statistics for reporting
        ils_stats = {
            'perturbed_conformers': needs_perturbation.sum().item() if needs_perturbation is not None else 0,
            'max_no_improvement': current_state['max_no_improvement'],
            'improvement_threshold': current_state['improvement_threshold'],
            'significantly_improved_conformers': significant_improvement.sum().item() if significant_improvement is not None else 0
        }
        
        # Batch processing stats
        batch_stats = {
            'current_batch': current_batch + 1,  # 1-indexed for display
            'total_batches': drug_conf.total_batches,
            'batch_step_count': current_state['batch_step_count'],
            'steps_per_batch': steps_per_batch,
            'total_grid_points': drug_conf.total_points,
            'active_conformers': drug_conf.n_active_conformers,
            'processed_batches': len(all_batches_best_scores)
        }
        
        return jsonify({
            'message': 'Optimization steps completed',
            'steps_performed': 10,
            'protein_pdb': current_state['pdb_data'],
            'conformer_pdbs': conformer_pdbs,
            'scores': conformer_scores[:n_top].detach().cpu().numpy().tolist() if n_top > 0 else [],
            'is_complete': is_complete,
            'is_stopped': is_stopped,
            'current_step': current_state['step_count'],
            'total_steps': current_state['total_steps'],
            'view_current_conformations': VIEW_CURRENT_CONFORMATIONS,
            'view_all_conformers': VIEW_ALL_CONFORMERS,
            'num_conformers_displayed': n_top,
            'total_conformers_available': len(conformer_coords) if conformer_coords is not None else 0,
            'ils_stats': ils_stats,
            'batch_stats': batch_stats,
            'batch_switched': should_switch_batch
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