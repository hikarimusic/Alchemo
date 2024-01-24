from flask import Flask, jsonify, request
from flask_cors import CORS
import random

app = Flask(__name__)
CORS(app)

def parse_pdb_data(pdb_data):
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
            
            base_structure.append([record_type, atom_serial_number, atom_name, alt_loc, residue_name, chain_id, 
                                   residue_seq_number, insertion_code, x, y, z, occupancy, temp_factor, 
                                   segment_id, element_symbol, charge])

    return base_structure

def generate_frames(base_structure, num_frames, displacement=0.1):
    frames = []
    for frame in range(num_frames):
        frame_data = "HEADER    TEST PDB\n"
        for atom in base_structure:
            # Displacing only the x, y, z coordinates
            displaced_x = atom[8] + random.uniform(-displacement, displacement)
            displaced_y = atom[9] + random.uniform(-displacement, displacement)
            displaced_z = atom[10] + random.uniform(-displacement, displacement)

            frame_data += "ATOM  {:>5} {:<4}{:1}{:>3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}      {:<4}{:>2}{:2}\n".format(
                atom[1], atom[2], atom[3], atom[4], atom[5], atom[6], atom[7], displaced_x, displaced_y, 
                displaced_z, atom[11], atom[12], atom[13], atom[14], atom[15]
            )
        frame_data += "TER\nEND\n"
        frames.append(frame_data)
        
    return frames

@app.route('/generateProteinStructure', methods=['POST'])
def generate_protein_structure():
    data = request.json
    protein_sequence = data['sequence']
    print("HELLO", protein_sequence)
    # Code to generate PDB data from protein sequence
    pdb_data = "generated PDB data for protein"
    return jsonify(pdbData=pdb_data)

@app.route('/generateMoleculeStructure', methods=['POST'])
def generate_molecule_structure():
    data = request.json
    molecule_sequence = data['sequence']
    # Code to generate PDB data from molecule sequence
    pdb_data = "generated PDB data for molecule"
    return jsonify(pdbData=pdb_data)

@app.route('/animateProtein', methods=['POST'])
def animate_protein():
    content = request.json
    pdb_data = content['pdbData']
    base_structure = parse_pdb_data(pdb_data)
    frames = generate_frames(base_structure, 10)  # Generate 10 frames
    return jsonify({'frames': frames})

if __name__ == '__main__':
    app.run(debug=True)
