import numpy as np
from scipy.spatial.transform import Rotation as R
import re

def transform_points(points, point_a_index, plane_points_indices, side_chain_indices=[], side_remove_indices=[]):
    # Step 1: Rotate Points B, C, and D to the same plane (z=0)
    p_a, p_b, p_c, p_d = points[[point_a_index] + plane_points_indices]
    normal_vector = np.cross(p_c - p_b, p_d - p_b)
    normal_vector /= np.linalg.norm(normal_vector)  # Normalize the vector

    z_axis = np.array([0, 0, 1])
    rotation_axis = np.cross(normal_vector, z_axis)
    rotation_angle = np.arccos(np.clip(np.dot(normal_vector, z_axis), -1.0, 1.0))
    rotation_matrix = R.from_rotvec(rotation_angle * rotation_axis / np.linalg.norm(rotation_axis)).as_matrix()
    rotated_points = np.dot(points, rotation_matrix.T)
    average_z = np.mean([rotated_points[i][2] for i in plane_points_indices])
    rotated_points[:, 2] -= average_z
    # print("rotate plane")
    # print(rotated_points)

    # Step 2: Rotate around Z-axis to align x-coordinates of C and D
    _, p_c, p_d = rotated_points[plane_points_indices]  # Re-fetch the transformed points C and D
    angle_to_x_axis = np.arctan2(p_d[0] - p_c[0], p_d[1] - p_c[1])
    final_rotation = R.from_euler('z', angle_to_x_axis).as_matrix()
    rotated_points = np.dot(rotated_points, final_rotation.T)
    # print("rotate C=O")
    # print(rotated_points)

    # Step 3: Rotate Point A around the axis between B and C
    p_a, p_b, p_c = rotated_points[[point_a_index] + plane_points_indices[:2]]
    axis = p_c - p_b
    axis /= np.linalg.norm(axis)  # Normalize the axis
    projection = p_b + np.dot(p_a - p_b, axis) * axis
    angle = np.arcsin(p_a[2]/np.linalg.norm(p_a - projection)) + np.pi
    if np.cross(p_a - p_b, p_c - p_b)[2] < 0:
        angle = np.pi - angle
    rotation_matrix_a = R.from_rotvec(angle * axis).as_matrix()
    rotated_points[point_a_index] = np.dot(rotation_matrix_a, (p_a - projection).T).T + projection
    nh_index = [i for i in range(4, points.shape[0]-1)]
    p_nh = rotated_points[nh_index]
    rotated_points[nh_index] = np.dot(rotation_matrix_a, (p_nh - projection).T).T + projection
    # print("rotate N1")
    # print(rotated_points)

    # Step 4: Rotate last point
    p_a, p_b, p_c = rotated_points[[-1] + plane_points_indices[:2]]
    axis = p_c - p_b
    axis /= np.linalg.norm(axis)  # Normalize the axis
    projection = p_b + np.dot(p_a - p_b, axis) * axis
    angle = np.arcsin(p_a[2]/np.linalg.norm(p_a - projection))
    if np.cross(p_a - p_b, p_c - p_b)[2] < 0:
        angle = np.pi - angle
    rotation_matrix_a = R.from_rotvec(angle * axis).as_matrix()
    rotated_points[-1] = np.dot(rotation_matrix_a, (p_a - projection).T).T + projection
    # print("rotate N2")
    # print(rotated_points)

    # Step 5: Translate Point A to the origin
    translation = -rotated_points[point_a_index]
    rotated_points = rotated_points + translation
    # print("Translate")
    # print(rotated_points)

    # Step 6: Inverse
    _, p_c, p_d = rotated_points[plane_points_indices]  # Re-fetch the transformed points C and D
    if p_d[1] > p_c[1]:  # Invert y and z coordinates if y-coordinate of D is greater than C
        inversion_matrix = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
        rotated_points = np.dot(rotated_points, inversion_matrix)
    print("Inverse")
    print(rotated_points)

    # Step 7: Rotate Side Chain
    if (len(side_chain_indices)==0):
        return rotated_points
    p_c, p_b, p_a = rotated_points[[1] + side_chain_indices]
    axis = p_c - p_b
    axis /= np.linalg.norm(axis)  # Normalize the axis
    projection = p_b + np.dot(p_a - p_b, axis) * axis
    # angle = np.arcsin((p_c[0]-p_b[0])/np.linalg.norm(p_a - projection)) + np.pi
    normal = np.cross(p_a - p_b, p_c - p_b)
    normal /= np.linalg.norm(normal)
    angle = np.arccos(np.dot(normal, np.array([1., 0., 0.]))) + np.pi
    if np.cross(p_a - p_b, p_c - p_b)[1] < 0:
        angle = - angle
    rotation_matrix_a = R.from_rotvec(angle * axis).as_matrix()
    move_idx = [i for i in range(4, points.shape[0]-1) if i not in side_remove_indices]
    rotated_points[move_idx] = np.dot(rotation_matrix_a, (rotated_points[move_idx] - projection).T).T + projection
    print("rotate side")
    print(rotated_points)

    return rotated_points

def modify_pdb_data(pdb_data, transformed_points):
    lines = pdb_data.split('\n')
    modified_structure = []
    atom_serial_number = 1

    for i, line in enumerate(lines):
        if line.startswith("ATOM"):
            # Extract basic information
            record_type = line[0:6].strip()
            residue_name = line[17:20].strip()
            chain_id = line[21].strip()
            insertion_code = line[26].strip()
            occupancy = float(line[54:60].strip())
            temp_factor = float(line[60:66].strip())
            segment_id = line[72:76].strip()
            element_symbol = line[76:78].strip()
            charge = line[78:80].strip()

            # Get transformed coordinates
            x, y, z = transformed_points[i-1]

            # Change residue sequence number
            residue_seq_number = 1 if i < len(lines) - 2 else 2

            # Construct the modified line
            modified_line = f"ATOM  {atom_serial_number:>5} {line[12:16]:<4}{line[16]:1}{residue_name:>3} {chain_id:1}{residue_seq_number:>4}{insertion_code:1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{temp_factor:>6.2f}      {segment_id:<4}{element_symbol:>2}{charge:2}"
            modified_structure.append(modified_line)

            atom_serial_number += 1

    return '\n'.join(modified_structure)


def print_formatted(array):
    for row in array:
        formatted_row = ' '.join(f'{x:7.3f}' for x in row)
        print(formatted_row)


res = ""


name = ">A"
data = """
ATOM    268  N   ALA A  18       3.147   2.159 -20.048  1.00  0.00           N  
ATOM    269  CA  ALA A  18       3.133   3.490 -19.458  1.00  0.00           C  
ATOM    270  C   ALA A  18       2.871   4.551 -20.522  1.00  0.00           C  
ATOM    271  O   ALA A  18       3.741   4.841 -21.345  1.00  0.00           O  
ATOM    272  CB  ALA A  18       4.442   3.766 -18.735  1.00  0.00           C  
ATOM    273  H   ALA A  18       3.990   1.662 -20.097  1.00  0.00           H  
ATOM    274  HA  ALA A  18       2.334   3.521 -18.732  1.00  0.00           H  
ATOM    275  HB1 ALA A  18       4.380   4.718 -18.230  1.00  0.00           H  
ATOM    276  HB2 ALA A  18       5.251   3.790 -19.452  1.00  0.00           H  
ATOM    277  HB3 ALA A  18       4.625   2.986 -18.011  1.00  0.00           H  
ATOM    278  N   THR A  19       1.661   5.108 -20.510  1.00  0.00           N
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 6], side_remove_indices=[5, 6])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">C"
data = """
ATOM    114  N   CYS A   7       4.071  -5.657  -5.445  1.00  0.00           N  
ATOM    115  CA  CYS A   7       4.063  -4.732  -6.572  1.00  0.00           C  
ATOM    116  C   CYS A   7       5.473  -4.337  -7.012  1.00  0.00           C  
ATOM    117  O   CYS A   7       5.636  -3.630  -8.001  1.00  0.00           O  
ATOM    118  CB  CYS A   7       3.259  -3.485  -6.205  1.00  0.00           C  
ATOM    119  SG  CYS A   7       2.965  -2.354  -7.583  1.00  0.00           S  
ATOM    120  H   CYS A   7       4.106  -5.302  -4.531  1.00  0.00           H  
ATOM    121  HA  CYS A   7       3.572  -5.228  -7.395  1.00  0.00           H  
ATOM    122  HB2 CYS A   7       2.296  -3.788  -5.823  1.00  0.00           H  
ATOM    123  HB3 CYS A   7       3.788  -2.938  -5.436  1.00  0.00           H  
ATOM    124  HG  CYS A   7       4.063  -2.337  -8.332  1.00  0.00           H  
ATOM    125  N   LEU A   8       6.490  -4.793  -6.296  1.00  0.00           N  
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[6, 7])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">D"
data = """
ATOM     83  N   ASP A   5       5.017  -7.221  -0.172  1.00  0.00           N  
ATOM     84  CA  ASP A   5       4.956  -6.177  -1.197  1.00  0.00           C  
ATOM     85  C   ASP A   5       5.039  -6.777  -2.596  1.00  0.00           C  
ATOM     86  O   ASP A   5       6.119  -6.931  -3.161  1.00  0.00           O  
ATOM     87  CB  ASP A   5       6.074  -5.147  -0.998  1.00  0.00           C  
ATOM     88  CG  ASP A   5       5.965  -4.421   0.328  1.00  0.00           C  
ATOM     89  OD1 ASP A   5       4.938  -3.744   0.557  1.00  0.00           O  
ATOM     90  OD2 ASP A   5       6.902  -4.528   1.149  1.00  0.00           O  
ATOM     91  H   ASP A   5       5.829  -7.307   0.371  1.00  0.00           H  
ATOM     92  HA  ASP A   5       4.004  -5.677  -1.096  1.00  0.00           H  
ATOM     93  HB2 ASP A   5       7.029  -5.649  -1.035  1.00  0.00           H  
ATOM     94  HB3 ASP A   5       6.025  -4.416  -1.792  1.00  0.00           H  
ATOM     95  N   LEU A   6       3.880  -7.102  -3.149  1.00  0.00           N  
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[8, 9])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">E"
data = """
ATOM    166  N   GLU A  11      10.194  -8.290  -6.211  1.00  0.00           N  
ATOM    167  CA  GLU A  11       9.341  -9.465  -6.157  1.00  0.00           C  
ATOM    168  C   GLU A  11       9.045  -9.985  -7.561  1.00  0.00           C  
ATOM    169  O   GLU A  11       8.230 -10.888  -7.747  1.00  0.00           O  
ATOM    170  CB  GLU A  11       8.042  -9.131  -5.427  1.00  0.00           C  
ATOM    171  CG  GLU A  11       8.193  -9.105  -3.918  1.00  0.00           C  
ATOM    172  CD  GLU A  11       8.750 -10.404  -3.375  1.00  0.00           C  
ATOM    173  OE1 GLU A  11       8.036 -11.429  -3.411  1.00  0.00           O  
ATOM    174  OE2 GLU A  11       9.907 -10.407  -2.912  1.00  0.00           O  
ATOM    175  H   GLU A  11       9.786  -7.420  -6.405  1.00  0.00           H  
ATOM    176  HA  GLU A  11       9.868 -10.231  -5.606  1.00  0.00           H  
ATOM    177  HB2 GLU A  11       7.706  -8.155  -5.749  1.00  0.00           H  
ATOM    178  HB3 GLU A  11       7.294  -9.864  -5.683  1.00  0.00           H  
ATOM    179  HG2 GLU A  11       8.863  -8.301  -3.650  1.00  0.00           H  
ATOM    180  HG3 GLU A  11       7.224  -8.929  -3.474  1.00  0.00           H  
ATOM    181  N   GLY A  12       9.715  -9.409  -8.549  1.00  0.00           N 
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[9, 10])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"

name = ">F"
data = """
ATOM    551  N   PHE A  38       9.748   9.557  -7.966  1.00  0.00           N  
ATOM    552  CA  PHE A  38       8.331   9.287  -8.144  1.00  0.00           C  
ATOM    553  C   PHE A  38       7.968   9.437  -9.617  1.00  0.00           C  
ATOM    554  O   PHE A  38       8.714  10.055 -10.382  1.00  0.00           O  
ATOM    555  CB  PHE A  38       7.474  10.221  -7.275  1.00  0.00           C  
ATOM    556  CG  PHE A  38       7.721  11.686  -7.510  1.00  0.00           C  
ATOM    557  CD1 PHE A  38       8.683  12.366  -6.779  1.00  0.00           C  
ATOM    558  CD2 PHE A  38       6.988  12.385  -8.456  1.00  0.00           C  
ATOM    559  CE1 PHE A  38       8.909  13.712  -6.988  1.00  0.00           C  
ATOM    560  CE2 PHE A  38       7.210  13.731  -8.670  1.00  0.00           C  
ATOM    561  CZ  PHE A  38       8.172  14.395  -7.935  1.00  0.00           C  
ATOM    562  H   PHE A  38      10.199  10.134  -8.617  1.00  0.00           H  
ATOM    563  HA  PHE A  38       8.151   8.264  -7.846  1.00  0.00           H  
ATOM    564  HB2 PHE A  38       6.432  10.026  -7.477  1.00  0.00           H  
ATOM    565  HB3 PHE A  38       7.676  10.014  -6.234  1.00  0.00           H  
ATOM    566  HD1 PHE A  38       9.260  11.833  -6.037  1.00  0.00           H  
ATOM    567  HD2 PHE A  38       6.236  11.865  -9.033  1.00  0.00           H  
ATOM    568  HE1 PHE A  38       9.661  14.230  -6.412  1.00  0.00           H  
ATOM    569  HE2 PHE A  38       6.633  14.263  -9.411  1.00  0.00           H  
ATOM    570  HZ  PHE A  38       8.347  15.448  -8.100  1.00  0.00           H  
ATOM    571  N   TRP A  39       6.841   8.870 -10.018  1.00  0.00           N 
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[11, 12])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">G"
data = """
ATOM    181  N   GLY A  12       9.715  -9.409  -8.549  1.00  0.00           N  
ATOM    182  CA  GLY A  12       9.558  -9.859  -9.918  1.00  0.00           C  
ATOM    183  C   GLY A  12       8.331  -9.286 -10.597  1.00  0.00           C  
ATOM    184  O   GLY A  12       7.730  -9.941 -11.452  1.00  0.00           O  
ATOM    185  H   GLY A  12      10.332  -8.673  -8.347  1.00  0.00           H  
ATOM    186  HA2 GLY A  12      10.434  -9.568 -10.482  1.00  0.00           H  
ATOM    187  HA3 GLY A  12       9.486 -10.935  -9.922  1.00  0.00           H  
ATOM    188  N   SER A  13       7.941  -8.081 -10.210  1.00  0.00           N  
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">H"
data = """
ATOM    363  N   HIS A  25      12.595  10.712 -15.027  1.00  0.00           N  
ATOM    364  CA  HIS A  25      11.652  10.389 -13.962  1.00  0.00           C  
ATOM    365  C   HIS A  25      11.104   8.992 -14.231  1.00  0.00           C  
ATOM    366  O   HIS A  25       9.967   8.837 -14.677  1.00  0.00           O  
ATOM    367  CB  HIS A  25      10.487  11.389 -13.902  1.00  0.00           C  
ATOM    368  CG  HIS A  25      10.887  12.822 -13.718  1.00  0.00           C  
ATOM    369  ND1 HIS A  25      10.885  13.738 -14.749  1.00  0.00           N  
ATOM    370  CD2 HIS A  25      11.259  13.507 -12.611  1.00  0.00           C  
ATOM    371  CE1 HIS A  25      11.233  14.920 -14.285  1.00  0.00           C  
ATOM    372  NE2 HIS A  25      11.466  14.810 -12.989  1.00  0.00           N  
ATOM    373  H   HIS A  25      12.252  10.860 -15.941  1.00  0.00           H  
ATOM    374  HA  HIS A  25      12.183  10.389 -13.021  1.00  0.00           H  
ATOM    375  HB2 HIS A  25       9.926  11.326 -14.822  1.00  0.00           H  
ATOM    376  HB3 HIS A  25       9.840  11.121 -13.080  1.00  0.00           H  
ATOM    377  HD2 HIS A  25      11.372  13.102 -11.616  1.00  0.00           H  
ATOM    378  HE1 HIS A  25      11.312  15.828 -14.865  1.00  0.00           H  
ATOM    379  HE2 HIS A  25      11.545  15.569 -12.371  1.00  0.00           H  
ATOM    380  N   PRO A  26      11.916   7.956 -13.981  1.00  0.00           N   
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[10, 11])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">I"
data = """
ATOM     64  N   ILE A   4       2.838  -9.533   1.579  1.00  0.00           N  
ATOM     65  CA  ILE A   4       4.154  -9.068   1.177  1.00  0.00           C  
ATOM     66  C   ILE A   4       4.006  -8.060   0.037  1.00  0.00           C  
ATOM     67  O   ILE A   4       2.973  -8.035  -0.626  1.00  0.00           O  
ATOM     68  CB  ILE A   4       5.049 -10.259   0.746  1.00  0.00           C  
ATOM     69  CG1 ILE A   4       6.487  -9.803   0.480  1.00  0.00           C  
ATOM     70  CG2 ILE A   4       4.471 -10.949  -0.481  1.00  0.00           C  
ATOM     71  CD1 ILE A   4       7.423 -10.937   0.121  1.00  0.00           C  
ATOM     72  H   ILE A   4       2.368 -10.185   1.008  1.00  0.00           H  
ATOM     73  HA  ILE A   4       4.613  -8.578   2.026  1.00  0.00           H  
ATOM     74  HB  ILE A   4       5.056 -10.976   1.554  1.00  0.00           H  
ATOM     75 HG12 ILE A   4       6.489  -9.102  -0.340  1.00  0.00           H  
ATOM     76 HG13 ILE A   4       6.876  -9.318   1.364  1.00  0.00           H  
ATOM     77 HG21 ILE A   4       4.395 -10.237  -1.291  1.00  0.00           H  
ATOM     78 HG22 ILE A   4       3.489 -11.337  -0.249  1.00  0.00           H  
ATOM     79 HG23 ILE A   4       5.119 -11.761  -0.776  1.00  0.00           H  
ATOM     80 HD11 ILE A   4       8.413 -10.546  -0.054  1.00  0.00           H  
ATOM     81 HD12 ILE A   4       7.065 -11.428  -0.772  1.00  0.00           H  
ATOM     82 HD13 ILE A   4       7.455 -11.648   0.934  1.00  0.00           H  
ATOM     83  N   ASP A   5       5.017  -7.221  -0.172  1.00  0.00           N
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[8, 9])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">K"
data = """
ATOM    341  N   LYS A  24      14.038  10.288 -17.196  1.00  0.00           N  
ATOM    342  CA  LYS A  24      14.775  10.825 -16.061  1.00  0.00           C  
ATOM    343  C   LYS A  24      13.900  10.814 -14.812  1.00  0.00           C  
ATOM    344  O   LYS A  24      14.392  10.911 -13.686  1.00  0.00           O  
ATOM    345  CB  LYS A  24      15.289  12.234 -16.370  1.00  0.00           C  
ATOM    346  CG  LYS A  24      16.410  12.241 -17.400  1.00  0.00           C  
ATOM    347  CD  LYS A  24      16.919  13.644 -17.688  1.00  0.00           C  
ATOM    348  CE  LYS A  24      15.888  14.479 -18.428  1.00  0.00           C  
ATOM    349  NZ  LYS A  24      16.414  15.827 -18.774  1.00  0.00           N  
ATOM    350  H   LYS A  24      13.353  10.850 -17.631  1.00  0.00           H  
ATOM    351  HA  LYS A  24      15.624  10.177 -15.890  1.00  0.00           H  
ATOM    352  HB2 LYS A  24      14.472  12.830 -16.750  1.00  0.00           H  
ATOM    353  HB3 LYS A  24      15.660  12.681 -15.460  1.00  0.00           H  
ATOM    354  HG2 LYS A  24      17.229  11.645 -17.026  1.00  0.00           H  
ATOM    355  HG3 LYS A  24      16.040  11.807 -18.319  1.00  0.00           H  
ATOM    356  HD2 LYS A  24      17.152  14.130 -16.752  1.00  0.00           H  
ATOM    357  HD3 LYS A  24      17.813  13.574 -18.290  1.00  0.00           H  
ATOM    358  HE2 LYS A  24      15.611  13.964 -19.336  1.00  0.00           H  
ATOM    359  HE3 LYS A  24      15.018  14.591 -17.799  1.00  0.00           H  
ATOM    360  HZ1 LYS A  24      16.742  16.314 -17.915  1.00  0.00           H  
ATOM    361  HZ2 LYS A  24      15.665  16.401 -19.220  1.00  0.00           H  
ATOM    362  HZ3 LYS A  24      17.214  15.741 -19.440  1.00  0.00           H  
ATOM    363  N   HIS A  25      12.595  10.712 -15.027  1.00  0.00           N 
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[9, 10])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">L"
data = """
ATOM     95  N   LEU A   6       3.880  -7.102  -3.149  1.00  0.00           N  
ATOM     96  CA  LEU A   6       3.793  -7.841  -4.408  1.00  0.00           C  
ATOM     97  C   LEU A   6       4.047  -6.964  -5.632  1.00  0.00           C  
ATOM     98  O   LEU A   6       4.234  -7.478  -6.736  1.00  0.00           O  
ATOM     99  CB  LEU A   6       2.408  -8.476  -4.527  1.00  0.00           C  
ATOM    100  CG  LEU A   6       2.048  -9.469  -3.424  1.00  0.00           C  
ATOM    101  CD1 LEU A   6       0.544  -9.674  -3.368  1.00  0.00           C  
ATOM    102  CD2 LEU A   6       2.757 -10.794  -3.652  1.00  0.00           C  
ATOM    103  H   LEU A   6       3.049  -6.861  -2.681  1.00  0.00           H  
ATOM    104  HA  LEU A   6       4.529  -8.620  -4.381  1.00  0.00           H  
ATOM    105  HB2 LEU A   6       1.670  -7.686  -4.529  1.00  0.00           H  
ATOM    106  HB3 LEU A   6       2.356  -8.994  -5.474  1.00  0.00           H  
ATOM    107  HG  LEU A   6       2.368  -9.074  -2.471  1.00  0.00           H  
ATOM    108 HD11 LEU A   6       0.196 -10.059  -4.313  1.00  0.00           H  
ATOM    109 HD12 LEU A   6       0.060  -8.729  -3.162  1.00  0.00           H  
ATOM    110 HD13 LEU A   6       0.306 -10.376  -2.582  1.00  0.00           H  
ATOM    111 HD21 LEU A   6       2.484 -11.490  -2.874  1.00  0.00           H  
ATOM    112 HD22 LEU A   6       3.826 -10.637  -3.637  1.00  0.00           H  
ATOM    113 HD23 LEU A   6       2.468 -11.195  -4.613  1.00  0.00           H  
ATOM    114  N   CYS A   7       4.071  -5.657  -5.445  1.00  0.00           N  
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[8, 9])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">M"
data = """
ATOM      1  N   MET A   1       3.853 -12.003   8.216  1.00  0.00           N  
ATOM      2  CA  MET A   1       2.824 -13.067   8.165  1.00  0.00           C  
ATOM      3  C   MET A   1       1.859 -12.792   7.018  1.00  0.00           C  
ATOM      4  O   MET A   1       0.858 -12.094   7.186  1.00  0.00           O  
ATOM      5  CB  MET A   1       2.079 -13.169   9.507  1.00  0.00           C  
ATOM      6  CG  MET A   1       1.544 -11.843  10.035  1.00  0.00           C  
ATOM      7  SD  MET A   1       0.709 -12.011  11.626  1.00  0.00           S  
ATOM      8  CE  MET A   1      -0.653 -13.089  11.181  1.00  0.00           C  
ATOM      9  H   MET A   1       3.398 -11.069   8.321  1.00  0.00           H  
ATOM     10  HA  MET A   1       3.326 -14.003   7.971  1.00  0.00           H  
ATOM     11  HB2 MET A   1       1.244 -13.843   9.389  1.00  0.00           H  
ATOM     12  HB3 MET A   1       2.754 -13.576  10.246  1.00  0.00           H  
ATOM     13  HG2 MET A   1       2.370 -11.156  10.150  1.00  0.00           H  
ATOM     14  HG3 MET A   1       0.844 -11.442   9.317  1.00  0.00           H  
ATOM     15  HE1 MET A   1      -0.263 -14.016  10.787  1.00  0.00           H  
ATOM     16  HE2 MET A   1      -1.263 -12.607  10.432  1.00  0.00           H  
ATOM     17  HE3 MET A   1      -1.252 -13.294  12.056  1.00  0.00           H  
ATOM     18  N   ARG A   2       2.186 -13.337   5.845  1.00  0.00           N  
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[8, 9])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">N"
data = """
ATOM    423  N   ASN A  29      10.937   3.435 -10.701  1.00  0.00           N  
ATOM    424  CA  ASN A  29      11.006   3.843  -9.297  1.00  0.00           C  
ATOM    425  C   ASN A  29       9.786   3.347  -8.525  1.00  0.00           C  
ATOM    426  O   ASN A  29       9.858   3.111  -7.319  1.00  0.00           O  
ATOM    427  CB  ASN A  29      11.103   5.371  -9.172  1.00  0.00           C  
ATOM    428  CG  ASN A  29       9.746   6.065  -9.169  1.00  0.00           C  
ATOM    429  OD1 ASN A  29       9.136   6.255  -8.122  1.00  0.00           O  
ATOM    430  ND2 ASN A  29       9.273   6.460 -10.338  1.00  0.00           N  
ATOM    431  H   ASN A  29      10.846   4.123 -11.395  1.00  0.00           H  
ATOM    432  HA  ASN A  29      11.892   3.400  -8.868  1.00  0.00           H  
ATOM    433  HB2 ASN A  29      11.609   5.618  -8.252  1.00  0.00           H  
ATOM    434  HB3 ASN A  29      11.677   5.754 -10.003  1.00  0.00           H  
ATOM    435 HD21 ASN A  29       9.817   6.289 -11.143  1.00  0.00           H  
ATOM    436 HD22 ASN A  29       8.403   6.905 -10.357  1.00  0.00           H  
ATOM    437  N   ILE A  30       8.676   3.176  -9.239  1.00  0.00           N 
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[8, 9])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">P"
data = """
ATOM    380  N   PRO A  26      11.916   7.956 -13.981  1.00  0.00           N  
ATOM    381  CA  PRO A  26      11.639   6.602 -14.463  1.00  0.00           C  
ATOM    382  C   PRO A  26      10.571   5.871 -13.653  1.00  0.00           C  
ATOM    383  O   PRO A  26      10.675   5.761 -12.434  1.00  0.00           O  
ATOM    384  CB  PRO A  26      12.994   5.912 -14.319  1.00  0.00           C  
ATOM    385  CG  PRO A  26      13.650   6.594 -13.169  1.00  0.00           C  
ATOM    386  CD  PRO A  26      13.161   8.016 -13.187  1.00  0.00           C  
ATOM    387  HA  PRO A  26      11.351   6.609 -15.503  1.00  0.00           H  
ATOM    388  HB2 PRO A  26      12.846   4.860 -14.125  1.00  0.00           H  
ATOM    389  HB3 PRO A  26      13.564   6.040 -15.228  1.00  0.00           H  
ATOM    390  HG2 PRO A  26      13.362   6.113 -12.246  1.00  0.00           H  
ATOM    391  HG3 PRO A  26      14.723   6.566 -13.289  1.00  0.00           H  
ATOM    392  HD2 PRO A  26      12.960   8.357 -12.183  1.00  0.00           H  
ATOM    393  HD3 PRO A  26      13.890   8.657 -13.663  1.00  0.00           H  
ATOM    394  N   PRO A  27       9.528   5.355 -14.334  1.00  0.00           N 
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"

name = ">Q"
data = """
ATOM    695  N   GLN A  47      -4.146   8.764 -17.925  1.00  0.00           N  
ATOM    696  CA  GLN A  47      -2.847   8.133 -17.769  1.00  0.00           C  
ATOM    697  C   GLN A  47      -3.010   6.709 -17.255  1.00  0.00           C  
ATOM    698  O   GLN A  47      -3.476   6.487 -16.137  1.00  0.00           O  
ATOM    699  CB  GLN A  47      -1.971   8.966 -16.827  1.00  0.00           C  
ATOM    700  CG  GLN A  47      -2.649   9.325 -15.514  1.00  0.00           C  
ATOM    701  CD  GLN A  47      -1.888  10.373 -14.730  1.00  0.00           C  
ATOM    702  OE1 GLN A  47      -1.022  10.057 -13.921  1.00  0.00           O  
ATOM    703  NE2 GLN A  47      -2.205  11.634 -14.973  1.00  0.00           N  
ATOM    704  H   GLN A  47      -4.717   8.883 -17.134  1.00  0.00           H  
ATOM    705  HA  GLN A  47      -2.379   8.098 -18.742  1.00  0.00           H  
ATOM    706  HB2 GLN A  47      -1.073   8.409 -16.604  1.00  0.00           H  
ATOM    707  HB3 GLN A  47      -1.699   9.884 -17.327  1.00  0.00           H  
ATOM    708  HG2 GLN A  47      -3.638   9.705 -15.724  1.00  0.00           H  
ATOM    709  HG3 GLN A  47      -2.728   8.433 -14.910  1.00  0.00           H  
ATOM    710 HE21 GLN A  47      -2.905  11.817 -15.634  1.00  0.00           H  
ATOM    711 HE22 GLN A  47      -1.723  12.333 -14.485  1.00  0.00           H  
ATOM    712  N   GLU A  48      -2.652   5.742 -18.083  1.00  0.00           N
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[9, 10])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"

name = ">R"
data = """
ATOM    888  N   ARG A  58      -2.199 -13.511  -4.808  1.00  0.00           N  
ATOM    889  CA  ARG A  58      -3.466 -13.302  -4.132  1.00  0.00           C  
ATOM    890  C   ARG A  58      -3.442 -11.945  -3.442  1.00  0.00           C  
ATOM    891  O   ARG A  58      -2.748 -11.762  -2.444  1.00  0.00           O  
ATOM    892  CB  ARG A  58      -3.715 -14.424  -3.116  1.00  0.00           C  
ATOM    893  CG  ARG A  58      -5.186 -14.671  -2.809  1.00  0.00           C  
ATOM    894  CD  ARG A  58      -5.815 -13.542  -2.008  1.00  0.00           C  
ATOM    895  NE  ARG A  58      -7.257 -13.720  -1.861  1.00  0.00           N  
ATOM    896  CZ  ARG A  58      -7.914 -13.584  -0.710  1.00  0.00           C  
ATOM    897  NH1 ARG A  58      -7.258 -13.291   0.407  1.00  0.00           N  
ATOM    898  NH2 ARG A  58      -9.230 -13.742  -0.684  1.00  0.00           N  
ATOM    899  H   ARG A  58      -1.411 -13.773  -4.280  1.00  0.00           H  
ATOM    900  HA  ARG A  58      -4.251 -13.308  -4.875  1.00  0.00           H  
ATOM    901  HB2 ARG A  58      -3.295 -15.342  -3.502  1.00  0.00           H  
ATOM    902  HB3 ARG A  58      -3.215 -14.172  -2.193  1.00  0.00           H  
ATOM    903  HG2 ARG A  58      -5.723 -14.774  -3.740  1.00  0.00           H  
ATOM    904  HG3 ARG A  58      -5.272 -15.588  -2.243  1.00  0.00           H  
ATOM    905  HD2 ARG A  58      -5.361 -13.515  -1.029  1.00  0.00           H  
ATOM    906  HD3 ARG A  58      -5.628 -12.607  -2.518  1.00  0.00           H  
ATOM    907  HE  ARG A  58      -7.765 -13.955  -2.674  1.00  0.00           H  
ATOM    908 HH11 ARG A  58      -6.263 -13.169   0.393  1.00  0.00           H  
ATOM    909 HH12 ARG A  58      -7.754 -13.198   1.281  1.00  0.00           H  
ATOM    910 HH21 ARG A  58      -9.729 -13.965  -1.535  1.00  0.00           H  
ATOM    911 HH22 ARG A  58      -9.740 -13.648   0.184  1.00  0.00           H  
ATOM    912  N   ILE A  59      -4.188 -10.999  -3.990  1.00  0.00           N  
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[11, 12])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">S"
data = """
ATOM    144  N   SER A   9       9.601  -4.762  -8.201  1.00  0.00           N  
ATOM    145  CA  SER A   9      10.412  -5.540  -9.127  1.00  0.00           C  
ATOM    146  C   SER A   9      11.303  -6.492  -8.343  1.00  0.00           C  
ATOM    147  O   SER A   9      11.765  -7.506  -8.861  1.00  0.00           O  
ATOM    148  CB  SER A   9      11.263  -4.618 -10.005  1.00  0.00           C  
ATOM    149  OG  SER A   9      11.962  -5.353 -10.996  1.00  0.00           O  
ATOM    150  H   SER A   9       9.825  -3.815  -8.048  1.00  0.00           H  
ATOM    151  HA  SER A   9       9.748  -6.115  -9.753  1.00  0.00           H  
ATOM    152  HB2 SER A   9      10.622  -3.900 -10.496  1.00  0.00           H  
ATOM    153  HB3 SER A   9      11.979  -4.097  -9.388  1.00  0.00           H  
ATOM    154  HG  SER A   9      11.334  -5.918 -11.481  1.00  0.00           H  
ATOM    155  N   SER A  10      11.525  -6.167  -7.076  1.00  0.00           N 
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[6, 7])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">T"
data = """
ATOM    278  N   THR A  19       1.661   5.108 -20.510  1.00  0.00           N  
ATOM    279  CA  THR A  19       1.268   6.131 -21.476  1.00  0.00           C  
ATOM    280  C   THR A  19       2.233   7.313 -21.437  1.00  0.00           C  
ATOM    281  O   THR A  19       2.784   7.720 -22.460  1.00  0.00           O  
ATOM    282  CB  THR A  19      -0.161   6.635 -21.199  1.00  0.00           C  
ATOM    283  OG1 THR A  19      -1.039   5.525 -20.990  1.00  0.00           O  
ATOM    284  CG2 THR A  19      -0.676   7.480 -22.356  1.00  0.00           C  
ATOM    285  H   THR A  19       1.007   4.809 -19.835  1.00  0.00           H  
ATOM    286  HA  THR A  19       1.290   5.688 -22.461  1.00  0.00           H  
ATOM    287  HB  THR A  19      -0.146   7.245 -20.308  1.00  0.00           H  
ATOM    288  HG1 THR A  19      -0.543   4.700 -21.059  1.00  0.00           H  
ATOM    289 HG21 THR A  19      -1.669   7.838 -22.127  1.00  0.00           H  
ATOM    290 HG22 THR A  19      -0.708   6.881 -23.254  1.00  0.00           H  
ATOM    291 HG23 THR A  19      -0.017   8.322 -22.509  1.00  0.00           H  
ATOM    292  N   SER A  20       2.434   7.853 -20.249  1.00  0.00           N 
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[7, 8])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">V"
data = """
ATOM    214  N   VAL A  15       6.732  -2.888 -14.688  1.00  0.00           N  
ATOM    215  CA  VAL A  15       5.686  -2.035 -15.234  1.00  0.00           C  
ATOM    216  C   VAL A  15       5.743  -1.996 -16.762  1.00  0.00           C  
ATOM    217  O   VAL A  15       6.813  -2.107 -17.362  1.00  0.00           O  
ATOM    218  CB  VAL A  15       5.764  -0.597 -14.671  1.00  0.00           C  
ATOM    219  CG1 VAL A  15       5.355  -0.578 -13.205  1.00  0.00           C  
ATOM    220  CG2 VAL A  15       7.162  -0.019 -14.839  1.00  0.00           C  
ATOM    221  H   VAL A  15       7.668  -2.613 -14.778  1.00  0.00           H  
ATOM    222  HA  VAL A  15       4.737  -2.457 -14.939  1.00  0.00           H  
ATOM    223  HB  VAL A  15       5.071   0.021 -15.223  1.00  0.00           H  
ATOM    224 HG11 VAL A  15       5.396   0.436 -12.834  1.00  0.00           H  
ATOM    225 HG12 VAL A  15       6.032  -1.197 -12.635  1.00  0.00           H  
ATOM    226 HG13 VAL A  15       4.349  -0.957 -13.106  1.00  0.00           H  
ATOM    227 HG21 VAL A  15       7.413   0.019 -15.887  1.00  0.00           H  
ATOM    228 HG22 VAL A  15       7.875  -0.645 -14.321  1.00  0.00           H  
ATOM    229 HG23 VAL A  15       7.191   0.978 -14.424  1.00  0.00           H  
ATOM    230  N   ILE A  16       4.578  -1.860 -17.372  1.00  0.00           N 
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[7, 8])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">W"
data = """
ATOM   1404  N   TRP A  88     -12.161  -3.632 -11.787  1.00  0.00           N  
ATOM   1405  CA  TRP A  88     -12.061  -4.140 -10.431  1.00  0.00           C  
ATOM   1406  C   TRP A  88     -12.609  -3.126  -9.438  1.00  0.00           C  
ATOM   1407  O   TRP A  88     -13.656  -3.344  -8.829  1.00  0.00           O  
ATOM   1408  CB  TRP A  88     -10.600  -4.476 -10.104  1.00  0.00           C  
ATOM   1409  CG  TRP A  88     -10.425  -5.296  -8.860  1.00  0.00           C  
ATOM   1410  CD1 TRP A  88     -11.402  -5.916  -8.134  1.00  0.00           C  
ATOM   1411  CD2 TRP A  88      -9.187  -5.602  -8.208  1.00  0.00           C  
ATOM   1412  NE1 TRP A  88     -10.849  -6.581  -7.069  1.00  0.00           N  
ATOM   1413  CE2 TRP A  88      -9.491  -6.406  -7.093  1.00  0.00           C  
ATOM   1414  CE3 TRP A  88      -7.853  -5.274  -8.459  1.00  0.00           C  
ATOM   1415  CZ2 TRP A  88      -8.508  -6.885  -6.231  1.00  0.00           C  
ATOM   1416  CZ3 TRP A  88      -6.878  -5.751  -7.603  1.00  0.00           C  
ATOM   1417  CH2 TRP A  88      -7.209  -6.550  -6.502  1.00  0.00           C  
ATOM   1418  H   TRP A  88     -11.460  -3.034 -12.129  1.00  0.00           H  
ATOM   1419  HA  TRP A  88     -12.651  -5.041 -10.370  1.00  0.00           H  
ATOM   1420  HB2 TRP A  88     -10.175  -5.030 -10.926  1.00  0.00           H  
ATOM   1421  HB3 TRP A  88     -10.049  -3.555  -9.976  1.00  0.00           H  
ATOM   1422  HD1 TRP A  88     -12.455  -5.875  -8.372  1.00  0.00           H  
ATOM   1423  HE1 TRP A  88     -11.350  -7.102  -6.398  1.00  0.00           H  
ATOM   1424  HE3 TRP A  88      -7.577  -4.658  -9.303  1.00  0.00           H  
ATOM   1425  HZ2 TRP A  88      -8.747  -7.503  -5.379  1.00  0.00           H  
ATOM   1426  HZ3 TRP A  88      -5.840  -5.507  -7.782  1.00  0.00           H  
ATOM   1427  HH2 TRP A  88      -6.416  -6.900  -5.860  1.00  0.00           H  
ATOM   1428  N   ILE A  89     -11.900  -2.020  -9.285  1.00  0.00           N
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[14, 15])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


name = ">Y"
data = """
ATOM   1785  N   TYR A 113      -2.007  -9.290 -13.930  1.00  0.00           N  
ATOM   1786  CA  TYR A 113      -1.769  -7.964 -14.484  1.00  0.00           C  
ATOM   1787  C   TYR A 113      -2.772  -6.972 -13.901  1.00  0.00           C  
ATOM   1788  O   TYR A 113      -3.932  -7.319 -13.674  1.00  0.00           O  
ATOM   1789  CB  TYR A 113      -1.904  -7.987 -16.015  1.00  0.00           C  
ATOM   1790  CG  TYR A 113      -1.004  -8.993 -16.706  1.00  0.00           C  
ATOM   1791  CD1 TYR A 113       0.328  -8.702 -16.962  1.00  0.00           C  
ATOM   1792  CD2 TYR A 113      -1.490 -10.236 -17.097  1.00  0.00           C  
ATOM   1793  CE1 TYR A 113       1.153  -9.619 -17.588  1.00  0.00           C  
ATOM   1794  CE2 TYR A 113      -0.671 -11.159 -17.721  1.00  0.00           C  
ATOM   1795  CZ  TYR A 113       0.648 -10.845 -17.963  1.00  0.00           C  
ATOM   1796  OH  TYR A 113       1.473 -11.762 -18.579  1.00  0.00           O  
ATOM   1797  H   TYR A 113      -2.935  -9.606 -13.854  1.00  0.00           H  
ATOM   1798  HA  TYR A 113      -0.768  -7.659 -14.216  1.00  0.00           H  
ATOM   1799  HB2 TYR A 113      -2.923  -8.227 -16.274  1.00  0.00           H  
ATOM   1800  HB3 TYR A 113      -1.664  -7.007 -16.402  1.00  0.00           H  
ATOM   1801  HD1 TYR A 113       0.721  -7.741 -16.670  1.00  0.00           H  
ATOM   1802  HD2 TYR A 113      -2.525 -10.479 -16.907  1.00  0.00           H  
ATOM   1803  HE1 TYR A 113       2.186  -9.373 -17.777  1.00  0.00           H  
ATOM   1804  HE2 TYR A 113      -1.068 -12.118 -18.018  1.00  0.00           H  
ATOM   1805  HH  TYR A 113       0.988 -12.210 -19.287  1.00  0.00           H  
ATOM   1806  N   LEU A 114      -2.329  -5.750 -13.643  1.00  0.00           N 
"""
regex = r'ATOM\s+\d+\s+\w+\s+\w+\s+\w+\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)'
matches = re.findall(regex, data)
points = np.array([[float(x), float(y), float(z)] for x, y, z in matches])
np.set_printoptions(precision=3, suppress=True)
transformed_points = transform_points(points, point_a_index=0, plane_points_indices=[1, 2, 3],  side_chain_indices=[4, 5], side_remove_indices=[12, 13])
modified_pdb = modify_pdb_data(data, transformed_points)
print(name)
print(modified_pdb)
res += name + "\n" + modified_pdb + "\n"


with open("amino_acids.txt", "w") as pdb_file:
    pdb_file.write(res)