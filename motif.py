import os
from Bio import PDB
import numpy as np
import pandas as pd
import argparse
from Bio.PDB import PDBParser, PDBIO

def get_residue(pdb_path):
    """
    Get the residue coordinates list for pdb file
    """
    parser = PDB.PDBParser()
    structure = parser.get_structure('PDB', pdb_path)
    model = structure[0]
    for chain in model:
        for residue in chain:
            atom_coordinates = {}
            for atom in residue:
                atom1 = atom.get_name()
                atom_coordinates[atom1] = atom.get_coord()
    return atom_coordinates

def center_coordinates(coordinates):
    """
    Get the center coordinates of the residue
    """
    center = np.mean(coordinates, axis=0)
    return coordinates - center


def move_coord(motif_coords, pdb_path, nums):
    """
    Align the motif coordinates to the pdb coordinates
    """
    parser = PDB.PDBParser()
    structure = parser.get_structure('PDB', pdb_path)
    model = structure[0]
    for chain in model:
        for residue in chain:
            if residue.get_id()[1] == nums:
                atom_coordinates = {}
                for atom in residue:
                    atom1 = atom.get_name()
                    atom_coordinates[atom1] = atom.get_coord()
    pdb_coords = atom_coordinates
    overlap = ['N', 'CA', 'C', 'O']
    coords_dict1 = [pdb_coords[atom] for atom in overlap if atom in pdb_coords]
    coords_dict2 = [motif_coords[atom] for atom in overlap if atom in motif_coords]
    coords_dict1 = np.array(coords_dict1)
    coords_dict2 = np.array(coords_dict2)
    centroid1 = np.mean(coords_dict1, axis=0)
    centroid2 = np.mean(coords_dict2, axis=0)
    coords_dict1_translated = coords_dict1 - centroid1
    coords_dict2_translated = coords_dict2 - centroid2
    H = np.dot(coords_dict2_translated.T, coords_dict1_translated)
    U, _, Vt = np.linalg.svd(H)
    rotation_matrix = np.dot(Vt.T, U.T)
    transformed_coords = {}
    for atom, coord in motif_coords.items():
        # 平移
        coord_translated = coord - centroid2
        # 旋转
        coord_rotated = np.dot(rotation_matrix, coord_translated.T) + centroid1
        transformed_coords[atom] = coord_rotated

    return transformed_coords

def save_to_pdb(coords_dict, pdb_path, motif_path, nums):
    """
    Save the transformed coordinates to a new pdb file
    """
    parser = PDB.PDBParser()
    structure = parser.get_structure('PDB', motif_path)
    model = structure[0]
    for chain in model:
        for residue in chain:
            residue.id = (' ', nums, ' ')
            for atom_name, new_coord in coords_dict.items():
                atom = residue[atom_name]
                atom.coord = np.array(new_coord)
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save('temp.pdb')
    structure1 = parser.get_structure('PDB', pdb_path)
    model1 = structure1[0]
    for chain in model1:
        residues = list(chain)
        residue_to_delete = residues[nums - 1]
        chain.detach_child(residue_to_delete.get_id())
        residues = list(chain)
        residues.insert(nums -1, residue)
    
    aaa = list(chain)
    for i in aaa:
        chain.detach_child(i.get_id())

    for i in residues:
        print(i)
        chain.add(i)
    io = PDB.PDBIO()
    io.set_structure(structure1)
    io.save('motif.result.pdb')

def main():
    parser = argparse.ArgumentParser(description='Align motif to pdb')
    parser.add_argument('--motif_path', type=str, help='path to motif pdb file')
    parser.add_argument('--pdb_path', type=str, help='path to pdb file')
    parser.add_argument('--nums', type=int, help='residue number')
    args = parser.parse_args()
    motif_path = args.motif_path
    pdb_path = args.pdb_path
    nums = args.nums
    motif_coords = get_residue(motif_path)
    transformed = move_coord(motif_coords, pdb_path, nums)
    print(transformed)
    save_to_pdb(transformed, pdb_path, motif_path,nums)

if __name__ == '__main__':
    main()


