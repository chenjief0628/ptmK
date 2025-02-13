from Bio.PDB import PDBParser, is_aa
import numpy as np
import os
import pandas as pd
from Bio.PDB import DSSP
import re
import freesasa
import argparse


print("calculate activity distance and SASA............")

def get_lys_residues(model):
    lys_residues = {}
    for chain in model:
        for residue in chain:
            if residue.get_resname() == 'LYS':
                res_id = residue.get_id()
                key = (chain.id, res_id[1], res_id[2].strip())
                if 'CA' in residue:
                    lys_residues[key] = residue['CA'].get_coord()
                else:
                    print("There is no CA atom in " + chain.id + " " + res_id[1])
    return lys_residues


def calculate_centroid(model, residue_ids):
    ca_coords = []
    for chain in model:
        chain_id = chain.get_id()
    
    for residue_id in residue_ids:
        residue = model[chain_id][(' ', int(residue_id) -1 , ' ')]
        # 获取 CA 原子
        ca = residue['CA']
        coord = ca.get_coord()
        ca_coords.append(coord)
        print(ca_coords)
    
    centroid = np.mean(ca_coords, axis=0)
    return centroid

def cal_dist(dict_k, zx):
    list_d = []
    for key in dict_k:
        c1 = dict_k[key]
        c2 = zx
        d = np.linalg.norm(c1 - c2)
        list_d.append(d)
    return list_d

def cal_sasa(pdb_file,model):
    aim_list = []
    for chain in model:
        for residue in chain:
            if residue.get_resname() == 'LYS':
                res_id = residue.get_id()
                aim = res_id[1]
                aim_list.append(aim)
            
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)
    result = result.residueAreas()['A']
    vlist = []
    for aim in aim_list:
        v = result[str(aim)].total
        vlist.append(v)
    return vlist

def main():
    parser = argparse.ArgumentParser(description='clean pdb file')
    parser.add_argument('pdb_file', help='PDB file')
    parser.add_argument('--residue_ids', help='residue ids')
    args = parser.parse_args()
    pdb_file = args.pdb_file
    residue_ids = args.residue_ids
    residue_ids = residue_ids.split(',')
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    model = structure[0]
    ds = calculate_centroid(model, residue_ids)
    dict_k = get_lys_residues(model)
    d = cal_dist(dict_k,ds)
    vl = cal_sasa(pdb_file, model)
    list_k = []
    for k in dict_k:
        num = k[1]
        list_k.append(num)
    score = [1/a * 0.7 + b/max(vl) * 0.3 for a, b in zip(d, vl)]
    print(score)
    data = {'K num': list_k, 'Activate distance': d, 'SASA area':vl, 'score' : score}
    df = pd.DataFrame(data)
    df.to_csv("cal_distance_sasa.csv", index = False)

main()