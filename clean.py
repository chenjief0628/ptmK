import os
from Bio import PDB
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import argparse


class NonAminoAcidRemover(PDB.Select):
    def accept_residue(self, residue):
        # 检查残基是否是氨基酸
        if residue.get_id()[0] == ' ':
            return 1  # 是氨基酸，保留
        return 0  # 不是氨基酸，删除


def process_ent_file(ent_file_path, output_pdb_path):
    # 1. 使用 BioPython 读取 .ent 文件，提取 A 链
    parser = PDB.PPBuilder()
    
    # 读取 .ent 文件并解析
    structure = PDB.PDBParser(QUIET=True).get_structure('protein', ent_file_path)
    model = structure[0]
    chains = []
    for model in structure:
        for chain in model:
            id = chain.get_id()
            chains.append(id)
    if len(chains) > 1:
        print("chain id are " + ",".join(chains) + " and first chain was selected.")
    else:
        print("chan id is " + chains[0])
    a_chain = model[chains[0]]
    
    with open("temp.pdb",'w') as f:
        io = PDB.PDBIO()
        io.set_structure(a_chain)
        io.save(f,NonAminoAcidRemover())
        
    fixer = PDBFixer(filename= "temp.pdb")
    
    # 修复 A 链的缺失原子
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)  # 添加氢原子，使用7.0的pH值
    
    # 3. 将修复后的结构保存为 PDB 文件
    # 保存修复后的 PDB 文件
    with open(output_pdb_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"A 链已修复并保存为 {output_pdb_path}")

def main():
    print("aaa")
    parser = argparse.ArgumentParser(description='clean pdb file')
    parser.add_argument('--pdb_path', type=str, help='path to pdb file')
    parser.add_argument('--out_path', type=str, help='path to output pdb file')
    args = parser.parse_args()
    pdb_path = args.pdb_path
    output_pdb_path = args.out_path
    process_ent_file(pdb_path, output_pdb_path)


main()